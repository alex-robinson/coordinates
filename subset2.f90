

module subset2

    use coordinates
!     use interp2D
    use ncio 

    implicit none

    type subset_class

        type(points_class) :: pts
        type(grid_class)   :: grid 
        type(map_class)    :: map_tosub_grid
        type(map_class)    :: map_tosub, map_fromsub
        integer            :: npts, factor
        logical            :: subset
        double precision, allocatable :: var2D(:,:) 
        integer, dimension(:,:), allocatable :: mask2D 
        logical, dimension(:,:), allocatable :: mask_pack  
    end type 

    interface subset_to_grid  
        module procedure subset_to_grid_double, subset_to_grid_integer
    end interface

    interface subset_to_points  
        module procedure subset_to_points_double, subset_to_points_integer
    end interface

    private
    public :: subset_class 
    public :: subset_init, subset_redefine
    public :: subset_to_grid, subset_to_points
    public :: subset_gen_mask 

contains 

    subroutine subset_init(sub,grid,npts,factor,max_neighbors,lat_lim,load)
        ! Determine the coordinates of the domain

        implicit none 

        type(subset_class) :: sub    ! The subset definition
        type(grid_class)   :: grid   ! The original grid 
        integer :: npts              ! Number of points to allocate to the subset
        integer :: factor            ! Resolution factor (should be >= 1)
        integer, optional :: max_neighbors ! Maximum number of neighbors to use for mapping
        integer :: max_neighbs(2)
        double precision, optional :: lat_lim
        double precision :: lat_limit 
        logical, optional :: load 
        logical :: load_map

        character(len=12) :: suffix, suffix2 
        double precision, allocatable, dimension(:) :: x, y

        sub%subset = .TRUE.
        sub%factor = factor  
        if (npts .le. 0) then 
            sub%subset = .FALSE. 
            sub%factor = 1 
        end if 

        ! Define a suffix to append to the grid name
        if (sub%factor .ge. 10) then
            write(suffix, "(a4,i2)") "-ghi", sub%factor
            write(suffix2,"(a4,i2)") "-phi", sub%factor
        else if (sub%factor .ge. 1) then
            write(suffix, "(a4,i1)") "-ghi", sub%factor 
            write(suffix2,"(a4,i1)") "-phi", sub%factor
        else
            write(*,*) "subset_define:: Error: factor must be greater than or equal to one."
            stop 
        end if 

        max_neighbs = [6,9] 
        if (present(max_neighbors)) max_neighbs(2) = max_neighbors

        lat_limit   = 4.d0 
        if (present(lat_lim)) lat_limit = lat_lim 

        load_map = .FALSE. 
        if (present(load)) load_map = load 

        if (sub%subset) then 
            ! Initialize the new grid with input grid characteristics
            ! but with the new resolution
            call grid_init(sub%grid,name=trim(grid%name)//trim(suffix),mtype=grid%mtype, &
                           units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                           x0=grid%G%x(1),dx=grid%G%dx/dble(sub%factor),nx=(grid%G%nx-1)*sub%factor+1, &
                           y0=grid%G%y(1),dy=grid%G%dy/dble(sub%factor),ny=(grid%G%ny-1)*sub%factor+1, &
                           lambda=grid%proj%lambda,phi=grid%proj%phi,alpha=grid%proj%alpha, &
                           x_e=grid%proj%x_e,y_n=grid%proj%y_n)
            
            ! Make sure the subset npts is consistent with the new grid
            sub%npts = npts 
            if (sub%npts .gt. sub%grid%npts) sub%npts = sub%grid%npts 

            ! Allocate temporary x,y vectors to store points of interest
            ! for new coordinates definition
            if (allocated(x)) deallocate(x)
            if (allocated(y)) deallocate(y)
            allocate(x(sub%npts),y(sub%npts))

            ! For now, assign dummy x/y values with correct length.
            ! The actual values will be determined later from 
            ! packing/unpacking at each time step
            x = 0.d0
            x(1:2) = [grid%G%x(1),grid%G%x(grid%G%nx)]
            y = 0.d0 
            y(1:2) = [grid%G%y(1),grid%G%y(grid%G%ny)]

            ! Initialize the points class, which is a subset of points of the grid class
            call points_init(sub%pts,trim(grid%name)//trim(suffix2),mtype=grid%mtype, &
                             units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                             x=x,y=y, &
                             lambda=grid%proj%lambda,phi=grid%proj%phi, &
                             alpha=grid%proj%alpha,x_e=grid%proj%x_e,y_n=grid%proj%y_n)

            ! Initialize mapping to subset grid from input grid, used for mask_pack generation
            if (sub%factor .gt. 1) then 
                call map_init(sub%map_tosub_grid,grid,sub%grid, &
                              max_neighbors=max_neighbs(1),lat_lim=lat_limit,fldr="maps",load=load_map)
            end if
        else
            ! Intialize the sub grid & pts with the input grid characteristics
            sub%grid = grid 
            sub%npts = sub%grid%npts 
            call grid_to_points(sub%grid,sub%pts)

        end if 

        ! Initialize the self-filling map 

        ! Allocate the packing mask for later use 
        if (allocated(sub%mask_pack)) deallocate(sub%mask_pack)
        allocate(sub%mask_pack(sub%grid%G%nx,sub%grid%G%ny))
        sub%mask_pack = .TRUE. 
        
        ! Initialize arrays to handle data and interpolation mask on the grid 
        ! (to avoid reallocating continuously) 
        if (allocated(sub%var2D)) deallocate(sub%var2D)
        allocate(sub%var2D(sub%grid%G%nx,sub%grid%G%ny)) 
        if (allocated(sub%mask2D)) deallocate(sub%mask2D)
        allocate(sub%mask2D(sub%grid%G%nx,sub%grid%G%ny)) 

        write(*,"(a,i10,1x,a1,1x,i10)") "subset:: subset_init :: Initialized subset, npts = ", &
                                       sub%npts, "/",sub%grid%npts
        return

    end subroutine subset_init 

    subroutine subset_redefine(sub,grid,mask_pack,max_neighbors,lat_lim,load)
        ! Re-determine the coordinates of the current domain,
        ! which may be a subset of points from a grid
        ! This function allows the surface calculations to adapt to a changing 
        ! topography (eg, set increased density of points where topography is steep)

        implicit none 

        type(subset_class) :: sub
        type(grid_class)   :: grid   ! The original grid 
        logical, dimension(:,:) :: mask_pack 
        integer, optional :: max_neighbors ! Maximum number of neighbors to use for mapping
        integer :: max_neighbs(2)
        double precision, optional :: lat_lim
        double precision :: lat_limit 
        logical, optional :: load
        logical :: load_map

        max_neighbs = [6,9] 
        if (present(max_neighbors)) max_neighbs(2) = max_neighbors

        lat_limit   = 4.d0 
        if (present(lat_lim)) lat_limit = lat_lim 

        load_map = .FALSE.
        if (present(load)) load_map = load

        if (sub%subset) then 

            ! Check that mask_pack is consistent with subset of npts
            if (count(mask_pack) .ne. sub%npts) then
                write(*,"(a)") "subset:: subset_redefine:: Error: packing mask must specify the same "// &
                           "number of points as defined in the subset."
                write(*,*) "subset npts =",sub%npts
                write(*,*) "mask_pack total = ",count(mask_pack)
                write(*,*) "Maybe the packing mask has not been initialized well?"
                stop 
            end if 

            write(*,*) "subset_redefine::", sub%grid%npts, sub%pts%npts 

            ! Get a new subset of points from the grid
            call grid_to_points(sub%grid,sub%pts,mask_pack=mask_pack,define=.FALSE.)
            
            ! Initialize to and fro mappings for subset pts and input grid
            ! (map generation can take some time - should be called infrequently)
            if (sub%factor .gt. 1) then 

                call map_init(sub%map_tosub,grid,sub%pts, &
                              max_neighbors=max_neighbs(1),lat_lim=lat_limit,fldr="maps",load=load_map) 

                call map_init(sub%map_fromsub,sub%pts,grid, &
                              max_neighbors=max_neighbs(2),lat_lim=lat_limit,fldr="maps",load=load_map)
            end if

        end if 

        return

    end subroutine subset_redefine 

    subroutine subset_to_grid_double(sub,var1D,var2D,mask_pack,map, &
                                     method,radius,fill,border,missing_value)
        ! This subroutine maps a subset of points (var1D) onto
        ! a 2D array (var2D) of resolution grid. 
        ! The subset should already be initialized. 

        implicit none 
 
        type(subset_class), intent(INOUT), target  :: sub 
        type(map_class),    intent(IN),    target, optional :: map 
        double precision, intent(OUT)   :: var2D(:,:)
        double precision, intent(IN)    :: var1D(:)
        logical, intent(IN)             :: mask_pack(:,:)

        double precision, allocatable :: var2Dtmp(:,:) 
        integer, allocatable          :: mask2D(:,:)
        character(len=*)           :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border

        type(map_class), pointer :: map_local 

        ! Assign a missing_value for use with mapping routine
        missing_val = -9999.d0
        if (present(missing_value)) missing_val = missing_value 

        if ( (sub%subset .and. sub%factor .gt. 1) .or. present(map)) then 

!             ! Consistency check 
!             if (count(mask_pack) .ne. sub%npts) then 
!                 write(*,*) "subset_to_points:: Error: total masked points not equal to npts."
!                 write(*,*) "count(mask_pack) =", count(mask_pack)
!                 write(*,*) "sub%npts =", sub%npts
!                 write(*,*) "Make sure mask_pack has been properly generated."
!                 stop 
!             end if 

            ! Determine map to use here 
            map_local => sub%map_fromsub 
            if (present(map)) map_local => map 

!             write(*,*) trim(map_local%name1), " => ",trim(map_local%name2)

            ! Step 1: Unpack the 1D variable onto its corresponding 
            ! predefined 2D grid.
            sub%var2D = missing_val 
            sub%var2D = unpack(var1D,mask_pack,sub%var2D)  ! Unpack the 1D vector 

            ! Step 2: Map the temporary 2D array to the desired 2D resolution
            allocate(mask2D(map_local%G%nx,map_local%G%ny))

            call map_field(map_local,"Mapped variable",var1D,var2D,mask2D,method=method, &
                           radius=radius,fill=fill,border=border,missing_value=missing_val)

        else if (sub%subset) then 

            ! Step 1: unpack the 1D vector into the 2D array using mask
            allocate(var2Dtmp(size(var2D,1),size(var2D,2)))
            var2Dtmp = missing_val 
            var2D = unpack(var1D,mask_pack,var2Dtmp)

        else

            ! Step 1: reshape the 1D vector into the 2D array
            var2D = reshape(var1D,[size(var2D,1),size(var2D,2)])

        end if 

!         if (present(fill) .and. method .ne. "nn") then 
!             if (fill) call fill_weighted(var2D,missing_val)
!         end if 
        
        return

    end subroutine subset_to_grid_double

    subroutine subset_to_grid_integer(sub,var1D,var2D,mask_pack,map, &
                                      radius,fill,border,missing_value)
        ! This subroutine maps a subset of points (var1D) onto
        ! a 2D array (var2D) of resolution grid. 
        ! The subset should already be initialized.

        implicit none 
 
        type(subset_class), intent(INOUT)  :: sub 
        type(map_class), intent(IN), optional :: map 
        integer, intent(OUT)   :: var2D(:,:)
        integer, intent(IN)    :: var1D(:)
        logical, intent(IN)    :: mask_pack(:,:)

        double precision, allocatable :: var2Dtmp(:,:)

        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border

        allocate(var2Dtmp(size(var2D,1),size(var2D,2)))

        call subset_to_grid_double(sub,dble(var1D),var2Dtmp,mask_pack,map, &
                                   "nn",radius,fill,border,missing_value)
        var2D = nint(var2Dtmp)

        return 

    end subroutine subset_to_grid_integer 

    subroutine subset_to_points_double(sub,var2D,var1D,mask_pack,map, &
                                       method,radius,fill,border,missing_value)
        ! This subroutine maps a 2D array (var2D) onto
        ! a subset of points (var1D) of resolution sub%grid. 
        ! The subset should already be initialized.

        implicit none 
 
        type(subset_class), intent(INOUT), target     :: sub 
        type(map_class),    intent(IN),    target, optional :: map 
        double precision, intent(IN)    :: var2D(:,:)
        double precision, intent(OUT)   :: var1D(:)
        logical, intent(IN)             :: mask_pack(:,:)

        character(len=*)           :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border 

        type(map_class), pointer :: map_local 

        if ( (sub%subset .and. sub%factor .gt. 1) .or. present(map)) then 

!             ! Consistency check 
!             if (count(mask_pack) .ne. sub%npts) then 
!                 write(*,*) "subset_to_points:: Error: total masked points not equal to npts."
!                 write(*,*) "count(mask_pack) =", count(mask_pack)
!                 write(*,*) "sub%npts =", sub%npts
!                 write(*,*) "Make sure mask_pack has been properly generated."
!                 stop 
!             end if 

            ! Determine map to use here 
            map_local => sub%map_tosub 
            if (present(map)) map_local => map 

!             write(*,*) trim(map_local%name1), " => ",trim(map_local%name2)
!             write(*,*) "size(var2D): ",size(var2D,1), size(var2D,2)
!             write(*,*) "sub%npts: ", sub%npts 
!             write(*,*) "map%npts: ", map_local%npts 

            ! Assign a missing_value for use with mapping routine
            missing_val = -9999.d0
            if (present(missing_value)) missing_val = missing_value 

!             ! Step 1: Prefill the predefined 2D array of the subset
!             sub%var2D = missing_val
            var1D = missing_val 

            ! Step 2: Map the 2D array to the temporary 2D array of the subset
            call map_field(map_local,"Mapped variable",var2D,var1D,method=method, &
                           radius=radius,fill=fill,border=border,missing_value=missing_val)

!             ! Step 3: Pack the 2D variable onto its corresponding predefined 1D points.
!             var1D = pack(sub%var2D,mask_pack)
        
        else if (sub%subset) then 

            ! Step 1: Pack the 2D variable onto its corresponding predefined 1D points.
            var1D = pack(var2D,mask_pack)
        
        else

            ! Step 1: Reshape the 2D variable into the 1D vector.
            var1D = reshape(var2D,[size(var1D)])

        end if 
        
        return

    end subroutine subset_to_points_double

    subroutine subset_to_points_integer(sub,var2D,var1D,mask_pack,map, &
                                        radius,fill,border,missing_value)
        ! This subroutine maps a 2D array (var2D) onto
        ! a subset of points (var1D) of resolution sub%grid. 
        ! The subset should already be initialized.
        ! Note: currently var1D and var2D are double precision! 

        implicit none 
 
        type(subset_class), intent(INOUT)  :: sub 
        type(map_class), intent(IN), optional :: map 
        integer, intent(IN)    :: var2D(:,:)
        integer, intent(OUT)   :: var1D(:)
        logical, intent(IN)    :: mask_pack(:,:)

        double precision, allocatable :: var1Dtmp(:) 

        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border 

        allocate(var1Dtmp(size(var1D)))

        call subset_to_points_double(sub,dble(var2D),var1Dtmp,mask_pack,map, &
                                     "nn",radius,fill,border,missing_value)
        var1D = nint(var1Dtmp)

        return 

    end subroutine subset_to_points_integer




    subroutine subset_gen_mask(mask_pack,var,npts,min_spacing,method,map)
        ! This routine can be used in a generic way to generate a mask
        ! based on an input field, a minimum spacing and a function (min/max)
        ! Note: it is highly relevant to subset, although it doesn't 
        ! depend on the subset_class since the mask doesn't have to
        ! be the subset%mask_pack variable.
        implicit none 

        logical, dimension(:,:), intent(INOUT) :: mask_pack 
        double precision, dimension(:,:), intent(IN) :: var 
        type(map_class), intent(IN), optional :: map
        integer :: npts, min_spacing 
        character(len=*) :: method 

        integer :: nx, ny, npts0, loc(2), i,j
        double precision, dimension(:,:), allocatable :: var1
        integer, dimension(:,:), allocatable :: mask1 
        logical :: coarse 
        double precision :: var_lim, var_step  

        nx = size(mask_pack,1)
        ny = size(mask_pack,2)

        ! === Step 1: Get the field dimensions correct

        ! Allocate temporary variable fields
        allocate(var1(nx,ny),mask1(nx,ny))

        if (size(mask_pack) .eq. size(var)) then
            ! Field is the same size as the mask
 
            var1 = var 

        else 
            ! Mapping is needed to adjust field dimensions to fit mask

            if (.not. (present(map))) then 
                write(*,*) "subset_gen_mask:: ", &
                "Error: map must be provided as an argument if the field "//&
                "provided does not match the mask."
                write(*,*) "var dims: ", size(var,1), size(var,2)
                write(*,*) "mask_pack dims: ", nx, ny 
                stop 
            end if 

            call map_field(map,"subset mask",var,var1,mask1,method="radius",fill=.TRUE.)

        end if 

        ! Step 2: Check the method (min or max)
        if ( .not. (trim(method) .eq. "min" .or. trim(method) .eq. "max") )  then 
            write(*,*) "subset_gen_mask:: ", &
            "Error: method can only be 'min' or 'max', not "//trim(method)
            stop
        end if 

        ! Take the negative of field if the method is min
        if (trim(method) .eq. "min") var1 = -var1 

        ! === Step 2: Decide which indices of the mask 
        ! === should be true based on user criteria.

        if ( nx*ny .eq. npts ) then 
            ! Total points is equal to the grid size,
            ! therefore all indices are counted.
            mask_pack = .TRUE. 

        else 
            ! Initially set all points to .FALSE.
            mask_pack = .FALSE.
            npts0     = 0

            ! Fill in mask_pack evenly-spaced at desired minimum resolution
            if (min_spacing .gt. 0) then

                do i = 1,nx,min_spacing
                do j = 1,ny,min_spacing
                    if (count(mask_pack) .ge. npts) exit 
                    mask_pack(i,j) = .TRUE. 
                end do 
                end do
            end if 

            ! Now fill in the remainder of the mask with the locations of field maximums
            npts0 = count(mask_pack)
            coarse = .TRUE.
            var_lim  = maxval(var1)
            var_step = var_lim*0.05 

            do while(npts0 < npts)
                where(mask_pack) var1 = minval(var1)   ! Eliminate points that are already included
                
                if (coarse) then 
                    var_lim = var_lim - var_step 
                    if (count(var1 .ge. var_lim)+count(mask_pack) .lt. npts) then 
                        where (var1 .ge. var_lim) mask_pack = .TRUE. 
                    else
                        coarse = .FALSE. 
                    end if 
                else 
                    loc = maxloc(var1)              ! Determine location of maximum value
                    mask_pack(loc(1),loc(2)) = .TRUE. ! Add this point to mask 
                end if 
                npts0 = count(mask_pack)          ! Recount masked points
                !write(*,*) "npts0 =", npts0 
            end do 

        end if 

        write(*,*) "mask_pack total = ",count(mask_pack)," / " , (nx*ny)

        return 

    end subroutine subset_gen_mask


end module subset2


