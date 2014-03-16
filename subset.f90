

module subset

    use coordinates
    use ncio 

    implicit none

    type subset_class

        type(points_class) :: pts
        type(grid_class)   :: grid 
        type(map_class)    :: map_tosub, map_fromsub 
        integer            :: npts 
        logical            :: subset
        logical, dimension(:,:), allocatable :: mask_pack  
    end type 

    private
    public :: subset_class 
    public :: subset_init, subset_redefine
    public :: subset_to_grid, subset_to_points

contains 

    subroutine subset_init(sub,grid,factor,npts,max_neighbors,lat_lim)
        ! Determine the coordinates of the domain

        implicit none 

        type(subset_class) :: sub    ! The subset definition
        type(grid_class)   :: grid   ! The original grid 
        integer :: npts              ! Number of points to allocate to the subset
        integer :: factor            ! Resolution factor (should be >= 1)
        integer, optional :: max_neighbors ! Maximum number of neighbors to use for mapping
        integer :: max_neighbs 
        double precision, optional :: lat_lim
        double precision :: lat_limit 

        character(len=12) :: suffix 
        double precision, allocatable, dimension(:) :: x, y

        ! Define a suffix to append to the grid name
        if (factor .ge. 10) then
            write(suffix,"(a3,i2)") "-hi", factor 
        else if (factor .ge. 1) then
            write(suffix,"(a3,i1)") "-hi", factor 
        else
            write(*,*) "subset_define:: Error: factor must be greater than or equal to one."
            stop 
        end if 

        max_neighbs = 10 
        if (present(max_neighbors)) max_neighbs = max_neighbors

        lat_limit   = 4.d0 
        if (present(lat_lim)) lat_limit = lat_lim 

        ! Initialize the new grid with input grid characteristics
        ! but with the new resolution
        call grid_init(sub%grid,name=trim(grid%name)//trim(suffix),mtype=grid%mtype, &
                       units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                       x0=grid%G%x(1),dx=grid%G%dx/dble(factor),nx=(grid%G%nx-1)*factor+1, &
                       y0=grid%G%y(1),dy=grid%G%dy/dble(factor),ny=(grid%G%ny-1)*factor+1, &
!                        x0=grid%G%x(1),dx=dx,nx=( grid%G%x(grid%G%nx)-grid%G%x(1) ) / dx, &
!                        y0=grid%G%y(1),dy=dy,ny=( grid%G%y(grid%G%ny)-grid%G%y(1) ) / dy, &
                       lambda=grid%proj%lambda,phi=grid%proj%phi,alpha=grid%proj%alpha, &
                       x_e=grid%proj%x_e,y_n=grid%proj%y_n)

        ! Allocate the packing mask for later use 
        if (allocated(sub%mask_pack)) deallocate(sub%mask_pack)
        allocate(sub%mask_pack(sub%grid%G%nx,sub%grid%G%ny))
        sub%mask_pack = .TRUE. 

        ! Make sure the subset npts is consistent with the new grid
        sub%npts = npts 
        if (sub%npts .gt. sub%grid%npts .or. sub%npts .le. 0) sub%npts = sub%grid%npts

        if (sub%npts .eq. sub%grid%npts) then 

            sub%subset = .FALSE.
            call grid_to_points(sub%grid,sub%pts)

        else

            sub%subset = .TRUE. 

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
            call points_init(sub%pts,trim(grid%name)//trim(suffix),mtype=grid%mtype, &
                             units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                             x=x,y=y, &
                             lambda=grid%proj%lambda,phi=grid%proj%phi, &
                             alpha=grid%proj%alpha,x_e=grid%proj%x_e,y_n=grid%proj%y_n)

            ! Initialize to and fro mappings for subset grid and input grid
            ! (intial map generation can take some time)
            call map_init(sub%map_tosub,grid,sub%grid, &
                          max_neighbors=max_neighbs,lat_lim=lat_limit,fldr="maps",load=.TRUE.)

            call map_init(sub%map_fromsub,sub%grid,grid, &
                          max_neighbors=max_neighbs,lat_lim=lat_limit,fldr="maps",load=.TRUE.)

        end if 
        
        write(*,"(a,i10,1x,a1,1x,i10)") "subset:: subset_init :: Initialized subset, npts = ", &
                                       sub%npts, "/",sub%grid%npts
        return

    end subroutine subset_init 

    subroutine subset_redefine(sub,mask_pack)
        ! Re-determine the coordinates of the current domain,
        ! which may be a subset of points from a grid
        ! This function allows the surface calculations to adapt to a changing 
        ! topography (eg, set increased density of points where topography is steep)

        implicit none 

        type(subset_class) :: sub
        logical, dimension(:,:) :: mask_pack 

        ! Check that mask_pack is consistent with subset of npts
        if (count(mask_pack) .ne. sub%npts) then
            write(*,"(a)") "subset:: subset_redefine:: Error: packing mask must specify the same "// &
                       "number of points as defined in the subset."
            write(*,*) "subset npts =",sub%npts
            write(*,*) "mask_pack total = ",count(mask_pack)
            stop 
        end if 

        write(*,*) "subset_redefine::"
        write(*,*) sub%grid%npts, sub%pts%npts 

        ! Get a new subset of points from the grid
        call grid_to_points(sub%grid,sub%pts,mask_pack=mask_pack,define=.FALSE.)

        return

    end subroutine subset_redefine 

    subroutine subset_to_grid(sub,var1D,var2D,mask_pack,map, &
                                    method,radius,fill,border,missing_value)
        ! This subroutine maps a subset of points (var1D) onto
        ! a 2D array (var2D) of resolution grid. 
        ! The subset should already be initialized.
        ! Note: currently var1D and var2D are double precision! 

        implicit none 
 
        type(subset_class), intent(IN)  :: sub 
        type(map_class), intent(IN), optional :: map 
        double precision, intent(OUT)   :: var2D(:,:)
        double precision, intent(IN)    :: var1D(:)
        logical, intent(IN)             :: mask_pack(:,:)

        double precision, allocatable :: var2Dtmp(:,:) 
        integer, allocatable          :: mask2D(:,:)

        character(len=*)           :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border

        type(map_class) :: map_local 

        ! Determine map to use here 
        map_local = sub%map_fromsub 
        if (present(map)) map_local = map 

        ! Assign a missing_value for use with mapping routine
        missing_val = -9999.d0
        if (present(missing_value)) missing_val = missing_value 

        write(*,*) trim(map_local%name1), " => ",trim(map_local%name2)

        ! Step 1: Unpack the 1D variable onto its corresponding 
        ! predefined 2D grid.
        call grid_allocate(sub%grid,var2Dtmp)        ! Allocate 2D array
        var2Dtmp = missing_val                       ! Prefill with missing_value
        var2Dtmp = unpack(var1D,mask_pack,var2Dtmp)  ! Unpack the 1D vector 
        
        ! Step 2: Map the temporary 2D array to the desired 2D resolution
        if (allocated(mask2D)) deallocate(mask2D)
        allocate(mask2D(map_local%G%nx,map_local%G%ny))

        call map_field(map_local,"Mapped variable",var2Dtmp,var2D,mask2D,method=method, &
                       radius=radius,fill=fill,border=border,missing_value=missing_val)

        return

    end subroutine subset_to_grid

    subroutine subset_to_points(sub,var2D,var1D,mask_pack,map, &
                                      method,radius,fill,border,missing_value)
        ! This subroutine maps a 2D array (var2D) onto
        ! a subset of points (var1D) of resolution sub%grid. 
        ! The subset should already be initialized.
        ! Note: currently var1D and var2D are double precision! 

        implicit none 
 
        type(subset_class), intent(IN)  :: sub 
        type(map_class), intent(IN), optional :: map 
        double precision, intent(IN)    :: var2D(:,:)
        double precision, intent(OUT)   :: var1D(:)
        logical, intent(IN)             :: mask_pack(:,:)

        double precision, allocatable :: var2Dtmp(:,:) 
        integer, allocatable          :: mask2D(:,:)

        character(len=*)           :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val
        logical, optional :: fill, border 

        type(map_class) :: map_local 

        ! Determine map to use here 
        map_local = sub%map_tosub 
        if (present(map)) map_local = map 

        ! Assign a missing_value for use with mapping routine
        missing_val = -9999.d0
        if (present(missing_value)) missing_val = missing_value 

        ! Step 1: Allocate a predefined 2D array.
        call grid_allocate(sub%grid,var2Dtmp)        ! Allocate 2D array
        var2Dtmp = missing_val                       ! Prefill with missing_value

        ! Step 2: Map the 2D array to the temporary 2D array of the subset
        if (allocated(mask2D)) deallocate(mask2D)
        allocate(mask2D(map_local%G%nx,map_local%G%ny))
        call map_field(map_local,"Mapped variable",var2D,var2Dtmp,mask2D,method=method, &
                       radius=radius,fill=fill,border=border,missing_value=missing_val)

        ! Step 3: Pack the 2D variable onto its corresponding predefined 1D points.
        var1D = pack(var2Dtmp,mask_pack)
        
        return

    end subroutine subset_to_points


end module subset


