

module subset

    use coordinates
    use ncio 

    implicit none

    type subset_class

        type(points_class) :: pts
        type(grid_class)   :: grid 
        type(map_class)    :: map_tosub, map_fromsub 
        integer            :: npts 
    end type 

    private
    public :: subset_class 
    public :: subset_init, subset_points_redefine
    public :: subset_points_to_grid, subset_grid_to_points

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
            write(suffix,"(a3,i2)") "_hi", factor 
        else if (factor .ge. 1) then
            write(suffix,"(a3,i1)") "_hi", factor 
        else
            write(*,*) "subset_define:: Error: factor must be greater than one."
            stop 
        end if 

        max_neighbs = 10 
        if (present(max_neighbors)) max_neighbs = max_neighbors

        lat_limit   = 4.d0 
        if (present(lat_lim)) lat_limit = lat_lim 

!         ! First intialize a grid definition if it is needed
!         write(*,*) "Grid file: "//trim(coord_info(1))

!         ! Get x and y axis values of desired grid
!         nx = nc_size(coord_info(1),coord_info(2))
!         ny = nc_size(coord_info(1),coord_info(3))
!         if (allocated(x_axis)) deallocate(x_axis)
!         if (allocated(y_axis)) deallocate(y_axis)
!         allocate(x_axis(nx),y_axis(ny))
!         call nc_read(coord_info(1),coord_info(2),x_axis)
!         call nc_read(coord_info(1),coord_info(3),y_axis)

!         ! Initialize the grid with the correct axis values from the file
!         call grid_init(grid,filename=coord_def,x=x_axis,y=y_axis)
        
        ! Assume that the original grid is already defined. The subset
        ! will be generated from a higher resolution version of the same
        ! grid 

        ! First, initialize the new grid with input grid characteristics
        ! but with the new resolution
        call grid_init(sub%grid,name=trim(grid%name)//trim(suffix),mtype=grid%mtype, &
                       units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                       x0=grid%G%x(1),dx=grid%G%dx/factor,nx=(grid%G%nx-1)*factor+1, &
                       y0=grid%G%y(1),dy=grid%G%dy/factor,ny=(grid%G%ny-1)*factor+1, &
                       lambda=grid%proj%lambda,phi=grid%proj%phi,alpha=grid%proj%alpha, &
                       x_e=grid%proj%x_e,y_n=grid%proj%y_n)

        ! Make sure the subset npts is consistent with the new grid
        sub%npts = npts 
        if (sub%npts .gt. sub%grid%npts .or. sub%npts .lt. 0) sub%npts = sub%grid%npts

        ! Allocate temporary x,y vectors to store points of interest
        ! for new coordinates definition
        if (allocated(x)) deallocate(x)
        if (allocated(y)) deallocate(y)
        allocate(x(npts),y(npts))

        ! For now, assign dummy x/y values with correct length.
        ! The actual values will be determined later from 
        ! packing/unpacking at each time step
        x = 0.d0
        x(1:2) = [grid%G%x(1),grid%G%x(grid%G%nx)]
        y = 0.d0 
        y(1:2) = [grid%G%y(1),grid%G%y(grid%G%ny)]

        call points_init(sub%pts,trim(grid%name)//trim(suffix),mtype=grid%mtype, &
                         units=grid%units,planet=grid%planet%name,lon180=grid%is_lon180, &
                         x=x,y=y, &
                         lambda=grid%proj%lambda,phi=grid%proj%phi, &
                         alpha=grid%proj%alpha,x_e=grid%proj%x_e,y_n=grid%proj%y_n)

        call map_init(sub%map_tosub,grid,sub%grid, &
                      max_neighbors=max_neighbs,lat_lim=lat_limit,fldr="maps",load=.TRUE.)

        call map_init(sub%map_fromsub,sub%grid,grid, &
                      max_neighbors=max_neighbs,lat_lim=lat_limit,fldr="maps",load=.TRUE.)

        return

    end subroutine subset_init 

    subroutine subset_points_redefine(pts,grid,mask_pack,is_grid)
        ! Re-determine the coordinates of the current domain,
        ! which may be a subset of points from a grid
        ! This function allows the surface calculations to adapt to a changing 
        ! topography (eg, set increased density of points where topography is steep)

        implicit none 

        type(points_class) :: pts
        type(grid_class)   :: grid
        logical            :: is_grid 
        logical, dimension(:,:) :: mask_pack 

        ! Get a new subset of points from the grid
        call grid_to_points(grid,pts,mask_pack)

        return

    end subroutine subset_points_redefine 

    subroutine subset_points_to_grid(grid,var1D,var2D,mask_pack,map, &
                                      method,radius,missing_value)
        ! This function will return a 2D array (double precision!) of values
        ! obtained from a 1D array that was either packed directly from the 
        ! 2D array, or the 1D array should be interpolated onto the 2D grid
        ! using map_field
        ! Output array var2D should have size of grid, or grid_new if present.

        implicit none 
 
        type(grid_class)    :: grid
        type(map_class), optional :: map   
        double precision :: var1D(:), var2D(:,:)
        logical          :: mask_pack(:,:)

        double precision, allocatable :: var2Dtmp(:,:) 
        integer, allocatable          :: mask2D(:,:)

        character(len=*), optional :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val

        ! Assign a missing_value for use with mapping routine
        missing_val = -9999.d0
        if (present(missing_value)) missing_val = missing_value 

        ! Step 1: Unpack the 1D variable onto its corresponding 
        ! predefined 2D grid.
        call grid_allocate(grid,var2Dtmp)      ! Allocate 2D array
        var2Dtmp = missing_val                 ! Prefill with missing_value
        var2Dtmp = unpack(var1D,mask_pack,var2Dtmp)  ! Unpack the 1D vector 
        
        if (present(map)) then
            ! Map the temporary 2D array to the desired 2D resolution

            if (allocated(mask2D)) deallocate(mask2D)
            allocate(mask2D(map%G%nx,map%G%ny))
            call map_field(map,"Mapped variable",var2Dtmp,var2D,mask2D, &
                           method=method,radius=radius,missing_value=missing_val)

        else
            ! Directly fill the output array with the unpacked values
            
            where (var2Dtmp .ne. missing_value) var2D = var2Dtmp

        end if 
        
        return

    end subroutine subset_points_to_grid

    subroutine subset_grid_to_points(pts,var2D,var1D,mask_pack,map, &
                                      method,radius,missing_value)
        ! This function will return a 1D array (double precision!) of values
        ! obtained from a 2D array that was either packed directly from the 
        ! 2D array, or the 2D array should be interpolated onto the 1D grid
        ! using map_field
        ! Output array var2D should have size of grid, or grid_new if present.

        implicit none 
 
        type(points_class)  :: pts
        type(map_class), optional :: map   
        double precision :: var2D(:,:), var1D(:)
        logical          :: mask_pack(:,:)

        double precision, allocatable :: var1Dtmp(:) 
        integer, allocatable          :: mask1D(:)

        character(len=*), optional :: method
        double precision, optional :: radius, missing_value 
        double precision :: missing_val

        ! Assign a missing_value for use with mapping routine
        missing_val = -9999.d0
        if (present(missing_value)) missing_val = missing_value 

        ! Step 1: Pack the 2D variable onto its corresponding 
        ! predefined 1D points.
        call points_allocate(pts,var1Dtmp)    ! Allocate 1D array
        var1Dtmp = pack(var2D,mask_pack)      ! Pack the 2D vector 
        
        if (present(map)) then
            ! Map the temporary 1D array to the desired 1D resolution

            if (allocated(mask1D)) deallocate(mask1D)
            allocate(mask1D(map%npts))

            call map_field(map,"Mapped variable",var1Dtmp,var1D,mask1D, &
                           method=method,radius=radius,missing_value=missing_val)

        else
            ! Directly fill the output array with the unpacked values
            
            where (var1Dtmp .ne. missing_value) var1D = var1Dtmp

        end if 
        
        return

    end subroutine subset_grid_to_points


end module subset


