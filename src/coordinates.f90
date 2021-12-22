!> ####################################################################
!! **Module**     : coordinates \n
!! **Author**     : Alex Robinson \n 
!! **Purpose**    : This module defines the coordinate classes to handle
!!                  grids and points in a given reference space. This module can
!!                  be used to interpolate and project between different
!!                  grids and points. Coordinates are represented by lon, lat at 
!!                  each point, and cartesian x,y values for stereographic grids.
!!                  Mapping classes store interpolation information for
!!                  'fast mapping' between different grids and sets of points, 
!!                  using projection routines from the oblimap2 package.
!! ####################################################################
module coordinates

    use coord_constants 

    use oblimap_projection_module
    use planet 
    use ncio

!     use gaussian_filter  
!     use interp2D 
    
!     use mod_toms526 

    implicit none 

    type points_class 

        character (len=128) :: name     ! name of this domain (world, GRL, etc)
        character (len=128) :: mtype    ! map type: latlon, cartesian, stereographic, etc (following cf conventions)
        character (len=128) :: units    ! units of the axes
        type(planet_class)  :: planet   ! which planet are we on?

        ! Projection parameters
        logical :: is_cartesian, is_projection
        logical :: is_lon180
        type(projection_class) :: proj

        ! Points information
        integer :: npts
        real(dp), allocatable, dimension(:)   :: x, y, lon, lat, area
        real(dp), allocatable, dimension(:)   :: dx, dy 
        integer,  allocatable, dimension(:)   :: border
        real(dp) :: xy_conv 

    end type 

    type grid_axis_class
        integer :: nx, ny
        real(dp) :: x0, dx, y0, dy
        real(dp), allocatable, dimension(:) :: x, y    
    end type

    type grid_class 

        character (len=128) :: name 
        character (len=128) :: mtype     ! latlon, cartesian, stereographic, etc
        character (len=128) :: units
        type(planet_class)  :: planet 

        ! Projection parameters
        logical :: is_cartesian, is_projection
        logical :: is_lon180
        type(projection_class) :: proj

        ! Points information
        integer :: npts
        real(dp), allocatable, dimension(:,:) :: x, y, lon, lat, area
        integer,  allocatable, dimension(:,:) :: border
        real(dp) :: xy_conv

        ! Grid axes information 
        type(grid_axis_class) :: G
        
    end type 

    interface grid_init
        module procedure grid_init_from_opts, grid_init_from_par
        module procedure grid_init_from_grid, grid_init_from_points
    end interface 
    
    interface grid_allocate 
        module procedure grid_allocate_integer, grid_allocate_float
        module procedure grid_allocate_double,  grid_allocate_logical
    end interface

    interface points_init
        module procedure points_init_from_opts, points_init_from_par
        module procedure points_init_from_points, points_init_from_grid
        module procedure points_init_from_file 
    end interface 

    interface points_allocate 
        module procedure points_allocate_integer, points_allocate_float
        module procedure points_allocate_double,  points_allocate_logical
    end interface

    interface compare_coord 
        module procedure compare_coord_grid_grid, compare_coord_points_points
    end interface 

    private
    public :: points_class, points_init, points_allocate, points_write, points_print 
    public :: grid_axis_class, grid_class, grid_init, grid_allocate, grid_write, grid_print
    public :: grid_to_points, points_to_grid
    public :: compare_coord
    public :: pts_which_nearest 

contains

    subroutine axis_init(x,x0,dx)

        implicit none 

        real(dp) :: x(:)
        real(dp), optional :: x0, dx
        real(dp) :: dx_tmp 
        integer :: i, nx  

        nx = size(x) 

        do i = 1, nx 
            x(i) = dble(i-1)
        end do 

        dx_tmp = 1.d0 
        if (present(dx)) dx_tmp = dx 
        
        x = x*dx_tmp  

        if (present(x0)) then 
            x = x + x0 
        else
            x = x + (-(nx-1.d0)/2.d0*dx_tmp)
        end if 

        return 
    end subroutine axis_init

    subroutine grid_area(grid)
        ! Calculate the area of the cell each grid point represents

        implicit none 

        type(grid_class) :: grid 
        real(dp) :: x1, x2, y1, y2 
        real(dp) :: vertx(4), verty(4), vertlon(4), vertlat(4)
        integer  :: i, j, q 
        
        do j = 1, grid%G%ny 
            do i = 1, grid%G%nx 

                ! Get vertices of cell around current grid point
                if (i .eq. 1) then 
                    x1 =  grid%G%x(i) - (grid%G%x(i+1)-grid%G%x(i)) / 2.0_dp
                    x2 = (grid%G%x(i) + grid%G%x(i+1)) / 2.0_dp
                else if (i .eq. grid%G%nx) then 
                    x1 = (grid%G%x(i) + grid%G%x(i-1)) / 2.0_dp
                    x2 =  grid%G%x(i) + (grid%G%x(i)-grid%G%x(i-1)) / 2.0_dp
                else
                    x1 = (grid%G%x(i) + grid%G%x(i-1)) / 2.0_dp
                    x2 = (grid%G%x(i) + grid%G%x(i+1)) / 2.0_dp
                end if 

                if (j .eq. 1) then 
                    y1 =  grid%G%y(j) - (grid%G%y(j+1)-grid%G%y(j)) / 2.0_dp
                    y2 = (grid%G%y(j) + grid%G%y(j+1)) / 2.0_dp
                else if (j .eq. grid%G%ny) then 
                    y1 = (grid%G%y(j) + grid%G%y(j-1)) / 2.0_dp
                    y2 =  grid%G%y(j) + (grid%G%y(j)-grid%G%y(j-1)) / 2.0_dp
                else
                    y1 = (grid%G%y(j) + grid%G%y(j-1)) / 2.0_dp
                    y2 = (grid%G%y(j) + grid%G%y(j+1)) / 2.0_dp
                end if 

                ! Convert grid points to meters as necessary
                vertx = [x1,x1,x2,x2] * grid%xy_conv 
                verty = [y1,y2,y2,y1] * grid%xy_conv 
                
                if (.not. grid%is_projection .and. grid%is_cartesian) then 
                    ! Calculate cartesian area
                    grid%area(i,j) = cartesian_area(vertx,verty)

                else if (grid%is_projection) then 
                    ! Calculate projected area 

                    do q = 1, 4
!                         call inverse_oblique_sg_projection(vertx(q),verty(q), &
!                                                            vertlon(q),vertlat(q),grid%proj)
                        call oblimap_projection_inverse(vertx(q),verty(q), &
                                                        vertlon(q),vertlat(q),grid%proj)
                    end do 

                    ! Calculate latlon area
                    grid%area(i,j) = planet_area(grid%planet%a,grid%planet%f,vertlon,vertlat)

                else 
                    ! Calculate latlon area 
                    grid%area(i,j) = planet_area(grid%planet%a,grid%planet%f,vertx,verty)

                end if 

            end do 
        end do 

        ! Return grid area in same units as axes
        if (.not. grid%xy_conv .eq. 0.d0) grid%area = grid%area / (grid%xy_conv*grid%xy_conv) 

        return

    end subroutine grid_area

    subroutine grid_init_from_grid(grid,grid0,name,x,y,x0,dx,nx,y0,dy,ny)
        ! Initialize a new grid based on an old grid but
        ! with new x/y coordinates 

        implicit none

        type(grid_class)   :: grid, grid0
        character(len=*)   :: name  
        integer, optional  :: nx, ny
        real(dp), optional :: x(:), y(:), x0, dx, y0, dy
        
        call grid_init_from_opts(grid,name,mtype=grid0%mtype,units=grid0%units, &
                                 planet=grid0%planet%name,lon180=grid0%is_lon180, &
                                 x=x,y=y,x0=x0,dx=dx,nx=nx,y0=y0,dy=dy,ny=ny, &
                                 lambda=grid0%proj%lambda,phi=grid0%proj%phi, &
                                 alpha=grid0%proj%alpha,x_e=grid0%proj%x_e,y_n=grid0%proj%y_n)

        return

    end subroutine grid_init_from_grid

    subroutine grid_init_from_points(grid,pts0,name,x,y,x0,dx,nx,y0,dy,ny)
        ! Initialize a new grid based on an old grid but
        ! with new x/y coordinates 

        implicit none

        type(grid_class)   :: grid
        type(points_class) :: pts0
        character(len=*)   :: name  
        integer, optional  :: nx, ny
        real(dp), optional :: x(:), y(:), x0, dx, y0, dy
        
        call grid_init_from_opts(grid,name,mtype=pts0%mtype,units=pts0%units, &
                                 planet=pts0%planet%name,lon180=pts0%is_lon180, &
                                 x=x,y=y,x0=x0,dx=dx,nx=nx,y0=y0,dy=dy,ny=ny, &
                                 lambda=pts0%proj%lambda,phi=pts0%proj%phi, &
                                 alpha=pts0%proj%alpha,x_e=pts0%proj%x_e,y_n=pts0%proj%y_n)

        return

    end subroutine grid_init_from_points

    subroutine grid_init_from_par(grid,filename,x,y,x0,dx,nx,y0,dy,ny)

        implicit none

        type(grid_class)   :: grid 
        character(len=*)   :: filename 
        integer, optional  :: nx, ny
        real(dp), optional :: x(:), y(:), x0, dx, y0, dy

        character(len=256)   :: name, mtype, units
        character(len=256)   :: planet  
        logical              :: lon180
        real(dp)             :: lambda, phi, alpha, x_e, y_n 

        namelist /map/ name, mtype, units, planet, lon180, &
                       lambda, phi, alpha, x_e, y_n  

        open(7,file=trim(filename))
        read(7,nml=map)
        close(7)

        call grid_init_from_opts(grid,name,mtype,units,planet,lon180, &
                       x,y,x0,dx,nx,y0,dy,ny,lambda,phi,alpha,x_e,y_n)

        return

    end subroutine grid_init_from_par

    subroutine grid_init_from_opts(grid,name,mtype,units,planet,lon180, &
                                   x,y,x0,dx,nx,y0,dy,ny,lambda,phi,alpha,x_e,y_n)

        implicit none 

        type(grid_class)   :: grid 
        type(points_class) :: pts
        character(len=*)   :: name, mtype, units
        character(len=*), optional :: planet  
        logical, optional  :: lon180
        integer, optional  :: nx, ny 
        real(dp), optional :: x(:), y(:), x0, dx, y0, dy 
        real(dp), optional :: lambda, phi, alpha, x_e, y_n 
        real(dp) :: tmp1, tmp2
        integer :: i, j 
        real(dp), allocatable :: tmpvec1(:), tmpvec2(:) 

        ! Initially deallocate grid arrays just in case 
        if (allocated(grid%G%x))     deallocate(grid%G%x)
        if (allocated(grid%G%y))     deallocate(grid%G%y)
        if (allocated(grid%x))       deallocate(grid%x)
        if (allocated(grid%y))       deallocate(grid%y)
        if (allocated(grid%lon))     deallocate(grid%lon)
        if (allocated(grid%lat))     deallocate(grid%lat)

        ! Set up the grid axis info
        ! Note: axis values can represent cartesian or latlon dimensions
        ! depending on map type (mtype)
        
        ! x-axis
        if (present(x)) then 
            grid%G%nx = size(x) 
            allocate(grid%G%x(grid%G%nx))
            grid%G%x  = x 
        else 
            grid%G%nx = nx 
            allocate(grid%G%x(grid%G%nx))
            call axis_init(grid%G%x,x0=x0,dx=dx)
        end if 

        ! y-axis
        if (present(y)) then 
            grid%G%ny = size(y) 
            allocate(grid%G%y(grid%G%ny))
            grid%G%y  = y
        else 
            grid%G%ny = ny
            allocate(grid%G%y(grid%G%ny))
            call axis_init(grid%G%y,x0=y0,dx=dy)
        end if

        ! How many points does grid contain?
        grid%npts = grid%G%nx * grid%G%ny 

        ! Determine dx, dy if not provided 
        if (present(dx)) then
            grid%G%dx = dx 
        else
            grid%G%dx = sum(grid%G%x(2:grid%G%nx)-grid%G%x(1:grid%G%nx-1))/(grid%G%nx-1)
        end if 

        if (present(dy)) then
            grid%G%dy = dy 
        else
            grid%G%dy = sum(grid%G%y(2:grid%G%ny)-grid%G%y(1:grid%G%ny-1))/(grid%G%ny-1)
        end if 
        
        ! Allocate and generate 2D point sets (x,y)
        ! Note x,y represent cartesian values or latlon 
        ! depending on mtype
        ! lon,lat allocated as duplicates in latlon case, for consistency
        allocate(grid%x(grid%G%nx,grid%G%ny))
        allocate(grid%y(grid%G%nx,grid%G%ny))
        allocate(grid%lon(grid%G%nx,grid%G%ny))
        allocate(grid%lat(grid%G%nx,grid%G%ny))

        ! Store axis values in 2D arrays too
        do i = 1, grid%G%nx 
            grid%y(i,:) = grid%G%y 
        end do 

        do j = 1, grid%G%ny 
            grid%x(:,j) = grid%G%x 
        end do 

        ! Initialize data as points, then convert back to grid using axis info
        call points_init_from_opts(pts,name,mtype,units,planet,lon180,reshape(grid%x,[grid%npts]),&
                         reshape(grid%y,[grid%npts]),lambda=lambda,phi=phi,alpha=alpha,x_e=x_e,y_n=y_n)
        call points_to_grid(pts,grid)

        ! Calculate grid cell areas 
        if (allocated(grid%area))   deallocate(grid%area)
        allocate(grid%area(grid%G%nx,grid%G%ny))
        call grid_area(grid) 
        
        ! Assign border points
        if (allocated(grid%border)) deallocate(grid%border)
        allocate(grid%border(grid%G%nx,grid%G%ny))
        grid%border              = 0
        grid%border(1,:)         = 1 
        grid%border(:,1)         = 1 
        grid%border(grid%G%nx,:) = 1 
        grid%border(:,grid%G%ny) = 1 

        ! Adjust border points if global latlon grid is being used 
        select case(trim(grid%mtype))

            case("latitude_longitude","latlon","gaussian")
                ! latitude_longitude grid 

                ! First check longitude distances
                tmp1 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(1),grid%G%x(2),grid%G%y(1))
                tmp2 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(1),grid%G%x(grid%G%nx),grid%G%y(1))
                
                ! If distance between first and last longitude is less than 
                ! twice the distance of the first and second longitude,
                ! then grid wraps around the globe longitudinally
                if (tmp2 .le. tmp1*2.0_dp) then 
                    grid%border(1,:)         = 0
                    grid%border(grid%G%nx,:) = 0 
                end if 

                ! Next check latitude distances

                ! First check lower latitude
                tmp1 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(1),grid%G%x(1),grid%G%y(2))
                tmp2 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(1),grid%G%x(1),-90.0_dp)
                if (tmp2 .le. tmp1*2.0_dp) grid%border(:,1)         = 0

                ! Next check upper latitude
                tmp1 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(grid%G%ny),grid%G%x(1),grid%G%y(grid%G%ny-1))
                tmp2 = planet_distance(grid%planet%a,grid%planet%f,grid%G%x(1),grid%G%y(grid%G%ny),grid%G%x(1),90.0_dp)
                if (tmp2 .le. tmp1*2.0_dp) grid%border(:,grid%G%ny) = 0 

            case DEFAULT 
                ! Projection and/or cartesian grid 

                ! Pass - no adjustment needed 
        end select 

        ! Make sure final lon/lat values are in desired range
        ! (default range is 0=>360)
        if ( .not. grid%is_cartesian ) then 
            if ( grid%is_lon180 ) &  ! Then make range -180=>180
                where( grid%G%x .gt. 180.0_dp ) grid%G%x = grid%G%x - 360.0_dp 

                ! Also for plotting programs, we need to shift matrix, so that
                ! axis goes from -180 to 180

        end if 

        ! Finally check that the proj%alpha used is reasonable 
        if ( grid%is_projection ) then 
            tmp1 = optimal_alpha(grid%planet%R,grid%G%nx,grid%G%ny, &
                                  grid%G%dx*grid%xy_conv,grid%G%dy*grid%xy_conv)

            write(*,"(a,f6.1)") "Optimal alpha (degrees) = ", tmp1 
        end if 

        ! Output a summary of grid axis information
        call grid_print(grid)

        return 

    end subroutine grid_init_from_opts

    subroutine points_init_from_points(pts,pts0,name,filename,x,y,dx,dy,latlon,skip)

        implicit none

        type(points_class) :: pts, pts0 
        character(len=*)   :: name 
        character(len=*), optional :: filename 
        real(dp), optional :: x(:), y(:), dx(:), dy(:)
        logical, optional  :: latlon
        integer, optional  :: skip 

        if (present(filename)) then 
            call points_init_from_file(pts,name,mtype=pts0%mtype,units=pts0%units, &
                                       filename=trim(filename),planet=pts0%planet%name, &
                                       lon180=pts0%is_lon180, &
                                       lambda=pts0%proj%lambda,phi=pts0%proj%phi, &
                                       alpha=pts0%proj%alpha,x_e=pts0%proj%x_e,y_n=pts0%proj%y_n, &
                                       latlon=latlon,skip=skip)

        else
            call points_init_from_opts(pts,name,mtype=pts0%mtype,units=pts0%units, &
                                     planet=pts0%planet%name,lon180=pts0%is_lon180,x=x,y=y,dx=dx,dy=dy, &
                                     lambda=pts0%proj%lambda,phi=pts0%proj%phi, &
                                     alpha=pts0%proj%alpha,x_e=pts0%proj%x_e,y_n=pts0%proj%y_n, &
                                     latlon=latlon)
        end if 

        return

    end subroutine points_init_from_points

    subroutine points_init_from_grid(pts,grid0,name,filename,x,y,dx,dy,latlon,skip)

        implicit none

        type(points_class) :: pts
        type(grid_class)   :: grid0 
        character(len=*)   :: name  
        character(len=*), optional :: filename
        real(dp), optional :: x(:), y(:), dx(:), dy(:)
        logical, optional  :: latlon 
        integer, optional  :: skip 

        if (present(filename)) then 
            call points_init_from_file(pts,name,mtype=grid0%mtype,units=grid0%units, &
                                       filename=trim(filename),planet=grid0%planet%name, &
                                       dx=grid0%G%dx,dy=grid0%G%dy,lon180=grid0%is_lon180, &
                                       lambda=grid0%proj%lambda,phi=grid0%proj%phi, &
                                       alpha=grid0%proj%alpha,x_e=grid0%proj%x_e,y_n=grid0%proj%y_n, &
                                       latlon=latlon,skip=skip)

        else
            call points_init_from_opts(pts,name,mtype=grid0%mtype,units=grid0%units, &
                                       planet=grid0%planet%name,lon180=grid0%is_lon180,x=x,y=y,dx=dx,dy=dy, &
                                       lambda=grid0%proj%lambda,phi=grid0%proj%phi, &
                                       alpha=grid0%proj%alpha,x_e=grid0%proj%x_e,y_n=grid0%proj%y_n, &
                                       latlon=latlon)
        end if

        return

    end subroutine points_init_from_grid

    subroutine points_init_from_par(pts,filename,x,y,dx,dy,latlon)

        implicit none

        type(points_class) :: pts 
        character(len=*)   :: filename 
        real(dp), optional :: x(:), y(:)
        real(dp), optional :: dx(:), dy(:) 
        logical, optional  :: latlon 

        character(len=256)   :: name, mtype, units
        character(len=256)   :: planet  
        logical              :: lon180
        real(dp)             :: lambda, phi, alpha, x_e, y_n 

        namelist /map/ name, mtype, units, planet, lon180, &
                       lambda, phi, alpha, x_e, y_n  

        open(7,file=trim(filename))
        read(7,nml=map)
        close(7)

        call points_init_from_opts(pts,name,mtype,units,planet,lon180, &
                         x,y,dx,dy,lambda,phi,alpha,x_e,y_n,latlon)

        return

    end subroutine points_init_from_par

    subroutine points_init_from_file(pts,name,mtype,units,filename,planet,dx,dy,lon180, &
                                     lambda,phi,alpha,x_e,y_n,latlon,skip)


        implicit none

        type(points_class) :: pts 
        character(len=*) :: name, mtype, units
        character(len=*) :: filename  
        character(len=*), optional :: planet
        real(dp), optional :: dx, dy 
        logical,  optional :: lon180
        real(dp), optional :: lambda, phi, alpha, x_e, y_n 
        logical,  optional :: latlon 
        integer,  optional :: skip 

        integer :: i, io, n, nskip   
        integer, parameter :: nmax = 1000000
        character(len=10)  :: tmp_char 
        real(dp) :: x_in(nmax), y_in(nmax)
        real(dp), allocatable :: x(:), y(:)
        real(dp), allocatable :: dx_vec(:), dy_vec(:)

        nskip = 0 
        if (present(skip)) nskip = skip 

        open(7,file=trim(filename))
        if (nskip .gt. 0) then 
            do i = 1, nskip 
                read(7,*) tmp_char 
            end do 
        end if 
        do i = 1, nmax
            read(7,*,iostat=io) x_in(i), y_in(i)
            if (io .lt. 0) exit 
        end do 
        close(7)

        n = i-1
        allocate(x(n),y(n))

        x = x_in(1:n)
        y = y_in(1:n)

        write(*,*) "points_init_from_file: ", minval(x), maxval(x)
        write(*,*) "points_init_from_file: ", minval(y), maxval(y)
        
        ! Set dx, dy according to input option 
        allocate(dx_vec(n),dy_vec(n))
        dx_vec = 0.0_dp
        dy_vec = 0.0_dp 
        if (present(dx)) dx_vec = dx 
        if (present(dy)) dy_vec = dy 
         
        call points_init_from_opts(pts,name,mtype,units,planet,lon180, &
                         x,y,dx_vec,dy_vec,lambda,phi,alpha,x_e,y_n,latlon)

        write(*,*) "points_init_from_file: ", minval(pts%lon), maxval(pts%lon)
        write(*,*) "points_init_from_file: ", minval(pts%lat), maxval(pts%lat)
        
        return

    end subroutine points_init_from_file

    subroutine points_init_from_opts(pts,name,mtype,units,planet,lon180,x,y,dx,dy, &
                                     lambda,phi,alpha,x_e,y_n,latlon,verbose)

        use oblimap_projection_module 

        implicit none 

        type(points_class) :: pts 
        real(dp)           :: x(:), y(:)
        real(dp), optional :: dx(:), dy(:) 
        character(len=*) :: name, mtype, units 
        character(len=*), optional :: planet
        character(len=256) :: planet_name
        logical, optional  :: lon180
        real(dp), optional :: lambda, phi, alpha, x_e, y_n 
        logical, optional :: latlon 
        logical, optional :: verbose 
        
        ! Local variables 
        logical :: latlon_in 
        logical :: verbose_in 
        integer :: i, nborder 
        integer, allocatable :: tmpi(:)

        ! Assign basic dataset info 
        pts%name  = trim(name)
        pts%mtype = trim(mtype)
        pts%units = trim(units)
        
        verbose_in = .TRUE. 
        if (present(verbose)) verbose_in = verbose 

        ! Also define the standard cf name
        select case(trim(pts%mtype))
            case("latitude_longitude","latlon","gaussian")
                pts%is_cartesian  = .FALSE. 
                pts%is_projection = .FALSE. 
            case("stereographic","polar_stereographic","lambert_azimuthal_equal_area")
                pts%is_cartesian  = .TRUE. 
                pts%is_projection = .TRUE. 
            case("cartesian")
                pts%is_cartesian  = .TRUE. 
                pts%is_projection = .FALSE. 
            case DEFAULT 
                write(*,"(a7,a20,a)")  &
                    "coord::","points_init: ","error: map type not allowed:"//trim(pts%mtype)
                write(*,*) "    map type must be one of the following: "
                write(*,*) "        latitude_longitude"
                write(*,*) "        latlon"
                write(*,*) "        gaussian"
                write(*,*) "        cartesian"
                write(*,*) "        stereographic"
                write(*,*) "        polar_stereographic"
                write(*,*) "        lambert_azimuthal_equal_area"
                write(*,*) 
                stop 
        end select 

        ! Make sure we can convert the units of the points as needed
        select case(trim(pts%units))
            case("kilometers","km")
                pts%xy_conv = 1.d3  ! To convert from km => m (for cartesian grids)
            case DEFAULT
                pts%xy_conv = 1.d0 
        end select

        ! Check latlon range (0=>360 or -180=>180)
        pts%is_lon180 = .FALSE.
        if (present(lon180)) pts%is_lon180 = lon180 

        ! Make sure x and y vectors have the same length
        if (size(x) .ne. size(y)) then
            write(*,"(a7,a20,a)") "coord::","points_init: ", &
                       "error: x and y points must have the same length."
            stop
        end if

        ! Check whether input points are xy values or latlon values (for projected grid)
        latlon_in = .FALSE.
        if (present(latlon)) latlon_in = latlon

        if (latlon_in .and. (pts%is_cartesian .and. .not. pts%is_projection)) then 
            write(*,*) "points_init:: error: x/y input values can only &
                       &be latlon values for projected grids or for latlon grids."
            stop 
        end if 

        ! Assign point information
        pts%npts = size(x)

        ! Reallocate point vectors
        if (allocated(pts%x))      deallocate(pts%x)
        if (allocated(pts%y))      deallocate(pts%y)
        if (allocated(pts%lon))    deallocate(pts%lon)
        if (allocated(pts%lat))    deallocate(pts%lat)
        if (allocated(pts%border)) deallocate(pts%border)
        if (allocated(pts%area))   deallocate(pts%area)
        allocate(pts%x(pts%npts),pts%y(pts%npts))
        allocate(pts%lon(pts%npts),pts%lat(pts%npts))
        allocate(pts%border(pts%npts))
        allocate(pts%area(pts%npts))

        if (allocated(pts%dx))      deallocate(pts%dx)
        if (allocated(pts%dy))      deallocate(pts%dy)
        allocate(pts%dx(pts%npts),pts%dy(pts%npts))

        ! Careful here: if the argument x = pts%x, the reallocation
        ! causes a memory gap sometimes. Is there a check for this?
        ! (ajr, 2013-09-28)
        pts%x = x 
        pts%y = y 

        ! Set dx/dy to zero for now, without additional information 
        ! (if it is a grid, info will be provided in grid_init)
        pts%dx = 0.0_dp 
        pts%dy = 0.0_dp 
        if (present(dx)) pts%dx = dx 
        if (present(dy)) pts%dy = dy 

        ! For now set points cell areas to 1
        ! (if it is a grid, area will be calculated afterwards in grid_init)
        pts%area   = 1.0_dp

        ! If it is just a set of points, none are considered border points
        pts%border = 0

        ! Initialize planetary reference information
        planet_name = "WGS84"
        if (present(planet)) planet_name = trim(planet)
        call planet_init(pts%planet,planet_name)

        ! If cartesian grid is a projection, calculate the corresponding latlon points
        ! Note: xy points not needed for a latlon grid except during mapping
        if ( pts%is_projection ) then
            ! Points are projected, find coordinates (latlon or xy)
            
            ! Initialize projection information
            call projection_init(pts%proj,trim(pts%mtype),pts%planet, &
                                 lambda,phi,alpha,x_e,y_n)

            if (latlon_in) then 
                ! Starting with latlon points 

                pts%lon = pts%x 
                pts%lat = pts%y 

                do i = 1, pts%npts       
                    call oblimap_projection(pts%lon(i),pts%lat(i),pts%x(i),pts%y(i),pts%proj)
                    pts%x(i) = pts%x(i)/pts%xy_conv
                    pts%y(i) = pts%y(i)/pts%xy_conv
                end do

            else 
                ! Starting with xy points 

                do i = 1, pts%npts       
                    call oblimap_projection_inverse(pts%x(i)*pts%xy_conv,pts%y(i)*pts%xy_conv, &
                                                    pts%lon(i),pts%lat(i),pts%proj)
                end do

            end if 

        else if ( .not. pts%is_cartesian ) then 
            ! Points are defined in latlon space

            pts%lon = pts%x 
            pts%lat = pts%y 

            ! Initialize projection information with planet information
            ! (to ensure parameters a and f are defined)
            call projection_init(pts%proj,"Undefined",pts%planet)

        end if 

        ! Make sure final lon/lat values are in desired range
        ! (default range is 0=>360)
        if ( pts%is_projection .or. .not. pts%is_cartesian ) then 
            if ( pts%is_lon180 ) &  ! Then make range -180=>180
                where( pts%lon .gt. 180.0_dp ) pts%lon = pts%lon - 360.0_dp 
        end if 

        ! Print a summary of the set of points if desired
        if (verbose_in) call points_print(pts)

        return 

    end subroutine points_init_from_opts

    subroutine grid_to_points(grid,pts,mask_pack,define)
        ! Converts a grid class to points class
        ! (ie, 2D grid => 1D vector of points)

        implicit none 

        type(grid_class)   :: grid 
        type(points_class) :: pts 
        logical, optional  :: mask_pack(:,:)
        integer, allocatable :: maski(:,:)
        logical, optional  :: define 
        logical            :: define_fields 

        ! Are we defining a new set of points?
        define_fields = .TRUE.
        if (present(define)) define_fields = define 

        ! Assign points constants
        if (define_fields) then 
            pts%name          = trim(grid%name) 
            pts%mtype         = trim(grid%mtype)
            pts%units         = trim(grid%units)
            pts%npts          = grid%npts
            pts%is_cartesian  = grid%is_cartesian 
            pts%is_projection = grid%is_projection 
            pts%is_lon180     = grid%is_lon180
            pts%planet        = grid%planet
            pts%proj          = grid%proj
            pts%xy_conv       = grid%xy_conv 
        end if 

        ! Make a mask for packing (since mask_pack is optional!)
        if (allocated(maski)) deallocate(maski)
        allocate(maski(grid%G%nx,grid%G%ny))
        maski = 1 
        if (present(mask_pack)) then
            where(.not. mask_pack) maski = 0
        end if 

        ! Make sure mask is consistent with desired pts output
        if (sum(maski) .ne. pts%npts) then
            write(*,*) "grid_to_points:: Error, "// &
                       "total masked values not equal to npts."
            write(*,*) "count(mask_pack) npts:",sum(maski), pts%npts 
            write(*,*) "Grid name:   "//trim(grid%name)
            write(*,*) "Points name: "//trim(pts%name)
            stop 
        end if 

        ! Deallocate all points fields
        if (define_fields) then 
            if (allocated(pts%x))      deallocate(pts%x)
            if (allocated(pts%y))      deallocate(pts%y)
            if (allocated(pts%area))   deallocate(pts%area)
            if (allocated(pts%dx))     deallocate(pts%dx)
            if (allocated(pts%dy))     deallocate(pts%dy)
            if (allocated(pts%border)) deallocate(pts%border)
            if (allocated(pts%lon))    deallocate(pts%lon)
            if (allocated(pts%lat))    deallocate(pts%lat)
        end if 

        ! Store x,y points
        if (define_fields) allocate(pts%x(pts%npts),pts%y(pts%npts))
        pts%x = pack(grid%x,maski==1)
        pts%y = pack(grid%y,maski==1)

        ! Store area 
        if (define_fields) allocate(pts%area(pts%npts))
        pts%area = pack(grid%area,maski==1)

        ! Store border 
        if (define_fields) allocate(pts%border(pts%npts))
        pts%border = pack(grid%border,maski==1)

        ! Store dx, dy 
        if (define_fields) allocate(pts%dx(pts%npts),pts%dy(pts%npts))
        pts%dx = grid%G%dx 
        pts%dy = grid%G%dy 

        ! If lon,lat points exist on grid, store them too
        if (allocated(grid%lon) .and. allocated(grid%lat)) then 
            if (define_fields) allocate(pts%lon(pts%npts),pts%lat(pts%npts))
            pts%lon = pack(grid%lon,maski==1)
            pts%lat = pack(grid%lat,maski==1)
        end if 

        return

    end subroutine grid_to_points

    subroutine points_to_grid(pts,grid,mask_pack,define)
        ! Converts pts class to a grid class
        ! (ie 1D vector => 2D grid points),
        ! assuming that grid axis info (grid%G) has 
        ! already been defined.

        implicit none 

        type(grid_class)   :: grid 
        type(points_class) :: pts 
        logical, optional  :: mask_pack(:,:)
        integer, allocatable :: maski(:,:)
        logical, optional  :: define 
        logical            :: define_fields 

        ! Are we defining a new set of points?
        define_fields = .TRUE.
        if (present(define)) define_fields = define 

        ! Assign grid constants
        grid%name          = trim(pts%name) 
        grid%mtype         = trim(pts%mtype)
        grid%units         = trim(pts%units)
        grid%npts          = pts%npts
        grid%is_cartesian  = pts%is_cartesian 
        grid%is_projection = pts%is_projection 
        grid%is_lon180     = pts%is_lon180
        grid%planet        = pts%planet 
        grid%proj          = pts%proj
        grid%xy_conv       = pts%xy_conv 

        ! Make a mask for unpacking
        if (allocated(maski)) deallocate(maski)
        allocate(maski(grid%G%nx,grid%G%ny))
        maski = 1 
        if (present(mask_pack)) then
            where(.not. mask_pack) maski = 0
        end if 

        ! Make sure mask is consistent with desired pts output
        if (sum(maski) .ne. pts%npts) then
            write(*,*) "points_to_grid:: Error, "// &
                       "total masked values not equal to npts, sum(mask) npts:", &
                       sum(maski), pts%npts 
            write(*,*) "Grid name:   "//trim(grid%name)
            write(*,*) "Points name: "//trim(pts%name)
            stop 
        end if 

        ! Reallocate all grid fields 
        if (define_fields .and. allocated(grid%x))      deallocate(grid%x)
        if (define_fields .and. allocated(grid%y))      deallocate(grid%y)
        if (define_fields .and. allocated(grid%area))   deallocate(grid%area)
        if (define_fields .and. allocated(grid%border)) deallocate(grid%border)
        if (define_fields .and. allocated(grid%lon))    deallocate(grid%lon)
        if (define_fields .and. allocated(grid%lat))    deallocate(grid%lat)
        
        if (define_fields) then
            allocate(grid%x(grid%G%nx,grid%G%ny),grid%y(grid%G%nx,grid%G%ny))
            grid%x = MISSING_VALUE_DEFAULT
            grid%y = MISSING_VALUE_DEFAULT
        end if 

        grid%x = unpack(pts%x,maski==1,grid%x)
        grid%y = unpack(pts%y,maski==1,grid%y)
        
        ! Store area 
        if (define_fields) then
            allocate(grid%area(grid%G%nx,grid%G%ny))
            grid%area = MISSING_VALUE_DEFAULT 
        end if 
        grid%area = unpack(pts%area,maski==1,grid%area)

        ! Store border 
        if (define_fields) then 
            allocate(grid%border(grid%G%nx,grid%G%ny))
            grid%border = MISSING_VALUE_DEFAULT
        end if 
        grid%border = unpack(pts%border,maski==1,grid%border)

        if (allocated(pts%lon) .and. allocated(pts%lat)) then 
            if (define_fields) then
                allocate(grid%lon(grid%G%nx,grid%G%ny),grid%lat(grid%G%nx,grid%G%ny))
                grid%lon = MISSING_VALUE_DEFAULT 
                grid%lat = MISSING_VALUE_DEFAULT 
            end if 
            grid%lon = unpack(pts%lon,maski==1,grid%lon)
            grid%lat = unpack(pts%lat,maski==1,grid%lat)
        end if 

        return

    end subroutine points_to_grid

    subroutine grid_allocate_integer(grid,var,value)

        implicit none 

        type(grid_class) :: grid 
        integer, allocatable :: var(:,:)
        integer, optional :: value 

        if (allocated(var)) deallocate(var)
        allocate(var(grid%G%nx,grid%G%ny))
        if (present(value)) var = value 

        return 

    end subroutine grid_allocate_integer

    subroutine grid_allocate_double(grid,var,value)

        implicit none 

        type(grid_class) :: grid 
        real(dp), allocatable :: var(:,:)
        real(dp), optional :: value 

        if (allocated(var)) deallocate(var)
        allocate(var(grid%G%nx,grid%G%ny))
        if (present(value)) var = value 
        
        return 

    end subroutine grid_allocate_double

    subroutine grid_allocate_float(grid,var,value)

        implicit none 

        type(grid_class) :: grid 
        real(4), allocatable :: var(:,:)
        real(4), optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(grid%G%nx,grid%G%ny))
        if (present(value)) var = value 
        
        return 

    end subroutine grid_allocate_float

    subroutine grid_allocate_logical(grid,var,value)

        implicit none 

        type(grid_class) :: grid 
        logical, allocatable :: var(:,:)
        logical, optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(grid%G%nx,grid%G%ny))
        if (present(value)) var = value 
        
        return 

    end subroutine grid_allocate_logical

    subroutine points_allocate_integer(points,var,value)

        implicit none 

        type(points_class) :: points 
        integer, allocatable :: var(:)
        integer, optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(points%npts))
        if (present(value)) var = value 
        
        return 

    end subroutine points_allocate_integer

    subroutine points_allocate_double(points,var,value)

        implicit none 

        type(points_class) :: points 
        real(dp), allocatable :: var(:)
        real(dp), optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(points%npts))
        if (present(value)) var = value 
        
        return 

    end subroutine points_allocate_double

    subroutine points_allocate_float(points,var,value)

        implicit none 

        type(points_class) :: points 
        real(4), allocatable :: var(:)
        real(4), optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(points%npts))
        if (present(value)) var = value 
        
        return 

    end subroutine points_allocate_float

    subroutine points_allocate_logical(points,var,value)

        implicit none 

        type(points_class) :: points 
        logical, allocatable :: var(:)
        logical, optional :: value 
        
        if (allocated(var)) deallocate(var)
        allocate(var(points%npts))
        if (present(value)) var = value 
        
        return 

    end subroutine points_allocate_logical

    function compare_coord_points_points(pts1,pts2) result(same_map)

        implicit none 

        type(points_class) :: pts1, pts2 
        logical :: same_map 

        ! Check if both maps use the same projection  
        if (pts1%is_projection .and. pts2%is_projection        .and. &
            trim(pts1%mtype)       .eq. trim(pts2%mtype)       .and. &
            same_projection(pts1%proj,pts2%proj) ) then 
            
            ! Both maps come from the same projection
            same_map = .TRUE. 

        else if (pts1%is_cartesian .and. .not. pts1%is_projection .and. &
                 pts2%is_cartesian .and. .not. pts2%is_projection) then 
            ! Both maps are generic cartesian grids (assume origin is the same)
            same_map = .TRUE. 

        else if (.not. pts1%is_cartesian .and. .not. pts2%is_cartesian) then 
            ! Both maps are latlon maps
            same_map = .TRUE. 

        else
            same_map = .FALSE.

        end if 

        return 

    end function compare_coord_points_points

    function compare_coord_grid_grid(grid1,grid2) result(same_map)

        implicit none 

        type(grid_class) :: grid1, grid2
        logical :: same_map 

        ! Check if both maps use the same projection  
        if (grid1%is_projection .and. grid2%is_projection        .and. &
            trim(grid1%mtype)       .eq. trim(grid2%mtype)       .and. &
            same_projection(grid1%proj,grid2%proj) ) then 
            
            ! Both maps come from the same projection
            same_map = .TRUE. 

        else if (grid1%is_cartesian .and. .not. grid1%is_projection .and. &
                 grid2%is_cartesian .and. .not. grid2%is_projection) then 
            ! Both maps are generic cartesian grids (assume origin is the same)
            same_map = .TRUE. 

        else if (.not. grid1%is_cartesian .and. .not. grid2%is_cartesian) then 
            ! Both maps are latlon maps
            same_map = .TRUE. 

        else
            same_map = .FALSE.

        end if 

        return 

    end function compare_coord_grid_grid

    subroutine grid_write(grid,fnm,xnm,ynm,create)

        implicit none 
        type(grid_class) :: grid 
        character(len=*) :: fnm,xnm,ynm
        logical :: create  

        ! Local variables 
        character(len=24) :: xunits, yunits 

        if (trim(grid%units) .eq. "degrees") then 
            xunits = "degrees_east"
            yunits = "degrees_north"
        else 
            xunits = trim(grid%units)
            yunits = trim(grid%units)
        end if 

        ! Create the netcdf file if desired
        if (create) then 
            call nc_create(fnm)
        
            ! Add grid axis variables to netcdf file
            call nc_write_dim(fnm,xnm,x=grid%G%x,units=xunits)
            call nc_write_dim(fnm,ynm,x=grid%G%y,units=yunits)
            
            if (grid%is_projection) then 
                call nc_write_attr(fnm,xnm,"standard_name","projection_x_coordinate")
                call nc_write_attr(fnm,ynm,"standard_name","projection_y_coordinate")
                ! call nc_write_attr(fnm,xnm,"_CoordinateAxisType","GeoX")
                ! call nc_write_attr(fnm,ynm,"_CoordinateAxisType","GeoY")
            end if 

        end if 

        if (grid%is_projection) then 
            ! Add coordinate reference system information 
            call nc_write_map(fnm,grid%mtype,grid%proj%lambda,phi=grid%proj%phi, &
                              alpha=grid%proj%alpha,x_e=grid%proj%x_e,y_n=grid%proj%y_n, &
                              is_sphere=grid%planet%is_sphere,semi_major_axis=grid%planet%a, &
                              inverse_flattening=1.d0/grid%planet%f)
        end if 

        if (grid%is_projection .or. grid%is_cartesian) then 
            call nc_write(fnm,"x2D",grid%x,dim1=xnm,dim2=ynm,grid_mapping="crs")
            call nc_write_attr(fnm,"x2D","units",grid%units)
            call nc_write(fnm,"y2D",grid%y,dim1=xnm,dim2=ynm,grid_mapping="crs")
            call nc_write_attr(fnm,"y2D","units",grid%units)
            
        end if 

        if ( grid%is_projection ) then 
            ! Write 2D arrays of lon/lat values that correspond to current projection
            ! (else there are no lon/lat values, or they are the 1D axes)
            call nc_write(fnm,"lon2D",grid%lon,dim1=xnm,dim2=ynm,grid_mapping="crs")
            !call nc_write_attr(fnm,"lon2D","_CoordinateAxisType","Lon")
            call nc_write_attr(fnm,"lon2D","units","degrees_east")
                
            call nc_write(fnm,"lat2D",grid%lat,dim1=xnm,dim2=ynm,grid_mapping="crs")
            !call nc_write_attr(fnm,"lat2D","_CoordinateAxisType","Lat")
            call nc_write_attr(fnm,"lat2D","units","degrees_north")
            
        end if 

        if ( grid%is_projection ) then
            call nc_write(fnm,"area",  grid%area,  dim1=xnm,dim2=ynm,grid_mapping="crs")
            call nc_write(fnm,"border",grid%border,dim1=xnm,dim2=ynm,grid_mapping="crs")
        else 
            call nc_write(fnm,"area",  grid%area,  dim1=xnm,dim2=ynm)
            call nc_write(fnm,"border",grid%border,dim1=xnm,dim2=ynm)
        end if

        ! Add coordinates information as needed
        if ( grid%is_projection ) then 
            call nc_write_attr(fnm,"x2D",   "coordinates","lat2D lon2D")
            call nc_write_attr(fnm,"y2D",   "coordinates","lat2D lon2D")
            call nc_write_attr(fnm,"area",  "coordinates","lat2D lon2D")
            call nc_write_attr(fnm,"border","coordinates","lat2D lon2D")
        end if 

        return
        
    end subroutine grid_write

    subroutine points_write(pts,fnm,xnm,ynm,create)

        implicit none 
        type(points_class) :: pts 
        character(len=*) :: fnm,xnm,ynm
        logical :: create  

        ! Local variables 
        character(len=24) :: xunits, yunits 

        if (trim(pts%units) .eq. "degrees") then 
            xunits = "degrees_east"
            yunits = "degrees_north"
        else 
            xunits = trim(pts%units)
            yunits = trim(pts%units)
        end if 

        ! Create the netcdf file if desired
        if (create) then 
            call nc_create(fnm)
        
            ! Add axis variable to netcdf file
            call nc_write_dim(fnm,"point",x=1,nx=pts%npts,units="n")
        end if 

        ! Add projection information if needed
        if (pts%is_projection) then 
            call nc_write_map(fnm,pts%mtype,pts%proj%lambda,phi=pts%proj%phi,&
                            alpha=pts%proj%alpha,x_e=pts%proj%x_e,y_n=pts%proj%y_n, &
                            is_sphere=pts%planet%is_sphere,semi_major_axis=pts%planet%a, &
                            inverse_flattening=1.d0/pts%planet%f)
        end if 

        if (pts%is_projection .or. pts%is_cartesian) then 
            call nc_write(fnm,xnm,pts%x,dim1="point",grid_mapping=pts%mtype)
            call nc_write_attr(fnm,xnm,"units",xunits)
            call nc_write(fnm,ynm,pts%y,dim1="point",grid_mapping=pts%mtype)
            call nc_write_attr(fnm,ynm,"units",yunits)
        end if 

        if (pts%is_projection) then 
            call nc_write_attr(fnm,xnm,"standard_name","projection_x_coordinate")
            call nc_write_attr(fnm,ynm,"standard_name","projection_y_coordinate")
            ! call nc_write_attr(fnm,xnm,"_CoordinateAxisType","GeoX")
            ! call nc_write_attr(fnm,ynm,"_CoordinateAxisType","GeoY")
        end if 

        if ( pts%is_projection .or. (.not. pts%is_cartesian) ) then  
            call nc_write(fnm,"lon",pts%lon,dim1="point",grid_mapping=pts%mtype)
            !call nc_write_attr(fnm,"lon","_CoordinateAxisType","Lon")

            call nc_write(fnm,"lat",pts%lat,dim1="point",grid_mapping=pts%mtype)
            !call nc_write_attr(fnm,"lat","_CoordinateAxisType","Lat")

        end if 

        ! if ( pts%is_projection ) then
        !     call nc_write(fnm,"area",  pts%area,  dim1=xnm,dim2=ynm,pts_mapping="crs")
        !     call nc_write(fnm,"border",pts%border,dim1=xnm,dim2=ynm,pts_mapping="crs")
        ! else 
        !     call nc_write(fnm,"area",  pts%area,  dim1=xnm,dim2=ynm)
        !     call nc_write(fnm,"border",pts%border,dim1=xnm,dim2=ynm)
        ! end if

        return

    end subroutine points_write

    subroutine grid_print(grid)

        implicit none 

        type(grid_class), intent(IN) :: grid
        type(points_class) :: pts 

        write(*,*) "== Grid summary   ====="
        write(*,"(a16,i6,i6)")     "nx,ny = ",grid%G%nx, grid%G%ny 
        write(*,"(a16,2g12.5)")    "dx,dy = ",grid%G%dx, grid%G%dy
        write(*,"(a16,2g12.5)") "range(x-axis) = ",minval(grid%G%x),maxval(grid%G%x)
        write(*,"(a16,2g12.5)") "range(y-axis) = ",minval(grid%G%y),maxval(grid%G%y)
        write(*,*) 

        return

    end subroutine grid_print

    subroutine points_print(pts)

        implicit none 

        type(points_class), intent(IN) :: pts 

        write(*,*) "== Points summary ========================="
        write(*,"(a24,a)") "Data set: ",trim(pts%name)
        write(*,"(a24,a)") "Map type: ",trim(pts%mtype)
        write(*,"(a24,a)") "Coordinate units: ",trim(pts%units)
        write(*,"(a16,i8)")           "Total points = ",pts%npts
        if (pts%is_cartesian) then 
            write(*,"(a16,2g12.5)")    "range(x) = ",minval(pts%x),maxval(pts%x)
            write(*,"(a16,2g12.5)")    "range(y) = ",minval(pts%y),maxval(pts%y)
        end if 
        if (pts%is_projection .or. .not. pts%is_cartesian) then 
            write(*,"(a16,2g12.5)")    "range(lon) = ",minval(pts%lon),maxval(pts%lon)
            write(*,"(a16,2g12.5)")    "range(lat) = ",minval(pts%lat),maxval(pts%lat)
        end if 
        
        if (pts%is_projection) then 
!             write(*,*) "Projection information"
            write(*,"(a16,a)")         "proj method = ", trim(pts%proj%method)
            write(*,"(a16,g12.5)")     "          a = ", pts%proj%a
            write(*,"(a16,g12.5)")     "          f = ", pts%proj%f
            write(*,"(a16,g12.5)")     "     lambda = ", pts%proj%lambda
            write(*,"(a16,g12.5)")     "        phi = ", pts%proj%phi
            write(*,"(a16,g12.5)")     "      alpha = ", pts%proj%alpha
            write(*,"(a16,g12.5)")     "        x_e = ", pts%proj%x_e
            write(*,"(a16,g12.5)")     "        y_n = ", pts%proj%y_n
        end if 
        
!         write(*,*)

        return

    end subroutine points_print
    
    function pts_which_nearest(pts,x,y,latlon) result(idx)

        implicit none 

        type(points_class), intent(IN) :: pts 
        double precision,   intent(IN) :: x, y 
        logical, optional  :: latlon 
        integer            :: idx 

        ! Local variables
        logical               :: is_latlon
        real(dp), allocatable :: dist(:)
        integer :: i 

        ! Determine if input location is in units of latlon or xy
        is_latlon = .FALSE.
        if (present(latlon)) is_latlon = latlon 

        if (is_latlon .and. (pts%is_cartesian .and. .not. pts%is_projection)) then 
            write(*,*) "pts_which_nearest:: error: "// &
                       "x and y cannot by latlon values unless coordinates have latlon defined."
            stop
        end if 

        allocate(dist(pts%npts))

        if (is_latlon) then 
            ! Calculate distances using planet distances (in meters)
            do i = 1, pts%npts
                dist(i) = planet_distance(pts%planet%a,pts%planet%f,x,y,pts%lon(i),pts%lat(i))
            end do 

        else
            ! Calculate distances using cartesian distances (in meters)
            do i = 1, pts%npts
                dist(i) = cartesian_distance(x,y,pts%x(i)*pts%xy_conv,pts%y(i)*pts%xy_conv)
            end do 

        end if 

        ! Find the index of the nearest point
        idx = minloc(dist,dim=1)

        return 

    end function pts_which_nearest

    subroutine progress(j,ntot)
        ! Progress bar, thanks to:
        ! http://thelazycatholic.wordpress.com/2010/08/19/progress-bar-in-fortran/

        implicit none
        integer :: j,k,ntot
        character(len=18) :: bar="\r???% |          |"
        real(dp), parameter :: perc_out = 5.d0 
        real(dp) :: perc
        perc = dble(j)/dble(ntot)*100.d0

        if (mod(perc,perc_out).le.0.1d-3) then 
            ! updates the fraction of calculation done
            write(unit=bar(2:4),fmt="(i3)") 10*j
            do k = 1, j
            bar(7+k:7+k)="*"
            enddo

            ! print the progress bar.
            write(*,'(a)',advance='no') bar

            if (perc==100) write(*,*)
        end if 

        return

    end subroutine progress

end module coordinates

