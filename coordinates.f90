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

    use oblimap_projection_module
    use planet 
    use ncio 

    implicit none 

    !! real(dp) definition and some internal constants
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    type points_class 

        character (len=128) :: name     ! name of this domain (world, GRL, etc)
        character (len=128) :: mtype    ! map type: latlon, cartesian, stereographic, etc
        character (len=128) :: units    ! units of the axes
        type(planet_class)  :: planet   ! which planet are we on?

        ! Projection parameters
        logical :: is_cartesian, is_projection
        logical :: is_lon180
        type(projection_class) :: proj

        ! Points information
        integer :: npts
        real(dp), allocatable, dimension(:)   :: x, y, lon, lat, area
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

        type(grid_axis_class) :: G
        integer :: npts
        real(dp), allocatable, dimension(:,:) :: x, y, lon, lat, area
        integer,  allocatable, dimension(:,:) :: border
        real(dp) :: xy_conv

    end type 

    type map_helper_class 
        real(dp), dimension(:), allocatable   :: var, weight
        real(dp), dimension(:), allocatable   :: var_tmp 
        integer,  dimension(:), allocatable   :: mask2
        logical,  dimension(:), allocatable   :: maskp 
    end type 

    type map_class
        character (len=128) :: name1, name2     ! Names of coordinate set1 and set2
        character (len=128) :: mtype
        character (len=128) :: units 
        type(planet_class)  :: planet 

        ! Projection parameters
        logical  :: is_cartesian, is_projection, is_same_map
        logical  :: is_lon180
        type(projection_class) :: proj

        ! Grid related variables
        logical :: is_grid
        type(grid_axis_class) :: G    ! Used if is_grid==.TRUE. 
 
        integer :: nbnd               ! Number of boundary points 

        ! Vectors of x, y, lon and lat points (will be allocated to size npts)
        integer :: npts
        real(dp), allocatable, dimension(:) :: x, y, lon, lat
        real(dp) :: xy_conv

        ! Neighbor arrays (will be allocated to size npts,nmax)
        integer :: nmax 
        integer,  dimension(:,:), allocatable :: i, quadrant, border
        real(sp), dimension(:,:), allocatable :: dist, weight

        ! Helper arrays for map_field
        ! (preallocating here saves time during mapping)
!         type(map_helper_class) :: tmp   

    end type

    interface grid_init
        module procedure grid_init_from_opts, grid_init_from_par
        module procedure grid_init_from_grid
    end interface 
    
    interface grid_allocate 
        module procedure grid_allocate_integer, grid_allocate_float
        module procedure grid_allocate_double,  grid_allocate_logical
    end interface

    interface points_init
        module procedure points_init_from_opts, points_init_from_par
    end interface 

    interface points_allocate 
        module procedure points_allocate_integer, points_allocate_float
        module procedure points_allocate_double,  points_allocate_logical
    end interface

    interface map_init 
        module procedure map_init_points_points, map_init_grid_grid 
        module procedure map_init_grid_points,   map_init_points_grid
    end interface

    interface map_field 
        module procedure map_field_grid_grid_double,   map_field_points_points_double
        module procedure map_field_grid_points_double, map_field_points_grid_double
        module procedure map_field_grid_grid_integer,  map_field_points_points_integer
        module procedure map_field_grid_points_integer, map_field_points_grid_integer
    end interface

    interface compare_map 
        module procedure compare_map_grid_grid, compare_map_points_points 
        module procedure compare_map_map_grid 
    end interface

    private
    public :: points_class, points_init, points_allocate, points_write, points_print 
    public :: grid_class, grid_init, grid_allocate, grid_write, grid_print
    public :: grid_to_points, points_to_grid
    public :: compare_map
    public :: map_class, map_init, map_field, map_print

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
        real(dp) :: vertx(4), verty(4) 
        integer  :: i, j 
        
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
                vertx = (/x1,x1,x2,x2/) * grid%xy_conv 
                verty = (/y1,y2,y2,y1/) * grid%xy_conv 
                
                if (grid%is_cartesian) then 
                    grid%area(i,j) = cartesian_area(vertx,verty)
                else
                    grid%area(i,j) = planet_area(grid%planet%a,grid%planet%f,vertx,verty)
                end if 

            end do 
        end do 

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
        call points_init_from_opts(pts,name,mtype,units,planet,lon180,reshape(grid%x,(/grid%npts/)),&
                         reshape(grid%y,(/grid%npts/)),lambda,phi,alpha,x_e,y_n)
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
        if ( trim(grid%mtype) .eq. "latlon" ) then 

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

        end if 

        ! Make sure final lon/lat values are in desired range
        ! (default range is 0=>360)
        if ( .not. grid%is_cartesian ) then 
            if ( grid%is_lon180 ) &  ! Then make range -180=>180
                where( grid%G%x .gt. 180.0_dp ) grid%G%x = grid%G%x - 360.0_dp 

                ! Also for plotting programs, we need to shift matrix, so that
                ! axis goes from -180 to 180

        end if 

        ! Output a summary of grid axis information
        call grid_print(grid)

        return 

    end subroutine grid_init_from_opts

    subroutine points_init_from_points(pts,pts0,name,x,y)

        implicit none

        type(points_class) :: pts, pts0 
        character(len=*)   :: name  
        real(dp), optional :: x(:), y(:)

        call points_init_from_opts(pts,name,mtype=pts0%mtype,units=pts0%units, &
                                 planet=pts0%planet%name,lon180=pts0%is_lon180,x=x,y=y, &
                                 lambda=pts0%proj%lambda,phi=pts0%proj%phi, &
                                 alpha=pts0%proj%alpha,x_e=pts0%proj%x_e,y_n=pts0%proj%y_n)

        return

    end subroutine points_init_from_points

    subroutine points_init_from_par(pts,filename,x,y)

        implicit none

        type(points_class) :: pts 
        character(len=*)   :: filename 
        real(dp), optional :: x(:), y(:)

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
                         x,y,lambda,phi,alpha,x_e,y_n)

        return

    end subroutine points_init_from_par

    subroutine points_init_from_opts(pts,name,mtype,units,planet,lon180,x,y, &
                                     lambda,phi,alpha,x_e,y_n)

        use oblimap_projection_module 

        implicit none 

        type(points_class) :: pts 
        real(dp)           :: x(:), y(:)
        character(len=*) :: name, mtype, units 
        character(len=*), optional :: planet
        character(len=256) :: planet_name
        logical, optional  :: lon180
        real(dp), optional :: lambda, phi, alpha, x_e, y_n 
        integer :: i, nborder 
        integer, allocatable :: tmpi(:)

        ! Assign basic dataset info 
        pts%name  = trim(name)
        pts%mtype = trim(mtype)
        pts%units = trim(units)
        
        ! Define map characteristics
        select case(trim(pts%mtype))
            case("latlon")
                pts%is_cartesian  = .FALSE. 
                pts%is_projection = .FALSE. 
            case("stereographic","polar stereographic")
                pts%is_cartesian  = .TRUE. 
                pts%is_projection = .TRUE. 
            case("cartesian")
                pts%is_cartesian  = .TRUE. 
                pts%is_projection = .FALSE. 
            case DEFAULT 
                write(*,"(a7,a20,a)")  &
                    "coord::","points_init: ","error: map type not allowed:"//trim(pts%mtype)
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
        allocate(pts%border(pts%npts))
        allocate(pts%area(pts%npts))

        ! Careful here: if the argument x = pts%x, the reallocation
        ! causes a memory gap sometimes. Is there a check for this?
        ! (ajr, 2013-09-28)
        pts%x = x 
        pts%y = y 

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
            ! Points are projected from latlon space, find latlon coordinates
            allocate(pts%lon(pts%npts),pts%lat(pts%npts))

            ! Initialize projection information
            call projection_init(pts%proj,"stereographic",pts%planet, &
                                 lambda,phi,alpha,x_e,y_n)

            do i = 1, pts%npts       
                call inverse_oblique_sg_projection(pts%x(i)*pts%xy_conv,pts%y(i)*pts%xy_conv, &
                                                   pts%lon(i),pts%lat(i),pts%proj)
            end do

        else if ( .not. pts%is_cartesian ) then 
            ! Points are defined in latlon space
            allocate(pts%lon(pts%npts),pts%lat(pts%npts))

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

        ! Print a summary of the set of points 
        call points_print(pts)

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

    subroutine map_init_grid_grid(map,grid1,grid2,max_neighbors,lat_lim,fldr,load,save)
        ! Generate mapping weights from grid1 to grid2

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid1, grid2 
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        ! Initialize map grid axis info (since mapping to grid2)
        map%is_grid = .TRUE. 
!         map%G%nx    = grid2%G%nx 
!         map%G%ny    = grid2%G%ny 
!         if (allocated(map%G%x)) deallocate(map%G%x)
!         if (allocated(map%G%y)) deallocate(map%G%y)
!         allocate(map%G%x(map%G%nx),map%G%y(map%G%ny))
!         map%G%x = grid2%G%x 
!         map%G%y = grid2%G%y 
        map%G = grid2%G 

        call grid_to_points(grid1,pts1)
        call grid_to_points(grid2,pts2)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)

        return 

    end subroutine map_init_grid_grid

    subroutine map_init_grid_points(map,grid1,pts2,max_neighbors,lat_lim,fldr,load,save)
        ! Generate mapping weights from grid to set of points

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid1 
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        map%is_grid = .FALSE. 
        call grid_to_points(grid1,pts1)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)

        return 

    end subroutine map_init_grid_points

    subroutine map_init_points_grid(map,pts1,grid2,max_neighbors,lat_lim,fldr,load,save)
        ! Generate mapping weights from set of points to grid

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        ! Initialize map grid axis info (since mapping to grid2)
        map%is_grid = .TRUE. 
!         map%G%nx    = grid2%G%nx 
!         map%G%ny    = grid2%G%ny 
!         if (allocated(map%G%x)) deallocate(map%G%x)
!         if (allocated(map%G%y)) deallocate(map%G%y)
!         allocate(map%G%x(map%G%nx),map%G%y(map%G%ny))
!         map%G%x = grid2%G%x 
!         map%G%y = grid2%G%y 
        map%G = grid2%G 

        ! Convert grid2 to points for map initialization
        call grid_to_points(grid2,pts2)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)

        return 

    end subroutine map_init_points_grid 

    subroutine map_init_points_points(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)
        ! Generate mapping weights from set of points to another set of points

        implicit none 

        type(points_class) :: pts1, pts2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        map%is_grid = .FALSE. 
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)

        return 

    end subroutine map_init_points_points

    subroutine map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,fldr,load,save)
        ! Generate mapping weights between sets of pointss

        implicit none 

        type(points_class) :: pts1, pts2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated
        logical :: load_file, save_file, fldr_exists, file_exists 
        character(len=256) :: mapfldr 

        ! Load file if it exists by default
        load_file = .TRUE. 
        if (present(load)) load_file = load 

        ! Save generated map to file by default
        save_file = .TRUE. 
        if (present(save)) save_file = save 

        mapfldr = "maps"
        if (present(fldr)) mapfldr = trim(fldr)

        ! Assign map constant information
        map%name1         = trim(pts1%name) 
        map%name2         = trim(pts2%name)
        map%mtype         = trim(pts2%mtype)
        map%units         = trim(pts2%units)
        map%is_projection = pts2%is_projection 
        map%is_cartesian  = pts2%is_cartesian
        map%is_lon180     = pts2%is_lon180
        map%planet        = pts2%planet
        map%proj          = pts2%proj 
        map%npts          = pts2%npts
        map%nmax          = max_neighbors 
        map%xy_conv       = pts2%xy_conv 

        ! Check if the same map is defined for both sets of points
        map%is_same_map = compare_map(pts1,pts2)

        ! Note: do not assign max distance here, save all distances
        ! up until the maximum number of neighbors
        ! Later, when loading map, let use choose max_distance
        ! In this way, less recalculation of maps will be needed
        ! when the max_distance changes.

        ! Determine if file matching these characteristics exists
        inquire(file=map_filename(map,mapfldr),exist=file_exists)

        !! Now load map information from file if exists and is desired
        !! or else calculate weights and store in file. 
        if ( load_file .and. file_exists ) then 

            ! Read map weights and info from file
            call map_read(map,mapfldr)

        else

            ! Reallocate and assign map point arrays
            if(allocated(map%x)) deallocate(map%x)
            if(allocated(map%y)) deallocate(map%y)
            if(allocated(map%lon)) deallocate(map%lon)
            if(allocated(map%lat)) deallocate(map%lat)
            allocate(map%x(map%npts),map%y(map%npts))
            allocate(map%lon(map%npts),map%lat(map%npts))
            
            map%x   = pts2%x
            map%y   = pts2%y 
            map%lon = pts2%lon
            map%lat = pts2%lat

            ! Reallocate map neighborhood arrays
            if (allocated(map%i)) deallocate(map%i)
            if (allocated(map%dist)) deallocate(map%dist)
            if (allocated(map%weight)) deallocate(map%weight)
            if (allocated(map%quadrant)) deallocate(map%quadrant)
            if (allocated(map%border)) deallocate(map%border)
            allocate(map%i(map%npts,map%nmax))
            allocate(map%dist(map%npts,map%nmax))
            allocate(map%weight(map%npts,map%nmax))
            allocate(map%quadrant(map%npts,map%nmax))
            allocate(map%border(map%npts,map%nmax))

            ! Calculate map weights (time consuming!)
            call map_calc_weights(map,pts1,pts2,lat_lim,2.0_dp)

            ! Write new map to file
            if (save_file) call map_write(map,mapfldr) 

        end if 

        ! Print map summary
        call map_print(map)

        return

    end subroutine map_init_internal

    subroutine map_calc_weights(map,pts1,pts2,lat_lim,shepard_exponent)
        implicit none 

        type(map_class) :: map 
        type(points_class),   intent(IN)  :: pts1, pts2 
        real(dp),             intent(IN)  :: shepard_exponent
        real(dp), optional :: lat_lim 

        real(dp), parameter :: DIST_ZERO_OFFSET = 1.0_dp  ! Change dist of zero to 1 m
        integer :: i, i1, kc, k
        real(dp) :: x, y, lon, lat
        real(dp) :: dist, lat_limit

        real :: start, finish

        map%i        = ERR_IND 
        map%dist     = ERR_DIST  
        map%quadrant = 0 
        map%border   = 0 

        lat_limit = 5.0_dp 
        if (present(lat_lim)) lat_limit = lat_lim 
        write(*,"(a,i12,a,f6.2)") "Total points to calculate=",pts2%npts,"  lat_lim=",lat_limit

        ! For each grid point in the new grid,
        ! Find points within a rough radius,
        ! calculate the distance to the current point
        ! and store in map
        call cpu_time(start)
        do i = 1, pts2%npts
                
            ! Get current xy and latlon coordinates
            x   = pts2%x(i)*pts2%xy_conv
            y   = pts2%y(i)*pts2%xy_conv
            lon = pts2%lon(i)
            lat = pts2%lat(i)

            ! Get distance in meters to current point on grid2
            ! for each point on grid1
            do i1 = 1, pts1%npts

                if ( dabs(pts1%lat(i1)-lat) .le. lat_limit ) then
 
                    if (map%is_same_map .and. map%is_cartesian) then
                        ! Use cartesian values to determine distance
                        dist = cartesian_distance(x,y,pts1%x(i1)*pts1%xy_conv,pts1%y(i1)*pts1%xy_conv)

                    else
                        ! Use planetary (latlon) values
                        !dist = spherical_distance(map%planet%a,map%planet%f,lon,lat,pts1%lon(i1),pts1%lat(i1))
                        dist = planet_distance(map%planet%a,map%planet%f,lon,lat,pts1%lon(i1),pts1%lat(i1))

                    end if 

                    ! Make sure no zero distances exist!
                    if (dist .lt. DIST_ZERO_OFFSET) dist = DIST_ZERO_OFFSET

                    do kc = 1, map%nmax
                        if (dist .lt. map%dist(i,kc)) exit
                    end do 

                    if (kc .le. map%nmax) then 

                        if (kc .le. map%nmax-1) then 
                            map%dist(i,kc:map%nmax)     = cshift(map%dist(i,kc:map%nmax),-1)
                            map%i(i,kc:map%nmax)        = cshift(map%i(i,kc:map%nmax),-1)
                            map%quadrant(i,kc:map%nmax) = cshift(map%quadrant(i,kc:map%nmax),-1)
                            map%border(i,kc:map%nmax)   = cshift(map%border(i,kc:map%nmax),-1)
                        end if 

                        map%dist(i,kc)   = dist 
                        map%i(i,kc)      = i1 
                        map%border(i,kc) = pts1%border(i1)

                        ! Get quadrants of neighbors
                        if (map%is_same_map .and. map%is_cartesian) then
                            ! Use cartesian points to determine quadrants
                            map%quadrant(i,kc) = quadrant_cartesian(x,y,pts1%x(i1)*pts1%xy_conv,pts1%y(i1)*pts1%xy_conv)
                        else
                            ! Use planetary (latlon) points
                            map%quadrant(i,kc) = quadrant_latlon(lon,lat,pts1%lon(i1),pts1%lat(i1))
                        end if 

                    end if 
                end if 

            end do

            ! Output every 1000 rows to check progress
            if (mod(i,1000)==0) write(*,*) "  ",i, " / ",pts2%npts,"   : ",map%dist(i,1)
        end do

        call cpu_time(finish)
        write(*,"(a,a,f7.2)") "map_calc_weights:: "//trim(map%name1)//" => "//trim(map%name2)//": ", &
                              "Calculation time (min.) =", (finish-start)/60.0_dp

        ! Also calculate shephard weights
        map%weight = 1.0_dp / (map%dist**shepard_exponent)

        return
    end subroutine map_calc_weights

    function compare_map_points_points(pts1,pts2) result(same_map)

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

    end function compare_map_points_points

    function compare_map_grid_grid(grid1,grid2) result(same_map)

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

    end function compare_map_grid_grid

    function compare_map_map_grid(map1,grid2) result(same_map)

        implicit none 

        type(map_class)  :: map1 
        type(grid_class) :: grid2
        logical :: same_map 

        ! Check if both maps use the same projection  
        if (map1%is_projection .and. grid2%is_projection        .and. &
            trim(map1%mtype)       .eq. trim(grid2%mtype)       .and. &
            same_projection(map1%proj,grid2%proj) ) then 
            
            ! Both maps come from the same projection
            same_map = .TRUE. 

        else if (map1%is_cartesian .and. .not. map1%is_projection .and. &
                 grid2%is_cartesian .and. .not. grid2%is_projection) then 
            ! Both maps are generic cartesian grids (assume origin is the same)
            same_map = .TRUE. 

        else if (.not. map1%is_cartesian .and. .not. grid2%is_cartesian) then 
            ! Both maps are latlon maps
            same_map = .TRUE. 

        else
            same_map = .FALSE.

        end if 

        return 

    end function compare_map_map_grid


    subroutine map_field_grid_grid_integer(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        integer, dimension(:,:), intent(IN)   :: var1
        integer, dimension(:,:), intent(OUT)  :: var2
        integer, dimension(:,:),  intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN), optional  :: mask_pack 
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable   :: mask_pack_vec 
        integer :: nx2, ny2, npts2, npts1 

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 
        npts1 = size(var1,1)*size(var1,2)

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,reshape(dble(var1),[npts1]),var2_vec,mask2_vec, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(int(var2_vec), [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        return

    end subroutine map_field_grid_grid_integer

    subroutine map_field_grid_points_integer(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        integer, dimension(:,:), intent(IN)   :: var1
        integer, dimension(:), intent(OUT)    :: var2
        integer, dimension(:),  intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border  
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec

        integer :: npts1 

        npts1 = size(var1,1)*size(var1,2)
        allocate(var2_vec(size(var2)))

        call map_field_points_points_double(map,name,reshape(dble(var1),[npts1]),var2_vec,mask2, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack)

        var2 = int(var2_vec)

        return

    end subroutine map_field_grid_points_integer

    subroutine map_field_points_grid_integer(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        integer, dimension(:), intent(IN)     :: var1
        integer, dimension(:,:), intent(INOUT)  :: var2
        integer, dimension(:,:),  intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN),  optional :: mask_pack 
        
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable   :: mask_pack_vec 
        integer :: nx2, ny2, npts2

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, (/npts2 /))
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,dble(var1),var2_vec,mask2_vec, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(int(var2_vec),[nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        return

    end subroutine map_field_points_grid_integer

    subroutine map_field_points_points_integer(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)          :: map 
        integer, dimension(:), intent(IN)    :: var1
        integer, dimension(:), intent(INOUT) :: var2
        integer, dimension(:),  intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border  
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec

        allocate(var2_vec(size(var2)))

        call map_field_points_points_double(map,name,dble(var1),var2_vec,mask2, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack)

        var2 = int(var2_vec)

        return

    end subroutine map_field_points_points_integer

    subroutine map_field_grid_grid_double(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack,grid0)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:,:), intent(IN)  :: var1
        real(dp), dimension(:,:), intent(INOUT) :: var2
        integer, dimension(:,:),  intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN),  optional :: mask_pack 
        type(grid_class), optional :: grid0 

        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable   :: mask_pack_vec 
        integer :: nx2, ny2, npts2, npts1 

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 
        npts1 = size(var1,1)*size(var1,2)

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,reshape(var1,[npts1]),var2_vec,mask2_vec, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack_vec)
    
        var2  = reshape(var2_vec, [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        return

    end subroutine map_field_grid_grid_double

    subroutine map_field_grid_points_double(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:,:), intent(IN)  :: var1
        real(dp), dimension(:), intent(INOUT) :: var2
        integer, dimension(:),  intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border  
        real(dp) :: shephard_exponent

        integer :: npts1 

        npts1 = size(var1,1)*size(var1,2)

        call map_field_points_points_double(map,name,reshape(var1,[npts1]),var2,mask2, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack)

        return

    end subroutine map_field_grid_points_double

    subroutine map_field_points_grid_double(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:), intent(IN)    :: var1
        real(dp), dimension(:,:), intent(INOUT) :: var2
        integer,  dimension(:,:), intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN),  optional :: mask_pack 
        
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shephard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable :: mask_pack_vec 
        integer :: nx2, ny2, npts2

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,var1,var2_vec,mask2_vec, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(var2_vec, [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        return

    end subroutine map_field_points_grid_double

    subroutine map_field_points_points_double1(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)
        ! Methods include "radius", "nn" (nearest neighbor) and "quadrant"
        
        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:), intent(IN)    :: var1
        real(dp), dimension(:), intent(INOUT) :: var2
        integer,  dimension(:), intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        logical,  dimension(:), allocatable   :: maskp 
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border 
        logical            :: fill_pts, fill_border
        real(dp) :: shephard_exponent
        real(dp) :: max_distance, missing_val 

        type map_local_type
            integer :: npts 
            real(dp), dimension(:,:), allocatable :: var
            integer,  dimension(:,:), allocatable :: i, quadrant, border  
            real(sp), dimension(:,:), allocatable :: dist, weight
            logical,  dimension(:), allocatable   :: is_border 
            real(sp), dimension(:), allocatable   :: num, denom 
            real(sp), dimension(:), allocatable   :: var2, field
            integer,  dimension(:), allocatable   :: mask2, fieldm 
        end type 
        real(dp), dimension(:), allocatable :: tmp1 

        type(map_local_type) :: maplocal

        integer :: i, k, q, j, ntot, check  
        logical :: found 

        ! Set neighborhood radius to very large value (to include all neighbors)
        ! or to radius specified by user
        max_distance = 1E7_dp
        if (present(radius)) max_distance = radius 

        ! Set grid missing value by default or that that specified by user
        missing_val  = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_val 

        ! By default, grid points with missing values will not be filled in
        fill_pts = .TRUE. 
        if (present(fill)) fill_pts = fill 

        ! By default, border points will not be filled in 
        fill_border = .FALSE. 
        if (present(border)) fill_border = border 

        ! By default, all var2 points are interpolated
        allocate(maskp(size(var2)))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = mask_pack 

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2 = missing_val 

        ! Check to make sure map is the right size
        if (maxval(map%i(1,:)) .gt. size(var1,1)) then
            write(*,*) "Problem! Indices get too big."
            write(*,*) maxval(map%i(1,:)), size(var1,1)
        end if 

        ! ## Eliminate unwanted neighbors ##

        ! Allocate and populated local map
        maplocal%npts = count(maskp)  
        allocate(maplocal%i(maplocal%npts,map%nmax))
        allocate(maplocal%dist(maplocal%npts,map%nmax))
        allocate(maplocal%weight(maplocal%npts,map%nmax))
        allocate(maplocal%var(maplocal%npts,map%nmax))
        allocate(maplocal%quadrant(maplocal%npts,map%nmax))
        allocate(maplocal%border(maplocal%npts,map%nmax))

        allocate(maplocal%is_border(maplocal%npts))
        allocate(maplocal%num(maplocal%npts))
        allocate(maplocal%denom(maplocal%npts))

        do i = 1, map%nmax 
            maplocal%i(:,i)        = pack(map%i(:,i),maskp)
            maplocal%dist(:,i)     = pack(map%dist(:,i),maskp)
            maplocal%weight(:,i)   = pack(map%weight(:,i),maskp)
            maplocal%quadrant(:,i) = pack(map%quadrant(:,i),maskp)
            maplocal%border(:,i)   = pack(map%border(:,i),maskp)
        end do 

        ! Eliminate missing indices
!         maplocal%i        = map%i
!         maplocal%dist     = map%dist         
        where (maplocal%i .eq. ERR_IND) 
            maplocal%i    = 1               ! Set dummy accessible index > 0
            maplocal%dist = max_distance 
        end where 

!         ! Eliminate points that shouldn't be interpolated 
!         do i = 1, map%nmax 
!             where (.not. maskp) maplocal%dist(:,i) = max_distance 
!         end do 
        
!         ! Eliminate neighbors outside of distance limit (in meters)
!         where(maplocal%dist .gt. max_distance) maplocal%dist = max_distance

        ! Eliminate border points 
        maplocal%is_border = .FALSE. 
        if ( (.not. fill_border) ) then 
            where ( sum(maplocal%border,dim=2) .gt. 0.25_dp ) maplocal%is_border = .TRUE. 
            do i = 1, map%nmax
                where (maplocal%is_border) maplocal%dist(:,i) = max_distance 
            end do 
        end if 

        ! Eliminate extra neighbors depending on interpolation method
        if (method .eq. "nn") then 
            maplocal%dist(:,2:map%nmax) = max_distance 
        else if (method .eq. "quadrant") then 
            call quadrant_search(maplocal%dist,maplocal%quadrant,max_distance)
        end if 

        write(*,*) "maskp: ",count(maskp), size(maplocal%var,1), size(maplocal%var,2)

        ! Populate new local var1 with neighbors and weights
!         allocate(tmp1(size(var1,1)))
        maplocal%var = missing_val 
        do i = 1, map%nmax
            tmp1 = var1(maplocal%i(:,i))
            maplocal%var(:,i) = pack(tmp1,maskp)
!             maplocal%var(:,i) = pack(var1(maplocal%i(:,i)),maskp)
            write(*,*) i, "tmp1: ",minval(tmp1),maxval(tmp1)
            write(*,*) i, "var: ",minval(maplocal%var(:,i)),maxval(maplocal%var(:,i))
            
        end do  

!         maplocal%weight = map%weight 
        where (maplocal%dist .ge. max_distance .or. maplocal%var .eq. missing_val) 
            maplocal%weight = 0.d0 
            maplocal%var    = 1.d0 
        end where 

!         write(*,*) "Ready to interpolate..."
!         write(*,*) "var:        ", minval(maplocal%var),maxval(maplocal%var)
!         write(*,*) "weight:     ", minval(maplocal%weight),maxval(maplocal%weight)
!         write(*,*) "var*weight: ", minval(maplocal%var*maplocal%weight), &
!                                    maxval(maplocal%var*maplocal%weight)
!         stop 

        write(*,*) "count maskp: ",count(maskp)
        write(*,*) "range weight: ",minval(maplocal%weight), maxval(maplocal%weight)
        write(*,*) "range var:    ",minval(maplocal%var), maxval(maplocal%var)
        write(*,*) "dim : ", size(sum( maplocal%var*maplocal%weight, dim=2 ),1)

!         maplocal%num   = 0.d0 
!         maplocal%denom = 0.d0 
!         do i = 1, map%nmax 
!             maplocal%num   = maplocal%num   + maplocal%var(:,i)*maplocal%weight(:,i)
!             maplocal%denom = maplocal%denom + maplocal%weight(:,i)
!         end do
!         maplocal%num   = sum( maplocal%var*maplocal%weight, dim=2 ) 
!         maplocal%denom = sum( maplocal%weight, dim=2 ) 
        maplocal%num   = sum( maplocal%var, dim=2 ) 
        maplocal%denom = 1.d0 !sum( maplocal%weight, dim=2 ) 

        where(abs(maplocal%num)   .lt. 1d-20) maplocal%num   = 0.d0 
        where(abs(maplocal%denom) .lt. 1d-20) maplocal%denom = 0.d0 
             
!         write(*,*) "num:   ",minval(maplocal%num,maplocal%denom .gt. 0.d0), maxval(maplocal%num,maplocal%denom .gt. 0.d0)
!         write(*,*) "denom: ",minval(maplocal%denom,maplocal%denom .gt. 0.d0), maxval(maplocal%denom,maplocal%denom .gt. 0.d0)

        allocate(maplocal%var2(maplocal%npts))
        allocate(maplocal%mask2(maplocal%npts))
        allocate(maplocal%field(size(var2,1)))
        allocate(maplocal%fieldm(size(var2,1)))
        maplocal%field  = missing_val  
        maplocal%fieldm = nint(missing_val) 

        maplocal%var2  = pack(var2,maskp)
        maplocal%mask2 = 0 
        where (maplocal%denom .gt. 0.d0) 
            maplocal%var2  = maplocal%num / maplocal%denom
            maplocal%mask2 = 1 
        end where 

        var2  = unpack(maplocal%var2, maskp, maplocal%field  ) 
        if (present(mask2)) mask2 = unpack(maplocal%mask2,maskp, maplocal%fieldm )
        
!         if (maxval(maplocal%var2) > 3000d0) then 
!             write(*,*) "Weighted interp finished: "
!             write(*,*) trim(name)," :", minval(var1), maxval(var1)
!             write(*,*) trim(name)," :", minval(maplocal%var2), maxval(maplocal%var2)
!             write(*,*) trim(name)," :", minval(var2), maxval(var2)
!             stop
!         end if  

!         ! Perform weighted interpolations 
!         do i = 1, map%npts 

!             ntot = count(maplocal%weight .gt. 0.d0)

!             if ( ntot .ge. 1) then 

!                 ! Calculate the weighted average
!                 var2(i)  = weighted_ave(var1(maplocal%i(i,:)),dble(maplocal%weight(i,:)))
!                 mask2(i) = 1

!             else
!                 ! If no neighbors exist, field not mapped here.
!                 mask2(i) = 0 

!             end if 

!             write(*,*) i 
!         end do 

!         write(*,*) "Done allocating."
!         stop 

!         ! Fill missing points with nearest neighbor if desired
!         ! Note, will not necessarily fill ALL points, if 
!         ! no neighbor within nmax can be found without a missing value
!         if ( fill_pts .and. var2(i) .eq. missing_val) then  
!             do k = 1, map%nmax 
!                 if (var1(map%i(i,k)) .ne. missing_val) then
!                     var2(i)  = var1(map%i(i,k))
!                     mask2(i) = 2 
!                     exit 
!                 end if 
!             end do 
!         end if 

        return 

    end subroutine map_field_points_points_double1

    subroutine quadrant_search(dist,quadrant,max_distance)
        ! Set distances to max distance if a neighbor in a given quadrant
        ! has been found 

        implicit none 

        real(sp), dimension(:,:), intent(INOUT) :: dist 
        integer,  dimension(:,:), intent(IN)    :: quadrant 
!         logical,  dimension(:), intent(IN)      :: mask 
        real(dp), intent(IN)                    :: max_distance 
        integer :: i, q, k 
        logical :: found 

        ! For quadrant method, limit the number of neighbors to 
        ! 4 points in different quadrants

        do i = 1, size(dist,1)

!             if (mask(i)) then 

                do q = 1, 4 
                    found = .FALSE. 
                    do k = 1, size(dist,2)
                        if (dist(i,k) .lt. max_distance .and. quadrant(i,k) .eq. q) then 
                            if (found) then 
                                dist(i,k) = max_distance 
                            else
                                found = .TRUE. 
                            end if 
                        end if 
                    end do 
                end do 

!             end if

        end do 

        return

    end subroutine quadrant_search

    subroutine map_field_points_points_double(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)
        ! Methods include "radius", "nn" (nearest neighbor) and "quadrant"
        
        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:), intent(IN)    :: var1
        real(dp), dimension(:), intent(INOUT) :: var2
        integer,  dimension(:), intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        logical,  dimension(:), allocatable   :: maskp 
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border 
        logical            :: fill_pts, fill_border
        real(dp) :: shephard_exponent
        real(dp) :: max_distance, missing_val 
        real(dp), dimension(:), allocatable   :: weight_neighb, v_neighb
        real(dp), dimension(:), allocatable   :: v_neighb_tmp 
        integer,  dimension(:), allocatable   :: mask2_local
        integer :: i, k, q, j, ntot, check  
        logical :: found 

        ! Set neighborhood radius to very large value (to include all neighbors)
        ! or to radius specified by user
        max_distance = 1E7_dp
        if (present(radius)) max_distance = radius 

        ! Set grid missing value by default or that that specified by user
        missing_val  = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value 

        ! By default, fill in grid points with missing values
        fill_pts = .TRUE. 
        if (present(fill)) fill_pts = fill 

        ! By default, border points will not be filled in 
        fill_border = .FALSE. 
        if (present(border)) fill_border = border 

        ! By default, all var2 points are interpolated
        allocate(maskp(size(var2)))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = mask_pack 

        allocate(v_neighb(map%nmax),weight_neighb(map%nmax),v_neighb_tmp(map%nmax))
        
        ! Initialize mask to show which points have been mapped
        allocate(mask2_local(map%npts))
        mask2_local = 0 

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2 = missing_val 

!         ! Check to make sure map is the right size
!         if (maxval(map%i(1,:)) .gt. size(var1,1)) then
!             write(*,*) "Problem! Indices get too big."
!             write(*,*) maxval(map%i(1,:)), size(var1,1)
!         end if 

        do i = 1, map%npts 

            if (maskp(i)) then ! Only perform calculations for packing mask points

            ! Initialize neighbor variable values from input
            v_neighb = missing_val 
            do k = 1, map%nmax 
                if (map%i(i,k) .gt. 0) v_neighb(k) = var1(map%i(i,k))
            end do 

            ! Eliminate neighbors outside of distance limit (in meters)
            where(map%dist(i,:) .gt. max_distance) v_neighb = missing_val 

            ! Skip remaining calculations if no neighbors are found
            check = count(.not. v_neighb .eq. missing_val) 
            if (check .gt. 0) then 

                ! If method == nn (nearest neighbor), limit neighbors to 1
                if (method .eq. "nn") then
!                     v_neighb(2:map%nmax) = missing_val

                    found = .FALSE. 
                    do k = 1, map%nmax 
                        if (v_neighb(k) .ne. missing_val) then 
                            if (found) then 
                                v_neighb(k) = missing_val 
                            else
                                found = .TRUE. 
                            end if 
                        end if 
                    end do

                else if (method .eq. "quadrant") then 
                    ! For quadrant method, limit the number of neighbors to 
                    ! 4 points in different quadrants
                    do q = 1, 4
                        found = .FALSE. 
                        do k = 1, map%nmax 
                            if (v_neighb(k) .ne. missing_val .and. map%quadrant(i,k) .eq. q) then 
                                if (found) then 
                                    v_neighb(k) = missing_val 
                                else
                                    found = .TRUE. 
                                end if 
                            end if 
                        end do 
                    end do 

                end if 

                ! Check number of neighbors available for calculations
                ntot = count(v_neighb .ne. missing_val)

                ! Check if a large fraction of neighbors are border points
                ! (if so, do not interpolate here)
                if ( (.not. fill_border) .and. ntot .gt. 0) then 
                    if ( sum(map%border(i,1:ntot))/dble(ntot) .gt. 0.25_dp ) ntot = 0
                end if 

                ! Fill in temp neighbors with valid values when available
                v_neighb_tmp  = v_neighb 
                v_neighb      = missing_val 
                weight_neighb = 0.0_dp 

                ! Reinitialize temp neighbor and weight values so that they appear in order (1:ntot)
                if (ntot .gt. 0) then 
                    q = 0 
                    do k = 1, map%nmax 
                        if (v_neighb_tmp(k) .ne. missing_val) then
                            q = q+1
                            weight_neighb(q) = map%weight(i,k)
                            v_neighb(q)      = var1(map%i(i,k))
                        end if 
                    end do 
                end if           


                if ( ntot .gt. 1) then 

                    ! Calculate the weighted average
                    var2(i)  = weighted_ave(v_neighb(1:ntot),weight_neighb(1:ntot))
                    mask2_local(i) = 1

                else if (ntot .eq. 1) then
                    var2(i)  = v_neighb(1)
                    mask2_local(i) = 1 

                else
                    ! If no neighbors exist, field not mapped here.
                    mask2_local(i) = 0  

                end if 

                ! Fill missing points with nearest neighbor if desired
                ! Note, will not necessarily fill ALL points, if 
                ! no neighbor within nmax can be found without a missing value
                if ( fill_pts .and. var2(i) .eq. missing_val) then 
                    do k = 1, map%nmax  
                        if (var1(map%i(i,k)) .ne. missing_val) then
                            var2(i)  = var1(map%i(i,k))
                            mask2_local(i) = 2 
                            exit 
                        end if
                    end do 
                end if 

            end if ! End of neighbor checking if-statement 
            end if ! End packing mask check if-statement 

        end do 

        where (dabs(var2) .lt. 1d-20) var2 = 0.d0 

        !write(*,*) "Mapped field: "//trim(name)
!         if (count(var2 .eq. missing_val) .gt. 0) &
!             write(*,*) "   **missing points remaining: ", count(var2 .eq. missing_val)

!         if (count(mask2_local .eq. 0) .gt. 0 .and. .not. fill_pts) then 
!             write(*,*) "Warning, array contains non-interpolated points."
!             write(*,*) "Ensure that it was already properly intialized with data."
!         end if 
        
        ! If interpolation mask available, send to output
        if (present(mask2)) mask2 = mask2_local 

        return
    end subroutine map_field_points_points_double

    function n_quadrants(quadrant)

        implicit none 

        integer, intent(IN) :: quadrant(:) 
        integer :: n_quadrants, q, k 

        n_quadrants = 0 
        do q = 1, 4
            do k = 1, size(quadrant)
                if (quadrant(k) .eq. q) then 
                    n_quadrants = n_quadrants + 1
                    exit 
                end if 
            end do 
        end do 

        return

    end function n_quadrants


    subroutine points_density(density,pts,max_distance)
        ! Get the density of neighbors
        ! for each point in a set within a given radius
        ! NOT TESTED, 2013-10-07, ajr 

        implicit none 

        real(dp) :: density(:), max_distance
        type(points_class) :: pts 

        real(dp), dimension(:), allocatable :: dist

        integer :: i, i1 

        allocate(dist(pts%npts))

        do i = 1, pts%npts 

            dist = ERR_DIST

            do i1 = 1, pts%npts 
                dist(i1) = planet_distance(pts%planet%a,pts%planet%f,pts%lon(i),pts%lat(i),pts%lon(i1),pts%lat(i1))
                density(i) = count(dist .le. max_distance) / (max_distance*max_distance*pi)
            end do 
        end do 

        return 

    end subroutine points_density

    function map_filename(map,fldr)
        ! Output the standard map filename with input folder name
        implicit none 

        type(map_class),  intent(IN) :: map 
        character(len=*), intent(IN) :: fldr 
        character(len=256) :: map_filename
        character(len=2) :: char2

        if (map%nmax .ge. 10) then
            write(char2,"(i2)") map%nmax
        else
            write(char2,"(i1,i1)") 0, map%nmax 
        end if 

        map_filename = trim(fldr)//"/map_"//trim(map%name1)//"_"//trim(map%name2)//"_"//trim(char2)//".nc"

        return
    end function map_filename

    subroutine map_write(map,fldr)

        implicit none 

        type(map_class),  intent(IN) :: map 
        character(len=*), intent(IN) :: fldr 
        character(len=256) :: fnm 
        character(len=128) :: dim1, dim2 

        fnm = map_filename(map,fldr)

        ! Create the netcdf file and the dimension variables
        call nc_create(fnm)
        call nc_write_attr(fnm,"title","Mapping "//trim(map%name1)//" => "//trim(map%name2))

        ! Write generic dimensions
        call nc_write_dim(fnm,"point",    x=1,nx=map%npts,units="n")
        call nc_write_dim(fnm,"neighbor", x=1,nx=map%nmax,units="n")
        call nc_write_dim(fnm,"parameter",x=1,units="")
        call nc_write_dim(fnm,"planetpar",x=1,nx=3,units="")
        call nc_write_dim(fnm,"projpar",  x=1,nx=5,units="")

        ! Write grid/vector specific dimensions and variables
        if (map%is_grid) then
            ! Write variables in a gridded format

            if (trim(map%mtype) .eq. "latlon") then 
                dim1 = "lon"
                dim2 = "lat"
                call nc_write_dim(fnm,dim1,x=map%G%x)
                call nc_write_dim(fnm,dim2,x=map%G%y)
            else 
                dim1 = "xc"
                dim2 = "yc" 
                call nc_write_dim(fnm,dim1,x=map%G%x,units=trim(map%units))
                call nc_write_dim(fnm,dim2,x=map%G%y,units=trim(map%units))
            end if 

            call nc_write(fnm,"x2D",  map%x,  dim1=dim1,dim2=dim2)
            call nc_write(fnm,"y2D",  map%y,  dim1=dim1,dim2=dim2)
            call nc_write(fnm,"lon2D",map%lon,dim1=dim1,dim2=dim2)
            call nc_write(fnm,"lat2D",map%lat,dim1=dim1,dim2=dim2)

            ! Warning about size 
            ! size of nx*ny*nmax < 1579220 breaks on the cluster !!!!
            if (map%G%nx*map%G%ny*map%nmax > 1500000) then 
                write(*,*) "Warning: map size is very large (nx*ny*nmax): ",map%G%nx*map%G%ny*map%nmax 
                write(*,*) "  It may cause a segmentation fault. If so, reduce the number of neighbors."
            end if 

            ! Write the map information (in grid format)
            call nc_write(fnm,"i",       map%i,       dim1=dim1,dim2=dim2,dim3="neighbor")
            call nc_write(fnm,"dist",    map%dist,    dim1=dim1,dim2=dim2,dim3="neighbor")
            call nc_write(fnm,"weight",  map%weight,  dim1=dim1,dim2=dim2,dim3="neighbor")
            call nc_write(fnm,"quadrant",map%quadrant,dim1=dim1,dim2=dim2,dim3="neighbor")
            call nc_write(fnm,"border",  map%border,  dim1=dim1,dim2=dim2,dim3="neighbor")

            ! Write grid specific parameters
            call nc_write(fnm,"nx",map%G%nx,dim1="parameter")
            call nc_write(fnm,"ny",map%G%ny,dim1="parameter")

        else
            ! Write variables in a vector format
            call nc_write(fnm,"x",       map%x,       dim1="point",dim2="neighbor")
            call nc_write(fnm,"y",       map%y,       dim1="point",dim2="neighbor")
            call nc_write(fnm,"lon",     map%lon,     dim1="point",dim2="neighbor")
            call nc_write(fnm,"lat",     map%lat,     dim1="point",dim2="neighbor")
            call nc_write(fnm,"i",       map%i,       dim1="point",dim2="neighbor")
            call nc_write(fnm,"dist",    map%dist,    dim1="point",dim2="neighbor")
            call nc_write(fnm,"weight",  map%weight,  dim1="point",dim2="neighbor")
            call nc_write(fnm,"quadrant",map%quadrant,dim1="point",dim2="neighbor")
            call nc_write(fnm,"border",  map%border,  dim1="point",dim2="neighbor")

        end if 

        ! Write generic map parameters
        call nc_write(fnm,"mtype",        map%mtype)
        call nc_write(fnm,"units",        map%units)
        call nc_write(fnm,"is_cartesian", map%is_cartesian, dim1="parameter")
        call nc_write(fnm,"is_projection",map%is_projection,dim1="parameter")
        call nc_write(fnm,"is_lon180",    map%is_lon180,    dim1="parameter")
        call nc_write(fnm,"is_same_map",  map%is_same_map,  dim1="parameter")
        call nc_write(fnm,"is_grid",      map%is_grid,      dim1="parameter")
        call nc_write(fnm,"npts",         map%npts,         dim1="parameter")  
        call nc_write(fnm,"nmax",         map%nmax,         dim1="parameter")        
        call nc_write(fnm,"xy_conv",      map%xy_conv,      dim1="parameter") 

        if (map%is_projection .or. .not. map%is_cartesian) then 
            call nc_write(fnm,"planet_name",map%planet%name)
            call nc_write(fnm,"planet_info",[map%planet%a, map%planet%f, map%planet%R], &
                          dim1="planetpar")
        end if 

        if (map%is_projection) then 
            call nc_write(fnm,"proj_name",map%proj%name)
            call nc_write(fnm,"proj_info",[map%proj%lambda, map%proj%phi, map%proj%alpha, &
                                 map%proj%x_e, map%proj%y_n], dim1="projpar")    
        end if 
        
        write(*,*) "Map written to file: "//trim(fnm)

        return 
    end subroutine map_write

    subroutine map_read(map,fldr)

        implicit none 

        type(map_class),  intent(INOUT) :: map 
        character(len=*), intent(IN)    :: fldr 
        character(len=256) :: fnm  

        integer,          dimension(:,:,:), allocatable :: tmpi
        real(dp), dimension(:,:,:), allocatable :: tmpd
        real(dp) :: tmp5(5) 

        fnm = map_filename(map,fldr)

        ! Deallocate all map fields initially
        if (allocated(map%G%x))      deallocate(map%G%x)
        if (allocated(map%G%y))      deallocate(map%G%y)
        if (allocated(map%x))        deallocate(map%x) 
        if (allocated(map%y))        deallocate(map%y) 
        if (allocated(map%lon))      deallocate(map%lon) 
        if (allocated(map%lat))      deallocate(map%lat) 
        if (allocated(map%i))        deallocate(map%i) 
        if (allocated(map%dist))     deallocate(map%dist) 
        if (allocated(map%weight))   deallocate(map%weight) 
        if (allocated(map%quadrant)) deallocate(map%quadrant) 
        if (allocated(map%border))   deallocate(map%border) 

        ! Read generic map parameters
        call nc_read(fnm,"mtype",        map%mtype        )
        call nc_read(fnm,"units",        map%units        )
        call nc_read(fnm,"is_cartesian", map%is_cartesian )
        call nc_read(fnm,"is_projection",map%is_projection)
        call nc_read(fnm,"is_lon180",    map%is_lon180    )
        call nc_read(fnm,"is_same_map",  map%is_same_map  )
        call nc_read(fnm,"is_grid",      map%is_grid      )
        call nc_read(fnm,"npts",         map%npts         ) 
        call nc_read(fnm,"nmax",         map%nmax         )        
        call nc_read(fnm,"xy_conv",      map%xy_conv      )        
        
        if (map%is_projection .or. .not. map%is_cartesian) then 
            call nc_read(fnm,"planet_name",map%planet%name)
            call nc_read(fnm,"planet_info",tmp5(1:3))
            map%planet%a = tmp5(1)
            map%planet%f = tmp5(2)
            map%planet%R = tmp5(3)
        end if 

        if (map%is_projection) then 
            call nc_read(fnm,"proj_name",map%proj%name)
            call nc_read(fnm,"proj_info",tmp5(1:5))
            map%proj%lambda = tmp5(1)
            map%proj%phi    = tmp5(2)
            map%proj%alpha  = tmp5(3) 
            map%proj%x_e    = tmp5(4) 
            map%proj%y_n    = tmp5(5) 
        end if 

        ! Allocate map fields 
        allocate(map%x(map%npts),map%y(map%npts))
        allocate(map%lon(map%npts),map%lat(map%npts))
        allocate(map%i(map%npts,map%nmax),map%dist(map%npts,map%nmax))
        allocate(map%weight(map%npts,map%nmax),map%quadrant(map%npts,map%nmax))
        allocate(map%border(map%npts,map%nmax))

        ! Read grid/vector specific dimensions and variables
        if (map%is_grid) then
            ! Read variables in a gridded format

            call nc_read(fnm,"nx",map%G%nx)
            call nc_read(fnm,"ny",map%G%ny)

            ! Allocate grid fields
            allocate(map%G%x(map%G%nx))
            allocate(map%G%y(map%G%ny))

            if (trim(map%mtype) .eq. "latlon") then 
                call nc_read(fnm,"lon",map%G%x)
                call nc_read(fnm,"lat",map%G%y)
            else 
                call nc_read(fnm,"xc",map%G%x)
                call nc_read(fnm,"yc",map%G%y)
            end if 

            if (allocated(tmpd)) deallocate(tmpd)
            allocate(tmpd(map%G%nx,map%G%ny,1))

            call nc_read(fnm,"x2D",tmpd(:,:,1))
            map%x   = reshape(tmpd,(/map%npts/))
            call nc_read(fnm,"y2D",tmpd(:,:,1))
            map%y   = reshape(tmpd,(/map%npts/))
            call nc_read(fnm,"lon2D",tmpd(:,:,1))
            map%lon = reshape(tmpd,(/map%npts/))
            call nc_read(fnm,"lat2D",tmpd(:,:,1))
            map%lat = reshape(tmpd,(/map%npts/))
            
            if (allocated(tmpi)) deallocate(tmpi)
            if (allocated(tmpd)) deallocate(tmpd)
            allocate(tmpi(map%G%nx,map%G%ny,map%nmax),tmpd(map%G%nx,map%G%ny,map%nmax))

            call nc_read(fnm,"i",tmpi)
            map%i = reshape(tmpi,(/map%npts,map%nmax/))
            call nc_read(fnm,"dist",tmpd)
            map%dist = reshape(tmpd,(/map%npts,map%nmax/))
            call nc_read(fnm,"weight",tmpd)
            map%weight = reshape(tmpd,(/map%npts,map%nmax/))
            call nc_read(fnm,"quadrant",tmpd)
            map%quadrant = reshape(tmpd,(/map%npts,map%nmax/))
            call nc_read(fnm,"border",tmpi)
            map%border = reshape(tmpi,(/map%npts,map%nmax/))
            
        else
            ! Read variables in a vector format

            call nc_read(fnm,       "x",map%x)
            call nc_read(fnm,       "y",map%y)
            call nc_read(fnm,     "lon",map%lon)
            call nc_read(fnm,     "lat",map%lat)
            call nc_read(fnm,       "i",map%i)
            call nc_read(fnm,    "dist",map%dist)
            call nc_read(fnm,  "weight",map%weight)
            call nc_read(fnm,"quadrant",map%quadrant)
            call nc_read(fnm,  "border",map%border)

        end if 

        write(*,*) "Map read from file: "//trim(fnm)

        return 
    end subroutine map_read

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

    subroutine grid_write(grid,fnm,xnm,ynm,create)

        implicit none 
        type(grid_class) :: grid 
        character(len=*) :: fnm,xnm,ynm
        logical :: create  

        ! Create the netcdf file if desired
        if (create) then 
            call nc_create(fnm)
        
            ! Add grid axis variables to netcdf file
            call nc_write_dim(fnm,xnm,x=grid%G%x,units=grid%units)
            call nc_write_dim(fnm,ynm,x=grid%G%y,units=grid%units)
        end if 

        ! Add projection information if needed
        if (grid%is_projection) &
            call nc_write_map(fnm,grid%mtype,grid%proj%lambda,phi=grid%proj%phi, &
                              x_e=grid%proj%x_e,y_n=grid%proj%y_n)

        if (grid%is_projection .or. grid%is_cartesian) then 
            call nc_write(fnm,"x2D",grid%x,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
            call nc_write(fnm,"y2D",grid%y,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
        end if 
        if (.not. (grid%is_cartesian .and. .not. grid%is_projection)) then 
            call nc_write(fnm,"lon2D",grid%lon,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
            call nc_write(fnm,"lat2D",grid%lat,dim1=xnm,dim2=ynm,grid_mapping=grid%name)
        end if 

        call nc_write(fnm,"area",  grid%area,  dim1=xnm,dim2=ynm,grid_mapping=grid%name)
        call nc_write(fnm,"border",grid%border,dim1=xnm,dim2=ynm,grid_mapping=grid%name)

        return
    end subroutine grid_write

    subroutine points_write(pts,fnm,xnm,ynm,create)

        implicit none 
        type(points_class) :: pts 
        character(len=*) :: fnm,xnm,ynm
        logical :: create  

        ! Create the netcdf file if desired
        if (create) then 
            call nc_create(fnm)
        
            ! Add axis variable to netcdf file
            call nc_write_dim(fnm,"point",x=1,nx=pts%npts,units="n")
        end if 

        ! Add projection information if needed
        if (pts%is_projection) &
            call nc_write_map(fnm,pts%mtype,pts%proj%lambda,phi=pts%proj%phi,&
                              x_e=pts%proj%x_e,y_n=pts%proj%y_n)

        if (pts%is_projection .or. pts%is_cartesian) then 
            call nc_write(fnm,xnm,pts%x,dim1="point",grid_mapping=pts%name)
            call nc_write(fnm,ynm,pts%y,dim1="point",grid_mapping=pts%name)
        end if 
        if (.not. (pts%is_cartesian .and. .not. pts%is_projection)) then 
            call nc_write(fnm,"lon",pts%lon,dim1="point",grid_mapping=pts%name)
            call nc_write(fnm,"lat",pts%lat,dim1="point",grid_mapping=pts%name)
        end if 

!         call nc_write(fnm,"area",  pts%area,  dim1="point",grid_mapping=pts%name)
!         call nc_write(fnm,"border",pts%border,dim1="point",grid_mapping=pts%name)

        return
    end subroutine points_write

    subroutine grid_print(grid)

        implicit none 

        type(grid_class), intent(IN) :: grid
        type(points_class) :: pts 

        write(*,*) "== Grid summary   ========================="
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
            write(*,*) "Projection information"
            write(*,"(a16,g12.5)")     "a = ",     pts%proj%a
            write(*,"(a16,g12.5)")     "f = ",     pts%proj%f
            write(*,"(a16,g12.5)")     "lambda = ",pts%proj%lambda
            write(*,"(a16,g12.5)")     "phi = ",   pts%proj%phi
            write(*,"(a16,g12.5)")     "alpha = ", pts%proj%alpha
            write(*,"(a16,g12.5)")     "x_e = ",   pts%proj%x_e
            write(*,"(a16,g12.5)")     "y_n = ",   pts%proj%y_n
        end if 
        
        write(*,*)

        return
    end subroutine points_print
    
    subroutine map_print(map)

        implicit none 

        type(map_class), intent(IN) :: map 

        write(*,*) "== Mapping summary ========================="
        write(*,*) trim(map%name1)," => ",trim(map%name2)
        write(*,"(a16,g12.5)")     "             a = ",map%planet%a
        write(*,"(a16,g12.5)")     "             f = ",map%planet%f
        write(*,"(a16,i8,i6)")    "     npts,nmax = ",map%npts,map%nmax 
        if (map%is_grid) then
            write(*,*) "  Grid axis information"
            write(*,"(a16,i8,i6)") "         nx,ny = ",map%G%nx, map%G%ny 
        end if 

        write(*,"(a,g11.4,a)") " ** Size in memory ~", &
                   ( (map%npts*map%nmax)*2.d0*4.d0 + (map%npts*map%nmax)*3.d0*4.d0 + &
                     (map%npts)*2.d0*8.d0 )  *7.6294d-6,"Mb"
        ! 2 real arrays (4 bytes per value): dist, weight
        ! 3 integer array (4 bytes per value): i, quadrant, border
        ! 4 real arrays (8 bytes per value): x, y, lon, lat

        write(*,*)
        write(*,*) 

        return
    end subroutine map_print

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

