module coordinates_mapping 
    
    use coordinates 

    use oblimap_projection_module
    use planet 
    use polygons 
    use ncio

    use gaussian_filter  
    use interp2D 
    use mod_toms526 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = -9999.0_dp 
    
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    type pt_wts_class 
        integer :: n  
        integer,  allocatable :: i(:)               ! Length of pts1/grid1 neighbors
        real(sp), allocatable :: x(:), y(:)
        real(sp), allocatable :: dist(:), weight(:), area(:)
        integer,  allocatable :: quadrant(:), border(:) 
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

        ! Neighbor info (allocate to size npts)
        type(pt_wts_class), allocatable :: map(:)

    end type

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
        module procedure compare_map_map_grid 
    end interface

    private 
    public :: compare_map
    public :: map_class, map_init, map_field, map_print
    public :: map_field_conservative_map1

contains 

    subroutine map_init_grid_grid(map,grid1,grid2,max_neighbors,lat_lim,dist_max,fldr,load,save)
        ! Generate mapping weights from grid1 to grid2

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid1, grid2 
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        real(dp), optional :: dist_max           ! maximum distance within which to include neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        ! Initialize map grid axis info (since mapping to grid2)
        map%is_grid = .TRUE. 
        map%G       = grid2%G 

        call grid_to_points(grid1,pts1)
        call grid_to_points(grid2,pts2)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)

        return 

    end subroutine map_init_grid_grid

    subroutine map_init_grid_points(map,grid1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)
        ! Generate mapping weights from grid to set of points

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid1 
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        real(dp), optional :: dist_max           ! maximum distance within which to include neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        map%is_grid = .FALSE. 
        call grid_to_points(grid1,pts1)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)

        return 

    end subroutine map_init_grid_points

    subroutine map_init_points_grid(map,pts1,grid2,max_neighbors,lat_lim,dist_max,fldr,load,save)
        ! Generate mapping weights from set of points to grid

        implicit none 

        type(points_class) :: pts1, pts2
        type(grid_class)   :: grid2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        real(dp), optional :: dist_max           ! maximum distance within which to include neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        ! Initialize map grid axis info (since mapping to grid2)
        map%is_grid = .TRUE.  
        map%G       = grid2%G 

        ! Convert grid2 to points for map initialization
        call grid_to_points(grid2,pts2)
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)

        return 

    end subroutine map_init_points_grid 

    subroutine map_init_points_points(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)
        ! Generate mapping weights from set of points to another set of points

        implicit none 

        type(points_class) :: pts1, pts2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        real(dp), optional :: dist_max           ! maximum distance within which to include neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated

        map%is_grid = .FALSE. 
        call map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)

        return 

    end subroutine map_init_points_points

    subroutine map_init_internal(map,pts1,pts2,max_neighbors,lat_lim,dist_max,fldr,load,save)
        ! Generate mapping weights between sets of pointss

        implicit none 

        type(points_class) :: pts1, pts2
        type(map_class) :: map 
        integer :: max_neighbors                 ! maximum number of neighbors to allocate 
        real(dp), optional :: lat_lim            ! Latitude limit to search for neighbors
        real(dp), optional :: dist_max           ! maximum distance within which to include neighbors
        character(len=*), optional :: fldr       ! Directory in which to save/load map
        logical, optional :: load                ! Whether loading is desired if map exists already
        logical, optional :: save                ! Whether to save map to file after it's generated
        logical :: load_file, save_file, fldr_exists, file_exists 
        character(len=256) :: mapfldr 

        integer :: i 
        integer, allocatable :: ii(:) 
        logical :: pts_is_latlon  

        pts_is_latlon = (.not. pts1%is_cartesian .and. .not. pts2%is_cartesian)

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
        map%is_same_map = compare_coord(pts1,pts2)

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

            allocate(map%map(map%npts))

            ! Calculate map weights (time consuming!)
            call map_calc_weights(map,pts1,pts2,2.0_dp,lat_lim,dist_max)

            ! If the grids are compatible, calculate area weighting too
            if (same_projection(pts1%proj,pts2%proj) .and. .not. pts_is_latlon) then 
                ! Only valid for cartesian grids on the same projection right now 

                write(*,*) "Calculating conservative weights ..."
                do i = 1, pts2%npts
                    if (map%map(i)%dist(1) .lt. ERR_DIST) then 
                        ii = map%map(i)%i 
                        map%map(i)%area = calc_weights_interpconserv1(x=real(pts1%x(ii)*pts1%xy_conv), &
                                            y=real(pts1%y(ii)*pts1%xy_conv), &
                                            dx=real(pts1%dx(ii)*pts1%xy_conv),dy=real(pts1%dy(ii)*pts1%xy_conv), &
                                            xout=real(pts2%x(i)*pts2%xy_conv),yout=real(pts2%y(i)*pts2%xy_conv), &
                                            dxout=real(pts2%dx(i)*pts2%xy_conv),dyout=real(pts2%dy(i)*pts2%xy_conv))
                        map%map(i)%area = map%map(i)%area / (pts1%xy_conv*pts1%xy_conv)   ! Convert back to axis-units of pts1
                        
!                         map%map(i)%area = calc_weights_interpconserv1(x=real(pts1%x(ii)), &
!                                             y=real(pts1%y(ii)), &
!                                             dx=real(pts1%dx(ii)),dy=real(pts1%dy(ii)), &
!                                             xout=real(pts2%x(i)),yout=real(pts2%y(i)), &
!                                             dxout=real(pts2%dx(i)),dyout=real(pts2%dy(i)))
!                         map%map(i)%area = map%map(i)%area
                    
                    end if 
                    
                    ! Check progress
                    if (mod(i,100)==0) write(*,"(a,i10,a3,i12,a5,2g12.3)")  &
                                "  ",i, " / ",pts2%npts,"   : ", sum(map%map(i)%area), pts2%dx(i)*pts2%dy(i)
                end do 

            end if 


            ! Write new map to file
            if (save_file) call map_write(map,mapfldr) 

        end if 

        ! Print map summary
        call map_print(map)

        return

    end subroutine map_init_internal

    subroutine map_allocate_map(mp,n)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mp 
        integer, intent(IN) :: n  

        mp%n = n  

        ! First deallocate if needed 
        if (allocated(mp%i))        deallocate(mp%i) 
        if (allocated(mp%x))        deallocate(mp%x) 
        if (allocated(mp%y))        deallocate(mp%y)  
        if (allocated(mp%dist))     deallocate(mp%dist) 
        if (allocated(mp%weight))   deallocate(mp%weight) 
        if (allocated(mp%area))     deallocate(mp%area) 
        if (allocated(mp%quadrant)) deallocate(mp%quadrant) 
        if (allocated(mp%border))   deallocate(mp%border) 
        
        allocate(mp%i(n))
        allocate(mp%x(n),mp%y(n))
        allocate(mp%dist(n), mp%weight(n), mp%area(n))
        allocate(mp%quadrant(n), mp%border(n))

        return 

    end subroutine map_allocate_map

    subroutine map_calc_weights(map,pts1,pts2,shepard_exponent,lat_lim,dist_max)
        implicit none 

        type(map_class) :: map 
        type(points_class),   intent(IN)  :: pts1, pts2 
        real(dp),             intent(IN)  :: shepard_exponent
        real(dp),             intent(IN), optional :: lat_lim 
        real(dp),             intent(IN), optional :: dist_max 

        real(dp), parameter :: DIST_ZERO_OFFSET = 1.0_dp  ! Change dist of zero to 1 m
        integer :: i, i1, kc, k, n, n_tot  
        real(dp) :: x, y, lon, lat
        real(dp) :: dist, lat_limit, dist_maximum 

        integer, parameter :: map_nmax = 100   ! No more than 1000 neighbors 
        real(dp) :: map_x(map_nmax), map_y(map_nmax), map_dist(map_nmax)
        integer  :: map_i(map_nmax), map_quadrant(map_nmax), map_border(map_nmax)

        real :: start, finish

        n_tot = 0 

        map%i        = ERR_IND 
        map%dist     = ERR_DIST  
        map%quadrant = 0 
        map%border   = 0 

        ! Limit neighborhood to search 
        lat_limit = 5.0_dp 
        if (present(lat_lim)) lat_limit = lat_lim

        dist_maximum = ERR_DIST 
        if (present(dist_max)) dist_maximum = dist_max 

        write(*,"(a,i12,a,f6.2,a,g12.3)") "Total points to calculate=",pts2%npts, &
                    "  lat_lim=",lat_limit, "  dist_max=",dist_maximum 

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

            ! Reset temp map values 
            map_i        = ERR_IND
            map_x        = 0.0_dp 
            map_y        = 0.0_dp 
            map_dist     = ERR_DIST 
            map_quadrant = 0 
            map_border   = 0 

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

                    do kc = 1, map_nmax
                        if (dist .lt. map_dist(kc)) exit
                    end do 

                    if (kc .le. map%nmax .and. dist .lt. dist_maximum) then 

                        if (kc .le. map%nmax-1) then 
                            map_i(kc:map_nmax)        = cshift(map_i(kc:map_nmax),-1)
                            map_x(kc:map_nmax)        = cshift(map_x(kc:map_nmax),-1)
                            map_y(kc:map_nmax)        = cshift(map_y(kc:map_nmax),-1)
                            map_dist(kc:map_nmax)     = cshift(map_dist(kc:map_nmax),-1)
                            map_quadrant(kc:map_nmax) = cshift(map_quadrant(kc:map_nmax),-1)
                            map_border(kc:map_nmax)   = cshift(map_border(kc:map_nmax),-1)
                        end if 

                        map_i(kc)      = i1 
                        map_x(kc)      = pts1%x(i1)
                        map_y(kc)      = pts1%y(i1) 
                        map_dist(kc)   = dist 
                        map_border(kc) = pts1%border(i1)

                        ! Get quadrants of neighbors
                        if (map%is_same_map .and. map%is_cartesian) then
                            ! Use cartesian points to determine quadrants
                            map_quadrant(kc) = quadrant_cartesian(x,y,pts1%x(i1)*pts1%xy_conv,pts1%y(i1)*pts1%xy_conv)
                        else
                            ! Use planetary (latlon) points
                            map_quadrant(kc) = quadrant_latlon(lon,lat,pts1%lon(i1),pts1%lat(i1))
                        end if 

                    end if 
                end if 

            end do

            ! Store in main map now 
            map%i(i,:)        = map_i(1:map%nmax)
            map%dist(i,:)     = map_dist(1:map%nmax)
            map%weight(i,:)   = 1.0_dp / (map%dist(i,:)**shepard_exponent)
            map%quadrant(i,:) = map_quadrant(1:map%nmax)
            map%border(i,:)   = map_border(1:map%nmax)

            ! Store in new main map now 
            n=count(map_i .ne. ERR_IND)
            if (n .gt. 0) then 
                n = min(n,map%nmax)   ! Limit neighbors to max_neighbors specified if needed
                call map_allocate_map(map%map(i),n=n)
                map%map(i)%i        = map_i(1:n)
                map%map(i)%x        = map_x(1:n)
                map%map(i)%y        = map_y(1:n)
                map%map(i)%dist     = map_dist(1:n)
                map%map(i)%weight   = 1.0_dp / (map%map(i)%dist**shepard_exponent)
                map%map(i)%quadrant = map_quadrant(1:n)
                map%map(i)%border   = map_border(1:n) 
                map%map(i)%area     = 0.0_dp 

            else 
                ! Store filler index
                n = 1 
                call map_allocate_map(map%map(i),n=n)
                map%map(i)%i        = ERR_IND 
                map%map(i)%x        = 0.0_dp 
                map%map(i)%y        = 0.0_dp 
                map%map(i)%dist     = ERR_DIST
                map%map(i)%weight   = 0.0_dp 
                map%map(i)%quadrant = 1
                map%map(i)%border   = 0
                map%map(i)%area     = 0.0_dp
 
            end if 

            n_tot = n_tot + n 

            ! Output every 1000 rows to check progress
            if (mod(i,1000)==0) write(*,*) "  ",i, " / ",pts2%npts,"   : ",map%dist(i,1), " : ", n
        end do

        call cpu_time(finish)
        write(*,"(a,a,f7.2)") "map_calc_weights:: "//trim(map%name1)//" => "//trim(map%name2)//": ", &
                              "Calculation time (min.) =", (finish-start)/60.0_dp

        write(*,*) "map sizes (old, new): ", size(map%i,1)*size(map%i,2), n_tot 

        return

    end subroutine map_calc_weights

    function calc_weights_interpconserv1(x,y,dx,dy,xout,yout,dxout,dyout,latlon) result(area)
        ! Calculate 1st order conservative interpolation for a 
        ! point given a vector of its neighbors (x,y)

        real(sp), intent(IN) :: x(:), y(:), dx(:), dy(:) 
        real(sp), intent(IN) :: xout, yout, dxout, dyout 
        logical,  intent(IN), optional :: latlon  
        real(sp) :: area(size(x,1))

        ! Local variables
        logical  :: is_latlon 
        type(polygon)       :: pol  
!         integer, parameter :: nx = 31, ny = 31, npts = nx*ny
        integer :: nx, ny, npts 
        real(sp) :: x1, y1
        integer  :: npts_in  
        integer :: i, j, now  
        real(sp) :: missing_val
        real(sp) :: area_target 

        is_latlon = .FALSE. 
        if (present(latlon)) is_latlon = latlon 

        if (is_latlon) then 
            ! Latlon is currently not yet supported, set area to 0 
            area = 0.d0 

        else 
            ! Calculate cartesian area of points inside of target point 

            ! Generate polygon representing boundaries of target point 
            pol = create_polygon(real([xout-dxout/2.d0,xout-dxout/2.d0,xout+dxout/2.d0,xout+dxout/2.d0]), &
                                 real([yout-dyout/2.d0,yout+dyout/2.d0,yout+dyout/2.d0,yout-dyout/2.d0]))

            ! Loop over source points and get the area of each source
            ! polygon that is inside of the target polygon 
            ! - Save the area (absolute area, not fraction)
            area_target = dxout*dyout 
            area        = 0.d0 

            do now = 1, size(x) 
                
                nx = max(5,int(dx(now)/dxout))
                ny = max(5,int(dy(now)/dyout))
                npts = nx*ny 

                npts_in   = 0

                do j = 1, ny 
                    do i = 1, nx 
                        x1 = (x(now)-dx(now)/2.d0) + (dx(now))*dble(i-1)/dble(nx) + 0.5d0*1.d0/dble(nx)
                        y1 = (y(now)-dy(now)/2.d0) + (dy(now))*dble(j-1)/dble(ny) + 0.5d0*1.d0/dble(ny) 
                        if (point_in_polygon(real(x1),real(y1),pol)) npts_in = npts_in+1
                    end do
                end do 
                
                area(now) = dble(npts_in)/dble(npts) * dx(now)*dy(now) 

                ! If the source points area adds up to the target cell area, exit loop
                if (sum(area) .ge. area_target) exit 
            end do 

        end if 


        return 

    end function calc_weights_interpconserv1

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
            call nc_write_attr(fnm,"proj_name","method",map%proj%method)
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


    ! === CONSERVATIVE FIELD MAPPING SUBROUTINES === 

    subroutine map_field_conservative_map1(mp,varname,var1,var2,fill,missing_value,mask_pack)
        ! Map a field from grid1 to grid2 using a predefined map of neighbor weights
        ! generated using `map_conserv_init` 

        implicit none 

        type(pt_wts_class), intent(IN)    :: mp(:)          ! Map values grid1=>grid2
        character(len=*),   intent(IN)    :: varname        ! Name of the variable being mapped
        double precision,   intent(IN)    :: var1(:,:)      ! Input variable
        double precision,   intent(INOUT) :: var2(:,:)      ! Output variable
        logical,  optional :: fill
        double precision, intent(IN), optional :: missing_value  ! Points not included in mapping
        logical,          intent(IN), optional :: mask_pack(:,:)

        ! Local variables
        integer :: npts1
        integer :: npts2
        logical          :: fill_pts
        double precision :: missing_val 
        integer :: i, j 
        logical, allocatable  :: maskp(:)
        integer, allocatable :: ii(:) 
        real(dp), allocatable :: area(:)

        real(dp), allocatable :: var1_vec(:), var2_vec(:) 

        npts1 = size(var1,1)*size(var1,2)
        npts2 = size(var2,1)*size(var2,2)

        ! By default, fill in target grid points with missing values
        fill_pts = .TRUE. 
        if (present(fill)) fill_pts = fill 

        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! By default, all var2 points are interpolated
        allocate(maskp(npts2))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = reshape(mask_pack,[npts2])

        ! Store var1 in vector format
        allocate(var1_vec(npts1)) 
        var1_vec = reshape(var1,[npts1])

        ! Store var2 in vector format
        allocate(var2_vec(npts2))
        var2_vec = reshape(var2,[npts2])

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2_vec = missing_val 

        ! Loop over target points
        do i = 1, npts2

            if (maskp(i)) then 
                ! Only interpolate for desired target points 

                ii = mp(i)%i  
                area = mp(i)%area 
                where (var1_vec(ii) .eq. missing_val) area = 0.d0 

                if (sum(area) .gt. 0.d0) then 
                    ! If an interpolation point was found, calculate interpolation 

                    var2_vec(i) = sum(var1_vec(ii)*area,mask=area.gt.0.d0) &
                                  / sum(area,mask=area.gt.0.d0)

                end if 

            end if 

        end do 

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

        if (fill_pts) then 
            call fill_nearest(var2,missing_value=missing_val)
        end if 

        return 

    end subroutine map_field_conservative_map1 

    function calc_grid_total(x,y,var,xlim,ylim) result(tot)

        implicit none 

        double precision, intent(IN)  :: x(:), y(:), var(:,:)
        double precision, intent(IN)  :: xlim(2), ylim(2)   ! Extent over which to calculate total
        double precision :: tot 

        double precision :: dx, dy, d_extra 
        double precision :: weight(size(var,1),size(var,2))
        double precision :: xwt, ywt 
        integer :: i, j 

        ! Determine the resolution of the grid 
        dx = abs(x(2)-x(1)) 
        dy = abs(y(2)-y(1))

        ! Initially set weight to zero 
        weight = 0.d0 

        ! Find weights over the whole domain 
        do j = 1, size(y)
            do i = 1, size(x)

                if (x(i)-dx/2.d0 .le. xlim(1) .and. x(i)+dx/2.d0 .gt. xlim(1)) then 
                    xwt = ( (x(i)+dx/2.d0) - xlim(1) ) / dx 
                else if (x(i)-dx/2.d0 .lt. xlim(2) .and. x(i)+dx/2.d0 .ge. xlim(2)) then 
                    xwt = ( xlim(2) - (x(i)-dx/2.d0) ) / dx 
                else if (x(i)-dx/2.d0 .ge. xlim(1) .and. x(i)+dx/2.d0 .le. xlim(2)) then 
                    xwt = 1.d0 
                else 
                    xwt = 0.d0 
                end if 

                if (y(j)-dy/2.d0 .le. ylim(1) .and. y(j)+dy/2.d0 .gt. ylim(1)) then 
                    ywt = ( (y(j)+dy/2.d0) - ylim(1) ) / dy 
                else if (y(j)-dy/2.d0 .lt. ylim(2) .and. y(j)+dy/2.d0 .ge. ylim(2)) then 
                    ywt = ( ylim(2) - (y(j)-dy/2.d0) ) / dy 
                else if (y(j)-dy/2.d0 .ge. ylim(1) .and. y(j)+dy/2.d0 .le. ylim(2)) then 
                    ywt = 1.d0 
                else 
                    ywt = 0.d0 
                end if 
                
                if (xwt*ywt .gt. 0.d0) weight(i,j) = xwt*ywt 
            end do 
        end do 

        ! Check that weights are within range 0:1 
        if (minval(weight) .lt. 0.d0 .or. maxval(weight) .gt. 1.d0) then 
            write(*,*) "calc_grid_total:: error: weights less than 0 or greater than 1!"
            write(*,*) "weight: ", minval(weight), maxval(weight)
            stop 
        end if 

!         do i = 1, size(x)
!             write(*,"(5f10.2)") x(i), weight(i,1:4)
!         end do 

!         do j = 1, size(y)
!             write(*,"(5f10.2)") y(j), weight(1:4,j)
!         end do 
        
        ! Calculate grid total 
        tot = sum(var*weight*dx*dy)

        return 

    end function calc_grid_total

    ! === FIELD MAPPING SUBROUTINES === 


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
        real(dp) :: shepard_exponent

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
        
        var2  = reshape(nint(var2_vec), [nx2,ny2])
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
        real(dp) :: shepard_exponent

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
        real(dp) :: shepard_exponent

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
        
        var2  = reshape(nint(var2_vec),[nx2,ny2])
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
        real(dp) :: shepard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec

        allocate(var2_vec(size(var2)))

        call map_field_points_points_double(map,name,dble(var1),var2_vec,mask2, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack)

        var2 = nint(var2_vec)

        return

    end subroutine map_field_points_points_integer

    subroutine map_field_grid_grid_double(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack,sigma)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:,:), intent(IN)  :: var1
        real(dp), dimension(:,:), intent(INOUT) :: var2
        integer, dimension(:,:),  intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN),  optional :: mask_pack 
        real(dp), optional :: sigma 

        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shepard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable   :: mask_pack_vec 
        integer :: nx2, ny2, npts2, npts1 
        character(len=24) :: method_local 
        real(dp), dimension(:,:), allocatable :: var2tmp 

        method_local = trim(method)
        if (method .eq. "nng") method_local = "nn" 
        
        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 
        npts1 = size(var1,1)*size(var1,2)

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,reshape(var1,[npts1]),var2_vec,mask2_vec, &
                                     method_local,radius,fill,border,missing_value, &
                                     mask_pack_vec)
    
        var2  = reshape(var2_vec, [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        if (method .eq. "nng") then 
            if (.not. present(sigma)) then 
                write(*,*) "map_field:: error: method 'nng' requires &
                           &that the optional sigma argument be specified."
                stop 
            end if 

            if (present(fill)) then 
                if (fill) call fill_nearest(var2,missing_value=missing_value)
            end if 

            allocate(var2tmp(nx2,ny2))
            var2tmp = var2 
            call filter_gaussian(input=var2tmp,output=var2,sigma=sigma,dx=map%G%dx,&
                        mask=reshape(mask_pack_vec,[nx2,ny2]) .and. var2tmp .ne. missing_value)
        
        end if 

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
        real(dp) :: shepard_exponent

        integer :: npts1 

        npts1 = size(var1,1)*size(var1,2)

        call map_field_points_points_double(map,name,reshape(var1,[npts1]),var2,mask2, &
                                     method,radius,fill,border,missing_value, &
                                     mask_pack)

        return

    end subroutine map_field_grid_points_double

    subroutine map_field_points_grid_double(map,name,var1,var2,mask2,method,radius,fill,border, &
                                            missing_value,mask_pack,sigma)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(dp), dimension(:), intent(IN)    :: var1
        real(dp), dimension(:,:), intent(INOUT) :: var2
        integer,  dimension(:,:), intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN),  optional :: mask_pack 
        real(dp), optional :: sigma 

        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp) :: shepard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable :: mask_pack_vec 
        integer :: nx2, ny2, npts2
        character(len=24) :: method_local 
        real(dp), dimension(:,:), allocatable :: var2tmp 

        method_local = trim(method)
        if (method .eq. "nng") method_local = "nn" 
        
        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,var1,var2_vec,mask2_vec, &
                                     method_local,radius,fill,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(var2_vec, [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        if (method .eq. "nng") then 
            if (.not. present(sigma)) then 
                write(*,*) "map_field:: error: method 'nng' requires &
                           &that the optional sigma argument be specified."
                stop 
            end if 

            if (present(fill)) then 
                if (fill) call fill_nearest(var2,missing_value=missing_value)
            end if 
            
            allocate(var2tmp(nx2,ny2))
            var2tmp = var2 
            call filter_gaussian(input=var2tmp,output=var2,sigma=sigma,dx=map%G%dx,&
                        mask=reshape(mask_pack_vec,[nx2,ny2]) .and. var2tmp .ne. missing_value)
        
        end if 

        return

    end subroutine map_field_points_grid_double

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
        real(dp) :: shepard_exponent
        real(dp) :: max_distance, missing_val 
        real(dp), dimension(:), allocatable   :: weight_neighb, v_neighb
        real(dp), dimension(:), allocatable   :: dist_neighb
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
        allocate(dist_neighb(map%nmax))

        ! Initialize mask to show which points have been mapped
        allocate(mask2_local(map%npts))
        mask2_local = 0 

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2 = missing_val 

        ! Loop over the new grid points and perform mapping
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
                if (trim(method) .eq. "nn") then
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

                else if (trim(method) .eq. "quadrant") then 
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
                dist_neighb   = max_distance

                ! Reinitialize temp neighbor and weight values so that they appear in order (1:ntot)
                if (ntot .gt. 0) then 
                    q = 0 
                    do k = 1, map%nmax 
                        if (v_neighb_tmp(k) .ne. missing_val) then
                            q = q+1
                            weight_neighb(q) = map%weight(i,k)
                            dist_neighb(q)   = map%dist(i,k)
                            v_neighb(q)      = var1(map%i(i,k))
                        end if 
                    end do 
                end if           


                if ( ntot .gt. 1) then 

                    ! Calculate the weighted average (using distance weighting)
                    var2(i)  = weighted_ave(v_neighb(1:ntot),weight_neighb(1:ntot))
!                   var2(i)  = weighted_ave_shepard(v_neighb(1:ntot),dist_neighb(1:ntot),shephard_exponent=2.d0)
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
                if ( fill_pts .and. (var2(i) .eq. missing_val)) then
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

        where (dabs(var2) .lt. 1d-12) var2 = 0.d0 

!         write(*,*) "Mapped field: "//trim(name)

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



    ! === HELPER SUBROUTINES === 



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





    ! EXTRA - UNUSED (BUGGY??)
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
        real(dp) :: shepard_exponent
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

end module coordinates_mapping 