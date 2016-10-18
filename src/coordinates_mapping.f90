module coordinates_mapping 
    
    use coordinates 

    use oblimap_projection_module
    use planet 
    use polygons 
    use ncio
    use index

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
        integer,  allocatable :: nn(:)                        ! Only used for combining multiple pt_wts_class objects
        integer,  allocatable :: i(:), quadrant(:), border(:) ! Length of pts1/grid1 neighbors
        real(sp), allocatable :: x(:), y(:), dist(:), weight(:), area(:)
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
        integer  :: map_i(map_nmax), map_quadrant(map_nmax), map_border(map_nmax)
        real(dp) :: map_x(map_nmax), map_y(map_nmax), map_dist(map_nmax)
        
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
            map_quadrant = 0 
            map_border   = 0 
            map_x        = 0.0_dp 
            map_y        = 0.0_dp 
            map_dist     = ERR_DIST 
            
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
                map%map(i)%i        = 1 
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
        ! point given a vector of its neighbors (x,y) with corresponding 
        ! resolutions (dx,dy)
        ! Works only for regular (square) grid points

        real(sp), intent(IN) :: x(:), y(:), dx(:), dy(:) 
        real(sp), intent(IN) :: xout, yout, dxout, dyout 
        logical,  intent(IN), optional :: latlon  
        real(sp) :: area(size(x,1))

        ! Local variables
        logical  :: is_latlon 
        type(polygon) :: pol  
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

        type(pt_wts_class) :: mp_vec 
        integer :: n_vec 

        ! Pack the neighborhood map into vector format for writing 
        write(*,*) "Packing mp_vec..."
        call pack_neighbors(map%map,mp_vec)
        n_vec = size(mp_vec%i)

        ! Determine the filename automatically
        fnm = map_filename(map,fldr)

        ! Create the netcdf file and the dimension variables
        call nc_create(fnm)
        call nc_write_attr(fnm,"title","Mapping "//trim(map%name1)//" => "//trim(map%name2))

        ! Write generic dimensions
        call nc_write_dim(fnm,"point",       x=1,nx=map%npts,units="n")
        call nc_write_dim(fnm,"parameter",   x=1,units="")
        call nc_write_dim(fnm,"planetpar",   x=1,nx=3,units="")
        call nc_write_dim(fnm,"projpar",     x=1,nx=5,units="")
        call nc_write_dim(fnm,"neighbor_vec",x=1,nx=n_vec,units="n")
        
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

            ! Write grid specific parameters
            call nc_write(fnm,"nx",map%G%nx,dim1="parameter")
            call nc_write(fnm,"ny",map%G%ny,dim1="parameter")

        else
            ! Write variables in a vector format
            call nc_write(fnm,"x",       map%x,       dim1="point")
            call nc_write(fnm,"y",       map%y,       dim1="point")
            call nc_write(fnm,"lon",     map%lon,     dim1="point")
            call nc_write(fnm,"lat",     map%lat,     dim1="point")

        end if 

        write(*,*) "Writing mp_vec..."
        write(*,*) "size(mp_vec%nn): ", size(mp_vec%nn), minval(mp_vec%nn), maxval(mp_vec%nn)
        write(*,*) "size(mp_vec%i):  ", size(mp_vec%i),  minval(mp_vec%i),  maxval(mp_vec%i)
        
        ! Write neighborhood information 
        call nc_write(fnm,"mp_vec_nn",      mp_vec%nn,      dim1="point")
        call nc_write(fnm,"mp_vec_i",       mp_vec%i,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_x",       mp_vec%x,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_y",       mp_vec%y,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_dist",    mp_vec%dist,    dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_weight",  mp_vec%weight,  dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_quadrant",mp_vec%quadrant,dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_border",  mp_vec%border,  dim1="neighbor_vec")


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

        type(pt_wts_class) :: mp_vec 
        integer :: n_vec 

        ! Determine filename automatically 
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

        ! Allocate map neighborhood fields 
        ! Allocate the nn vector to the length of points being mapped to 
        allocate(mp_vec%nn(map%npts))
        
        n_vec = nc_size(fnm,"mp_vec_i")

        ! Allocate remaining vectors to correct length
        allocate(mp_vec%i(n_vec))
        allocate(mp_vec%quadrant(n_vec))
        allocate(mp_vec%border(n_vec))
        allocate(mp_vec%x(n_vec))
        allocate(mp_vec%y(n_vec))
        allocate(mp_vec%dist(n_vec))
        allocate(mp_vec%weight(n_vec))
        allocate(mp_vec%area(n_vec))
        
        ! Load neighborhood vectors 
        call nc_read(fnm,"mp_vec_i",        mp_vec%i)
        call nc_read(fnm,"mp_vec_quadrant", mp_vec%quadrant)
        call nc_read(fnm,"mp_vec_border",   mp_vec%border)
        call nc_read(fnm,"mp_vec_x",        mp_vec%x)
        call nc_read(fnm,"mp_vec_y",        mp_vec%y)
        call nc_read(fnm,"mp_vec_dist",     mp_vec%dist)
        call nc_read(fnm,"mp_vec_weight",   mp_vec%weight)
        call nc_read(fnm,"mp_vec_area",     mp_vec%area)
        
        ! Unpack neighbor vectors 
        call unpack_neighbors(mp_vec,map%map)


!         ! Allocate map fields 
!         allocate(map%x(map%npts),map%y(map%npts))
!         allocate(map%lon(map%npts),map%lat(map%npts))
!         allocate(map%i(map%npts,map%nmax),map%dist(map%npts,map%nmax))
!         allocate(map%weight(map%npts,map%nmax),map%quadrant(map%npts,map%nmax))
!         allocate(map%border(map%npts,map%nmax))

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

        else
            ! Read variables in a vector format

            call nc_read(fnm,       "x",map%x)
            call nc_read(fnm,       "y",map%y)
            call nc_read(fnm,     "lon",map%lon)
            call nc_read(fnm,     "lat",map%lat)

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
        ! Methods include "radius", "nn" (nearest neighbor), "quadrant"
        ! See `map_field_conservative_map` for 1st order conservative mapping
        
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
        integer,  dimension(:), allocatable   :: mask2_local
        integer :: i, k, q, j, ntot, check, n1  
        logical :: found 

        type(pt_wts_class)    :: map_now 
        real(dp), allocatable :: map_now_var(:) 
        integer,  allocatable :: ii(:) 

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

        ! Initialize mask to show which points have been mapped
        allocate(mask2_local(map%npts))
        mask2_local = 0 

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2 = missing_val 

        ! Loop over the new grid points and perform mapping
        do i = 1, map%npts 

            if (maskp(i)) then ! Only perform calculations for packing mask points

                ! Get current map and size of neighborhood
                map_now = map%map(i) 
                n1      = size(map_now%i)

                ! Get current variable values 
                map_now_var = var1(map_now%i)

                ! Eliminate neighbors outside of distance limit
                where(map_now%dist .gt. max_distance) map_now_var = missing_val

                ! Skip remaining calculations if no neighbors are found
                check = count(.not. map_now_var .eq. missing_val) 
                if (check .gt. 0) then 

                    ! If method == nn (nearest neighbor), limit neighbors to 1
                    if (trim(method) .eq. "nn") then

                        found = .FALSE. 
                        do k = 1, n1 
                            if (map_now_var(k) .ne. missing_val) then 
                                if (found) then 
                                    map_now_var(k) = missing_val 
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
                            do k = 1, n1 
                                if (map_now_var(k) .ne. missing_val .and. map_now%quadrant(k) .eq. q) then 
                                    if (found) then 
                                        map_now_var(k) = missing_val 
                                    else
                                        found = .TRUE. 
                                    end if 
                                end if 
                            end do 
                        end do 

                    end if 

                    ! Check number of neighbors available for calculations
                    call which(map_now_var .ne. missing_val,ii)
                    ntot = size(ii,1)

                    ! Check if a large fraction of neighbors are border points
                    ! (if so, do not interpolate here)
                    if ( (.not. fill_border) .and. ntot .gt. 0) then 
                        if ( sum(map_now%border(ii))/dble(ntot) .gt. 0.25_dp ) ntot = 0
                    end if 

                    ! Reset weights according to missing values 
                    where (map_now_var .eq. missing_val) 
                        map_now%weight = 0.0_dp 
                        map_now%dist   = ERR_DIST
                        map_now%area   = 0.0_dp 
                    end where 

                    if ( ntot .gt. 1) then 

                        ! Calculate the weighted average (using distance weighting)
                        var2(i)        = weighted_ave(map_now_var(ii),dble(map_now%weight(ii)))
    !                   var2(i)        = weighted_ave_shepard(map_now_var(ii),map_now%dist(ii),shephard_exponent=2.d0)
                        mask2_local(i) = 1

                    else if (ntot .eq. 1) then
                        call which(map_now%weight .ge. maxval(map_now%weight)*0.9999_dp,ii)
                        var2(i)        = map_now_var(ii(1))
                        mask2_local(i) = 1 

                    else
                        ! If no neighbors exist, field not mapped here.
                        mask2_local(i) = 0  

                    end if 

                    ! Fill missing points with nearest neighbor if desired
                    ! Note, will not necessarily fill ALL points, if 
                    ! no neighbor can be found without a missing value
                    if ( fill_pts .and. (var2(i) .eq. missing_val)) then
                        do k = 1, n1  
                            if (var1(map_now%i(k)) .ne. missing_val) then
                                var2(i)  = var1(map_now%i(k))
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

    subroutine pack_neighbors(map,map_vec)

        implicit none 

        type(pt_wts_class),     intent(IN)  :: map(:) 
        type(pt_wts_class), intent(OUT) :: map_vec 

        integer :: ntot, i, k, k1

        ! Deallocate map_vec vectors 
        if (allocated(map_vec%nn))       deallocate(map_vec%nn)
        if (allocated(map_vec%i))        deallocate(map_vec%i)
        if (allocated(map_vec%quadrant)) deallocate(map_vec%quadrant)
        if (allocated(map_vec%border))   deallocate(map_vec%border)
        if (allocated(map_vec%x))        deallocate(map_vec%x)
        if (allocated(map_vec%y))        deallocate(map_vec%y)
        if (allocated(map_vec%dist))     deallocate(map_vec%dist)
        if (allocated(map_vec%weight))   deallocate(map_vec%weight)
        if (allocated(map_vec%area))     deallocate(map_vec%area)
        
        ! Allocate the nn vector to the length of points being mapped to 
        allocate(map_vec%nn(size(map)))

        ! Determine length of storage vectors 
        ntot = 0 

        do i = 1, size(map)
            map_vec%nn(i) = size(map(i)%i,1)
            ntot = ntot + map_vec%nn(i)
        end do 

        ! Allocate remaining vectors
        call map_allocate_map(map_vec,ntot)
        
        k = 1 
        do i = 1, size(map)

            k1 = k + map_vec%nn(i) - 1

            map_vec%i(k:k1)        = map(i)%i 
            map_vec%quadrant(k:k1) = map(i)%quadrant 
            map_vec%border(k:k1)   = map(i)%border 
            map_vec%x(k:k1)        = map(i)%x 
            map_vec%y(k:k1)        = map(i)%y 
            map_vec%dist(k:k1)     = map(i)%dist 
            map_vec%weight(k:k1)   = map(i)%weight 
            map_vec%area(k:k1)     = map(i)%area 
            
            k = k1 + 1
        end do 

        return 

    end subroutine pack_neighbors 

    subroutine unpack_neighbors(map_vec,map)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: map_vec 
        type(pt_wts_class), intent(OUT)   :: map(:) 
        
        integer :: i, k, k1

        k = 1 
        do i = 1, size(map)

            ! Allocate current map object to right length  
            call map_allocate_map(map(i),map_vec%nn(i))

            k1 = k + map_vec%nn(i) - 1
            
            map(i)%i        = map_vec%i(k:k1)
            map(i)%quadrant = map_vec%quadrant(k:k1)
            map(i)%border   = map_vec%border(k:k1)
            map(i)%x        = map_vec%x(k:k1)
            map(i)%y        = map_vec%y(k:k1)
            map(i)%dist     = map_vec%dist(k:k1)
            map(i)%weight   = map_vec%weight(k:k1)
            map(i)%area     = map_vec%area(k:k1)
            
            k = k1 + 1
        end do 

        ! Deallocate map_vec vectors 
        if (allocated(map_vec%nn))       deallocate(map_vec%nn)
        if (allocated(map_vec%i))        deallocate(map_vec%i)
        if (allocated(map_vec%quadrant)) deallocate(map_vec%quadrant)
        if (allocated(map_vec%border))   deallocate(map_vec%border)
        if (allocated(map_vec%x))        deallocate(map_vec%x)
        if (allocated(map_vec%y))        deallocate(map_vec%y)
        if (allocated(map_vec%dist))     deallocate(map_vec%dist)
        if (allocated(map_vec%weight))   deallocate(map_vec%weight)
        if (allocated(map_vec%area))     deallocate(map_vec%area)
        
        return 

    end subroutine unpack_neighbors 

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




end module coordinates_mapping 