module coordinates_mapping 

    !$ use omp_lib
    use coord_constants 
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

    type pt_wts_class 
        integer :: n
        integer,  allocatable :: nn(:)                        ! Only used for combining multiple pt_wts_class objects
        integer,  allocatable :: i(:), quadrant(:), border(:) ! Length of pts1/grid1 neighbors
        real(sp), allocatable :: x(:), y(:), dist(:), weight(:), area(:) 

        ! For bilinear interpolation only:
        integer,  allocatable :: iquad(:) 
        real(sp), allocatable :: alpha1(:), alpha2(:), alpha3(:)  
    end type 

    
    type map_class
        character (len=128) :: name1, name2     ! Names of coordinate set1 (source) and set2 (target)
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

        ! Maximum number of neighbors to be stored for each point
        integer :: nmax 
        integer,  allocatable :: nn(:)

        ! Neighbor info (allocate to size npts)
        type(pt_wts_class), allocatable :: map(:)

    end type

    interface map_init 
        module procedure map_init_points_points, map_init_grid_grid 
        module procedure map_init_grid_points,   map_init_points_grid
    end interface

    interface map_field 
        module procedure map_field_grid_grid_double,    map_field_points_points_double
        module procedure map_field_grid_points_double,  map_field_points_grid_double
        module procedure map_field_grid_grid_float,     map_field_points_points_float
        module procedure map_field_grid_points_float,   map_field_points_grid_float
        module procedure map_field_grid_grid_integer,   map_field_points_points_integer
        module procedure map_field_grid_points_integer, map_field_points_grid_integer
    end interface

    interface compare_map  
        module procedure compare_map_map_grid 
    end interface

    private 
    public :: compare_map
    public :: map_class, map_init, map_field, map_print
    public :: pt_wts_class 

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
        ! Generate mapping weights between sets of points

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

        integer :: i, k  
        real(dp), allocatable :: xxnow(:), yynow(:), dxxnow(:), dyynow(:) 
        logical :: pts_is_latlon, pts_same_proj  

        pts_is_latlon = (.not. pts1%is_cartesian .and. .not. pts2%is_cartesian)
        pts_same_proj = same_projection(pts1%proj,pts2%proj)

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

            ! Allocate the neighborhood object
            allocate(map%map(map%npts))
            if (allocated(map%nn)) deallocate(map%nn)
            allocate(map%nn(map%npts))

            ! Calculate map weights (time consuming!)
            call map_calc_weights(map,pts1,pts2,2.0_dp,lat_lim,dist_max)

            ! If the grids are compatible, calculate area weighting too
            if ((pts1%is_projection .or. .not. pts1%is_cartesian) .and. pts2%is_projection ) then 
                ! Only valid for projected source and target grids right now 
                ! If pts1 and pts2 are not the same projection, conservative weights will 
                ! not be fully accurate, because of distorted dx/dy values.

                write(*,*) "Calculating conservative weights ..."
                do i = 1, pts2%npts
                    if (map%map(i)%dist(1) .lt. ERR_DIST) then 

                        if (allocated(xxnow))  deallocate(xxnow)
                        if (allocated(yynow))  deallocate(yynow)
                        if (allocated(dxxnow)) deallocate(dxxnow)
                        if (allocated(dyynow)) deallocate(dyynow)
                        k = size(map%map(i)%i,1)
                        allocate(xxnow(k),yynow(k),dxxnow(k),dyynow(k))

                        ! Determine current point on target grid 
                        ! xnow/ynow, dxnow/dynow are in [m]
                        if (pts_same_proj) then 
                            xxnow  = pts1%x(map%map(i)%i)*pts1%xy_conv
                            yynow  = pts1%y(map%map(i)%i)*pts1%xy_conv
                            dxxnow = pts1%dx(map%map(i)%i)*pts1%xy_conv
                            dyynow = pts1%dy(map%map(i)%i)*pts1%xy_conv
                        else 
                            do k = 1, size(map%map(i)%i)
                                call oblimap_projection(pts1%lon(map%map(i)%i(k)),pts1%lat(map%map(i)%i(k)), &
                                                        xxnow(k),yynow(k),pts2%proj)
                            end do 
                            dxxnow = pts1%dx(map%map(i)%i)*pts1%xy_conv
                            dyynow = pts1%dy(map%map(i)%i)*pts1%xy_conv
                        end if 

                        map%map(i)%area = calc_weights_interpconserv1( &
                                            x=real(xxnow),y=real(yynow), &
                                            dx=real(dxxnow),dy=real(dyynow), &
                                            xout=real(pts2%x(i)*pts2%xy_conv), &
                                            yout=real(pts2%y(i)*pts2%xy_conv), &
                                            dxout=real(pts2%dx(i)*pts2%xy_conv), &
                                            dyout=real(pts2%dy(i)*pts2%xy_conv))
                        map%map(i)%area = map%map(i)%area / (pts2%xy_conv*pts2%xy_conv)   ! Convert back to axis-units of pts2
                        
                    end if 
                    
                    ! Check progress
                    if (mod(i,100)==0) write(*,"(a,i10,a3,i12,a5,2g12.3)")  &
                                "  ",i, " / ",pts2%npts,"   : ", sum(map%map(i)%area), pts2%area(i) ! pts2%dx(i)*pts2%dy(i)
                end do 

            else 

                write(*,*) "Note: no conservative weights calculated since source and target grids &
                           &are on different projections..."

            end if 


            ! Write new map to file
            if (save_file) call map_write(map,mapfldr) 

        end if 

        ! Print map summary
        call map_print(map)

        return

    end subroutine map_init_internal

    subroutine map_allocate_map(mp,n,nblin)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mp 
        integer, intent(IN) :: n
        integer, intent(IN) :: nblin 
        mp%n = n  

        ! First deallocate if needed 
        call map_deallocate_map(mp) 

        allocate(mp%i(n),mp%quadrant(n), mp%border(n))
        allocate(mp%x(n),mp%y(n))
        allocate(mp%dist(n), mp%weight(n), mp%area(n))
        
        allocate(mp%iquad(nblin*4))
        allocate(mp%alpha1(nblin))
        allocate(mp%alpha2(nblin))
        allocate(mp%alpha3(nblin))
        
        return 

    end subroutine map_allocate_map

    subroutine map_deallocate_map(mp)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mp  

        ! First deallocate if needed 
        if (allocated(mp%i))        deallocate(mp%i) 
        if (allocated(mp%quadrant)) deallocate(mp%quadrant) 
        if (allocated(mp%border))   deallocate(mp%border) 
        if (allocated(mp%x))        deallocate(mp%x) 
        if (allocated(mp%y))        deallocate(mp%y)  
        if (allocated(mp%dist))     deallocate(mp%dist) 
        if (allocated(mp%weight))   deallocate(mp%weight) 
        if (allocated(mp%area))     deallocate(mp%area) 
        
        if (allocated(mp%iquad))    deallocate(mp%iquad)
        if (allocated(mp%alpha1))   deallocate(mp%alpha1)
        if (allocated(mp%alpha2))   deallocate(mp%alpha2)
        if (allocated(mp%alpha3))   deallocate(mp%alpha3)
        
        return 

    end subroutine map_deallocate_map

    subroutine map_calc_weights(map,pts1,pts2,shepard_exponent,lat_lim,dist_max)
        implicit none 

        type(map_class) :: map 
        type(points_class),   intent(IN)  :: pts1, pts2 
        real(dp),             intent(IN)  :: shepard_exponent
        real(dp),             intent(IN), optional :: lat_lim 
        real(dp),             intent(IN), optional :: dist_max   ! in units of pts2 axes!!

        real(dp), parameter :: DIST_ZERO_OFFSET = 1.0_dp  ! Change dist of zero to 1 m
        integer :: i, i1, kc, k, n 
        real(dp) :: x, y, lon, lat
        real(dp) :: dist, lat_limit, dist_maximum 
        real(dp) :: lon180, lon180_1

        integer, parameter :: map_nmax = 1000   ! No more than 1000 neighbors 
        type(pt_wts_class) :: mp_all 

        real :: start, finish

        ! Limit neighborhood to search 
        lat_limit = 5.0_dp 
        if (present(lat_lim)) lat_limit = lat_lim

        ! Distances are handled in units of target points
        dist_maximum = ERR_DIST 
        if (present(dist_max)) dist_maximum = dist_max

        write(*,"(a,i12,a,f6.2,a,g12.3)") "Total points to calculate=",pts2%npts, &
                    "  lat_lim=",lat_limit, "  dist_max=",dist_maximum

        ! Allocate the map to contain all neighbors within distance
        call map_allocate_map(mp_all,map_nmax,nblin=1)

        ! For each grid point in the new grid,
        ! Find points within a rough radius,
        ! calculate the distance to the current point
        ! and store in map
        call cpu_time(start)

        ! Note: for older versions of the ifort compiler (eg v12),
        ! the allocatable object mp_all appears not to be allocated inside
        ! of the processing loop below. The solution for such older compilers 
        ! is to use the following line:
        !
        ! !$omp parallel do private(i,x,y,lon,lat,i1,dist,kc,n) firstprivate(mp_all)
        !
        ! while for newer compilers, the following line works well:
        !
        ! !$omp parallel do private(i,x,y,lon,lat,mp_all,i1,dist,kc,n)
        !
        ! However, the solution for older compilers does not work for newer versions
        ! of ifort (eg ifort v18). Therefore, the default is to keep the
        ! newer solution activated, and revert to not using openmp with the
        ! older compilers (or activate the old solution as needed).
        
        !$omp parallel do private(i,x,y,lon,lat,lon180,lon180_1,mp_all,i1,dist,kc,n)
        do i = 1, pts2%npts
            
            !!$ print *,omp_get_thread_num(),'/',omp_get_num_threads()

            ! Get current xy and latlon coordinates
            x   = pts2%x(i)*pts2%xy_conv
            y   = pts2%y(i)*pts2%xy_conv
            lon = pts2%lon(i)
            lat = pts2%lat(i)

            ! Reset temp map values 
            mp_all%i        = ERR_IND
            mp_all%quadrant = 0 
            mp_all%border   = 0 
            mp_all%x        = 0.0_dp 
            mp_all%y        = 0.0_dp 
            mp_all%dist     = ERR_DIST 
            
            ! longitude between -180 and 180
            if (lon .ge. 180_dp) then
                lon180 = lon - 360_dp
            else if (lon .lt. -180_dp) then
                lon180 = lon + 360_dp
            else
                lon180 = lon
            end if

            ! Get distance in meters to current point on grid2
            ! for each point on grid1
            do i1 = 1, pts1%npts

                ! longitude between -180  and 180
                if (pts1%lon(i1) .ge. 180_dp) then
                    lon180_1 = pts1%lon(i1) - 360_dp
                else if (pts1%lon(i1) .lt. -180_dp) then
                    lon180_1 = pts1%lon(i1) + 360_dp
                else
                    lon180_1 = pts1%lon(i1)
                end if

                !if ( dabs(pts1%lat(i1)-lat).le.lat_limit .and.  (dabs(lon180_1-lon180).le.lat_limit .or. dabs(lon180_1-lon180).ge.(360_dp-lat_lim)) ) then
                if ( dabs(pts1%lat(i1)-lat).le.lat_limit) then

                    if (map%is_same_map .and. map%is_cartesian) then
                        ! Use cartesian values to determine distance
                        dist = cartesian_distance(x,y,pts1%x(i1)*pts1%xy_conv,pts1%y(i1)*pts1%xy_conv)

                    else
                        ! Use planetary (latlon) values
                        !dist = spherical_distance(map%planet%a,map%planet%f,lon,lat,pts1%lon(i1),pts1%lat(i1))
                        dist = planet_distance(map%planet%a,map%planet%f,lon,lat,pts1%lon(i1),pts1%lat(i1))

                    end if 

                    ! Convert distance to units of target grid 
                    ! Note: latlon grids use meters as distance unit (ie, pts2%xy_conv=1)
                    dist = dist / pts2%xy_conv 

                    ! Make sure no zero distances exist!
                    if (dist .lt. DIST_ZERO_OFFSET) dist = DIST_ZERO_OFFSET

                    ! Find location to insert current neighbor to
                    ! keep distances sorted in ascending order 
                    do kc = 1, map_nmax
                        if (dist .lt. mp_all%dist(kc)) exit
                    end do 

                    if (kc .le. map_nmax .and. dist .lt. dist_maximum) then 

                        ! Shift all other neighbors with larger distances right
                        ! to make room for new neighbor
                        if (kc .le. map_nmax-1) then 
                            mp_all%i(kc:map_nmax)        = cshift(mp_all%i(kc:map_nmax),-1)
                            mp_all%x(kc:map_nmax)        = cshift(mp_all%x(kc:map_nmax),-1)
                            mp_all%y(kc:map_nmax)        = cshift(mp_all%y(kc:map_nmax),-1)
                            mp_all%dist(kc:map_nmax)     = cshift(mp_all%dist(kc:map_nmax),-1)
                            mp_all%quadrant(kc:map_nmax) = cshift(mp_all%quadrant(kc:map_nmax),-1)
                            mp_all%border(kc:map_nmax)   = cshift(mp_all%border(kc:map_nmax),-1)
                        end if 

                        ! Store new neighbor in current index 
                        mp_all%i(kc)      = i1 
                        mp_all%x(kc)      = pts1%x(i1)
                        mp_all%y(kc)      = pts1%y(i1) 
                        mp_all%dist(kc)   = dist 
                        mp_all%border(kc) = pts1%border(i1)

                        ! Get quadrants of neighbors
                        if (map%is_same_map .and. map%is_cartesian) then
                            ! Use cartesian points to determine quadrants
                            mp_all%quadrant(kc) = quadrant_cartesian(x,y,pts1%x(i1)*pts1%xy_conv,pts1%y(i1)*pts1%xy_conv)
                        else
                            ! Use planetary (latlon) points
                            mp_all%quadrant(kc) = quadrant_latlon(lon,lat,pts1%lon(i1),pts1%lat(i1))
                        end if 

                        ! if (lat.lt.-87.5) then
                        !     write(*,*)
                        !     write(*,*) lat,pts1%lat(i1)
                        !     write(*,*) lon,pts1%lon(i1)
                        !     write(*,*) dist,mp_all%quadrant(kc)
                        ! endif

                    end if 
                end if 

            end do

            ! Store in neighborhood map now 
            n=count(mp_all%i .ne. ERR_IND)
            if (n .gt. 0) then 
                n = min(n,map%nmax)   ! Limit neighbors to max_neighbors specified if needed
                map%nn(i) = n         ! Store neighbor count
                call map_allocate_map(map%map(i),n=n,nblin=1)
                map%map(i)%i        = mp_all%i(1:n)
                map%map(i)%x        = mp_all%x(1:n)
                map%map(i)%y        = mp_all%y(1:n)
                map%map(i)%dist     = mp_all%dist(1:n)
                map%map(i)%weight   = 1.0_dp / (map%map(i)%dist**shepard_exponent)
                map%map(i)%quadrant = mp_all%quadrant(1:n)
                map%map(i)%border   = mp_all%border(1:n) 
                map%map(i)%area     = 0.0_dp 

            else 
                ! Store filler index
                n = 1 
                map%nn(i) = n 
                call map_allocate_map(map%map(i),n=n,nblin=1)
                map%map(i)%i        = 1 
                map%map(i)%x        = 0.0_dp 
                map%map(i)%y        = 0.0_dp 
                map%map(i)%dist     = ERR_DIST
                map%map(i)%weight   = 0.0_dp 
                map%map(i)%quadrant = 1
                map%map(i)%border   = 0
                map%map(i)%area     = 0.0_dp
 
            end if 

            ! Perform additional calculations to facilitate bilinear interpolation 
            call map_calc_weights_bilin(mp_all,x,y,lon,lat,pts1,map%planet, &
                    use_cartesian=(map%is_same_map .and. map%is_cartesian) )

            ! Store in current point's map 
            map%map(i)%iquad  = mp_all%iquad 
            map%map(i)%alpha1 = mp_all%alpha1 
            map%map(i)%alpha2 = mp_all%alpha2 
            map%map(i)%alpha3 = mp_all%alpha3 
            
            ! Output every 1000 rows to check progress
            if (mod(i,1000)==0) write(*,"(a,i10,a3,i12,a5,g12.3)")  &
                                    "  ",i, " / ",pts2%npts,"   : ",map%map(i)%dist(1)
        end do
        !$omp end parallel do

        call cpu_time(finish)
        write(*,"(a,a,f7.2)") "map_calc_weights:: "//trim(map%name1)//" => "//trim(map%name2)//": ", &
                              "Calculation time (min.) =", (finish-start)/60.0_dp

        return

    end subroutine map_calc_weights

    subroutine map_calc_weights_bilin(mnow,x,y,lon,lat,pts1,planet,use_cartesian)
        ! Pre-calculate the weighting and indices to make
        ! bilinear interpolation online faster 
        ! Note: assumes all variable data will be available 
        ! from the nearest four neighbors from four quadrants, ie, no missing values
        ! otherwise, the interpolated value will be missing too. 
        ! Note: also assumed that input points are on a regular grid, this should be
        ! managed by the user - otherwise interpolation will be incorrect. 

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mnow 
        real(dp),              intent(IN) :: x, y, lon, lat     ! Current location on target grid
        type(points_class),    intent(IN) :: pts1               ! Source points
        type(planet_class),    intent(IN) :: planet             ! Planet parameters of map
        logical,               intent(IN) :: use_cartesian 

        ! Local variables 
        integer :: q, k, n1, i1, i2, i3, i4, ntot  
        real(dp) :: xy_conv    
        real(dp) :: dx1, dx1_tot 
        real(dp) :: dx2, dx2_tot 
        real(dp) :: dy1, dy1_tot 
        real(dp) :: dy2, dy2_tot 
        real(dp) :: dx, dx_tot, dy, dy_tot 
        real(dp) :: lat_mid_1, lat_mid_0   
        real(dp) :: y_mid_1, y_mid_0 

        xy_conv = pts1%xy_conv 

        ! Initially set all 4 iquad indices to missing, 
        ! and bilin weights to missing 
        mnow%iquad  = ERR_IND
        mnow%alpha1 = 0.0_dp
        mnow%alpha2 = 0.0_dp
        mnow%alpha3 = 0.0_dp
        
        ! Get size of neighborhood 
        n1=count(mnow%i .ne. ERR_IND)

        if (n1 .ge. 4) then 
            ! At least four neighbors are needed for these calculations 

            ! Store the indices of four quadrants
            do q = 1, 4
                do k = 1, n1
                    if (mnow%quadrant(k) .eq. q) then
                        ! Store first index of this quadrant that we find,
                        ! which should also be the closest since neighbors are sorted,
                        ! then skip the rest 
                        mnow%iquad(q) = mnow%i(k)
                        exit 
                    end if
                end do
            end do

            ! Check number of neighbors available for calculations
            ntot = count(mnow%iquad .ne. ERR_IND)
            !write(*,*) ntot 

            if (ntot.eq.4) then
                ! All four quadrant neighbors exist 
                ! calculate bilin weights 

                ! Get indices of relevant neighbors
                i1 = mnow%iquad(1)    ! Quadrant 1 (above-right of point)
                i2 = mnow%iquad(2)    ! Quadrant 2 (above-left  of point)
                i3 = mnow%iquad(3)    ! Quadrant 3 (below-left  of point)
                i4 = mnow%iquad(4)    ! Quadrant 4 (below-right of point) 

                if (use_cartesian) then
                    ! Use cartesian values to determine distance
                    ! Note: extra interpolation to y_mid_ points is not necessary
                    ! on a Cartesian grid, but the form is kept to be consistent
                    ! with the latlon grid algorithm 

                    y_mid_1 = 0.5_dp*(pts1%y(i2)+pts1%y(i1))
                    y_mid_0 = 0.5_dp*(pts1%y(i3)+pts1%y(i4))
                    
                    dx1     = cartesian_distance(                 x,   y_mid_1*xy_conv, &
                                                 pts1%x(i2)*xy_conv,pts1%y(i2)*xy_conv)
                    dx1_tot = cartesian_distance(pts1%x(i1)*xy_conv,pts1%y(i1)*xy_conv, &
                                                 pts1%x(i2)*xy_conv,pts1%y(i2)*xy_conv)
                    
                    dx2     = cartesian_distance(                 x,   y_mid_0*xy_conv, &
                                                 pts1%x(i3)*xy_conv,pts1%y(i3)*xy_conv)
                    dx2_tot = cartesian_distance(pts1%x(i4)*xy_conv,pts1%y(i4)*xy_conv, &
                                                 pts1%x(i3)*xy_conv,pts1%y(i3)*xy_conv)
                    
                    dy1     = cartesian_distance(                 x,   y, &
                                                                  x,   y_mid_0*xy_conv)
                    dy1_tot = cartesian_distance(                 x,   y_mid_1*xy_conv, &
                                                                  x,   y_mid_0*xy_conv)
                    
                else
                    ! Use planetary (latlon) values
                    ! Note: it is not strictly correct to average the latitudes
                    ! evenly - this should be done proportionally to the distance,
                    ! but to the first order it's ok. 

                    lat_mid_1 = 0.5_dp*(pts1%lat(i2)+pts1%lat(i1))
                    lat_mid_0 = 0.5_dp*(pts1%lat(i3)+pts1%lat(i4))
                    
                    dx1     = planet_distance(planet%a,planet%f,         lon,  lat_mid_1, &
                                                                pts1%lon(i2),pts1%lat(i2))
                    dx1_tot = planet_distance(planet%a,planet%f,pts1%lon(i1),pts1%lat(i1), &
                                                                pts1%lon(i2),pts1%lat(i2))
                    
                    dx2     = planet_distance(planet%a,planet%f,         lon,  lat_mid_0, &
                                                                pts1%lon(i3),pts1%lat(i3))
                    dx2_tot = planet_distance(planet%a,planet%f,pts1%lon(i4),pts1%lat(i4), &
                                                                pts1%lon(i3),pts1%lat(i3))
                    
                    dy1     = planet_distance(planet%a,planet%f,         lon,        lat, &
                                                                         lon,   lat_mid_0)
                    dy1_tot = planet_distance(planet%a,planet%f,         lon,   lat_mid_1, &
                                                                         lon,   lat_mid_0)
                    
                end if 

                ! if (lat.lt.-87.5) then
                !     write(*,*)
                !     write(*,*) lat,lon
                !     write(*,*) pts1%lat(i1),pts1%lon(i1)
                !     write(*,*) pts1%lat(i2),pts1%lon(i2)
                !     write(*,*) pts1%lat(i3),pts1%lon(i3)
                !     write(*,*) pts1%lat(i4),pts1%lon(i4)
                !     write(*,*) dx1,dx1_tot
                !     write(*,*) dx2,dx2_tot
                !     write(*,*) dy1,dy1_tot
                ! endif

                if (dx1_tot .gt. 0.0) mnow%alpha1 = dx1 / dx1_tot 
                if (dx2_tot .gt. 0.0) mnow%alpha2 = dx2 / dx2_tot
                if (dy1_tot .gt. 0.0) mnow%alpha3 = dy1 / dy1_tot

            end if

        end if 

        return 

    end subroutine map_calc_weights_bilin

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
                
                ! Define number of subgrid points to use for area calculation
                nx = max(5,int(dx(now)/dxout))
                ny = max(5,int(dy(now)/dyout))
                npts = nx*ny 

!                 write(*,*) now, nx, ny, npts 
                
                npts_in   = 0

                ! Determine how many subgrid points are inside of the target point boundaries
                do j = 1, ny 
                    do i = 1, nx 
                        x1 = (x(now)-dx(now)/2.d0) + (dx(now))*dble(i-1)/dble(nx) + 0.5d0*1.d0/dble(nx)
                        y1 = (y(now)-dy(now)/2.d0) + (dy(now))*dble(j-1)/dble(ny) + 0.5d0*1.d0/dble(ny) 
                        if (point_in_polygon(real(x1),real(y1),pol)) npts_in = npts_in+1
                    end do
                end do 
                
                ! Get area of target point that current neighbor point contains
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
        character(len=10) :: char10

        write(char10,"(i10)") map%nmax
        map_filename = trim(fldr)//"/map_"//trim(map%name1)//"_"//trim(map%name2)// &
                            "_"//trim(adjustl(char10))//".nc"

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
        call nc_write_dim(fnm,"point_4",     x=1,nx=map%npts*4,units="n")
        
        ! Write grid/vector specific dimensions and variables
        if (map%is_grid) then
            ! Write variables in a gridded format

            select case(trim(map%mtype))

                case("latitude_longitude","latlon","gaussian")
                    ! latitude_longitude 

                    dim1 = "lon"
                    dim2 = "lat"
                    call nc_write_dim(fnm,dim1,x=map%G%x)
                    call nc_write_dim(fnm,dim2,x=map%G%y)
                
                case DEFAULT 
                    ! Projection and/or cartesian 

                    dim1 = "xc"
                    dim2 = "yc" 
                    call nc_write_dim(fnm,dim1,x=map%G%x,units=trim(map%units))
                    call nc_write_dim(fnm,dim2,x=map%G%y,units=trim(map%units))
            end select  

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

        ! Write neighborhood information 
        call nc_write(fnm,"mp_vec_nn",      mp_vec%nn,      dim1="point")
        call nc_write(fnm,"mp_vec_i",       mp_vec%i,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_quadrant",mp_vec%quadrant,dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_border",  mp_vec%border,  dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_x",       mp_vec%x,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_y",       mp_vec%y,       dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_dist",    mp_vec%dist,    dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_weight",  mp_vec%weight,  dim1="neighbor_vec")
        call nc_write(fnm,"mp_vec_area",    mp_vec%area,    dim1="neighbor_vec")
        
        call nc_write(fnm,"mp_vec_iquad",   mp_vec%iquad,   dim1="point_4")
        call nc_write(fnm,"mp_vec_alpha1",  mp_vec%alpha1,  dim1="point")
        call nc_write(fnm,"mp_vec_alpha2",  mp_vec%alpha2,  dim1="point")
        call nc_write(fnm,"mp_vec_alpha3",  mp_vec%alpha3,  dim1="point")
        
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
        
        ! Determine mp_vector length
        n_vec = nc_size(fnm,"neighbor_vec")

        ! Allocate remaining vectors to correct length
        call map_allocate_map(mp_vec,n_vec,nblin=map%npts)

        ! Load neighborhood vectors 
        call nc_read(fnm,"mp_vec_nn",       mp_vec%nn)
        call nc_read(fnm,"mp_vec_i",        mp_vec%i)
        call nc_read(fnm,"mp_vec_quadrant", mp_vec%quadrant)
        call nc_read(fnm,"mp_vec_border",   mp_vec%border)
        call nc_read(fnm,"mp_vec_x",        mp_vec%x)
        call nc_read(fnm,"mp_vec_y",        mp_vec%y)
        call nc_read(fnm,"mp_vec_dist",     mp_vec%dist)
        call nc_read(fnm,"mp_vec_weight",   mp_vec%weight)
        call nc_read(fnm,"mp_vec_area",     mp_vec%area)
        
        call nc_read(fnm,"mp_vec_iquad",    mp_vec%iquad)
        call nc_read(fnm,"mp_vec_alpha1",   mp_vec%alpha1)
        call nc_read(fnm,"mp_vec_alpha2",   mp_vec%alpha2)
        call nc_read(fnm,"mp_vec_alpha3",   mp_vec%alpha3)
        
        ! Allocate map%map 
        if (allocated(map%map)) deallocate(map%map)
        allocate(map%map(map%npts))
            
        ! Store neighbor counts
        if (allocated(map%nn)) deallocate(map%nn)
        allocate(map%nn(map%npts))
        map%nn = mp_vec%nn 

        ! Unpack neighbor vectors 
        call unpack_neighbors(mp_vec,map%map)

        ! Read grid/vector specific dimensions and variables
        if (map%is_grid) then
            ! Read variables in a gridded format

            call nc_read(fnm,"nx",map%G%nx)
            call nc_read(fnm,"ny",map%G%ny)

            ! Allocate grid fields
            allocate(map%G%x(map%G%nx))
            allocate(map%G%y(map%G%ny))

            select case(trim(map%mtype))

                case("latitude_longitude","latlon","gaussian")
                    ! latitude_longitude 

                    call nc_read(fnm,"lon",map%G%x)
                    call nc_read(fnm,"lat",map%G%y)

                case DEFAULT 
                    ! Projection and/or cartesian 

                    call nc_read(fnm,"xc",map%G%x)
                    call nc_read(fnm,"yc",map%G%y)

            end select 

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
        integer :: i, n_vec 

        n_vec = 0
        do i = 1, map%npts 
            n_vec = n_vec + size(map%map(i)%i)
        end do 

        write(*,*) "== Mapping summary ========================="
        write(*,*) trim(map%name1)," => ",trim(map%name2)
        if (map%is_grid) then
!             write(*,*) "  Grid axis information"
            write(*,"(a17,i11,i11)") "nx,ny = ", map%G%nx, map%G%ny 
        end if 
        write(*,"(a17,i11,i11)")     "npts, nvec = ", map%npts, n_vec 
        

        write(*,"(a17,g11.4,a3)")   "Size in memory ~ ", &
            ( 4.d0*(map%npts*8.d0) &
            + 2.d0*(n_vec*4.d0)    &
            + 5.d0*(n_vec*8.d0)    &
            + 7.d0*(map%npts*4.d0) ) *7.6294d-6, " Mb"      
        ! 4 real arrays (8 bytes per value): x, y, lon, lat
        ! 3 integer array (4 bytes per value): i, quadrant, border
        ! 5 real arrays (8 bytes per value): x, y, dist, weight, area
        ! 7 real arays (4 bytes per value): iquad, alpha1, alpha2, alpha3
        write(*,*)
        write(*,*) 

        return
    end subroutine map_print

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
                                     method,radius,border,missing_value, &
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
                                     method,radius,border,missing_value, &
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
                                     method,radius,border,missing_value, &
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
                                     method,radius,border,missing_value, &
                                     mask_pack)

        var2 = nint(var2_vec)

        return

    end subroutine map_field_points_points_integer

    subroutine map_field_grid_grid_float(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack,sigma)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(sp), dimension(:,:), intent(IN)   :: var1
        real(sp), dimension(:,:), intent(OUT)  :: var2
        integer,  dimension(:,:), intent(OUT), optional :: mask2
        logical,  dimension(:,:), intent(IN), optional  :: mask_pack 
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border
        real(dp), optional :: sigma 
        real(dp) :: shepard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec
        integer,  dimension(:), allocatable   :: mask2_vec
        logical,  dimension(:), allocatable   :: mask_pack_vec 
        integer :: nx2, ny2, npts2, npts1 
        character(len=24) :: method_local
        real(dp), dimension(:,:), allocatable :: var2dp 
        real(dp) :: missing_val 

        method_local = trim(method)
        if (method .eq. "nng") method_local = "nn"

        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 
        npts1 = size(var1,1)*size(var1,2)

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,reshape(dble(var1),[npts1]),var2_vec,mask2_vec, &
                                     method_local,radius,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(real(var2_vec,sp), [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        if (method .eq. "nng") then
            if (.not. present(sigma)) then
                write(*,*) "map_field:: error: method 'nng' requires &
                           &that the optional sigma argument be specified."
                stop
            end if
            
            allocate(var2dp(nx2,ny2))
            var2dp = dble(var2)

            if (present(fill)) then
                if (fill) call fill_nearest(var2dp,missing_value=missing_val)
            end if

            call filter_gaussian(var=var2dp,sigma=sigma,dx=map%G%dx,&
                        mask=reshape(mask_pack_vec,[nx2,ny2]) .and. var2 .ne. missing_val)

            var2 = real(var2dp,sp)

        end if

        return

    end subroutine map_field_grid_grid_float

    subroutine map_field_grid_points_float(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(sp), dimension(:,:), intent(IN)   :: var1
        real(sp), dimension(:), intent(OUT)    :: var2
        integer,  dimension(:), intent(OUT), optional :: mask2
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
                                     method,radius,border,missing_value, &
                                     mask_pack)

        var2 = real(var2_vec,sp)

        return

    end subroutine map_field_grid_points_float

    subroutine map_field_points_grid_float(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)           :: map 
        real(sp), dimension(:), intent(IN)     :: var1
        real(sp), dimension(:,:), intent(INOUT)  :: var2
        integer,  dimension(:,:),  intent(OUT), optional :: mask2
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
                                     method,radius,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(real(var2_vec,sp),[nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        return

    end subroutine map_field_points_grid_float

    subroutine map_field_points_points_float(map,name,var1,var2,mask2,method,radius,fill,border,missing_value,mask_pack)

        implicit none 

        type(map_class), intent(IN)          :: map 
        real(sp), dimension(:), intent(IN)    :: var1
        real(sp), dimension(:), intent(INOUT) :: var2
        integer, dimension(:),  intent(OUT), optional :: mask2
        logical,  dimension(:), intent(IN),  optional :: mask_pack 
        
        character(len=*) :: name, method
        real(dp), optional :: radius, missing_value 
        logical,  optional :: fill, border  
        real(dp) :: shepard_exponent

        real(dp), dimension(:), allocatable   :: var2_vec

        allocate(var2_vec(size(var2)))

        call map_field_points_points_double(map,name,dble(var1),var2_vec,mask2, &
                                     method,radius,border,missing_value, &
                                     mask_pack)

        var2 = real(var2_vec,sp)

        return

    end subroutine map_field_points_points_float

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
        real(dp) :: missing_val 
        integer :: nx2, ny2, npts2, npts1 
        character(len=24) :: method_local 

        method_local = trim(method)
        if (method .eq. "nng") method_local = "nn" 
        
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value 

        nx2   = size(var2,1)
        ny2   = size(var2,2)
        npts2  = nx2*ny2 
        npts1 = size(var1,1)*size(var1,2)

        allocate(var2_vec(npts2),mask2_vec(npts2),mask_pack_vec(npts2))
        var2_vec = reshape(var2, [npts2])
        mask_pack_vec = .TRUE. 
        if (present(mask_pack)) mask_pack_vec = reshape(mask_pack,[npts2])

        call map_field_points_points_double(map,name,reshape(var1,[npts1]),var2_vec,mask2_vec, &
                                     method_local,radius,border,missing_val, &
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
                if (fill) call fill_nearest(var2,missing_value=missing_val)
            end if 
             
            call filter_gaussian(var=var2,sigma=sigma,dx=map%G%dx,&
                        mask=reshape(mask_pack_vec,[nx2,ny2]) .and. var2 .ne. missing_val)
        
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
                                     method,radius,border,missing_value, &
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
                                     method_local,radius,border,missing_value, &
                                     mask_pack_vec)
        
        var2  = reshape(var2_vec, [nx2,ny2])
        if (present(mask2)) mask2 = reshape(mask2_vec,[nx2,ny2])

        ! Fill missing values if requested
        if (present(fill)) then 
            !if (fill) call fill_nearest(var2,missing_value=missing_value)               ! Old routine  < 2020.06.18
            if (fill) call fill_nearest(var2,missing_value,fill_value=missing_value, &   ! New routine >= 2020.06.18
                                                fill_dist=1d20,n=5,dx=map%G%dx)    
        end if 
        
        ! Apply additional Gaussian smoothing if method==nng
        if (method .eq. "nng") then 
            if (.not. present(sigma)) then 
                write(*,*) "map_field:: error: method 'nng' requires &
                           &that the optional sigma argument be specified."
                stop 
            end if 

            call filter_gaussian(var=var2,sigma=sigma,dx=map%G%dx,&
                        mask=reshape(mask_pack_vec,[nx2,ny2]) .and. var2 .ne. missing_value)
        
        end if 

        return

    end subroutine map_field_points_grid_double

    subroutine map_field_points_points_double(map,name,var1,var2,mask2,method,radius,border,missing_value,mask_pack)
        ! Methods include "radius", "nn" (nearest neighbor), "quadrant"
        ! See `map_field_conservative_map1` for 1st order conservative mapping
        ! This is `map_field_internal` basically and is 
        ! a wrapper for calling specific methods 

        implicit none 

        type(map_class),  intent(IN)    :: map
        character(len=*), intent(IN)    :: name 
        real(dp),         intent(IN)    :: var1(:)
        real(dp),         intent(INOUT) :: var2(:)
        integer,          intent(OUT), optional :: mask2(:)
        character(len=*), intent(IN) :: method

        real(dp), intent(IN),  optional :: radius
        logical,  intent(IN),  optional :: border 
        real(dp), intent(IN),  optional :: missing_value 
        logical,  intent(IN),  optional :: mask_pack(:) 

        ! Local variables 
        real(dp) :: max_distance, missing_val  
        logical, allocatable :: maskp(:)
        integer, allocatable :: mask2_local(:) 

        ! Set neighborhood radius to very large value (to include all neighbors)
        ! or to radius specified by user
        max_distance = 1E7_dp
        if (present(radius)) max_distance = radius 

        ! Set grid missing value by default or that that specified by user
        missing_val  = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value 

        ! By default, all var2 points are interpolated
        allocate(maskp(size(var2)))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = mask_pack 
        
        ! Initialize mask to show which points have been mapped
        allocate(mask2_local(map%npts))
        mask2_local = 0 
        
        ! Apply method of choice 
        select case(trim(method))

            case("nn","nearest")

                call map_field_nearest(map,maskp,var1,var2,mask2_local,max_distance,missing_val)

            case("bilinear")

                call map_field_bilinear(map,maskp,var1,var2,mask2_local,max_distance,missing_val)

            case("quadrant")

                call map_field_radius(map,maskp,var1,var2,mask2_local,max_distance,missing_val, &
                                      is_quadrant=.TRUE.,border=border)

            case("radius","shepard")

                call map_field_radius(map,maskp,var1,var2,mask2_local,max_distance,missing_val, &
                                      is_quadrant=.FALSE.,border=border)

            case DEFAULT 

                write(*,*) "map_field:: error: method not recognized: "//trim(method)
                stop 

        end select 

        ! Post-processing 

        ! Avoid underflow errors
        where (dabs(var2) .lt. 1d-12) var2 = 0.d0 

        ! If interpolation mask available, send to output
        if (present(mask2)) mask2 = mask2_local 

        return 

    end subroutine map_field_points_points_double

    subroutine map_field_nearest(map,maskp,var1,var2,mask2,max_distance,missing_val)

        implicit none 


        type(map_class), intent(IN)  :: map       ! Map information
        logical,         intent(IN)  :: maskp(:)  ! Mask for which points should be interpolated 
        real(dp),        intent(IN)  :: var1(:)   ! Source variable 
        real(dp),        intent(OUT) :: var2(:)   ! Target grid variable 
        integer,         intent(OUT) :: mask2(:)  ! Interpolation mask
        real(dp),        intent(IN)  :: max_distance  
        real(dp),        intent(IN)  :: missing_val 

        ! Local variables 
        integer :: i, k, nmax, n  
        real(dp), allocatable :: mp_var(:) 

        ! Allocate a storage array to hold neighbors
        ! (as large as the largest available neighborhood)
        nmax = maxval(map%nn)
        allocate(mp_var(nmax))

        ! Loop over the new grid points and perform mapping
        do i = 1, map%npts 

            if (maskp(i)) then ! Only perform calculations for packing mask points

                ! Get size of neighborhood 
                n = size(map%map(i)%i)
                
                ! Store neighborhood variable values 
                mp_var       = missing_val 
                mp_var(1:n)  = var1(map%map(i)%i)

                ! Find nearest available neighbor, save it if valid
                k = minloc(map%map(i)%dist,mask=mp_var(1:n).ne.missing_val,dim=1)
                if (k .gt. 0) then 
                    if (map%map(i)%dist(k) .lt. max_distance) then
                        var2(i)  = mp_var(k) 
                        mask2(i) = 1
                    end if 
                end if 

            end if 

        end do 
                        
        return 

    end subroutine map_field_nearest

    subroutine map_field_bilinear(map,maskp,var1,var2,mask2,max_distance,missing_val)

        implicit none 


        type(map_class), intent(IN)  :: map       ! Map information
        logical,         intent(IN)  :: maskp(:)  ! Mask for which points should be interpolated 
        real(dp),        intent(IN)  :: var1(:)   ! Source variable 
        real(dp),        intent(OUT) :: var2(:)   ! Target grid variable 
        integer,         intent(OUT) :: mask2(:)  ! Interpolation mask
        real(dp),        intent(IN)  :: max_distance  
        real(dp),        intent(IN)  :: missing_val 

        ! Local variables 
        integer  :: i, ntot, n, k
        real(dp) :: mp_var(4)   ! Only four neighbors needed

        ! Loop over the new grid points and perform mapping
        do i = 1, map%npts 

            if (maskp(i)) then ! Only perform calculations for packing mask points

                if (count(map%map(i)%iquad.ne.ERR_IND).eq.4) then  
                    ! All four neighbor locations are available, proceed... 

                    ! Store the four neighbor values based on map indices
                    mp_var = var1(map%map(i)%iquad)
                    
                    ! How many neighbors are actually available
                    ntot = count(mp_var.ne.missing_val)

                    if (ntot.eq.4) then 
                        ! All four variable values are available, proceed...

                        var2(i)  = calc_bilin_point(mp_var(1),mp_var(2), &
                                                    mp_var(3),mp_var(4), &
                                                    map%map(i)%alpha1(1), &
                                                    map%map(i)%alpha2(1), &
                                                    map%map(i)%alpha3(1))
                        
                        mask2(i) = 1

                    end if 

                else ! nearest neighbour

                    ! Get size of neighborhood
                    n = size(map%map(i)%i)

                    ! Store neighborhood variable values
                    mp_var       = missing_val
                    mp_var(1:n)  = var1(map%map(i)%i)

                    ! Find nearest available neighbor, save it if valid
                    k = minloc(map%map(i)%dist,mask=mp_var(1:n).ne.missing_val,dim=1)
                    if (map%map(i)%dist(k) .lt. max_distance) then
                        var2(i)  = mp_var(k)
                        mask2(i) = 1
                    end if

                end if 
                
            end if 

        end do 
                        
        return 

    end subroutine map_field_bilinear
    
    subroutine map_field_radius(map,maskp,var1,var2,mask2,max_distance,missing_val,is_quadrant,border)
        ! Map a field using distance weighting
        ! (optionally limit neighbors to four neighbors in different quadrants - is_quadrant=.TRUE.)

        implicit none 

        type(map_class), intent(IN)  :: map       ! Map information
        logical,         intent(IN)  :: maskp(:)  ! Mask for which points should be interpolated 
        real(dp),        intent(IN)  :: var1(:)   ! Source variable 
        real(dp),        intent(OUT) :: var2(:)   ! Target grid variable 
        integer,         intent(OUT) :: mask2(:)  ! Interpolation mask
        real(dp),        intent(IN)  :: max_distance  
        real(dp),        intent(IN)  :: missing_val 
        logical,         intent(IN)  :: is_quadrant 
        logical,         intent(IN), optional :: border 

        ! Local variables
        logical :: fill_border  
        integer :: i, k, q, nmax, n, ntot
        logical :: found     
        real(dp), allocatable :: mp_var(:) 
        type(pt_wts_class)    :: mp_now

        ! By default, border points will not be filled in 
        fill_border = .FALSE. 
        if (present(border)) fill_border = border 

        ! Allocate a storage array to hold neighbors
        ! (as large as the largest available neighborhood)
        nmax = maxval(map%nn)
        allocate(mp_var(nmax))
        call map_allocate_map(mp_now,nmax,nblin=1)

        ! Loop over the new grid points and perform mapping
        do i = 1, map%npts 

            if (maskp(i)) then ! Only perform calculations for packing mask points

                ! Distance weighted interpolation 

                ! Reset current map values 
                mp_var          = missing_val 
                mp_now%dist     = missing_val 
                mp_now%weight   = 0.0_dp 
                mp_now%quadrant = 1
                mp_now%border   = 0

                ! Store current map

                ! Get size of neighborhood 
                n = size(map%map(i)%i)
                
                ! Store neighborhood variable values 
                mp_var(1:n)          = var1(map%map(i)%i)
                mp_now%dist(1:n)     = map%map(i)%dist
                mp_now%weight(1:n)   = map%map(i)%weight
                mp_now%quadrant(1:n) = map%map(i)%quadrant
                mp_now%border(1:n)   = map%map(i)%border

                if (is_quadrant) then 
                    ! For quadrant method, limit the number of neighbors to 
                    ! 4 points in different quadrants from neighborhood
                    ! with data available  
                    do q = 1, 4
                        found = .FALSE. 
                        do k = 1, n 
                            if (mp_var(k) .ne. missing_val .and. mp_now%quadrant(k) .eq. q) then 
                                if (found) then 
                                    mp_var(k) = missing_val 
                                else
                                    found = .TRUE. 
                                end if 
                            end if 
                        end do 
                    end do 

                end if  

                ! Eliminate neighbors outside of distance limit
                where(mp_now%dist .gt. max_distance) mp_var = missing_val

                ! Set all missing values to maximum distance
                where(mp_var .eq. missing_val)  mp_now%dist = ERR_DIST

                ! Check number of neighbors available for calculations
                ntot = count(mp_var .ne. missing_val)

                ! Check if a large fraction of neighbors are border points
                ! (if so, do not interpolate here)
                if ( (.not. fill_border) .and. ntot .gt. 0) then 
                    if ( sum(mp_now%border,mask=mp_var.ne.missing_val)/dble(ntot) .gt. 0.25_dp ) &
                            ntot = 0
                end if 

                ! Apply appropriate interpolation calculation
                if ( ntot .gt. 1) then 

                    ! Calculate the weighted average (using distance weighting)
                    var2(i)        = weighted_ave(mp_var,real(mp_now%weight,dp),mask=mp_var .ne. missing_val)
!                     var2(i)        = weighted_ave_shepard(mp_var,mp_now%dist,shepard_exponent=2.d0, &
!                                                           mask=mp_var .ne. missing_val)
                    mask2(i) = 1

                else if (ntot .eq. 1) then
                    k = minloc(mp_now%dist,1)
                    var2(i)  = mp_var(k)
                    mask2(i) = 1 
                else
                    ! If no neighbors exist, field not mapped here.
                    mask2(i) = 0  

                end if 

            end if 

        end do 
                        
        return 

    end subroutine map_field_radius


    ! === HELPER SUBROUTINES === 

    function calc_bilin_point(z1,z2,z3,z4,alpha1,alpha2,alpha3) result(var)
        ! Helper function to calculate the bilinear interpolate
        ! given the four quadrant neighbor values and the weights alpha* 

        implicit none 

        real(dp), intent(IN) :: z1, z2, z3, z4
        real(sp), intent(IN) :: alpha1, alpha2, alpha3
        real(dp)             :: var 

        ! Local variables 
        real(dp) :: p0, p1 

        p1  = z2 + alpha1*(z1-z2)
        p0  = z3 + alpha2*(z4-z3)
        
        var = p0 + alpha3*(p1-p0)

        return 

    end function calc_bilin_point 

    subroutine pack_neighbors(map,map_vec)

        implicit none 

        type(pt_wts_class), intent(IN)  :: map(:) 
        type(pt_wts_class), intent(OUT) :: map_vec 

        integer :: ntot, npt, i, k, k1

        ! Determine how many points there are 
        npt = size(map)
        
        ! Allocate the nn vector to the length of points being mapped to 
        if (allocated(map_vec%nn)) deallocate(map_vec%nn)
        allocate(map_vec%nn(npt))

        ! Determine length of storage vectors 
        ntot = 0 

        do i = 1, size(map)
            map_vec%nn(i) = size(map(i)%i,1)
            ntot = ntot + map_vec%nn(i)
        end do 

        ! Allocate remaining vectors
        call map_allocate_map(map_vec,ntot,nblin=npt)
        
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

        ! Store bilinear indices and weights 
        k = 1 
        do i = 1, size(map)
            k1 = k + 4 - 1 
            map_vec%iquad(k:k1)    = map(i)%iquad(1:4) 
            map_vec%alpha1(i)      = map(i)%alpha1(1) 
            map_vec%alpha2(i)      = map(i)%alpha2(1) 
            map_vec%alpha3(i)      = map(i)%alpha3(1) 
            k = k1 + 1 
        end do 

        return 

    end subroutine pack_neighbors 

    subroutine unpack_neighbors(map_vec,map)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: map_vec 
        type(pt_wts_class), intent(OUT)   :: map(:) 
        
        integer :: i, k, k1, npt 

        k = 1 
        do i = 1, size(map)

            ! Allocate current map object to right length  
            call map_allocate_map(map(i),map_vec%nn(i),nblin=1)

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

        ! Upack bilin map info 
        k = 1
        do i = 1, size(map)
            k1 = k + 4 - 1
            map(i)%iquad(1:4) = map_vec%iquad(k:k1)
            map(i)%alpha1(1)  = map_vec%alpha1(i)
            map(i)%alpha2(1)  = map_vec%alpha2(i)
            map(i)%alpha3(1)  = map_vec%alpha3(i)
            k = k1 + 1 
        end do 


        ! Deallocate map_vec vectors 
        call map_deallocate_map(map_vec) 

        return 

    end subroutine unpack_neighbors 

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
