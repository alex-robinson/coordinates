module interp2D_conservative
    
    use coordinates 
    use coordinates_mapping 

    use planet 
    use oblimap_projection_module
    use polygons 
    use geodesic 
    use index 
    use interp2D 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT

    type pt_wts_class 
        integer :: nx, ny  
        integer,  allocatable :: i(:), j(:)
        real(dp), allocatable :: dist(:,:), weight(:,:), area(:,:)
    end type 

    type map_conserv_class
        type(grid_class) :: grid1, grid2 

        ! Neighbor info
        type(pt_wts_class), allocatable :: map(:,:)

        ! Projection info 
        logical  :: is_cartesian, is_projection, is_same_map
        logical  :: is_lon180
        logical :: is_grid

    end type

    interface map_field_conservative
        module procedure :: map_field_conservative_map
        module procedure :: map_field_conservative_grid
    end interface

    private 
    public :: map_conserv_class 
    public :: map_conservative_init
    public :: map_field_conservative
    public :: map_field_conservative_smooth

contains 

    function gen_fraction_mask(grid1,grid2,varname,var1) result(var2)

        implicit none 

        type(grid_class), intent(IN)  :: grid1  ! Original grid information
        type(grid_class), intent(IN)  :: grid2  ! New grid 
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        logical,          intent(IN)  :: var1(:,:)      ! Input variable
        double precision              :: var2(grid2%G%nx,grid2%G%ny)      ! Output variable
        

        return 

    end function gen_fraction_mask 

    function gen_stat_mask(grid1,grid2,varname,var1,stat) result(var2)

        implicit none 

        type(grid_class), intent(IN)  :: grid1  ! Original grid information
        type(grid_class), intent(IN)  :: grid2  ! New grid 
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        double precision, intent(IN)  :: var1(:,:)      ! Input variable
        character(len=*), intent(IN)  :: stat           ! Statistic to calculate (mean,min,max,sd)
        double precision              :: var2(grid2%G%nx,grid2%G%ny)      ! Output variable
        

        return 

    end function gen_stat_mask 

    subroutine map_field_conservative_smooth(map,map_hilo,varname,var1,var2,fill,missing_value,mask_pack)
        ! Map a low resolution field to a high resolution grid, but impose local smoothing
        ! that retains conservative properties of field

        implicit none 

        type(map_conserv_class), intent(IN) :: map        ! Map grid1=>grid2
        type(map_conserv_class), intent(IN) :: map_hilo   ! Map grid2=>grid1 
        
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        double precision, intent(IN)  :: var1(:,:)      ! Input variable
        double precision, intent(OUT) :: var2(:,:)      ! Output variable
        logical,  optional :: fill
        double precision, intent(IN), optional :: missing_value  ! Points not included in mapping
        logical,          intent(IN), optional :: mask_pack(:,:)

        ! Local variables
        logical          :: fill_pts
        double precision :: missing_val 
        logical, allocatable :: maskp(:,:) 
        integer :: i, j 
        integer, allocatable :: ii(:), jj(:) 
        real(dp) :: area(size(map%grid1%x,1),size(map%grid1%x,2))

        real(dp) :: low1(size(map%grid1%x,1),size(map%grid1%x,2))
        real(dp) :: low2(size(map%grid1%x,1),size(map%grid1%x,2))
        real(dp) :: high0(size(map%grid2%x,1),size(map%grid2%x,2))
        real(dp) :: high(size(map%grid2%x,1),size(map%grid2%x,2))
        real(dp) :: high_corr(size(map%grid2%x,1),size(map%grid2%x,2))
        integer :: k, n_iter, n_border, nxh, nyh  
        real(dp) :: target_val, current_val, err_percent 

!         ! Generate conservation map for going from high to low resolution 
!         call map_conservative_init(map_hilo,map%grid2,map%grid1)
        
        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! Determine number of "border" grid points in high resolution field
        n_border = int( map%grid1%G%dx / map%grid2%G%dx )
        nxh      = map%grid2%G%nx 
        nyh      = map%grid2%G%ny 

        ! Calculate the target conservation value for the whole domain 
        target_val = sum(var1*map%grid1%G%dx*map%grid1%G%dy)

        ! Step 0: interpolate the field from low resolution to high resolution
        ! conservatively (results in blocky map) 
        call map_field_conservative(map,varname,var1,high0,fill,missing_val,mask_pack)
        high = high0 

        n_iter = 1 

        ! Begin smoothing iteration
        do k = 1, n_iter 

            ! Step 1: re-interpolate back to low resolution
            low1 = interp_nearest(x=map%grid2%G%x,y=map%grid2%G%y,z=high, &
                                  xout=map%grid1%G%x,yout=map%grid1%G%y, &
                                  missing_value=missing_val)

            ! Step 2: use smooth interpolation to create a 
            ! new high res field with interpolation points 
            ! exactly on the values of the array `low1`
            high = interp_bilinear(x=map%grid1%G%x,y=map%grid1%G%y,z=low1, &
                                   xout=map%grid2%G%x,yout=map%grid2%G%y, &
                                   missing_value=missing_val)

            ! Step 3: integrate over high resolution points that fall into the
            ! coarse resolution grid. This is the target conservation value
            call map_field_conservative(map_hilo,varname,high,low2,missing_value=missing_val)

            ! Step 4: Find conservation error (low2-low) and create a high res version of that.
            ! Note, you compare against the original low here, not low1 from above.
            call map_field_conservative(map,varname,low2-var1,high_corr,fill,missing_val,mask_pack)

            ! Step 5: Apply the correction 
            high = high - high_corr

            ! Assign border values from blocky map to avoid artefacts 
            high(1:n_border,:)           = high0(1:n_border,:)
            high((nxh-n_border+1):nxh,:) = high0((nxh-n_border+1):nxh,:)
            high(:,1:n_border)           = high0(:,1:n_border)
            high(:,(nyh-n_border+1):nyh) = high0(:,(nyh-n_border+1):nyh)

            ! Eliminate correction for border points 
            high_corr(1:n_border,:)           = 0.d0
            high_corr((nxh-n_border+1):nxh,:) = 0.d0
            high_corr(:,1:n_border)           = 0.d0
            high_corr(:,(nyh-n_border+1):nyh) = 0.d0

            ! Calculate the current conservation value for the whole domain,
            ! Stop if stopping criterion is reached 
            current_val = sum(high*map%grid2%G%dx*map%grid2%G%dy)
            err_percent = (current_val-target_val) / target_val * 100.d0 

            write(*,"(a,i4,4g10.1)") "map_field_conservative_smooth:: k, err_percent: ", k, &
                        target_val, current_val, err_percent, maxval(high_corr)

            ! If error is less than 1%, exit the iterative loop 
            if (abs(err_percent) .lt. 1.d0) exit 

        end do 

        ! After iterations store smooth high resolution grid for output 
        var2 = high 

            
        return 

    end subroutine map_field_conservative_smooth

    subroutine map_field_conservative_grid(grid1,grid2,varname,var1,var2,fill,missing_value,mask_pack)
        ! Map a field from grid1 to grid2 without defining any map variable beforehand.
        ! This is somewhat slower than performing the map initialization and 
        ! field mapping in one step, but it is more consistent with the general
        ! two-step approach (map_init, then map_field)

        implicit none 

        type(grid_class), intent(IN)  :: grid1, grid2   ! Input and output grid information
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        double precision, intent(IN)  :: var1(:,:)      ! Input variable
        double precision, intent(OUT) :: var2(:,:)      ! Output variable
        logical,  optional :: fill
        double precision, intent(IN), optional :: missing_value  ! Points not included in mapping
        logical,          intent(IN), optional :: mask_pack(:,:)

        ! Local variables
        type(map_conserv_class) :: map      ! Map grid1=>grid2
        
        ! Define the map locally 
        call map_conservative_init(map,grid1,grid2)

        ! Interpolate the field 
        call map_field_conservative_map(map,varname,var1,var2,fill,missing_value,mask_pack)

        return 

    end subroutine map_field_conservative_grid

    subroutine map_field_conservative_map(map,varname,var1,var2,fill,missing_value,mask_pack)
        ! Map a field from grid1 to grid2 using a predefined map of neighbor weights
        ! generated using `map_conserv_init` 

        implicit none 

        type(map_conserv_class), intent(IN) :: map      ! Map grid1=>grid2
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        double precision, intent(IN)  :: var1(:,:)      ! Input variable
        double precision, intent(OUT) :: var2(:,:)      ! Output variable
        logical,  optional :: fill
        double precision, intent(IN), optional :: missing_value  ! Points not included in mapping
        logical,          intent(IN), optional :: mask_pack(:,:)

        ! Local variables
        logical          :: fill_pts
        double precision :: missing_val 
        logical, allocatable :: maskp(:,:) 
        integer :: i, j 
        integer, allocatable :: ii(:), jj(:) 
        real(dp) :: area(size(map%grid1%x,1),size(map%grid1%x,2))

        ! Check if the grids are compatible for this mapping routine
        if (.not. same_projection(map%grid1%proj,map%grid2%proj)) then 
            write(*,*) "map_field_conservative:: error:  &
                       &Currently this subroutine only interpolates between grids &
                       &on the same projection. Try again."
            stop 
        end if 

        ! By default, fill in target grid points with missing values
        fill_pts = .TRUE. 
        if (present(fill)) fill_pts = fill 

        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! By default, all var2 points are interpolated
        allocate(maskp(size(var2,1),size(var2,2)))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = mask_pack 

        ! If fill is desired, initialize output points to missing values
        if (fill_pts) var2 = missing_val 

        ! A loop to get started 
        do j = 1, map%grid2%G%ny
            do i = 1, map%grid2%G%nx 

                if (maskp(i,j)) then 
                    ! Only interpolate for desired target points 

                    ii = map%map(i,j)%i 
                    jj = map%map(i,j)%j 

                    area = 0.d0 
                    area(ii,jj) = map%map(i,j)%area 

                    if (sum(area(ii,jj)) .gt. 0.d0) then 
                        ! If an interpolation point was found, calculate interpolation 

                        var2(i,j) = sum(var1(ii,jj)*area(ii,jj), &
                            mask=area(ii,jj).gt.0.d0 .and. var1(ii,jj).ne.missing_val) &
                                      / sum(area(ii,jj), &
                            mask=area(ii,jj).gt.0.d0 .and. var1(ii,jj).ne.missing_val)

                    end if 

                end if 

            end do 
        end do 

        if (fill_pts) then 
            call fill_nearest(var2,missing_value=missing_val)
        end if 

        return 

    end subroutine map_field_conservative_map 

    subroutine map_conservative_init(map,grid1,grid2)
        ! Initialize a conservative weighting map to go from grid1 to grid2

        implicit none 

        type(map_conserv_class), intent(OUT) :: map 
        type(grid_class),        intent(IN)  :: grid1, grid2 

        ! Local variables          
        double precision :: x1, y1, x2, y2 
        integer :: i, j 
        integer, allocatable :: ii(:), jj(:) 
        real(dp) :: area(size(grid1%x,1),size(grid1%x,2))
        logical :: is_hi2lo 

        write(*,*) "Mapping weights: ", trim(grid1%name), " => ", trim(grid2%name)

        ! Determine direction of mapping 
        is_hi2lo = .TRUE. 
        if (grid2%G%dx .lt. grid1%G%dx) is_hi2lo = .FALSE. 

        ! Assign grids to map for later use 
        map%grid1 = grid1 
        map%grid2 = grid2 

        ! Allocate the map object to the size of the target grid
        if (allocated(map%map)) deallocate(map%map)
        allocate(map%map(grid2%G%nx,grid2%G%ny))
                    
        ! Loop over target grid, calculate area of points within each target grid point
        do j = 1, grid2%G%ny
            if (mod(j,10).eq.0) write(*,*) "j / ny: ", j, grid2%G%ny 

            do i = 1, grid2%G%nx 

                ! Get different range of neighbor indices depending on which direction
                ! of interpolation is being performed (to minimize neighborhood size)

                if (is_hi2lo) then
                    ! From high to low resolution 

                    call which(grid1%G%x .ge. grid2%x(i,j)-grid2%G%dx .and. grid1%G%x .le. grid2%x(i,j)+grid2%G%dx,ii)
                    call which(grid1%G%y .ge. grid2%y(i,j)-grid2%G%dy .and. grid1%G%y .le. grid2%y(i,j)+grid2%G%dy,jj)
                 
                else
                    ! From low to high resolution 

                    call which(grid1%G%x .ge. grid2%x(i,j)-grid1%G%dx .and. grid1%G%x .le. grid2%x(i,j)+grid1%G%dx,ii)
                    call which(grid1%G%y .ge. grid2%y(i,j)-grid1%G%dy .and. grid1%G%y .le. grid2%y(i,j)+grid1%G%dy,jj)
                 
                end if 

                area = 0.d0 
                area(ii,jj) = interpconserv1_weights(x=grid1%x(ii,jj),y=grid1%y(ii,jj), &
                                                  dx=grid1%G%dx,dy=grid1%G%dy, &
                                                  xout=grid2%x(i,j),yout=grid2%y(i,j), &
                                                  dxout=grid2%G%dx,dyout=grid2%G%dy)

                if (sum(area(ii,jj)) .gt. 0.d0) then 
                    ! If an interpolation point was found, calculate interpolation 

                    ! Allocate one map point to indicate no neighbors are available 
                    call map_allocate_map(map%map(i,j),nx=size(ii),ny=size(jj))
                    map%map(i,j)%i = ii 
                    map%map(i,j)%j = jj  

                    ! Store the area values in the map
                    map%map(i,j)%area = area(ii,jj) 

!                     write(*,*) i, j, sum(area(ii,jj)), sum(map%map(i,j)%area)

                    ! Set the others to zero 
                    map%map(i,j)%dist   = 0.d0 
                    map%map(i,j)%weight = 0.d0 

                else 
                    ! Allocate one map point to indicate no neighbors are available 
                    call map_allocate_map(map%map(i,j),nx=2,ny=2)
                    map%map(i,j)%i = [1,2]    ! Fill values with area set to zero!
                    map%map(i,j)%j = [1,2]    ! Fill values with area set to zero! 

                    ! Set the others to zero 
                    map%map(i,j)%area   = 0.d0 
                    map%map(i,j)%dist   = 0.d0 
                    map%map(i,j)%weight = 0.d0 

                end if 

            end do 
        end do 

        write(*,*) "Mapping completed."
        write(*,*) 

        return 

    end subroutine map_conservative_init 

    subroutine map_allocate_map(mp,nx,ny)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mp 
        integer, intent(IN) :: nx, ny  

        mp%nx = nx 
        mp%ny = ny 

        ! First deallocate if needed 
        if (allocated(mp%i)) deallocate(mp%i) 
        if (allocated(mp%j)) deallocate(mp%j) 
        if (allocated(mp%dist)) deallocate(mp%dist) 
        if (allocated(mp%weight)) deallocate(mp%weight) 
        if (allocated(mp%area)) deallocate(mp%area) 
        
        allocate(mp%i(nx),mp%j(ny))
        allocate(mp%dist(nx,ny), mp%weight(nx,ny), mp%area(nx,ny))

        return 

    end subroutine map_allocate_map

    function interpconserv1_weights(x,y,dx,dy,xout,yout,dxout,dyout,latlon,missing_value) result(area)
        ! Calculate 1st order conservative interpolation for a 
        ! point given its nearest neighbors  
        ! Better for going from high res to low res (hi2lo)
        real(dp), intent(IN) :: x(:,:), y(:,:), dx, dy 
        real(dp), intent(IN) :: xout, yout, dxout, dyout 
        logical,  intent(IN), optional :: latlon
        real(dp), intent(IN), optional :: missing_value  
        real(dp) :: area(size(x,1),size(x,2))

        ! Local variables
        logical  :: is_latlon 
        type(polygon)       :: pol  
        integer, parameter :: nx = 10, ny = 10, npts = nx*ny
        real(dp) :: x1, y1
        integer  :: npts_in 
        integer :: i, j, now  
        real(dp) :: missing_val
        real(dp) :: x_vec(size(x,1)*size(x,2)), y_vec(size(x,1)*size(x,2)) 
        real(dp) :: area_vec(size(x,1)*size(x,2))
        real(dp) :: area_target 

        is_latlon = .FALSE. 
        if (present(latlon)) is_latlon = latlon 

        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! Generate polygon representing boundaries of target point 
        pol = create_polygon(real([xout-dxout/2.d0,xout-dxout/2.d0,xout+dxout/2.d0,xout+dxout/2.d0]), &
                             real([yout-dyout/2.d0,yout+dyout/2.d0,yout+dyout/2.d0,yout-dyout/2.d0]))

        x_vec = reshape(x,[size(x_vec)])
        y_vec = reshape(y,[size(y_vec)])
        
        ! Loop over source points and get the area of each source
        ! polygon that is inside of the target polygon 
        ! - Save the area (absolute area, not fraction)
        area_target = dxout*dyout 
        area_vec    = 0.d0 

        do now = 1, size(x_vec) 
         
            npts_in = 0
            do i = 1, nx
            do j = 1, ny 
                x1 = (x_vec(now)-dx/2.d0) + (dx)*dble(i-1)/dble(nx) + (dx)*0.5d0/dble(nx-1) 
                y1 = (y_vec(now)-dy/2.d0) + (dy)*dble(j-1)/dble(ny) + (dy)*0.5d0/dble(ny-1) 
                if (point_in_polygon(real(x1),real(y1),pol)) npts_in = npts_in+1

            end do
            end do 
            
            area_vec(now) = dble(npts_in)/dble(npts) * dx*dy 

            ! If the source points area adds up to the target cell area, exit loop
            if (sum(area_vec) .ge. area_target) exit 
        end do 

        ! Return 2D area 
        area = reshape(area_vec,[size(x,1),size(x,2)])

        return 

    end function interpconserv1_weights

end module interp2D_conservative 
