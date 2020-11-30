module interp2D_conservative
    
    use coord_constants 
    use coordinates 

    use planet 
    use oblimap_projection_module
    use polygons 
    use geodesic 
    use index 
    use interp2D 

    implicit none 

    type pt_wts_class 
        integer :: n  
        integer,  allocatable :: x(:), y(:), i(:)               ! Length of pts1/grid1 neighbors
        real(dp), allocatable :: dist(:), weight(:), area(:)
    end type 
! ** Now also defined in coordinates_mapping - this module is 
!    meant to stand alone though. 

    type map_conserv_class
!         type(grid_class) :: grid1, grid2 

        ! Neighbor info
        type(pt_wts_class), allocatable :: map(:)   ! Length of pts2/grid2

        ! Projection info 
        logical  :: is_cartesian, is_projection, is_same_map
        logical  :: is_lon180
        logical  :: is_grid

    end type

    interface map_field_conservative
        module procedure map_field_conservative_map
        module procedure map_field_conservative_grid
    end interface

    private 
    public :: map_conserv_class 
    public :: map_conservative_init
    public :: map_field_conservative
    public :: map_field_conservative_smooth
    !public :: calc_grid_total

contains 

    function gen_fraction_mask(grid1,grid2,varname,var1) result(var2)

        implicit none 

        type(grid_class), intent(IN)  :: grid1  ! Original grid information
        type(grid_class), intent(IN)  :: grid2  ! New grid 
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        logical,          intent(IN)  :: var1(:,:)      ! Input variable
        double precision              :: var2(grid2%G%nx,grid2%G%ny)      ! Output variable
        
        var2 = 0.d0 

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
        
        var2 = 0.d0 

        return 

    end function gen_stat_mask 

    subroutine map_field_conservative_smooth(map,map_hilo,grid1,grid2,varname,var1,var2,fill,missing_value,mask_pack)
        ! Map a low resolution field to a high resolution grid, but impose local smoothing
        ! that retains conservative properties of field

        implicit none 

        type(map_conserv_class), intent(IN) :: map        ! Map grid1=>grid2
        type(map_conserv_class), intent(IN) :: map_hilo   ! Map grid2=>grid1 
        type(grid_class),        intent(IN) :: grid1, grid2 ! Grids 

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
        real(dp) :: area(size(grid1%x,1),size(grid1%x,2))

        real(dp) :: low1(size(grid1%x,1),size(grid1%x,2))
        real(dp) :: low2(size(grid1%x,1),size(grid1%x,2))
        real(dp) :: low_corr(size(grid1%x,1),size(grid1%x,2))
        real(dp) :: high0(size(grid2%x,1),size(grid2%x,2))
        real(dp) :: high(size(grid2%x,1),size(grid2%x,2))
        real(dp) :: high_corr(size(grid2%x,1),size(grid2%x,2))
        integer :: k, n_iter, n_border, nxh, nyh  
        real(dp) :: target_val, current_val, err_percent, norm_val 
        real(dp) :: xlim(2), ylim(2) 

!         ! Generate conservation map for going from high to low resolution 
!         call map_conservative_init(map_hilo,grid2,grid1)
        
        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! Determine number of "border" grid points in high resolution field
        n_border = int( grid1%G%dx / grid2%G%dx )
        nxh      = grid2%G%nx 
        nyh      = grid2%G%ny 

        ! Calculate the target conservation value for the whole domain 
        xlim = [minval(grid1%G%x+grid1%G%dx/2.d0),maxval(grid1%G%x-grid1%G%dx/2.d0)]
        ylim = [minval(grid1%G%y+grid1%G%dy/2.d0),maxval(grid1%G%y-grid1%G%dy/2.d0)]
        target_val = calc_grid_total(grid1%G%x,grid1%G%y,var1,xlim=xlim,ylim=ylim)
        
        ! Normalization of error vals by sd(var1)
        norm_val = sqrt(sum(var1,mask=var1.ne.missing_val)**2) / dble(count(var1.ne.missing_val))

        ! Step 0: interpolate the field from low resolution to high resolution
        ! conservatively (results in blocky map) 
        call map_field_conservative(map,varname,var1,high0,fill,missing_val,mask_pack)
        high     = high0 
        low1     = var1 
        low_corr = 0.d0 

        n_iter = 50 

        ! Begin smoothing iteration
        do k = 1, n_iter 

            ! Step 1: re-interpolate back to low resolution
            low1 = low1 - low_corr 

            ! Step 2: use smooth interpolation to create a 
            ! new high res field with interpolation points 
            ! exactly on the values of the array `low1`
            high = interp_bilinear(x=grid1%G%x,y=grid1%G%y,z=low1, &
                                   xout=grid2%G%x,yout=grid2%G%y, &
                                   missing_value=missing_val)

            ! Step 3: integrate over high resolution points that fall into the
            ! coarse resolution grid. This is the target conservation value
            call map_field_conservative(map_hilo,varname,high,low2,missing_value=missing_val)

            ! Step 4: Find conservation error by comparing new low res field to original.
            low_corr = (low2-var1)

            ! Eliminate correction for border points
            low_corr(1,:) = 0.d0 
            low_corr(:,1) = 0.d0 
            low_corr(grid1%G%nx,:) = 0.d0 
            low_corr(:,grid1%G%ny) = 0.d0 
            
            ! Calculate the current conservation value for the whole domain
            current_val = calc_grid_total(grid2%G%x,grid2%G%y,high,xlim=xlim,ylim=ylim)
            err_percent = (current_val-target_val) / target_val * 100.d0 

            write(*,"(a,i4,5g11.2)") "map_field_conservative_smooth:: k, err_percent: ", k, &
                        target_val, current_val, err_percent, maxval(low_corr), norm_val

            ! Stop if stopping criterion is reached (eg, less than 1% of variance of field)
            if (maxval(low_corr) .lt. 0.01d0*norm_val) exit 
        end do 

        ! Apply the last correction field to the high resolution field before output
        call map_field_conservative(map,varname,low_corr,high_corr,fill,missing_val,mask_pack)
        
        ! After iterations store smooth high resolution grid for output 
        var2 = high - high_corr 

            
        return 

    end subroutine map_field_conservative_smooth

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

        ! Calculate grid total 
        tot = sum(var*weight*dx*dy)

        return 

    end function calc_grid_total

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
        character(len=*), intent(IN)    :: varname        ! Name of the variable being mapped
        double precision, intent(IN)    :: var1(:,:)      ! Input variable
        double precision, intent(INOUT) :: var2(:,:)      ! Output variable
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

                ii = map%map(i)%i  

                area = map%map(i)%area 
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

    end subroutine map_field_conservative_map 

    subroutine map_conservative_init(map,grid1,grid2)
        ! Initialize a conservative weighting map to go from grid1 to grid2

        implicit none 

        type(map_conserv_class), intent(OUT) :: map 
        type(grid_class),        intent(IN)  :: grid1, grid2 

        ! Local variables   
        integer :: npts2      
        double precision :: x1, y1, x2, y2 
        integer :: i, j 
        integer,  allocatable :: ii(:), jj(:)
        real(dp), allocatable :: area(:)
        logical :: is_hi2lo 

        type(points_class) :: pts1, pts2 
        real(dp) :: dx1, dy1, dx2, dy2 

        write(*,*) "Mapping weights: ", trim(grid1%name), " => ", trim(grid2%name)

        npts2 = grid2%G%nx*grid2%G%ny 
        dx1   = grid1%G%dx 
        dy1   = grid1%G%dy 
        dx2   = grid2%G%dx 
        dy2   = grid2%G%dy 
        
        ! Check if the grids are compatible for this mapping routine
        if (.not. same_projection(grid1%proj,grid2%proj)) then 
            write(*,*) "map_conservative_init:: error:  &
                       &Currently this subroutine only interpolates between grids/points &
                       &on the same projection. Try again."
            stop 
        end if 

        ! Determine direction of mapping 
        is_hi2lo = .TRUE. 
        if (grid2%G%dx .lt. grid1%G%dx) is_hi2lo = .FALSE. 

        ! Get points from grids for calculations 
        call grid_to_points(grid1,pts1)
        call grid_to_points(grid2,pts2)

        ! Allocate the map object to the size of the target grid
        if (allocated(map%map)) deallocate(map%map)
        allocate(map%map(pts2%npts))
        
        ! Allocate area to size of pts1
        allocate(area(pts1%npts))

        ! Loop over target grid, calculate area of points within each target grid point
        do i = 1, pts2%npts

            ! Get different range of neighbor indices depending on which direction
            ! of interpolation is being performed (to minimize neighborhood size)

            if (is_hi2lo) then
                ! From high to low resolution 

                call which(pts1%x .ge. pts2%x(i)-dx2 .and. pts1%x .le. pts2%x(i)+dx2 .and. & 
                           pts1%y .ge. pts2%y(i)-dy2 .and. pts1%y .le. pts2%y(i)+dy2, ii)
            else
                ! From low to high resolution 

                call which(pts1%x .ge. pts2%x(i)-dx1 .and. pts1%x .le. pts2%x(i)+dx1 .and. & 
                           pts1%y .ge. pts2%y(i)-dy1 .and. pts1%y .le. pts2%y(i)+dy1, ii)
            end if 

            area = 0.d0 
            area(ii) = interpconserv1_weights(x=pts1%x(ii),y=pts1%y(ii),dx=dx1,dy=dy1, &
                                              xout=pts2%x(i),yout=pts2%y(i),dxout=dx2,dyout=dy2)

            if (sum(area(ii)) .gt. 0.d0) then 
                ! If an interpolation point was found, calculate interpolation 

                ! Allocate map points for non-zero areas
                call which(area .gt. 0.d0, jj)
                call map_allocate_map(map%map(i),n=size(jj))
                map%map(i)%i = jj 

                ! Store the area values in the map
                map%map(i)%area = area(jj) 

                ! Set the other values to zero for now
                map%map(i)%dist   = 0.d0 
                map%map(i)%weight = 0.d0 

            else 
                ! Allocate one map point to indicate no neighbors are available 
                call map_allocate_map(map%map(i),n=1)
                map%map(i)%i = [1]
                map%map(i)%area = 0.d0 ! Fill values with area set to zero! 
                
                ! Set the other values to zero 
                map%map(i)%dist   = 0.d0 
                map%map(i)%weight = 0.d0 

            end if 

            ! Check progress
            if (mod(i,1000)==0) write(*,"(a,i10,a3,i12,a5,2g12.3)")  &
                                "  ",i, " / ",npts2,"   : ", sum(map%map(i)%area), dx2*dy2
        end do 

        write(*,*) "Mapping completed."
        write(*,*) 

        return 

    end subroutine map_conservative_init 

    subroutine map_allocate_map(mp,n)

        implicit none 

        type(pt_wts_class), intent(INOUT) :: mp 
        integer, intent(IN) :: n  

        mp%n = n  

        ! First deallocate if needed 
        if (allocated(mp%i)) deallocate(mp%i) 
        if (allocated(mp%dist)) deallocate(mp%dist) 
        if (allocated(mp%weight)) deallocate(mp%weight) 
        if (allocated(mp%area)) deallocate(mp%area) 
        
        allocate(mp%i(n))
        allocate(mp%dist(n), mp%weight(n), mp%area(n))

        return 

    end subroutine map_allocate_map

    function interpconserv1_weights(x,y,dx,dy,xout,yout,dxout,dyout,latlon) result(area)
        ! Calculate 1st order conservative interpolation for a 
        ! point given a vector of its neighbors (x,y)

        real(dp), intent(IN) :: x(:), y(:), dx, dy 
        real(dp), intent(IN) :: xout, yout, dxout, dyout 
        logical,  intent(IN), optional :: latlon
        real(dp) :: area(size(x,1))

        ! Local variables
        logical  :: is_latlon 
        type(polygon)       :: pol  
!         integer, parameter :: nx = 31, ny = 31, npts = nx*ny
        integer :: nx, ny, npts 
        real(dp) :: x1, y1
        integer  :: npts_in  
        integer :: i, j, now  
        real(dp) :: missing_val
        real(dp) :: area_target 

        is_latlon = .FALSE. 
        if (present(latlon)) is_latlon = latlon 

        ! Generate polygon representing boundaries of target point 
        pol = create_polygon(real([xout-dxout/2.d0,xout-dxout/2.d0,xout+dxout/2.d0,xout+dxout/2.d0]), &
                             real([yout-dyout/2.d0,yout+dyout/2.d0,yout+dyout/2.d0,yout-dyout/2.d0]))

        ! Loop over source points and get the area of each source
        ! polygon that is inside of the target polygon 
        ! - Save the area (absolute area, not fraction)
        area_target = dxout*dyout 
        area        = 0.d0 

        nx = max(1,int(dx/dxout))
        ny = max(1,int(dy/dyout))
        npts = nx*ny 

        do now = 1, size(x) 
            
            npts_in   = 0

            do j = 1, ny 
                do i = 1, nx 
                    x1 = (x(now)-dx/2.d0) + (dx)*dble(i-1)/dble(nx) + 0.5d0*1.d0/dble(nx)
                    y1 = (y(now)-dy/2.d0) + (dy)*dble(j-1)/dble(ny) + 0.5d0*1.d0/dble(ny) 
                    if (point_in_polygon(real(x1),real(y1),pol)) npts_in = npts_in+1
                end do
            end do 
            
            area(now) = dble(npts_in)/dble(npts) * dx*dy 

            ! If the source points area adds up to the target cell area, exit loop
            if (sum(area) .ge. area_target) exit 
        end do 

        return 

    end function interpconserv1_weights

end module interp2D_conservative 
