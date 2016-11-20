module coordinates_mapping_conservative
    ! This module will eventually replace interp2D_conservative, 
    ! to handle generic cases of conservative mapping within the
    ! standard `map_class` framework. 
    
    use coord_constants
    use coordinates 
    use coordinates_mapping 
    use interp2D 
    use index
    use planet

    implicit none 

    private
    public :: map_field_conservative_map1
    public :: calc_grid_total 

contains 

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
 
                area = mp(i)%area 
                where (var1_vec(mp(i)%i) .eq. missing_val) area = 0.d0 

                write(*,*) i, minval(mp(i)%i), maxval(mp(i)%i), sum(mp(i)%area), sum(area)

                if (sum(area) .gt. 0.d0) then 
                    ! If an interpolation point was found, calculate interpolation 

                    var2_vec(i) = sum(var1_vec(mp(i)%i)*area,mask=area.gt.0.d0) &
                                  / sum(area,mask=area.gt.0.d0)

                end if 

            end if 

        end do 

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

!         if (fill_pts) then 
!             call fill_nearest(var2,missing_value=missing_val)
!         end if 

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

        ! Calculate grid total 
        tot = sum(var*weight*dx*dy)

        return 

    end function calc_grid_total






end module coordinates_mapping_conservative