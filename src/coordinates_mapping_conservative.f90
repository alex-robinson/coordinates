module coordinates_mapping_conservative

    use coordinates 
    use coordinates_mapping 
    use interp2D 
    use index
    use planet

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = -9999.0_dp 
    
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    private
    public :: map_field_conservative_map1

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








end module coordinates_mapping_conservative