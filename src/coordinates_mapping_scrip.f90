module coordinates_mapping_scrip
    ! Define and perform mapping using the 
    ! SCRIP file format for storing mapping
    ! weights and neighbors.

    use coord_constants
    use coordinates
    use ncio 
    use index 

    implicit none 

    type map_scrip_class

        character(len=256) :: src_name 
        character(len=256) :: dst_name 
        character(len=512) :: map_fname 

        ! ========================================================
        ! Variables below are defined to be consistent with 
        ! a SCRIP format netcdf file
        ! ========================================================
        
        integer :: src_grid_size
        integer :: dst_grid_size 
        integer :: dst_grid_corners
        integer :: src_grid_rank
        integer :: dst_grid_rank
        integer :: num_links
        integer :: num_wgts

        integer, allocatable :: src_grid_dims(:) 
        integer, allocatable :: dst_grid_dims(:) 
        real(8), allocatable :: src_grid_center_lat(:) 
        real(8), allocatable :: dst_grid_center_lat(:) 
        real(8), allocatable :: src_grid_center_lon(:) 
        real(8), allocatable :: dst_grid_center_lon(:) 
        real(8), allocatable :: dst_grid_corner_lat(:,:)
        real(8), allocatable :: dst_grid_corner_lon(:,:)
        integer, allocatable :: src_grid_imask(:) 
        integer, allocatable :: dst_grid_imask(:) 
        real(8), allocatable :: src_grid_area(:) 
        real(8), allocatable :: dst_grid_area(:)
        real(8), allocatable :: src_grid_frac(:) 
        real(8), allocatable :: dst_grid_frac(:)
        integer, allocatable :: src_address(:)
        integer, allocatable :: dst_address(:)
        real(8), allocatable :: remap_matrix(:,:) 

    end type 

    interface map_scrip_field
        module procedure map_scrip_field_double
        module procedure map_scrip_field_float
        module procedure map_scrip_field_integer 
    end interface 

    private 
    public :: map_scrip_class 
    public :: map_scrip_field
    public :: map_scrip_init
    public :: map_scrip_load 


contains 
    
    subroutine map_scrip_field_integer(map,var_name,var1,var2,method,fill,missing_value,mask_pack)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class), intent(IN), target :: map 
        character(len=*),      intent(IN)    :: var_name 
        integer,               intent(IN)    :: var1(:,:) 
        integer,               intent(INOUT) :: var2(:,:) 
        character(len=*),      intent(IN)    :: method
        logical,               intent(IN), optional :: fill            ! Fill cells with no available values?
        double precision,      intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,               intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate

        ! Local variables 
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field(map,var_name,var1dp,var2dp,method,fill,missing_value,mask_pack)

        var2 = int(var2dp) 

        return 

    end subroutine map_scrip_field_integer

    subroutine map_scrip_field_float(map,var_name,var1,var2,method,fill,missing_value,mask_pack)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class), intent(IN), target :: map 
        character(len=*),      intent(IN)    :: var_name 
        real(sp),              intent(IN)    :: var1(:,:) 
        real(sp),              intent(INOUT) :: var2(:,:) 
        character(len=*),      intent(IN)    :: method
        logical,               intent(IN), optional :: fill            ! Fill cells with no available values?
        double precision,      intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,               intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate

        ! Local variables 
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field(map,var_name,var1dp,var2dp,method,fill,missing_value,mask_pack)

        var2 = real(var2dp,sp) 

        return 

    end subroutine map_scrip_field_float

    subroutine map_scrip_field_double(map,var_name,var1,var2,method,fill,missing_value,mask_pack)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class), intent(IN), target :: map 
        character(len=*),      intent(IN)    :: var_name 
        real(8),               intent(IN)    :: var1(:,:) 
        real(8),               intent(INOUT) :: var2(:,:) 
        character(len=*),      intent(IN)    :: method
        logical,               intent(IN), optional :: fill            ! Fill cells with no available values?
        double precision,      intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,               intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate

        ! Local variables 
        integer :: n, k, npts1, npts2         
        logical          :: fill_pts
        double precision :: missing_val 
        logical          :: fixed_values  
        logical, allocatable  :: maskp(:)
        real(dp), allocatable :: area(:)
        integer :: i, j, j1, j2  

        real(dp), allocatable, target :: var1_vec(:)
        real(dp), allocatable :: var2_vec(:) 
        real(dp) :: area_tot, pt_ave, pt_var   
        integer  :: npt_now, num_links_now 

        real(dp), pointer :: var1_now(:) 
        real(dp), pointer :: wts1_now(:) 
        real(dp) :: wts1_tot 

        npts1 = size(var1,1)*size(var1,2)
        npts2 = size(var2,1)*size(var2,2)

        ! Confirm that source (var1) and destination (var2)
        ! arrays match map. 
        if (npts1 .ne. map%src_grid_size) then 
            write(*,*) "map_scrip_field:: Error: source array and map size do not match."
            write(*,*) "size(var1): ", size(var1,1), size(var1,2), " = ", npts1
            write(*,*) "map%src_grid_size = ", map%src_grid_size 
            stop 
        end if 
        if (npts2 .ne. map%dst_grid_size) then 
            write(*,*) "map_scrip_field:: Error: dest. array and map size do not match."
            write(*,*) "size(var2): ", size(var2,1), size(var2,2), " = ", npts2
            write(*,*) "map%dst_grid_size = ", map%dst_grid_size 
            stop 
        end if 
        
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

        j1 = 0 
        j2 = 0 

        ! Loop over target points
        do k = 1, npts2 

            if (maskp(k)) then 
                ! Only interpolate for desired target points 
        
                ! Find the range of link indices that correspond 
                ! to the current point k, ie, such that:
                ! map%dst_address(j1:j2) == k 
                ! Note: dst_address can be expected to be sorted 
                ! in ascending order.
                j1 = j2+1 

                ! Check index associated with this address. If 
                ! it is greater than the current index k, it 
                ! means this point has no interpolation links,
                ! so skip this iteration of the main loop.
                if (map%dst_address(j1) .gt. k) then 
                    j1 = j1-1 
                    cycle 
                end if 

                ! Given j1 is the start of the addresses associated 
                ! with the current index k, find the upper range 
                ! such that map%dst_address(j1:j2) == k and it 
                ! covers all addresses equal to k.
                do j = j1, map%num_links
                    if (map%dst_address(j) .eq. map%dst_address(j1) ) then 
                        j2 = j 
                    else 
                        exit 
                    end if 
                end do 

                ! Determine the number of links 
                num_links_now = j2-j1+1

                ! Define pointers to range of relevant input data and weights
                nullify(var1_now)
                nullify(wts1_now) 
                allocate(var1_now(num_links_now))
                allocate(wts1_now(num_links_now))

                ! Assign data and weights to pointers
                var1_now = var1_vec(map%src_address(j1:j2))
                wts1_now = map%remap_matrix(1,j1:j2)

                ! Calculate the total weight associated with this point,
                ! accounting for missing values in the source array.
                wts1_tot = sum(wts1_now,mask=var1_now .ne. missing_val)

                if (wts1_tot .gt. 0.0d0) then 
                    ! Interpolation data found, proceed to interpolate this point
                    
                    var2_vec(k) = 0.0d0 

                    select case(trim(method))

                        case("mean")
                            ! Calculate the area-weighted mean 

                            var2_vec(k) = sum((wts1_now/wts1_tot)*var1_now,mask=var1_now .ne. missing_val)

                        case("count")
                            ! Choose the most frequently occurring value, weighted by area

                            var2_vec(k) = maxcount(var1_now,wts1_now,missing_val)

                        case("stdev")
                            ! Calculate the weighted standard deviation 
                            ! using unbiased estimator correction 

                            npt_now = count(var1_now .ne. missing_val)

                            if (npt_now .gt. 2) then
                                ! Only calculate stdev for 2 or more input points

                                pt_ave      = sum((wts1_now/wts1_tot)*var1_now,mask=var1_now .ne. missing_val)
                                var2_vec(k) = (npt_now/(npt_now - 1.0)) &
                                               * sum((wts1_now/wts1_tot)*(var1_now-pt_ave)**2, & 
                                                                            mask=var1_now .ne. missing_val)
                                var2_vec(k) = sqrt(var2_vec(k))
                                
                            else
                                ! Otherwise assume standard deviation is zero 
                                var2_vec(k) = 0.0d0 

                            end if 

                        case DEFAULT 

                            write(*,*) "map_scrip_field:: Error: interpolation method not recognized."
                            write(*,*) "method = ", trim(method) 
                            stop 

                    end select 

                end if 

            end if 

        end do 

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

        return 

    end subroutine map_scrip_field_double

    subroutine map_scrip_init(map,src_name,dst_name,fldr,src_nc,load)
        ! Use cdo to generate scrip map based on grid 
        ! definitions. 

        ! 1. Assume that grid description text files already exist
        !    for each grid. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map     ! map object to be initialized
        character(len=*), intent(IN) :: src_name        ! Source grid name
        character(len=*), intent(IN) :: dst_name        ! Dest./target grid name
        character(len=*), intent(IN) :: fldr            ! Folder where grid desciptions can be found
        character(len=*), intent(IN) :: src_nc          ! Path to source netcdf file containing grid/variables (needed by cdo)
        logical,          intent(IN), optional :: load  ! Load map from file if available? 

        ! Local variables 
        character(len=512)  :: fnm1
        character(len=512)  :: fnm2
        character(len=512)  :: fnm_map 
        character(len=2048) :: cdo_cmd
        logical :: load_map 
        logical :: map_exists  
        logical :: cdo_success 

        ! Determine whether map file should be loaded if available 
        load_map = .TRUE. 
        if (present(load)) load_map = load 

        ! Step 1: call cdo to generate mapping weights in a scrip file 

        ! Generate grid description filenames 
        fnm1 = trim(fldr)//"/"//"grid_"//trim(src_name)//".txt"
        fnm2 = trim(fldr)//"/"//"grid_"//trim(dst_name)//".txt"

        ! Determine map filename from grid names and folder 
        fnm_map = gen_map_filename(src_name,dst_name,fldr)
        
        ! Check if scrip weights file already exists  
        inquire(file=trim(fnm_map),exist=map_exists)

        if ( (.not. map_exists) .or. (.not. load_map) ) then 
            ! If no map exists yet, or loading is not desired, 
            ! then call cdo to generate a new map file. 

            ! Define cdo command to generate mapping weights from 
            ! src grid (fnm1) to dest grid (fnm2) using example netcdf 
            ! grid file (src_nc) and storing weights in map file (fnm_map).
            ! cdo command output is redirected to a file '.tmpcdoout'.
            cdo_cmd = "cdo gencon,"//trim(fnm2)//" -setgrid,"//trim(fnm1)// &
                    " "//trim(src_nc)//" "//trim(fnm_map)//" &> .tmpcdoout"

            write(*,*) "cdo command: "
            write(*,*) trim(cdo_cmd) 

            write(*,"(a)",advance='no') "Calling via system call... "
            call system(cdo_cmd)
            write(*,*) "done." 

            ! Check if scrip weights file was written 
            inquire(file=trim(fnm_map),exist=cdo_success)

            if (.not. cdo_success) then 
                write(*,*) "map_scrip_init:: Error: scrip map file was not written. &
                & This may mean that the system call to cdo was unsucessful. Check the &
                &cdo log file '.tmpcdoout'."
                stop 
            end if 

        end if 

        ! Step 2: load map weights and initialize map_scrip_class object 
        call map_scrip_load(map,src_name,dst_name,fldr)

        return 

    end subroutine map_scrip_init

    subroutine map_scrip_load(map,src_name,dst_name,fldr)
        ! Load a map_scrip_class object into memory
        ! from a netcdf file. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map
        character(len=*), intent(IN) :: src_name
        character(len=*), intent(IN) :: dst_name 
        character(len=*), intent(IN) :: fldr  
        
        ! Local variables 
        integer, allocatable :: dims(:) 
        character(len=56), allocatable :: dim_names(:) 

        ! Define map names 
        map%src_name = trim(src_name) 
        map%dst_name = trim(dst_name) 

        ! Determine filename from grid names and folder 
        map%map_fname = gen_map_filename(src_name,dst_name,fldr)
        
        !write(*,*) "Loading SCRIP map from file: "//trim(filename) 
        !write(*,*) "" 

        call nc_dims(map%map_fname,"src_grid_center_lat",dim_names,dims)
        map%src_grid_size = dims(1) 
        call nc_dims(map%map_fname,"dst_grid_corner_lat",dim_names,dims)
        map%dst_grid_corners = dims(1) 
        map%dst_grid_size    = dims(2) 
        call nc_dims(map%map_fname,"src_grid_dims",dim_names,dims)
        map%src_grid_rank = dims(1) 
        call nc_dims(map%map_fname,"dst_grid_dims",dim_names,dims)
        map%dst_grid_rank = dims(1) 
        call nc_dims(map%map_fname,"remap_matrix",dim_names,dims)
        map%num_wgts  = dims(1)
        map%num_links = dims(2) 

        ! write(*,*) "src_grid_size:    ", map%src_grid_size 
        ! write(*,*) "dst_grid_size:    ", map%dst_grid_size 
        ! write(*,*) "dst_grid_corners: ", map%dst_grid_corners 
        ! write(*,*) "src_grid_rank:    ", map%src_grid_rank 
        ! write(*,*) "dst_grid_rank:    ", map%dst_grid_rank 
        ! write(*,*) "num_links:        ", map%num_links 
        ! write(*,*) "num_wgts:         ", map%num_wgts 
        
        ! Allocate map_scrip to match dimensions 
        call map_scrip_alloc(map)

        ! Load map from file 
        call nc_read(map%map_fname,"src_grid_dims",map%src_grid_dims)
        call nc_read(map%map_fname,"dst_grid_dims",map%dst_grid_dims)
        call nc_read(map%map_fname,"src_grid_center_lat",map%src_grid_center_lat)
        call nc_read(map%map_fname,"dst_grid_center_lat",map%dst_grid_center_lat)
        call nc_read(map%map_fname,"src_grid_center_lon",map%src_grid_center_lon)
        call nc_read(map%map_fname,"dst_grid_center_lon",map%dst_grid_center_lon)
        call nc_read(map%map_fname,"dst_grid_corner_lat",map%dst_grid_corner_lat)
        call nc_read(map%map_fname,"dst_grid_corner_lon",map%dst_grid_corner_lon)
        call nc_read(map%map_fname,"src_grid_imask",map%src_grid_imask)
        call nc_read(map%map_fname,"dst_grid_imask",map%dst_grid_imask)
        call nc_read(map%map_fname,"src_grid_area",map%src_grid_area)
        call nc_read(map%map_fname,"dst_grid_area",map%dst_grid_area)
        call nc_read(map%map_fname,"src_grid_frac",map%src_grid_frac)
        call nc_read(map%map_fname,"dst_grid_frac",map%dst_grid_frac)
        call nc_read(map%map_fname,"src_address",map%src_address)
        call nc_read(map%map_fname,"dst_address",map%dst_address)
        call nc_read(map%map_fname,"remap_matrix",map%remap_matrix)

        !write(*,*) "range(remap_matrix): ", minval(map%remap_matrix), maxval(map%remap_matrix)
        !write(*,*) "Loaded SCRIP weights."

        ! Summary print line
        !write(*,*) "Loaded SCRIP map from file: "//trim(filename) 
        write(*,*) "Loaded SCRIP map: "//trim(map%src_name)//" => "//trim(map%dst_name) 

        return 

    end subroutine map_scrip_load

    subroutine map_scrip_alloc(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 

        ! First deallocate everything for safety 

        call map_scrip_dealloc(map) 

        ! Proceed to allocation 

        allocate(map%src_grid_dims(map%src_grid_rank))
        allocate(map%dst_grid_dims(map%dst_grid_rank))
        allocate(map%src_grid_center_lat(map%src_grid_size))
        allocate(map%dst_grid_center_lat(map%dst_grid_size))
        allocate(map%src_grid_center_lon(map%src_grid_size))
        allocate(map%dst_grid_center_lon(map%dst_grid_size))
        allocate(map%dst_grid_corner_lat(map%dst_grid_corners,map%dst_grid_size))
        allocate(map%dst_grid_corner_lon(map%dst_grid_corners,map%dst_grid_size))
        allocate(map%src_grid_imask(map%src_grid_size))
        allocate(map%dst_grid_imask(map%dst_grid_size))
        allocate(map%src_grid_area(map%src_grid_size))
        allocate(map%dst_grid_area(map%dst_grid_size))
        allocate(map%src_grid_frac(map%src_grid_size))
        allocate(map%dst_grid_frac(map%dst_grid_size))
        allocate(map%src_address(map%num_links))
        allocate(map%dst_address(map%num_links))
        allocate(map%remap_matrix(map%num_wgts,map%num_links))

        return 

    end subroutine map_scrip_alloc

    subroutine map_scrip_dealloc(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 
        
        if (allocated(map%src_grid_dims))       deallocate(map%src_grid_dims)
        if (allocated(map%dst_grid_dims))       deallocate(map%dst_grid_dims)
        if (allocated(map%src_grid_center_lat)) deallocate(map%src_grid_center_lat)
        if (allocated(map%dst_grid_center_lat)) deallocate(map%dst_grid_center_lat)
        if (allocated(map%src_grid_center_lon)) deallocate(map%src_grid_center_lon)
        if (allocated(map%dst_grid_center_lon)) deallocate(map%dst_grid_center_lon)
        if (allocated(map%dst_grid_corner_lat)) deallocate(map%dst_grid_corner_lat)
        if (allocated(map%dst_grid_corner_lon)) deallocate(map%dst_grid_corner_lon)
        if (allocated(map%src_grid_imask))      deallocate(map%src_grid_imask)
        if (allocated(map%dst_grid_imask))      deallocate(map%dst_grid_imask)
        if (allocated(map%src_grid_area))       deallocate(map%src_grid_area)
        if (allocated(map%dst_grid_area))       deallocate(map%dst_grid_area)
        if (allocated(map%src_grid_frac))       deallocate(map%src_grid_frac)
        if (allocated(map%dst_grid_frac))       deallocate(map%dst_grid_frac)
        if (allocated(map%src_address))         deallocate(map%src_address)
        if (allocated(map%dst_address))         deallocate(map%dst_address)
        if (allocated(map%remap_matrix))        deallocate(map%remap_matrix)
        
        return 

    end subroutine map_scrip_dealloc

    function gen_map_filename(src_name,dst_name,fldr) result(filename)
        ! Output the standard map filename with input folder name
        implicit none 

        character(len=*), intent(IN) :: src_name
        character(len=*), intent(IN) :: dst_name 
        character(len=*), intent(IN) :: fldr 
        character(len=256) :: filename

        filename = trim(fldr)//"/scrip_"//trim(src_name)//"_"//trim(dst_name)//".nc"

        return

    end function gen_map_filename

end module coordinates_mapping_scrip