module coordinates_mapping_scrip
    ! Define and perform mapping using the 
    ! SCRIP file format for storing mapping
    ! weights and neighbors.

    use coord_constants
    use coordinates
    use grid_to_cdo 
    use ncio 
    use index 
    use interp2D 
    use gaussian_filter 

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
    public :: map_scrip_init_from_griddesc
    public :: map_scrip_load 
    public :: map_scrip_end

contains 
    
    subroutine map_scrip_field_integer(map,var_name,var1,var2,method,reset,missing_value, &
                                        mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        integer,                intent(IN)    :: var1(:,:) 
        integer,                intent(INOUT) :: var2(:,:) 
        character(len=*),       intent(IN)    :: method
        logical,                intent(IN), optional :: reset           ! Fill cells with no available values?
        double precision,       intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN), optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN), optional :: filt_method     ! Method to use for filtering
        double precision,       intent(IN), optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN), optional :: verbose         ! Print information
        
        ! Local variables 
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field(map,var_name,var1dp,var2dp,method,reset,missing_value, &
                                            mask_pack,fill_method,filt_method,filt_par)

        var2 = int(var2dp) 

        return 

    end subroutine map_scrip_field_integer

    subroutine map_scrip_field_float(map,var_name,var1,var2,method,reset,missing_value, &
                                            mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        real(sp),               intent(IN)    :: var1(:,:) 
        real(sp),               intent(INOUT) :: var2(:,:) 
        character(len=*),       intent(IN)    :: method
        logical,                intent(IN), optional :: reset           ! Reset var2 initially to missing_value?
        double precision,       intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN), optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN), optional :: filt_method     ! Method to use for filtering
        double precision,       intent(IN), optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN), optional :: verbose         ! Print information
        
        ! Local variables 
        real(dp), allocatable :: var1dp(:,:) 
        real(dp), allocatable :: var2dp(:,:) 
        
        allocate(var1dp(size(var1,1),size(var1,2)))
        allocate(var2dp(size(var2,1),size(var2,2)))
        
        var1dp = real(var1,dp)
        var2dp = real(var2,dp)
        
        call map_scrip_field(map,var_name,var1dp,var2dp,method,reset,missing_value, &
                                            mask_pack,fill_method,filt_method,filt_par)

        var2 = real(var2dp,sp) 

        return 

    end subroutine map_scrip_field_float

    subroutine map_scrip_field_double(map,var_name,var1,var2,method,reset,missing_value, &
                                                        mask_pack,fill_method,filt_method,filt_par,verbose)
        ! Map a variable field var1 from a src_grid to variable field var2 on dst_grid 

        ! Note: method='mean' is analogous to the method normalize_opt='fracarea' 
        ! desribed in the SCRIP documention (Fig. 2.4 in scripusers.pdf). The 
        ! other methods normalize_opt=['destarea','none'] have not been implemented.

        implicit none 

        type(map_scrip_class),  intent(IN), target :: map 
        character(len=*),       intent(IN)    :: var_name 
        real(8),                intent(IN)    :: var1(:,:) 
        real(8),                intent(INOUT) :: var2(:,:) 
        character(len=*),       intent(IN)    :: method
        logical,                intent(IN), optional :: reset           ! Reset var2 initially to missing_value?
        double precision,       intent(IN), optional :: missing_value   ! Points not included in mapping
        logical,                intent(IN), optional :: mask_pack(:,:)  ! Mask for where to interpolate
        character(len=*),       intent(IN), optional :: fill_method     ! Method to fill in remaining missing values
        character(len=*),       intent(IN), optional :: filt_method     ! Method to use for filtering
        double precision,       intent(IN), optional :: filt_par(:)     ! gaussian=[sigma,dx]; poisson=[tol]
        logical,                intent(IN), optional :: verbose         ! Print information
        
        ! Local variables 
        integer :: n, k, npts1, npts2         
        logical          :: reset_pts
        double precision :: missing_val 
        logical          :: verbose_out
        logical, allocatable  :: maskp(:)
        real(dp), allocatable :: area(:)
        integer :: i, j, j1, j2  

!cmw        real(dp), allocatable, target :: var1_vec(:)
        real(dp), allocatable :: var1_vec(:)
        real(dp), allocatable :: var2_vec(:) 
        real(dp) :: area_tot, pt_ave, pt_var   
        integer  :: npt_now, num_links_now 

!cmw        real(dp), pointer :: var1_now(:) 
!cmw        real(dp), pointer :: wts1_now(:) 
        real(dp), allocatable :: var1_now(:) 
        real(dp), allocatable :: wts1_now(:) 
        real(dp) :: wts1_tot 

        logical, allocatable  :: mask2(:,:) 

        integer :: npts_apply 
        real(dp) :: mean2, mean2b           ! Check mean before/after filtering

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
        
        ! By default, reset target grid points to missing values initially
        reset_pts = .TRUE. 
        if (present(reset)) reset_pts = reset 

        ! By defualt missing value is the coordinates package default value
        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! By default, not verbose output 
        verbose_out = .FALSE. 
        if (present(verbose)) verbose_out = verbose 

        ! By default, all var2 points are interpolated
        allocate(maskp(npts2))
        maskp = .TRUE. 
        if (present(mask_pack)) maskp = reshape(mask_pack,[npts2])

        ! Count total points to be applied 
        npts_apply = count(maskp)
        
        ! Store var1 in vector format
        allocate(var1_vec(npts1)) 
        var1_vec = reshape(var1,[npts1])

        ! Store var2 in vector format
        allocate(var2_vec(npts2))
        var2_vec = reshape(var2,[npts2])

        ! Reset output points to missing values
        if (reset_pts) var2_vec = missing_val 

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
!cmw                nullify(var1_now)
!cmw                nullify(wts1_now) 
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

                ! cmw
                deallocate(var1_now)
                deallocate(wts1_now)

            end if 

        end do 

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

        ! Allocate mask if needed 
        if (present(fill_method) .or. present(filt_method)) then

            allocate(mask2(size(var2,1),size(var2,2)))
            mask2 = reshape(maskp,[size(var2,1),size(var2,2)])

        end if 

        ! === Filling ===
        ! Fill in remaining missing values 

        if (present(fill_method)) then 

            select case(trim(fill_method))

                case("weighted")

                    call fill_weighted(var2,missing_val,n=6,mask=mask2)

                case DEFAULT ! eg "none"

                    ! Pass - no filling applied 

            end select

        end if 

        ! === Filtering ===
        ! Now perform filtering (smoothing) steps if
        ! the right arguments have been provided. 
        
        if (present(filt_method)) then 

            ! Calculate grid average before filtering 
            if (verbose_out .and. npts_apply .gt. 0) then 
                mean2 = sum(var2,mask=mask2) / real(npts_apply,dp)
            end if 

            select case(trim(filt_method))

                case("gaussian")

                    call filter_gaussian(var2,sigma=filt_par(1),dx=filt_par(2),mask=mask2)
                    
                case("gaussian-fast")

                    call filter_gaussian_fast(var2,sigma=filt_par(1),dx=filt_par(2),mask=mask2)
        

                case("poisson")

                    call filter_poisson(var2,mask=mask2,tol=filt_par(1), &
                                    missing_value=missing_val,wrapx=.FALSE.,verbose=.FALSE.)

                case DEFAULT ! == "none"

                    ! Pass - no filtering applied 

            end select 

            ! Calculate grid average after filtering 
            if (verbose_out .and. npts_apply .gt. 0) then 
                mean2b = sum(var2,mask=mask2) / real(npts_apply,dp)
            end if 

            if (verbose_out) then 
                ! Print summary of filtering 
                write(*,"(4a,2g14.5)") var_name, " - ",filt_method, ": mean[orig,filtered]: ", mean2, mean2b
            end if 

        end if 


        return 

    end subroutine map_scrip_field_double

    subroutine map_scrip_init(mps,grid1,grid2,fldr,load,clean)
        ! Generate mapping weights from grid1 to grid2

        implicit none 

        type(map_scrip_class), intent(INOUT) :: mps 
        type(grid_class), intent(IN)    :: grid1
        type(grid_class), intent(IN)    :: grid2 
        character(len=*), intent(IN), optional :: fldr      ! Directory in which to save/load map
        logical,          intent(IN), optional :: load      ! Whether loading is desired if map exists already
        logical,          intent(IN), optional :: clean     ! Whether to delete intermediate grid desc / grid files

        ! Local variables 
        logical :: load_file, fldr_exists, file_exists 
        character(len=256) :: mapfldr 
        character(len=512) :: src_nc 
        character(len=12)  :: xnm, ynm 
        character(len=512) :: filename, cmd 

        ! Load file if it exists by default
        load_file = .TRUE. 
        if (present(load)) load_file = load 

        mapfldr = "maps"
        if (present(fldr)) mapfldr = trim(fldr)

if (.FALSE.) then
    ! To do - add to scrip file? 

        ! ! Assign map constant information
        ! mps%name1         = trim(grid1%name) 
        ! mps%name2         = trim(grid2%name)
        ! mps%mtype         = trim(grid2%mtype)
        ! mps%units         = trim(grid2%units)
        ! mps%is_projection = grid2%is_projection 
        ! mps%is_cartesian  = grid2%is_cartesian
        ! mps%is_lon180     = grid2%is_lon180
        ! mps%planet        = grid2%planet
        ! mps%proj          = grid2%proj 
        ! mps%npts          = grid2%npts
        ! mps%nmax          = max_neighbors 
        ! mps%xy_conv       = grid2%xy_conv 

        ! ! Check if the same map is defined for both sets of points
        ! mps%is_same_map = compare_coord(grid1,grid2)
end if 

        ! Note: do not assign max distance here, save all distances
        ! up until the maximum number of neighbors
        ! Later, when loading map, let use choose max_distance
        ! In this way, less recalculation of maps will be needed
        ! when the max_distance changes.

        ! Determine if file matching these characteristics exists
        inquire(file=gen_map_filename(grid1%name,grid2%name,mapfldr),exist=file_exists)

        !! Now load map information from file if exists and is desired
        !! or else calculate weights and store in file. 
        if ( load_file .and. file_exists ) then 

            ! Read map from file
            call map_scrip_load(mps,grid1%name,grid2%name,mapfldr)

        else
            
            ! == Write grid description files to mapfldr

            call grid_cdo_write_desc_short(grid1,fldr=mapfldr) 
            call grid_cdo_write_desc_short(grid2,fldr=mapfldr) 

            ! ! Testing: overwrite projection grid desc with explicit cells
            ! if (grid1%is_projection) then 
            !     ! call grid_cdo_write_desc_explicit_proj(grid1%lon,grid1%lat,grid1%name,mapfldr)
            
            !     src_nc = trim(mapfldr)//"/grid_"//trim(grid1%name)//".nc"
            !     if ( (.not. grid1%is_projection) .and. (.not. grid1%is_cartesian) ) then 
            !         xnm = "lon"
            !         ynm = "lat"
            !     else 
            !         xnm = "xc"
            !         ynm = "yc"
            !     end if 

            !     call grid_write(grid1,fnm=src_nc,xnm=xnm,ynm=ynm,create=.TRUE.)

            !     call grid_cdo_write_desc_via_cdo(grid1%name,mapfldr,src_nc)

            ! end if 
            ! if (grid2%is_projection) then 
            !     ! call grid_cdo_write_desc_explicit_proj(grid2%lon,grid2%lat,grid2%name,mapfldr)
                
            !     src_nc = trim(mapfldr)//"/grid_"//trim(grid2%name)//".nc"
            !     if ( (.not. grid2%is_projection) .and. (.not. grid2%is_cartesian) ) then 
            !         xnm = "lon"
            !         ynm = "lat"
            !     else 
            !         xnm = "xc"
            !         ynm = "yc"
            !     end if 

            !     call grid_write(grid2,fnm=src_nc,xnm=xnm,ynm=ynm,create=.TRUE.)

            !     call grid_cdo_write_desc_via_cdo(grid2%name,mapfldr,src_nc)

            ! end if 

            
            ! == Generate source-grid file for use with `cdo gencon` call

            src_nc = trim(mapfldr)//"/grid_"//trim(grid1%name)//".nc"
            if ( (.not. grid1%is_projection) .and. (.not. grid1%is_cartesian) ) then 
                xnm = "lon"
                ynm = "lat"
            else 
                xnm = "xc"
                ynm = "yc"
            end if 

            call grid_write(grid1,fnm=src_nc,xnm=xnm,ynm=ynm,create=.TRUE.)


            ! == Generate the SCRIP map via a cdo call:

            call map_scrip_init_from_griddesc(mps,grid1%name,grid2%name,mapfldr,src_nc,load=.FALSE.)


            ! ==  Delete intermediate files if desired

            if (clean) then 

                ! Remove source grid file 
                cmd = "rm "//trim(src_nc)
                write(*,*) trim(cmd) 

                write(*,"(a)",advance='no') "Calling via system call... "
                call system(cmd)
                write(*,*) "done." 

                ! Remove grid description files 
                filename = trim(mapfldr)//"/"//"grid_"//trim(grid1%name)//".txt"
                cmd = "rm "//trim(filename)
                write(*,*) trim(cmd) 

                write(*,"(a)",advance='no') "Calling via system call... "
                call system(cmd)
                write(*,*) "done." 

                filename = trim(mapfldr)//"/"//"grid_"//trim(grid2%name)//".txt"
                cmd = "rm "//trim(filename)
                write(*,*) trim(cmd) 

                write(*,"(a)",advance='no') "Calling via system call... "
                call system(cmd)
                write(*,*) "done." 

            end if 

        end if 

        return 

    end subroutine map_scrip_init

    subroutine map_scrip_init_from_griddesc(map,src_name,dst_name,fldr,src_nc,load)
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

            ! Call cdo command via system call
            call call_system_cdo(cdo_cmd)

        end if 

        ! Step 2: load map weights and initialize map_scrip_class object 
        call map_scrip_load(map,src_name,dst_name,fldr)

        return 

    end subroutine map_scrip_init_from_griddesc

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
        call nc_dims(map%map_fname,"dst_grid_center_lat",dim_names,dims)
        map%dst_grid_size = dims(1) 
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lat")) then 
          call nc_dims(map%map_fname,"dst_grid_corner_lat",dim_names,dims)
          map%dst_grid_corners = dims(1) 
        else
          map%dst_grid_corners = 0
        endif
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
        ! Note: it seems dst_grid_corner_lat and dst_grid_corner_lon
        ! do not exist in all map files generated by cdo. This may be 
        ! related to the type of grid, but it is unclear. Nonetheless,
        ! these variables are so far not used in the map_field routine above.
        ! Below they are only read-in if available. 
        call nc_read(map%map_fname,"src_grid_dims",map%src_grid_dims)
        call nc_read(map%map_fname,"dst_grid_dims",map%dst_grid_dims)
        call nc_read(map%map_fname,"src_grid_center_lat",map%src_grid_center_lat)
        call nc_read(map%map_fname,"dst_grid_center_lat",map%dst_grid_center_lat)
        call nc_read(map%map_fname,"src_grid_center_lon",map%src_grid_center_lon)
        call nc_read(map%map_fname,"dst_grid_center_lon",map%dst_grid_center_lon)
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lat")) then 
            call nc_read(map%map_fname,"dst_grid_corner_lat",map%dst_grid_corner_lat)
        end if 
        if (nc_exists_var(map%map_fname,"dst_grid_corner_lon")) then 
            call nc_read(map%map_fname,"dst_grid_corner_lon",map%dst_grid_corner_lon)
        end if
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

    subroutine map_scrip_end(map)
        ! Allocate arrays in map_scrip_class. This
        ! routine assumes that size parameters within
        ! the object are already specified. 

        implicit none 

        type(map_scrip_class), intent(INOUT) :: map 
        
        ! Simply deallocate the map object 
        call map_scrip_dealloc(map) 

        return 

    end subroutine map_scrip_end

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
        if (map%dst_grid_corners>0) then
          allocate(map%dst_grid_corner_lat(map%dst_grid_corners,map%dst_grid_size))
          allocate(map%dst_grid_corner_lon(map%dst_grid_corners,map%dst_grid_size))
        endif
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
