module coordinates_mapping_scrip
    ! Define and perform mapping using the 
    ! SCRIP file format for storing mapping
    ! weights and neighbors.

    use coord_constants
    use ncio 

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


contains 
    
    subroutine map_scrip_field(map,var_name,var1,var2,method,fill,missing_value,mask_pack)
        ! Load a map_scrip_class object into memory
        ! from a netcdf file. 

        implicit none 

        type(map_scrip_class), intent(IN)    :: map 
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
        integer :: i, j 

        real(dp), allocatable :: var1_vec(:), var2_vec(:) 
        real(dp) :: area_tot, pt_ave, pt_var   
        integer  :: npt_now 


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
        
        ! First reset new interpolation address to zero 
        do n = 1, npts2 
            if (maskp(n)) then 
                do k = 1, map%num_links 
                    if (map%dst_address(k) .eq. n) then 
                        var2_vec(n) = 0.0d0
                        exit 
                    end if 
                end do 
            end if 
        end do 

        select case(trim(method))
        
            case ("fracarea")

                do n = 1,map%num_links
                    if (maskp(map%dst_address(n)) .and. var1_vec(map%src_address(n)) .ne. missing_val) then 
                        var2_vec(map%dst_address(n)) = var2_vec(map%dst_address(n)) +      &
                                      map%remap_matrix(1,n)*var1_vec(map%src_address(n))
                    end if
                end do

            case ("destarea")

                do n = 1, map%num_links
                    if (maskp(map%dst_address(n)) .and. var1_vec(map%src_address(n)) .ne. missing_val) then 
                        var2_vec(map%dst_address(n)) = var2_vec(map%dst_address(n)) +        &
                                        (map%remap_matrix(1,n)*var1_vec(map%src_address(n)))/ &
                                        (map%dst_grid_frac(map%dst_address(n)))
                    end if 
                end do

            case ("none")

                do n = 1, map%num_links
                    if (maskp(map%dst_address(n)) .and. var1_vec(map%src_address(n)) .ne. missing_val) then 
                        var2_vec(map%dst_address(n)) = var2_vec(map%dst_address(n)) +        &
                                       (map%remap_matrix(1,n)*var1_vec(map%src_address(n)))/ &
                                       (map%dst_grid_area(map%dst_address(n))*map%dst_grid_frac(map%dst_address(n)))
                    end if
                end do

            case DEFAULT 

                write(*,*) "map_scrip_field:: Error: method not recongized."
                write(*,*) "method = ", trim(method)
                stop 

        end select

        ! Send back to 2D array 
        var2 = reshape(var2_vec,[size(var2,1),size(var2,2)])

        return 

    end subroutine map_scrip_field

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

        write(*,*) "src_grid_size:    ", map%src_grid_size 
        write(*,*) "dst_grid_size:    ", map%dst_grid_size 
        write(*,*) "dst_grid_corners: ", map%dst_grid_corners 
        write(*,*) "src_grid_rank:    ", map%src_grid_rank 
        write(*,*) "dst_grid_rank:    ", map%dst_grid_rank 
        write(*,*) "num_links:        ", map%num_links 
        write(*,*) "num_wgts:         ", map%num_wgts 
        
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