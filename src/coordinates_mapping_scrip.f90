module coordinates_mapping_scrip
    ! Define and perform mapping using the 
    ! SCRIP file format for storing mapping
    ! weights and neighbors.

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

    
    subroutine map_scrip_field(map,var_name,dst,src,,dst_mask,method)
        ! Load a map_scrip_class object into memory
        ! from a netcdf file. 

        implicit none 

        type(map_scrip_class), intent(IN)  :: map 
        character(len=*),      intent(IN)  :: var_name 
        real(8),               intent(OUT) :: dst(:,:) 
        real(8),               intent(IN)  :: src(:,:) 
        character(len=*),      intent(IN), optional :: method

        ! Local variables 
        integer :: n 
        character(len=56) :: normalize_opt
        real(8), allocatable :: dst_array(:) 
        real(8), allocatable :: src_array(:) 
        
        allocate(dst_array(map%dst_grid_size))
        allocate(src_array(map%src_grid_size))
        
        dst_array = 0.0
        src_array = reshape(src,[map%src_grid_size])

        select case(trim(normalize_opt))
        
            case ('fracarea')

                do n = 1,map%num_links
                    dst_array(map%dst_address(n)) = dst_array(map%dst_address(n)) +      &
                                    map%remap_matrix(1,n)*src_array(map%src_address(n))
                end do

            case ('destarea')

                do n = 1, map%num_links
                    dst_array(map%dst_address(n)) = dst_array(map%dst_address(n)) +        &
                                    (map%remap_matrix(1,n)*src_array(map%src_address(n)))/ &
                                    (map%dst_grid_frac(map%dst_address(n)))
                end do

            case ('none')

                do n = 1, map%num_links
                    dst_array(map%dst_address(n)) = dst_array(map%dst_address(n)) +        &
                                    (map%remap_matrix(1,n)*src_array(map%src_address(n)))/ &
                                    (map%dst_grid_area(map%dst_address(n))*map%dst_grid_frac(map%dst_address(n)))
                end do

            case DEFAULT 

                write(*,*) "map_scrip_field:: Error: normalize_opt case not recongized."
                write(*,*) "normalize_opt = ", normalize_opt
                stop 

        end select

        dst = reshape(dst_array,[size(dst,1),size(dst,2)])
        
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