module coord 
    ! This module wraps all submodules of the coord library
    ! By inserting `use coord` into a program, the user should
    ! have access to all public functions and subroutines defined
    ! in the modules included below. 

!     use ncio 

    use index 
    use interp1D 
    use interp2D 
    use interp_time 
    use polygons 
    use planet 
    use geodesic 
    use oblimap_projection_module

    use coordinates     ! grid_class, points_class and init methods    
    use coordinates_mapping
    use coordinates_mapping_conservative 
    use coordinates_mapping_scrip
    use grid_to_cdo 

    use subset2 
!     use subset1 
    use grid_gen 

!     use interp2D_conservative
!     use loess 
!     use gaussian_filter 
!     use mod_toms526

    implicit none 

    real(4), parameter :: COORD_VERSION = 1.0 

contains 

    function get_coord_version() result(ver)

        implicit none 

        real(4) :: ver 

        ver = COORD_VERSION

        return 

    end function 

end module coord
