!> ####################################################################
!! **Module**     : coordinates_sigma \n
!! **Author**     : Alex Robinson \n 
!! **Purpose**    : This module defines the coordinate classes to handle
!!                  vector and point transformations to/from a cartesian
!!                  and sigma (0 to 1) coordinate system.
!! ####################################################################
module coordinates_sigma 

    use coord_constants 

    implicit none 

    type points_sigma_class 

        character (len=128) :: name     ! name of this domain (sigma_std,sigma_exp, etc)
        character (len=128) :: units    ! units of the axes
        
        ! Projection parameters
        logical :: is_dz_const

        ! Points information
        integer :: npts
        real(dp), allocatable, dimension(:)   :: z, zeta    ! cartesian, sigma
        real(dp), allocatable, dimension(:)   :: dz, dzeta  
        integer,  allocatable, dimension(:)   :: border
        real(dp) :: xy_conv 

    end type 

    private
    public :: points_sigma_class 
    public :: points_sigma_init
    public :: points_sigma_write, points_sigma_print

contains 



    subroutine points_sigma_init(pts,name,units)

        implicit none 

        type(points_sigma_class) :: pts 
        character(len=*) :: name, units 


        ! Make sure we can convert the units of the points as needed
        select case(trim(pts%units))
            case("kilometers","km")
                pts%xy_conv = 1.d3  ! To convert from km => m
            case DEFAULT
                pts%xy_conv = 1.d0 
        end select

        ! Reallocate point vectors
        if (allocated(pts%z))      deallocate(pts%z)
        if (allocated(pts%zeta))   deallocate(pts%zeta)
        if (allocated(pts%dz))     deallocate(pts%dz)
        if (allocated(pts%dzeta))  deallocate(pts%dzeta)

        allocate(pts%z(pts%npts), pts%zeta(pts%npts))
        allocate(pts%dz(pts%npts),pts%dzeta(pts%npts))


        dzeta_c = 1.0_dp/real(KCMAX,dp)

        if (DEFORM >= eps) then

            flag_aa_nonzero = .true.   ! non-equidistant grid

            aa = DEFORM
            ea = exp(aa)

            kc=0
            zeta_c(kc)         = 0.0_dp
            eaz_c(kc)          = 1.0_dp
            eaz_c_quotient(kc) = 0.0_dp

            do kc=1, KCMAX-1
                zeta_c(kc) = kc*dzeta_c
                eaz_c(kc)  = exp(aa*zeta_c(kc))
                eaz_c_quotient(kc) = (eaz_c(kc)-1.0_dp)/(ea-1.0_dp)
            end do

            kc=KCMAX
            zeta_c(kc)         = 1.0_dp
            eaz_c(kc)          = exp(aa)
            eaz_c_quotient(kc) = 1.0_dp

        else

            flag_aa_nonzero = .false.   ! equidistant grid

            aa = 0.0_dp
            ea = 1.0_dp

            kc=0
            zeta_c(kc)         = 0.0_dp
            eaz_c(kc)          = 1.0_dp
            eaz_c_quotient(kc) = 0.0_dp

            do kc=1, KCMAX-1
                zeta_c(kc) = kc*dzeta_c
                eaz_c(kc)  = 1.0_dp
                eaz_c_quotient(kc) = zeta_c(kc)
            end do

            kc=KCMAX
            zeta_c(kc)         = 1.0_dp
            eaz_c(kc)          = 1.0_dp
            eaz_c_quotient(kc) = 1.0_dp

        end if


        return 

    end subroutine points_sigma_init 

    subroutine points_sigma_update(pts,z0,z1)

        implicit none 

        type(points_sigma_class) :: pts 
        real(dp) :: z0, z1 



        return 

    end subroutine points_sigma_update 

end module coordinates_sigma 
