! File name: oblimap_projection_module.f90
!
! Copyright (C) 2013 Thomas Reerink & Michael Kliphuis.
!
! This file is distributed under the terms of the
! GNU General Public License.
!
! This file is part of OBLIMAP 2.0
!
! OBLIMAP's scientific documentation and its first open source
! release (see the supplement) is published at:
! http://www.geosci-model-dev.net/3/13/2010/gmd-3-13-2010.html
!
! OBLIMAP is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! OBLIMAP is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with OBLIMAP. If not, see <http://www.gnu.org/licenses/>.
!
!
! OBLIMAP is maintained by:
!
! Thomas Reerink
! Institute for Marine and Atmospheric Research Utrecht (IMAU)
! Utrecht University
! Princetonplein 5
! 3584 CC Utrecht
! The Netherlands
!
! email: tjreerink@gmail.com
!

!! ####################################################################
!! MODIFICATIONS
!! ajr, 2013-09-17: encapsulated funcionality further
!! ajr, 2015-03-05: Added additional projection routines: LEA, etc.
!! 
!! ####################################################################

MODULE oblimap_projection_module

    use coord_constants
    use planet 

    implicit none 
   
    type projection_class

        character(len=256) :: name, method  

        ! Projection parameters
        real(dp) :: lambda, lambda_M 
        real(dp) :: phi, phi_M 
        real(dp) :: alpha, alpha_stereographic
        real(dp) :: x_e, y_n 

        ! Planet parameters 
        character(len=256) :: planet_name 
        logical            :: is_sphere 

        ! Error tolerance level
        real(dp)      :: ebs = 1.e-8_dp

        ! Geometric constants (spherical planet)
        real(dp)      :: R, earth_radius
        
        ! Geometric constants (ellipsoid planet, wG84)
        real(dp)      :: a
        real(dp)      :: e
        real(dp)      :: f
        real(dp)      :: am
        real(dp)      :: akm
        real(dp)      :: chi_M
        
        real(dp)      :: q_M      
        real(dp)      :: q_polar  
        real(dp)      :: beta_M   
        real(dp)      :: R_q_polar
        real(dp)      :: D 

    end type 

    private
    public :: projection_class, projection_init, same_projection
    public :: optimal_alpha
    public :: oblimap_projection, oblimap_projection_inverse

contains
  
  subroutine projection_init(proj,name,planet,lambda,phi,alpha,x_e,y_n,method)
    
    implicit none
    
    type(projection_class) :: proj 
    character(len=*)       :: name
    type(planet_class)     :: planet 
    real(dp), optional     :: lambda, phi, alpha 
    real(dp), optional     :: x_e, y_n 
    character(len=*), optional :: method 
    character(len=256) :: proj_method 

    proj%name = trim(name) 

    ! Save planet info locally for ease of use 
    proj%f   = planet%f 
    proj%a   = planet%a 
    proj%e   = planet%e 
    proj%R   = planet%R 
    proj%earth_radius = planet%R 

    proj%planet_name = trim(planet%name)
    proj%is_sphere   = planet%is_sphere 

    ! Determine best default project method based on planet information
    ! and map type

    if (trim(proj%name) .eq. "stereographic" .or. &
        trim(proj%name) .eq. "polar_stereographic") then 

        if (proj%is_sphere) then 
            proj%method = "oblique_sg_projection"
        else 
            proj%method = "oblique_sg_projection_ellipsoid_snyder" 
        end if 

    else if (trim(proj%name) .eq. "lambert_azimuthal_equal_area" .or. &
             trim(proj%name) .eq. "lambert_azimuthal_equal_area") then 
        
        if (proj%is_sphere) then 
            proj%method = "oblique_laea_projection_snyder"
        else 
            proj%method = "oblique_laea_projection_ellipsoid_snyder" 
        end if 
        
    else
        proj%method = "Undefined"
    end if 

    ! Determine actual method to use based on user input if given
    if (present(method)) proj%method = trim(method)

!     proj%R            = 6.371221E6_dp ! Radius of the planet (for now it's Earth)
!     proj%earth_radius = 6.371221E6_dp ! redundant def of Earth radius (why not use generic R?)

!     proj%a      = 6378137._dp    ! equatorial ellipsoid radius,   a in Snyder, WGS84
!     proj%e      = 0.081819191_dp ! eccentricity of the ellipsoid, e in Snyder, WGS84, see Snyder p. 13
!     proj%f      = 1.0_dp/298.257223563_dp  ! Flattening of the ellipsoid
    
    ! Default dummy projection values
    proj%phi    = 0.0_dp 
    proj%lambda = 0.0_dp 
    proj%alpha  = 0.0_dp 
    proj%x_e    = 0.0_dp
    proj%y_n    = 0.0_dp 

    ! Get projection center if provided
    ! Also add additional projection rotation parameters if they exist
    ! (currently not used, but need to be stored in netcdf files)
    if (present(phi))    proj%phi    = phi
    if (present(lambda)) proj%lambda = lambda 
    if (present(x_e))    proj%x_e    = x_e 
    if (present(y_n))    proj%y_n    = y_n 

    if (trim(proj%method) .ne. "Undefined") then 
        ! latitude_longitude grid assumed.

        ! === Determine alpha ================================

        if ( trim(proj%name) .eq. "polar_stereographic" ) then 
            ! Use alpha determined from phi directly 
            ! If desired phi = 90 or -90, make sure alpha is non-zero

            proj%alpha = max(90.0_dp - dabs(proj%phi),0.0001_dp) 

        else if (present(alpha)) then 
            ! Use user-determined alpha

            if ( alpha .gt. 0.0_dp ) then 

                proj%alpha  = alpha 

            else 

                write(*,*) "oblimap_projection_module::"
                write(*,*) "alpha must be greater than zero, alpha = ", alpha
                stop 
            
            end if 

        else 
            ! Use optimal alpha 

            write(*,*) "oblimap_projection_module::"
            write(*,*) "    Optimal alpha calculation not active."
            write(*,*) "    Please specify a positive value for alpha."
            write(*,*)
            stop

        end if 

    end if 

    proj%phi_M               = degrees_to_radians * proj%phi
    proj%lambda_M            = degrees_to_radians * proj%lambda
    proj%alpha_stereographic = degrees_to_radians * proj%alpha

    ! See equations (14-15) and (21-27) on page 160 in Snyder (1987), akm corresponds with 2a*k0*m_1 in Snyder:
    proj%am     = proj%a * (COS(proj%phi_M) / DSQRT(1._dp - (proj%e * SIN(proj%phi_M))**2))
    proj%akm    = (1._dp + COS(proj%alpha_stereographic)) * proj%am
    ! See equations (3-1a) on page 160 in Snyder (1987),  chi_M corresponds with chi_1 in Snyder:
    if (dabs(proj%phi) .lt. 90.0_dp) then 
      proj%chi_M  = 2._dp * ATAN(DSQRT(((1._dp +       SIN(proj%phi_M)) / (1._dp -       SIN(proj%phi_M))) * &
                           ((1._dp - proj%e * SIN(proj%phi_M)) / (1._dp + proj%e * SIN(proj%phi_M)))**(proj%e))) - 0.5_dp * pi
    else if (proj%phi .eq. 90.0_dp) then 
      proj%chi_M = pi/2.0_dp   ! North Pole
    else
      proj%chi_M = -pi/2.0_dp  ! South Pole
    end if

    ! See equation (3-12) on page 187 in Snyder (1987):
    proj%q_M   = (1._dp - proj%e**2) * ((SIN(proj%phi_M) / &
                   (1._dp - (proj%e * SIN(proj%phi_M))**2)) - (1._dp / (2._dp * proj%e)) * &
                   LOG((1._dp - proj%e * SIN(proj%phi_M)) / (1._dp + proj%e * SIN(proj%phi_M)))) 
    ! See equation (3-12) on page 187 in Snyder (1987):
    proj%q_polar   = (1._dp - proj%e**2) * ((1._dp / (1._dp - proj%e**2)) - (1._dp / (2._dp * proj%e)) * &
                      LOG((1._dp - proj%e) / (1._dp + proj%e)))
    
    ! See equation (3-11) on page 187 in Snyder (1987):
    proj%beta_M    = ASIN(proj%q_M / proj%q_polar)
    ! See equation (3-13) on page 187 in Snyder (1987):
    proj%R_q_polar = proj%a * DSQRT(0.5_dp * proj%q_polar)
    ! See equation (24-20) on page 187 in Snyder (1987):
    proj%D         = proj%am / (proj%R_q_polar * COS(proj%phi_M))

    ! if (trim(proj%method) .ne. "Undefined") call projection_print(proj)
    
    return

  end subroutine projection_init

  subroutine projection_print(proj)

    implicit none 

    type(projection_class) :: proj 


        character(len=256) :: name, method  


        ! Planet parameters 
        character(len=256) :: planet_name 
        logical            :: is_sphere 

        ! Error tolerance level
        real(dp)      :: ebs = 1.e-8_dp

        ! Geometric constants (spherical planet)
        real(dp)      :: R, earth_radius
    
    write(*,*) "== Projection summary ========================="
    write(*,"(a20,a)")     "name = ", trim(proj%name)
    write(*,"(a20,a)") "method = ", trim(proj%method)
    write(*,"(a20,a)") "planet = ", trim(proj%planet_name)
    write(*,"(a20,l2)") "is sphere? = ", proj%is_sphere
    write(*,"(a20,g15.3)") "lambda = ", proj%lambda
    write(*,"(a20,g15.3)") "phi = ", proj%phi
    write(*,"(a20,g15.3)") "alpha = ", proj%alpha
    write(*,"(a20,g15.3)") "x_e = ", proj%x_e
    write(*,"(a20,g15.3)") "y_n = ", proj%y_n
    write(*,*) 
    write(*,"(a20,g15.3)") "lambda_M = ", proj%lambda_M
    write(*,"(a20,g15.3)") "phi_M = ", proj%phi_M
    write(*,"(a20,g15.3)") "alpha_ster = ", proj%alpha_stereographic
    write(*,*)
    write(*,"(a20,g15.3)") "R = ", proj%R 
    write(*,"(a20,g15.3)") "earth_radius = ", proj%earth_radius 
    write(*,"(a20,g15.3)") "a = ", proj%a 
    write(*,"(a20,g15.3)") "e = ", proj%e 
    write(*,"(a20,g15.3)") "f = ", proj%f 
    write(*,"(a20,g15.3)") "am = ", proj%am
    write(*,"(a20,g15.3)") "akm = ", proj%akm 
    write(*,"(a20,g15.3)") "chi_M = ", proj%chi_M 
    write(*,"(a20,g15.3)") "q_M = ", proj%q_M 
    write(*,"(a20,g15.3)") "q_polar = ", proj%q_polar 
    write(*,"(a20,g15.3)") "beta_M = ", proj%beta_M 
    write(*,"(a20,g15.3)") "R_q_polar = ", proj%R_q_polar 
    write(*,"(a20,g15.3)") "D = ", proj%D 
    
    return 

  end subroutine projection_print

  function same_projection(proj1,proj2) result(same_proj)

    implicit none 
    type(projection_class), intent(IN) :: proj1, proj2 
    logical  :: same_proj 
    real(dp) :: eps = 1.e-8_dp

    same_proj = .FALSE. 

    if (trim(proj1%method) .eq. trim(proj2%method) .and. &
        proj1%lambda .eq. proj2%lambda    .and. &
        proj1%phi    .eq. proj2%phi       .and. &
        proj1%alpha  .eq. proj2%alpha     .and. &
        proj1%x_e    .eq. proj2%x_e       .and. &
        proj1%y_n    .eq. proj2%y_n       .and. &
        dabs(proj1%f-proj2%f) .le. eps    .and. &
        dabs(proj1%a-proj2%a) .le. eps    .and. &
        dabs(proj1%e-proj2%e) .le. eps    .and. &
        dabs(proj1%R-proj2%R) .le. eps     ) then

            same_proj = .TRUE. 

    end if 

    return 

  end function same_projection

  function optimal_alpha(R,nx,ny,dx,dy) result(alpha)
    ! Given the dimensions and resolution of a projected grid,
    ! determine the optimal angle alpha (degrees) for the 
    ! oblique projection 

    implicit none 

    real(dp) :: R, dx, dy, alpha 
    integer  :: nx, ny 
    real(dp) :: val 

    val = 1.d0/R * sqrt( 1.d0/(2.d0*pi) * nx*ny*abs(dx*dy) )
    alpha = asin(val) * radians_to_degrees

    return

  end function optimal_alpha

  subroutine oblimap_projection(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)

    implicit none 

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    select case(trim(proj%method))

        case("oblique_sg_projection","oblique_sg_projection_snyder")
            ! These choices use the same forward projection routine 

            call oblique_sg_projection(lambda,phi,x_IM_P_prime,y_IM_P_prime,proj)

        case("oblique_laea_projection_snyder")

            call oblique_laea_projection_snyder(lambda,phi,x_IM_P_prime,y_IM_P_prime,proj)

        case("oblique_sg_projection_ellipsoid_snyder")

            call oblique_sg_projection_ellipsoid_snyder(lambda,phi,x_IM_P_prime,y_IM_P_prime,proj)

        case("oblique_laea_projection_ellipsoid_snyder")

            call oblique_laea_projection_ellipsoid_snyder(lambda,phi,x_IM_P_prime,y_IM_P_prime,proj)

        case DEFAULT 
            write(*,*) "oblimap_projection:: error: projection method not recognized: "// &
                       trim(proj%method)
            write(*,*) "Only the following options are possible: "
            write(*,*) "oblique_sg_projection"
            write(*,*) "oblique_sg_projection_snyder"
            write(*,*) "oblique_laea_projection_snyder"
            write(*,*) "oblique_sg_projection_ellipsoid_snyder"
            write(*,*) "oblique_laea_projection_ellipsoid_snyder"
            write(*,*) 

            stop 

    end select 

    return 

  end subroutine oblimap_projection

  subroutine oblimap_projection_inverse(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)

    implicit none 

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    select case(trim(proj%method))

        case("oblique_sg_projection")

            call inverse_oblique_sg_projection(x_IM_P_prime,y_IM_P_prime,lambda_P,phi_P,proj)

        case("oblique_sg_projection_snyder")

            call inverse_oblique_sg_projection_snyder(x_IM_P_prime,y_IM_P_prime,lambda_P,phi_P,proj)

        case("oblique_laea_projection_snyder")

            call inverse_oblique_laea_projection_snyder(x_IM_P_prime,y_IM_P_prime,lambda_P,phi_P,proj)

        case("oblique_sg_projection_ellipsoid_snyder")

            call inverse_oblique_sg_projection_ellipsoid_snyder(x_IM_P_prime,y_IM_P_prime,lambda_P,phi_P,proj)

        case("oblique_laea_projection_ellipsoid_snyder")

            call inverse_oblique_laea_projection_ellipsoid_snyder(x_IM_P_prime,y_IM_P_prime,lambda_P,phi_P,proj)

        case DEFAULT 
            write(*,*) "oblimap_projection:: error: projection method not recognized: "// &
                       trim(proj%method)
            write(*,*) "Only the following options are possible: "
            write(*,*) "oblique_sg_projection"
            write(*,*) "oblique_sg_projection_snyder"
            write(*,*) "oblique_laea_projection_snyder"
            write(*,*) "oblique_sg_projection_ellipsoid_snyder"
            write(*,*) "oblique_laea_projection_ellipsoid_snyder"
            write(*,*) 

            stop 

    end select 

    return 

  end subroutine oblimap_projection_inverse

  SUBROUTINE oblique_sg_projection(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate 
    ! system, with coordinates (x,y).
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: t_P_prime
    
    ! For North and South Pole: proj%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = degrees_to_radians * phi       
    lambda_P = degrees_to_radians * lambda
    
    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(proj%alpha_stereographic)) / &
                 (1._dp + COS(phi_P) * COS(proj%phi_M) * COS(lambda_P - proj%lambda_M) + SIN(phi_P) * SIN(proj%phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  proj%R * (COS(phi_P) * SIN(lambda_P - proj%lambda_M)) * t_P_prime
    y_IM_P_prime =  proj%R * (SIN(phi_P) * COS(proj%phi_M) - &
                     (COS(phi_P) * SIN(proj%phi_M)) * COS(lambda_P - proj%lambda_M)) * t_P_prime
    
    return

  END SUBROUTINE oblique_sg_projection

  SUBROUTINE inverse_oblique_sg_projection(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with an inverse oblique stereographic projection the 
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime
    REAL(dp)              :: y_3D_P_prime
    REAL(dp)              :: z_3D_P_prime
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P
    REAL(dp)              :: y_3D_P
    REAL(dp)              :: z_3D_P
    
    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = proj%R * COS(proj%alpha_stereographic) * COS(proj%lambda_M) * &
                     COS(proj%phi_M) - SIN(proj%lambda_M) * &
                     x_IM_P_prime - COS(proj%lambda_M) * SIN(proj%phi_M) * y_IM_P_prime
    y_3D_P_prime = proj%R * COS(proj%alpha_stereographic) * SIN(proj%lambda_M) * &
                     COS(proj%phi_M) + COS(proj%lambda_M) * &
                     x_IM_P_prime - SIN(proj%lambda_M) * SIN(proj%phi_M) * y_IM_P_prime
    z_3D_P_prime = proj%R * COS(proj%alpha_stereographic) * &
                     SIN(proj%phi_M)                                 + &
                     COS(proj%phi_M) * y_IM_P_prime
    
    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(proj%lambda_M) * COS(proj%phi_M) * x_3D_P_prime  +  SIN(proj%lambda_M) * &
         COS(proj%phi_M) * y_3D_P_prime  +  SIN(proj%phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * proj%R**2 + 2._dp * proj%R * a) / &
           (proj%R**2 + 2._dp * proj%R * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  proj%R * COS(proj%lambda_M) * COS(proj%phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  proj%R * SIN(proj%lambda_M) * COS(proj%phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  proj%R *                      SIN(proj%phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
     lambda_P = 180._dp + radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
     lambda_P =           radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 360._dp + radians_to_degrees * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
     lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
     lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
     lambda_P =   0._dp
    END IF

   ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
   IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
    phi_P = radians_to_degrees * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2)) 
   ELSE IF(z_3D_P >  0._dp) THEN
    phi_P =   90._dp
   ELSE IF(z_3D_P <  0._dp) THEN
    phi_P =  -90._dp
   END IF

   return

  END SUBROUTINE inverse_oblique_sg_projection

    SUBROUTINE inverse_oblique_sg_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection the 
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C     ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (20-18) on page 159 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)

    ! In case point P coincides with M (see condition at the first line of page  159 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = radians_to_degrees * proj%lambda_M
     phi_P    = radians_to_degrees * proj%phi_M

    ELSE 

     ! See equation (21-15) on page 159 Snyder (1987), because the denominator is always positive this ATAN doesn't 
     ! need a correction like note 2 on page ix in Snyder (1987):
     angle_C  = 2._dp * ATAN(rho / ((1._dp + COS(proj%alpha_stereographic)) * proj%earth_radius))
     
     ! See equation (20-14) on page 158 Snyder (1987):
     if (rho /= 0._dp) then 
         phi_P    = radians_to_degrees * &
         ( ASIN(COS(angle_C)*SIN(proj%phi_M) + ((y_IM_P_prime*SIN(angle_C)*COS(proj%phi_M)) / rho)) )
     else 
         phi_P    = 1._dp
     end if 

     ! See equation (20-15) on page 159 Snyder (1987):
     numerator   = x_IM_P_prime * SIN(angle_C)
     denumerator = rho * COS(proj%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(proj%phi_M) * SIN(angle_C)
     lambda_P    = radians_to_degrees * (proj%lambda_M + arctangens_quotient(numerator, denumerator))
    
    END IF 

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

  END SUBROUTINE inverse_oblique_sg_projection_snyder



  SUBROUTINE oblique_laea_projection_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection the 
    ! longitude-latitude coordinates which coincide with the GCM grid points to the rectangular IM 
    ! coordinate system, with coordinates (x,y).
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: t_P_prime
    
    ! For North and South Pole: proj%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = degrees_to_radians * phi       
    lambda_P = degrees_to_radians * lambda

    ! See equation (21-4) on page 185 of Snyder (1987):
    t_P_prime = SQRT(2._dp / (1._dp + COS(phi_P) * COS(proj%phi_M) * COS(lambda_P - proj%lambda_M) + &
                            SIN(phi_P) * SIN(proj%phi_M)))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010), page 185 of Snyder (1987):
    x_IM_P_prime =  proj%earth_radius * (COS(phi_P)*SIN(lambda_P - proj%lambda_M)) * t_P_prime
    y_IM_P_prime =  proj%earth_radius *  &
        (SIN(phi_P)*COS(proj%phi_M) - (COS(phi_P)*SIN(proj%phi_M)) * COS(lambda_P-proj%lambda_M))*t_P_prime
  END SUBROUTINE oblique_laea_projection_snyder


  SUBROUTINE inverse_oblique_laea_projection_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection 
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C              ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (20-18) on page 187 Snyder (1987):
    rho      = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)


    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN

     lambda_P = radians_to_degrees * proj%lambda_M
     phi_P    = radians_to_degrees * proj%phi_M

    ELSE 

     ! See equation (24-16) on page 187 Snyder (1987):
     angle_C  = 2._dp * ASIN(rho / (2._dp * proj%earth_radius))
     
     ! See equation (20-14) on page 186 Snyder (1987):
     if (rho /= 0._dp) then 
         phi_P    = radians_to_degrees* &
          ( ASIN(COS(angle_C)*SIN(proj%phi_M) + ((y_IM_P_prime * SIN(angle_C) * COS(proj%phi_M)) / rho)) )
     else 
         phi_P    = 1._dp
     end if 
     
     ! See equation (20-15) on page 186 Snyder (1987):
     numerator   = x_IM_P_prime * SIN(angle_C)
     denumerator = rho * COS(proj%phi_M) * COS(angle_C) - y_IM_P_prime * SIN(proj%phi_M) * SIN(angle_C)
     lambda_P    = radians_to_degrees * (proj%lambda_M + arctangens_quotient(numerator, denumerator))
    
    END IF 

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    
  END SUBROUTINE inverse_oblique_laea_projection_snyder


! ===== NEW ========
! obtained from source published in: http://www.geosci-model-dev-discuss.net/gmd-2016-124/

  SUBROUTINE oblique_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with Snyder's oblique stereographic projection for the ellipsoid
    ! the the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y).
    !
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395

    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda        ! in degrees
    REAL(dp), INTENT(IN)  :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT) :: y_IM_P_prime  ! in meter

    ! Local variables:
    REAL(dp)              :: phi_P         ! in radians,  phi    in Snyder (1987)
    REAL(dp)              :: lambda_P      ! in radians,  lambda in Snyder (1987)
    REAL(dp)              :: chi_P         ! in radians,  chi    in Snyder (1987)
    REAL(dp)              :: A

    ! For North and South Pole: proj%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IF(proj%name .eq. "polar_stereographic") THEN
     ! The polar case is excepted from the oblique formula's, see page 161 Snyder (1987)

     ! Output: x_IM_P_prime, y_IM_P_prime
     CALL polar_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ELSE
     ! The oblique case, see page 160 Snyder (1987)

     ! Convert longitude-latitude coordinates to radians:
     phi_P    = degrees_to_radians * phi
     lambda_P = degrees_to_radians * lambda

     ! See equations (3-1a) and (21-27) on page 160 in Snyder (1987):
     chi_P = 2._dp * ATAN(SQRT(((1._dp +       SIN(phi_P)) / (1._dp -       SIN(phi_P))) * &
                    ((1._dp - proj%e * SIN(phi_P)) / (1._dp + proj%e * SIN(phi_P)))**(proj%e))) - 0.5_dp * pi
     A     = proj%akm / (COS(proj%chi_M) * (1._dp + SIN(proj%chi_M) * SIN(chi_P) +  &
                         COS(proj%chi_M) * COS(chi_P) * COS(lambda_P - proj%lambda_M)))


     ! See equations (21-24) and (21-25) on page 160 in Snyder (1987):
     x_IM_P_prime =  A * COS(chi_P) * SIN(lambda_P - proj%lambda_M)
     y_IM_P_prime =  A * (COS(proj%chi_M) * SIN(chi_P) - SIN(proj%chi_M) * COS(chi_P) * COS(lambda_P - proj%lambda_M))
    END IF
  END SUBROUTINE oblique_sg_projection_ellipsoid_snyder



  SUBROUTINE polar_sg_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with Snyder's polar stereographic projection for the ellipsoid
    ! the the longitude-latitude coordinates which coincide with the GCM grid points to
    ! the rectangular IM coordinate system, with coordinates (x,y). See Snyder (1987) p. 161.
    !
    ! The examples of Snyder (1987) at p. 314-315 with the international ellipsoid are used to
    ! validate this forward SG projection on the ellipsoid for the polar aspect.
    !
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    
    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda        ! in degrees
    REAL(dp), INTENT(IN)  :: phi           ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT) :: y_IM_P_prime  ! in meter

    ! Local variables:
    REAL(dp)              :: phi_P         ! in radians
    REAL(dp)              :: lambda_P      ! in radians
    REAL(dp)              :: phi_C         ! in radians,  phi_c  in Snyder (1987), the standard parallel
    REAL(dp)              :: t_P           ! 
    REAL(dp)              :: t_C           ! 
    REAL(dp)              :: m_C           ! 
    REAL(dp)              :: rho           ! 
    REAL(dp)              :: pf            ! polar factor: -1.0 for SP and +1.0 for NP
    REAL(dp)              :: k0            ! 

    !IF(proj%phi_M == -90.0_dp * degrees_to_radians) THEN
    if(proj%phi_M .lt. 0.0_dp) then 
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    phi_C    = (degrees_to_radians * 90.0_dp - proj%alpha_stereographic) * pf

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = degrees_to_radians * phi
    lambda_P = degrees_to_radians * lambda

    IF(proj%alpha_stereographic == 0.0_dp) THEN
     ! This variant is only considered for the case that proj%alpha_stereographic = 0, i.e. k0 = 1, for which the else-option does not work).
     ! The POLAR ASPECT WITH KNOWN k0 case:
    !t_P = TAN(pi / 4.0_dp - phi_P * pf / 2.0_dp) / ( (1.0_dp - proj%e * SIN(phi_P * pf)) / (1.0_dp + proj%e * SIN(phi_P * pf)) )**(proj%e / 2.0_dp)                      ! (15-9)  on page 161 in Snyder (1987)
     t_P = ( ((1.0_dp - SIN(phi_P * pf)) / (1.0_dp + SIN(phi_P * pf))) *  &
             ((1.0_dp + proj%e * SIN(phi_P * pf)) / (1.0_dp - proj%e * SIN(phi_P * pf)))**proj%e )**0.5_dp    ! (15-9a) on page 161 in Snyder (1987)
     k0  = 0.5_dp * (1.0_dp + COS(proj%alpha_stereographic))                                                                                                        ! in fact it will be always 1 because it is only used for proj%alpha_stereographic = 0
     rho = 2.0_dp * proj%a * k0 * t_P / ( (1.0_dp + proj%e)**(1.0_dp + proj%e) * (1.0_dp - proj%e)**(1.0_dp - proj%e) )**0.5_dp                                                 ! (21-33) on page 161 in Snyder (1987)
    ELSE
     ! The POLAR ASPECT WITH KNOWN STANDARD PARALLEL NOT AT POLE case:
    !t_P   = TAN(pi / 4.0_dp - phi_P * pf / 2.0_dp) / ( (1.0_dp - proj%e * SIN(phi_P * pf)) / (1.0_dp + proj%e * SIN(phi_P * pf)) )**(proj%e / 2.0_dp)                    ! (15-9)  on page 161 in Snyder (1987)
     t_P   = ( ((1.0_dp - SIN(phi_P * pf)) / (1.0_dp + SIN(phi_P * pf))) * &
               ((1.0_dp + proj%e * SIN(phi_P * pf)) / (1.0_dp - proj%e * SIN(phi_P * pf)))**proj%e )**0.5_dp  ! (15-9a) on page 161 in Snyder (1987)
    !t_C   = TAN(pi / 4.0_dp - phi_C * pf / 2.0_dp) / ( (1.0_dp - proj%e * SIN(phi_C * pf)) / (1.0_dp + proj%e * SIN(phi_C * pf)) )**(proj%e / 2.0_dp)                    ! (15-9)  on page 161 in Snyder (1987)
     t_C   = ( ((1.0_dp - SIN(phi_C * pf)) / (1.0_dp + SIN(phi_C * pf))) * &
               ((1.0_dp + proj%e * SIN(phi_C * pf)) / (1.0_dp - proj%e * SIN(phi_C * pf)))**proj%e )**0.5_dp  ! (15-9a) on page 161 in Snyder (1987)
     m_C   = COS(phi_C * pf) / (1.0_dp - proj%e**2 * SIN(phi_C * pf)**2)**0.5_dp                                                                                    ! (14-15) on page 160 in Snyder (1987)
     rho   = proj%a * m_C * t_P / t_C                                                                                                                               ! (21-34) on page 161 in Snyder (1987)
    END IF

    x_IM_P_prime   =   rho * SIN(lambda_P * pf - proj%lambda_M * pf) * pf                                                                                           ! (21-30) on page 161 in Snyder (1987)
    y_IM_P_prime   = - rho * COS(lambda_P * pf - proj%lambda_M * pf) * pf                                                                                           ! (21-31) on page 161 in Snyder (1987)
  END SUBROUTINE polar_sg_projection_ellipsoid_snyder

           
  SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with Snyder's inverse oblique stereographic projection for the ellipsoid
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    
    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C       ! in radians
    REAL(dp)              :: chi_P         ! in radians,  chi in Snyder (1987)
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator

    IF(proj%name .eq. "polar_stereographic") THEN
     ! The inverse polar case is excepted from the inverse oblique formula's, see page 161 Snyder (1987)

     ! Output: lambda_P, phi_P
     call inverse_polar_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ELSE
     ! The inverse oblique case, see page 160 Snyder (1987)

     ! See equation (20-18) on page 162 Snyder (1987):
     rho     = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)

     ! In case point P coincides with M (see condition at the first line of page  162 Snyder (1987):
     IF(rho == 0._dp) THEN

      lambda_P = radians_to_degrees * proj%lambda_M
      phi_P    = radians_to_degrees * proj%phi_M

     ELSE 

      ! See equation (21-38) on page 162 Snyder (1987):
      angle_C = 2._dp * ATAN(rho * COS(proj%chi_M) / proj%akm)

      ! See equations (21-37) on page 161 in Snyder (1987):
      chi_P   = ASIN(COS(angle_C) * SIN(proj%chi_M) + y_IM_P_prime * SIN(angle_C) * COS(proj%chi_M) / rho)

      ! See equation (3-5) on page 162 instead of equation (3-4) on page 161 Snyder (1987):
      phi_P = radians_to_degrees * (chi_P + &
              (proj%e**2 / 2._dp + 5._dp * proj%e**4 / 24._dp +          &
                 proj%e**6 / 12._dp  +   13._dp * proj%e**8 /    360._dp) * SIN(2._dp * chi_P) + &
              (                 7._dp * proj%e**4 / 48._dp +             &
                 29._dp * proj%e**6 / 240._dp +  811._dp * proj%e**8 /  11520._dp) * SIN(4._dp * chi_P) + &
              (   7._dp * proj%e**6 / 120._dp +   81._dp * proj%e**8 /   1120._dp) * SIN(6._dp * chi_P) + &
              (        4279._dp * proj%e**8 / 161280._dp) * SIN(8._dp * chi_P))

      ! See equation (21-36) on page 161 Snyder (1987):
      numerator   = x_IM_P_prime * SIN(angle_C)
      denumerator = rho * COS(proj%chi_M) * COS(angle_C) - y_IM_P_prime * SIN(proj%chi_M) * SIN(angle_C)
      lambda_P    = radians_to_degrees * (proj%lambda_M + arctangens_quotient(numerator, denumerator))

     END IF 

     ! Our choice is to return lambda in the 0-360 degree range:
     IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp

    END IF

  END SUBROUTINE inverse_oblique_sg_projection_ellipsoid_snyder


  SUBROUTINE inverse_polar_sg_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with Snyder's inverse polar stereographic projection for the ellipsoid
    ! the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees. See Snyder (1987) p. 162.
    !
    ! The examples of Snyder (1987) at p. 317-318 with the international ellipsoid are used to
    ! validate this inverse SG projection on the ellipsoid for the polar aspect.
    !
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    
    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees,  lambda in Snyder (1987)
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees,  phi    in Snyder (1987)

    ! Local variables:
    REAL(dp)              :: chi_P         ! in radians,  chi    in Snyder (1987)
    REAL(dp)              :: phi_C         ! in radians,  phi_c  in Snyder (1987), the standard parallel
    REAL(dp)              :: t_P
    REAL(dp)              :: t_C
    REAL(dp)              :: m_C
    REAL(dp)              :: rho
    REAL(dp)              :: pf            ! polar factor: -1.0 for SP and +1.0 for NP
    REAL(dp)              :: k0

    ! IF(proj%phi_M == -90.0_dp * degrees_to_radians) THEN
    if (proj%phi_M .lt. 0.0_dp) then 
     pf = -1.0_dp                                       ! The polar factor for the SP
    ELSE
     pf =  1.0_dp                                       ! The polar factor for the NP
    END IF

    phi_C = (degrees_to_radians * 90.0_dp - proj%alpha_stereographic) * pf

    rho   = SQRT(x_IM_P_prime**2 + y_IM_P_prime**2)                                                                                                           ! (20-18) on page 162 in Snyder (1987)

    IF(proj%alpha_stereographic == 0.0_dp) THEN
     ! This variant is only considered for the case that proj%alpha_stereographic = 0, i.e. k0 = 1, for which the else-option does not work).
     ! The POLAR ASPECT WITH KNOWN k0 case:
     k0  = 0.5_dp * (1.0_dp + COS(proj%alpha_stereographic))                                                                                                    ! in fact it will be always 1 because it is only used for proj%alpha_stereographic = 0
     t_P = rho * ( (1.0_dp + proj%e)**(1.0_dp + proj%e) * (1.0_dp - proj%e)**(1.0_dp - proj%e) )**0.5_dp / (2.0_dp * proj%a * k0)                                            ! (21-39) on page 162 in Snyder (1987)
    ELSE
     ! The POLAR ASPECT WITH KNOWN STANDARD PARALLEL NOT AT POLE case:
    !t_C = TAN(pi / 4.0_dp - phi_C * pf / 2.0_dp) / ( (1.0_dp - proj%e * SIN(phi_C * pf)) / (1.0_dp + proj%e * SIN(phi_C * pf)) )**(proj%e / 2.0_dp)                   ! (15-9)  on page 161 in Snyder (1987)
     t_C = ( ((1.0_dp - SIN(phi_C * pf)) / (1.0_dp + SIN(phi_C * pf))) *  &
                ((1.0_dp + proj%e * SIN(phi_C * pf)) / (1.0_dp - proj%e * SIN(phi_C * pf)))**proj%e )**0.5_dp ! (15-9a) on page 161 in Snyder (1987)
     m_C = COS(phi_C * pf) / (1.0_dp - proj%e**2 * SIN(phi_C * pf)**2)**0.5_dp                                                                                   ! (14-15) on page 160 in Snyder (1987)
     t_P = rho * t_C / (proj%a * m_C)                                                                                                                            ! (21-40) on page 162 in Snyder (1987)
    END IF

    ! Note: Eventually replace ATAN with arctangens_quotient(numerator = rho * t_C, denumerator = proj%a * m_C) :
    chi_P = pi / 2.0_dp - 2.0_dp * ATAN(t_P)                                                                                                                ! (7-13)  on page 162 in Snyder (1987)

    ! See equation (3-5) on page 162 instead of equation (7-9) on page 162 Snyder (1987):
    phi_P = radians_to_degrees * (chi_P + &
            (proj%e**2 / 2._dp + 5._dp * proj%e**4 / 24._dp +      &
                proj%e**6 / 12._dp  +   13._dp * proj%e**8 /    360._dp) * SIN(2._dp * chi_P) + &
            (                 7._dp * proj%e**4 / 48._dp +         & 
                29._dp * proj%e**6 / 240._dp +  811._dp * proj%e**8 /  11520._dp) * SIN(4._dp * chi_P) + &
            (        7._dp * proj%e**6 / 120._dp +   81._dp * proj%e**8 /   1120._dp) * SIN(6._dp * chi_P) + &
            (          4279._dp * proj%e**8 / 161280._dp) * SIN(8._dp * chi_P)) *pf                 ! (3-5)   on page 162 in Snyder (1987)
    lambda_P = radians_to_degrees * (proj%lambda_M * pf + arctangens_quotient(x_IM_P_prime * pf, - y_IM_P_prime * pf)) *pf                                     ! (20-16) on page 162 in Snyder (1987)

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp)    lambda_P = lambda_P + 360._dp

  END SUBROUTINE inverse_polar_sg_projection_ellipsoid_snyder


! ===== END NEW ========

  SUBROUTINE oblique_laea_projection_ellipsoid_snyder(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with Snyder's oblique Lambert azimuthal equal-area projection for 
    ! the ellipsoid the longitude-latitude coordinates which coincide with the GCM grid points to 
    ! the rectangular IM coordinate system, with coordinates (x,y).
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    
    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: lambda
    REAL(dp), INTENT(IN)  :: phi

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_P_prime
    REAL(dp), INTENT(OUT) :: y_IM_P_prime

    ! Local variables:
    REAL(dp)              :: phi_P
    REAL(dp)              :: lambda_P
    REAL(dp)              :: q_P      ! q in Snyder (1987)
    REAL(dp)              :: beta
    REAL(dp)              :: B
    
    ! For North and South Pole: proj%lambda_M = 0._dp, to generate the correct IM coordinate 
    ! system, see the oblimap_configuration_module and see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = degrees_to_radians * phi       
    lambda_P = degrees_to_radians * lambda
    
    ! See equation (3-12) on page 187 in Snyder (1987):
    q_P = (1._dp - proj%e**2) * ((SIN(phi_P) / (1._dp - (proj%e * SIN(phi_P))**2)) - &
        (1._dp / (2._dp*proj%e)) * LOG((1._dp - proj%e*SIN(phi_P)) / (1._dp + proj%e * SIN(phi_P)))) 
    ! See equation (3-11) on page 187 in Snyder (1987):
    beta = ASIN(q_P / proj%q_polar)
    ! See equation (24-19) on page 187 in Snyder (1987):
    B = proj%R_q_polar* &
        SQRT(2._dp / (1._dp + SIN(proj%beta_M)*SIN(beta) + COS(proj%beta_M)*COS(beta) * COS(lambda_P-proj%lambda_M)))

    ! See equation (24-17) and (24-18) on page 187 in Snyder (1987):
    x_IM_P_prime = B * proj%D * COS(beta) * SIN(lambda_P - proj%lambda_M)
    y_IM_P_prime = (B / proj%D) * &
        (COS(proj%beta_M)*SIN(beta) - SIN(proj%beta_M)*COS(beta) * COS(lambda_P-proj%lambda_M))
  END SUBROUTINE oblique_laea_projection_ellipsoid_snyder



  SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder(x_IM_P_prime, y_IM_P_prime, lambda_P, phi_P, proj)
    ! This subroutine projects with Snyder's inverse oblique Lambert azimuthal equal-area projection for 
    ! the ellipsoid the (x,y) coordinates which coincide with the IM grid points to the longitude-latitude 
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    ! 
    ! For more information about M, proj%alpha_stereographic, the center of projection and the used 
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    ! and
    !  Snyder (1987), map projections: A working manual, http://pubs.er.usgs.gov/usgspubs/pp/pp1395
    
    IMPLICIT NONE

    type(projection_class), INTENT(IN) :: proj 

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime
    REAL(dp), INTENT(IN)  :: y_IM_P_prime

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P
    REAL(dp), INTENT(OUT) :: phi_P

    ! Local variables:
    REAL(dp)              :: rho
    REAL(dp)              :: angle_C           ! In radians
    REAL(dp)              :: beta              ! In radians
    REAL(dp)              :: numerator
    REAL(dp)              :: denumerator
    
    ! See equation (24-28) on page 189 Snyder (1987):
    rho      = SQRT((x_IM_P_prime / proj%D)**2 + (y_IM_P_prime / proj%D)**2)

    ! In case point P coincides with M (see the condition down equation (20-14) on page 186 Snyder (1987):
    IF(rho == 0._dp) THEN
     lambda_P = radians_to_degrees * proj%lambda_M
     phi_P    = radians_to_degrees * proj%phi_M

    ELSE 

     ! See equation (24-29) on page 189 Snyder (1987):
     angle_C  = 2._dp * ASIN(rho / (2._dp * proj%R_q_polar))
     
     ! See equation (24-30) on page 189 Snyder (1987):
     beta = ASIN(COS(angle_C)*SIN(proj%beta_M) + (proj%D*y_IM_P_prime*SIN(angle_C)*COS(proj%beta_M) / rho)) 

     ! See equation (3-18) on page 189 instead of equation (3-16) on page 188 Snyder (1987):
     phi_P = radians_to_degrees * (beta + &
             (proj%e**2 / 3._dp + 31._dp * proj%e**4 / 180._dp + 517._dp * proj%e**6 /  5040._dp) * SIN(2._dp * beta) + &
             (                 23._dp * proj%e**4 / 360._dp + 251._dp * proj%e**6 /  3780._dp) * SIN(4._dp * beta) + &
             (                                             761._dp * proj%e**6 / 45360._dp) * SIN(6._dp * beta)) 
     
     ! See equation (20-26) on page 188 Snyder (1987):
     numerator   = x_IM_P_prime * SIN(angle_C)
     denumerator = proj%D * rho * COS(proj%beta_M) * COS(angle_C) - proj%D**2 * y_IM_P_prime * SIN(proj%beta_M) * SIN(angle_C)
     lambda_P    = radians_to_degrees * (proj%lambda_M + arctangens_quotient(numerator, denumerator))
    
    END IF 

    ! Our choice is to return lambda in the 0-360 degree range:
    IF(lambda_P < 0._dp) lambda_P = lambda_P + 360._dp
    

  END SUBROUTINE inverse_oblique_laea_projection_ellipsoid_snyder

  FUNCTION arctangens_quotient(numerator, denumerator) RESULT(angle)

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: numerator  
    REAL(dp), INTENT(IN)  :: denumerator

    ! Result variables:
    REAL(dp)              :: angle       ! In radians

    ! Local variables:
    REAL(dp)              :: quadrant_correction
    
    ! See note 2 on page ix in Snyder (1987), to distinguish between the quadrants:
    quadrant_correction = 0._dp
    IF(denumerator <  0._dp) quadrant_correction =          pi
    IF(denumerator == 0._dp) quadrant_correction = 0.5_dp * pi
    IF(numerator   <  0._dp) quadrant_correction = - quadrant_correction
    
    IF(denumerator == 0._dp) THEN
     IF(numerator == 0._dp) THEN
      ! The angle is indetermined, usually zero is taken:
      angle = 0._dp
     ELSE
      angle = 0.5_dp * pi * (numerator / ABS(numerator))
     END IF
    ELSE
     angle = ATAN(numerator / denumerator) + quadrant_correction
    END IF

    return 

  END FUNCTION arctangens_quotient



  SUBROUTINE rotation_projection(x_IM, y_IM, x_IM_prime, y_IM_prime, theta)
    ! This subroutine transforms the 2D coordinates which are projected with a rotation from
    ! a rectangular IM coordinate system (x,y) to another rectangular IM prime coordinate
    ! system (x',y')

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM
    REAL(dp), INTENT(IN)  :: y_IM
    REAL(dp), INTENT(IN)  :: theta 
    
    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM_prime
    REAL(dp), INTENT(OUT) :: y_IM_prime

    ! See for a derivation: http://www.youtube.com/watch?v=h11ljFJeaLo
    x_IM_prime =   x_IM * COS(theta) + y_IM * SIN(theta)
    y_IM_prime = - x_IM * SIN(theta) + y_IM * COS(theta)

    return

  END SUBROUTINE rotation_projection



  SUBROUTINE inverse_rotation_projection(x_IM_prime, y_IM_prime, x_IM, y_IM, theta)
    ! This subroutine transforms the 2D coordinates which are projected with a inverse rotation
    ! from a rectangular IM prime coordinate system (x',y') to another rectangular IM coordinate
    ! system (x,y)

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_prime
    REAL(dp), INTENT(IN)  :: y_IM_prime
    REAL(dp), INTENT(IN)  :: theta 

    ! Output variables:
    REAL(dp), INTENT(OUT) :: x_IM
    REAL(dp), INTENT(OUT) :: y_IM

    ! This equations are derived by taking a linear combination of the rotation projection equations.
    x_IM = x_IM_prime * COS(theta) - y_IM_prime * SIN(theta)
    y_IM = x_IM_prime * SIN(theta) + y_IM_prime * COS(theta)

    return 

  END SUBROUTINE inverse_rotation_projection

END MODULE oblimap_projection_module
