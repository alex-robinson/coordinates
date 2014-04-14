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

!! ########################################### !!
!! NOW MODIFIED TO BE MORE ENCAPSULATED        !!
!! ajr, 2013-09-17                             !!
!! ########################################### !!

MODULE oblimap_projection_module

    use planet 

    implicit none 

    !! real(dp) definition
    integer, parameter :: dp    = kind(1.0d0)

    ! Mathematical constants
    real(dp), parameter  :: pi  = 2._dp*acos(0._dp)
    real(dp), parameter  :: degrees_to_radians = pi / 180._dp  ! Conversion factor between radians and degrees
    real(dp), parameter  :: radians_to_degrees = 180._dp / pi  ! Conversion factor between degrees and radians
        
    type projection_class

        character(len=256) :: name 

        ! Projection parameters
        real(dp) :: lambda, lambda_M 
        real(dp) :: phi, phi_M 
        real(dp) :: alpha, alpha_stereographic
        real(dp) :: x_e, y_n 

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
    public :: pi, degrees_to_radians, radians_to_degrees
    public :: projection_class, projection_init, same_projection
    public :: oblique_sg_projection, inverse_oblique_sg_projection

CONTAINS

  SUBROUTINE projection_init(proj,name,planet,lambda,phi,alpha,x_e,y_n)
    
    IMPLICIT NONE

    type(projection_class) :: proj 
    character(len=*)       :: name
    type(planet_class)     :: planet 
    real(dp), optional     :: lambda, phi, alpha 
    real(dp), optional     :: x_e, y_n 

    proj%name = trim(name) 

    proj%f   = planet%f 
    proj%a   = planet%a 
    proj%e   = planet%e 
    proj%R   = planet%R 
    proj%earth_radius = planet%R 

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

    if ( present(alpha) ) then    

        if ( alpha .gt. 0.0_dp ) then 

            proj%alpha  = alpha 

        else
            write(*,*) "oblimap_projection_module::"
            write(*,*) "    Optimal alpha calculation not yet developed."
            write(*,*) "    Instead, please specify a positive value for alpha."
            write(*,*)
            stop

!             ! Check for optimal alpha
!             test1 = grid%nx*grid%ny*(grid%dx*grid%xy_conv)*(grid%dy*grid%xy_conv)
!             test2 = 2.0_dp * pi * grid%P%R**2
!             if (test1 .lt. test2) then  ! Projected grid is smaller than 1 hemisphere
!                 alpha = asin( (1.0_dp / grid%P%R) sqrt((1/(2.0_dp*pi)*test1)) )
!             else                        ! Projected grid is equal to 1 hemisphere
!                 gamma = 45.0_dp * degrees_to_radians
!                 alpha = 2.0_dp * atan( sqrt(1.0_dp/(2.0_dp*pi))*tan(gamma) )  ! alpha = 43.5 deg
!                 alpha = alpha * radians_to_degrees
!             end if

        end if 

    end if

    proj%phi_M               = degrees_to_radians * proj%phi
    proj%lambda_M            = degrees_to_radians * proj%lambda
    proj%alpha_stereographic = degrees_to_radians * proj%alpha

    ! See equations (14-15) and (21-27) on page 160 in Snyder (1987), akm corresponds with 2a*k0*m_1 in Snyder:
    proj%am     = proj%a * (COS(proj%phi_M) / DSQRT(1._dp - (proj%e * SIN(proj%phi_M))**2))
    proj%akm    = (1._dp + COS(proj%alpha_stereographic)) * proj%am
    ! See equations (3-1a) on page 160 in Snyder (1987),  chi_M corresponds with chi_1 in Snyder:
    proj%chi_M  = 0.0_dp 
    if (dabs(proj%phi) .lt. 90.0_dp) &
      proj%chi_M  = 2._dp * ATAN(DSQRT(((1._dp +       SIN(proj%phi_M)) / (1._dp -       SIN(proj%phi_M))) * &
                           ((1._dp - proj%e * SIN(proj%phi_M)) / (1._dp + proj%e * SIN(proj%phi_M)))**(proj%e))) - 0.5_dp * pi

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

    return

  END SUBROUTINE projection_init 

  function same_projection(proj1,proj2) result(same_proj)

    implicit none 
    type(projection_class), intent(IN) :: proj1, proj2 
    logical  :: same_proj 
    real(dp) :: eps = 1.e-8_dp

    same_proj = .FALSE. 

    if (proj1%lambda .eq. proj2%lambda    .and. &
        proj1%phi    .eq. proj2%phi       .and. &
        proj1%alpha  .eq. proj2%alpha     .and. &
        proj1%x_e    .eq. proj2%x_e       .and. &
        proj1%y_n    .eq. proj2%y_n       .and. &
        dabs(proj1%f-proj2%f) .le. eps    .and. &
        dabs(proj1%a-proj2%a) .le. eps    .and. &
        dabs(proj1%e-proj2%e) .le. eps    .and. &
        dabs(proj1%R-proj2%R) .le. eps         ) then

            same_proj = .TRUE. 

    end if 

    return 

    end function same_projection

  SUBROUTINE oblique_sg_projection(lambda, phi, x_IM_P_prime, y_IM_P_prime, proj)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate 
    ! system, with coordinates (x,y).
    ! 
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
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
    ! For more information about M, C%alpha_stereographic, the center of projection and the used 
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


  FUNCTION arctanges_quotient(numerator, denumerator) RESULT(angle)

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

  END FUNCTION arctanges_quotient



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
