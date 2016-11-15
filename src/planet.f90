module planet

    use coord_constants
    use geodesic

    implicit none 

    type planet_class
        character(len=256) :: name 
        logical :: is_sphere
        real(dp) :: a, e, f
        real(dp) :: R      ! in case of sphere (equivalent to a=R, f=1)
    end type 

    interface cartesian_distance
        module procedure cartesian_distance_float
        module procedure cartesian_distance_dble 
    end interface 

    interface planet_distance
        module procedure planet_distance_float
        module procedure planet_distance_dble 
    end interface 
    
    private
    public :: planet_class, planet_init, planet_print 
    public :: planet_distance, cartesian_distance, weighted_ave_shepard, weighted_ave
    public :: planet_area, cartesian_area 
    public :: spherical_distance, quadrant_latlon, quadrant_cartesian

contains

    subroutine planet_init(now,name)

        implicit none 

        type(planet_class) :: now
        character(len=*) name  

        now%name = trim(name)
        
        now%a         = 1.0_dp 
        now%f         = 1.0_dp 
        now%e         = 0.0_dp
        now%is_sphere = .TRUE. 

        select case(trim(now%name))
            case("WGS84")
                now%a = 6378137.0_dp             ! Equatorial ellipsoid radius,   a in Snyder, WGS84
                now%f = 1.0_dp/298.257223563_dp  ! Flattening of the ellipsoid
                now%R = 6.371221E6_dp            ! Radius of round Earth
                now%e = 0.081819191_dp           ! Eccentricity of the ellipsoid, e in Snyder, WGS84, see Snyder p. 13
                now%is_sphere = .FALSE.         

            case("Spherical Earth")
                now%R = 6.371221E6_dp            ! Radius of round Earth
                now%a = now%R                    ! a is the same as the radius
                now%f = 1.0_dp                   ! flattening is zero 
                now%e = 0.0_dp                   ! Elliptical parameter is zero 

            case DEFAULT
                
                write(*,"(a,a)") "planet:: planet_init: ", &
                                       "error: planet type not yet defined: "//trim(now%name)
                write(*,"(a)") "  Available choices are: 'WGS84', 'Spherical Earth' "
                write(*,*)
                stop 

        end select

        !call planet_print(now) 

        return

    end subroutine planet_init 

    subroutine planet_print(now)

        implicit none 

        type(planet_class), intent(IN) :: now 

        write(*,*) "== Planet parameters =="
        write(*,*) "a = ", now%a 
        write(*,*) "f = ", now%f 
        write(*,*) "e = ", now%e 
        write(*,*) "R = ", now%R
        write(*,*) "is sphere? ", now%is_sphere 

        return

    end subroutine planet_print 

    function planet_distance_float(a,f,lon1,lat1,lon2,lat2) result(dist)

        implicit none 

        real(dp), intent(IN) :: a, f   ! Most likely from planet object (with dp precision)
        real(4),  intent(IN) :: lon1, lat1, lon2, lat2
        real(4) :: dist
        real(8) :: s12
        real(8) :: azi1, azi2, a12, m12, MM12, MM21, SS12
        integer omask
        
        call invers(a, f, dble(lat1), dble(lon1), dble(lat2), dble(lon2), &
                    s12, azi1, azi2, omask, a12, m12, MM12, MM21, SS12)

        dist = real(s12)    ! in meters
!         dist = (/ s12, azi1 /)    ! in meters, degrees

        return

    end function planet_distance_float

    function planet_distance_dble(a,f,lon1,lat1,lon2,lat2) result(dist)

        implicit none 

        real(8), intent(IN) :: a, f, lon1, lat1, lon2, lat2
        real(8) :: dist
        real(8) :: s12
        real(8) :: azi1, azi2, a12, m12, MM12, MM21, SS12
        integer omask
        
        call invers(a, f, lat1, lon1, lat2, lon2, &
                    s12, azi1, azi2, omask, a12, m12, MM12, MM21, SS12)

        dist = s12    ! in meters
!         dist = (/ s12, azi1 /)    ! in meters, degrees

        return

    end function planet_distance_dble

    elemental function cartesian_distance_float(x1,y1,x2,y2) result(dist)

        implicit none 

        real(4), intent(IN) :: x1, y1, x2, y2
        real(4) :: dist 

        dist = sqrt( (y2-y1)**2 + (x2-x1)**2 )

        return

    end function cartesian_distance_float  

    elemental function cartesian_distance_dble(x1,y1,x2,y2) result(dist)

        implicit none 

        real(8), intent(IN) :: x1, y1, x2, y2
        real(8) :: dist 

        dist = dsqrt( (y2-y1)**2 + (x2-x1)**2 )

        return

    end function cartesian_distance_dble 

    function spherical_distance(a,f,lon1,lat1,lon2,lat2)
        ! Based on python code here:
        ! http://www.johndcook.com/python_longitude_latitude.html

        implicit none 

        real(dp) :: a, f, lon1, lat1, lon2, lat2
        real(dp) :: spherical_distance
        real(dp) :: phi1, phi2, theta1, theta2
        real(dp) :: tmp

        ! phi = 90 - latitude
        phi1 = (90.0_dp - lat1)*degrees_to_radians
        phi2 = (90.0_dp - lat2)*degrees_to_radians
            
        ! theta = longitude
        theta1 = lon1*degrees_to_radians
        theta2 = lon2*degrees_to_radians
            
        ! Compute spherical distance from spherical coordinates.
            
        ! For two locations in spherical coordinates 
        ! (1, theta, phi) and (1, theta, phi)
        ! cosine( arc length ) = 
        !    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
        ! distance = rho * arc length
        
        tmp =  sin(phi1)*sin(phi2)*cos(theta1-theta2) + cos(phi1)*cos(phi2) 
        spherical_distance = acos( tmp )*a

        return

    end function spherical_distance

    function planet_area(a,f,lons,lats)

        implicit none 

        real(dp) :: a, f, lons(:), lats(:)
        real(dp) :: planet_area
        real(dp) :: AA, PP 

        if (size(lons) .ne. size(lats)) then 
            write(*,"(a,a)") "planet:: planet_area: ", &
                             "error: lons and lats must be the same length."
            write(*,*)       "nlon, nlat: ",size(lons),size(lats)
            stop 
        end if 

        if (size(lons) .le. 2) then 
            write(*,"(a,a)") "planet:: planet_area: ", &
                             "error: more than 2 points must be given to calculate area."
            stop 
        end if 

        call area(a, f, lats, lons, size(lons), AA, PP)

        ! Return the area of the lon/lat region
        planet_area = dabs(AA) 

        return 
    
    end function planet_area 

    function cartesian_area(x,y)
        
        implicit none 

        real(dp) :: x(:), y(:)
        real(dp) :: cartesian_area 
        real(dp), allocatable :: tmp(:)
        integer  :: i, n 

        if (size(x) .ne. size(y)) then 
            write(*,"(a,a)") "planet:: cartesian_area: ", &
                             "error: x and y must be the same length."
            write(*,*)       "nx, ny: ",size(x),size(y)
            stop 
        end if 

        if (size(x) .le. 2) then 
            write(*,"(a,a)") "planet:: cartesian_area: ", &
                             "error: more than 2 points must be given to calculate area."
            stop 
        end if 

        n = size(x)
        allocate(tmp(n))

        ! Calculate the area of an irregular polygon
        do i = 1, n 
            if (i .lt. n) then 
                tmp(i) = x(i)*y(i+1)-y(i)*x(i+1)
            else
                tmp(i) = x(i)*y(1)-y(i)*x(1)
            end if 
        end do  


        cartesian_area = dabs(sum(tmp)) / 2.0_dp 

        return 
    
    end function cartesian_area 

    elemental function shepard_weight(dist,shepard_exponent)
        implicit none 

        real(dp), intent(IN) :: dist 
        real(dp), intent(IN), optional :: shepard_exponent
        real(dp) :: shep_e
        real(dp) :: shepard_weight

        shep_e = 2.0_dp 
        if (present(shepard_exponent)) shep_e = shepard_exponent

        shepard_weight = 1.0_dp / (dist**shep_e)

        return 
    end function shepard_weight 


    function weighted_ave(var,weight,mask)

        implicit none

        real(dp) :: var(:), weight(:)
        logical,  optional :: mask(:) 
        real(dp) :: weighted_ave
        real(dp) :: numerator, denominator
        integer :: i, n

        logical :: msk(size(var)) 

        msk = .TRUE. 
        if (present(mask)) msk = mask 

        numerator            = sum( var*weight, mask=msk ) 
        denominator          = sum( weight,     mask=msk ) 

        if (denominator .gt. 0.0_dp) then
            weighted_ave = numerator / denominator
!             weighted_ave = min(dabs(weighted_ave),1d20)
        else 
            write(*,*) "weighted_ave:: No weights provided."
            weighted_ave = -9999.0_dp
        end if 

        return

    end function weighted_ave

    function weighted_ave_shepard(var,dist,shepard_exponent,mask)

        implicit none

        real(dp) :: var(:), dist(:)
        real(dp), optional :: shepard_exponent
        logical,  optional :: mask(:) 
        real(dp) :: shep_e, weight, numerator, denominator
        real(dp) :: weighted_ave_shepard
        integer :: i, n

        logical :: msk(size(var)) 

        msk = .TRUE. 
        if (present(mask)) msk = mask 

        shep_e = 2.0_dp 
        if (present(shepard_exponent)) shep_e = shepard_exponent

        numerator            = 0.0_dp
        denominator          = 0.0_dp
        weight               = 0.0_dp 

        n = size(var)

        do i = 1, n

            if (msk(i)) then   ! If not a missing point 
                ! See denominator in equation (2.17) and equation (2.19) in Reerink et al. (2010):
                weight = 1.0_dp / (dist(i)**shep_e)
                numerator   = numerator + var(i) * weight
                denominator = denominator + weight
            end if 

        end do 

        if (denominator .gt. 0.0_dp) then 
            weighted_ave_shepard = numerator / denominator
        else 
            weighted_ave_shepard = 0.0_dp 
        end if 

        return

    end function weighted_ave_shepard

    function quadrant_latlon(lon1, lat1, lon2, lat2)
        ! Determing the quadrant in which a 'projected point' (lon2,lat2)
        ! is situated relative to another (lon1,lat1).
        !   quadrants:
        !   II  |   I
        !       |     
        !  -----|-----
        !       |     
        !  III  |  IV 
        !
        implicit none 

        real(dp), intent(IN)  :: lon1, lat1, lon2, lat2
        integer :: quadrant, quadrant_latlon

        if( ((lon2 - lon1) .le.  180.0_dp .and. (lon2 - lon1) .gt. 0.0_dp) .or. &
            ((lon2 - lon1) .le. -180.0_dp)) then
            if(lat2 .gt. lat1) then
                quadrant = 1
            else
                quadrant = 4
            end if
        else
            if(lat2 .gt. lat1) then
                quadrant = 2
            else
                quadrant = 3
            end if
        end if

        quadrant_latlon = quadrant

        return

    end function quadrant_latlon

    function quadrant_cartesian(x1,y1,x2,y2)
        ! Determing the quadrant in which a point (x2,y2) is situated 
        ! relative to another grid point (x1,y1).
        !   quadrants:
        !   II  |   I
        !       |     
        !  -----|-----
        !       |     
        !  III  |  IV 
        !
        implicit none 

        real(dp), intent(IN)  :: x1,y1,x2,y2
        integer :: quadrant, quadrant_cartesian

        if (      x2 >  x1 .AND. y2 >= y1) then
            quadrant = 1
        else if ( x2 <= x1 .AND. y2 >  y1) then
            quadrant = 2
        else if ( x2 <  x1 .AND. y2 <= y1) then
            quadrant = 3
        else if ( x2 >= x1 .AND. y2 <= y1) then
            quadrant = 4
        else
            write(*,"(a,a)") "planet:: map_quadrant: ","error: quadrant not found."
            stop 
        end if

        quadrant_cartesian = quadrant

        return

    end function quadrant_cartesian

end module planet