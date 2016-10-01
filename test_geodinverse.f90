


program geodinverse

    use planet 

    implicit none 

    double precision :: a, f, dist

    ! WGS84 values
    a = 6378137d0
    f = 1d0/298.257223563d0

    ! Default fortran test case for geodinverse.for
    ! in legacy GeographicLib-1.32 code:
    ! http://geographiclib.sourceforge.net/html/Fortran/
    dist = planet_distance(a,f,30.d0, 0.d0, 29.5d0, 179.5d0)

    write(*,*) "lat,lon 1:", " 30.0 N","   0.0 E"
    write(*,*) "lat,lon 2:", " 29.5 N"," 179.5 E"
    write(*,*) "dist =", dist, "m" 

end program