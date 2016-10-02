module interp2D_conservative
    
    use coordinates 
    use coordinates_mapping 

    use polygons 
    use geodesic 
    
    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    private 

contains 

    function interpconserv1(x,y,z,xout,yout) result(zout)
        ! Calculate 1st order conservative interpolation for a 
        ! point given its nearest neighbors  
        real(dp) :: x(:), y(:), z(:)
        real(dp) :: xout, yout, zout 


        return 

    end function interpconserv1


end module interp2D_conservative 
