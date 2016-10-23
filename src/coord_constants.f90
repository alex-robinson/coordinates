module coord_constants
    
    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the coord library (sp,dp)
    integer,  parameter :: prec = sp 


    ! Missing value and aliases
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large) and error index 
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    public 

contains 

    function get_coord_precision() result(coord_prec)

        implicit none 

        integer :: coord_prec 

        coord_prec = kind(prec)

        return 

    end function get_coord_precision

end module coord_constants 