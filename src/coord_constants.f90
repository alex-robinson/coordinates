module coord_constants
    
    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    
    ! Missing value and aliases
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large) and error index 
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 

    public 

end module coord_constants 