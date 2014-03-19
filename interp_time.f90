
module interp_time
    
    use interp1D 
    
    implicit none 

    private
    public :: convert_monthly1D_daily1D

contains


    function convert_monthly1D_daily1D(dat_mon,dat_m0,dat_m13,nd) result(dat_day)
        
        implicit none

        double precision, dimension(:), intent(IN)   :: dat_mon
        double precision,               intent(IN)   :: dat_m0, dat_m13 
        double precision, dimension(nd)  :: dat_day
        double precision, dimension(:), allocatable  :: dat_mon14

        integer :: nx, ny, nd, i, j
        double precision, dimension(:), allocatable :: x, xout 

        if (nd .ne. 360) then 
            write(*,*) "convert_monthly1D_daily1D:: ", &
            "Error: currently this routine only works with 360-day years."
            stop 
        end if 

        allocate(x(14),xout(nd)) 
        allocate(dat_mon14(14))

        dat_mon14(1)    = dat_m0
        dat_mon14(14)   = dat_m13 
        dat_mon14(2:13) = dat_mon 

        ! Define x-values for converting to daily data from months using splines
        do i = 1, 14
            x(i) = dble(i-1)*30-15
        end do 
        do i = 1, nd
            xout(i) = dble(i)
        end do 

        dat_day = interp_spline(x,dat_mon14,xout)

        return 

    end function convert_monthly1D_daily1D


end module interp_time 
