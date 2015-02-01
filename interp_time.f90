
module interp_time
    
    use interp1D 
    
    implicit none 

    interface time_average 
        module procedure time_average_2D 
        module procedure time_average_3D 
        module procedure time_average_4D 
    end interface 

    private
    public :: time_average 
    public :: convert_monthly_daily_1D, convert_monthly_daily_3D
    public :: convert_monthly_daily_4D

contains

    subroutine convert_monthly_daily_4D(mon4D,day4D,days)

        implicit none 

        double precision, dimension(:,:,:,:) :: mon4D, day4D 
        double precision, dimension(:,:), allocatable :: mon0, mon13
        integer, dimension(:), optional :: days
        integer :: nx, ny, nm, nd, nk 
        integer :: i, j, k 

        nx = size(mon4D,1)
        ny = size(mon4D,2)
        nm = size(mon4D,3)
        nk = size(mon4D,4)

        nd = size(day4D,3)

        if (allocated(mon0))  deallocate(mon0)
        if (allocated(mon13)) deallocate(mon13) 
        allocate(mon0(nx,ny),mon13(nx,ny))

        do k = 1, nk 

            write(*,*) "monthly=>daily : ",k,"/",nk

            if (k .eq. 1) then 
                mon0  = mon4D(:,:,12,1)
                mon13 = mon4D(:,:,1,k+1)
            else if (k .eq. nk) then 
                mon0  = mon4D(:,:,12,k-1)
                mon13 = mon4D(:,:,1,k)
            else
                mon0  = mon4D(:,:,12,k-1)
                mon13 = mon4D(:,:,1,k+1)
            end if 

            do i = 1, nx 
            do j = 1, ny 
                day4D(i,j,:,k) = convert_monthly_daily_1D(mon4D(i,j,:,k), &
                                 mon0(i,j),mon13(i,j),nd,days)
            end do 
            end do 

        end do 

        return 

    end subroutine convert_monthly_daily_4D

    subroutine convert_monthly_daily_3D(mon3D,day3D,days)

        implicit none 

        double precision, dimension(:,:,:) :: mon3D, day3D 
        double precision, dimension(:,:), allocatable :: mon0, mon13
        integer, dimension(:), optional :: days 
        integer :: nx, ny, nm, nd, nk 
        integer :: i, j, k 

        nx = size(mon3D,1)
        ny = size(mon3D,2)
        nm = size(mon3D,3)

        nd = size(day3D,3)

        if (allocated(mon0))  deallocate(mon0)
        if (allocated(mon13)) deallocate(mon13) 
        allocate(mon0(nx,ny),mon13(nx,ny))

        mon0  = mon3D(:,:,12)
        mon13 = mon3D(:,:,1)

        do i = 1, nx 
        do j = 1, ny 
            day3D(i,j,:) = convert_monthly_daily_1D(mon3D(i,j,:), &
                             mon0(i,j),mon13(i,j),nd,days)
        end do 
        end do 

        return 

    end subroutine convert_monthly_daily_3D

    function convert_monthly_daily_1D(dat_mon,dat_m0,dat_m13,nd,days) result(dat_day)
        
        implicit none

        double precision, dimension(:), intent(IN)   :: dat_mon
        double precision,               intent(IN)   :: dat_m0, dat_m13
        integer :: nd 
        integer,          dimension(:), intent(IN), optional :: days 
        double precision, dimension(nd)  :: dat_day
        double precision, dimension(:), allocatable  :: dat_mon14
        double precision, dimension(:), allocatable :: x, xout 
        integer :: i 

        allocate(x(14),xout(nd)) 
        allocate(dat_mon14(14))

        ! Define x-values for converting to daily data from months using splines
        do i = 1, 14
            x(i) = dble(i-1)*30-15
        end do 
        
        if (present(days)) then 
            xout = days 
        else 
            do i = 1, nd
                xout(i) = dble(i)
            end do 
        end if 

        if (.not. present(days) .and. nd .ne. 360) then 
            write(*,*) "convert_monthly_daily_1D:: ", &
            "Error: This routine only works with 360-day years &
             &unless specific days are provided."
            stop 
        else if (nd .ne. size(xout)) then 
            write(*,*) "convert_monthly_daily_1D:: ", &
            "Error: Number of days given and size of fields do not match."
            stop 
        end if 

        if (minval(xout) .lt. x(1) .or. maxval(xout) .gt. x(14)) then 
            write(*,*) "convert_monthly_days_1D:: ", &
            "Error: currently this routine only works with 360-day years."
            stop 
        end if 

        dat_mon14(1)    = dat_m0
        dat_mon14(14)   = dat_m13 
        dat_mon14(2:13) = dat_mon 

        
        dat_day = interp_spline(x,dat_mon14,xout)

        return 

    end function convert_monthly_daily_1D

    function time_average_2D(var2D) result(var1D)
        implicit none

        double precision, dimension(:,:) :: var2D 
        double precision, dimension(size(var2D,1)) :: var1D 
        integer :: nx, ny, nk, i, j 

        nx = size(var2D,1) 
        nk = size(var2D,2)

        do i = 1, nx 
            var1D(i) = sum(var2D(i,:)) / dble(nk)
        end do 

        return
    end function time_average_2D 

    function time_average_3D(var3D) result(var2D)
        implicit none

        double precision, dimension(:,:,:) :: var3D 
        double precision, dimension(size(var3D,1),size(var3D,2)) :: var2D 
        integer :: nx, ny, nk, i, j 

        nx = size(var3D,1)
        ny = size(var3D,2) 
        nk = size(var3D,3)

        do i = 1, nx 
            do j = 1, ny 
                var2D(i,j) = sum(var3D(i,j,:)) / dble(nk)
            end do 
        end do 

        return
    end function time_average_3D 

    function time_average_4D(var4D) result(var3D)
        implicit none

        double precision, dimension(:,:,:,:) :: var4D 
        double precision, dimension(size(var4D,1),size(var4D,2),size(var4D,3)) :: var3D 
        integer :: nx, ny, nm, nk, i, j, m

        nx = size(var4D,1)
        ny = size(var4D,2) 
        nm = size(var4D,3)
        nk = size(var4D,4)

        do i = 1, nx 
            do j = 1, ny
                do m = 1, nm  
                    var3D(i,j,m) = sum(var4D(i,j,m,:)) / dble(nk)
                end do 
            end do 
        end do 

        return
    end function time_average_4D 

end module interp_time 
