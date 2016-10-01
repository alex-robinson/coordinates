

program test_loess

    use loess 

    implicit none 


    double precision, allocatable :: x(:), y(:), yfit(:) 
    integer :: tmp(20)
    integer :: ny = 2014-1880+1
    character(len=1000) :: tmpc 
    integer :: k 

    allocate(x(ny),y(ny),yfit(ny))

    open(10,file="data/GLB.Ts+dSST.txt")

    ! Read header 
    do k=1, 8
        read(10,*) tmpc 
    end do 

    ! Read data (annual mean => col 14)
    do k = 1, ny 
        read(10,*) tmp 
        x(k) = dble(tmp(1)) 
        y(k) = dble(tmp(14))*0.1d0 
    end do 

    close(10) 

    write(*,*) "years: ", minval(x), maxval(x)
    write(*,*) "dT:    ", minval(y), maxval(y)

    y(10:15) = -9999.0 
    y(50:60) = -9999.0 

    yfit = loess_smooth(x,y,L=5.d0,missing_value=-9999.d0)

    write(*,*) "dTfit: ", minval(yfit), maxval(yfit)

    open(11,file="output/loess_check.txt")

    write(11,"(3a14)") "year", "dT", "dTfit"

    do k = 1, ny 
        write(11,"(i14,2f14.2)") int(x(k)), y(k), yfit(k) 
    end do 

end program test_loess 


