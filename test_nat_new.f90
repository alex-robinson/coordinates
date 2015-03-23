

program test_nat_new

    implicit none 

    real(4), allocatable  :: x(:), y(:), z(:)   ! Original field
    real(4), allocatable  :: x1(:), y1(:), z1(:), z1x(:), z1y(:)  ! Interpolant locations
    real(4), allocatable  :: x2(:), y2(:), z2(:)                  ! Interpolant locations (new)

    integer :: nstore, npts, nx2, ny2, npts2 
    integer :: i, j, k 
    real(4) :: xmin, dx, ymin, dy 
    character(len=10) :: tmp 

    npts = 5000
    allocate(x(npts),y(npts),z(npts))
    allocate(x1(npts),y1(npts),z1(npts),z1x(npts),z1y(npts))

    nx2 = 20
    ny2 = 20 
    npts2 = nx2*ny2 
    allocate(x2(npts2),y2(npts2),z2(npts2))

    ! Read in data from file
    nstore = 13 
    open(nstore,file="fort.13",status="old")
    read(nstore,*) tmp 
    read(nstore,*) tmp 
    do i = 1, npts 
        read(nstore,*) x(i), y(i), z(i), z1(i), z1x(i), z1y(i)
    end do 
    close(nstore)

    ! Prepare new array of points for interpolation 
    xmin = minval(x) + (maxval(x)-minval(x))*0.1
    dx   = (maxval(x)-minval(x))*0.8 / nx2 
    ymin = minval(y) + (maxval(y)-minval(y))*0.1
    dy   = (maxval(y)-minval(y))*0.8 / ny2 

    write(*,*) "xmin,dx = ", xmin, dx
    write(*,*) "ymin,dy = ", ymin, dy 

    k = 0
    do i = 1, nx2 
    do j = 1, ny2 
        k = k+1
        x2(k) = xmin + (i-1)*dx 
        y2(k) = ymin + (j-1)*dy 
    end do 
    end do 
    
    write(*,*) "Old domain"
    write(*,*) "x :", minval(x), maxval(x)
    write(*,*) "y :", minval(y), maxval(y)
    write(*,*) "z :", minval(z), maxval(z)

    ! Interpolate all the points 
    do k = 1, npts2
        call tile4_natint_pt(x,y,z,x2(k),y2(k),zout=z2(k))
    end do 

    write(*,*) "New domain"
    write(*,*) "x2:", minval(x2), maxval(x2)
    write(*,*) "y2:", minval(y2), maxval(y2)
    write(*,*) "z2:", minval(z2), maxval(z2)
    
    open(nstore,file="fort.13_natnew",status="unknown")
    do i = 1, npts2 
        write(nstore,"(3g15.5)") x2(i), y2(i), z2(i)
    end do 
    close(nstore)

    
contains


    subroutine tile4_natint_pt(x,y,z,xout,yout,zout)

        implicit none 

        real(4), intent(IN)  :: x(:), y(:), z(:)   ! Original field
        real(4), intent(IN)  :: xout, yout         ! Interpolant location
        real(4), intent(OUT) :: zout               ! Interpolant value 

        integer :: npts 

        ! Internal variables needed for NNBRHG routine
        integer, parameter :: jcns   = 4      ! Not sure - related to "constraints"
        integer, parameter :: ltop   = 10000  ! Buffer space needed by routines
        integer, parameter :: knbmax = 50     ! Working array size
        integer, parameter :: nfmax  = 8 
        real(4), parameter :: epscn  = 0.01 
        real(4), parameter :: epspt  = 0.01 

        real(4), allocatable :: cn(:,:), jaddr(:)
        real(4), allocatable :: pt(:,:), val(:)
        integer, allocatable :: naddr(:), l(:)
        real(4), allocatable :: sbarea(:), ptoff(:,:), valoff(:)
        real(4), allocatable :: delsba(:,:)
        real(4) :: area, xpt, ypt, z0, zx0, zy0, zh, zxh, zyh
        integer :: nfree, nstart, nptsin, lptr, lbase
        integer :: iflag, igarb, nindex, knb, icase

        npts = size(x)

        allocate(cn(3,jcns))
        allocate(jaddr(jcns))
        allocate(pt(2,npts),val(npts))
        allocate(naddr(npts))
        allocate(l(ltop))
        allocate(sbarea(knbmax))
        allocate(ptoff(2,knbmax))
        allocate(valoff(knbmax))
        allocate(delsba(2,knbmax))

        ! Fill in the local variables 
        pt(1,:) = x 
        pt(2,:) = y 
        val     = z 
        xpt     = xout 
        ypt     = yout 
        
        ! Constraints array 
        CALL SETREC(cn,jcns,iflag,minval(x)-0.5,maxval(x)+0.5,minval(y)-0.5,maxval(y)+0.5) !    XMIN,XMAX,YMIN,YMAX)

        ! Call tile4m to generate the tesselation of input points
        CALL tile4m(cn,jcns,jaddr,pt,npts,naddr,nfree,nstart,nptsin,  &
                    l,ltop,lptr,lbase,epscn,epspt,iflag,igarb)

        if (minval(abs([0,4,5]-iflag)) .ne. 0) then 
            write(*,*) "tile4m:: ", xpt, ypt, iflag 
            stop 
        end if 

        ! Call NNBRHG to find the natural neighbor C1 interpolant (ZH)
        ! and C0 interpolant (Z0) at the point XPT, YPT 
        nstart = 1 

        call NNBRHG(cn,jcns,jaddr,pt,npts,naddr,nfree,nstart,                  &
                    nptsin,l,ltop,lptr,lbase,epscn,epspt,iflag,igarb,xpt,ypt,  &
                    nindex,area,sbarea,delsba,ptoff,knbmax,knb,icase,          &
                    val,z0,zx0,zy0,zh,zxh,zyh) 

!         write(*,*) "tile4_natint_pt:: "
!         write(*,*) minval(pt(1,:)), maxval(pt(1,:))
!         write(*,*) minval(pt(2,:)), maxval(pt(2,:))
!         write(*,*) minval(val),     maxval(val)
        write(*,*) xpt, ypt, z0, zh, zh-val(nindex), nindex, iflag, igarb
        
        ! Return the output of interest
        zout = zh 

        return

    end subroutine tile4_natint_pt 

end program
