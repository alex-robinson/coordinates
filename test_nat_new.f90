

program test_nat_new

    implicit none 









contains


    subroutine tile4_natint_pt(x,y,z,xout,yout,zout)

        implicit none 

        real(4), intent(IN)  :: x(:), y(:), z(:)   ! Original field
        real(4), intent(IN)  :: xout, yout         ! Interpolant location
        real(4), intent(OUT) :: zout               ! Interpolant value 

        integer :: npts 

        ! Internal variables needed for NNBRHG routine
        integer, parameter :: jcns   = 4     ! Not sure
        integer, parameter :: ltop   = 6000  ! Not sure 
        integer, parameter :: knbmax = 50    ! Working array size
        integer, parameter :: nfmax  = 8 
        real(4), parameter :: epscn  = 0.0 
        real(4), parameter :: epspt  = 0.0 

        real(4), allocatable :: cn(:,:), jaddr(:)
        real(4), allocatable :: pt(:,:), val(:)
        integer, allocatable :: naddr(:), l(:)
        real(4), allocatable :: sbarea(:), ptoff(:,:), valoff(:), grad(:,:)
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
        allocate(grad(2,npts))
        allocate(delsba(2,knbmax))

        ! Fill in the local variables 
        pt(1,:) = x 
        pt(2,:) = y 
        val     = z 
        xpt     = xout 
        ypt     = yout 

        ! Call NNBRHG to find the natural neighbor C1 interpolant (ZH)
        ! and C0 interpolant (Z0) at the point XPT, YPT 

        call NNBRHG(cn,jcns,jaddr,pt,npts,naddr,nfree,nstart,                  &
                    nptsin,l,ltop,lptr,lbase,epscn,epspt,iflag,igarb,xpt,ypt,  &
                    nindex,area,sbarea,delsba,ptoff,knbmax,knb,icase,          &
                    val,z0,zx0,zy0,zh,zxh,zyh) 

        return

        end subroutine tile4_natint_pt 

end program