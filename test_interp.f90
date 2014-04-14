

program test

    use interp1D
    use interp2D
    use interp_time
    use coordinates 
    use ncio 

    implicit none

    integer :: t

    double precision, parameter :: missing_value = -9999.d0 

    type(grid_class) :: grid, gridlo, gridhi 
    double precision, dimension(:), allocatable   :: x, y
    double precision, dimension(:,:), allocatable :: var, varlo, varhi
    character(len=256) :: outfldr,file_input, file_inputhi
    character(len=256) :: file_outlo, file_outhi
    integer :: i, j 

    outfldr = "output/interp/"
    file_input   = trim(outfldr)//"GRL-50KM_TOPO.nc"
    file_inputhi = trim(outfldr)//"GRL-20KM_TOPO.nc"
    file_outhi   = trim(outfldr)//"GRL-20KM_TOPO1.nc"

    allocate(x(37),y(61))
    call nc_read(file_input,"xc",x)
    call nc_read(file_input,"yc",y)
    call grid_init(grid,name="GRL50KM",mtype="stereographic",units="km",lon180=.TRUE.,&
                   lambda=-40d0,phi=72.d0,alpha=7.5d0,x=x,y=y)

    deallocate(x)
    deallocate(y)
    allocate(x(90),y(150))
    call nc_read(file_inputhi,"xc",x)
    call nc_read(file_inputhi,"yc",y)
    call grid_init(gridhi,name="GRL20KM",mtype="stereographic",units="km",lon180=.TRUE., &
                   lambda=-40d0,phi=72.d0,alpha=7.5d0,x=x,y=y)

    ! Allocate grid variables for testing 
    call grid_allocate(grid,var)
    call grid_allocate(gridhi,varhi)
    
    ! Load an example low resolution variable
    call nc_read(file_input,"zs",var)
    write(*,*) "zs range: ",minval(var), maxval(var)

    ! Interpolation one row 
    varhi = missing_value 

!     write(*,*) 
!     do i = 5, gridhi%G%nx-5
!         do j = 1, grid%G%nx
!             if ()
!             varhi(i,:) = interp_linear(grid%G%x,var(,gridhi%G%x) 
!     end do 

!     i = 30 
!     j = 74 
!     varhi(:,j) = interp_linear(grid%G%x,var(:,i),gridhi%G%x) 

    ! Output test interpolation
!     call grid_write(gridhi,file_outhi,xnm="xc",ynm="yc",create=.TRUE.)
!     call nc_write(file_outhi,"zs",varhi(:,j),dim1="xc")

!     varhi(:,j) = interp_spline(grid%G%x,var(:,i),gridhi%G%x) 
!     call nc_write(file_outhi,"zs_sp",varhi(:,j),dim1="xc")
    

    ! Test bilinear interpolation
    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)

    call grid_write(gridhi,file_outhi,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write(file_outhi,"zs",varhi,dim1="xc",dim2="yc")

end program test