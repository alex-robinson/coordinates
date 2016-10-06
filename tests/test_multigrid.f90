

program test

    use coordinates 
    use grid_gen 

    use ncio 

    implicit none

    integer :: t

    double precision, parameter :: missing_value = -9999.d0 

    type(grid_class) :: grid, gridlo, gridhi 
    double precision, dimension(:), allocatable   :: x, y
    double precision, dimension(:,:), allocatable :: var, varlo, varhi
    integer, dimension(:,:), allocatable :: mask, masklo, maskhi 
    character(len=256) :: outfldr,file_input, file_inputhi
    character(len=256) :: file_outlo, file_outhi
    integer :: i, j 

    outfldr = "output/interp/"
    file_input   = trim(outfldr)//"GRL-50KM_TOPO.nc"
    file_inputhi = trim(outfldr)//"GRL-20KM_TOPO.nc"
    file_outhi   = trim(outfldr)//"GRL-20KM_TOPO1.nc"
    file_outlo   = trim(outfldr)//"GRL-50KM_TOPO1.nc"

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

    ! Check if the grids are defined as the same map 
    write(*,*) "Same grid: ", compare_coord(grid,gridhi)


end program test