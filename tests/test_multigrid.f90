

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

    type(multigrid_class) :: mgrid 

    outfldr = "output/interp/"
    file_input   = trim(outfldr)//"GRL-5KM_TOPO.nc"
    file_inputhi = trim(outfldr)//"GRL-20KM_TOPO.nc"
    file_outhi   = trim(outfldr)//"GRL-20KM_TOPO1.nc"
    file_outlo   = trim(outfldr)//"GRL-50KM_TOPO1.nc"

    call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    ! Test multigrid initialization 
    call multigrid_init(mgrid,grid,dx=[20.d0,50.d0,100.d0])

end program test