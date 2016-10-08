

program test

    use coordinates 
    use interp2D_conservative
    use grid_gen 

    use ncio 

    implicit none

    integer :: t

    double precision, parameter :: mv = -9999.d0 

    type(grid_class) :: grid
    double precision, dimension(:,:), allocatable :: var, var1, var2 
    integer, dimension(:,:), allocatable :: mask
    character(len=256) :: outfldr, file_input
    character(len=256) :: file_out_suffix, file_out 
    character(len=56)  :: varname 
    integer :: q 

    type(multigrid_class) :: mgrid 

    outfldr = "output/multigrid/"
    file_input      = "data/GRL-5KM_TOPO-B13.nc"
    file_out_suffix = "_TOPO.nc"
    varname         = "zb"

    call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=360,dy=5.d0,ny=600, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)
    call grid_allocate(grid,var1)
    call nc_read(file_input,varname,var1)

    call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    ! Test multigrid initialization 
    call multigrid_init(mgrid,grid,dx=[10.d0,20.d0,50.d0,100.d0])


    ! Load test data 
    call grid_allocate(grid,var)
!     call nc_read(file_input,varname,var)
    var(1:grid%G%nx-1,1:grid%G%ny-1) = var1
    var(:,grid%G%ny) = var(:,grid%G%ny-1)
    var(grid%G%nx,:) = var(grid%G%nx-1,:)

    do q = 1, mgrid%n_grids 

        ! Interp field to subgrid
        call grid_allocate(mgrid%grid(q),var1)
        call map_field_conservative(grid,mgrid%grid(q),varname,var,var1)

        ! Write to file 
        file_out = trim(outfldr)//trim(mgrid%grid(q)%name)//trim(file_out_suffix)
        call grid_write(mgrid%grid(q),file_out,xnm="xc",ynm="yc",create=.TRUE.)

        call nc_write(file_out,varname,var1,dim1="xc",dim2="yc")


        write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
                sum(var*grid%G%dx*grid%G%dy), &
                sum(var1*mgrid%grid(q)%G%dx*mgrid%grid(q)%G%dy), & 
                100*(sum(var1*mgrid%grid(q)%G%dx*mgrid%grid(q)%G%dy) - sum(var*grid%G%dx*grid%G%dy)) &
                        / sum(var*mgrid%grid(q)%G%dx*mgrid%grid(q)%G%dy)                  

    end do 

end program test