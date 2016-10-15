

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

    type(multigrid_class)   :: mgrid

    double precision :: target_val, current_val, err_percent 
    double precision :: xlim(2), ylim(2) 

    outfldr = "output/multigrid/"
    file_input      = "data/GRL-5KM_TOPO-B13.nc"
    file_out_suffix = "_TOPO.nc"
    varname         = "zb"

    call grid_init(grid,name="GRL-5KM",mtype="stereographic",units="kilometers", &
                               lon180=.TRUE.,dx=5.d0,nx=361,dy=5.d0,ny=601, &
                               lambda=-40.d0,phi=72.d0,alpha=8.4d0)

    ! Test multigrid initialization 
!     call multigrid_init(mgrid,grid,dx=[10.d0,20.d0,50.d0,100.d0])
    call multigrid_init(mgrid,grid,dx=[50.d0])
    
    ! Load test data 
    call grid_allocate(grid,var)
    call nc_read(file_input,varname,var)

!     xlim = [minval(grid%G%x),maxval(grid%G%y)]
!     ylim = [minval(grid%G%y),maxval(grid%G%y)]
    xlim = [ -850.d0, 850.d0]
    ylim = [-1450.d0,1450.d0]
    
!     target_val = sum(var*grid%G%dx*grid%G%dy)
    target_val = calc_grid_total(grid%G%x,grid%G%y,var,xlim=xlim,ylim=ylim)

    ! Loop over multigrids
    do q = 1, mgrid%n_grids 

!         ! Generate maps (hi-lo,lo-hi)
!         call map_conservative_init(mgrid%map(q),grid,mgrid%grid(q))
!         call map_conservative_init(mlohi,mgrid%grid(q),grid)

        ! Interp field to subgrid   
        call grid_allocate(mgrid%grid(q),var1)
        call map_field_conservative(mgrid%map_to_mgrid(q),varname,var,var1)

        ! Write to file 
        file_out = trim(outfldr)//trim(mgrid%grid(q)%name)//trim(file_out_suffix)
        call grid_write(mgrid%grid(q),file_out,xnm="xc",ynm="yc",create=.TRUE.)
        call nc_write(file_out,varname,var1,dim1="xc",dim2="yc")

!         current_val = sum(var1*mgrid%grid(q)%G%dx*mgrid%grid(q)%G%dy)
        current_val = calc_grid_total(mgrid%grid(q)%G%x,mgrid%grid(q)%G%y,var1,xlim=xlim,ylim=ylim)
        err_percent = 100.d0 * (current_val-target_val) / target_val
        write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
                target_val, current_val, err_percent                  

        ! Remap to high resolution 
        call grid_allocate(grid,var2)
!         call map_field_conservative(mgrid%map_from_mgrid(q),varname,var1,var2)
        call map_field_conservative_smooth(mgrid%map_from_mgrid(q),mgrid%map_to_mgrid(q),varname,var1,var2)

        file_out = trim(outfldr)//trim(mgrid%grid(q)%name)//"_5KM"//trim(file_out_suffix)
        call grid_write(grid,file_out,xnm="xc",ynm="yc",create=.TRUE.)
        call nc_write(file_out,varname,var2,dim1="xc",dim2="yc")

!         current_val = sum(var2*grid%G%dx*grid%G%dy)
        current_val = calc_grid_total(mgrid%grid(q)%G%x,mgrid%grid(q)%G%y,var1,xlim=xlim,ylim=ylim)
        err_percent = 100.d0 * (current_val-target_val) / target_val
        write(*,"(a,3g12.4)") "mass comparison (hi, con-hi, % diff): ", &
                target_val, current_val, err_percent                  

    end do 

end program test