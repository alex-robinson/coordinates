

program test
 
    use coord 
    use ncio 
    use coordinates_mapping_conservative
!     use interp2D_conservative

    implicit none

    integer :: t

    double precision, parameter :: mv = -9999.d0 

    type(grid_class) :: grid
    double precision, dimension(:,:), allocatable :: var, var1, var2 
    integer, dimension(:,:), allocatable :: mask
    character(len=256) :: outfldr, file_input
    character(len=256) :: file_out_suffix, file_out 
    character(len=56)  :: varnames(4), varname 
    integer :: q, k 

    type(multigrid_class)   :: mgrid

    double precision :: target_val, current_val, err_percent 
    double precision :: xlim(2), ylim(2) 

    outfldr = "output/multigrid/"
    file_input      = "data/GRL-5KM_TOPO-RTOPO-2.0.1.nc"
    file_out_suffix = "_TOPO.nc"
    
    call grid_init(grid,name="GRL-5KM",mtype="polar_stereographic",units="kilometers", &
                        lon180=.TRUE.,x0=-720.d0,dx=5.0d0,nx=337,y0=-3450.d0,dy=5.0d0,ny=577, &
                        lambda=-45.d0,phi=70.d0)

    ! Test multigrid initialization 
!     call multigrid_init(mgrid,grid,dx=[10.d0,20.d0,50.d0,100.d0])
    call multigrid_init(mgrid,grid,dx=[20.d0]) 
    call grid_allocate(grid,var)
    
    xlim = [ -720.d0, 960.d0]
    ylim = [-3450.d0, 570.d0]
    
    varnames = ["z_srf","z_bed","H_ice","mask "]

    do k = 1, size(varnames)

        varname = varnames(k)
        write(*,*) k, trim(varname)

        ! Load test data
        call nc_read(file_input,varname,var)
        target_val = calc_grid_total(grid%G%x,grid%G%y,var,xlim=xlim,ylim=ylim)

        ! Loop over multigrids
        do q = 1, mgrid%n 

            ! Interp field to subgrid   
            call grid_allocate(mgrid%grid(q),var1)

            if (trim(varname) .ne. "mask") then 
!                 call map_field_conservative(mgrid%map_to(q),varname,var,var1)
                call map_field_conservative_map1(mgrid%map_to(q)%map,varname,var,var1,method="mean")

                current_val = calc_grid_total(mgrid%grid(q)%G%x,mgrid%grid(q)%G%y,var1,xlim=xlim,ylim=ylim)
                err_percent = 100.d0 * (current_val-target_val) / target_val
                write(*,"(a,3g12.4)") "mass comparison (hi, con, % diff): ", &
                        target_val, current_val, err_percent                  

            else 
!                 var1 = interp_nearest(x=grid%G%x,y=grid%G%y,z=var, &
!                                       xout=mgrid%grid(q)%G%x,yout=mgrid%grid(q)%G%y)

                call map_field_conservative_map1(mgrid%map_to(q)%map,varname,var,var1,method="count")

            end if 

            ! Write to file 
            file_out = trim(outfldr)//trim(mgrid%grid(q)%name)//trim(file_out_suffix)
            if (k.eq.1) call grid_write(mgrid%grid(q),file_out,xnm="xc",ynm="yc",create=.TRUE.)
            call nc_write(file_out,varname,var1,dim1="xc",dim2="yc")


            ! Remap to high resolution 
            call grid_allocate(grid,var2)

            if (trim(varname) .ne. "mask") then 
!                call map_field_conservative(mgrid%map_from(q),varname,var1,var2)

!                 call map_field_conservative_smooth(mgrid%map_from(q),mgrid%map_to(q), &
!                             mgrid%grid(q),grid,varname,var1,var2)
                
                call map_field_conservative_map1(mgrid%map_from(q)%map,varname,var1,var2,method="mean")

                current_val = calc_grid_total(mgrid%grid(q)%G%x,mgrid%grid(q)%G%y,var1,xlim=xlim,ylim=ylim)
                err_percent = 100.d0 * (current_val-target_val) / target_val
                write(*,"(a,3g12.4)") "mass comparison (hi, con-hi, % diff): ", &
                        target_val, current_val, err_percent                  

            else 

!                 var2 = interp_nearest(x=mgrid%grid(q)%G%x,y=mgrid%grid(q)%G%y,z=var1, &
!                                       xout=grid%G%x,yout=grid%G%y)

                call map_field_conservative_map1(mgrid%map_from(q)%map,varname,var1,var2,method="count")

            end if 

            file_out = trim(outfldr)//trim(mgrid%grid(q)%name)//"_5KM"//trim(file_out_suffix)
            if (k.eq.1) call grid_write(grid,file_out,xnm="xc",ynm="yc",create=.TRUE.)
            call nc_write(file_out,varname,var2,dim1="xc",dim2="yc")

                
        end do 
    end do 

end program test