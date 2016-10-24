
!! TO COMPILE : 
!! gfortran -fcheck=all -I/opt/local/include -o test_coord.x ../ncio/ncio3.f90 oblimap_projection_module.f90 geodesic.f90 planet.f90 coordinates.f90 test_coordinates.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!! or
!! ifort -g -I/home/robinson/apps/netcdf/netcdf/include -o test_coord.x ../ncio/ncio3.f90 oblimap_projection_module.f90 geodesic.f90 planet.f90 coordinates.f90 test_coordinates.f90 -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

program test_coordinates

    use coord 
    use ncio 
    
    implicit none

    type(grid_class) :: rembo_grl_lo, rembo_grl_hi, ecmwf
    type(map_class)  :: map1, map2, map3, map4

    double precision, dimension(:,:), allocatable :: zs1,   zs2, zs3
    integer,          dimension(:,:), allocatable :: mask1, mask2, mask3

    !!! ECMWF latlon grid 

    ! ECMWF 2-deg latlon 
    call grid_init(ecmwf,name="GRL-LATLON05",mtype="latlon",units="degrees",lon180=.TRUE., &
                     x0=-100.d0,dx=0.5d0,nx=231,y0=50.d0,dy=0.5d0,ny=81)

    call grid_write(ecmwf,"maps/grid_ecmwf.nc",xnm="lon",ynm="lat",create=.TRUE.)

    ! Oblique stereographic
    call grid_init(rembo_grl_lo,name="GRL-20KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                     x0=-750.d0,dx=20.d0,nx=76,y0=-1500.d0,dy=20.d0,ny=151, &
                     lambda=320.d0,phi=72.d0,alpha=7.5d0,x_e=0.d0,y_n=0.d0)

    call grid_write(rembo_grl_lo,"maps/grid_rembo-grl-lo.nc",xnm="xc",ynm="yc",create=.TRUE.)

    call grid_init(rembo_grl_hi,name="GRL-10KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                     x0=-750.d0,dx=10.d0,nx=151,y0=-1500.d0,dy=10.d0,ny=301, &
                     lambda=320.d0,phi=72.d0,alpha=7.5d0,x_e=0.d0,y_n=0.d0)

    call grid_write(rembo_grl_hi,"maps/grid_rembo-grl-hi.nc",xnm="xc",ynm="yc",create=.TRUE.)


    !! Begin mapping testing
    write(*,*) " == MAPPING =="
    write(*,*) 

    call map_init(map1,ecmwf,rembo_grl_lo,20,"maps",.TRUE.)
    call map_init(map2,rembo_grl_lo,ecmwf,20,"maps",.TRUE.)
    call map_init(map3,rembo_grl_lo,rembo_grl_hi,12,"maps",.TRUE.)
    call map_init(map4,rembo_grl_hi,rembo_grl_lo,16,"maps",.TRUE.)

    allocate( zs1(ecmwf%G%nx,ecmwf%G%ny), mask1(ecmwf%G%nx,ecmwf%G%ny) )
    allocate( zs2(rembo_grl_lo%G%nx,rembo_grl_lo%G%ny), mask2(rembo_grl_lo%G%nx,rembo_grl_lo%G%ny) )
    allocate( zs3(rembo_grl_hi%G%nx,rembo_grl_hi%G%ny), mask3(rembo_grl_hi%G%nx,rembo_grl_hi%G%ny) )
    call nc_read("data/ECMWF/NEW/ERA-INTERIM-GRL-invariant_historical_mon_197901-201212.nc",zs1,"z")
    zs1 = zs1 / 9.8d0 
    call nc_write("maps/grid_ecmwf.nc",zs1,"zs",dim1="lon",dim2="lat")

    call map_field(map1,"zs",zs1,zs2,mask2,"shepard",125.d3)
    call nc_write("maps/grid_rembo-grl-lo.nc",zs2,"zs",dim1="xc",dim2="yc")
    call nc_write("maps/grid_rembo-grl-lo.nc",mask2,"mask",dim1="xc",dim2="yc")

    call map_field(map2,"zs",zs2,zs1,mask1,"shepard",50.d3)
    call nc_write("maps/grid_ecmwf.nc",zs1,"zs1",dim1="lon",dim2="lat")
    call nc_write("maps/grid_ecmwf.nc",mask1,"mask1",dim1="lon",dim2="lat")

    call map_field(map3,"zs",zs2,zs3,mask3,"shepard",20.d3)
    call nc_write("maps/grid_rembo-grl-hi.nc",zs3,"zs",dim1="xc",dim2="yc")
    call nc_write("maps/grid_rembo-grl-hi.nc",mask3,"mask",dim1="xc",dim2="yc")

    call map_field(map4,"zs",zs3,zs2,mask2,"shepard",20.d3)
    call nc_write("maps/grid_rembo-grl-lo.nc",zs2,"zs2",dim1="xc",dim2="yc")
    call nc_write("maps/grid_rembo-grl-lo.nc",mask2,"mask2",dim1="xc",dim2="yc")

end program test_coordinates