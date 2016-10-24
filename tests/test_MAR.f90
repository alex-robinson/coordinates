
!! TO COMPILE : 
!! gfortran -fcheck=all -I/opt/local/include -o test_MAR.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_MAR.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!! or
!! ifort -g -I/home/robinson/apps/netcdf/netcdf/include -o test_MAR.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_MAR.f90 -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

program test_MAR

    use coord 
    use ncio 

    implicit none

    type(grid_class) :: gMAR, gSICO
    type(map_class)  :: mMAR_SICO

    double precision, dimension(:,:), allocatable :: var_in,  var_out
    integer,          dimension(:,:), allocatable :: mask_out

    character(len=256) :: file_input, file_output, file_test, var_name
    integer :: t, ttot

    ! =======================================================================
    !
    ! Step 1: Define global input grid and load data that will be used here
    !
    ! =======================================================================

    ! Define file names for input and output of global grids
    file_input     = "/Users/robinson/wrk_local/data/MAR/ICE.1979-2006_ydaymean_noleapyear.nc"
    file_output    = "/Users/robinson/wrk_local/data/MAR/ICE.1979-2006_sicogrid_cc.nc"
    file_test      = "/Users/robinson/wrk_local/data/MAR/grid_MAR.nc"

    ! Define MAR grid and input variable field
    call grid_init(gMAR,name="MAR-25KM",mtype="stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-750.d0,dx=25.d0,nx=58,y0=-1200.d0,dy=25.d0,ny=108, &
                   lambda=320.d0,phi=72.d0,alpha=7.5d0)
    call grid_allocate(gMAR, var_in)

    call grid_write(gMAR,file_test,xnm="xc",ynm="yc",create=.TRUE.)

    ! Define SICOPOLIS grid and output variable field
    call grid_init(gSICO,name="SICO-20KM",mtype="polar_stereographic",units="kilometers",lon180=.TRUE., &
                   x0=-800.d0,dx=20.d0,nx=76,y0=-3400.d0,dy=20.d0,ny=151, &
                   lambda=-39.d0,phi=90.d0,alpha=7.5d0)
    call grid_allocate(gSICO, var_out)
    call grid_allocate(gSICO, mask_out)

    ! Output new grid to file
    call nc_create(file_output)
    call nc_write_dim(file_output,"xc",  x=gSICO%G%x,units="kilometers")
    call nc_write_dim(file_output,"yc",  x=gSICO%G%y,units="kilometers")
    call nc_write_dim(file_output,"time",x=1.d0,dx=1.d0,nx=365,units="years",calendar="365_day")
    
    call grid_write(gSICO,file_output,xnm="xc",ynm="yc",create=.FALSE.)
    

    ! =======================================================================
    !
    ! Step 2: Map the fields 
    !
    ! =======================================================================
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    ! max_neighbors is the maximum neighbors to be stored for each point
    ! lat_lim is the range of latitudes relative to a given point to check neighbor distances (to speed things up)
    call map_init(mMAR_SICO,gMAR,gSICO,max_neighbors=20,lat_lim=4.0d0,fldr="maps",load=.TRUE.)

    ! Map each field back to the SICO domain using the radius method

    ttot = 365
    var_name = "TT"

    do t = 1, ttot 

        ! Read current timestep of data from MAR grid
        call nc_read(file_input,var_name,var_in,start=[1,1,1,t],count=[gMAR%G%nx,gMAR%G%ny,1,1])
        
        ! Map the MAR data to the SICO grid (distance-weighted averaging within 50km radius)
        call map_field(mMAR_SICO,var_name,var_in,var_out,mask_out,"shepard",50.d3,fill=.TRUE.)

        ! Write new SICO gridded data to grid file 
        call nc_write(file_output,var_name,var_out,dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,t],count=[gSICO%G%nx,gSICO%G%ny,1] )
    end do 


end program test_MAR