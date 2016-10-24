
program test_etopo

    use coord 
    use ncio 
    
    implicit none

    type(grid_class) :: grid0, grid1, grid2, grid3, grid4
    double precision, allocatable :: var0(:,:), var1(:,:) 
    double precision :: tmplon(720), tmplat(360)
    
    character(len=256) :: file0
    double precision :: table(4,4)
    integer :: i 

    write(*,*) 

    ! =======================================================================
    !
    ! Define global input grid and load input data
    !
    ! =======================================================================

    file0  = "data/ETOPO1-ICE-050deg.nc"

    ! Define input latlon grid
    call nc_read(file0,"lon",tmplon)
    call nc_read(file0,"lat",tmplat)
    call grid_init(grid0,name="GLOBAL-050deg",mtype="latlon",units="degrees",lon180=.TRUE., &
                     x=tmplon,y=tmplat)

    ! Allocate arrays of the size of the global grid
    call grid_allocate(grid0, var0)

    ! Load original data
    call nc_read(file0,"z",var0)

    ! =========================================================
    !
    ! PROJECTED GRID DEFINITIONS
    !
    ! =========================================================

    call grid_init(grid1,name="ANT-40KM",mtype="polar_stereographic",units="kilometers", &
                   lon180=.TRUE.,dx=40.d0,nx=156,dy=40.d0,ny=146,lambda=0.d0,phi=-71.d0)

    call grid_init(grid2,name="ANT-20KM",mtype="polar_stereographic",units="kilometers", &
                   lon180=.TRUE.,dx=20.d0,nx=311,dy=20.d0,ny=291,lambda=0.d0,phi=-71.d0)

    call grid_init(grid3,name="ANT-10KM",mtype="polar_stereographic",units="kilometers", &
                   lon180=.TRUE.,dx=10.d0,nx=621,dy=10.d0,ny=581,lambda=0.d0,phi=-71.d0)

!     call grid_init(grid4,name="ANT-5KM",mtype="polar_stereographic",units="kilometers", &
!                    lon180=.TRUE.,dx=5.d0,nx=1241,dy=5.d0,ny=1161,lambda=0.d0,phi=-71.d0)


    ! =========================================================
    !
    ! TEST MAPPING SPEEDS
    !
    ! =========================================================
    table = 0.d0 

    table(1,:) = test_mapping(var0,grid0,grid1,niter=100,lat_lim=0.5d0)
    table(2,:) = test_mapping(var0,grid0,grid2,niter=100,lat_lim=0.5d0)
    table(3,:) = test_mapping(var0,grid0,grid3,niter=100,lat_lim=0.5d0)
!     table(4,:) = test_mapping(var0,grid0,grid4,niter=100,lat_lim=0.5d0)

    write(*,*)
    write(*,*) "Mapping/interpolation summary"
    write(*,"(5a14)") "npts", "map_time", "n_iter", "interp_time", "   (min.)"

    do i = 1, 4 
        write(*,"(f14.1,f14.3,f14.1,f14.5)") table(i,:)
    end do 


    write(*,*) 
    write(*,*) "Done!"
    write(*,*) 

contains

    function test_mapping(var0,grid0,grid1,niter,lat_lim) result(table)

        implicit none 

        double precision :: var0(:,:)
        type(grid_class) :: grid0, grid1
        integer          :: niter 
        double precision :: lat_lim 
        double precision :: table(4)

        ! Local variables
        type(map_class)    :: map
        double precision, allocatable :: var1(:,:)
        integer,          allocatable :: mask1(:,:)

        double precision   :: start, finish 
        double precision   :: mapping_time, interp_time 
        character(len=256) :: file1
        integer            :: i 

        ! Allocate output arrays 
        call grid_allocate(grid1,var1)
        call grid_allocate(grid1,mask1)
        
        ! First, perform mapping once, record time
        write(*,*) "" 
        write(*,*) "Mapping ..."
        write(*,*) "" 

        call cpu_time(start)

        ! Create a map object for grid0=>grid1 mapping
        call map_init(map,grid0,grid1,max_neighbors=20,lat_lim=lat_lim,dist_max=200d0, &
                        fldr="maps",load=.TRUE.)
    
        call cpu_time(finish)
        mapping_time = (finish-start)/60.0d0


        ! Next perform interpolation niter times
        write(*,*) "" 
        write(*,*) "Interpolating ..."
        write(*,*) "" 

        call cpu_time(start)
        
        do i = 1, niter 

            call map_field(map,"zs",var0,var1,mask1,method="radius")

        end do 

        call cpu_time(finish)
        interp_time = (finish-start)/60.0d0

        write(*,*) "Summary: "//trim(map%name1)//" => "//trim(map%name2)
        write(*,"(a,f10.2,i4,f10.5)") "Mapping time, niter, interp. time (min.): ", &
                                         mapping_time, niter, interp_time 

        ! Write results to netcdf file 
        call map_field(map,"zs",var0,var1,mask1,method="radius")

        file1  = "output/"//trim(grid1%name)//"_ETOPO1-ICE-050deg.nc"
        call grid_write(grid1,file1,xnm="xc",ynm="yc",create=.TRUE.)
        call nc_write(file1,"zs",var1,dim1="xc",dim2="yc")

        table = [dble(grid1%npts),mapping_time,dble(niter),interp_time]

        return 


    end function test_mapping 


end program test_etopo

