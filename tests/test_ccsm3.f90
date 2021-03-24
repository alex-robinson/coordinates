
!! TO COMPILE : 
!! gfortran -fcheck=all -I/opt/local/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!! or
!! ifort -g -I/home/robinson/apps/netcdf/netcdf/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

program test_ccsm3

    use coord 
    use ncio 
    
    implicit none

    type(grid_class) :: gCCSM3, gREG
    type(map_class)  :: mCCSM3_REG, mREG_CCSM3

    type vars_type 
        character (len=256) :: name 
        integer             :: nx, ny 
        double precision    :: lambda, phi, alpha
        double precision, dimension(:,:), allocatable :: Ts, MB, Hs
        integer,          dimension(:,:), allocatable :: mask

    end type 

    type(vars_type) :: CCSM3a, CCSM3b, REG

    character(len=256) :: file_input, file_gCCSM3a, file_GCCSM3b, file_gREG
    character(len=256) :: file_new 

    double precision :: tmplon(128), tmplat(64)

    integer :: q 

    type(map_scrip_class) :: mps1, mps2

    ! ==== SCRIPS testing ========

    ! call map_scrip_load(mps1,"ccsm3","ANT-20KM","grids")
    ! call map_scrip_load(mps2,"ANT-20KM","ccsm3","grids")

    ! ============================

    ! =======================================================================
    !
    ! Step 1: Define global input grid and load data that will be used here
    !
    ! =======================================================================

    ! Define file names for input and output of global grids
    file_input     = "data/ccsm_example_dec_feb_pd.nc"
    file_gCCSM3a   = "output/ccsm3/grid_CCSM3-T42a.nc"
    file_gCCSM3b   = "output/ccsm3/grid_CCSM3-T42b_quadrant.nc"
    
    ! CCSM3 T42 latlon grid
    call nc_read(file_input,"lon",tmplon)
    call nc_read(file_input,"lat",tmplat)
    call grid_init(gCCSM3,name="CCSM3-T42",mtype="gaussian",units="degrees",lon180=.FALSE., &
                     x=tmplon,y=tmplat)
!     call grid_init(gCCSM3,filename="maps/grid_CCSM3-T42.txt",x=tmplon,y=tmplat)

    ! Write grid information to file
    call grid_write(gCCSM3,file_gCCSM3a,xnm="lon",ynm="lat",create=.TRUE.)
    
    ! Allocate arrays of the size of the global grid
    call grid_allocate(gCCSM3, CCSM3a%Ts)
    call grid_allocate(gCCSM3, CCSM3a%MB)
    call grid_allocate(gCCSM3, CCSM3a%Hs)
    call grid_allocate(gCCSM3, CCSM3a%mask)

    call grid_allocate(gCCSM3, CCSM3b%Ts)
    call grid_allocate(gCCSM3, CCSM3b%MB)
    call grid_allocate(gCCSM3, CCSM3b%Hs)
    call grid_allocate(gCCSM3, CCSM3b%mask)

    ! Load original GCM data
    call nc_read(file_input,"TS",CCSM3a%Ts)
    call nc_read(file_input,"PRECSC",CCSM3a%MB)
    call nc_read(file_input,"PRECSL",CCSM3a%Hs)
    CCSM3a%MB = (CCSM3a%MB + CCSM3a%Hs) * (1.000d0/0.910d0) * 31556926.d0  ! mie/a
    call nc_read(file_input,"PHIS",CCSM3a%Hs)
    CCSM3a%Hs = CCSM3a%Hs / 9.81d0

    CCSM3b%Ts = CCSM3a%Ts 
    CCSM3b%MB = CCSM3a%MB 
    CCSM3b%Hs = CCSM3a%Hs 
    
    ! Write original data to grid file 
    call nc_write(file_gCCSM3a,"Ts",CCSM3a%Ts,dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3a,"MB",CCSM3a%MB,dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3a,"Hs",CCSM3a%Hs,dim1="lon",dim2="lat")

    ! =======================================================================
    !
    ! Step 2: Generate region specific grids
    !
    ! =======================================================================

    ! REG%name = "ANT-20KM"
    ! REG%nx   = 281
    ! REG%ny   = 281 
    ! REG%lambda =   0.d0 
    ! REG%phi    = -90.d0
    ! REG%alpha  =  19.d0 

    ! REG%name = "GRL-20KM"
    ! REG%nx   = 76
    ! REG%ny   = 151 
    ! REG%lambda = 320.d0 
    ! REG%phi    = 72.d0
    ! REG%alpha  = 7.5d0

    ! REG%name = "GRL-20KM"
    ! REG%nx   = 76
    ! REG%ny   = 151 
    ! REG%lambda = -45.d0 
    ! REG%phi    =  70.d0
    ! REG%alpha  = 7.5d0

    REG%name = "NH-40KM"
    REG%nx   = 313
    REG%ny   = 313 
    REG%lambda = -44.d0 
    REG%phi    =  70.d0
    REG%alpha  = 0.0d0

!     REG%name = "HIM-20KM"
!     REG%nx   = 200
!     REG%ny   = 200 
!     REG%lambda = 90.d0 
!     REG%phi    = 32.d0
!     REG%alpha  = 14.5d0

    ! Oblique stereographic grid for region of interest
    ! call grid_init(gREG,name=REG%name,mtype="stereographic",units="km",lon180=.FALSE., &
    !                  dx=20.d0,nx=REG%nx,dy=20.d0,ny=REG%ny, &
    !                  lambda=REG%lambda,phi=REG%phi,alpha=REG%alpha)
    call grid_init(gREG,name=REG%name,mtype="polar_stereographic",units="km",lon180=.FALSE., &
                     dx=40.d0,nx=REG%nx,dy=40.d0,ny=REG%ny, &
                     lambda=REG%lambda,phi=REG%phi,alpha=REG%alpha)

    ! call grid_init(gREG,name=REG%name,mtype="polar_stereographic",units="km",lon180=.TRUE., &
    !                  x0=-720.d0,dx=20.d0,nx=REG%nx,y0=-3450.0d0,dy=20.d0,ny=REG%ny, &
    !                  lambda=REG%lambda,phi=REG%phi,alpha=REG%alpha)

    ! Allocate arrays of the size of the regional grid
    call grid_allocate(gREG, REG%Ts)
    call grid_allocate(gREG, REG%MB)
    call grid_allocate(gREG, REG%Hs)
    call grid_allocate(gREG, REG%mask)
    REG%mask = 0    

    ! =======================================================================
    !
    ! Step 3: Map the fields 
    !
    ! =======================================================================
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Initialize 'to' and 'fro' mappings
    ! max_neighbors is the maximum neighbors to be stored for each point
    ! lat_lim is the range of latitudes relative to a given point to check neighbor distances (to speed things up)
    call map_init(mCCSM3_REG,gCCSM3,gREG,max_neighbors=6, lat_lim=3.0d0,dist_max=1000d0,fldr="maps",load=.TRUE.)
    call map_init(mREG_CCSM3,gREG,gCCSM3,max_neighbors=30,lat_lim=2.0d0,dist_max=1000d3,fldr="maps",load=.TRUE.)

    ! === Quadrant method ===
    file_gCCSM3b   = "output/ccsm3/grid_CCSM3-T42b_quadrant.nc"
    file_gREG      = "output/ccsm3/grid_"//trim(REG%name)//"_quadrant.nc"
    
    ! Map each field to the regional domain using the quadrant method (no max_distance required here)
    call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="quadrant")
    call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="quadrant")
    call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="quadrant") 

    CCSM3b%Ts = CCSM3a%Ts 
    CCSM3b%MB = CCSM3a%MB 
    CCSM3b%Hs = CCSM3a%Hs 
    call map_field(mREG_CCSM3,"Ts",REG%Ts,CCSM3b%Ts,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"MB",REG%MB,CCSM3b%MB,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"Hs",REG%Hs,CCSM3b%Hs,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)

    ! Write new regional data to grid file
    call grid_write(gREG,file_gREG,xnm="xc",ynm="yc",create=.TRUE.) 
    call nc_write(file_gREG,"Ts",  REG%Ts,  dim1="xc",dim2="yc",grid_mapping="crs") 
    call nc_write(file_gREG,"MB",  REG%MB,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"Hs",  REG%Hs,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"mask",REG%mask,dim1="xc",dim2="yc",grid_mapping="crs")
    
    ! Write new CCSM3 data to grid file 
    call grid_write(gCCSM3,file_gCCSM3b,xnm="lon",ynm="lat",create=.TRUE.)
    call nc_write(file_gCCSM3b,"Ts",  CCSM3b%Ts,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"MB",  CCSM3b%MB,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"Hs",  CCSM3b%Hs,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"mask",CCSM3b%mask,dim1="lon",dim2="lat")

    ! Calculate statistics concerning remapping
    ! (as in Table 3 of Reerink et al, 2010)
    write(*,*) "=== Quadrant method ==="
    call grid_stats("Ts",CCSM3a%Ts,CCSM3b%Ts,CCSM3b%mask)
    call grid_stats("MB",CCSM3a%MB,CCSM3b%MB,CCSM3b%mask)
    call grid_stats("Hs",CCSM3a%Hs,CCSM3b%Hs,CCSM3b%mask)

    ! === Bilinear method ===
    file_gCCSM3b   = "output/ccsm3/grid_CCSM3-T42b_bilinear.nc"
    file_gREG      = "output/ccsm3/grid_"//trim(REG%name)//"_bilinear.nc"
    
    call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="bilinear")
    call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="bilinear")
    call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="bilinear")

    CCSM3b%Ts = CCSM3a%Ts 
    CCSM3b%MB = CCSM3a%MB 
    CCSM3b%Hs = CCSM3a%Hs 
    call map_field(mREG_CCSM3,"Ts",REG%Ts,CCSM3b%Ts,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"MB",REG%MB,CCSM3b%MB,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"Hs",REG%Hs,CCSM3b%Hs,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)

    ! Write new regional data to grid file
    call grid_write(gREG,file_gREG,xnm="xc",ynm="yc",create=.TRUE.) 
    call nc_write(file_gREG,"Ts",  REG%Ts,  dim1="xc",dim2="yc",grid_mapping="crs") 
    call nc_write(file_gREG,"MB",  REG%MB,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"Hs",  REG%Hs,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"mask",REG%mask,dim1="xc",dim2="yc",grid_mapping="crs")

    ! Write new CCSM3 data to grid file 
    call grid_write(gCCSM3,file_gCCSM3b,xnm="lon",ynm="lat",create=.TRUE.)
    call nc_write(file_gCCSM3b,"Ts",  CCSM3b%Ts,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"MB",  CCSM3b%MB,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"Hs",  CCSM3b%Hs,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"mask",CCSM3b%mask,dim1="lon",dim2="lat")

    ! Calculate statistics concerning remapping
    ! (as in Table 3 of Reerink et al, 2010)
    write(*,*) "=== Bilinear method ==="
    call grid_stats("Ts",CCSM3a%Ts,CCSM3b%Ts,CCSM3b%mask)
    call grid_stats("MB",CCSM3a%MB,CCSM3b%MB,CCSM3b%mask)
    call grid_stats("Hs",CCSM3a%Hs,CCSM3b%Hs,CCSM3b%mask)

    
    ! === Nearest-neighbor method ===
    file_gCCSM3b   = "output/ccsm3/grid_CCSM3-T42b_nn.nc"
    file_gREG      = "output/ccsm3/grid_"//trim(REG%name)//"_nn.nc"

    call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="nn")
    call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="nn")
    call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="nn")

    CCSM3b%Ts = CCSM3a%Ts 
    CCSM3b%MB = CCSM3a%MB 
    CCSM3b%Hs = CCSM3a%Hs 
    call map_field(mREG_CCSM3,"Ts",REG%Ts,CCSM3b%Ts,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"MB",REG%MB,CCSM3b%MB,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"Hs",REG%Hs,CCSM3b%Hs,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)
    
    ! Write new regional data to grid file
    call grid_write(gREG,file_gREG,xnm="xc",ynm="yc",create=.TRUE.) 
    call nc_write(file_gREG,"Ts",  REG%Ts,  dim1="xc",dim2="yc",grid_mapping="crs") 
    call nc_write(file_gREG,"MB",  REG%MB,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"Hs",  REG%Hs,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"mask",REG%mask,dim1="xc",dim2="yc",grid_mapping="crs")

    ! Write new CCSM3 data to grid file 
    call grid_write(gCCSM3,file_gCCSM3b,xnm="lon",ynm="lat",create=.TRUE.)
    call nc_write(file_gCCSM3b,"Ts",  CCSM3b%Ts,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"MB",  CCSM3b%MB,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"Hs",  CCSM3b%Hs,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"mask",CCSM3b%mask,dim1="lon",dim2="lat")

    ! (as in Table 3 of Reerink et al, 2010)
    write(*,*) "=== Nearest-neighbor method ==="
    call grid_stats("Ts",CCSM3a%Ts,CCSM3b%Ts,CCSM3b%mask)
    call grid_stats("MB",CCSM3a%MB,CCSM3b%MB,CCSM3b%mask)
    call grid_stats("Hs",CCSM3a%Hs,CCSM3b%Hs,CCSM3b%mask)

    
    ! === SCRIP method ===

    file_gCCSM3b   = "output/ccsm3/grid_CCSM3-T42b_scrip.nc"
    file_gREG      = "output/ccsm3/grid_"//trim(REG%name)//"_scrip.nc"
    
    write(*,*) "=== SCRIP method ==="
    
    call map_scrip_init(mps1,gCCSM3,gREG,fldr="maps",load=.FALSE.,clean=.FALSE.,method="con")
    call map_scrip_init(mps2,gREG,gCCSM3,fldr="maps",load=.FALSE.,clean=.FALSE.,method="con")

    call map_scrip_field(mps1,"Ts",  CCSM3a%Ts,  REG%Ts,  REG%mask, method="mean")
    call map_scrip_field(mps1,"MB",  CCSM3a%MB,  REG%MB,  REG%mask, method="mean")
    call map_scrip_field(mps1,"Hs",  CCSM3a%Hs,  REG%Hs,  REG%mask, method="mean")

    CCSM3b%Ts = CCSM3a%Ts 
    CCSM3b%MB = CCSM3a%MB 
    CCSM3b%Hs = CCSM3a%Hs 
    call map_scrip_field(mps2,"Ts",  REG%Ts,  CCSM3b%Ts,  CCSM3b%mask, method="mean",reset=.FALSE.)
    call map_scrip_field(mps2,"MB",  REG%MB,  CCSM3b%MB,  CCSM3b%mask, method="mean",reset=.FALSE.)
    call map_scrip_field(mps2,"Hs",  REG%Hs,  CCSM3b%Hs,  CCSM3b%mask, method="mean",reset=.FALSE.)

    ! Write new regional data to grid file
    call grid_write(gREG,file_gREG,xnm="xc",ynm="yc",create=.TRUE.) 
    call nc_write(file_gREG,"Ts",  REG%Ts,  dim1="xc",dim2="yc",grid_mapping="crs") 
    call nc_write(file_gREG,"MB",  REG%MB,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"Hs",  REG%Hs,  dim1="xc",dim2="yc",grid_mapping="crs")
    call nc_write(file_gREG,"mask",REG%mask,dim1="xc",dim2="yc",grid_mapping="crs")

    ! Write new CCSM3 data to grid file 
    call grid_write(gCCSM3,file_gCCSM3b,xnm="lon",ynm="lat",create=.TRUE.)
    call nc_write(file_gCCSM3b,"Ts",  CCSM3b%Ts,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"MB",  CCSM3b%MB,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"Hs",  CCSM3b%Hs,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"mask",CCSM3b%mask,dim1="lon",dim2="lat")

    ! (as in Table 3 of Reerink et al, 2010)
    write(*,*) "=== SCRIP-mean method ==="
    call grid_stats("Ts",CCSM3a%Ts,CCSM3b%Ts,CCSM3b%mask)
    call grid_stats("MB",CCSM3a%MB,CCSM3b%MB,CCSM3b%mask)
    call grid_stats("Hs",CCSM3a%Hs,CCSM3b%Hs,CCSM3b%mask)

    
contains

    subroutine grid_stats(name,var1,var2,mask2)
        implicit none 

        character(len=*) :: name
        double precision :: var1(:,:), var2(:,:)
        integer          :: mask2(:,:)

        double precision, allocatable :: err(:)
        integer :: i, j, n

        double precision :: fld_ave, fld_range(2), MAE, AE_SD, RRD

        allocate(err(count(mask2 .eq. 1)))

        n = 0 
        do i = 1, size(var1,1)
            do j = 1, size(var1,2)
                if (mask2(i,j) .eq. 1) then 
                    n = n+1
                    err(n) = var2(i,j)-var1(i,j)
                end if 
            end do 
        end do 

        fld_ave = sum(var2*mask2) / dble(n) 
        fld_range(1) = minval(var2,mask2 .eq. 1)
        fld_range(2) = maxval(var2,mask2 .eq. 1)
        MAE = sum(dabs(err)) / dble(n)
        AE_SD = dsqrt( sum( (dabs(err) - (sum(dabs(err))/n))**2 ) / (n-1) )
        RRD   = MAE / (fld_range(2)-fld_range(1)) *100.d0

        write(*,"(a,3f10.1,3f12.4)") name, fld_range, fld_ave, MAE, AE_SD*2.d0, RRD

        return 

    end subroutine grid_stats

end program test_ccsm3