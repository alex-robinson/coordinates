
!! TO COMPILE : 
!! gfortran -fcheck=all -I/opt/local/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!! or
!! ifort -g -I/home/robinson/apps/netcdf/netcdf/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

program test_ccsm3_diffusion

    use coord 
    use ncio 
    use gaussian_filter 

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

    type gauss_type 
        real(4), allocatable :: input(:,:), output(:,:), grad(:,:), kappa(:,:)
        logical, allocatable :: mask(:,:)
        real, dimension(:, :), allocatable :: kernel
        integer :: ntot, nbin, nnow 
        real(4) :: sigmak, sigma, grad_now 

    end type 

    type(gauss_type) :: gtest 
    integer :: i 
    integer :: nt 
    real(4) :: dt, kappa, gradmax  

    ! =======================================================================
    !
    ! Step 1: Define global input grid and load data that will be used here
    !
    ! =======================================================================

    ! Define file names for input and output of global grids
    file_input     = "data/ccsm_example_dec_feb_pd.nc"
    file_gCCSM3a   = "output/ccsm3diff/grid_CCSM3-T42a.nc"
    file_gCCSM3b   = "output/ccsm3diff/grid_CCSM3-T42b_nn.nc"
    
    ! CCSM3 T42 latlon grid
    call nc_read(file_input,"lon",tmplon)
    call nc_read(file_input,"lat",tmplat)
    call grid_init(gCCSM3,name="CCSM3-T42",mtype="latlon",units="degrees",lon180=.FALSE., &
                     x=tmplon,y=tmplat)
!     call grid_init(gCCSM3,filename="maps/grid_CCSM3-T42.txt",x=tmplon,y=tmplat)

    ! Write grid information to files
    call grid_write(gCCSM3,file_gCCSM3a,xnm="lon",ynm="lat",create=.TRUE.)
    call grid_write(gCCSM3,file_gCCSM3b,xnm="lon",ynm="lat",create=.TRUE.)

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

!     REG%name = "ANT-20KM"
!     REG%nx   = 281
!     REG%ny   = 281 
!     REG%lambda =   0.d0 
!     REG%phi    = -90.d0
!     REG%alpha  =  19.d0 

    REG%name = "GRL-20KM"
    REG%nx   = 76
    REG%ny   = 151 
    REG%lambda = 320.d0 
    REG%phi    = 72.d0
    REG%alpha  = 7.5d0

!     REG%name = "HIM-20KM"
!     REG%nx   = 200
!     REG%ny   = 200 
!     REG%lambda = 90.d0 
!     REG%phi    = 32.d0
!     REG%alpha  = 14.5d0

    ! Oblique stereographic grid for region of interest
    call grid_init(gREG,name=REG%name,mtype="stereographic",units="kilometers",lon180=.FALSE., &
                     dx=20.d0,nx=REG%nx,dy=20.d0,ny=REG%ny, &
                     lambda=REG%lambda,phi=REG%phi,alpha=REG%alpha)

    ! Output grid to file
    file_gREG = "output/ccsm3diff/grid_"//trim(REG%name)//"_nn.nc"
    call grid_write(gREG,file_gREG,xnm="xc",ynm="yc",create=.TRUE.)

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
    call map_init(mCCSM3_REG,gCCSM3,gREG,max_neighbors=10,lat_lim=4.0d0,fldr="maps",load=.TRUE.)
    call map_init(mREG_CCSM3,gREG,gCCSM3,max_neighbors=10,lat_lim=2.0d0,fldr="maps",load=.TRUE.)

    ! Map each field to the regional domain using the quadrant method (no max_distance required here)
    call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="nn")
    call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="nn")
    call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="nn")
!     call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="nng",sigma=80.d0)
!     call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="nng",sigma=80.d0)
!     call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="nng",sigma=80.d0)

    ! Write new regional data to grid file 
    call nc_write(file_gREG,"Ts",  REG%Ts,  dim1="xc",dim2="yc")
    call nc_write(file_gREG,"MB",  REG%MB,  dim1="xc",dim2="yc")
    call nc_write(file_gREG,"Hs",  REG%Hs,  dim1="xc",dim2="yc")
    call nc_write(file_gREG,"mask",REG%mask,dim1="xc",dim2="yc")

    ! ======================================================================
    ! ======================================================================

    ! Testing variable sigmas for Gaussian filter 
    call grid_allocate(gREG,gtest%input)
    call grid_allocate(gREG,gtest%output)
    call grid_allocate(gREG,gtest%mask)
    call grid_allocate(gREG,gtest%grad)
    call grid_allocate(gREG,gtest%kappa)
            
    ! Map field with nearest neighbor algorithm
    call map_field(mCCSM3_REG,"Ts",CCSM3a%Ts,REG%Ts,method="nn")
    gtest%input = real(REG%Ts)
    gtest%grad = hgrad(gtest%input,dx=real(gREG%G%dx),norm=.TRUE.)
    call nc_write(file_gREG,"dTs",  gtest%grad,  dim1="xc",dim2="yc")

    call map_field(mCCSM3_REG,"MB",CCSM3a%MB,REG%MB,method="nn")
    gtest%input = real(REG%MB)
    gtest%grad = hgrad(gtest%input,dx=real(gREG%G%dx),norm=.TRUE.)
    call nc_write(file_gREG,"dMB",  gtest%grad,  dim1="xc",dim2="yc")

    call map_field(mCCSM3_REG,"Hs",CCSM3a%Hs,REG%Hs,method="nn")
    gtest%input = real(REG%Hs)
    gtest%grad = hgrad(gtest%input,dx=real(gREG%G%dx),norm=.TRUE.)
!     call nc_write(file_gREG,"dHs",  gtest%grad,  dim1="xc",dim2="yc")


    ! Output grid to file
    dt = 1.0 
    nt = 1000 
    kappa = 30.0 

    write(*,*) "Timestep = ", dt 
    write(*,*) "Timestep max = ", &
            diff2D_timestep(real(gREG%G%dx*gREG%xy_conv),real(gREG%G%dy*gREG%xy_conv),kappa)

    file_new = "output/ccsm3diff/"//trim(REG%name)//"_diffuse.nc"
    call grid_write(gREG,file_new,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write_dim(file_new,"time",x=0.0,dx=dt,nx=1,units="seconds",unlimited=.TRUE.)

    call diffuse2D(gtest%input,gtest%grad,dx=real(gREG%G%dx))
    write(*,*) "stdev = ", stdev_2D(gtest%grad,mv=-9999.0)
    gradmax = maxval(abs(gtest%grad))
    write(*,*) "gradmax = ", gradmax 

    do i = 1, nt 
        call diffuse2D(gtest%input,gtest%grad,dx=real(gREG%G%dx))
!         gtest%kappa = kappa*(abs(gtest%grad)/maxval(abs(gtest%grad)))
        gtest%kappa = kappa 
        gradmax     = maxval(abs(gtest%grad))

        write(*,"(i6,3g15.3)") i, minval(gtest%grad), maxval(gtest%grad),stdev_2D(gtest%grad,mv=-9999.0)
        gtest%output = gtest%input 

        if ( i-1 .gt. 0) then 
            where(abs(gtest%grad) .gt. gradmax*0.85) &
                gtest%output = gtest%input + dt*gtest%kappa*gtest%grad 
        end if 

        call nc_write(file_new,"time",i*dt,dim1="time",start=[i],count=[1])
        call nc_write(file_new,"dHs2",  gtest%grad,  dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,i],count=[gREG%G%nx,gREG%G%ny,1])
        call nc_write(file_new,"Hs2",  gtest%output,  dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,i],count=[gREG%G%nx,gREG%G%ny,1])
        call nc_write(file_new,"Hs2_Hs",  gtest%output-REG%Hs,  dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,i],count=[gREG%G%nx,gREG%G%ny,1])
        gtest%input = gtest%output 

        if (gradmax .lt. 0.02) exit 
    end do 

    stop 

    gtest%sigma  = 200.d0 
    gtest%sigmak =  10.d0 
    gtest%ntot   = (gtest%sigma / gtest%sigmak)**2
    gtest%nbin   = gtest%ntot/10 

    call gaussian_kernel(sigma=gtest%sigmak/real(gREG%G%dx),kernel=gtest%kernel)

    gtest%nnow = 9
    gtest%grad_now = maxval(gtest%grad)*gtest%nnow/gtest%nbin 

    do i = 1, gtest%ntot 

        if (mod(i,gtest%nbin) .eq. 0) then 
            gtest%nnow = max(0,gtest%nnow-1)
            gtest%grad_now = maxval(gtest%grad)* gtest%nnow / gtest%nbin 
        end if
        gtest%mask = .FALSE.
        where(gtest%grad .ge. gtest%grad_now) gtest%mask = .TRUE. 

        write(*,*) i, mod(i,gtest%nbin), gtest%nnow, gtest%grad_now, count(gtest%mask)

        call convolve(gtest%input, gtest%kernel, gtest%output, mask=gtest%mask)
        gtest%input = gtest%output 

    end do 

    call nc_write(file_gREG,"Hs_sigma100", gtest%output,  dim1="xc",dim2="yc")

    gtest%grad = hgrad(gtest%output,dx=real(gREG%G%dx),norm=.TRUE.)
    call nc_write(file_gREG,"dHs_sigma100", gtest%grad,  dim1="xc",dim2="yc")

    stop 

    ! ======================================================================
    ! ======================================================================


    ! Map each field back to the CCSM3 domain using the radius method
    call map_field(mREG_CCSM3,"Ts",REG%Ts,CCSM3b%Ts,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"MB",REG%MB,CCSM3b%MB,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
    call map_field(mREG_CCSM3,"Hs",REG%Hs,CCSM3b%Hs,CCSM3b%mask,"shepard",125.d3,fill=.FALSE.)
!     call map_field(mREG_CCSM3,"Ts",REG%Ts,CCSM3b%Ts,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)
!     call map_field(mREG_CCSM3,"MB",REG%MB,CCSM3b%MB,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)
!     call map_field(mREG_CCSM3,"Hs",REG%Hs,CCSM3b%Hs,CCSM3b%mask,"nn",125.d3,fill=.FALSE.)

    ! Write new CCSM3 data to grid file 
    call nc_write(file_gCCSM3b,"Ts",  CCSM3b%Ts,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"MB",  CCSM3b%MB,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"Hs",  CCSM3b%Hs,  dim1="lon",dim2="lat")
    call nc_write(file_gCCSM3b,"mask",CCSM3b%mask,dim1="lon",dim2="lat")

    ! Calculate statistics concerning remapping
    ! (as in Table 3 of Reerink et al, 2010)
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


    function stdev_2D(x,mv) result(sdx)
        
        implicit none 

        real(4) :: x(:,:), mv 
        real(4) :: sdx 
        real(4) :: ave, totsq 
        integer :: ntot 

        ntot = count(x.ne.mv)

        sdx   = -1.0
        if (ntot .gt.2) then 
            ave   = sum(x,mask=x.ne.mv)/max(1,ntot)
            totsq = sum( (x-ave)*(x-ave) )
            sdx   = sqrt(totsq/(ntot-1))
        end if 

        return

    end function stdev_2D 

    function stdev_vec(x,mv) result(sdx)
        
        implicit none 

        real(4) :: x(:), mv 
        real(4) :: sdx 
        real(4) :: ave, totsq 
        integer :: ntot

        ntot = count(x.ne.mv)

        sdx   = -1.0
        if (ntot .gt.2) then
            ave   = sum(x,mask=x.ne.mv)/max(1,ntot)
            totsq = sum( (x-ave)*(x-ave) )
            sdx   = sqrt(totsq/(ntot-1))
        end if 

        return

    end function stdev_vec 


end program test_ccsm3_diffusion