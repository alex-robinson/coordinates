
program test_etopo

    use ncio 
    use coordinates
    use mod_toms526 

    implicit none

    type(grid_class) :: grid0, grid1 
    type(map_class)  :: map_g0g1, map_g1g0 

    type vars_type 
        character (len=256) :: name 
        integer             :: nx, ny 
        double precision    :: lambda, phi, alpha
        double precision, dimension(:,:), allocatable :: zs
        integer,          dimension(:,:), allocatable :: mask
    end type 

    type(vars_type) :: vars0, vars1, vars0b

    character(len=256) :: file0, file1, file0b 

    double precision :: tmplon(721), tmplat(361)
    character(len=256) :: gridname1, gridname0, outfldr 
    
    write(*,*) 

    ! =========================================================
    !
    ! USER DEFINITIONS
    !
    ! =========================================================
 
    gridname1 = "NH-40KM"
    gridname0 = "GLOBAL-050deg"
    outfldr   = "output/"


    ! Generate file names 
    file0  = "data/ETOPO1-ICE-050deg.nc"
    file1  = trim(outfldr)//trim(gridname1)//"_zs.nc"
    file0b = trim(outfldr)//trim(gridname0)//"_zs.nc"
    
    ! =========================================================
    !
    ! OUTPUT GRID DEFINITION
    !
    ! =========================================================

    select case(trim(gridname1))
        case("NH-40KM")
            call grid_init(grid1,name="NH-40KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=281,dy=40.d0,ny=281, &
                           lambda=0.d0,phi=90.d0,alpha=44.7d0)

        case("NH-20KM")
            call grid_init(grid1,name="NH-20KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=20.d0,nx=561,dy=20.d0,ny=561, &
                           lambda=0.d0,phi=90.d0,alpha=44.7d0)

        case DEFAULT
            write(*,*) "test_etopo:: error: grid name not recognized: "//trim(gridname1)
            stop 

    end select

    ! Allocate data arrays
    call grid_allocate(grid1, vars1%zs)
    call grid_allocate(grid1, vars1%mask)

    ! Write grid definition to output file
    call grid_write(grid1,file1,xnm="xc",ynm="yc",create=.TRUE.)

    ! =======================================================================
    !
    ! Define global input grid and load input data
    !
    ! =======================================================================

    ! Define input latlon grid
    call nc_read(file0,"lon",tmplon)
    call nc_read(file0,"lat",tmplat)
    call grid_init(grid0,name=gridname0,mtype="latlon",units="degrees",lon180=.TRUE., &
                     x=tmplon,y=tmplat)

    ! Allocate arrays of the size of the global grid
    call grid_allocate(grid0, vars0%zs)
    call grid_allocate(grid0, vars0%mask)

    ! Allocate arrays to check reinterpolation error 
    call grid_allocate(grid0, vars0b%zs)
    call grid_allocate(grid0, vars0b%mask)

    ! Load original GCM data
    call nc_read(file0,"z",vars0%zs)
    vars0b%zs = vars0%zs 

    ! Write reinterpolation grid to second output file
    call grid_write(grid0,file0b,xnm="xc",ynm="yc",create=.TRUE.)

    ! =======================================================================
    !
    ! Map the fields 
    !
    ! =======================================================================
    write(*,*) 
    write(*,*) " === MAPPING === "
    write(*,*) 

    ! Map each field to the regional domain using the quadrant method (no max_distance required here)
    call map_init(map_g0g1,grid0,grid1,max_neighbors=4,lat_lim=2.0d0,fldr="maps",load=.TRUE.)
    call map_field(map_g0g1,"zs",vars0%zs,vars1%zs,vars1%mask,method="quadrant")
    call nc_write(file1,"zs",vars1%zs,dim1="xc",dim2="yc")
    call nc_write(file1,"mask",vars1%mask,dim1="xc",dim2="yc")

    ! Map each field back to the original domain using the radius method
    call map_init(map_g1g0,grid1,grid0,max_neighbors=6,lat_lim=2.0d0,fldr="maps",load=.TRUE.)
    call map_field(map_g1g0,"zs",vars1%zs,vars0b%zs,vars0b%mask,"shepard",125.d3,fill=.FALSE.)
    call nc_write(file0b,"zs",  vars0b%zs,dim1="xc",dim2="yc")
    call nc_write(file0b,"mask",vars0b%mask,dim1="xc",dim2="yc")

    ! Calculate statistics concerning remapping
    ! (as in Table 3 of Reerink et al, 2010)
    call grid_stats("zs",vars0%zs,vars0b%zs,vars0%mask)


    stop 

    ! =======================================================================
    !
    ! Testing Akima routines
    !
    ! =======================================================================
!     write(*,*) 
!     write(*,*) " === Akima === "
!     write(*,*) 

!     write(*,*) "in  x: ", minval(gCCSM3%x),   maxval(gCCSM3%x) 
!     write(*,*) "in  y: ", minval(gCCSM3%y),   maxval(gCCSM3%y) 
!     write(*,*) "in  z: ", minval(CCSM3a%Ts), maxval(CCSM3a%Ts) 
!     write(*,*) "out x: ", minval(gREG%lon),   maxval(gREG%lon) 
!     write(*,*) "out y: ", minval(gREG%lat),   maxval(gREG%lat) 
    
!     call interp2D_akima(gCCSM3%x,gCCSM3%y,CCSM3a%Ts,xout=gREG%lon,yout=gREG%lat,zout=REG%Ts)

!     write(*,*) "out z: ", minval(REG%Ts), maxval(REG%Ts) 
    
!     call nc_write(file_gREG,"Ts_akima",  REG%Ts,  dim1="xc",dim2="yc")

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

end program test_etopo

