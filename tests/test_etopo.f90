
program test_etopo

    use coord 
    use ncio 
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

    double precision :: tmplon(720), tmplat(360)
    character(len=256) :: gridname1, gridname0, outfldr 
    
    type akima_test_type 
        real(8), allocatable :: x(:), y(:), z(:) 
    end type 

    type(akima_test_type) :: ak_in, ak_out 
    integer :: i 

    write(*,*) 

    ! =======================================================================
    !
    ! Testing Akima routines
    !
    ! =======================================================================
    write(*,*) 
    write(*,*) " === Akima === "
    write(*,*) 

    allocate(ak_in%x(101),ak_in%y(101),ak_in%z(101))
    allocate(ak_out%x(20),ak_out%y(20),ak_out%z(20))
    
    call random_number(ak_in%x)
    call random_number(ak_in%y)
    ak_in%z = sin(ak_in%x*2.d0*3.14159)*cos(ak_in%y*2.d0*3.14159)

    write(*,*) "ak_in%x: ", minval(ak_in%x), maxval(ak_in%x)
    write(*,*) "ak_in%y: ", minval(ak_in%y), maxval(ak_in%y)
    write(*,*) "ak_in%z: ", minval(ak_in%z), maxval(ak_in%z)
    
    do i = 1, 20 
        ak_out%x(i) = dble(i-1) / 20.d0 
        ak_out%y(i) = dble(i-1) / 20.d0 
    end do 

    call random_number(ak_out%x)
    call random_number(ak_out%y)
    
    write(*,*) "ak_out%x: ", minval(ak_out%x), maxval(ak_out%x)
    write(*,*) "ak_out%y: ", minval(ak_out%y), maxval(ak_out%y)
    
    call interp2D_akima(ak_in%x,ak_in%y,ak_in%z, &
                        xout=ak_out%x,yout=ak_out%y,zout=ak_out%z)

    write(*,*) "ak_out%z: ", minval(ak_out%z), maxval(ak_out%z)
    
    
    call write_ascii("output/akima_in.txt",ak_in%x,ak_in%y,ak_in%z)
    call write_ascii("output/akima_out.txt",ak_out%x,ak_out%y,ak_out%z)

    stop 

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
        case("NH-320KM")
            call grid_init(grid1,name="NH-320KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=320.d0,nx=28,dy=320.d0,ny=26, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)
        
        case("NH-160KM")
            call grid_init(grid1,name="NH-160KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=160.d0,nx=56,dy=160.d0,ny=52, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)
        
        case("NH-80KM")
            call grid_init(grid1,name="NH-80KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=80.d0,nx=112,dy=80.d0,ny=104, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        case("NH-40KM")
            call grid_init(grid1,name="NH-40KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=40.d0,nx=224,dy=40.d0,ny=208, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        case("NH-20KM")
            call grid_init(grid1,name="NH-20KM",mtype="stereographic",units="kilometers", &
                           lon180=.TRUE.,dx=20.d0,nx=448,dy=20.d0,ny=416, &
                           lambda=-53.d0,phi=78.d0,alpha=32.7d0)

        case DEFAULT
            write(*,*) "test_etopo:: error: grid name not recognized: "//trim(gridname1)
            stop 

    end select

    ! Allocate data arrays
    call grid_allocate(grid1, vars1%zs)
    call grid_allocate(grid1, vars1%mask)

    ! Write grid definition to output file
    call grid_write(grid1,file1,xnm="xc",ynm="yc",create=.TRUE.)

!     stop 

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
    call grid_stats("zs",vars0%zs,vars0b%zs,vars0b%mask)

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

        if (count(mask2 .eq. 1) .gt. 0) then 
            fld_ave = sum(var2*mask2) / dble(n) 
            fld_range(1) = minval(var2,mask2 .eq. 1)
            fld_range(2) = maxval(var2,mask2 .eq. 1)
            MAE = sum(dabs(err)) / dble(n)
            AE_SD = dsqrt( sum( (dabs(err) - (sum(dabs(err))/n))**2 ) / (n-1) )
            RRD   = MAE / (fld_range(2)-fld_range(1)) *100.d0

            write(*,"(a,3f10.1,3f12.4)") name, fld_range, fld_ave, MAE, AE_SD*2.d0, RRD
        else
            write(*,*) "Interpolation mask empty - no interpolation performed?"
        end if 

        return 

    end subroutine grid_stats 

    subroutine write_ascii(filename,x,y,z)
        implicit none 

        character(len=*), intent(IN) :: filename 
        double precision, intent(IN) :: x(:), y(:), z(:) 
        integer :: i 

        open(10,file=trim(filename))
        write(10,"(3a12)") "x", "y", "z" 
        do i = 1, size(x)
            write(10,"(3f12.2)") x(i), y(i), z(i)
        end do 
        close(10)

        return 
    end subroutine write_ascii 

end program test_etopo

