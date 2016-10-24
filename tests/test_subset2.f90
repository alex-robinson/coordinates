
!! TO COMPILE : 
!! gfortran -fcheck=all -I/opt/local/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/opt/local/lib -lnetcdff -lnetcdf
!! or
!! ifort -g -I/home/robinson/apps/netcdf/netcdf/include -o test_ccsm3.x ../ncio/ncio3.f90 geodesic.f90 planet.f90 projection_oblimap2.f90 coordinates.f90 test_ccsm3.f90 -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

program test_subset

    use coord 
    use ncio 
    
    implicit none

    type(grid_class)   :: grid, gridlo
    type(points_class) :: ptshi
    type(subset_class) :: sub 

    type(map_class)  :: mgridhi_gridlo, mgridlo_gridhi

    type vars_type 
        character (len=256) :: name, file 
        integer             :: nx, ny 
        double precision    :: lambda, phi, alpha
        double precision, dimension(:), allocatable  :: x, y 
        double precision, dimension(:,:), allocatable :: lon, lat, zs, dzs, diff
        integer,          dimension(:,:), allocatable :: mask
    end type 

    type vars1D_type 
        double precision, dimension(:), allocatable :: zs, dzs
        integer,          dimension(:), allocatable :: mask
    end type 

    type(vars_type)   :: set, setlo, set1, set2
    type(vars1D_type) :: sethi 

    logical, dimension(:,:), allocatable :: mask_pack 

    character(len=256) :: file_input, outfldr, file_outhi, file_outlo 
    character(len=256) :: file_out1, file_out2
    integer :: i, j 
    logical :: load 

    double precision, parameter :: missing_value = -9999.d0 

    outfldr = "output/subsetting2/"
    file_input = trim(outfldr)//"GRL-20KM_TOPO.nc"
    file_outhi = trim(outfldr)//"GRL-20KM_TOPO_hi.nc"
    file_outlo = trim(outfldr)//"GRL-20KM_TOPO_lo.nc"
    file_out1  = trim(outfldr)//"GRL-20KM_TOPO1.nc" 
    file_out2  = trim(outfldr)//"GRL-20KM_TOPO2.nc" 
        
    ! Load stored maps from memory?
    load = .TRUE. 

    ! =======================================================================
    !
    ! Step 1: Define main grid and load data that will be used here
    !
    ! =======================================================================

    set%nx     = nc_size(file_input,"xc")
    set%ny     = nc_size(file_input,"yc")
    set%lambda = 320.d0
    set%phi    =  72.d0
    set%alpha  =   7.5d0 

    allocate(set%x(set%nx))
    allocate(set%y(set%ny))
    call nc_read(file_input,"xc",set%x)
    call nc_read(file_input,"yc",set%y)
    call grid_init(grid,name="GRL20KM",mtype="stereographic",units="km", &
                   lambda=set%lambda,phi=set%phi,alpha=set%alpha,x=set%x,y=set%y)

    call grid_allocate(grid, set%lon)
    call grid_allocate(grid, set%lat)
    call grid_allocate(grid, set%zs)
    call grid_allocate(grid, set%dzs)
    call grid_allocate(grid, set%mask)
    call grid_allocate(grid, set%diff)

    call nc_read(file_input,"lon2D",set%lon)
    call nc_read(file_input,"lat2D",set%lat)
    call nc_read(file_input,"zs",set%zs)
    call nc_read(file_input,"mask",set%mask)

    ! Get horizontal gradient of elevation
    set%dzs = 0.d0 
    do i = 2,grid%G%nx-1
    do j = 2,grid%G%ny-1  
        set%dzs(i,j) = dabs(0.5d0*(set%zs(i+1,j)-set%zs(i-1,j))/grid%G%dx &
                            + 0.5d0*(set%zs(i,j+1)-set%zs(i,j-1))/grid%G%dy)
    end do 
    end do 


    ! Output our basic diagnostic file too 
    set1 = set 
    call grid_write(grid,file_out1,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write(file_out1,"zs",set1%zs,dim1="xc",dim2="yc")
    call nc_write(file_out1,"dzs",set1%dzs,dim1="xc",dim2="yc")
    call nc_write(file_out1,"zs2",(set1%zs*1e-3)**3,dim1="xc",dim2="yc")

    ! =======================================================================
    !
    ! Step 2: Define subset and map gridded data
    !
    ! =======================================================================

    call subset_init(sub,grid,factor=3,npts=grid%npts*3,load=load)
!     call subset_init(sub,grid,factor=2,npts=0,load=load)
!     call subset_init(sub,grid,factor=1,npts=1000,load=load)
    
    ! Specify the local points object with subset information (redundant, but convenient)
    ptshi = sub%pts 
    
    ! Generate a mask for the subset of points based on elevation gradient
    ! and redefine the subset based on this mask 
    call subset_gen_mask(sub%mask_pack,set%dzs,ptshi%npts, &
                         min_spacing=floor(40.d0/sub%grid%G%dx),method="max",map=sub%map_tosub_grid)
    write(*,*) "ptshi%npts =",ptshi%npts 
    call subset_redefine(sub,grid,sub%mask_pack,load=load)
    ptshi = sub%pts 

    ! Check the mask
    call grid_write(sub%grid,"output/subsetting2/mask_pack.nc",xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write("output/subsetting2/mask_pack.nc","mask_pack",sub%mask_pack,dim1="xc",dim2="yc") 

    ! Map the original gridded data to the 1D subset
    call points_allocate(ptshi, sethi%zs,missing_value)
    call points_allocate(ptshi, sethi%dzs,missing_value)
    call points_allocate(ptshi, sethi%mask,nint(missing_value))

    call subset_to_points(sub,set%zs,sethi%zs,sub%mask_pack,method="radius")
    call subset_to_points(sub,set%dzs,sethi%dzs,sub%mask_pack,method="radius")
    call subset_to_points(sub,set%mask,sethi%mask,sub%mask_pack)

    call points_write(ptshi,file_outhi,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write(file_outhi,"zs",sethi%zs,dim1="point")
    call nc_write(file_outhi,"dzs",sethi%dzs,dim1="point")
    call nc_write(file_outhi,"zs2",(sethi%zs*1e-3)**3,dim1="point")
    call nc_write(file_outhi,"mask",sethi%mask,dim1="point")
    
    ! Map the calculated data back to the original grid 
    call subset_to_grid(sub,sethi%zs,set1%zs,sub%mask_pack,method="radius")
    call subset_to_grid(sub,sethi%dzs,set1%dzs,sub%mask_pack,method="radius")
    call subset_to_grid(sub,sethi%mask,set1%mask,sub%mask_pack)

    call nc_write(file_out1,"zs",set1%zs,dim1="xc",dim2="yc")
    call nc_write(file_out1,"dzs",set1%dzs,dim1="xc",dim2="yc")
    call nc_write(file_out1,"mask",set1%mask,dim1="xc",dim2="yc")

    set1%diff = missing_value
    where(.not. set1%zs .eq. missing_value) set1%diff = set1%zs-set%zs
    call nc_write(file_out1,"zs_diff",set1%diff,dim1="xc",dim2="yc")

    ! Check how well filling works
    call subset_to_grid(sub,sethi%zs,set1%zs,sub%mask_pack,method="radius",fill=.TRUE.)
    set1%diff = missing_value
    where(.not. set1%zs .eq. missing_value) set1%diff = set1%zs-set%zs
    call nc_write(file_out1,"zs_diff_fill",set1%diff,dim1="xc",dim2="yc")

    ! =======================================================================
    !
    ! Step 3: Map pts between another grid
    !
    ! =======================================================================

    ! Initialize the low resolution grid
    call grid_init(gridlo,name="GRL50KM",mtype=grid%mtype,units=grid%units, &
                   lambda=grid%proj%lambda,phi=grid%proj%phi,alpha=grid%proj%alpha, &
                   x0=grid%G%x(1),dx=50.d0,nx=floor((grid%G%nx-1)*(grid%G%dx/50.d0)+1), &
                   y0=grid%G%y(1),dy=50.d0,ny=floor((grid%G%ny-1)*(grid%G%dy/50.d0)+1))
    call grid_write(gridlo,file_outlo,xnm="xc",ynm="yc",create=.TRUE.)  

    ! Initialize maps to/from low resolution grid
    call map_init(mgridhi_gridlo,sub%pts,gridlo,max_neighbors=9,lat_lim=2.d0,load=load)
    call map_init(mgridlo_gridhi,gridlo,sub%pts,max_neighbors=6,lat_lim=2.d0,load=load)
    
    ! Map the calculated data to the low resolution grid 
    call grid_allocate(gridlo,setlo%zs,  missing_value)
    call grid_allocate(gridlo,setlo%dzs, missing_value)
    call grid_allocate(gridlo,setlo%mask,nint(missing_value))

    do i = 1, 100 ! Loop for performance testing 

    call subset_to_grid(sub,sethi%zs,  setlo%zs,  sub%mask_pack,map=mgridhi_gridlo,method="radius",border=.TRUE.)
    call subset_to_grid(sub,sethi%dzs, setlo%dzs, sub%mask_pack,map=mgridhi_gridlo,method="radius",border=.TRUE.)
    call subset_to_grid(sub,sethi%mask,setlo%mask,sub%mask_pack,map=mgridhi_gridlo,border=.TRUE.)

    ! Write low resolution data
    write(*,*) "min/max zs: ",minval(setlo%zs),maxval(setlo%zs)
    call nc_write(file_outlo,"zs",  setlo%zs,  dim1="xc",dim2="yc",missing_value=missing_value)
    call nc_write(file_outlo,"dzs", setlo%dzs, dim1="xc",dim2="yc",missing_value=missing_value)
    call nc_write(file_outlo,"mask",setlo%mask,dim1="xc",dim2="yc",missing_value=nint(missing_value))

    ! Remap to hi resolution
    call subset_to_points(sub,setlo%zs,sethi%zs,sub%mask_pack,map=mgridlo_gridhi,method="radius",border=.TRUE.)
    call subset_to_points(sub,setlo%dzs,sethi%dzs,sub%mask_pack,map=mgridlo_gridhi,method="radius",border=.TRUE.)
    call subset_to_points(sub,setlo%mask,sethi%mask,sub%mask_pack,map=mgridlo_gridhi,border=.TRUE.)

    ! Remap to original grid
    set2 = set 
    call subset_to_grid(sub,sethi%zs,  set2%zs,  sub%mask_pack,method="radius")
    call subset_to_grid(sub,sethi%dzs, set2%dzs, sub%mask_pack,method="radius")
    call subset_to_grid(sub,sethi%mask,set2%mask,sub%mask_pack)

    end do ! End performance testing loop 

    ! Output our basic diagnostic file too 
    call grid_write(grid,file_out2,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write(file_out2,"zs",set2%zs,dim1="xc",dim2="yc")
    call nc_write(file_out2,"dzs",set2%dzs,dim1="xc",dim2="yc")
    call nc_write(file_out2,"mask",set2%mask,dim1="xc",dim2="yc")

contains

    subroutine gen_subset_mask(mask_pack,dzs,dx,dx_min,npts,map)

        implicit none 

        logical, dimension(:,:), intent(INOUT) :: mask_pack 
        double precision, dimension(:,:), intent(IN) :: dzs 
        type(map_class), intent(IN), optional :: map
        double precision :: dx, dx_min
        integer :: npts 

        integer :: nx, ny, by, npts0, loc(2), i,j
        double precision, dimension(:,:), allocatable :: dzs1
        integer, dimension(:,:), allocatable :: mask1 
        logical :: coarse 
        double precision :: dzs_lim, dzs_step  

        nx = size(mask_pack,1)
        ny = size(mask_pack,2)

        if (size(mask_pack) .eq. size(dzs)) then 
            mask_pack = .TRUE. 

        else 
        
            allocate(dzs1(nx,ny),mask1(nx,ny))

            call map_field(map,"Subset mask",dzs,dzs1,mask1,method="radius",fill=.TRUE.)

            ! Initially set all points to .FALSE.
            mask_pack = .FALSE.

            ! Fill in mask_pack evenly-spaced at desired minimum resolution
            if (dx_min .gt. 0.d0) then 
                by = floor( dx_min / dx ) 
                do i = 1,nx,by
                do j = 1,ny,by 
                    mask_pack(i,j) = .TRUE. 
                end do 
                end do
            end if 

            ! Now fill in the remainder of the mask with the highest elevation gradients
            npts0 = count(mask_pack)
            coarse = .TRUE.
            dzs_lim  = maxval(dzs1)
            dzs_step = dzs_lim*0.05 

            do while(npts0 < npts)
                where(mask_pack) dzs1 = 0.d0       ! Eliminate points that are already included
                
                if (coarse) then 
                    dzs_lim = dzs_lim - dzs_step 
                    if (count(dzs1 .ge. dzs_lim)+count(mask_pack) .lt. npts) then 
                        where (dzs1 .ge. dzs_lim) mask_pack = .TRUE. 
                    else
                        coarse = .FALSE. 
                    end if 
                else 
                    loc = maxloc(dzs1)                 ! Determine location of highest slope
                    mask_pack(loc(1),loc(2)) = .TRUE. ! Add this point to mask 
                end if 
                npts0 = count(mask_pack)          ! Recount masked points
                !write(*,*) "npts0 =", npts0 
            end do 

        end if 

        write(*,*) "mask_pack total = ",count(mask_pack)," / " , (nx*ny)

        return 

    end subroutine gen_subset_mask 


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

end program test_subset