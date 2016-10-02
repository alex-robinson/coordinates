

program test

    use interp1D
    use interp2D
    use interp_time
    use coordinates 
    use coordinates_mapping 
    use ncio 

    use gaussian_filter 

    implicit none

    integer :: t

    double precision, parameter :: missing_value = -9999.d0 

    type(grid_class) :: grid, gridlo, gridhi 
    double precision, dimension(:), allocatable   :: x, y
    double precision, dimension(:,:), allocatable :: var, varlo, varhi
    integer, dimension(:,:), allocatable :: mask, masklo, maskhi 
    character(len=256) :: outfldr,file_input, file_inputhi
    character(len=256) :: file_outlo, file_outhi
    integer :: i, j 

    real(4), dimension(:,:), allocatable :: vargauss
    real, dimension(:, :), allocatable :: kernel

    outfldr = "output/interp/"
    file_input   = trim(outfldr)//"GRL-50KM_TOPO.nc"
    file_inputhi = trim(outfldr)//"GRL-20KM_TOPO.nc"
    file_outhi   = trim(outfldr)//"GRL-20KM_TOPO1.nc"

    allocate(x(37),y(61))
    call nc_read(file_input,"xc",x)
    call nc_read(file_input,"yc",y)
    call grid_init(grid,name="GRL50KM",mtype="stereographic",units="km",lon180=.TRUE.,&
                   lambda=-40d0,phi=72.d0,alpha=7.5d0,x=x,y=y)

    deallocate(x)
    deallocate(y)
    allocate(x(90),y(150))
    call nc_read(file_inputhi,"xc",x)
    call nc_read(file_inputhi,"yc",y)
    call grid_init(gridhi,name="GRL20KM",mtype="stereographic",units="km",lon180=.TRUE., &
                   lambda=-40d0,phi=72.d0,alpha=7.5d0,x=x,y=y)

    ! Check if the grids are defined as the same map 
    write(*,*) "Same grid: ", compare_map(grid,gridhi)

    ! Allocate grid variables for testing 
    call grid_allocate(grid,var)
    call grid_allocate(gridhi,varhi)
    
    call grid_allocate(grid,mask)
    call grid_allocate(gridhi,maskhi)
    
    ! Load an example low resolution variable
    call nc_read(file_input,"zs",var)
    write(*,*) "zs range: ",minval(var), maxval(var)

    call nc_read(file_input,"mask",mask)
    write(*,*) "mask range: ",minval(mask), maxval(mask)

    ! Add some missing data 
    var(10:15,10:15) = missing_value
    var(30:37,:)     = missing_value 
    var(:,55:61)     = missing_value 

    ! Test bilinear interpolation
    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    
    call grid_write(gridhi,file_outhi,xnm="xc",ynm="yc",create=.TRUE.)
    call nc_write(file_outhi,"zs",varhi,dim1="xc",dim2="yc")

    ! Test missing value filling routines
    call fill_weighted(varhi,missing_value,n=1)
    call nc_write(file_outhi,"zs_filled1",varhi,dim1="xc",dim2="yc")

    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    call fill_weighted(varhi,missing_value,n=2)
    call nc_write(file_outhi,"zs_filled2",varhi,dim1="xc",dim2="yc")

    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    call fill_weighted(varhi,missing_value,n=3)
    call nc_write(file_outhi,"zs_filled3",varhi,dim1="xc",dim2="yc")

    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    call fill_weighted(varhi,missing_value,n=4)
    call nc_write(file_outhi,"zs_filled4",varhi,dim1="xc",dim2="yc")

    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    call fill_nearest(varhi,missing_value,n=4)
    call nc_write(file_outhi,"zs_nearest",varhi,dim1="xc",dim2="yc")

    ! Add some missing data 
    mask(10:15,10:15) = int(missing_value)
    mask(30:37,:)     = int(missing_value)
    mask(:,55:61)     = int(missing_value)

    maskhi = interp_nearest(grid%G%x,grid%G%y,mask,gridhi%G%x,gridhi%G%y,int(missing_value))
    call nc_write(file_outhi,"mask",maskhi,dim1="xc",dim2="yc")
    
    maskhi = interp_nearest(grid%G%x,grid%G%y,mask,gridhi%G%x,gridhi%G%y,nint(missing_value))
    call fill_nearest(maskhi,nint(missing_value),n=4)
    call nc_write(file_outhi,"mask_filled",maskhi,dim1="xc",dim2="yc")



    ! Test Gaussian smoothing

    ! Test bilinear interpolation
    varhi = interp_bilinear(grid%G%x,grid%G%y,var,gridhi%G%x,gridhi%G%y,missing_value)
    
    maskhi = 0
    where(varhi.eq.missing_value) maskhi = 1

!     call run_gaussian_filter(sigma=2.0, truncate=4.0, kx=5, ky=5, kernel, &
!                              nx, ny, input, output, mask, has_mask)
    
    call grid_allocate(gridhi,vargauss)
    
    call gaussian_kernel(sigma=0.1,kernel=kernel)
    call convolve(real(varhi), kernel, vargauss, mask=maskhi.eq.1)

    call nc_write(file_outhi,"zs_gauss",vargauss,dim1="xc",dim2="yc")

end program test