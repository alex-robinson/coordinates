


program test_climber
    
    use coordinates
    use ncio 

    use gaussian_filter 

    implicit none


    character(len=256) :: file_in, file_out 
    real(4), allocatable :: var_in(:,:), var_out(:,:), mask(:,:)  
    integer :: nx, ny 

    real, dimension(:, :), allocatable :: kernel

    file_in  = "output/interp/NH-40KM_present.nc"
    file_out = "output/interp/NH-40KM_present.nc"

    nx = nc_size(file_in,"xc")
    ny = nc_size(file_in,"yc")

    allocate(var_in(nx,ny))
    allocate(var_out(nx,ny))
    allocate(mask(nx,ny))

    call nc_read(file_in,"t2m_ann",var_in)

    call gaussian_kernel(sigma=2.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"t2m_ann_gauss_2",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=5.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"t2m_ann_gauss_5",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=10.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"t2m_ann_gauss_10",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=20.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"t2m_ann_gauss_20",var_out,dim1="xc",dim2="yc")

    call nc_read(file_in,"pr_ann",var_in)

    call gaussian_kernel(sigma=2.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"pr_ann_gauss_2",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=5.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"pr_ann_gauss_5",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=10.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"pr_ann_gauss_10",var_out,dim1="xc",dim2="yc")

    call gaussian_kernel(sigma=20.0,kernel=kernel)
    call convolve(var_in, kernel, var_out)
    call nc_write(file_out,"pr_ann_gauss_20",var_out,dim1="xc",dim2="yc")

end program 

