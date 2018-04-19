


program test_climber
    
    use coord 
    use ncio 
    
    use gaussian_filter 

    implicit none


    character(len=256) :: file_in, file_out 
    real(4), allocatable :: var_in(:,:), var_out(:,:), mask(:,:)  
    integer :: nx, ny, nz, k 
    real(4) :: dx 

    real, dimension(:, :), allocatable :: kernel

    file_in  = "output/interp/NH-40KM_present.nc"
    file_out = "output/interp/NH-40KM_present.nc"
    dx = 40.0 
    
    nx = nc_size(file_in,"xc")
    ny = nc_size(file_in,"yc")

    allocate(var_in(nx,ny))
    allocate(var_out(nx,ny))
    allocate(mask(nx,ny))

    ! ===Climate fields === 

    if (.TRUE.) then 

        call nc_read(file_in,"t2m_ann",var_in)
        var_out = var_in
        call filter_gaussian(var=var_out,sigma=240.0,dx=dx,mask=var_in.lt.0.0)
        call nc_write(file_out,"t2m_ann_gaussb",var_out,dim1="xc",dim2="yc")

        ! Equivalent to above
        var_out = var_in 
        do k = 1, 16
            call filter_gaussian(var=var_out,sigma=60.0,dx=dx,mask=var_in.lt.0.0)
        end do 
        call nc_write(file_out,"t2m_ann_gaussc",var_out,dim1="xc",dim2="yc")

    end if 

    if (.FALSE.) then 

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

    end if 

    ! === Ocean field ===

    if (.FALSE.) then 

    file_in  = "output/interp/NH-40KM_present_ocean.nc"
    file_out = "output/interp/NH-40KM_present_ocean.nc"

    nz = nc_size(file_in,"z_ocn")

    call gaussian_kernel(sigma=2.0,kernel=kernel)
    do k = 1, nz 
        call nc_read(file_in,"to_ann",var_in,start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
        call convolve(var_in, kernel, var_out)
        call nc_write(file_out,"to_ann_gauss_2",var_out,dim1="xc",dim2="yc",dim3="z_ocn", &
                      start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
    end do 

    call gaussian_kernel(sigma=5.0,kernel=kernel)
    do k = 1, nz 
        call nc_read(file_in,"to_ann",var_in,start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
        call convolve(var_in, kernel, var_out)
        call nc_write(file_out,"to_ann_gauss_5",var_out,dim1="xc",dim2="yc",dim3="z_ocn", &
                      start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
    end do 

    call gaussian_kernel(sigma=10.0,kernel=kernel)
    do k = 1, nz 
        call nc_read(file_in,"to_ann",var_in,start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
        call convolve(var_in, kernel, var_out)
        call nc_write(file_out,"to_ann_gauss_10",var_out,dim1="xc",dim2="yc",dim3="z_ocn", &
                      start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
    end do 

    call gaussian_kernel(sigma=20.0,kernel=kernel)
    do k = 1, nz 
        call nc_read(file_in,"to_ann",var_in,start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
        call convolve(var_in, kernel, var_out)
        call nc_write(file_out,"to_ann_gauss_20",var_out,dim1="xc",dim2="yc",dim3="z_ocn", &
                      start=[1,1,k],count=[size(var_in,1),size(var_in,2),1])
    end do 

    end if 



end program 

