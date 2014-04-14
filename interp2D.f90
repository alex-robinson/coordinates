
module interp2D
    
    implicit none 

    !! real(dp) definition and some internal constants
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    private
    public :: interp_bilinear

contains


    function interp_bilinear(x,y,z,xout,yout,missing_value,mask) result(zout)
        ! Find closest x-indices and closest y-indices on original
        ! grid (assume these and next indices will bracket our point)
        ! Perform weighted interpolation 
        implicit none 

        real(dp), dimension(:) :: x, y, xout, yout 
        real(dp), dimension(:,:) :: z
        real(dp), optional :: missing_value 
        logical, dimension(:,:), optional :: mask 
        real(dp), dimension(size(xout,1),size(yout,1)) :: zout 
        logical,  dimension(size(xout,1),size(yout,1)) :: mask_interp 
        real(dp) :: missing_val 

        integer, dimension(size(xout,1)) :: x_idx
        integer, dimension(size(yout,1)) :: y_idx
        
        integer :: i, j, i1, j1  
        integer :: nx, ny, nx1, ny1 
        double precision :: alpha1, alpha2, p0, p1 

        nx = size(x,1)
        ny = size(y,1)

        nx1 = size(xout,1)
        ny1 = size(yout,1)

        ! Determine which points we are interested in interpolating
        mask_interp = .TRUE. 
        if (present(mask)) mask_interp = mask 

        ! Determine missing value if present 
        missing_val = MISSING_VALUE_DEFAULT
        if (present(missing_value)) missing_val = missing_value

        ! Get x-indices corresponding to nearest neighbor
        ! greater-than-equal-to x-value of interest
        do i1 = 1, nx1 

            if (xout(i1) .le. x(1)) then 
                x_idx(i1) = -1

            else if (xout(i1) .ge. (x(nx))) then 
                x_idx(i1) = -2 
            else 

                do i = 1, nx 
                    if (x(i) .ge. xout(i1)) exit 
                end do 

                x_idx(i1) = i 
            end if 

        end do 

        ! Get y-indices corresponding to nearest neighbor
        ! greater-than-equal-to y-value of interest
        do j1 = 1, ny1 

            if (yout(j1) .le. y(1)) then 
                y_idx(j1) = -1

            else if (yout(j1) .ge. (y(ny))) then 
                y_idx(j1) = -2 
            else 

                do j = 1, ny 
                    if (y(j) .ge. yout(j1)) exit 
                end do 

                y_idx(j1) = j 
            end if 

        end do 


        ! Now loop over output grid points and perform
        ! bilinear interpolation where desired 
        zout = missing_val 
        do i1 = 1, nx1 
        do j1 = 1, ny1 

            ! Only interpolate points of interest 
            if (mask_interp(i1,j1)) then 

                i = x_idx(i1)
                j = y_idx(j1) 

                ! Only interpolate points inside the original grid (for now)
                if (i .gt. 0 .and. i-1 .gt. 0 .and. j .gt. 0 .and. j-1 .gt. 0) then 
                    
                    alpha1 = (xout(i1) - x(i-1)) / (x(i)-x(i-1))
                    p0 = z(i-1,j-1) + alpha1*(z(i,j-1)-z(i-1,j-1))
                    p1 = z(i-1,j)   + alpha1*(z(i,j)-z(i-1,j))
                    
                    alpha2 = (yout(j1) - y(j-1)) / (y(j)-y(j-1))
                    zout(i1,j1) = p0 + alpha2*(p1-p0)

                end if 

            end if 

        end do 
        end do
        

        return 

    end function interp_bilinear




    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine :  t o h i r e s
    ! Author     :  Alex Robinson
    ! Purpose    :  interpolate a cartesian grid to a higher resolution
    !               Note: assumes lower res is a multiple of higher res!
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    subroutine tohires(in,out,ratio,m)

        implicit none

        real(dp) :: in(:,:), out(:,:)
        integer :: nxl, nyl, nx, ny 
        integer :: nin, i, j, inow, jnow, ih, jh, jh2, q, m
        real(dp) :: ratio, f, frac

        nxl = size(in,1)
        nyl = size(in,2)
        nx  = size(out,1)
        ny  = size(out,2) 

        nin = ceiling(1/ratio) - 1 
        frac = 1.d0 / (nin+1)
                
        out = 0.d0

        ! Loop through lo-res grid
        do inow = 1, nxl-1
            ih = 1 + ceiling((inow-1)/ratio)
          
            do jnow = 1, nyl
                jh = 1 + ceiling((jnow-1)/ratio)

                ! Save matching points horizontally
                out(jh,ih) = in(jnow,inow)
                out(jh,ih+nin+1) = in(jnow,inow+1)

                ! Fill in gaps horizontally
                do q = 1, nin
                    f = q*frac
                    out(jh,ih+q) = (1-f)*in(jnow,inow) + f*in(jnow,inow+1)
                end do

            end do
        end do

        ! Loop through hi-res grid
        do ih= 1, nx
            do jnow= 1, nyl-1

                jh  = 1 + ceiling((jnow-1)/ratio)
                jh2 = 1 + ceiling(jnow/ratio)

                ! Fill in gaps vertically
                do q = 1, nin
                    f = q*frac
                    out(jh+q,ih) = (1-f)*out(jh,ih) + f*out(jh2,ih)
                end do

            end do
        end do      

        if (m .eq. 1) then
            write(*,*) "Code to fix mask!" ! Not yet needed
        end if  

        return
    end subroutine tohires

end module interp2D