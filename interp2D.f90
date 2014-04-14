
module interp2D
    
    implicit none 

    !! real(dp) definition and some internal constants
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    private
    public :: interp_bilinear_grid

contains


    function interp_bilinear_grid(x,y,z,xout,yout) result(zout)
        ! Perform bilinear interpolation from a regular grid
        ! to a regular grid. This means, linear interpolation
        ! row-by-row, then column by column (not as general as
        ! pure bilinear interpolation from random nearest neighbors)

        implicit none 

        real(dp), dimension(:,:) :: x, y, z 
        real(dp), dimension(:,:) :: xout, yout 
        real(dp), dimension(size(xout,1),size(xout,2)) :: zout 


        return 

    end function interp_bilinear_grid




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