
module interp2D
    
    implicit none 

    !! real(dp) definition and some internal constants
    integer,  parameter :: dp  = kind(1.0d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: ERR_DIST = 1E8_dp 
    integer,  parameter :: ERR_IND  = -1 
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    interface interp_bilinear 
        module procedure interp_bilinear_dble 
    end interface

    private
    public :: interp_bilinear
    public :: fill_weighted

contains


    function interp_bilinear_dble(x,y,z,xout,yout,missing_value,mask) result(zout)
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
                    
                    ! Only interpolate points with all neighbors available
                    if (count([z(i-1,j),z(i,j),z(i,j-1),z(i-1,j-1)] .eq. missing_val) .eq. 0) then
                        alpha1 = (xout(i1) - x(i-1)) / (x(i)-x(i-1))
                        p0 = z(i-1,j-1) + alpha1*(z(i,j-1)-z(i-1,j-1))
                        p1 = z(i-1,j)   + alpha1*(z(i,j)-z(i-1,j))
                        
                        alpha2 = (yout(j1) - y(j-1)) / (y(j)-y(j-1))
                        zout(i1,j1) = p0 + alpha2*(p1-p0)
                    end if 

                end if 

            end if 

        end do 
        end do
        

        return 

    end function interp_bilinear_dble

!     subroutine fill_bilinear_dble(z,missing_value)

!         implicit none 

!         real(dp), intent(INOUT) :: z(:,:)
!         real(dp), intent(IN) :: missing_value
!         integer :: i, j, q, nx, ny 
!         integer :: x_idx(2), y_idx(2)
!         real(dp) :: neighb(4), weights(4) 

!         nx = size(z,1)
!         ny = size(z,2)

!         do i = 1, nx 
!             do j = 1, ny 

!                 if (z(i,j) .eq. missing_value) then 
!                     x_idx = find_neighbors_1D(z(:,j),i,missing_value)
!                     y_idx = find_neighbors_1D(z(i,:),j,missing_value)

!                     ! At least one of four neighbors found
!                     if (count(x_idx .lt. 0) .eq. 0 .and. count(y_idx .lt. 0) .lt. 4) then 

!                     weights = 0.d0 
!                     if (x_idx(1) .lt. 0 .and. x_idx(2)) then 
!                         neighb(1:2)  = 0.d0 
!                         weights(1:2) = 0.d0 
!                     else 
!                         neighb(1)  = z(x_idx(1),j)
!                         weights(1) = 
!                     ! ## Now handle different interpolation cases depending on available neighbors

!                     8
!                     end if 


!                 end if 

!             end do 
!         end do 

!         return 

!     end subroutine fill_bilinear_dble 


    function find_neighbors_1D(x,i,missing) result(ineighb)

        implicit none 

        double precision, intent(IN) :: x(:), missing 
        integer, intent(IN) :: i
        integer :: ineighb(2) 
        integer :: j, step  


        ! Find negative neighbors 
        ineighb(1) = -1 
        if (i .gt. 1) then 
            do j = i-1, 1, -1
                if (.not. x(j) .eq. missing) then 
                    ineighb(1) = j 
                    exit
                end if  
            end do 
        end if 

        ! Find positive neighbors 
        ineighb(2) = -1 
        if (i .lt. size(x)) then 
            do j = i+1, size(x)
                if (.not. x(j) .eq. missing) then 
                    ineighb(2) = j
                    exit
                end if 
            end do 
        end if 



        return 

    end function find_neighbors_1D

    subroutine fill_weighted(var,missing_value,fill_value,nr)
        implicit none 
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j
        integer, parameter :: qmax = 40 ! Iterations 
!         integer, parameter :: nr   = 1  ! Neighbor radius
        integer :: nr 
        double precision, dimension (2*nr+1,2*nr+1) :: neighb, weight, weight0 
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        ! Define the default neighbor weighting 
        do i = -nr, nr 
        do j = -nr, nr 
            weight0(i,j) = dsqrt(dble(i)*dble(i) + dble(j)*dble(j))
        end do 
        end do 


        do q = 1, qmax 

            filled = missing_value 

            do i = 1+nr, nx-nr
                do j = 1+nr, ny-nr
                    
                    if (var(i,j) .eq. missing_value) then 

                        neighb = var(i-nr:i+nr,j-nr:j+nr)

                        weight = weight0
                        where (neighb .eq. missing_value) weight = 0.d0
                        wtot = sum(weight)

                        if (wtot/sum(weight0) .gt. 0.5d0) filled(i,j) = sum(neighb*weight)/wtot

                    end if 
                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            !write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var(1+nr+1:nx-nr-1,1+nr+1:ny-nr-1) .eq. missing_value) .eq. 0 ) exit 
            write(*,*) "Still missing... ", count(var(1+nr:nx-nr,1+nr:ny:nr) .eq. missing_value), &
                        " of ", nx*ny
        end do 

        ! Fill in boundaries too 
        do i = 1+nr, 1, -1
            where(var(i,:) .eq. missing_value) var(i,:) = var(i+1,:)
        end do 
        do i = nx-nr, nx
            where(var(i,:) .eq. missing_value) var(i,:) = var(i-1,:)
        end do 
        do j = 1+nr, 1, -1
            where(var(:,j) .eq. missing_value) var(:,j) = var(:,j+1)
        end do 
        do j = ny-nr, ny
            where(var(:,j) .eq. missing_value) var(:,j) = var(:,j-1)
        end do 


!         if (count(var .eq. missing_value) .gt. 0) then 
!             where(var .eq. missing_value) var = minval(var,mask=var .ne. missing_value)
!         end if 

        write(*,*) "Fill iterations: ",q 

        return
    end subroutine fill_weighted

    subroutine fill_weighted1(var,missing_value,fill_value,nr)

        implicit none 
        
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j
        integer, parameter :: qmax = 40 ! Iterations 
!         integer, parameter :: nr   = 1  ! Neighbor radius
        integer :: nr 
        double precision, dimension (2*nr+1,2*nr+1) :: neighb, weight, weight0 
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled

        logical, dimension(:,:), allocatable :: is_missing
        integer, dimension(:,:), allocatable :: n_missing
        integer :: di_lo, di_hi, dj_lo, dj_hi, ntot 

!         integer, dimension(:), allocatable :: nx_miss, ny_miss
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))
        allocate(is_missing(nx,ny),n_missing(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        ! Define the default neighbor weighting 
        do i = -nr, nr 
        do j = -nr, nr 
            weight0(i,j) = dsqrt(dble(i)*dble(i) + dble(j)*dble(j))
        end do 
        end do 

        do q = 1, 1 !qmax 

            filled = missing_value 

            is_missing = .FALSE. 
            where(var .eq. missing_value) is_missing = .TRUE. 

            ! Find out where missing neighbor values are aggregated
            n_missing = 0 
            do i = 1, nx
                do j = 1, ny
                    if (var(i,j) .eq. missing_value) then 
                        di_lo = nr - max(nr-( i-1),0)
                        di_hi = nr - max(nr-(nx-i),0)
                        dj_lo = nr - max(nr-( j-1),0)
                        dj_hi = nr - max(nr-(ny-j),0)

                        n_missing(i,j) = count(var(i-di_lo:i+di_hi,j-dj_lo:j+dj_hi) .eq. missing_value)

                        ! Neighborhood size
                        ntot = (di_hi+di_lo+1)*(dj_hi+dj_lo+1)

!                         if (i .eq. nx .and. j .eq. ny) then 
!                             write(*,*) i, j, n_missing(i,j), ntot, dble(n_missing(i,j))/dble(ntot)
!                             write(*,*) di_lo, di_hi, dj_lo, dj_hi 
!                             write(*,*) var(i-di_lo:i+di_hi,j-dj_lo:j+dj_hi)
!                             stop 
!                         end if 

                        ! Prioritizes interpolation of points with more of neighbors
                        if (dble(n_missing(i,j))/dble(ntot) .le. 1.0d0) then 

                            neighb = missing_value 
                            neighb(1:di_lo+1+di_hi,1:dj_lo+1+dj_hi) = var(i-di_lo:i+di_hi,j-dj_lo:j+dj_hi)
                            weight = 1.d0
                            where (neighb .eq. missing_value) weight = 0.d0

                            filled(i,j) = sum(neighb*weight)/sum(weight) 
                        end if 
                    end if 
                end do 
            end do 

            ! Now fill in original matrix with this pass's fill values
            where(filled .ne. missing_value) var = filled 

            !write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 
        end do

        write(*,*) "Fill iterations: ",q 

        return

    end subroutine fill_weighted1

    subroutine fill_weighted0(var,missing_value,fill_value,nr)
        implicit none 
        double precision, dimension(:,:) :: var 
        double precision :: missing_value 
        double precision, optional :: fill_value

        integer :: q, nx, ny, i, j
        integer, parameter :: qmax = 40 ! Iterations 
!         integer, parameter :: nr   = 1  ! Neighbor radius
        integer :: nr 
        double precision, dimension (2*nr+1,2*nr+1) :: neighb, weight, weight0 
        double precision :: wtot, mval 
        double precision, dimension(:,:), allocatable :: filled
        nx = size(var,1)
        ny = size(var,2) 

        allocate(filled(nx,ny))

        if (present(fill_value)) then
            where(var .eq. missing_value) var = fill_value 
        end if 

        ! Define the default neighbor weighting 
        do i = -nr, nr 
        do j = -nr, nr 
            weight0(i,j) = dsqrt(dble(i)*dble(i) + dble(j)*dble(j))
        end do 
        end do 


        do q = 1, qmax 

            filled = missing_value 

            do i = 1+nr, nx-nr
                do j = 1+nr, ny-nr
                    neighb = var(i-nr:i+nr,j-nr:j+nr)

                    weight = 0.d0 
                    where (neighb .ne. missing_value) weight = 1.d0
!                     weight = weight0
!                     where (neighb .eq. missing_value) weight = 0.d0
                    wtot = sum(weight)

                    if (wtot .gt. 0.d0) then 
                        mval = sum(neighb*weight)/wtot
                        where (neighb .eq. missing_value) neighb = mval 
                    end if 

                    filled(i-nr:i+nr,j-nr:j+nr) = neighb 

                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            !write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 

            ! Now reverse directions for better results
            filled = missing_value 

            do i = nx-nr, 1+nr, -1
                do j = ny-nr, 1+nr, -1
                    neighb = var(i-nr:i+nr,j-nr:j+nr)

                    weight = 0.d0 
                    where (neighb .ne. missing_value) weight = 1.d0
                    wtot = sum(weight)

                    if (wtot .gt. 0.d0) then 
                        mval = sum(neighb*weight)/wtot
                        where (neighb .eq. missing_value) neighb = mval 
                    end if 

                    filled(i-nr:i+nr,j-nr:j+nr) = neighb 

                end do 
            end do 

            where(filled .ne. missing_value) var = filled 

            !write(*,*) q," : Missing values: ", count(var .eq. missing_value)
            if ( count(var .eq. missing_value) .eq. 0 ) exit 
        end do 
        if (count(var .eq. missing_value) .gt. 0) then 
            where(var .eq. missing_value) var = minval(var,mask=var .ne. missing_value)
        end if 

        write(*,*) "Fill iterations: ",q 

        return
    end subroutine fill_weighted0

end module interp2D