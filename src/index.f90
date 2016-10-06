
module index

    implicit none 

    interface unique 
        module procedure unique_float, unique_dble
    end interface 

    private
    public :: which, unique 

contains 

    subroutine which(x,ind,stat)
        ! Analagous to R::which function
        ! Returns indices that match condition x==.TRUE.

        implicit none 

        logical :: x(:)
        integer, allocatable :: tmp(:), ind(:)
        integer, optional :: stat  
        integer :: n, i  

        n = count(x)
        allocate(tmp(n))
        tmp = 0 

        n = 0
        do i = 1, size(x) 
            if (x(i)) then 
                n = n+1
                tmp(n) = i 
            end if
        end do 

        if (present(stat)) stat = n 

        if (allocated(ind)) deallocate(ind)
        
        if (n .eq. 0) then 
            allocate(ind(1))
            ind = -1 
        else
            allocate(ind(n))
            ind = tmp(1:n)
        end if 
        
        return 

    end subroutine which 

    subroutine unique_float(xu,x)
        ! Return only the unique values of a vector
        ! http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran

        implicit none 

        real(4) :: x(:)          ! The input
        real(4) :: res(size(x))  ! The unique values
        real(4), allocatable :: xu(:)  ! The output 
        integer :: k                   ! The number of unique elements
        integer :: i, j
        real(4), parameter :: tol = 1d-5
        logical :: found 

        k = 1
        res(1) = x(1)
        do i=2,size(x)
            found = .FALSE.
            do j=1,k
                if (abs(res(j)-x(i)) .le. tol) then 
                   ! Found a match so start looking again
                   found = .TRUE. 
                   cycle 
                end if
            end do
            ! No match found so add it to the output
            if (.not. found) then 
                k = k + 1
                res(k) = x(i)
            end if 
        end do

        ! Store output in properly sized output vector
        if(allocated(xu)) deallocate(xu)
        allocate(xu(k))
        xu = res(1:k)

        return 

    end subroutine unique_float 

    subroutine unique_dble(xu,x)
        ! Return only the unique values of a vector
        ! http://rosettacode.org/wiki/Remove_duplicate_elements#Fortran

        implicit none 

        double precision :: x(:)          ! The input
        double precision :: res(size(x))  ! The unique values
        double precision, allocatable :: xu(:)  ! The output 
        integer :: k                   ! The number of unique elements
        integer :: i, j
        double precision, parameter :: tol = 1d-5
        logical :: found 

        k = 1
        res(1) = x(1)
        do i=2,size(x)
            found = .FALSE.
            do j=1,k
                if (abs(res(j)-x(i)) .le. tol) then 
                   ! Found a match so start looking again
                   found = .TRUE. 
                   cycle 
                end if
            end do
            ! No match found so add it to the output
            if (.not. found) then 
                k = k + 1
                res(k) = x(i)
            end if 
        end do

        ! Store output in properly sized output vector
        if(allocated(xu)) deallocate(xu)
        allocate(xu(k))
        xu = res(1:k)

        return 

    end subroutine unique_dble

end module index