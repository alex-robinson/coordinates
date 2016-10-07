module grid_gen 
    ! This module is designed to generate various resolution
    ! versions of the same grid, with interpolation between 
    ! them, including staggered grids

    use coordinates 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 

    type multigrid_class 
        integer :: n_grids
        type(grid_class), allocatable :: grid(:)
    end type 

    private 
    public :: multigrid_class
    public :: multigrid_init 
    public :: grid_gen_staggered

contains 

    subroutine multigrid_init(mgrid,grid,dx)
        ! Generate multiple grids of different resolutions
        ! (descending from high to low resolution)

        implicit none 

        type(multigrid_class), intent(OUT) :: mgrid   ! Multigrid 
        type(grid_class),      intent(IN)  :: grid    ! Base grid (highest resolution)
        real(dp), intent(IN) :: dx(:)   ! Sub-grid resolutions 

        ! Local variables 
        integer :: q, i, j  
        character(len=256) :: grid_name 
        real(dp), allocatable :: x(:), y(:) 
        real(dp) :: nx_dble(size(dx)), ny_dble(size(dy))
        integer  :: nx(size(dx)), ny(size(dy))
        character(len=12)  :: dx_str 

        ! How many sub-grids are desired? 
        mgrid%n_grids = size(dx,1)

        ! === Check consistency of desired grids ===
        ! 1. New grid resolutions result in integer multiple grid sizes
        !    of the original grid 

        nx_dble = (grid%G%nx-1)*grid%G%dx / dx
        nx      = int(nx_dble)
        if (sum(nx_dble-nx) .ne. 0.d0) then 

            write(*,*) "multigrid_init:: error: desired grid resolutions &
            &are not consistent."
            write(*,"(a,10f10.1)") "dx:      ", dx 
            write(*,"(a,10f10.1)") "nx:      ", nx_dble
            write(*,"(a,10f10.1)") "err(nx): ", nx_dble-nx
            stop 

        end if 


        ! Allocate sub-grids 
        if (allocated(mgrid%grid)) deallocate(mgrid%grid)
        allocate(mgrid%grid(mgrid%n_grids))


        ! Loop over dx values, generate new axis values and initialize
        ! subgrid objects
        do q = 1, mgrid%n_grids 

            write(dx_str,"(i5)") int(dx(q))
            dx_str = trim(adjustl(dx_str))
            grid_name = trim(grid%name)//"_"//trim(dx_str)//"KM"

            if (allocated(x)) deallocate(x)
            if (allocated(y)) deallocate(y)
            allocate(x(nx(q)),y(10))

            do i = 1, size(x)
                x(i) = minval(grid%G%x) + (i-1)*dx(q)
            end do 
            do j = 1, size(y)
                y(j) = minval(grid%G%y) + (j-1)*dx(q)
            end do 

            call grid_init(mgrid%grid(q),grid,name=trim(grid_name),x=x,y=y)

        end do 

        return 

    end subroutine multigrid_init 


    subroutine grid_gen_staggered()

        implicit none 


        return 

    end subroutine grid_gen_staggered


end module grid_gen 

