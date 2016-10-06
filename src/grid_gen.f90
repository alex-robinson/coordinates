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
    public :: multigrid_init 
    public :: grid_gen_staggered

contains 

    subroutine multigrid_init(mgrid,grid,dx)

        implicit none 

        type(multigrid_class ), intent(OUT) :: mgrid   ! Multigrid 
        type(grid_class ),      intent(IN)  :: grid    ! Base grid (highest resolution)
        real(dp), intent(IN) :: dx(:)   ! Sub-grid resolutions 

        return 

    end subroutine multigrid_init 


    subroutine grid_gen_staggered()

        implicit none 


        return 

    end subroutine grid_gen_staggered


end module grid_gen 

