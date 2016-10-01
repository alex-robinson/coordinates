module grid_gen 
    ! This module is designed to generate various resolution
    ! versions of the same grid, with interpolation between 
    ! them, including staggered grids

    use coordinates 

    implicit none 


    private 
    public :: grid_gen_staggered

contains 

    subroutine grid_gen_staggered()

        implicit none 


        return 

    end subroutine grid_gen_staggered


end module grid_gen 

