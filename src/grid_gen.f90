module grid_gen 
    ! This module is designed to generate various resolution
    ! versions of the same grid, with interpolation between them

    use coord_constants
    use coordinates 
    use coordinates_mapping
!     use interp2D_conservative 
    use ncio 

    implicit none 

    type multigrid_class 
        integer :: n
        type(grid_class), allocatable :: grid(:)
!         type(map_conserv_class), allocatable :: map_to(:), map_from(:)
        type(map_class), allocatable :: map_to(:), map_from(:)
    end type 

    private 
    public :: multigrid_class
    public :: multigrid_init 

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
        real(dp) :: nx_dble(size(dx)), ny_dble(size(dx))
        integer  :: nx(size(dx)), ny(size(dx))
        character(len=12)  :: dx_str 
        logical :: err 

        ! === Check consistency of desired grids ===
        ! 1. New grid resolutions result in integer multiple grid sizes
        !    of the original grid 

        nx_dble = (grid%G%nx-1)*grid%G%dx / dx + 1
        nx      = int(nx_dble)
        
        ny_dble = (grid%G%ny-1)*grid%G%dy / dx + 1
        ny      = int(ny_dble)

        err = .FALSE. 

        if (sum(nx_dble-nx) .ne. 0.d0) then 
            write(*,*) "multigrid_init:: error: desired grid resolutions &
            &are not consistent for the x-dimension."
            write(*,"(a,10f10.1)") "dx:      ", dx 
            write(*,"(a,10f10.1)") "nx:      ", nx_dble
            write(*,"(a,10f10.1)") "err(nx): ", nx_dble-nx
            err = .TRUE.
        end if 
        if (sum(ny_dble-ny) .ne. 0.d0) then 
            write(*,*) "multigrid_init:: error: desired grid resolutions &
            &are not consistent for the y-dimension."
            write(*,"(a,10f10.1)") "dx:      ", dx 
            write(*,"(a,10f10.1)") "ny:      ", ny_dble
            write(*,"(a,10f10.1)") "err(ny): ", ny_dble-ny
            err = .TRUE.
        end if 

        if (err) stop 

        ! == Sub-grid initialization == 

        ! How many sub-grids are desired? 
        mgrid%n = size(dx,1)

        ! Allocate sub-grids 
        if (allocated(mgrid%grid)) deallocate(mgrid%grid)
        allocate(mgrid%grid(mgrid%n))

        ! Allocate maps 
        if (allocated(mgrid%map_to)) deallocate(mgrid%map_to)
        allocate(mgrid%map_to(mgrid%n))

        if (allocated(mgrid%map_from)) deallocate(mgrid%map_from)
        allocate(mgrid%map_from(mgrid%n))

        ! Loop over dx values, generate new axis values and initialize
        ! subgrid objects
        do q = 1, mgrid%n 

            write(dx_str,"(i5)") int(dx(q))
            dx_str = trim(adjustl(dx_str))
            grid_name = trim(grid%name)//"_"//trim(dx_str)//"KM"

            if (allocated(x)) deallocate(x)
            if (allocated(y)) deallocate(y)
            allocate(x(nx(q)),y(ny(q)))

            do i = 1, size(x)
                x(i) = minval(grid%G%x) + (i-1)*dx(q)
            end do 
            do j = 1, size(y)
                y(j) = minval(grid%G%y) + (j-1)*dx(q)
            end do 

            call grid_init(mgrid%grid(q),grid,name=trim(grid_name),x=x,y=y)


            ! Generate conservative maps (to, from)
!             call map_conservative_init(mgrid%map_to(q),  grid,mgrid%grid(q))
!             call map_conservative_init(mgrid%map_from(q),mgrid%grid(q),grid)
            call map_init(mgrid%map_to(q),  grid,mgrid%grid(q),max_neighbors=20,lat_lim=1.d0)
            call map_init(mgrid%map_from(q),mgrid%grid(q),grid,max_neighbors=20,lat_lim=1.d0)

        end do 

        ! == Done == 

        write(*,*) "multigrid_init completed. Summary: "
        write(*,"(a5,3a8,4a10)") "Grid", "dx", "nx", "ny", "xmin", "xmax", "ymin", "ymax"
        write(*,"(i5,f8.1,2i8,4f10.2)") 1, grid%G%dx, grid%G%nx, grid%G%ny,  &
                minval(grid%G%x), maxval(grid%G%x), minval(grid%G%y), maxval(grid%G%y)
        do q = 1, mgrid%n 
            write(*,"(i5,f8.1,2i8,4f10.2)") q+1, mgrid%grid(q)%G%dx, mgrid%grid(q)%G%nx, mgrid%grid(q)%G%ny,  &
                minval(mgrid%grid(q)%G%x), maxval(mgrid%grid(q)%G%x), minval(mgrid%grid(q)%G%y), maxval(mgrid%grid(q)%G%y)
        end do 
        write(*,*) 

        return 

    end subroutine multigrid_init 

end module grid_gen 

