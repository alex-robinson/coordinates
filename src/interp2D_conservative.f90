module interp2D_conservative
    
    use coordinates 
    use coordinates_mapping 

    use planet 
    use oblimap_projection_module
    use polygons 
    use geodesic 
    use index 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)
    real(dp), parameter :: MISSING_VALUE_DEFAULT = -9999.0_dp 
    real(dp), parameter :: mv = MISSING_VALUE_DEFAULT

!     interface map_calc_conserv1_weights
!         module procedure map_calc_conserv1_weights_points_points
!         module procedure map_calc_conserv1_weights_grid_grid
!     end interface 

    private 
!     public :: map_calc_conserv1_weights
    public :: map_field_conservative

contains 

    subroutine map_field_conservative(grid1,grid2,varname,var1,var2,missing_value)

        implicit none 

        type(grid_class), intent(IN)  :: grid1  ! Original grid information
        type(grid_class), intent(IN)  :: grid2  ! New grid 
        character(len=*), intent(IN)  :: varname        ! Name of the variable being mapped
        double precision, intent(IN)  :: var1(:,:)      ! Input variable
        double precision, intent(OUT) :: var2(:,:)      ! Output variable
        double precision, intent(IN), optional :: missing_value  ! Points not included in mapping

        double precision :: missing_val 

        double precision :: x1, y1, x2, y2 
        integer :: i, j 
        integer, allocatable :: ii(:), jj(:) 

        real(dp) :: area(size(grid1%x,1),size(grid1%x,2))

        ! Helpful fields available from grid_class objects grid and grid1:
        ! grid%x    : 2D array of projected x-values [km]
        ! grid%y    : 2D array of projected y-values [km]
        ! grid%lon  : 2D array of lon-values [degrees]
        ! grid%lat  : 2D array of lat-values [degrees]
        ! grid%area : 2D array of grid-cell areas [m]
        ! grid%G%x  : vector of x-values that defines the x-axis [km]
        ! grid%G%y  : vector of y-values that defines the y-axis [km]
        ! grid%G%nx : length of x-axis 
        ! grid%G%ny : length of y-axis 

        ! Check if the grids are compatible for this mapping routine
        if (.not. same_projection(grid1%proj,grid2%proj)) then 
            write(*,*) "map_field_conservative:: error:  &
                       &Currently this subroutine only interpolates between grids &
                       &on the same projection. Try again."
            stop 
        end if 

        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! Initially set all values to missing
        var2 = missing_val

        ! A loop to get started 
        do j = 1, grid2%G%ny
            write(*,*) "j / ny: ", j, grid2%G%ny 

            do i = 1, grid2%G%nx 

                ! == TO DO == 

                call which(grid1%G%x .ge. grid2%x(i,j)-grid2%G%dx .and. grid1%G%x .le. grid2%x(i,j)+grid2%G%dx,ii)
                call which(grid1%G%y .ge. grid2%y(i,j)-grid2%G%dy .and. grid1%G%y .le. grid2%y(i,j)+grid2%G%dy,jj)
                
                area(ii,jj) = interpconserv1_weights(x=grid1%x(ii,jj),y=grid1%y(ii,jj),dx=grid1%G%dx,dy=grid1%G%dy, &
                                              xout=grid2%x(i,j),yout=grid2%y(i,j), &
                                              dxout=grid2%G%dx,dyout=grid2%G%dy)

                if (sum(area(ii,jj)) .gt. 0.d0) then 

                    var2(i,j) = sum(var1(ii,jj)*area(ii,jj),mask=area(ii,jj).gt.0.d0) &
                                  / sum(area(ii,jj),mask=area(ii,jj).gt.0.d0)

                else 
                    var2(i,j) = missing_val 

                end if 

            end do 
        end do 

        return 

    end subroutine map_field_conservative 

    function interpconserv1_weights(x,y,dx,dy,xout,yout,dxout,dyout,latlon,missing_value) result(area)
        ! Calculate 1st order conservative interpolation for a 
        ! point given its nearest neighbors  
        real(dp), intent(IN) :: x(:,:), y(:,:), dx, dy 
        real(dp), intent(IN) :: xout, yout, dxout, dyout 
        logical,  intent(IN), optional :: latlon
        real(dp), intent(IN), optional :: missing_value  
        real(dp) :: area(size(x,1),size(x,2))

        ! Local variables
        logical  :: is_latlon 
        type(polygon)       :: pol  
        integer, parameter :: nx = 15, ny = 15, npts = nx*ny
        real(dp) :: x1, y1
        integer  :: npts_in 
        integer :: i, j, now  
        real(dp) :: missing_val
        real(dp) :: x_vec(size(x,1)*size(x,2)), y_vec(size(x,1)*size(x,2)) 
        real(dp) :: area_vec(size(x,1)*size(x,2))
        real(dp) :: area_target 

        is_latlon = .FALSE. 
        if (present(latlon)) is_latlon = latlon 

        missing_val = mv 
        if (present(missing_value)) missing_val = missing_value

        ! Generate polygon representing boundaries of target point 
        pol = create_polygon(real([xout-dxout/2.d0,xout-dxout/2.d0,xout+dxout/2.d0,xout+dxout/2.d0]), &
                             real([yout-dyout/2.d0,yout+dyout/2.d0,yout+dyout/2.d0,yout-dyout/2.d0]))

        x_vec = reshape(x,[size(x_vec)])
        y_vec = reshape(y,[size(y_vec)])
        
        ! Loop over source points and get the area of each source
        ! polygon that is inside of the target polygon 
        ! - Save the area (absolute area, not fraction)
        area_target = dxout*dyout 
        area_vec    = 0.d0 

        do now = 1, size(x_vec) 
         
            npts_in = 0
            do i = 1, nx
            do j = 1, ny 
                x1 = (x_vec(now)-dx/2.d0) + (dx)*dble(i-1)/dble(nx) + (dx)*0.5d0/dble(nx-1) 
                y1 = (y_vec(now)-dy/2.d0) + (dy)*dble(j-1)/dble(ny) + (dy)*0.5d0/dble(ny-1) 
                if (point_in_polygon(real(x1),real(y1),pol)) npts_in = npts_in+1

            end do
            end do 
            
            area_vec(now) = dble(npts_in)/dble(npts) * dx*dy 

            ! If the source points area adds up to the target cell area, exit loop
            if (sum(area_vec) .ge. area_target) exit 
        end do 

        ! Return 2D area 
        area = reshape(area_vec,[size(x,1),size(x,2)])

        return 

    end function interpconserv1_weights

!     subroutine map_calc_conserv1_weights_grid_grid(map,grid1)

!         implicit none 

!         type(map_class),  intent(INOUT) :: map
!         type(grid_class), intent(IN)    :: grid1   

!         type(points_class) :: pts1 
!         integer :: i, k, ntot  

!         ! Convert the grid to points for calculating weights more easily
!         call grid_to_points(grid1,pts1)

!         ! Update the map weights with area inside polygon calculations
!         do i = 1, map%npts  
 
!             do k = 1, map%nmax 
!                 if (map%i(i,k) .le. 0) exit 
!             end do 
!             ntot = k-1 

!             if (ntot .gt. 0) then
!                 ! Update weights if their are neighbors available 

!                 map%weight(i,1:ntot) = interpconserv1_weights(x=pts1%lon(map%i(i,1:ntot)),y=pts1%lat(map%i(i,1:ntot)), &
!                                                       dx=grid1%G%dx,dy=grid1%G%dy,xout=map%x(i),yout=map%y(i), &
!                                                       dxout=map%G%dx,dyout=map%G%dy)
!                 write(*,*) i, minval(map%weight(i,1:ntot)), maxval(map%weight(i,1:ntot))
!             end if 

!         end do 

!         write(*,*) "** Updated map weights to conservative (1st) order scheme."

!         return 

!     end subroutine map_calc_conserv1_weights_grid_grid

!     subroutine map_calc_conserv1_weights_points_points(map,pts1,dx1,dy1,dx2,dy2)

!         implicit none 

!         type(map_class),    intent(INOUT) :: map
!         type(points_class), intent(IN)    :: pts1 
!         real(dp), intent(IN) :: dx1, dy1, dx2, dy2  

!         integer :: i

!         ! Update the map weights with area inside polygon calculations
!         do i = 1, map%npts  

!             map%weight(i,:) = interpconserv1_weights(x=pts1%x(map%i(i,:)),y=pts1%y(map%i(i,:)), &
!                                                       dx=dx1,dy=dy1,xout=map%x(i),yout=map%y(i), &
!                                                       dxout=dx2,dyout=dy2)

!         end do 

!         write(*,*) "** Updated map weights to conservative (1st) order scheme."
        
!         return 

!     end subroutine map_calc_conserv1_weights_points_points

end module interp2D_conservative 
