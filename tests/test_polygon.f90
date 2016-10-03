
! Ray-casting algorithm and point_in_polygon test from Rosetta code:
! http://rosettacode.org/wiki/Ray-casting_algorithm#Fortran

program test_points
    
    use polygons 

    implicit none 

    integer :: npts 
    real(4), allocatable :: xx(:), yy(:)
    real(4) :: x, y 
    logical :: inside 

!     npts = 14 
!     allocate(xx(npts),yy(npts))

!     xx = [0.0,10.0,10.0,0.0,2.5,7.5,7.5,2.5,0.0,10.0,3.0,7.0,7.0,3.0]
!     yy = [0.0,0.0,10.0,10.0,2.5,2.5,7.5,7.5,5.0,5.0,0.0,0.0,10.0,10.0]
    
    npts = 4
    allocate(xx(npts),yy(npts))

    xx = [0.0,1.0,1.0,0.0]
    yy = [0.0,0.0,1.0,1.0]

    x = 0.1 
    y = 0.1 

    write(*,*) "Points inside polygon?"
    inside = point_in_polygon(x, y, xx, yy)
    write(*,*) x, y, inside

end program 