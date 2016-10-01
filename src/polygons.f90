
! Ray-casting algorithm and point_in_polygon test from Rosetta code:
! http://rosettacode.org/wiki/Ray-casting_algorithm#Fortran

module polygons 

    implicit none

    real, parameter :: eps = 0.00001

    type point
        real :: x, y
    end type point

    interface operator (-)
        module procedure pt_sub
    end interface

    interface len
        module procedure pt_len
    end interface
    
    type polygon
        type(point), dimension(:), allocatable :: points
        integer,     dimension(:), allocatable :: vertices
    end type polygon

    interface create_polygon
        module procedure create_polygon_vts
        module procedure create_polygon_novts
        module procedure create_polygon_vals
    end interface 

    interface point_in_polygon 
        module procedure point_is_inside_poly
        module procedure point_is_inside_points
        module procedure points_1D_is_inside_points
        module procedure points_2D_is_inside_points   
    end interface 

    private 
    public :: point, polygon  
    public :: create_polygon
    public :: point_in_polygon

contains
    
    ! ## POINT functions ##
    function pt_sub(a, b) result(c)
        type(point), intent(in) :: a, b
        type(point) :: c

        c = point(a%x - b%x, a%y - b%y)
    end function pt_sub
 
    function pt_len(a) result(l)
        type(point), intent(in) :: a
        real :: l

        l = sqrt((a%x)**2 + (a%y)**2)
    end function pt_len
 

    ! ## POLYGON functions ##
    function create_polygon_vts(pts, vt)
        type(polygon) :: create_polygon_vts
        type(point), dimension(:), intent(in) :: pts
        integer, dimension(:), intent(in) :: vt

        integer :: np, nv

        np = size(pts,1)
        nv = size(vt,1)

        allocate(create_polygon_vts%points(np), create_polygon_vts%vertices(nv))
        create_polygon_vts%points   = pts
        create_polygon_vts%vertices = vt

    end function create_polygon_vts
    
    function create_polygon_novts(pts)
        type(polygon) :: create_polygon_novts
        type(point), dimension(:), intent(in) :: pts
        integer, dimension(:), allocatable :: vt

        integer :: np, nv
        integer :: i, k 

        np = size(pts,1)
        nv = np*2

        ! Define vertices (loop around point indices)
        allocate(vt(nv)) 
        vt(1)  = 1
        vt(nv) = 1 
        
        do i = 1, np
            if (i > 1) then
                k = i*2-1 
                vt(k-1:k) = i 
            end if 
        end do 

        allocate(create_polygon_novts%points(np), create_polygon_novts%vertices(nv))
        create_polygon_novts%points   = pts
        create_polygon_novts%vertices = vt

    end function create_polygon_novts
    
    function create_polygon_vals(xx,yy)
        real(4), intent(IN) :: xx(:), yy(:)
        type(polygon) :: create_polygon_vals
        type(point), dimension(:), allocatable :: pts
        integer, dimension(:), allocatable :: vt

        integer :: np, nv
        integer :: i, k 

        np = size(xx)
        nv = np*2

        ! Define pts 
        allocate(pts(np))
        do i = 1, np 
            pts(i) = point(xx(i),yy(i))
        end do 

        ! Define vertices (loop around point indices)
        allocate(vt(nv)) 
        vt(1)  = 1
        vt(nv) = 1 
        
        do i = 1, np
            if (i > 1) then
                k = i*2-1 
                vt(k-1:k) = i 
            end if 
        end do 

        allocate(create_polygon_vals%points(np), create_polygon_vals%vertices(nv))
        create_polygon_vals%points   = pts
        create_polygon_vals%vertices = vt

    end function create_polygon_vals
 
    subroutine free_polygon(pol)
        type(polygon), intent(inout) :: pol

        deallocate(pol%points, pol%vertices)

    end subroutine free_polygon
 

 
    ! ## RAY CASTING functions ##

    function ray_intersects_seg(p0, a0, b0) result(intersect)
        type(point), intent(in) :: p0, a0, b0
        logical :: intersect

        type(point) :: a, b, p
        real :: m_red, m_blue

        p = p0
        ! let variable "a" be the point with smallest y coordinate
        if ( a0%y > b0%y ) then
            b = a0
            a = b0
        else
            a = a0
            b = b0
        end if

        if ( (p%y == a%y) .or. (p%y == b%y) ) p%y = p%y + eps
 
        intersect = .false.
 
        if ( (p%y > b%y) .or. (p%y < a%y) ) return
        if ( p%x > max(a%x, b%x) ) return

        if ( p%x < min(a%x, b%x) ) then
            intersect = .true.
        else
            if ( abs(a%x - b%x) > tiny(a%x) ) then
                m_red = (b%y - a%y) / (b%x - a%x)
            else
                m_red = huge(m_red)
            end if
            if ( abs(a%x - p%x) > tiny(a%x) ) then
                m_blue = (p%y - a%y) / (p%x - a%x)
            else
                m_blue = huge(m_blue)
            end if
            if ( m_blue >= m_red ) then
                intersect = .true.
            else
                intersect = .false.
            end if
        end if
 
    end function ray_intersects_seg
    
    ! ## Functions to test individual points inside of polygons ##
    function point_is_inside_poly(p, pol) result(inside)
        logical :: inside
        type(point),   intent(in) :: p
        type(polygon), intent(in) :: pol

        integer :: i, cnt, pa, pb

        cnt = 0
        do i = lbound(pol%vertices,1), ubound(pol%vertices,1), 2
            pa = pol%vertices(i)
            pb = pol%vertices(i+1)
            if ( ray_intersects_seg(p, pol%points(pa), pol%points(pb)) ) cnt = cnt + 1
        end do

        inside = .true.
        if ( mod(cnt, 2) == 0 ) inside = .false.

    end function point_is_inside_poly
 

    function point_is_inside_points(x, y, xx, yy) result(inside)
        real(4) :: x, y
        real(4) :: xx(:), yy(:) 
        logical :: inside

        type(point)   :: p
        type(polygon) :: pol

        integer :: i, cnt, pa, pb

        ! Create point from values
        p = point(x,y)

        ! Create polygon from border values
        pol = create_polygon(xx,yy)

        ! Calculate number of ray intersections
        cnt = 0
        do i = lbound(pol%vertices,1), ubound(pol%vertices,1), 2
            pa = pol%vertices(i)
            pb = pol%vertices(i+1)
            if ( ray_intersects_seg(p, pol%points(pa), pol%points(pb)) ) cnt = cnt + 1
        end do

        ! If count is even, point is outside polygon
        inside = .true.
        if ( mod(cnt, 2) == 0 ) inside = .false.

    end function point_is_inside_points
 

    function points_1D_is_inside_points(x, y, xx, yy) result(inside)
        real(4) :: x(:), y(:)
        real(4) :: xx(:), yy(:) 
        logical :: inside(size(x))

        integer :: i 

        do i = 1, size(x) 
            inside(i) = point_is_inside_points(x(i),y(i),xx,yy)
        end do 

        return 

    end function points_1D_is_inside_points
 
    function points_2D_is_inside_points(x, y, xx, yy) result(inside)
        real(4) :: x(:,:), y(:,:)
        real(4) :: xx(:), yy(:) 
        logical :: inside(size(x,1),size(x,2))

        integer :: i, j 

        do i = 1, size(x,1)
        do j = 1, size(x,2) 
            inside(i,j) = point_is_inside_points(x(i,j),y(i,j),xx,yy)
        end do 
        end do 
        
        return 

    end function points_2D_is_inside_points
 





end module polygons
