

program test_proj 

    use coord 
    use ncio 
    
    implicit none 

    type(points_class) :: pts 

    integer :: i 

    call points_init(pts,name="Bamber-corners",mtype="polar_stereographic",units="kilometers", &
                     lambda=-39.d0,phi=90.d0,alpha=19.d0,lon180=.TRUE., &
                     x=[-800.d0,-800.d0,700.d0,700.d0],y=[-3400.d0,-600.d0,-600.d0,-3400.d0])


    do i = 1, size(pts%x,1)
        write(*,"(2f10.3,a,2f10.3)") pts%x(i), pts%y(i), " => ", pts%lat(i), pts%lon(i)
    end do 

end program test_proj