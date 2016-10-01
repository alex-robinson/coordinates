

program test_nat 

!     implicit none 

!

! Code converted using TO_F90 by Alan Miller
! Date: 2015-03-17  Time: 14:20:49

!    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977
!    PROGRAM TESTC 4.3 DATED 26 MARCH 1981

!    COPYRIGHT (C) 1981 UNIVERSITY OF BATH


!    This program tests the TILE 4 routines
!        SCATRC
!        TILE4M
!        GRLD
!        RECT1G
!    and the routines that they call.

!    It constructs a tessellation of a random point pattern
!    on a rectangle, establishes the values of some test
!    functions at the data points and then reconstructs the
!    functions over a smaller rectangle using the C1 natural
!    neighbour interpolant.  The values of the original
!    functions and their reconstructed values are written
!    to a file without carriage control characters for
!    later plotting.  This can be supressed by setting the variable
!    NSTORE zero or negative.  Statistics on the accuracy of the
!    reconstructions are written to another file with carriage
!    control characters.

!2      INTEGER*2 L
EXTERNAL uni4m
DIMENSION cn(3,4),jaddr(4),pt(2,400),naddr(400),l(6000),  &
    val(400),sbarea(50),ptoff(2,50),valoff(50),grad(2,400),  &
    delsba(2,50),z(25,25),zx(25,25),zy(25,25)

DATA jcns,npts,ltop,knbmax,mm,nn/4,400,6000,50,25,25/
DATA xmin1,xmax1,ymin1,ymax1/-0.5,1.5,-0.5,1.5/
DATA xmin2,xmax2,ymin2,ymax2/0.0,1.0,0.0,1.0/
DATA mseed/97129/
DATA nfmax/8/
DATA epscn,epspt/0.0,0.0/
DATA nstart/0/
DATA nwrite,nstore/12,13/

!    Record NFMAX, MM, NN and the smaller rectangle
!    on NSTORE to be used for latter plotting if required.

IF(nstore <= 0) GO TO 8
WRITE(nstore,8000)nfmax,mm,nn
WRITE(nstore,8001)xmin2,xmax2,ymin2,ymax2

!    Calculate the values needed to scan the grid on the
!    smaller rectangle.

8   deltax = (xmax2-xmin2)/FLOAT(mm)
deltay = (ymax2-ymin2)/FLOAT(nn)
xnum = FLOAT(mm*nn)
xstart = xmin2-0.5*deltax
ystart = ymin2-0.5*deltay

!    Generate NPTS points on the larger rectangle and tessellate
!    them.

msd = mseed
CALL scatrc(cn,jcns,pt,npts,iflag,uni4m,msd,xmin1,xmax1, ymin1,ymax1)
IF(iflag == 0) GO TO 1
WRITE(nwrite,9000)iflag
STOP 612
1   CALL tile4m(cn,jcns,jaddr,pt,npts,naddr,nfree,nstart,nptsin,  &
    l,ltop,lptr,lbase,epscn,epspt,iflag,igarb)
IF(iflag == 0) GO TO 2
WRITE(nwrite,9000)iflag
STOP 613

!    Loop through the NFMAX functions.

2   DO   nfunc = 1,nfmax
  
!    Set the function values at the data sites in VAL.
  
  DO   npt = 1,npts
    val(npt) = func(pt(1,npt),pt(2,npt),nfunc)
  END DO
  
!    Call GRLD to estimate the gradients at the data sites.
  
  CALL grld(cn,jcns,pt,npts,naddr,l,ltop,iflag,val,sbarea,  &
      ptoff,valoff,knbmax,grad)
  IF(iflag == 0) GO TO 4
  WRITE(nwrite,9000)iflag
  STOP 614
  
!    Call RECT1G to compute the C1 interpolant for the function
!    at the nodes of the MM*NN grid on the smaller rectangle.
  4   CALL rect1g(cn,jcns,jaddr,pt,npts,naddr,nfree,nstart,nptsin,  &
      l,ltop,lptr,lbase,epscn,epspt,iflag,igarb,sbarea,delsba,  &
      ptoff,knbmax,val,grad,z,zx,zy,mm,mm,nn,xmin2,xmax2,ymin2, ymax2,zlo,zhi)
  IF(iflag == 0) GO TO 5
  WRITE(nwrite,9000)iflag
  STOP 615
  
!    Record the original function and its reconstruction on NSTORE
!    if required and compute error statistics.
  
  q = 0 

  5   rmse = 0.0
  abse = 0.0
  amxe = 0.0
  xpt = xstart
  DO   m = 1,mm
    ypt = ystart
    xpt = xpt+deltax
    DO   n = 1,nn
      ypt = ypt+deltay
      f = func(xpt,ypt,nfunc)
      fzdif = ABS(z(m,n)-f)
      rmse = rmse+fzdif*fzdif
      abse = abse+fzdif
      amxe = AMAX1(amxe,fzdif)
      IF(nstore <= 0) CYCLE
      WRITE(nstore,8001) xpt, ypt, f,z(m,n),zx(m,n),zy(m,n)

    END DO
  END DO
  
!    Write out the error statistics.
  
  rmse = SQRT(rmse/xnum)
  abse = abse/xnum
  WRITE(nwrite,9001)nfunc,rmse,abse,amxe
END DO

!    Run complete.

WRITE(nwrite,9002)
STOP


8000   FORMAT(3(i13,2X))
8001   FORMAT(6(e13.6,2X))
9000   FORMAT(///1H ,19(1H*)/1H ,16H error: iflag = ,i3/1H ,19(1H*)///)

9001   FORMAT(///1H ,10HFUNCTION: ,i3/1H ,16HRMS error:      ,e13.6/  &
    1H ,16HABSOLUTE error: ,e13.6/1H ,16HMAXIMUM error:  ,e13.6 ///)
9002   FORMAT(////1H ,22HTESTC 4.3 run complete//)

!***********************************************************************

contains 

    FUNCTION func(x,y,nfunc)
        !    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977
        !    FUNCTION FUNC 4.1 DATED 8 JANUARY 1981

        !    This routine supplies the value of one of eight
        !    test functions at the point (X,Y).

        IF(nfunc <= 0 .OR. nfunc > 8) STOP 616

        SELECT CASE ( nfunc )
          CASE (    1)
            
            !    Function 1: Cliff edge.
            func = 0.0
            IF(x > 0.5) func = 1.0

          CASE (    2)
            
            !    Function 2: V-shaped valley.
            func = ABS(x-y)

          CASE (    3)
            
            !    Function 3: Spherical quadratic (parabolic function).
            func = (x-0.25)**2+(y-0.25)**2

          CASE (    4)
            
            !    Function 4: Non-spherical quadratic (hyperbolic paraboloid).
            func = (x-0.5)**2-(y-0.5)**2

          CASE (    5)
            
            !    Function 5: Spherical bivariate normal density function.
            func = EXP(-8.0*((x-0.75)**2+(y-0.75)**2))

          CASE (    6)
            
            !    Function 6: Sum of two normal densities to give a bimodal
            !                density function.
            func = EXP(-16.0*((x-0.75)**2+(y-0.25)**2))
            func = func+EXP(-16.0*((x-0.25)**2+(y-0.75)**2))

          CASE (    7)
            
            !    Function 7: Fan of ripples.
            func = SIN(6.283185*x*(y+1.0))

          CASE (    8)
            
            !    Function 8: Concentric circular ripples.
            func = COS(12.566371*SQRT((x-0.25)**2+(y-0.25)**2))

        END SELECT

        RETURN

    END FUNCTION func


end program 