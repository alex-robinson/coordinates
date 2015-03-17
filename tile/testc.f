C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    PROGRAM TESTC 4.3 DATED 26 MARCH 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C 
C    This program tests the TILE 4 routines 
C        SCATRC 
C        TILE4M 
C        GRLD 
C        RECT1G 
C    and the routines that they call. 
C 
C    It constructs a tessellation of a random point pattern 
C    on a rectangle, establishes the values of some test 
C    functions at the data points and then reconstructs the 
C    functions over a smaller rectangle using the C1 natural 
C    neighbour interpolant.  The values of the original 
C    functions and their reconstructed values are written 
C    to a file without carriage control characters for 
C    later plotting.  This can be supressed by setting the variable 
C    NSTORE zero or negative.  Statistics on the accuracy of the 
C    reconstructions are written to another file with carriage 
C    control characters. 
C 
C2      INTEGER*2 L 
        EXTERNAL UNI4M 
        DIMENSION CN(3,4),JADDR(4),PT(2,400),NADDR(400),L(6000), 
     1    VAL(400),SBAREA(50),PTOFF(2,50),VALOFF(50),GRAD(2,400), 
     2    DELSBA(2,50),Z(25,25),ZX(25,25),ZY(25,25) 
C 
        DATA JCNS,NPTS,LTOP,KNBMAX,MM,NN/4,400,6000,50,25,25/ 
        DATA XMIN1,XMAX1,YMIN1,YMAX1/-0.5,1.5,-0.5,1.5/ 
        DATA XMIN2,XMAX2,YMIN2,YMAX2/0.0,1.0,0.0,1.0/ 
        DATA MSEED/97129/ 
        DATA NFMAX/8/ 
        DATA EPSCN,EPSPT/0.0,0.0/ 
        DATA NSTART/0/ 
        DATA NWRITE,NSTORE/12,13/ 
C 
C    Record NFMAX, MM, NN and the smaller rectangle 
C    on NSTORE to be used for latter plotting if required. 
C 
        IF(NSTORE.LE.0) GOTO 8 
        WRITE(NSTORE,8000)NFMAX,MM,NN 
        WRITE(NSTORE,8001)XMIN2,XMAX2,YMIN2,YMAX2 
C 
C    Calculate the values needed to scan the grid on the 
C    smaller rectangle. 
C 
    8   DELTAX = (XMAX2-XMIN2)/FLOAT(MM) 
        DELTAY = (YMAX2-YMIN2)/FLOAT(NN) 
        XNUM = FLOAT(MM*NN) 
        XSTART = XMIN2-0.5*DELTAX 
        YSTART = YMIN2-0.5*DELTAY 
C 
C    Generate NPTS points on the larger rectangle and tessellate 
C    them. 
C 
        MSD = MSEED 
        CALL SCATRC(CN,JCNS,PT,NPTS,IFLAG,UNI4M,MSD,XMIN1,XMAX1, 
     1    YMIN1,YMAX1) 
        IF(IFLAG.EQ.0) GOTO 1 
        WRITE(NWRITE,9000)IFLAG 
        STOP 612 
    1   CALL TILE4M(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB) 
        IF(IFLAG.EQ.0) GOTO 2 
        WRITE(NWRITE,9000)IFLAG 
        STOP 613 
C 
C    Loop through the NFMAX functions. 
C 
    2   DO 7  NFUNC = 1,NFMAX 
C 
C    Set the function values at the data sites in VAL. 
C 
        DO 3  NPT = 1,NPTS 
    3   VAL(NPT) = FUNC(PT(1,NPT),PT(2,NPT),NFUNC) 
C 
C    Call GRLD to estimate the gradients at the data sites. 
C 
        CALL GRLD(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VAL,SBAREA, 
     1    PTOFF,VALOFF,KNBMAX,GRAD) 
        IF(IFLAG.EQ.0) GOTO 4 
        WRITE(NWRITE,9000)IFLAG 
        STOP 614 
C 
C    Call RECT1G to compute the C1 interpolant for the function 
C    at the nodes of the MM*NN grid on the smaller rectangle. 
C 
    4   CALL RECT1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,SBAREA,DELSBA, 
     2    PTOFF,KNBMAX,VAL,GRAD,Z,ZX,ZY,MM,MM,NN,XMIN2,XMAX2,YMIN2, 
     3    YMAX2,ZLO,ZHI) 
        IF(IFLAG.EQ.0) GOTO 5 
        WRITE(NWRITE,9000)IFLAG 
        STOP 615 
C 
C    Record the original function and its reconstruction on NSTORE 
C    if required and compute error statistics. 
C 
    5   RMSE = 0.0 
        ABSE = 0.0 
        AMXE = 0.0 
        XPT = XSTART 
        DO 6  M = 1,MM 
        YPT = YSTART 
        XPT = XPT+DELTAX 
        DO 6  N = 1,NN 
        YPT = YPT+DELTAY 
        F = FUNC(XPT,YPT,NFUNC) 
        FZDIF = ABS(Z(M,N)-F) 
        RMSE = RMSE+FZDIF*FZDIF 
        ABSE = ABSE+FZDIF 
        AMXE = AMAX1(AMXE,FZDIF) 
        IF(NSTORE.LE.0) GOTO 6 
        WRITE(NSTORE,8001) XPT, YPT, F,Z(M,N),ZX(M,N),ZY(M,N) 
    6   CONTINUE 
C 
C    Write out the error statistics. 
C 
        RMSE = SQRT(RMSE/XNUM) 
        ABSE = ABSE/XNUM 
    7   WRITE(NWRITE,9001)NFUNC,RMSE,ABSE,AMXE 
C 
C    Run complete. 
C 
        WRITE(NWRITE,9002) 
        STOP 
 8000   FORMAT(3(I13,2X)) 
 8001   FORMAT(6(E13.6,2X)) 
 9000   FORMAT(///1H ,19(1H*)/1H ,16H ERROR: IFLAG = ,I3/1H ,19(1H*)///) 
 9001   FORMAT(///1H ,10HFUNCTION: ,I3/1H ,16HRMS ERROR:      ,E13.6/ 
     1    1H ,16HABSOLUTE ERROR: ,E13.6/1H ,16HMAXIMUM ERROR:  ,E13.6 
     2    ///) 
 9002   FORMAT(////1H ,22HTESTC 4.3 RUN COMPLETE//) 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    FUNCTION FUNC 4.1 DATED 8 JANUARY 1981 
C 
        FUNCTION FUNC(X,Y,NFUNC) 
C 
C    This routine supplies the value of one of eight 
C    test functions at the point (X,Y). 
C 
        IF(NFUNC.LE.0.OR.NFUNC.GT.8) STOP 616 
        GOTO (1,2,3,4,5,6,7,8),NFUNC 
C 
C    Function 1: Cliff edge. 
C 
    1   FUNC = 0.0 
        IF(X.GT.0.5) FUNC = 1.0 
        RETURN 
C 
C    Function 2: V-shaped valley. 
C 
    2   FUNC = ABS(X-Y) 
        RETURN 
C 
C    Function 3: Spherical quadratic (parabolic function). 
C 
    3   FUNC = (X-0.25)**2+(Y-0.25)**2 
        RETURN 
C 
C    Function 4: Non-spherical quadratic (hyperbolic paraboloid). 
C 
    4   FUNC = (X-0.5)**2-(Y-0.5)**2 
        RETURN 
C 
C    Function 5: Spherical bivariate normal density function. 
C 
    5   FUNC = EXP(-8.0*((X-0.75)**2+(Y-0.75)**2)) 
        RETURN 
C 
C    Function 6: Sum of two normal densities to give a bimodal 
C                density function. 
C 
    6   FUNC = EXP(-16.0*((X-0.75)**2+(Y-0.25)**2)) 
        FUNC = FUNC+EXP(-16.0*((X-0.25)**2+(Y-0.75)**2)) 
        RETURN 
C 
C    Function 7: Fan of ripples. 
C 
    7   FUNC = SIN(6.283185*X*(Y+1.0)) 
        RETURN 
C 
C    Function 8: Concentric circular ripples. 
C 
    8   FUNC = COS(12.566371*SQRT((X-0.25)**2+(Y-0.25)**2)) 
        RETURN 
        END 
