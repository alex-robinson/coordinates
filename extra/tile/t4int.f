C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    INTERPOLATION ROUTINES LAST MODIFIED 13 MARCH 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RECT1G 4.1 DATED 13 MARCH 1981 
C 
        SUBROUTINE RECT1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,SBAREA, 
     2    DELSBA,PTOFF,KNBMAX,VAL,GRAD,Z,ZX,ZY,MM,M,N, 
     3    XMIN,XMAX,YMIN,YMAX,ZLO,ZHI) 
C 
C    Evaluates the C1 natural neighbour interpolant and its gradient 
C    at all points on a rectangular grid. 
C 
C    User provides the TILE4 database; working arrays SBAREA(KNBMAX), 
C    DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), and KNBMAX to dimension them; 
C    arrays VAL(NPTS),GRAD(2,NPTS) holding data values and gradients at 
C    data sites, with GRAD usually set by a call to GRLD; 
C    arrays Z(MM,N),ZX(MM,N),ZY(MM,N) in which the routine 
C    returns value and components of gradient at the grid points, with 
C    MM as the true first dimension of these arrays, and M,N as the 
C    numbers of grid points in the X and Y directions respectively; 
C    XMIN,XMAX,YMIN,YMAX to specify the position of the grid; and ZLO, 
C    ZHI in which the routine returns the minimum and maximum values of 
C    the interpolant on the grid.  The position and spacing of the grid 
C    is defined as follows.  The rectangle ÝXMIN,XMAX¨*ÝYMIN,YMAX¨ is 
C    divided into M*N small rectangles each (XMAX-XMIN)/FLOAT(M) by 
C    (YMAX-YMIN)/FLOAT(N); the interpolant is evaluated at the centres 
C    of these.  Thus grid point (MSCAN,NSCAN) is at 
C 
C       (XMIN+(FLOAT(MSCAN)-0.5)*(XMAX-XMIN)/FLOAT(M), 
C          YMIN+(FLOAT(NSCAN)-0.5)*(YMAX-YMIN)/FLOAT(N)) 
C 
C    The grid need not lie strictly within the window; at grid points 
C    outside the window, ZX and ZY are set conventionally to zero, and 
C    Z is set to 2.0*ZLO-ZHI to allow such points to be detected without 
C    the need for an array of flags.  IFLAG is zero on successful 
C    return.  Nonzero values indicate error or nonstandard returns as 
C    follows. 
C 
C       4  Some grid point outside window 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    Note that heap overflow is not here a catastrophic error.  A return 
C    value of 4 for IFLAG is not an error in this context - it simply 
C    alerts the user to the fact that his grid is partly outside the 
C    window.  If grid points ever duplicate data sites, a value of 5 for 
C    IFLAG is generated internally by NNBR1G, but this again is not an 
C    error in this context; value and gradient are set correctly, and 
C    the flag is not passed back to the user. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS), 
     2    GRAD(2,NPTS),Z(MM,N),ZX(MM,N),ZY(MM,N) 
C 
C    Initialise values. 
C 
        ZLO = 1.0 
        ZHI = 0.0 
        IFLAG = 0 
        NS2 = NSTART 
        DELTAX = (XMAX-XMIN)/FLOAT(M) 
        DELTAY = (YMAX-YMIN)/FLOAT(N) 
        X0 = XMIN-0.5*DELTAX 
        Y0 = YMIN-0.5*DELTAY 
        XSCAN = X0 
C 
C     Scan the rectangular grid. 
C 
        DO 1  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        IF(NS2.EQ.0) NS2 = NSTART 
        NS1 = NS2 
        NS2 = 0 
        DO 1  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL NNBR1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NS1,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFL,IGARB,XSCAN,YSCAN,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL,GRAD, 
     3    Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
        write(*,"(2f8.2,i6,2f8.2)") XSCAN, YSCAN, NINDEX, Z1, 
     1  Z1-VAL(NINDEX)
C 
C    Check IFL and take appropriate action. 
C 
        IF(IFL.EQ.0.OR.IFL.EQ.5) GOTO 2 
        IFLAG = IFL 
        IF(IFL.EQ.4) GOTO 1 
        RETURN 
C 
C    We are dealing with a grid point inside the window. 
C 
    2   Z(MSCAN,NSCAN) = Z1 
        ZX(MSCAN,NSCAN) = ZX1 
        ZY(MSCAN,NSCAN) = ZY1 
        IF(ZLO.LE.ZHI) GOTO 3 
        ZLO = Z1 
        ZHI = Z1 
        GOTO 4 
    3   IF(ZLO.GT.Z1) ZLO = Z1 
        IF(ZHI.LT.Z1) ZHI = Z1 
    4   IF(NS2.EQ.0) NS2 = NINDEX 
        NS1 = NINDEX 
    1   CONTINUE 
C 
C    Deal with the flagging of points outside the window if any. 
C 
        IF(IFLAG.EQ.0) RETURN 
        ZVAL = 2.0*ZLO-ZHI 
        XSCAN = X0 
        DO 5  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        DO 5  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFL,XSCAN,YSCAN,JINDEX,TSTVAL) 
        IF(IFL.EQ.0) GOTO 5 
        Z(MSCAN,NSCAN) = ZVAL 
        ZX(MSCAN,NSCAN) = 0.0 
        ZY(MSCAN,NSCAN) = 0.0 
    5   CONTINUE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RECTHG 4.1 DATED 13 MARCH 1981 
C 
        SUBROUTINE RECTHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,SBAREA, 
     2    DELSBA,PTOFF,KNBMAX,VAL,Z,ZX,ZY,MM,M,N, 
     3    XMIN,XMAX,YMIN,YMAX,ZLO,ZHI) 
C 
C    Evaluates the C1 natural neighbour interpolant and its gradient 
C    at all points on a rectangular grid, with the gradients 
C    at the data sites all forced to be zero and GRAD 
C    accordingly not required.  Other comments are as for 
C    RECT1G. 
C 
C    User provides the TILE4 database; working arrays SBAREA(KNBMAX), 
C    DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), and KNBMAX to dimension them; 
C    array VAL(NPTS) holding values at data sites; 
C    arrays Z(MM,N),ZX(MM,N),ZY(MM,N) in which the routine 
C    returns value and components of gradient at the grid points, with 
C    MM as the true first dimension of these arrays, and M,N as the 
C    numbers of grid points in the X and Y directions respectively; 
C    XMIN,XMAX,YMIN,YMAX to specify the position of the grid; and ZLO, 
C    ZHI in which the routine returns the minimum and maximum values of 
C    the interpolant on the grid.  The position and spacing of the grid 
C    is defined as follows.  The rectangle ÝXMIN,XMAX¨*ÝYMIN,YMAX¨ is 
C    divided into M*N small rectangles each (XMAX-XMIN)/FLOAT(M) by 
C    (YMAX-YMIN)/FLOAT(N); the interpolant is evaluated at the centres 
C    of these.  Thus grid point (MSCAN,NSCAN) is at 
C 
C       (XMIN+(FLOAT(MSCAN)-0.5)*(XMAX-XMIN)/FLOAT(M), 
C          YMIN+(FLOAT(NSCAN)-0.5)*(YMAX-YMIN)/FLOAT(N)) 
C 
C    The grid need not lie strictly within the window; at grid points 
C    outside the window, ZX and ZY are set conventionally to zero, and 
C    Z is set to 2.0*ZLO-ZHI to allow such points to be detected without 
C    the need for an array of flags.  IFLAG is zero on successful 
C    return.  Nonzero values indicate error or nonstandard returns as 
C    follows. 
C 
C       4  Some grid point outside window 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    Note that heap overflow is not here a catastrophic error.  A return 
C    value of 4 for IFLAG is not an error in this context - it simply 
C    alerts the user to the fact that his grid is partly outside the 
C    window.  If grid points ever duplicate data sites, a value of 5 for 
C    IFLAG is generated internally by NNBRHG, but this again is not an 
C    error in this context; value and gradient are set correctly, and 
C    the flag is not passed back to the user. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS), 
     2    Z(MM,N),ZX(MM,N),ZY(MM,N) 
C 
C    Initialise values. 
C 
        ZLO = 1.0 
        ZHI = 0.0 
        IFLAG = 0 
        NS2 = NSTART 
        DELTAX = (XMAX-XMIN)/FLOAT(M) 
        DELTAY = (YMAX-YMIN)/FLOAT(N) 
        X0 = XMIN-0.5*DELTAX 
        Y0 = YMIN-0.5*DELTAY 
        XSCAN = X0 
C 
C     Scan the rectangular grid. 
C 
        DO 1  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        IF(NS2.EQ.0) NS2 = NSTART 
        NS1 = NS2 
        NS2 = 0 
        DO 1  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL NNBRHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NS1,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFL,IGARB,XSCAN,YSCAN,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL, 
     3    Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
C 
C    Check IFL and take appropriate action. 
C 
        IF(IFL.EQ.0.OR.IFL.EQ.5) GOTO 2 
        IFLAG = IFL 
        IF(IFL.EQ.4) GOTO 1 
        RETURN 
C 
C    We are dealing with a grid point inside the window. 
C 
    2   Z(MSCAN,NSCAN) = Z1 
        ZX(MSCAN,NSCAN) = ZX1 
        ZY(MSCAN,NSCAN) = ZY1 
        IF(ZLO.LE.ZHI) GOTO 3 
        ZLO = Z1 
        ZHI = Z1 
        GOTO 4 
    3   IF(ZLO.GT.Z1) ZLO = Z1 
        IF(ZHI.LT.Z1) ZHI = Z1 
    4   IF(NS2.EQ.0) NS2 = NINDEX 
        NS1 = NINDEX 
    1   CONTINUE 
C 
C    Deal with the flagging of points outside the window if any. 
C 
        IF(IFLAG.EQ.0) RETURN 
        ZVAL = 2.0*ZLO-ZHI 
        XSCAN = X0 
        DO 5  MSCAN = 1,M 
        XSCAN = XSCAN+DELTAX 
        YSCAN = Y0 
        DO 5  NSCAN = 1,N 
        YSCAN = YSCAN+DELTAY 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFL,XSCAN,YSCAN,JINDEX,TSTVAL) 
        IF(IFL.EQ.0) GOTO 5 
        Z(MSCAN,NSCAN) = ZVAL 
        ZX(MSCAN,NSCAN) = 0.0 
        ZY(MSCAN,NSCAN) = 0.0 
    5   CONTINUE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR1G 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE NNBR1G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,GRAD,Z0,ZX0,ZY0,Z1,ZX1,ZY1) 
C 
C    Calculates the C0 and C1 natural neighbour interpolants and 
C    their gradients at the point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolants are wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) as working space, 
C    and KNBMAX to dimension them, arrays VAL(NPTS),GRAD(2,NPTS) 
C    holding the data values and gradients at the data sites.  GRAD 
C    will usually have been loaded by a call to GRLD.  Z0 returns the 
C    C0 natural neighbour interpolant, (ZX0,ZY0) its gradient; Z1 
C    returns the C1 natural neighbour interpolant, (ZX1,ZY1) its 
C    gradient.  IFLAG is zero on successful return.  Nonzero values 
C    indicate error or nonstandard returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    "Outside" and "duplicate" are interpreted in terms of the 
C    tolerance values EPSCN,EPSPT.  If duplication occurs, NINDEX 
C    returns the index of the duplicated point.  Z0 returns the 
C    value at NINDEX, and (ZX0,ZY0) conventionally return (0.0,0.0). 
C    (ZX1,ZY1) return the gradient at NINDEX, and Z1 returns the 
C    value at (XPT,YPT) on the plane through Z0 at NINDEX with 
C    gradient (ZX1,ZY1).  Note that heap overflow is not here a 
C    catastrophic error.  ICASE returns the edge effect indicator 
C    from the call to TRYSBG, that is, 5 at the boundary and 0 
C    elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), 
     2    VAL(NPTS),GRAD(2,NPTS) 
C 
C    Call TRYSBG to calculate subareas and their gradients. 
C 
        CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   UPT = XPT-PT(1,NINDEX) 
        VPT = YPT-PT(2,NINDEX) 
        Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        ZX1 = GRAD(1,NINDEX) 
        ZY1 = GRAD(2,NINDEX) 
        Z1 = Z0+ZX1*UPT+ZY1*VPT 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        S1X = 0.0 
        S1Y = 0.0 
        S2X = 0.0 
        S2Y = 0.0 
        SM1X = 0.0 
        SM1Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
        TM1X = 0.0 
        TM1Y = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        GX = 0.0 
        GY = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
        GXX = 0.0 
        GXY = 0.0 
        GYX = 0.0 
        GYY = 0.0 
        EDGEXX = 0.0 
        EDGEXY = 0.0 
        EDGEYX = 0.0 
        EDGEYY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up and calculate values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT*UPT+VPT*VPT 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
        GRX = GRAD(1,N) 
        GRY = GRAD(2,N) 
        ZETA = ZPT+GRX*UPT+GRY*VPT 
C 
C    Accumulate main values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZETA/DPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        S1X = S1X+SX*DPT+S*UPT/DPT 
        S1Y = S1Y+SY*DPT+S*VPT/DPT 
        S2X = S2X+SX*DSQ 
        S2Y = S2Y+SY*DSQ 
        SM1X = SM1X+(SX-S*UPT/DSQ)/DPT 
        SM1Y = SM1Y+(SY-S*VPT/DSQ)/DPT 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
        TM1X = TM1X+(SX*ZETA+S*(GRX-ZETA*UPT/DSQ))/DPT 
        TM1Y = TM1Y+(SY*ZETA+S*(GRY-ZETA*VPT/DSQ))/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        GX = GX+S*GRX 
        GY = GY+S*GRY 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
        GXX = GXX+SX*GRX 
        GXY = GXY+SY*GRX 
        GYX = GYX+SX*GRY 
        GYY = GYY+SY*GRY 
        EDGEXX = EDGEXX+SX*UPT+S 
        EDGEXY = EDGEXY+SY*UPT 
        EDGEYX = EDGEYX+SX*VPT 
        EDGEYY = EDGEYY+SY*VPT+S 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        IF(ICASE.EQ.0) GOTO 5 
        S2X = S2X+2.0*EDGEX 
        S2Y = S2Y+2.0*EDGEY 
        CORR = (GX*EDGEX+GY*EDGEY)/S0 
        T0 = T0+CORR 
        T0X = T0X+(GXX*EDGEX+GX*EDGEXX+GYX*EDGEY+GY*EDGEYX-CORR*S0X)/S0 
        T0Y = T0Y+(GXY*EDGEX+GX*EDGEXY+GYY*EDGEY+GY*EDGEYY-CORR*S0Y)/S0 
    5   Z1DEN = S1*S0+S2*SM1 
        Z1 = (S1*T0+S2*TM1)/Z1DEN 
        ZX1 = S1X*T0+S1*T0X+S2X*TM1+S2*TM1X 
        ZX1 = (ZX1-Z1*(S1X*S0+S1*S0X+S2X*SM1+S2*SM1X))/Z1DEN 
        ZY1 = S1Y*T0+S1*T0Y+S2Y*TM1+S2*TM1Y 
        ZY1 = (ZY1-Z1*(S1Y*S0+S1*S0Y+S2Y*SM1+S2*SM1Y))/Z1DEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 DATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBRHG 4.1 DATED 10 MARCH 1981 
C 
        SUBROUTINE NNBRHG00(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,Z0,ZX0,ZY0,ZH,ZXH,ZYH) 
C
C    Initial data values: PT, NPTS, VAL(NPTS),GRAD(2,NPTS)
C    CN = 
C    Working arrays: SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), 
C    Working dimension: KNBMAX
C    Calculates the C0 natural neighbour interpolant and its gradient, 
C    and also the C1 natural neighbour interpolant and its gradient with 
C    all data site gradients forced to zero, at the point (XPT,YPT). 
C    This routine thus produces the same effect as a call to NNBR1G 
C    with zeroes entered into GRAD, but the array GRAD does not need 
C    to be passed to it and the calculation is more efficient.  The 
C    C1 interpolant value is returned as ZH, its gradient as (ZXH,ZYH). 
C    Other details are as for NNBR1G. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 
C 
C    Call TRYSBG to calculate subareas and their gradients. 
C 
       CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        ZH = Z0 
        ZXH = 0.0 
        ZYH = 0.0 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        S1X = 0.0 
        S1Y = 0.0 
        S2X = 0.0 
        S2Y = 0.0 
        SM1X = 0.0 
        SM1Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
        TM1X = 0.0 
        TM1Y = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        EDGEX = 0.0 
        EDGEY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT**2+VPT**2 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZPT/DPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        S1X = S1X+SX*DPT+S*UPT/DPT 
        S1Y = S1Y+SY*DPT+S*VPT/DPT 
        S2X = S2X+SX*DSQ 
        S2Y = S2Y+SY*DSQ 
        SM1X = SM1X+(SX-S*UPT/DSQ)/DPT 
        SM1Y = SM1Y+(SY-S*VPT/DSQ)/DPT 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
        TM1X = TM1X+(SX-S*UPT/DSQ)*ZPT/DPT 
        TM1Y = TM1Y+(SY-S*VPT/DSQ)*ZPT/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        IF(ICASE.EQ.0) GOTO 5 
        S2X = S2X+2.0*EDGEX 
        S2Y = S2Y+2.0*EDGEY 
    5   ZHDEN = S1*S0+S2*SM1 
        ZH = (S1*T0+S2*TM1)/ZHDEN 
        ZXH = S1X*T0+S1*T0X+S2X*TM1+S2*TM1X 
        ZXH = (ZXH-ZH*(S1X*S0+S1*S0X+S2X*SM1+S2*SM1X))/ZHDEN 
        ZYH = S1Y*T0+S1*T0Y+S2Y*TM1+S2*TM1Y 
        ZYH = (ZYH-ZH*(S1Y*S0+S1*S0Y+S2Y*SM1+S2*SM1Y))/ZHDEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR1 4.2 DATED 13 MARCH 1981 
C 
        SUBROUTINE NNBR1(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE, 
     3    VAL,GRAD,Z0,Z1) 
C 
C    Calculates the C0 and C1 natural neighbour interpolants at the 
C    point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolants are wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX) as working space, and KNBMAX to 
C    dimension them, arrays VAL(NPTS),GRAD(2,NPTS) holding the data 
C    values and gradients at the data sites.  GRAD will usually have 
C    been loaded by a call to GRLD.  Z0 returns the C0 natural 
C    neighbour interpolant, Z1 the C1 natural neighbour interpolant. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    or nonstandard returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNBMAX too small 
C 
C    "Outside" and "duplicate" are interpreted in terms of the 
C    tolerance values EPSCN,EPSPT.  If duplication occurs, NINDEX 
C    returns the index of the duplicated point, and Z0,Z1 return 
C    correctly.  Note that heap overflow is not here a catastrophic 
C    error.  ICASE returns the edge effect indicator value from the 
C    call to TRYSBA, that is, 5 at the boundary and 0 elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS),GRAD(2,NPTS) 
C 
C    Call TRYSBA to calculate subareas. 
C 
        CALL TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is zero or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   UPT = XPT-PT(1,NINDEX) 
        VPT = YPT-PT(2,NINDEX) 
        Z0 = VAL(NINDEX) 
        Z1 = Z0+GRAD(1,NINDEX)*UPT+GRAD(2,NINDEX)*VPT 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators for the non-edge case. 
C 
    1   S0 = 0.0 
        S1 = 0.0 
        S2 = 0.0 
        SM1 = 0.0 
        T0 = 0.0 
        TM1 = 0.0 
C 
C    If needed, initialise accumulators for edge correction. 
C 
        IF(ICASE.EQ.0) GOTO 3 
        GX = 0.0 
        GY = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up and calculate values 
C    needed to update the accumulators. 
C 
    3   LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 4  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 4 
        UPT = -PTOFF(1,K) 
        VPT = -PTOFF(2,K) 
        DSQ = UPT**2+VPT**2 
        DPT = SQRT(DSQ) 
        S = SBAREA(K) 
        ZPT = VAL(N) 
        GRX = GRAD(1,N) 
        GRY = GRAD(2,N) 
        ZETA = ZPT+GRX*UPT+GRY*VPT 
C 
C    Accumulate main values. 
C 
        S0 = S0+S 
        S1 = S1+S*DPT 
        S2 = S2+S*DSQ 
        SM1 = SM1+S/DPT 
        T0 = T0+S*ZPT 
        TM1 = TM1+S*ZETA/DPT 
C 
C    Accumulate edge corrector values if needed. 
C 
        IF(ICASE.EQ.0) GOTO 4 
        GX = GX+S*GRX 
        GY = GY+S*GRY 
        EDGEX = EDGEX+S*UPT 
        EDGEY = EDGEY+S*VPT 
    4   CONTINUE 
C 
C    Calculate values for return, making edge corrections if needed. 
C 
        Z0 = T0/S0 
        IF(ICASE.EQ.0) GOTO 5 
        CORR = (GX*EDGEX+GY*EDGEY)/S0 
        T0 = T0+CORR 
    5   Z1DEN = S1*S0+S2*SM1 
        Z1 = (S1*T0+S2*TM1)/Z1DEN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR0G 4.3 DATED 9 MARCH 1981 
C 
        SUBROUTINE NNBR0G(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,VAL, 
     3    Z0,ZX0,ZY0) 
C 
C    Calculates the C0 natural neighbour interpolant and its gradient 
C    at the point (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolant is wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) as working space, 
C    and KNBMAX to dimension them, and an array VAL(NPTS) holding the 
C    data values at the data sites.  Z0 returns the C0 natural 
C    neighbour interpolant, (ZX0,ZY0) its gradient.  IFLAG is zero on 
C    successful return.  Nonzero values indicate error or nonstandard 
C    returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If duplication occurs, NINDEX returns the 
C    index of the duplicated point.  Z0 returns the value at NINDEX. 
C    (ZX0,ZY0) conventionally return (0.0,0.0).  Note that heap overflow 
C    is not here a catastrophic error.  ICASE returns the edge effect 
C    indicator flag from the call to TRYSBG, that is, 5 at the boundary 
C    and 0 elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 
C 
C    Call TRYSBG to calculate subareas and gradients. 
C 
        CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX) 
        ZX0 = 0.0 
        ZY0 = 0.0 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators.  There is no need to worry about the edge case, 
C    because we do not try to compensate in the C0 natural neighbour 
C    interpolant. 
C 
    1   S0 = 0.0 
        T0 = 0.0 
        S0X = 0.0 
        S0Y = 0.0 
        T0X = 0.0 
        T0Y = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
        LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 3 
        S = SBAREA(K) 
        SX = DELSBA(1,K) 
        SY = DELSBA(2,K) 
        ZPT = VAL(N) 
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        T0 = T0+S*ZPT 
        S0X = S0X+SX 
        S0Y = S0Y+SY 
        T0X = T0X+SX*ZPT 
        T0Y = T0Y+SY*ZPT 
    3   CONTINUE 
C 
C    Calculate values for return. 
C 
        Z0 = T0/S0 
        ZX0 = (T0X-Z0*S0X)/S0 
        ZY0 = (T0Y-Z0*S0Y)/S0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NNBR0 4.2 DATED 12 MARCH 1981 
C 
        SUBROUTINE NNBR0(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE,VAL,Z0) 
C 
C    Calculates the C0 natural neighbour interpolant at the point 
C    (XPT,YPT). 
C 
C    User provides the TILE4 database, the coordinates (XPT,YPT) of the 
C    point at which the interpolant is wanted, a variable NINDEX in 
C    which the index of the (a) nearest data site is returned, arrays 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX) as working space, and KNBMAX to 
C    dimension them, and an array VAL(NPTS) holding the data values at 
C    the data sites.  Z0 returns the C0 natural neighbour interpolant. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    or nonstandard returns as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If duplication occurs, NINDEX returns the 
C    index of the duplicated point, and Z0 returns the value at NINDEX. 
C    Note that heap overflow is not here a catastrophic error.  ICASE 
C    returns the edge effect indicator flag from the call to TRYSBA, 
C    that is, 5 at the boundary and 0 elsewhere. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VAL(NPTS) 
C 
C    Call TRYSBA to calculate subareas. 
C 
        CALL TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, 
     2    AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Check IFLAG and return unless it is 0 or 5. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        IF(IFLAG.EQ.5) GOTO 2 
        RETURN 
C 
C    Deal with the case where we have hit an accepted point. 
C 
    2   Z0 = VAL(NINDEX) 
        RETURN 
C 
C    Deal with the general case.  Initialise to zero all the 
C    accumulators.  There is no need to worry about the edge case, 
C    because we do not try to compensate in the C0 natural neighbour 
C    interpolant. 
C 
    1   S0 = 0.0 
        T0 = 0.0 
C 
C    Pick up the contiguity list, and enter a loop through the 
C    neighbours of the trial point.  Pick up values 
C    needed to update the accumulators. 
C 
        LLO = LPTR+2 
        LHI = LLO-1+L(LLO-1) 
        K = 0 
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N.LE.0) GOTO 3 
        S = SBAREA(K) 
        ZPT = VAL(N) 
C 
C    Accumulate values. 
C 
        S0 = S0+S 
        T0 = T0+S*ZPT 
    3   CONTINUE 
C 
C    Calculate value for return. 
C 
        Z0 = T0/S0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GRLD 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE GRLD(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VAL, 
     1    SBAREA,PTOFF,VALOFF,KNBMAX,GRAD) 
C 
C    Loads the gradient at each data site into GRAD. 
C 
C    User provides the appropriate part of the TILE4 database, and 
C    the values at the accepted points (data sites) in VAL. 
C    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VALOFF(KNBMAX) are needed as 
C    working space. 
C 
C    If KNBMAX is too small, error return with IFLAG set to 10 occurs. 
C    Otherwise IFLAG is zero on return.  For accepted point NPT, the 
C    gradient is loaded into (GRAD(1,NPT),GRAD(2,NPT)). 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP),VAL(NPTS), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX),VALOFF(KNBMAX),GRAD(2,NPTS) 
C 
C    Scan the points. 
C 
        DO 1  NPT = 1,NPTS 
C 
C    Calculate subareas. 
C 
        CALL ACCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NPT,AREA,SBAREA, 
     1    PTOFF,KNBMAX,KNB) 
C 
C    Check IFLAG.  Return if 10, loop if 9. 
C 
        IF(IFLAG.EQ.0) GOTO 2 
        IF(IFLAG.EQ.9) GOTO 1 
        RETURN 
C 
C    Calculate gradient. 
C 
    2   CALL CURLYD(NPTS,NADDR,L,LTOP,IFLAG,VAL,NPT,SBAREA,PTOFF,VALOFF, 
     1    KNB,BETAX,BETAY,GAMMA,ICASE) 
C 
C    Check IFLAG to make sure nothing has gone wrong. 
C 
        IF(IFLAG.EQ.0) GOTO 3 
        STOP 766 
C 
C    Save gradient. 
C 
    3   GRAD(1,NPT) = BETAX 
        GRAD(2,NPT) = BETAY 
    1   CONTINUE 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CURLYD 4.4 DATED 13 JANUARY 1981 
C 
        SUBROUTINE CURLYD(NPTS,NADDR,L,LTOP,IFLAG,VAL,NINDEX, 
     1    SBAREA,PTOFF,VALOFF,KNB,BETAX,BETAY,GAMMA,ICASE) 
C 
C    Applies the discrete differentiator CURLY-D at NINDEX. 
C 
C    NPTS,NADDR,L,LTOP are from the TILE4 database.  VAL contains values 
C    at accepted points.  NINDEX is the index of the point (data site) 
C    at which CURLY-D is to be applied.  KNB is the number of its 
C    neighbours, and for K = 1,...,KNB the subtile area for the Kth 
C    neighbour is held in SBAREA(K) and the offset to it as a real or 
C    virtual point in (PTOFF(1,K),PTOFF(2,K)).  KNB,SBAREA,PTOFF can be 
C    loaded by a call to ACCSBA.  VALOFF(K) returns the real or virtual 
C    value offset at the Kth neighbour, (BETAX,BETAY) the gradient and 
C    GAMMA the spherical coefficient at NINDEX.  IFLAG returns the value 
C    zero unless NINDEX is not the index of an accepted point, in which 
C    case it returns the value 9.  ICASE returns a value indicating the 
C    extent of the edge effects, according to the following conventions. 
C 
C       0  All neighbours of NINDEX are points.  Equations for BETAX, 
C          BETAY,GAMMA assumed well-conditioned.  Problems could arise 
C          only if NINDEX is too close to another point.  This should 
C          have been trapped earlier. 
C 
C       1  Three or more neighbours of NINDEX, but not all, are points 
C          and equations are well-conditioned.  This is the normal edge 
C          case.  The calculation is similar to case 0 but with extra 
C          terms to compensate for the edge effects. 
C 
C       2  Two neighbours of NINDEX are points.  More precisely, the 
C          equations for BETAX,BETAY,GAMMA are ill-conditioned or 
C          indeterminate, but on putting GAMMA equal to zero and 
C          discarding an equation, well-conditioned equations for 
C          BETAX,BETAY result. 
C 
C       3  One neighbour of NINDEX is a point.  More precisely, the 
C          equations for BETAX,BETAY,GAMMA are ill-conditioned or 
C          indeterminate, and on putting GAMMA equal to zero and 
C          discarding an equation, the resultant equations for 
C          BETAX,BETAY are still ill-conditioned or indeterminate. 
C          A direction for beta is constructed conventionally, and an 
C          appropriate magnitude is then determined. 
C 
C       4  No neighbours of NINDEX are points.  This occurs only if 
C          NINDEX is the sole accepted point.  The case is included 
C          for completeness, but will not normally arise. 
C 
C2      INTEGER*2 L 
        DIMENSION NADDR(NPTS),L(LTOP),VAL(NPTS),SBAREA(KNB), 
     1    PTOFF(2,KNB),VALOFF(KNB) 
C 
C    Test value for conditioning. 
C 
        DATA EPS1,EPS2/1.0E-6,1.0E-6/ 
C 
C    Check that NINDEX is the index of an accepted point, and if so pick 
C    up LLO as the base of its contiguity list.  If not, return with 
C    IFLAG set to 9. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    ZINDEX as the value at NINDEX.  Initialise accumulator variables. 
C 
    2   LHI = LLO-1+L(LLO-1) 
        ZINDEX = VAL(NINDEX) 
        AREA = 0.0 
        EDGEX = 0.0 
        EDGEY = 0.0 
        SUMSQ = 0.0 
        Q = 0.0 
        HXX = 0.0 
        HXY = 0.0 
        HYY = 0.0 
        PX = 0.0 
        PY = 0.0 
C 
C    Scan the neighbours, counting how many are points and how many are 
C    constraints.  Load VALOFF for points, and accumulate all necessary 
C    values. 
C 
        NBR = 0 
        JBR = 0 
        K = 0 
        DO 3  LL = LLO,LHI 
        K = K+1 
        N = L(LL) 
        IF(N) 4,701,5 
  701   STOP 767 
    4   JBR = JBR+1 
        GOTO 3 
    5   NBR = NBR+1 
        Z = VAL(N)-ZINDEX 
        VALOFF(K) = Z 
        S = SBAREA(K) 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        AREA = AREA+S 
        EDGEX = EDGEX+S*U 
        EDGEY = EDGEY+S*V 
        SUMSQ = SUMSQ+S*(U*U+V*V) 
        Q = Q+S*Z 
        S = S/(U*U+V*V) 
        HXX = HXX+S*U*U 
        HXY = HXY+S*U*V 
        HYY = HYY+S*V*V 
        PX = PX+S*Z*U 
        PY = PY+S*Z*V 
    3   CONTINUE 
        IF(K.NE.KNB) STOP 770 
C 
C    Split off edge cases, and deal with case 0. 
C 
        IF(JBR.GT.0) GOTO 6 
        ICASE = 0 
        DET = HXX*HYY-HXY*HXY 
        BETAX = (HYY*PX-HXY*PY)/DET 
        BETAY = (HXX*PY-HXY*PX)/DET 
        GAMMA = Q/SUMSQ 
        IFLAG = 0 
        RETURN 
C 
C    Split off bad edge cases, and try for case 1. 
C 
    6   IF(NBR.LE.2) GOTO 7 
        ICASE = 1 
        HXXMOD = HXX-EDGEX*EDGEX/SUMSQ 
        HXYMOD = HXY-EDGEX*EDGEY/SUMSQ 
        HYYMOD = HYY-EDGEY*EDGEY/SUMSQ 
        DET = HXXMOD*HYYMOD-HXYMOD*HXYMOD 
        IF(DET/AREA.LT.EPS1) GOTO 7 
        PXMOD = PX-Q*EDGEX/SUMSQ 
        PYMOD = PY-Q*EDGEY/SUMSQ 
        BETAX = (HYYMOD*PXMOD-HXYMOD*PYMOD)/DET 
        BETAY = (HXXMOD*PYMOD-HXYMOD*PXMOD)/DET 
        GAMMA = (Q-EDGEX*BETAX-EDGEY*BETAY)/SUMSQ 
        GOTO 10 
C 
C    Try for case 2. 
C 
    7   IF(NBR.LE.1) GOTO 8 
        ICASE = 2 
        DET = HXX*HYY-HXY*HXY 
        IF(DET/AREA.LT.EPS2) GOTO 8 
        BETAX = (HYY*PX-HXY*PY)/DET 
        BETAY = (HXX*PY-HXY*PX)/DET 
        GAMMA = 0.0 
        GOTO 10 
C 
C    Case 3. 
C 
    8   IF(NBR.EQ.0) GOTO 9 
        ICASE = 3 
        BETAX = PX/AREA 
        BETAY = PY/AREA 
        GAMMA = 0.0 
        GOTO 10 
C 
C    Case 4. 
C 
    9   ICASE = 4 
        BETAX = 0.0 
        BETAY = 0.0 
        GAMMA = 0.0 
C 
C    Calculate and load virtual value offsets. 
C 
   10   LL = LLO-1 
        DO 11  K = 1,KNB 
        LL = LL+1 
        IF(L(LL).GE.0) GOTO 11 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        VALOFF(K) = BETAX*U+BETAY*V+GAMMA*(U*U+V*V) 
   11   CONTINUE 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYSBG 4.1 DATED 6 MARCH 1981 
C 
        SUBROUTINE TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    For a trial point (XPT,YPT) calculates the contiguities this point 
C    would have if inserted into the tessellation, the offsets to its 
C    neighbours, the area of its tile, and the areas and gradients 
C    of the areas of its subtiles. 
C 
C    The user supplies the TILE4 data structure, the coordinates 
C    (XPT,YPT) of the trial point, unset variables NINDEX,AREA,KNB, 
C    ICASE, and arrays SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX). 
C    NINDEX returns the index of the accepted point nearest to the trial 
C    point.  AREA returns the area of the tile of the trial point. 
C    KNB returns the number of neighbours.  ICASE returns 0 if all these 
C    are points, and 5 if some are constraints.  The data item 
C    containing the contiguity list is returned at LPTR in the heap L, 
C    and neighbours are referenced in the order in which they occur in 
C    this contiguity list.  SBAREA(K) for K = 1,...,KNB returns the area 
C    of the subtile for the Kth neighbour.  (DELSBA(1,K),DELSBA(2,K)) 
C    for K = 1,...,KNB return the components of the gradient of the 
C    area of the subtile for the Kth neighbour.  (PTOFF(1,K),PTOFF(2,K)) 
C    for K = 1,...,KNB return the offsets to the Kth neighbour as a real 
C    or virtual point.  IFLAG is zero on successful return.  Nonzero 
C    values indicate error returns as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If error return with IFLAG set to 5 occurs, 
C    NINDEX correctly returns the index of the duplicated point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Call TRYPT to find the contiguity list for the trial point. 
C 
        CALL TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Check IFLAG and return if it is nonzero. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Call LOCSBG for the contiguity list at LPTR constructed for the 
C    trial point (XPT,YPT), to work out the subareas and gradients. 
C 
    1   CALL LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LPTR,XPT,YPT, 
     1    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBG, and these values can 
C    be passed straight back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ACCSBG 4.1 DATED 6 MARCH 1981 
C 
        SUBROUTINE ACCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB) 
C 
C    Calculates the coordinate offsets of the neighbours of NINDEX, the 
C    areas of their subtiles of the tile of NINDEX, the gradients of 
C    these subtile areas, and the total area of the tile of NINDEX. 
C 
C    Note: this routine is provided as a tool for future use.  It is not 
C    called as part of the standard interpolation procedure in TILE 4. 
C 
C    The user supplies arrays SBAREA(KNBMAX),DELSBA(2,KNBMAX), 
C    PTOFF(2,KNBMAX).  KNB returns the number of neighbours of NINDEX. 
C    If this would exceed KNBMAX, error return with IFLAG set to 10 
C    occurs.  If NINDEX is not the index of an accepted point, error 
C    return with IFLAG set to 9 occurs.  Otherwise IFLAG is zero on 
C    successful return.  AREA returns the total area of the tile of 
C    NINDEX.  (PTOFF(1,K),PTOFF(2,K)) for K = 1,...,KNB return the 
C    coordinate offsets of the neighbours of NINDEX as real or virtual 
C    points, in the order in which they are encountered in its 
C    contiguity list.  SBAREA(K) for K = 1,...,KNB returns the area of 
C    the subtile associated with the Kth neighbour, (DELSBA(1,K), 
C    DELSBA(2,K)) the components of the gradient of its area. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check that NINDEX is the index of an accepted point.  If not, 
C    return with IFLAG set to 9.  If so, LOCN is the address of its 
C    data item in the heap. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN.GT.0) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    2   XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Call subroutine LOCSBG to calculate the subareas and their 
C    gradients. 
C 
        CALL LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XINDEX,YINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBG, and these values can 
C    just be passed back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCSBG 4.2 DATED 6 MARCH 1981 
C 
        SUBROUTINE LOCSBG(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XPT,YPT,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Calculates coordinate offsets to neighbours, tile area, subtile 
C    areas, and gradients of subtile areas for a point (XPT,YPT) whose 
C    data item in the heap is held at LOCN.  The data item may be either 
C    that of an accepted point, or that of a trial point prepared at 
C    location LPTR by a call to TRYPT.  The user should not normally 
C    call subroutine LOCSBG directly, but rather via a call to ACCSBG 
C    to find the subareas and gradients for an accepted point, or via 
C    a call to TRYSBG to find those for a trial point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX) 
C 
C    No check is made on LOCN.  If it is not in fact the address of a 
C    data item, chaos will result.  Assuming this does not happen, pick 
C    up LLO as the base of the contiguity list, and LHI as the top of 
C    it.  Check that the list is not too long. 
C 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        IF(LHI-LLO+1.LE.KNBMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Find the coordinate offsets of all the neighbours.  Initialise 
C    DELSBA to zero. 
C 
    1   KNB = 0 
        DO 2  LL = LLO,LHI 
        KNB = KNB+1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XPT,YPT,L(LL),U,V) 
        PTOFF(1,KNB) = U 
        PTOFF(2,KNB) = V 
        DELSBA(1,KNB) = 0.0 
    2   DELSBA(2,KNB) = 0.0 
C 
C    Initialise AREA to zero, and UCURR,VCURR,DCURR to the components 
C    and squared length of the last offset.  Initial and ICASE. 
C 
        AREA = 0.0 
        UCURR = PTOFF(1,KNB) 
        VCURR = PTOFF(2,KNB) 
        DCURR = UCURR**2+VCURR**2 
        LL = LLO-1 
        ICASE = 0 
C 
C    The main loop scans through the neighbours.  It is a DO loop on 
C    KCURR, with LL scanning the original contiguity list in parallel. 
C 
        DO 3  KCURR = 1,KNB 
        LL = LL+1 
C 
C    Save old values and pick up new ones.  Set ICASE to 5 if KCURRth 
C    neighbour is virtual. 
C 
        UPREV = UCURR 
        VPREV = VCURR 
        DPREV = DCURR 
        UCURR = PTOFF(1,KCURR) 
        VCURR = PTOFF(2,KCURR) 
        DCURR = UCURR**2+VCURR**2 
        N = L(LL) 
        IF(N.LE.0) ICASE = 5 
C 
C    The subarea for KCURR is calculated by breaking the subtile into 
C    triangles and accumulating their areas in SBA, which is initialised 
C    to zero.  All these triangles have as a common vertex the vertex of 
C    the tile clockwise from KCURR.  Find the offset of this hinge point 
C    as (UHINGE,VHINGE). 
C 
        SBA = 0.0 
        C = 0.5/(UPREV*VCURR-UCURR*VPREV) 
        UHINGE = (DPREV*VCURR-DCURR*VPREV)*C 
        VHINGE = (DCURR*UPREV-DPREV*UCURR)*C 
C 
C    The first such triangle has the KCURRth face of the tile itself as 
C    a side.  One end of this side is the hinge point.  Find the other 
C    end, as (UVTX,VVTX). 
C 
        KEXTR = KCURR+1 
        IF(KEXTR.GT.KNB) KEXTR = 1 
        UEXTR = PTOFF(1,KEXTR) 
        VEXTR = PTOFF(2,KEXTR) 
        DEXTR = UEXTR**2+VEXTR**2 
        C = 0.5/(UCURR*VEXTR-UEXTR*VCURR) 
        UVTX = (DCURR*VEXTR-DEXTR*VCURR)*C 
        VVTX = (DEXTR*UCURR-DCURR*UEXTR)*C 
C 
C    If the KCURRth neighbour is real, calculate the main term in the 
C    subtile gradient - the only one if all neighbours are real. 
C 
        IF(N.LE.0) GOTO 4 
        FACTOR = 0.5*SQRT(((UVTX-UHINGE)**2+(VVTX-VHINGE)**2)/DCURR) 
        DELSBA(1,KCURR) = DELSBA(1,KCURR)+FACTOR*(UVTX+UHINGE) 
        DELSBA(2,KCURR) = DELSBA(2,KCURR)+FACTOR*(VVTX+VHINGE) 
C 
C    Now enter a test-controlled inner loop to find the remaining 
C    vertices of the subtile in anticlockwise order.  These vertices are 
C    internal to the tile.  Each new vertex, its predecessor, and the 
C    hinge point together form a triangle, and these triangles give the 
C    desired decomposition of the subtile.  We begin by saving old 
C    values. 
C 
    4   KSCAN = KEXTR 
        USCAN = UEXTR 
        VSCAN = VEXTR 
        UOLD = UVTX 
        VOLD = VVTX 
C 
C    Attempt to find the next vertex.  The candidates for defining the 
C    next vertex are all the neighbours following KSCAN and 
C    preceding KCURR.  They are scanned in an inner test-controlled loop 
C    to find which one gives the tightest intercept on the perpendicular 
C    bisector of (UCURR,VCURR) and (USCAN,VSCAN).  The winner, or the 
C    last such in the case of ties, becomes KEXTR.  If there are no 
C    such candidates, KEXTR is not advanced from KSCAN, the attempt to 
C    find a new vertex fails, and the subtile has been exhausted.  In 
C    that case set (UVTX,VVTX) to (UHINGE,VHINGE) for the benefit of 
C    the gradient calculation if KCURR is virtual. 
C 
        K = KSCAN 
    5   K = K+1 
        IF(K.GT.KNB) K = 1 
        IF(K.EQ.KCURR) GOTO 6 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        D = (U-USCAN)*(VCURR-VSCAN)-(UCURR-USCAN)*(V-VSCAN) 
        IF(D.LE.0.0) GOTO 5 
        T = ((U-USCAN)*(U-UCURR)+(V-VSCAN)*(V-VCURR))/D 
        IF(KEXTR.EQ.KSCAN) GOTO 7 
        IF(TEXTR-T) 5,8,7 
    7   TEXTR = T 
    8   KEXTR = K 
        UEXTR = U 
        VEXTR = V 
        GOTO 5 
    6   IF(KEXTR.NE.KSCAN) GOTO 9 
        UVTX = UHINGE 
        VVTX = VHINGE 
        GOTO 10 
C 
C    The new vertex has been identified.  Find its offset as (UVTX,VVTX) 
C    and hence find the area of the triangle and add it to SBA. 
C 
    9   U2 = USCAN-UCURR 
        V2 = VSCAN-VCURR 
        U3 = UEXTR-UCURR 
        V3 = VEXTR-VCURR 
        S2 = U2*(USCAN+UCURR)+V2*(VSCAN+VCURR) 
        S3 = U3*(UEXTR+UCURR)+V3*(VEXTR+VCURR) 
        C = 0.5/(U2*V3-U3*V2) 
        UVTX = (S2*V3-S3*V2)*C 
        VVTX = (S3*U2-S2*U3)*C 
        SBA = SBA+0.5*((UOLD-UHINGE)*(VVTX-VHINGE)-(UVTX-UHINGE)* 
     1    (VOLD-VHINGE)) 
C 
C    If the KCURRth neighbour is virtual then an infinitesimal 
C    displacement of the point itself causes an infinitesimal 
C    displacement of the KCURRth neighbour, related to that of the 
C    point by reflexion in the boundary.  The boundary is left fixed, 
C    and rather than calculate the usual term in the gradient for 
C    the corresponding subtile and then cancel it, we have simply 
C    omitted it.  But two systems of compensating terms are needed, 
C    which cancel one another in pairs but are allocated to different 
C    sums and so have an effect.  Those allocated to the gradient of the 
C    subtile area of the KCURRth neighbour are, like that subtile area 
C    itself, omitted from the summations involved in interpolant 
C    calculations and are thus in practice lost, but the balancing terms 
C    paired with them are allocated to the gradients of the areas of 
C    those subtiles contiguous to that of KCURR.  Where such subtiles 
C    correspond to real points, the terms have an effect in the 
C    interpolant calculations. 
C 
C    Check if the KCURRth neighbour is virtual. 
C 
   10   IF(N.GT.0) GOTO 11 
C 
C    It is.  Find the gradient contribution arising from the 
C    KCURR - KSCAN common boundary w.r.t. displacement of KCURR. 
C 
        FACTOR = 0.5*SQRT(((UVTX-UOLD)**2+(VVTX-VOLD)**2)/ 
     1    ((USCAN-UCURR)**2+(VSCAN-VCURR)**2)) 
        GX = FACTOR*(UVTX+UOLD-2.0*UCURR) 
        GY = FACTOR*(VVTX+VOLD-2.0*VCURR) 
C 
C    Reflect it in the generating point - KCURR boundary to obtain the 
C    gradient contribution w.r.t. displacement of the generating 
C    point. 
C 
        FACTOR = 2.0*(GX*UCURR+GY*VCURR)/DCURR 
        GX = GX-FACTOR*UCURR 
        GY = GY-FACTOR*VCURR 
C 
C    Add it into the gradient for KCURR to keep the books straight 
C    even though we do not use this value, and subtract it from the 
C    gradient for KSCAN, which is an effective change unless KSCAN 
C    is itself virtual. 
C 
        DELSBA(1,KCURR) = DELSBA(1,KCURR)+GX 
        DELSBA(2,KCURR) = DELSBA(2,KCURR)+GY 
        DELSBA(1,KSCAN) = DELSBA(1,KSCAN)-GX 
        DELSBA(2,KSCAN) = DELSBA(2,KSCAN)-GY 
   11   IF(KEXTR.NE.KSCAN) GOTO 4 
C 
C    Store the subarea and add it to the area. 
C 
        SBAREA(KCURR) = SBA 
    3   AREA = AREA+SBA 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYSBA 4.4 DATED 6 MARCH 1981 
C 
        SUBROUTINE TRYSBA(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, 
     2    NINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    For a trial point (XPT,YPT) calculates the contiguities this point 
C    would have if inserted into the tessellation, the offsets to its 
C    neighbours, the area of its tile, and the areas of its subtiles. 
C 
C    The user supplies the TILE4 data structure, the coordinates 
C    (XPT,YPT) of the trial point, unset variables NINDEX,AREA,KNB, 
C    ICASE, and arrays SBAREA(KNBMAX),PTOFF(2,KNBMAX).  NINDEX returns 
C    the index of the accepted point nearest to the trial point.  AREA 
C    returns the area of the tile of the trial point.  KNB returns the 
C    number of neighbours.  ICASE returns 0 if all these are points, and 
C    5 if some are constraints.  The data item containing the contiguity 
C    list is returned at LPTR in the heap L, and neighbours are 
C    referenced in the order in which they occur in this contiguity 
C    list.  SBAREA(K) for K = 1,...,KNB returns the area of the subtile 
C    for the Kth neighbour.  (PTOFF(1,K),PTOFF(2,K)) for K = 1,...,KNB 
C    return the offsets to the Kth neighbour as a real or virtual point. 
C    IFLAG is zero on successful return.  Nonzero values indicate error 
C    returns, as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap  overflow 
C       8  NSTART not the index of an accepted point 
C       10 Transfer array overflow - KNB would exceed KNBMAX 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are interpreted in terms of the tolerance 
C    values EPSCN,EPSPT.  If error return with IFLAG set to 5 occurs, 
C    NINDEX correctly returns the index of the duplicated point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Call TRYPT to find the contiguity list for the trial point. 
C 
        CALL TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Check IFLAG and return if it is nonzero. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Call LOCSBA for the contiguity list at LPTR constructed for the 
C    trial point (XPT,YPT), to work out subareas. 
C 
    1   CALL LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LPTR,XPT,YPT, 
     1    AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBA, and these values can 
C    be passed straight back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ACCSBA 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE ACCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    AREA,SBAREA,PTOFF,KNBMAX,KNB) 
C 
C    Calculates the coordinate offsets of the neighbours of NINDEX, 
C    the areas of their subtiles of the tile of NINDEX, and the area 
C    of the tile of NINDEX. 
C 
C    The user supplies arrays SBAREA(KNBMAX),PTOFF(2,KNBMAX).  KNB 
C    returns the number of neighbours of NINDEX.  If this would exceed 
C    KNBMAX, error return with IFLAG set to 10 occurs.  If NINDEX is 
C    not the index of an accepted point, error return with IFLAG set to 
C    9 occurs.  Otherwise IFLAG is zero on successful return.  AREA 
C    returns the total area of the tile of NINDEX.  (PTOFF(1,K), 
C    PTOFF(2,K)) for K = 1,...,KNB return the coordinate offsets of 
C    the neighbours of NINDEX as actual or virtual points, in the order 
C    in which they are encountered in its contiguity list.  SBAREA(K) 
C    for K = 1,...,KNB returns the area of the subtile associated with 
C    the Kth neighbour. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check that NINDEX is the index of an accepted point.  If not, 
C    return with IFLAG set to 9.  If so, LOCN is the address of its 
C    data item in the heap. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN.GT.0) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Pick up (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    2   XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Call subroutine LOCSBA to calculate the subareas. 
C 
        CALL LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XINDEX,YINDEX,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    IFLAG will be 0 or 10 on return from LOCSBA, and these values 
C    can just be passed back. 
C 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCSBA 4.3 DATED 10 MARCH 1981 
C 
        SUBROUTINE LOCSBA(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,LOCN, 
     1    XPT,YPT,AREA,SBAREA,PTOFF,KNBMAX,KNB,ICASE) 
C 
C    Calculates coordinate offsets to neighbours, tile area, and 
C    subtile areas for a point (XPT,YPT) whose data item in the heap 
C    is held at LOCN.  The data item may be either that of an 
C    accepted point, or that of a trial point prepared at location 
C    LPTR by a call to TRYPT.  The user should not normally call 
C    subroutine LOCSBA directly, but rather via a call to ACCSBA to 
C    find the subareas for an accepted point, or via a call to TRYSBA 
C    to find those for a trial point. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    SBAREA(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    No check is made on LOCN.  If it is not in fact the address of a 
C    data item, chaos will result.  Assuming this does not happen, 
C    pick up LLO as the base of the contiguity list, and LHI as 
C    the top of it.  Check that the list is not too long. 
C 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        IF(LHI-LLO+1.LE.KNBMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Find the coordinate offsets of all the neighbours. 
C 
    1   KNB = 0 
        DO 4  LL = LLO,LHI 
        KNB = KNB+1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XPT,YPT,L(LL),U,V) 
        PTOFF(1,KNB) = U 
    4   PTOFF(2,KNB) = V 
C 
C    Initialise AREA to zero, and UCURR,VCURR,DCURR to the components 
C    and squared length of the last offset. 
C 
        AREA = 0.0 
        UCURR = PTOFF(1,KNB) 
        VCURR = PTOFF(2,KNB) 
        DCURR = UCURR*UCURR+VCURR*VCURR 
        LL = LLO-1 
        ICASE = 0 
C 
C    The main loop scans through the neighbours.  Save old values and 
C    pick up new ones. 
C 
        DO 5  KCURR = 1,KNB 
        LL = LL+1 
        UPREV = UCURR 
        VPREV = VCURR 
        DPREV = DCURR 
        UCURR = PTOFF(1,KCURR) 
        VCURR = PTOFF(2,KCURR) 
        DCURR = UCURR*UCURR+VCURR*VCURR 
        IF(L(LL).LE.0) ICASE = 5 
C 
C    The subarea for KCURR is calculated by breaking the subtile into 
C    triangles and accumulating their areas in SBA, which is initialised 
C    to zero.  All these triangles have as a common vertex the vertex of 
C    the tile clockwise from KCURR.  Find the offset of this hinge point 
C    as (UHINGE,VHINGE). 
C 
        SBA = 0.0 
        C = 0.5/(UPREV*VCURR-UCURR*VPREV) 
        UHINGE = (DPREV*VCURR-DCURR*VPREV)*C 
        VHINGE = (DCURR*UPREV-DPREV*UCURR)*C 
C 
C    The first such triangle has the Kth face of the tile itself as a 
C    side.  One end of this side is the hinge point.  Find the other 
C    end, as (UVTX,VVTX). 
C 
        KEXTR = KCURR+1 
        IF(KEXTR.GT.KNB) KEXTR = 1 
        UEXTR = PTOFF(1,KEXTR) 
        VEXTR = PTOFF(2,KEXTR) 
        DEXTR = UEXTR*UEXTR+VEXTR*VEXTR 
        C = 0.5/(UCURR*VEXTR-UEXTR*VCURR) 
        UVTX = (DCURR*VEXTR-DEXTR*VCURR)*C 
        VVTX = (DEXTR*UCURR-DCURR*UEXTR)*C 
C 
C    Now enter a test-controlled inner loop to find the remaining 
C    vertices of the subtile in anticlockwise order.  These vertices 
C    are internal to the tile.  Each new vertex, its predecessor, and 
C    the hinge point together form a triangle, and these triangles give 
C    the desired decomposition of the subtile.  We begin by saving 
C    old values. 
C 
    6   KSCAN = KEXTR 
        USCAN = UEXTR 
        VSCAN = VEXTR 
        UOLD = UVTX 
        VOLD = VVTX 
C 
C    Attempt to find the next vertex.  The candidates for defining the 
C    next vertex are the neighbours following KSCAN and 
C    preceding KCURR.  They are scanned in an inner test-controlled loop 
C    to find which one gives the tightest intercept on the perpendicular 
C    bisector of (UCURR,VCURR) and (USCAN,VSCAN).  The winner, or the 
C    last such in the case of ties, becomes KEXTR.  If there are no 
C    such candidates, KEXTR is not advanced from KSCAN, the attempt to 
C    find a new vertex fails, and the subtile has been exhausted. 
C 
        K = KSCAN 
    7   K = K+1 
        IF(K.GT.KNB) K = 1 
        IF(K.EQ.KCURR) GOTO 8 
        U = PTOFF(1,K) 
        V = PTOFF(2,K) 
        D = (U-USCAN)*(VCURR-VSCAN)-(UCURR-USCAN)*(V-VSCAN) 
        IF(D.LE.0.0) GOTO 7 
        T = ((U-USCAN)*(U-UCURR)+(V-VSCAN)*(V-VCURR))/D 
        IF(KEXTR.EQ.KSCAN) GOTO 9 
        IF(TEXTR-T) 7,10,9 
    9   TEXTR = T 
   10   KEXTR = K 
        UEXTR = U 
        VEXTR = V 
        GOTO 7 
    8   IF(KEXTR.EQ.KSCAN) GOTO 11 
C 
C    The new vertex has been identified.  Find its offset as (UVTX,VVTX) 
C    and hence find the area of the triangle and add it to SBA. 
C 
        U2 = USCAN-UCURR 
        V2 = VSCAN-VCURR 
        U3 = UEXTR-UCURR 
        V3 = VEXTR-VCURR 
        S2 = U2*(USCAN+UCURR)+V2*(VSCAN+VCURR) 
        S3 = U3*(UEXTR+UCURR)+V3*(VEXTR+VCURR) 
        C = 0.5/(U2*V3-U3*V2) 
        UVTX = (S2*V3-S3*V2)*C 
        VVTX = (S3*U2-S2*U3)*C 
        SBA = SBA+0.5*((UOLD-UHINGE)*(VVTX-VHINGE)-(UVTX-UHINGE)* 
     1    (VOLD-VHINGE)) 
        GOTO 6 
C 
C    Store the subarea and add it to the area. 
C 
   11   SBAREA(KCURR) = SBA 
    5   AREA = AREA+SBA 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
