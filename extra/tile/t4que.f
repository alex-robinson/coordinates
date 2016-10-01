C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    INTERROGATION ROUTINES LAST MODIFIED 3 APRIL 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE WINVTX 4.1 DATED 4 JULY 1979 
C 
        SUBROUTINE WINVTX(CN,JCNS,L,LTOP,IFLAG,VTX,KVTMAX,KVTX, 
     1    XMIN,XMAX,YMIN,YMAX) 
C 
C    Calculates the vertices of the window in clockwise cyclic order, 
C    the first vertex being that between the last and first constraints 
C    in the boundary list.  The coordinates of the Kth vertex are 
C    returned in (VTX(1,K),VTX(2,K)).  KVTX returns the number of 
C    vertices.  If this would exceed KVTMAX, which is the second 
C    dimension of the user-supplied array VTX, error return with IFLAG 
C    set to 10 occurs.  Otherwise IFLAG is zero on return.  XMIN,XMAX, 
C    YMIN,YMAX return the smallest rectangle containing the window. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP),VTX(2,KVTMAX) 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
        IF(L(2).LE.KVTMAX) GOTO 1 
        IFLAG = 10 
        RETURN 
C 
C    Calculate the vertex coordinates by solving all the relevant 
C    simultaneous equations. 
C 
    1   LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        JCN = -L(LHI) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        KVTX = 0 
        DO 2  LL = LLO,LHI 
        AJOLD = AJ 
        BJOLD = BJ 
        CJOLD = CJ 
        JCN = -L(LL) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        KVTX = KVTX+1 
        D = AJOLD*BJ-AJ*BJOLD 
        VTX(1,KVTX) = (BJOLD*CJ-BJ*CJOLD)/D 
    2   VTX(2,KVTX) = (CJOLD*AJ-CJ*AJOLD)/D 
C 
C    Calculate the extreme coordinate values. 
C 
        XMIN = VTX(1,1) 
        XMAX = XMIN 
        YMIN = VTX(2,1) 
        YMAX = YMIN 
        DO 3  K = 2,KVTX 
        XMIN = AMIN1(XMIN,VTX(1,K)) 
        XMAX = AMAX1(XMAX,VTX(1,K)) 
        YMIN = AMIN1(YMIN,VTX(2,K)) 
    3   YMAX = AMAX1(YMAX,VTX(2,K)) 
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
C    SUBROUTINE TILVTX 4.1 DATED 4 JULY 1979 
C 
        SUBROUTINE TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTX,KVTMAX,KVTX,RELAT) 
C 
C    Calculates the coordinates (if RELAT is .FALSE.) or coordinate 
C    offsets from NINDEX (if RELAT is .TRUE.) of the vertices of 
C    the tile of NINDEX.  The first vertex is that between the last 
C    and first neighbours of NINDEX, and the order is anticlockwise 
C    cyclic.  The coordinates, or coordinate offsets, of the Kth 
C    vertex are returned in (VTX(1,K),VTX(2,K)).  KVTX returns the 
C    number of vertices.  If this would exceed KVTMAX, which is the 
C    second dimension of the user-supplied array VTX, error return 
C    with IFLAG set to 10 occurs.  If NINDEX is not the index of an 
C    accepted point, error return with IFLAG set to 9 occurs. 
C    Otherwise IFLAG is zero on return. 
C 
C2      INTEGER*2 L 
        LOGICAL RELAT 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX) 
C 
C    Check that NINDEX is the index of an accepted point, and if so 
C    pick up LLO as the beginning of its contiguity list.  Return with 
C    IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
    2   IF(L(LLO-1).LE.KVTMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    3   LHI = LLO-1+L(LLO-1) 
        XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Find the offset of LHI, and the squared length of the offset. 
C 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LHI),U,V) 
        D = U*U+V*V 
C 
C    Loop through the neighbours of NINDEX, calculating vertex offsets 
C    as we go. 
C 
        KVTX = 0 
        DO 4  LL = LLO,LHI 
        UOLD = U 
        VOLD = V 
        DOLD = D 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LL),U,V) 
        D = U*U+V*V 
        KVTX = KVTX+1 
        C = 0.5/(UOLD*V-U*VOLD) 
        VTX(1,KVTX) = (DOLD*V-D*VOLD)*C 
    4   VTX(2,KVTX) = (D*UOLD-DOLD*U)*C 
C 
C    Modify to absolute coordinates if wanted. 
C 
        IF(RELAT) GOTO 5 
        DO 6  K = 1,KVTX 
        VTX(1,K) = VTX(1,K)+XINDEX 
    6   VTX(2,K) = VTX(2,K)+YINDEX 
C 
C    Set IFLAG to zero and return. 
C 
    5   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE OFFSET 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE OFFSET(CN,JCNS,PT,NPTS,XORIG,YORIG,NPT,UPT,VPT) 
C 
C    If NPT is positive, returns as (UPT,VPT) the offset of the 
C    position of NPT from (XORIG,YORIG).  If NPT is negative, 
C    returns as (UPT,VPT) the offset from (XORIG,YORIG) to its 
C    reflexion in constraint -NPT.  A zero value of NPT is trapped 
C    as a hard error, but otherwise no check on NPT is made. 
C 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
        IF(NPT) 1,701,2 
  701   STOP 742 
    1   JCN = -NPT 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        DJ = 2.0*(AJ*XORIG+BJ*YORIG+CN(3,JCN))/(AJ*AJ+BJ*BJ) 
        UPT = -AJ*DJ 
        VPT = -BJ*DJ 
        RETURN 
    2   UPT = PT(1,NPT)-XORIG 
        VPT = PT(2,NPT)-YORIG 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INFORM 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE INFORM(VTX,KVTX,P,AREA,X,Y,XX,XY,YY) 
C 
C    Returns size and shape parameters for a convex polygon. 
C 
C    (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX hold the vertices of a 
C    convex polygon in anticlockwise cyclic order.  P returns the 
C    perimeter, AREA the area.  For a uniform probability measure on 
C    the polygon, (X,Y) return the mean, (XX,XY,YY) the variances and 
C    covariance.  If the polygon is regarded as a thin uniform 
C    lamina of unit density, then AREA is its mass, (X,Y) its centroid 
C    (centre of mass), and (AREA*XX,AREA*XY,AREA*YY) its moments and 
C    product of inertia. 
C 
        DIMENSION VTX(2,KVTX) 
        IF(KVTX.LT.3) STOP 743 
C 
C    Calculate the perimeter. 
C 
        P = 0.0 
        U = VTX(1,KVTX) 
        V = VTX(2,KVTX) 
        DO 1  K = 1,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K) 
        V = VTX(2,K) 
        UDIFF = U-UOLD 
        VDIFF = V-VOLD 
    1   P = P+SQRT(UDIFF*UDIFF+VDIFF*VDIFF) 
C 
C    Calculate zeroth, first, and second moments about the first 
C    vertex by breaking the polygon into triangles with it as a 
C    common vertex. 
C 
        AREA = 0.0 
        X = 0.0 
        Y = 0.0 
        XX = 0.0 
        XY = 0.0 
        YY = 0.0 
        X1 = VTX(1,1) 
        Y1 = VTX(2,1) 
        U = VTX(1,2)-X1 
        V = VTX(2,2)-Y1 
        DO 2  K = 3,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K)-X1 
        V = VTX(2,K)-Y1 
C 
C    Find the area of the current triangle and add it to AREA. 
C 
        A = (UOLD*V-U*VOLD)/2.0 
        AREA = AREA+A 
C 
C    Accumulate mean relative to (X1,Y1) times area. 
C 
        X = X+A*(UOLD+U)/3.0 
        Y = Y+A*(VOLD+V)/3.0 
C 
C    Accumulate second moments about (X1,Y1) times area. 
C 
        XX = XX+A*(UOLD*UOLD+UOLD*U+U*U)/6.0 
        XY = XY+A*(UOLD*VOLD+(UOLD*V+U*VOLD)/2.0+U*V)/6.0 
    2   YY = YY+A*(VOLD*VOLD+VOLD*V+V*V)/6.0 
C 
C    Normalise, centre, and return. 
C 
        X = X/AREA 
        Y = Y/AREA 
        XX = XX/AREA-X*X 
        XY = XY/AREA-X*Y 
        YY = YY/AREA-Y*Y 
        X = X+X1 
        Y = Y+Y1 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GRAZE 4.3 DATED 3 APRIL 1981 
C 
        SUBROUTINE GRAZE(VTXOFF,KVTX,R,AREA,EATEN,ANGLE) 
C 
C    (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,...,KVTX give in anticlockwise 
C    order the coordinate offsets from a reference point of the vertices 
C    of a polygonal field, not necessarily convex or containing the 
C    reference point, although usually convex in applications. 
C    A goat is tethered to the reference point by a rope of length R. 
C    AREA returns the total area of the field.  EATEN returns the area 
C    of the grazed part of the field.  ANGLE returns the angle in 
C    radians over which the rope is an effective restraint. 
C 
        DIMENSION VTXOFF(2,KVTX) 
C 
C    Only the square of R is used.  Initialise AREA,EATEN,ANGLE. 
C 
        RSQ = R*R 
        AREA = 0.0 
        EATEN = 0.0 
        ANGLE = 0.0 
C 
C    Return immediately if the field is degenerate. 
C 
        IF(KVTX.LT.3) RETURN 
C 
C    Pick up the offset of the last vertex. 
C 
        U1 = VTXOFF(1,KVTX) 
        V1 = VTXOFF(2,KVTX) 
C 
C    Loop through the vertices, saving the offset of the previous vertex 
C    and picking up the offset of the current vertex. 
C 
C    Comments refer to the usual case where the field is starlike about 
C    the reference point, but all area and angle calculations are 
C    carried out with signed values which always accumulate correctly. 
C 
        DO 1  K = 1,KVTX 
        U0 = U1 
        V0 = V1 
        U1 = VTXOFF(1,K) 
        V1 = VTXOFF(2,K) 
C 
C    Find in terms of a barycentric coordinate T the points where the 
C    line through the previous and current vertices meets the circle of 
C    radius R centre the reference point.  First find the discriminant 
C    of the relevant quadratic.  If it is nonpositive a sector is eaten 
C    within the current triangle. 
C 
        X = U0*V1-U1*V0 
        U = U1-U0 
        V = V1-V0 
        D = U*U+V*V 
        H = RSQ*D-X*X 
        IF(H.LE.0.0) GOTO 2 
C 
C    Otherwise the roots are real and distinct.  Find them. 
C 
        H = SQRT(H) 
        P = U0*U+V0*V 
        TA = (-P-H)/D 
        TB = (-P+H)/D 
C 
C    The eaten portion of the current triangle is as follows. 
C 
C      0.0.LT.1.0.LE.TA.LT.TB  sector 
C      TA.LT.TB.LE.0.0.LT.1.0  sector 
C      0.0.LT.TA.LT.TB.LT.1.0  sector plus triangle plus sector 
C      TA.LE.0.0.LT.1.0.LE.TB  triangle 
C      TA.LE.0.0.LT.TB.LT.1.0  triangle plus sector 
C      0.0.LT.TA.LT.1.0.LE.TB  sector plus triangle 
C 
C    Test economically to find which case applies. 
C 
        IF(TB.LE.0.0) GOTO 2 
        IF(TB.GE.1.0) GOTO 3 
        IF(TA.LE.0.0) GOTO 4 
C 
C    Sector plus triangle plus sector. 
C 
        UA = U0+TA*U 
        VA = V0+TA*V 
        UB = U0+TB*U 
        VB = V0+TB*V 
        EATEN = EATEN+UA*VB-UB*VA 
        ANGLE = ANGLE+ATAN2(U0*VA-UA*V0,U0*UA+V0*VA) 
     1    +ATAN2(UB*V1-U1*VB,UB*U1+VB*V1) 
        GOTO 1 
C 
C    Triangle plus sector. 
C 
    4   UB = U0+TB*U 
        VB = V0+TB*V 
        EATEN = EATEN+U0*VB-UB*V0 
        ANGLE = ANGLE+ATAN2(UB*V1-U1*VB,UB*U1+VB*V1) 
        GOTO 1 
C 
C    More tests. 
C 
    3   IF(TA.LE.0.0) GOTO 5 
        IF(TA.GE.1.0) GOTO 2 
C 
C    Sector plus triangle. 
C 
        UA = U0+TA*U 
        VA = V0+TA*V 
        EATEN = EATEN+UA*V1-U1*VA 
        ANGLE = ANGLE+ATAN2(U0*VA-UA*V0,U0*UA+V0*VA) 
        GOTO 1 
C 
C    Triangle. 
C 
    5   EATEN = EATEN+X 
        GOTO 1 
C 
C    Sector. 
C 
    2   ANGLE = ANGLE+ATAN2(X,U0*U1+V0*V1) 
C 
C    Increment AREA. 
C 
    1   AREA = AREA+X 
C 
C    On completion of the loop, combine and normalise values and return. 
C 
        AREA = 0.5*AREA 
        EATEN = 0.5*(EATEN+RSQ*ANGLE) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDISC 4.2 DATED 19 NOVEMBER 1980 
C 
        SUBROUTINE INDISC(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,VTXOFF, 
     1    KVTMAX,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    For values of R in RVAL(M) for M = 1,...,MRVAL, calculates the 
C    area of the region of the window within R of accepted points, 
C    and returns it as EVAL(M), and also calculates the total angle 
C    subtended at accepted points by the curved boundary of this 
C    region, and returns it as AVAL(M).  ARETOT returns the total area 
C    of the window, and ANGTOT returns the number of accepted points 
C    times two pi.  IFLAG is zero on successful return.  The only error 
C    return is 
C 
C      10  Overflow in VTXOFF. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL),AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and rejecting if necessary. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SUBDIV 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SUBDIV(VTX,KVTX,PROP) 
C 
C    The KVTX vertices of a convex polygon are held in (either) cyclic 
C    order in (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX.  PROP(K) returns 
C    the proportion of the total area of the polygon contained within 
C    the polygon defined by the first K vertices.  PROP(1) and 
C    PROP(2) are always zero, and PROP(KVTX) is always unity. 
C 
        DIMENSION VTX(2,KVTX),PROP(KVTX) 
        IF(KVTX.LT.3) STOP 744 
        PROP(1) = 0.0 
        PROP(2) = 0.0 
        A = 0.0 
        X1 = VTX(1,1) 
        Y1 = VTX(2,1) 
        U = VTX(1,2)-X1 
        V = VTX(2,2)-Y1 
        DO 1  K = 3,KVTX 
        UOLD = U 
        VOLD = V 
        U = VTX(1,K)-X1 
        V = VTX(2,K)-Y1 
        A = A+UOLD*V-U*VOLD 
    1   PROP(K) = A 
        DO 2  K = 3,KVTX 
    2   PROP(K) = PROP(K)/A 
        PROP(KVTX) = 1.0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDICO 4.2 DATED 19 NOVEMBER 1980 
C 
        SUBROUTINE INDICO(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NCODE, 
     1    VTXOFF,KVTMAX,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT,NC) 
C 
C    Carries out the same calculation as subroutine INDISC, but 
C    accumulates contributions only from those points whose code 
C    values as looked up in NCODE are equal to NC.  ARETOT returns 
C    the total tile area associated with such points, and ANGTOT 
C    returns the number of such points times two pi. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP),NCODE(NPTS), 
     1    VTXOFF(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL),AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and rejecting if necessary. 
C 
        DO 2  NINDEX = 1,NPTS 
        IF(NCODE(NINDEX).NE.NC) GOTO 2 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE WINBD 4.1 DATED 27 SEPTEMBER 1979 
C 
        SUBROUTINE WINBD(CN,JCNS,L,LTOP,XMIN,XMAX,YMIN,YMAX) 
C 
C    XMIN,XMAX,YMIN,YMAX return the smallest rectangle in the coordinate 
C    directions which contains the window. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP) 
C 
C    Calculate the vertex coordinates, updating extreme values. 
C 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        JCN = -L(LHI) 
        AJOLD = CN(1,JCN) 
        BJOLD = CN(2,JCN) 
        CJOLD = CN(3,JCN) 
        JCN = -L(LLO) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        D = AJOLD*BJ-AJ*BJOLD 
        XMIN = (BJOLD*CJ-BJ*CJOLD)/D 
        YMIN = (CJOLD*AJ-CJ*AJOLD)/D 
        XMAX = XMIN 
        YMAX = YMIN 
        LLO = LLO+1 
        DO 1  LL = LLO,LHI 
        AJOLD = AJ 
        BJOLD = BJ 
        CJOLD = CJ 
        JCN = -L(LL) 
        AJ = CN(1,JCN) 
        BJ = CN(2,JCN) 
        CJ = CN(3,JCN) 
        D = AJOLD*BJ-AJ*BJOLD 
        X = (BJOLD*CJ-BJ*CJOLD)/D 
        Y = (CJOLD*AJ-CJ*AJOLD)/D 
        XMIN = AMIN1(XMIN,X) 
        XMAX = AMAX1(XMAX,X) 
        YMIN = AMIN1(YMIN,Y) 
    1   YMAX = AMAX1(YMAX,Y) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILVTE 4.1 DATED 21 OCTOBER 1980 
C 
        SUBROUTINE TILVTE(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTX,KVTMAX,KVTX,RELAT,IEDGE) 
C 
C    Calculates the coordinates (if RELAT is .FALSE.) or coordinate 
C    offsets from NINDEX (if RELAT is .TRUE.) of the vertices of 
C    the tile of NINDEX.  The first vertex is that between the last 
C    and first neighbours of NINDEX, and the order is anticlockwise 
C    cyclic.  The coordinates, or coordinate offsets, of the Kth 
C    vertex are returned in (VTX(1,K),VTX(2,K)).  KVTX returns the 
C    number of vertices.  If this would exceed KVTMAX, which is the 
C    second dimension of the user-supplied array VTX, error return 
C    with IFLAG set to 10 occurs.  If NINDEX is not the index of an 
C    accepted point, error return with IFLAG set to 9 occurs. 
C    Otherwise IFLAG is zero on return.  Thus this routine does 
C    exactly the same as TILVTX, which it closely resembles.  However 
C    it returns one additional argument IEDGE.  This returns -1 if 
C    NINDEX is contiguous to some constraint, and 0 or 1 if not.  The 
C    value 1 is returned if the tile of NINDEX could be affected by 
C    points outside the current window if it were enlarged, and 
C    otherwise 0 is returned. 
C 
C2      INTEGER*2 L 
        LOGICAL RELAT 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTX(2,KVTMAX) 
C 
C    Check that NINDEX is the index of an accepted point, and if so 
C    pick up LLO as the beginning of its contiguity list.  Return with 
C    IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 1 
        LLO = NADDR(NINDEX)+2 
        IF(LLO.GT.2) GOTO 2 
    1   IFLAG = 9 
        RETURN 
C 
C    Check that KVTMAX is large enough, and return with IFLAG set to 
C    10 if not. 
C 
    2   IF(L(LLO-1).LE.KVTMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
C 
C    Pick up LHI as the top of the contiguity list for NINDEX, and 
C    (XINDEX,YINDEX) as the coordinates of NINDEX. 
C 
    3   LHI = LLO-1+L(LLO-1) 
        XINDEX = PT(1,NINDEX) 
        YINDEX = PT(2,NINDEX) 
C 
C    Initialise IEDGE to zero. 
C 
        IEDGE = 0 
C 
C    Find the offset of LHI, and the squared length of the offset. 
C 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,L(LHI),U,V) 
        D = U*U+V*V 
C 
C    Loop through the neighbours of NINDEX, calculating vertex offsets 
C    as we go, and resetting IEDGE to -1 if any neighbour is a 
C    constraint. 
C 
        KVTX = 0 
        DO 4  LL = LLO,LHI 
        UOLD = U 
        VOLD = V 
        DOLD = D 
        N = L(LL) 
        IF(N.LT.0) IEDGE = -1 
        CALL OFFSET(CN,JCNS,PT,NPTS,XINDEX,YINDEX,N,U,V) 
        D = U*U+V*V 
        KVTX = KVTX+1 
        C = 0.5/(UOLD*V-U*VOLD) 
        VTX(1,KVTX) = (DOLD*V-D*VOLD)*C 
    4   VTX(2,KVTX) = (D*UOLD-DOLD*U)*C 
C 
C    Unless IEDGE is already -1, scan vertices against constraints.  If 
C    any vertex is nearer to a constraint than to NINDEX, reset IEDGE 
C    to 1 and drop out of the scan. 
C 
        IF(IEDGE.EQ.-1) GOTO 7 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        DO 5  K = 1,KVTX 
        UVTX = VTX(1,K) 
        VVTX = VTX(2,K) 
        DVTX = SQRT(UVTX**2+VVTX**2) 
        XVTX = XINDEX+UVTX 
        YVTX = YINDEX+VVTX 
        DO 5  LL = LLO,LHI 
        J = -L(LL) 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        D = -(AJ*XVTX+BJ*YVTX+CN(3,J))/SQRT(AJ**2+BJ**2) 
        IF(D.LT.DVTX) GOTO 6 
    5   CONTINUE 
        GOTO 7 
    6   IEDGE = 1 
C 
C    Modify to absolute coordinates if wanted. 
C 
    7   IF(RELAT) GOTO 8 
        DO 9  K = 1,KVTX 
        VTX(1,K) = VTX(1,K)+XINDEX 
    9   VTX(2,K) = VTX(2,K)+YINDEX 
C 
C    Set IFLAG to zero and return. 
C 
    8   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE MEAD 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE MEAD(VTXOFF,KVTX,AREA,ECCIR,ABCEN) 
C 
C    Calculates tile area, eccircularity, and abcentricity as defined 
C    by R. Mead, "A relationship between individual plant spacing and 
C    yield", Annals of Botany (New Series) 30 (1966) pp. 301-309. 
C 
C    The statistics are calculated for a tile the coordinates offsets 
C    of whose KVTX vertices are held in anticlockwise cyclic order 
C    in (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,...,KVTX.  These values 
C    will usually have been set up by a call to TILVTX with RELAT 
C    set to .TRUE..  AREA, ECCIR, ABCEN return the area, 
C    eccircularity, and abcentricity respectively. 
C 
        DIMENSION VTXOFF(2,KVTX) 
        IF(KVTX.LT.3) STOP 745 
C 
C    Calculate the area of the tile and the offset of its centroid. 
C 
        UCURR = VTXOFF(1,KVTX) 
        VCURR = VTXOFF(2,KVTX) 
        AREA = 0.0 
        UCENTR = 0.0 
        VCENTR = 0.0 
        DO 1  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = VTXOFF(1,K) 
        VCURR = VTXOFF(2,K) 
        A = 0.5*(UOLD*VCURR-UCURR*VOLD) 
        AREA = AREA+A 
        UCENTR = UCENTR+A*(UOLD+UCURR)/3.0 
    1   VCENTR = VCENTR+A*(VOLD+VCURR)/3.0 
        UCENTR = UCENTR/AREA 
        VCENTR = VCENTR/AREA 
        OFFCEN = SQRT(UCENTR**2+VCENTR**2) 
C 
C    Calculate the weighted sum of centroid-to-vertex distances, using 
C    the exterior angles at the vertices as weights; they must in fact 
C    add up to two pi, but we obtain the normalising constant by 
C    accumulating them. 
C 
        WTSUM = 0.0 
        DISSUM = 0.0 
        UNEW = VTXOFF(1,KVTX) 
        VNEW = VTXOFF(2,KVTX) 
        UCURR = VTXOFF(1,KVTX-1) 
        VCURR = VTXOFF(2,KVTX-1) 
        DO 2  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = UNEW 
        VCURR = VNEW 
        UNEW = VTXOFF(1,K) 
        VNEW = VTXOFF(2,K) 
        DIST = SQRT((UCURR-UCENTR)**2+(VCURR-VCENTR)**2) 
        PN = UNEW-UCURR 
        QN = VNEW-VCURR 
        PO = UOLD-UCURR 
        QO = VOLD-VCURR 
        ANG = ATAN2(PN*QO-PO*QN,-PO*PN-QO*QN) 
        WTSUM = WTSUM+ANG 
    2   DISSUM = DISSUM+ANG*DIST 
        DISSUM = DISSUM/WTSUM 
C 
C    Calculate the statistics. 
C 
        ECCIR = DISSUM*SQRT(WTSUM*0.5/AREA) 
        ABCEN = OFFCEN/DISSUM 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE NLIST 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE NLIST(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    NBR,PTOFF,KNBMAX,KNB) 
C 
C    Returns in NBR(K) the Kth neighbour of NINDEX and in (PTOFF(1,K), 
C    PTOFF(2,K)) the offset to it as a real or virtual point for 
C    K = 1,...,KNB, where KNB returns the number of neighbours. 
C    IFLAG returns 0 on normal return, 9 if NINDEX is not the index 
C    of an accepted point, and 10 if KNB exceeds KNBMAX. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    NBR(KNBMAX),PTOFF(2,KNBMAX) 
C 
C    Check NINDEX. 
C 
        IF(NINDEX.LT.1.OR.NINDEX.GT.NPTS) GOTO 1 
        LOCN = NADDR(NINDEX) 
        IF(LOCN) 1,701,2 
  701   STOP 746 
    1   IFLAG = 9 
        RETURN 
C 
C    Copy the contiguity list after checking its size. 
C 
    2   KNB = L(LOCN+1) 
        IF(KNB.LE.KNBMAX) GOTO 3 
        IFLAG = 10 
        RETURN 
    3   K = 0 
        LLO = LOCN+2 
        LHI = LLO-1+L(LLO-1) 
        DO 4  LL = LLO,LHI 
        K = K+1 
    4   NBR(K) = L(LL) 
C 
C    Calculate and save the offsets. 
C 
        XORIG = PT(1,NINDEX) 
        YORIG = PT(2,NINDEX) 
        DO 5  K = 1,KNB 
        N = NBR(K) 
        CALL OFFSET(CN,JCNS,PT,NPTS,XORIG,YORIG,N,U,V) 
        PTOFF(1,K) = U 
    5   PTOFF(2,K) = V 
C 
C    Set IFLAG to 0 and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE AREAS 4.1 DATED 27 OCTOBER 1980 
C 
        SUBROUTINE AREAS(VTXOFF,KVTX,AIN,AREA,AOUT) 
C 
C    Calculates tile area together with the area of the largest disc 
C    centred at the generating point and contained within the tile and 
C    the area of the smallest disc centred at the generating point 
C    and containing the tile.  These values are returned as AREA, AIN, 
C    and AOUT respectively.  (VTXOFF(1,K),VTXOFF(2,K)) for K = 1,..., 
C    KVTX hold the vertex offsets from the generating point for the 
C    KVTX vertices of the tile in anticlockwise cyclic order.  These 
C    values will usually have been loaded by a call to TILVTX or to 
C    TILVTE. 
C 
        DIMENSION VTXOFF(2,KVTX) 
C 
C    The calculations are accomplished trivially in one scan. 
C 
        UCURR = VTXOFF(1,KVTX) 
        VCURR = VTXOFF(2,KVTX) 
        AREA = 0.0 
        AOUT = 0.0 
        AIN = UCURR**2+VCURR**2 
        DO 1  K = 1,KVTX 
        UOLD = UCURR 
        VOLD = VCURR 
        UCURR = VTXOFF(1,K) 
        VCURR = VTXOFF(2,K) 
        A2 = UOLD*VCURR-UCURR*VOLD 
        AREA = AREA+A2 
        AOUT = AMAX1(AOUT,UCURR**2+VCURR**2) 
    1   AIN = AMIN1(AIN,(A2**2)/((UCURR-UOLD)**2+(VCURR-VOLD)**2)) 
        AREA = 0.5*AREA 
        AOUT = 3.141593*AOUT 
        AIN = 3.141593*AIN 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRIM 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE TRIM(CN,JCNS,L,LTOP,EPSCN,IFLAG,XPT,YPT,VTXOFF, 
     1    KVTMAX,KVTX,ITRIM,WORK) 
C 
C    The array VTXOFF(2,KVTMAX) holds on entry the vertex offsets from 
C    (XPT,YPT) of a polygon, usually the tile of that point, in 
C    anticlockwise order.  This is modified by TRIM to give the vertex 
C    offsets of the intersection of it with a reduced window moved in 
C    by EPSCN from the original.  The original number of vertices is 
C    KVTX.  This is modified to the new number, or conventionally set 
C    to 1 if the trimmed polygon has empty interior.  IFLAG is zero on 
C    successful return, and returns 10 if at any stage in the 
C    calculation KVTX would exceed KVTMAX, in which case the information 
C    in VTXOFF is lost.  ITRIM returns -1,0,1 as (XPT,YPT) is strictly 
C    within, on, or strictly outside the reduced window.  The array 
C    WORK(2,KVTMAX) must be supplied as workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),L(LTOP),VTXOFF(2,KVTMAX),WORK(2,KVTMAX) 
C 
C    Initialise ITRIM.  Convert offsets to absolute coordinates. 
C 
        ITRIM = -1 
        DO 1  K = 1,KVTX 
        VTXOFF(1,K) = VTXOFF(1,K)+XPT 
    1   VTXOFF(2,K) = VTXOFF(2,K)+YPT 
C 
C    Chop the polygon by each effective constraint in turn, picking them 
C    up from the boundary list in the heap. 
C 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        DO 2  LL = LLO,LHI 
        J = -L(LL) 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        CJ = CN(3,J)+EPSCN*SQRT(AJ**2+BJ**2) 
C 
C    Test (XPT,YPT) and reset ITRIM if necessary. 
C 
        W = AJ*XPT+BJ*YPT+CJ 
        IF(W.EQ.0.0) ITRIM = MAX0(ITRIM,0) 
        IF(W.GT.0.0) ITRIM = 1 
C 
C    Quick scan through the vertices to avoid unnecessary work if the 
C    chop either misses or annihilates the polygon. 
C 
        IF(KVTX.LT.3) GOTO 2 
        KIN = 0 
        KOUT = 0 
        DO 3  K = 1,KVTX 
        IF(AJ*VTXOFF(1,K)+BJ*VTXOFF(2,K)+CJ) 4,3,5 
    4   KIN = KIN+1 
        GOTO 3 
    5   KOUT = KOUT+1 
    3   CONTINUE 
        IF(KOUT.EQ.0) GOTO 2 
        IF(KIN.GT.0) GOTO 6 
        KVTX = 1 
        GOTO 2 
C 
C    Properly chopped case.  Copy to WORK, modifying as we go and 
C    checking space.  Then copy back to VTXOFF. 
C 
    6   KCOPY = 0 
        X = VTXOFF(1,KVTX) 
        Y = VTXOFF(2,KVTX) 
        W = AJ*X+BJ*Y+CJ 
        DO 7  K = 1,KVTX 
        XOLD = X 
        YOLD = Y 
        WOLD = W 
        X = VTXOFF(1,K) 
        Y = VTXOFF(2,K) 
        W = AJ*X+BJ*Y+CJ 
        IF(W.GT.0.0) GOTO 8 
        IF(W.EQ.0.0.OR.WOLD.LE.0.0) GOTO 9 
        IF(KCOPY.LT.KVTMAX) GOTO 10 
        IFLAG = 10 
        RETURN 
   10   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = (WOLD*X-W*XOLD)/(WOLD-W) 
        WORK(2,KCOPY) = (WOLD*Y-W*YOLD)/(WOLD-W) 
    9   IF(KCOPY.LT.KVTMAX) GOTO 11 
        IFLAG = 10 
        RETURN 
   11   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = X 
        WORK(2,KCOPY) = Y 
        GOTO 7 
    8   IF(WOLD.GE.0.0) GOTO 7 
        IF(KCOPY.LT.KVTMAX) GOTO 12 
        IFLAG = 10 
        RETURN 
   12   KCOPY = KCOPY+1 
        WORK(1,KCOPY) = (WOLD*X-W*XOLD)/(WOLD-W) 
        WORK(2,KCOPY) = (WOLD*Y-W*YOLD)/(WOLD-W) 
    7   CONTINUE 
        KVTX = KCOPY 
        DO 13  K = 1,KVTX 
        VTXOFF(1,K) = WORK(1,K) 
   13   VTXOFF(2,K) = WORK(2,K) 
    2   CONTINUE 
C 
C    Return to relative coordinates. 
C 
        DO 14  K = 1,KVTX 
        VTXOFF(1,K) = VTXOFF(1,K)-XPT 
   14   VTXOFF(2,K) = VTXOFF(2,K)-YPT 
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
C    SUBROUTINE INDITR 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE INDITR(CN,JCNS,PT,NPTS,NADDR,L,LTOP,EPSCN,IFLAG, 
     1    VTXOFF,KVTMAX,WORK,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    Carries out the same calculation as INDISC, but trims all 
C    tiles by moving the window boundary in by amount EPSCN. 
C    ARETOT returns the area within the trimmed window, and 
C    ANGTOT returns two pi times the number of points within 
C    the trimmed window.  WORK(2,KVTMAX) must be supplied as 
C    workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),WORK(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL), 
     2    AVAL(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
    1   AVAL(M) = 0.0 
        ARETOT = 0.0 
        ANGTOT = 0.0 
C 
C    Scan points, finding the tile for each and checking IFLAG. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Trim the tile. 
C 
        CALL TRIM(CN,JCNS,L,LTOP,EPSCN,IFLAG,PT(1,NINDEX), 
     1    PT(2,NINDEX),VTXOFF,KVTMAX,KVTX,ITRIM,WORK) 
        IF(IFLAG.EQ.10) RETURN 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
    3   AVAL(M) = AVAL(M)+ANGLE 
        ARETOT = ARETOT+AREA 
        IF(ITRIM.EQ.-1) ANGTOT = ANGTOT+1.0 
    2   CONTINUE 
C 
C    Normalise ANGTOT, set IFLAG to zero, and return. 
C 
        ANGTOT = ANGTOT*6.2831853 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE INDIVT 4.1 DATED 3 APRIL 1981 
C 
        SUBROUTINE INDIVT(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG, 
     1    VTXOFF,KVTMAX,WORK,RVAL,EVAL,AVAL,MRVAL,ARETOT,ANGTOT) 
C 
C    Carries out the same calculation as INDISC, but trims all 
C    tiles by moving the window boundary in by amount RVAL(M) 
C    when making the calculations for that value.  This is done 
C    progressively, so for this routine the values in RVAL must 
C    be in ascending order.  ARETOT(M) returns the area within 
C    the trimmed window, and ANGTOT(M) returns the number of 
C    points within the trimmed window times two pi.  Note that 
C    both these quantities depend on the values in RVAL, and 
C    so in this routine ARETOT and ANGTOT have to be arrays. 
C    WORK(2,KVTMAX) must be supplied as workspace. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP), 
     1    VTXOFF(2,KVTMAX),WORK(2,KVTMAX),RVAL(MRVAL),EVAL(MRVAL), 
     2    AVAL(MRVAL),ARETOT(MRVAL),ANGTOT(MRVAL) 
C 
C    Initialise values. 
C 
        DO 1  M = 1,MRVAL 
        EVAL(M) = 0.0 
        AVAL(M) = 0.0 
        ARETOT(M) = 0.0 
    1   ANGTOT(M) = 0.0 
C 
C    Scan points, finding the tile for each and checking IFLAG. 
C 
        DO 2  NINDEX = 1,NPTS 
        CALL TILVTX(CN,JCNS,PT,NPTS,NADDR,L,LTOP,IFLAG,NINDEX, 
     1    VTXOFF,KVTMAX,KVTX,.TRUE.) 
        IF(IFLAG.EQ.10) RETURN 
        IF(IFLAG.EQ.9) GOTO 2 
C 
C    Calculate and accumulate quantities for each value of R. 
C 
        DO 3  M = 1,MRVAL 
C 
C    Trim the tile. 
C 
        CALL TRIM(CN,JCNS,L,LTOP,RVAL(M),IFLAG,PT(1,NINDEX), 
     1    PT(2,NINDEX),VTXOFF,KVTMAX,KVTX,ITRIM,WORK) 
        IF(IFLAG.EQ.10) RETURN 
C 
C    Calculate the quantities. 
C 
        CALL GRAZE(VTXOFF,KVTX,RVAL(M),AREA,EATEN,ANGLE) 
        EVAL(M) = EVAL(M)+EATEN 
        AVAL(M) = AVAL(M)+ANGLE 
        ARETOT(M) = ARETOT(M)+AREA 
    3   IF(ITRIM.EQ.-1) ANGTOT(M) = ANGTOT(M)+6.2831853 
    2   CONTINUE 
C 
C    Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
