C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    GENERAL ROUTINES LAST MODIFIED 16 MAY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    FUNCTION UNI4M 4.4 DATED 22 DECEMBER 1980 
C 
        FUNCTION UNI4M(MSEED) 
C 
C    Returns a pseudorandom number from the uniform distribution 
C    on (0,1).  The value of MSEED is used to generate the number 
C    and then modified.  This allows repeated calls to the function to 
C    generate a sequence of different random numbers, all depending on 
C    the initial value of MSEED. 
C 
C    The generator used is: 
C 
C      MSEED=MOD(24298*MSEED+99991,199017) 
C 
C    broken up so as to be safe in 24 bit integers. 
C    This is the generator used on Texas Instruments TI58/59 
C    calculators. 
C 
C    The user will usually supply UNI4M via the argument PRNG to some 
C    of the tile routines. 
C 
C    This generator is a convenient and simple one which is provided 
C    so as to make programs self-contained and portable.  Its 
C    statistical properties are respectable, but by no means as good 
C    as those of more sophisticated generators, and its cycle length 
C    is relatively short.  It should NOT normally be used for 
C    statistical simulation, for which purpose a generator such as 
C    that provided in the NAG library is to be preferred.  This 
C    generator is supplied with TILE 4 so as to facilitate 
C    the provision of test and demonstration programs. 
C 
C    Note: beware the fact that this function has a side-effect. 
C    Since it is used via EXTERNAL declarations, even the most 
C    arrogant optimising compiler should be wary, but you never 
C    know. 
C 
        MSEED = MOD(23*MOD(32*MOD(33*MOD(MSEED,199017),199017),199017)+ 
     1      10*MSEED+99991,199017) 
        UNI4M = (FLOAT(MSEED)+0.5)/199017.0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SHOVEL 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SHOVEL(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,IGARB,NPRINT) 
C 
C    Prints out the current status of the TILE4 database on stream 
C    NPRINT with carriage control characters. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
        WRITE(NPRINT,9001) LTOP,LBASE,LPTR,IGARB 
        LLO = 3 
        LHI = LLO-1+L(LLO-1) 
        WRITE(NPRINT,9002) JCNS,(L(LL), LL = LLO,LHI) 
        WRITE(NPRINT,9003) 
        DO 1  JCN = 1,JCNS 
        IF(JADDR(JCN).EQ.0) GOTO 2 
        LLO = JADDR(JCN)+2 
        LHI = LLO-1+L(LLO-1) 
        WRITE(NPRINT,9004) JCN,CN(1,JCN),CN(2,JCN),CN(3,JCN), 
     1    (L(LL), LL=LLO,LHI) 
        GO TO 1 
    2   WRITE(NPRINT,9005) JCN,CN(1,JCN),CN(2,JCN),CN(3,JCN) 
    1   CONTINUE 
        WRITE(NPRINT,9006) NPTS,NPTSIN,NFREE,NSTART 
        DO 3  NPT = 1,NPTS 
        LLL = NADDR(NPT) 
        IF(LLL) 4,701,5 
  701   STOP 747 
    5   LLO = LLL+2 
        LHI = L(LLO-1)+LLO-1 
        WRITE(NPRINT,9007) NPT,PT(1,NPT),PT(2,NPT), 
     1    (L(LL), LL = LLO,LHI) 
        GOTO 3 
    4   LLL = -LLL 
        WRITE(NPRINT,9008) NPT,LLL 
    3   CONTINUE 
        RETURN 
 9001   FORMAT(1H //1H ,14HFORTRAN TILE 4//1H ,10HHEAP SIZE ,I6/ 
     1    1H ,18HBASE OF WORKSPACE ,I4/1H ,19HBASE OF FREE SPACE ,I6/ 
     2    1H ,37HNUMBER OF FORCED GARBAGE COLLECTIONS ,I4) 
 9002   FORMAT(//1H ,22HNUMBER OF CONSTRAINTS ,I4//1H , 
     1    15HBOUNDARY LIST: ,10I6,(/1H ,15X,10I6)) 
 9003   FORMAT(//1H ,27HCONSTRAINT CONTIGUITY LISTS/) 
 9004   FORMAT(1H ,I4,3H : ,3F10.4,3H : ,10I6,(/1H ,40X,10I6)) 
 9005   FORMAT(1H ,I4,3H : ,3F10.4,12H : REDUNDANT) 
 9006   FORMAT(//1H ,26HNUMBER OF POINT LOCATIONS ,I6/ 
     1    1H ,26HNUMBER OF ACCEPTED POINTS ,I6/ 
     2    1H ,26HFREE LOCATION CHAIN ENTRY ,I6/ 
     3    1H ,20HWALK STARTING POINT ,6X,I6/// 
     4    1H ,22HPOINT CONTIGUITY LISTS/) 
 9007   FORMAT(1H ,I6,3H : ,2F10.4,3H : ,10I6,(/1H ,32X,10I6)) 
 9008   FORMAT(1H ,I6,20H : FREE - POINTS TO ,I6) 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE RANPT 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE RANPT(VTX,KVTX,PROP,PRNG,MSEED,XPT,YPT) 
C 
C    The KVTX vertices of a convex polygon are held in (either) cyclic 
C    order in (VTX(1,K),VTX(2,K)) for K = 1,...,KVTX, and the cumulative 
C    proportional subareas of a triangular subdivision are held in 
C    PROP(K) for K = 3,...,KVTX.  PROP can be loaded by a call to 
C    SUBDIV.  PRNG should be a U(0,1) pseudorandom number generator 
C    returning a single-precision real value and using and updating 
C    an integer seed supplied as MSEED.  The routine RANPT returns as 
C    (XPT,YPT) the coordinates of a uniform pseudorandom point in 
C    the polygon.  Normally three calls to PRNG are made. 
C 
        EXTERNAL PRNG 
        DIMENSION VTX(2,KVTX),PROP(KVTX) 
        IF(KVTX.LT.3) STOP 750 
        DO 1  J = 1,10 
        R1 = PRNG(MSEED) 
        R2 = PRNG(MSEED) 
        R3 = PRNG(MSEED) 
        IF(R1+R2-1.0) 2,1,3 
    3   R1 = 1.0-R1 
        R2 = 1.0-R2 
    2   DO 4  K = 3,KVTX 
        IF(PROP(K).GT.R3) GOTO 5 
    4   CONTINUE 
        STOP 751 
    5   XPT = R1*VTX(1,K-1)+R2*VTX(1,K)+(1.0-R1-R2)*VTX(1,1) 
        YPT = R1*VTX(2,K-1)+R2*VTX(2,K)+(1.0-R1-R2)*VTX(2,1) 
        RETURN 
    1   CONTINUE 
        STOP 752 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SCATRC 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SCATRC(CN,JCNS,PT,NPTS,IFLAG,PRNG,MSEED, 
     1    XMIN,XMAX,YMIN,YMAX) 
C 
C    Sets up entries in CN to specify as the window the rectangle 
C 
C       (XMIN,XMAX) * (YMIN,YMAX) 
C 
C    and then uses the external pseudorandom number generator 
C    PRNG, seeded by MSEED which it, and hence SCATRC, modifies, 
C    to set up entries in PT for NPTS points in random positions. 
C    IFLAG is zero on successful return.  Nonzero values indicate 
C    error return, as follows. 
C 
C       2  JCNS less than 4 
C       3  XMIN.GE.XMAX or YMIN.GE.YMAX 
C 
        EXTERNAL PRNG 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
C 
C    Check dimensionality. 
C 
        IF(NPTS.LE.0) STOP 753 
C 
C    Call SETREC to establish the window. 
C 
        CALL SETREC(CN,JCNS,IFLAG,XMIN,XMAX,YMIN,YMAX) 
C 
C    Check IFLAG on return from SETREC. 
C 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Set up the points using PRNG. 
C 
    1   DO 2  N = 1,NPTS 
        PT(1,N) = XMIN+(XMAX-XMIN)*PRNG(MSEED) 
    2   PT(2,N) = YMIN+(YMAX-YMIN)*PRNG(MSEED) 
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
C    SUBROUTINE SETREC 4.1 DATED 5 AUGUST 1980 
C 
        SUBROUTINE SETREC(CN,JCNS,IFLAG,XMIN,XMAX,YMIN,YMAX) 
C 
C    Sets up entries in CN to specify as the window the rectangle 
C 
C       (XMIN,XMAX) * (YMIN,YMAX). 
C 
C    IFLAG is zero on successful return.  Nonzero values indicate 
C    errors, as follows. 
C 
C       2   JCNS less than 4 
C       3   XMIN.GE.XMAX or YMIN.GE.YMAX 
C 
        DIMENSION CN(3,JCNS) 
C 
C    Check dimensionality. 
C 
        IF(JCNS.GE.4) GOTO 1 
        IFLAG = 2 
        RETURN 
C 
C    Check window specification. 
C 
    1   IF(XMIN.LT.XMAX.AND.YMIN.LT.YMAX) GOTO 2 
        IFLAG = 3 
        RETURN 
C 
C    Set up the window. 
C 
    2   CN(1,1) = 1.0 
        CN(2,1) = 0.0 
        CN(3,1) = -XMAX 
        CN(1,2) = 0.0 
        CN(2,2) = -1.0 
        CN(3,2) = YMIN 
        CN(1,3) = -1.0 
        CN(2,3) = 0.0 
        CN(3,3) = XMIN 
        CN(1,4) = 0.0 
        CN(2,4) = 1.0 
        CN(3,4) = -YMAX 
        IF(JCNS.EQ.4) GOTO 3 
        DO 4  J = 5,JCNS 
        CN(1,J) = 1.0 
        CN(2,J) = 0.0 
    4   CN(3,J) = -(XMAX+1.0) 
C 
C    Set IFLAG to zero and return. 
C 
    3   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TDUMP 4.2 DATED 16 FEBRUARY 1981 
C 
        SUBROUTINE TDUMP(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,NDUMP) 
C 
C    This subroutine dumps the entire tessellation data structure 
C    to device NDUMP using unformatted write statements.  It may 
C    then be read back subsequently using subroutine TLOAD (q.v.). 
C    The user is recommended to perform a garbage collection by a 
C    call to the routine GARBAJ immediately prior to a call to 
C    TDUMP in order to minimise the quantity of output. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Dump the simple INTEGERs. 
C 
        WRITE(NDUMP) JCNS,NPTS,NFREE,NSTART,NPTSIN,LTOP,LPTR,LBASE 
C 
C    Dump the simple REALs. 
C 
        WRITE(NDUMP) EPSCN,EPSPT 
C 
C    Dump the INTEGER arrays other than L. 
C 
        WRITE(NDUMP) (JADDR(J),J = 1,JCNS),(NADDR(N),N = 1,NPTS) 
C 
C    Dump L up to LPTR-1. 
C 
        LPT1 = LPTR-1 
        WRITE(NDUMP) (L(LL),LL = 1,LPT1) 
C 
C    Dump the REAL arrays. 
C 
        WRITE(NDUMP) ((CN(I,J),I = 1,3),J = 1,JCNS),(PT(1,N),PT(2,N), 
     1    N = 1,NPTS) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TLOAD 4.3 DATED 16 FEBRUARY 1981 
C 
        SUBROUTINE TLOAD(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,NLOAD) 
C 
C    This subroutine loads the entire tessellation data structure 
C    from device NLOAD using unformatted read statements. 
C    The routine checks that the dimensions of the 
C    arrays CN, JADDR, PT, NADDR are the same 
C    as they were when the file on device NLOAD was created 
C    using routine TDUMP (q.v.) and that the array L is large enough. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Load the simple INTEGERs and check the array dimensions. 
C 
        READ(NLOAD) JCNS1,NPTS1,NFREE,NSTART,NPTSIN,LTOP1,LPTR,LBASE 
        IF(JCNS1.NE.JCNS.OR.NPTS1.NE.NPTS.OR.LPTR.GT.LTOP) STOP 754 
C 
C    Load the simple REALs. 
C 
        READ(NLOAD) EPSCN,EPSPT 
C 
C    Load the INTEGER arrays other than L. 
C 
        READ(NLOAD) (JADDR(J),J = 1,JCNS),(NADDR(N),N = 1,NPTS) 
C 
C    Load L up to LPTR-1 
C 
        LPT1 = LPTR-1 
        READ(NLOAD) (L(LL),LL = 1,LPT1) 
C 
C    Load the REAL arrays. 
C 
        READ(NLOAD) ((CN(I,J),I = 1,3),J = 1,JCNS),(PT(1,N),PT(2,N), 
     1    N = 1,NPTS) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CNLD 4.2  DATED 16 MAY 1981 
C 
        SUBROUTINE CNLD(CN,JCNS,VTX,KVTX) 
C 
C    Subroutine to compute the window constraints from the coordinates 
C    of the window vertices. 
C 
C    The window vertex coordinates should be loaded by the user into 
C    (VTX(1,K),VTX(2,K)) K = 1......KVTX in clockwise order 
C    round the window.  On return the array CN will contain the proper 
C    description of the window.  The only error detected by the routine 
C    is the setting of KVTX to a greater value than JCNS.  This causes 
C    a trap to a numbered STOP statement.  The user should take care 
C    that the vertices that he specifies do form a convex window.  If 
C    they do not this will usually be trapped when WINLD is subsequently 
C    called, but this is not inevitable, as the intersecting half 
C    planes defined by a non-convex set of vertices may still form 
C    a valid, though incorrect, window. 
C 
        DIMENSION CN(3,JCNS),VTX(2,KVTX) 
C 
C    Check for a DIMENSION error. 
C 
        IF(KVTX.GT.JCNS.OR.KVTX.LT.3) STOP 755 
C 
C    Set up the first vertex and keep track of the minimum and maximum 
C    vertex x coordinates. 
C 
        XOLD = VTX(1,1) 
        YOLD = VTX(2,1) 
        XMX = XOLD 
        XMN = XMX 
C 
C    Scan through the vertices setting up CN. 
C 
        DO 1 K = 1,KVTX 
        IF(K.NE.KVTX) GOTO 2 
        XNEW = VTX(1,1) 
        YNEW = VTX(2,1) 
        GOTO 3 
    2   K1 = K+1 
        XNEW = VTX(1,K1) 
        YNEW = VTX(2,K1) 
        IF(XNEW.GT.XMX) XMX = XNEW 
        IF(XNEW.LT.XMN) XMN = XNEW 
    3   A = YOLD-YNEW 
        B = XNEW-XOLD 
C 
C    Normalise the gradient terms. 
C 
        SCA = 1.0/SQRT(A*A+B*B) 
        A = A*SCA 
        B = B*SCA 
        CN(3,K) = -A*XOLD-B*YOLD 
        CN(1,K) = A 
        CN(2,K) = B 
        XOLD = XNEW 
    1   YOLD = YNEW 
C 
C    Fill any extra CN entries with a harmless constraint. 
C 
        IF(JCNS.EQ.KVTX) RETURN 
        XMX = XMN-2.0*XMX 
        KVT1 = KVTX+1 
        DO 4 K = KVT1,JCNS 
        CN(1,K) = 1.0 
        CN(2,K) = 0.0 
    4   CN(3,K) = XMX 
        RETURN 
        END 
C 
C 
C*********************************************************************** 
C 
