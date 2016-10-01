C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    MANIPULATION ROUTINES LAST MODIFIED 28 JULY 1981 
C 
C    COPYRIGHT (C) 1981 UNIVERSITY OF BATH 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILE 4.6 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TILE(CN,JCNS,JADDR,PT,NPTS,NADDR,L,LTOP, 
     1    LPTR,IFLAG,IGARB) 
C 
C    A simple master routine for constructing a tessellation 
C    within an arbitrary window. 
C 
C    The subroutine is compatible with that offered as the 
C    master subroutine in TILE 003.  The calling program is 
C    responsible for entering the constraints into CN and the 
C    points into PT, and for setting the dimensions JCNS,NPTS, 
C    LTOP; JCNS serves also as the number of constraints and 
C    NPTS as the number of points.  Information is returned 
C    through the standard data structure in the heap L (which 
C    is garbaged immediately before return) accessed through 
C    the address vectors JADDR,NADDR.  LPTR is the base of the 
C    free area in L, IGARB the number of garbage collections 
C    excluding the final one.  IFLAG is zero on successful 
C    return.  A nonzero value of IFLAG indicates an error 
C    condition, as follows. 
C 
C       1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C       2  Unbounded window - constraints inadequate 
C       3  Empty window - constraints inconsistent 
C       4  No accepted points within window 
C       5  Attempt to insert duplicate point 
C       6  Heap overflow 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Hard error for nonpositive dimensions. 
C 
        IF(JCNS.LE.0.OR.NPTS.LE.0.OR.LTOP.LE.0) STOP 701 
C 
C    Set up the window and check its validity. 
C 
        CALL WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Initialise other variables. 
C 
    1   CALL CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    LPTR,LBASE,IGARB) 
C 
C    Add the points, checking the validity of each, and 
C    occasionally resetting the starting-point of the 
C    nearest-neighbour walk to keep it near to the centroid 
C    of the accepted points. 
C 
        NRESET = 2 
        DO 2  N = 1,NPTS 
        XPT = PT(1,N) 
        YPT = PT(2,N) 
        CALL ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,0.0,0.0, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
        IF(IFLAG.EQ.7) STOP 702 
        IF(IFLAG.EQ.8) STOP 703 
        IF(IFLAG.NE.5.AND.IFLAG.NE.6) GOTO 3 
        RETURN 
    3   IF(IFLAG.NE.4) GOTO 4 
        NSAVE = -NADDR(NPTS) 
        NADDR(NPTS) = -N 
        NFREE = -NADDR(N) 
        NADDR(N) = -NSAVE 
        GOTO 2 
    4   IF(NPTSIN.LT.NRESET) GOTO 2 
        NRESET = 2*NRESET 
        CALL MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP, 
     1    IFLAG,XBAR,YBAR,NINDEX) 
        IF(IFLAG.EQ.4) STOP 704 
        IF(IFLAG.EQ.8) STOP 705 
        NSTART = NINDEX 
    2   CONTINUE 
C 
C    Final garbage collection and return. 
C 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1981 
C    SUBROUTINE WINLD 4.5 DATED 28 SEPTEMBER 1981 
C 
        SUBROUTINE WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
C 
C    Constructs the description of the window from the 
C    given constraints. 
C 
C    The subroutine checks the validity of the window.  If 
C    it is valid, then a data item is constructed as follows 
C    at the bottom of the heap L: L(1) is unset; L(2) is the 
C    number of effective constraints; succeeding locations in 
C    L, as many as are needed, hold the list of effective 
C    constraints in clockwise cyclic order round the boundary 
C    of the window, each pointer back to a constraint being 
C    flagged with a negative sign.  LBASE is set to the next 
C    free location in L above this data item.  Entries in JADDR 
C    are set to positive values for effective constraints and 
C    to zero for redundant ones.  On successful exit from 
C    the subroutine, IFLAG returns the value zero.  Four 
C    possible nonzero values may be returned, each denoting 
C    a catastrophic error condition, as follows. 
C 
C       1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C       2  Unbounded window - constraints inadequate 
C       3  Empty window - constraints inconsistent 
C       6  Heap overflow 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),L(LTOP) 
        LOGICAL NEWCLK,NEWACK 
C 
C    Check for dud values of LTOP and JCNS. 
C 
        IF(LTOP.LT.5) GOTO 66 
        IF(JCNS.LT.2) GOTO 62 
C 
C    Initialise JADDR to zero and check for dud constraints. 
C 
        DO 2 J = 1,JCNS 
        JADDR(J) = 0 
        IF(CN(1,J).EQ.0.0.AND.CN(2,J).EQ.0.0) GOTO 61 
    2   CONTINUE 
C 
C    Find the first effective constraint. 
C 
        JCOUNT = 0 
        LBASE = 3 
        DO 3 J = 1,JCNS 
        A = CN(1,J) 
        B = CN(2,J) 
        C = CN(3,J) 
        NEWCLK = .FALSE. 
        NEWACK = .FALSE. 
        DO 4 J1 = 1,JCNS 
        IF(JADDR(J1).NE.0.OR.J.EQ.J1) GOTO 4 
        CALL UU(A,B,C,CN(1,J1),CN(2,J1),CN(3,J1),U,IWIPE) 
        IF(IWIPE) 5,6,7 
C 
C    Check parallel constraints. 
C 
    6   AMDA = (A*CN(1,J1)+B*CN(2,J1))/(A*A+B*B) 
        IF(U) 27,27,28 
   27   IF(AMDA) 63,61,31 
   28   IF(AMDA) 4,61,32 
C 
C    Constraint J1 is made redundant by constraint J. 
C 
   32   JADDR(J1) = -1 
        GOTO 4 
C 
C    The wipe is anti-clockwise. 
C 
    5   IF(NEWACK) GOTO 9 
        UACK = U 
        NEWACK = .TRUE. 
        GOTO 4 
    9   IF(U.GT.UACK) UACK = U 
        GOTO 4 
C 
C    The wipe is clockwise. 
C 
    7   IF(NEWCLK) GOTO 10 
        UCLK = U 
        JNEXT = J1 
        NEWCLK = .TRUE. 
        GOTO 4 
   10   IF(U.GE.UCLK) GOTO 4 
        UCLK = U 
        JNEXT = J1 
    4   CONTINUE 
C 
C    Check for unbounded window. 
C 
        IF(.NOT.(NEWCLK.AND.NEWACK)) GOTO 62 
C 
C    Check if the wipes overlap.  If they do the constraint is 
C    ineffective, as nowhere along its length does it form part 
C    of the window boundary. 
C 
        IF(UCLK.GT.UACK) GOTO 26 
C 
C    The constraint is ineffective. 
C 
   31   JADDR(J) = -1 
        JCOUNT = JCOUNT+1 
    3   CONTINUE 
C 
C    If this point is reached no constraint is effective: the window 
C    is empty. 
C 
C 
C    The first effective constraint has been found.  It is J. 
C 
   26   JFIRST = J 
C 
C    Consider the next constraint in a clockwise direction. 
C 
   11   L(LBASE) = -J 
        JADDR(J) = 1 
        LBASE = LBASE+1 
        IF(LBASE.GT.LTOP) GOTO 66 
        JCOUNT = JCOUNT+1 
C 
C    Check for cycling round and round. 
C 
        IF(JCOUNT.GE.JCNS) STOP 706 
C 
C    Compute the anticlockwise wipe of constraint J on JNEXT. 
C 
        CALL UU(CN(1,JNEXT),CN(2,JNEXT),CN(3,JNEXT),A,B,C,U,IWIPE) 
        IF(IWIPE) 13,777,777 
  777   STOP 707 
   13   UACK = U 
C 
C    Find the next clockwise constraint from JNEXT. 
C 
        A = CN(1,JNEXT) 
        B = CN(2,JNEXT) 
        C = CN(3,JNEXT) 
        NEWCLK = .FALSE. 
        DO 14 J1 = 1,JCNS 
        IF(J1.EQ.JFIRST) GOTO 12 
        IF(JADDR(J1).NE.0.OR.J1.EQ.JNEXT) GOTO 14 
   12   CALL UU(A,B,C,CN(1,J1),CN(2,J1),CN(3,J1),U,IWIPE) 
        IF(IWIPE) 15,16,17 
C 
C    Check parallel constraints. 
C 
   16   AMDA = (A*CN(1,J1)+B*CN(2,J1))/(A*A+B*B) 
        IF(U) 127,127,128 
  127   IF(AMDA) 63,61,710 
  128   IF(AMDA) 14,61,14 
C 
C    Anti-clockwise wipe - check for inconsistency. 
C 
   15   IF(J1.EQ.JFIRST) GOTO 14 
        IF(U.LE.UACK) GOTO 14 
  710   STOP 710 
C 
C    Clockwise wipe. 
C 
   17   IF(NEWCLK) GOTO 20 
        JNXT2 = J1 
        UCLK = U 
        NEWCLK = .TRUE. 
        GOTO 14 
   20   IF(U.GE.UCLK) GOTO 14 
        UCLK = U 
        JNXT2 = J1 
   14   CONTINUE 
C 
C    Check for unbounded window. 
C 
        IF(.NOT.NEWCLK) GOTO 62 
C 
C    Check if the window has been closed. 
C 
        IF(JNXT2.EQ.JFIRST) GOTO 22 
C 
C    It has not.  Go back and add JNEXT to the list. 
        J = JNEXT 
        JNEXT = JNXT2 
        GOTO 11 
C 
C    It has.  Add JNEXT as the last constraint. 
C 
   22   L(LBASE) = -JNEXT 
        JADDR(JNEXT) = 1 
        LBASE = LBASE+1 
        IF(LBASE.GT.LTOP) GOTO 66 
        JCNSIN = LBASE-3 
        L(2) = JCNSIN 
C 
C    Entries in JADDR that are still 0 must be ineffective. 
C    Find any initial ineffective constraints and set them 
C    to zero as well. 
C 
        DO 23 J = 1,JCNS 
        IF(JADDR(J).GT.0) GOTO 60 
   23   JADDR(J) = 0 
        GOTO 63 
C 
C    Successful return. 
C 
   60   IFLAG = 0 
        RETURN 
C 
C    Error returns. 
C 
   61   IFLAG = 1 
        RETURN 
   62   IFLAG = 2 
        RETURN 
   63   IFLAG = 3 
        RETURN 
   66   IFLAG = 6 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE CLEAR 4.1 DATED 17 AUGUST 1977 
C 
        SUBROUTINE CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    LPTR,LBASE,IGARB) 
C 
C    Initialises or reinitialises that part of the data 
C    structure not set up by subroutine WINLD. 
C 
C    The subroutine constructs a chain through the 
C    locations in NADDR with NFREE pointing to the start 
C    of the chain.  NSTART, NPTSIN, and IGARB are 
C    initialised to zero, and LPTR to LBASE, which has 
C    already been set up.  The effect of a call to 
C    subroutine CLEAR is to wipe the window clear of 
C    points ready for the construction of a completely 
C    new tessellation within the same window. 
C 
        DIMENSION NADDR(NPTS) 
        DO 1  N = 1,NPTS 
    1   NADDR(N) = -(N+1) 
        NFREE = 1 
        NSTART = 0 
        NPTSIN = 0 
        LPTR = LBASE 
        IGARB = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE ADDPT 4.5 DATED 13 JANUARY 1981 
C 
        SUBROUTINE ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Builds the point (XPT,YPT) into the tessellation and 
C    returns its index as NINDEX. 
C 
C    The subroutine attempts to build the new point whose 
C    coordinates are (XPT,YPT) into the tessellation.  A 
C    check is made that the new point is adequately inside 
C    the window.  Then the current data structure is 
C    interrogated to determine whether there are already 
C    any accepted points; if not, all that needs to be done 
C    is some appropriate initialisation.  If points are 
C    already accepted, the (a) nearest neighbour among them 
C    to the new point is found, and a check is made that 
C    the new point is not too close to it.  A free index is 
C    reserved, and NINDEX set to its value; the corresponding 
C    entries in PT are set to the coordinates of the new 
C    point.  Finally, the contiguity list for the new point 
C    is constructed, and those of its neighbours are 
C    modified, so as to modify the tessellation to include 
C    it.  On successful completion of this procedure, IFLAG 
C    returns the value zero.  Nonzero values indicate 
C    failure to insert the new point, as follows. 
C 
C       4  Attempt to insert point outside window 
C       5  Attempt to insert duplicate point 
C       6  Heap overflow *** catastrophic error 
C       7  Address vector already full 
C       8  NSTART not the index of an accepted point 
C 
C    "Outside" and "duplicate" are each interpreted in terms 
C    of the tolerance values EPSCN,EPSPT in the argument 
C    list.  If there are no previously accepted points, 
C    return with IFLAG set to 8 cannot of course occur; 
C    in that case NSTART is set up to the index of the 
C    newly-accepted first point, and otherwise it is left 
C    unchanged. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Call subroutine SEEPT to check that (XPT,YPT) is far 
C    enough inside the window; return if not, with NINDEX 
C    giving the index of an unsatisfied constraint. 
C 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    If any points are already accepted, jump forward 
C    to the main part of the routine. 
C 
    1   IF(NPTSIN.GT.0) GOTO 2 
C 
C    This is the first point.  Check that NFREE defines a 
C    valid location for it, set NINDEX to this value, and 
C    copy the coordinates of the new point to the corresponding 
C    entries in PT.  Reset NPTSIN to 1. 
C 
        IF(NFREE.LE.NPTS) GOTO 3 
        IFLAG = 7 
        RETURN 
    3   NINDEX = NFREE 
        NFREE = -NADDR(NFREE) 
        PT(1,NINDEX) = XPT 
        PT(2,NINDEX) = YPT 
        NPTSIN = 1 
C 
C    Check that there is enough space available in the heap 
C    for initialisation to be possible. 
C 
        IF(LPTR+6*L(2)+1.LE.LTOP) GOTO 4 
        IFLAG = 6 
        RETURN 
C 
C    LPTR is the first free location in the heap.  Build up 
C    contiguity lists from that location; first, that for 
C    the new point NINDEX, 
C 
    4   NADDR(NINDEX) = LPTR 
        L(LPTR) = NINDEX 
        L(LPTR+1) = L(2) 
        LPTR = LPTR+2 
        LHI = 2+L(2) 
        DO 5  LL = 3,LHI 
        LC = LHI+3-LL 
        L(LPTR) = L(LC) 
    5   LPTR = LPTR+1 
C 
C    next, that for the first effective constraint, 
C 
        J = -L(3) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LHI) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(4) 
        LPTR = LPTR+5 
C 
C    next, that for the last effective constraint, 
C 
        J = -L(LHI) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LHI-1) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(3) 
        LPTR = LPTR+5 
C 
C    and finally, for the remaining effective constraints; 
C    since the window is bounded there are at least three 
C    effective constraints. 
C 
        LHI = LHI-1 
        DO 6  LL = 4,LHI 
        J = -L(LL) 
        JADDR(J) = LPTR 
        L(LPTR) = -J 
        L(LPTR+1) = 3 
        L(LPTR+2) = L(LL-1) 
        L(LPTR+3) = NINDEX 
        L(LPTR+4) = L(LL+1) 
    6   LPTR = LPTR+5 
C 
C    Set NSTART to NINDEX, IFLAG to zero, and return. 
C 
        NSTART = NINDEX 
        IFLAG = 0 
        RETURN 
C 
C    The main part of the routine deals with the insertion 
C    of a new point when there are already some accepted 
C    points.  Call subroutine LOCPT to find the (a) nearest 
C    neighbour of the new point (XPT,YPT).  If it is too 
C    close, return with the index of the neighbour in NINDEX. 
C 
    2   CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 7 
        RETURN 
C 
C    Check that NFREE defines a valid location for the new 
C    point, set NINDEX to this value having saved the index 
C    of the nearest neighbour as NEAR, and copy the coordinates 
C    of the new point to the corresponding entries in PT. 
C    Increment NPTSIN. 
C 
    7   IF(NFREE.LE.NPTS) GOTO 8 
        IFLAG = 7 
        RETURN 
    8   NEAR = NINDEX 
        NINDEX = NFREE 
        NFREE = -NADDR(NFREE) 
        PT(1,NINDEX) = XPT 
        PT(2,NINDEX) = YPT 
        NPTSIN = NPTSIN+1 
C 
C    Initialise NNEW, the first item to be placed on the 
C    contiguity list for NINDEX; it is the neighbour of 
C    NINDEX clockwise from NEAR.  Flag degeneracy if it occurs. 
C 
        LLO = NADDR(NEAR)+2 
        LHI = LLO-1+L(LLO-1) 
        XNEAR = PT(1,NEAR) 
        YNEAR = PT(2,NEAR) 
        NNEW = 0 
        DO 9  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XNEAR,YNEAR,N,T,IWIPE) 
        IF(IWIPE.LE.0) GOTO 9 
        IF(NNEW.EQ.0) GOTO 10 
        IF(TNEW-T) 9,51,10 
   51   IDGNEW = 1 
        IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 9 
   10   IDGNEW = 0 
        NNEW = N 
        TNEW = T 
    9   CONTINUE 
C 
C    Initialise NCURR to NEAR and LEFF to LTOP.  The 
C    contiguity list for NINDEX is constructed in a 
C    temporary position working down from LTOP. 
C 
        NCURR = NEAR 
        LEFF = LTOP 
C 
C    Enter the main loop.  First update NOLD and NCURR, 
C    and pass back the degeneracy flag. 
C 
   11   NOLD = NCURR 
        NCURR = NNEW 
        IDGCUR = IDGNEW 
C 
C    Is NCURR a point or a constraint? 
C 
        IF(NCURR) 12,701,13 
  701   STOP 712 
C 
C    Constraint.  (XCURR,YCURR) is the reflexion of NINDEX 
C    in it. 
C 
   12   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XPT+BJ*YPT+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XPT-AJ*DJ 
        YCURR = YPT-BJ*DJ 
C 
C    Pick up beginning and size of contiguity list, and 
C    garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = JADDR(JCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        IFLAG = 6 
        RETURN 
C 
C    Point, with coordinates (XCURR,YCURR). 
C 
   13   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
C 
C    Pick up beginning and size of contiguity list, and 
C    garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = NADDR(NCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 14 
        IFLAG = 6 
        RETURN 
C 
C    Find the end of the old contiguity list for NCURR, and flag 
C    this list as defunct.  Save the value of LPTR as the address 
C    of the new contiguity list for NCURR, flag the new 
C    list as active, and place NINDEX on it. 
C 
   14   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LCURR = LPTR 
        L(LPTR) = NCURR 
        L(LPTR+2) = NINDEX 
        LPTR = LPTR+3 
C 
C    Contiguity is a symmetric relation: find NOLD in 
C    the old contiguity list for NCURR. 
C 
        DO 15  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 16 
   15   CONTINUE 
        STOP 713 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for NINDEX after NCURR; it is the neighbour of 
C    NINDEX clockwise from NCURR.  Flag degeneracy if it 
C    occurs, and backspace LOLD if degeneracy occurred 
C    last time. 
C 
   16   NNEW = 0 
        IDGNEW = 0 
        LL = LOLD 
        IF(IDGCUR.EQ.0) GOTO 17 
        LOLD = LOLD-1 
        IF(LOLD.LT.LLO) LOLD = LHI 
   17   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.GT.0) GOTO 18 
        IF(NNEW) 19,17,19 
   18   IF(NNEW.EQ.0) GOTO 20 
        IF(TNEW-T) 19,52,20 
   52   IDGNEW = 1 
        GOTO 21 
   20   NNEW = N 
        TNEW = T 
        GOTO 17 
C 
C    Complete the new contiguity list for NCURR, and set the 
C    appropriate pointer to it from either JADDR or NADDR. 
C 
   19   L(LPTR) = NNEW 
        LPTR = LPTR+1 
   21   L(LPTR) = N 
        LPTR = LPTR+1 
        IF(LL.EQ.LOLD) GOTO 22 
        LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        GOTO 21 
   22   L(LCURR+1) = LPTR-LCURR-2 
        IF(NCURR) 23,703,24 
  703   STOP 714 
   23   JADDR(JCURR) = LCURR 
        GOTO 25 
   24   NADDR(NCURR) = LCURR 
C 
C    Add NCURR to the contiguity list for NINDEX. 
C 
   25   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NEAR to the contiguity list for 
C    NINDEX, then the main loop is completed. 
C 
        IF(NCURR.NE.NEAR) GOTO 11 
C 
C    Move the contiguity list for NINDEX down to the base of 
C    the free area, and set up pointers and size. 
C 
        NADDR(NINDEX) = LPTR 
        L(LPTR) = NINDEX 
        L(LPTR+1) = LTOP-LEFF 
        LPTR = LPTR+2 
        LEFF = LEFF+1 
        DO 26  LL = LEFF,LTOP 
        L(LPTR) = L(LL) 
   26   LPTR = LPTR+1 
C 
C    The new point has been built into the tessellation 
C    successfully.  Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SEEPT 4.2 DATED 21 OCTOBER 1980 
C 
        SUBROUTINE SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    JINDEX,TSTVAL) 
C 
C    Tests whether (XPT,YPT) is inside the window. 
C 
C    The subroutine evaluates the perpendicular distance 
C    from (XPT,YPT) to each effective constraint, flagged 
C    with a negative sign on the inwards side.  TSTVAL returns 
C    the largest such value, JINDEX the constraint where it 
C    occurs - arbitrary choice if equality.  IFLAG returns 
C    the value zero unless any such value equals or exceeds 
C    zero, in which case it returns the value 4.  If EPSCN 
C    is positive, then -EPSCN is used as the cutoff 
C    value instead of zero. 
C 
        DIMENSION CN(3,JCNS),JADDR(JCNS) 
C 
C    Examine the constraints to find TSTVAL and JINDEX. 
C 
        JINDEX = 0 
        DO 1  J = 1,JCNS 
        IF(JADDR(J).LE.0) GOTO 1 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        T = (AJ*XPT+BJ*YPT+CN(3,J))/SQRT(AJ*AJ+BJ*BJ) 
        IF(JINDEX.EQ.0) GOTO 2 
        IF(T.LE.TSTVAL) GOTO 1 
    2   TSTVAL = T 
        JINDEX = -J 
    1   CONTINUE 
C 
C    Find the cutoff value. 
C 
        CUTOFF = 0.0 
        IF(EPSCN.GT.0.0) CUTOFF = -EPSCN 
C 
C    Report if (XPT,YPT) does not lie in the window. 
C 
        IF(TSTVAL.LT.CUTOFF) GOTO 3 
        IFLAG = 4 
        RETURN 
    3   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE LOCPT 4.2 DATED 18 JULY 1979 
C 
        SUBROUTINE LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
C 
C    Finds the nearest neighbour of (XPT,YPT). 
C 
C    The subroutine finds the index of the accepted point 
C    closest to the trial point (XPT,YPT), which does 
C    not need to be inside the window.  NINDEX returns the 
C    index of the closest point, and TSTVAL its squared 
C    distance from the trial point - arbitrary choice if 
C    equality.  IFLAG returns the value zero unless TSTVAL 
C    is zero, in which case it returns the value 5.  If 
C    EPSPT is positive then 4.0*EPSPT*EPSPT is used as the 
C    cutoff value instead of zero.  If the starting-point 
C    for the walk by which the nearest neighbour is found 
C    is not the index of an accepted point, then IFLAG 
C    returns the value 8 and no other useful information 
C    is provided. 
C 
C2      INTEGER*2 L 
        DIMENSION PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Check that NSTART is the index of an accepted point. 
C 
        IF(NSTART.LE.0.OR.NSTART.GT.NPTS) GOTO 5 
        IF(NADDR(NSTART).GT.0) GOTO 1 
    5   IFLAG = 8 
        RETURN 
C 
C    Initialise NINDEX to NSTART and then replace successively 
C    by any nearer accepted point until the (a) nearest is 
C    found. 
C 
    1   NINDEX = NSTART 
        XDIFF = XPT-PT(1,NINDEX) 
        YDIFF = YPT-PT(2,NINDEX) 
        TSTVAL = XDIFF*XDIFF+YDIFF*YDIFF 
    2   LLO = NADDR(NINDEX)+2 
        LHI = LLO-1+L(LLO-1) 
        DO 3  LL = LLO,LHI 
        NEXT = L(LL) 
        IF(NEXT.LE.0) GOTO 3 
        XDIFF = XPT-PT(1,NEXT) 
        YDIFF = YPT-PT(2,NEXT) 
        T = XDIFF*XDIFF+YDIFF*YDIFF 
        IF(T.GE.TSTVAL) GOTO 3 
        NINDEX = NEXT 
        TSTVAL = T 
        GOTO 2 
    3   CONTINUE 
C 
C    Find the cutoff value. 
C 
        CUTOFF = 0.0 
        IF(EPSPT.GT.0.0) CUTOFF = 4.0*EPSPT*EPSPT 
C 
C    Test the trial point for closeness to its nearest 
C    neighbour. 
C 
        IF(TSTVAL.GT.CUTOFF) GOTO 4 
        IFLAG = 5 
        RETURN 
    4   IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE MIDPT 4.2 DATED 3 JULY 1979 
C 
        SUBROUTINE MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP, 
     1    IFLAG,XBAR,YBAR,NINDEX) 
C 
C    Finds the centroid of the accepted points and the index 
C    of the nearest accepted point to it. 
C 
C    The subroutine returns the coordinates of the centroid 
C    of the accepted points as (XBAR,YBAR), and the index of 
C    the (a) nearest accepted point to it as NINDEX - arbitrary 
C    choice if equality.  IFLAG returns the value zero on 
C    successful completion of the calculation.  If NSTART is 
C    not the index of an accepted point then (XBAR,YBAR) is 
C    returned correctly, but IFLAG returns the value 8 and 
C    no useful information is returned in NINDEX.  If there 
C    are no accepted points then no useful information is 
C    returned; IFLAG returns the value 4. 
C 
C2      INTEGER*2 L 
        DIMENSION PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Check that there are some accepted points. 
C 
        IF(NPTSIN.GT.0) GOTO 1 
        IFLAG = 4 
        RETURN 
C 
C    Find the centroid of the accepted points. 
C 
    1   XBAR = 0.0 
        YBAR = 0.0 
        DO 2  N = 1,NPTS 
        IF(NADDR(N).LE.0) GOTO 2 
        XBAR = XBAR+PT(1,N) 
        YBAR = YBAR+PT(2,N) 
    2   CONTINUE 
        XBAR = XBAR/FLOAT(NPTSIN) 
        YBAR = YBAR/FLOAT(NPTSIN) 
C 
C    Call subroutine LOCPT to find the index of the nearest 
C    accepted point. 
C 
        CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,0.0, 
     1    IFLAG,XBAR,YBAR,NINDEX,TSTVAL) 
C 
C    Suppress the possible value 5 for IFLAG returned 
C    by this subroutine - it does not matter if the 
C    centroid coincides with an accepted point. 
C 
        IF(IFLAG.EQ.5) IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TT 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TT(CN,JCNS,PT,NPTS,X0,Y0,X1,Y1,N,T,IWIPE) 
C 
C    Calculates intercepts on the perpendicular bisector of 
C    the line joining (X0,Y0) and (X1,Y1). 
C 
        DIMENSION CN(3,JCNS),PT(2,NPTS) 
        IWIPE = 1 
        IF(N) 2,701,1 
  701   STOP 715 
    1   XN = PT(1,N) 
        YN = PT(2,N) 
        D = (XN-X0)*(Y1-Y0)-(YN-Y0)*(X1-X0) 
        T = (XN-X0)*(XN-X1)+(YN-Y0)*(YN-Y1) 
    6   IF(D) 3,4,5 
    3   IWIPE = -1 
    5   T = T/D 
        RETURN 
    4   IWIPE = 0 
        RETURN 
    2   J = -N 
        AJ = CN(1,J) 
        BJ = CN(2,J) 
        D = AJ*(Y1-Y0)-BJ*(X1-X0) 
        T = -AJ*(X0+X1)-BJ*(Y0+Y1)-2.0*CN(3,J) 
        GOTO 6 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE UU 4.1 DATED 22 JULY 1976 
C 
        SUBROUTINE UU(A0,B0,C0,A1,B1,C1,U,IWIPE) 
C 
C    Calculates intercepts on a constraint. 
C 
        IWIPE = 1 
        D = B0*A1-A0*B1 
        U = -C1+C0*(A0*A1+B0*B1)/(A0*A0+B0*B0) 
        IF(D) 1,2,3 
    1   IWIPE = -1 
    3   U = U/D 
        RETURN 
    2   IWIPE = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE GARBAJ 4.1 DATED 17 AUGUST 1977 
C 
        SUBROUTINE GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR, 
     1    LPROT) 
C 
C    Does a garbage collection on the heap.  Free space is 
C    returned to the top end, and LPTR, which points to the 
C    base of the free area, is reset.  The area below LPROT 
C    is protected. 
C 
C2      INTEGER*2 L 
        DIMENSION JADDR(JCNS),NADDR(NPTS),L(LTOP) 
C 
C    Return if there is nothing to garbage. 
C 
        IF(LPTR.EQ.LPROT) RETURN 
C 
C    The old value of LPTR is saved so that completion of the 
C    scan can be detected.  LPTR is reset to LPROT, and LSCAN 
C    is initialised to the same value. 
C 
        LSAVE = LPTR 
        LPTR = LPROT 
        LSCAN = LPROT 
C 
C    On entry to the scanning loop, LSCAN points to the first 
C    cell of the current data item.  If the entry in this cell 
C    is zero, the data item is defunct, and is skipped. 
C    Otherwise the data item is active, being the contiguity 
C    list either of an effective constraint (negative entry) 
C    or of an accepted point (positive entry).  In either of 
C    these latter cases the appropriate pointer into the 
C    heap is updated. 
C 
    1   N = L(LSCAN) 
        LLO = LSCAN 
        LSCAN = LSCAN+L(LSCAN+1)+2 
        IF(N) 2,3,4 
    2   J = -N 
        JADDR(J) = LPTR 
        GOTO 5 
    4   NADDR(N) = LPTR 
C 
C    Copy down the data item. 
C 
    5   LHI = LSCAN-1 
        DO 6  LL = LLO,LHI 
        L(LPTR) = L(LL) 
    6   LPTR = LPTR+1 
C 
C    Loop if there is more to do. 
C 
    3   IF(LSCAN.LT.LSAVE) GOTO 1 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE SUBPT 4.5 DATED 13 JANUARY 1981 
C 
        SUBROUTINE SUBPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,IFLAG,IGARB, 
     2    XPT,YPT,NINDEX) 
C 
C    Removes from the tessellation the point whose index 
C    is NINDEX and returns its coordinates as (XPT,YPT). 
C 
C    The subroutine attempts to remove from the tessellation 
C    the point whose index is NINDEX.  A check is made that 
C    this index refers to a currently accepted point.  Then 
C    the current data structure is interrogated to determine 
C    whether it is the only one; if so, all that is needed 
C    to delete it is a simple reinitialisation procedure 
C    similar to that carried out by subroutine CLEAR, but 
C    preserving the value of IGARB and the pointer chain in 
C    NADDR.  If not, a substantial calculation is needed. 
C    What is done is to carry out a scratchpad calculation 
C    of the tessellation of the whole window determined by 
C    those points contiguous to NINDEX.  This calculation 
C    is done by logic identical to that employed to 
C    construct a tessellation by repeated calls to 
C    subroutine ADDPT; the only reason why ADDPT cannot be 
C    used to do it is the need to administer a slightly more 
C    complicated data structure to maintain the scratchpad 
C    calculation distinct from the main data structure. 
C    The objects whose data items include scratchpad 
C    sections at the completion of the scratchpad 
C    calculation are (a) all effective constraints, and (b) 
C    all points contiguous to the point being removed. 
C    The form of such an extended data item is as follows: 
C    (1) back-reference; (2) total length of following part 
C    of data item; (3) main contiguity list; (4) marker zero; 
C    (5) scratchpad contiguity list.  Because the back- 
C    reference and length are maintained as in the normal 
C    form of the data base, garbage collection is not 
C    affected.  Once the scratchpad calculation is 
C    complete, a merging operation is carried out whereby 
C    the main and scratchpad lists for objects contiguous 
C    to the point being deleted are merged to form new 
C    lists describing the pattern of contiguities after 
C    removal of the point being deleted.  Then the 
C    scratchpad lists for effective constraints not 
C    contiguous to the point being deleted are killed; the 
C    data item for the deleted point is killed; NINDEX 
C    is chained into the beginning of the free space in NADDR; 
C    and finally the count of accepted points is decremented by 
C    one.  On successful completion of this procedure, 
C    return takes place with IFLAG set to zero; nonzero 
C    values indicate failure to remove the specified point, 
C    as follows. 
C 
C       6  Heap overflow *** catastrophic error 
C       9  NINDEX not the index of an accepted point 
C 
C    The coordinates of the removed point are returned as	 
C    (XPT,YPT).  If the removal leaves no accepted points, 
C    NSTART is reset to zero.  If the point being removed 
C    happens to coincide with NSTART, then NSTART is reset 
C    to a neighbouring point.  Otherwise it is not checked. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Check that NINDEX is the index of an accepted point, 
C    and return with IFLAG set to 9 if not. 
C 
        IF(NINDEX.LE.0.OR.NINDEX.GT.NPTS) GOTO 100 
        IF(NADDR(NINDEX).GT.0) GOTO 1 
  100   IFLAG = 9 
        RETURN 
C 
C    Copy the coordinates of NINDEX to (XPT,YPT); these 
C    values are not needed again. 
C 
    1   XPT = PT(1,NINDEX) 
        YPT = PT(2,NINDEX) 
C 
C    If NINDEX is the only accepted point, all that is 
C    needed to remove it is to reinitialise NSTART, NPTSIN, 
C    and LPTR, and to chain NINDEX back into NADDR. 
C 
        IF(NPTSIN.GT.1) GOTO 2 
        NSTART = 0 
        NADDR(NINDEX) = -NFREE 
        NFREE = NINDEX 
        NPTSIN = 0 
        LPTR = LBASE 
        RETURN 
C 
C    If NINDEX is not the only accepted point, a lengthy 
C    calculation is needed to remove it.  First pick up 
C    the length of its contiguity list.  Also initialise 
C    to zero the number of points accepted in the 
C    scratchpad calculation. 
C 
    2   LLOX = NADDR(NINDEX)+2 
        LSIZEX = L(LLOX-1) 
        NINX = 0 
C 
C    The outer loop of the scratchpad calculation is a 
C    scan through the objects contiguous to NINDEX. 
C 
        DO 3  LKX = 1,LSIZEX 
        LLX = NADDR(NINDEX)+1+LKX 
        NMOD = L(LLX) 
C 
C    Nothing to do unless NMOD is a point. 
C 
        IF(NMOD) 3,701,4 
  701   STOP 716 
C 
C    If any points are already accepted in the scratchpad 
C    calculation, jump forward to the main part of the 
C    loop. 
C 
    4   IF(NINX.GT.0) GOTO 5 
C 
C    This section initialises the scratchpad on detection 
C    of the first point contiguous to NINDEX.  Pick up 
C    beginning and size of its contiguity list, and garbage 
C    if necessary, checking that enough space can be recovered 
C    for its extended list. 
C 
        LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+L(2)+2.LE.LTOP) GOTO 6 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+L(2)+2.LE.LTOP) GOTO 6 
        IFLAG = 6 
        RETURN 
C 
C    Find top of main list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and kill 
C    old list.  Copy old main list across, and add marker. 
C 
    6   LHI = LLO-1+LSIZE 
        L(LPTR) = NMOD 
        L(LPTR+1) = LSIZE+1+L(2) 
        NADDR(NMOD) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 7  LL = LLO,LHI 
        L(LPTR) = L(LL) 
    7   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    The scratchpad part of the extended list is the 
C    boundary list in anticlockwise order. 
C 
        LHIC = 2+L(2) 
        DO 8  LL = 3,LHIC 
        LC = LHIC+3-LL 
        L(LPTR) = L(LC) 
    8   LPTR = LPTR+1 
C 
C    Now form the extended lists for the effective 
C    constraints; the scratchpad list in each case consists 
C    of the two contiguous effective constraints separated 
C    by the point NMOD. 
C 
        DO 9  LC = 3,LHIC 
        J = -L(LC) 
C 
C    Pick up beginning and size of main list, and garbage if 
C    necessary, checking that enough space can be recovered. 
C 
        LLO = JADDR(J)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LTOP) GOTO 10 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(J)+2 
        IF(LPTR+LSIZE+5.LE.LTOP) GOTO 10 
        IFLAG = 6 
        RETURN 
C 
C    Find top of main list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and kill 
C    old list.  Copy old main list across, and add marker. 
C 
   10   LHI = LLO-1+LSIZE 
        L(LPTR) = -J 
        L(LPTR+1) = LSIZE+4 
        JADDR(J) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 11  LL = LLO,LHI 
        L(LPTR) = L(LL) 
   11   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    Add scratchpad part of extended list, dealing 
C    properly with cases LC.EQ.3 and LC.EQ.LHIC. 
C 
        LC1 = LC-1 
        IF(LC1.LT.3) LC1 = LHIC 
        L(LPTR) = L(LC1) 
        L(LPTR+1) = NMOD 
        LC1 = LC+1 
        IF(LC1.GT.LHIC) LC1 = 3 
        L(LPTR+2) = L(LC1) 
    9   LPTR = LPTR+3 
C 
C    Initialise NEAR to NMOD, reset NINX to 1, and if 
C    necessary reset NSTART; the scratchpad calculation for 
C    the first point contiguous to NINDEX is then complete. 
C 
        NEAR = NMOD 
        NINX = 1 
        IF(NSTART.EQ.NINDEX) NSTART = NMOD 
        GOTO 3 
C 
C    Points subsequent to the first are added into the 
C    scratchpad data-base as in subroutine ADDPT, but apart 
C    from the need to copy items of the main data-base at each 
C    step, the calculations are a little simpler.  This is 
C    because various acceptability checks have already been 
C    made and because there is no need for sophistication 
C    over locating the nearest neighbour.  First extract the 
C    coordinates of the new point. 
C 
    5   XMOD = PT(1,NMOD) 
        YMOD = PT(2,NMOD) 
C 
C    Find its nearest neighbour among the points already in 
C    the scratchpad calculation. 
C 
        XDIFF = XMOD-PT(1,NEAR) 
        YDIFF = YMOD-PT(2,NEAR) 
        D = XDIFF*XDIFF+YDIFF*YDIFF 
        LLOX = NADDR(NINDEX)+2 
        LKX1 = LKX-1 
        DO 12  LKXX = 1,LKX1 
        LLXX = LLOX-1+LKXX 
        NTRY = L(LLXX) 
        IF(NTRY) 12,702,13 
  702   STOP 717 
   13   XDIFF = XMOD-PT(1,NTRY) 
        YDIFF = YMOD-PT(2,NTRY) 
        DTRY = XDIFF*XDIFF+YDIFF*YDIFF 
        IF(DTRY.GE.D) GOTO 12 
        D = DTRY 
        NEAR = NTRY 
   12   CONTINUE 
C 
C    Initialise NNEW, the first item to be placed on the 
C    scratchpad contiguity list for NMOD; it is the neighbour 
C    of NMOD clockwise from NEAR.  Flag degeneracy if it 
C    occurs. 
C 
C    First pick up the beginning and end of the extended 
C    contiguity list for NEAR, and its coordinates. 
C 
        LLO = NADDR(NEAR)+2 
        LHI = LLO-1+L(LLO-1) 
        XNEAR = PT(1,NEAR) 
        YNEAR = PT(2,NEAR) 
C 
C    Initialise NNEW to zero, and locate the marker to give 
C    the beginning of the scratchpad contiguity list. 
C 
        NNEW = 0 
        DO 14  LL = LLO,LHI 
        N = L(LL) 
        IF(N)  14,15,14 
   14   CONTINUE 
        STOP 720 
   15   LLO = LL+1 
C 
C    Now find NNEW as in subroutine ADDPT. 
C 
        DO 16  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XMOD,YMOD,XNEAR,YNEAR,N,T, 
     1    IWIPE) 
        IF(IWIPE.LE.0) GOTO 16 
        IF(NNEW.EQ.0) GOTO 17 
        IF(TNEW-T) 16,18,17 
   18   IDGNEW = 1 
        IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 16 
   17   IDGNEW = 0 
        NNEW = N 
        TNEW = T 
   16   CONTINUE 
C 
C    Initialise NCURR to NEAR and LEFF to LTOP.  The 
C    scratchpad contiguity list for NMOD is constructed 
C    in a temporary position working down from LTOP. 
C 
        NCURR = NEAR 
        LEFF = LTOP 
C 
C    Enter the main inner loop of the scratchpad 
C    calculation, carried out to incorporate NMOD into 
C    the scratchpad tessellation.  First update NOLD and 
C    NCURR, and pass back the degeneracy flag. 
C 
   19   NOLD = NCURR 
        NCURR = NNEW 
        IDGCUR = IDGNEW 
C 
C    Is NCURR a constraint or a point? 
C 
        IF(NCURR) 20,704,21 
  704   STOP 721 
C 
C    Constraint.  (XCURR,YCURR) is the reflexion of NMOD 
C    in it. 
C 
   20   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XMOD+BJ*YMOD+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XMOD-AJ*DJ 
        YCURR = YMOD-BJ*DJ 
C 
C    Pick up beginning and size of extended contiguity list,	 
C    and garbage if necessary, checking that enough space can 
C    be recovered. 
C 
        LLO = JADDR(JCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        IFLAG = 6 
        RETURN 
C 
C    Point, with coordinates (XCURR,YCURR). 
C 
   21   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
C 
C    Pick up beginning and size of extended contiguity list, 
C    and garbage if necessary, checking that enough space 
C    can be recovered. 
C 
        LLO = NADDR(NCURR)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NCURR)+2 
        IF(LPTR+LSIZE+5.LE.LEFF) GOTO 22 
        IFLAG = 6 
        RETURN 
C 
C    Find the end of the old contiguity list for NCURR, and 
C    kill that list.  Save the value of LPTR as the address 
C    of the new contiguity list for NCURR.  Flag the new 
C    list as active.  Copy the old main list across, and the 
C    marker, and start the new scratchpad list with NMOD. 
C 
   22   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LCURR = LPTR 
        L(LPTR) = NCURR 
        LPTR = LPTR+2 
        DO 23  LL = LLO,LHI 
        N = L(LL) 
        L(LPTR) = N 
        LPTR = LPTR+1 
        IF(N.EQ.0) GOTO 24 
   23   CONTINUE 
        STOP 722 
   24   LLO = LL+1 
        L(LPTR) = NMOD 
        LPTR = LPTR+1 
C 
C    Contiguity is a symmetric relation: find NOLD in the 
C    scratchpad part of the old extended contiguity list 
C    for NCURR.  Resetting LLO has already located the start 
C    of the scratchpad. 
C 
        DO 25  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 26 
   25   CONTINUE 
        STOP 723 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for NMOD after NCURR; it is the neighbour of 
C    NMOD clockwise from NCURR.  Flag degeneracy if it occurs, 
C    and backspace LOLD if degeneracy occurred last time. 
C 
   26   NNEW = 0 
        IDGNEW = 0 
        LL = LOLD 
        IF(IDGCUR.EQ.0) GOTO 27 
        LOLD = LOLD-1 
        IF(LOLD.LT.LLO) LOLD = LHI 
   27   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XMOD,YMOD,XCURR,YCURR,N,T, 
     1    IWIPE) 
        IF(IWIPE.GT.0) GOTO 28 
        IF(NNEW) 29,27,29 
   28   IF(NNEW.EQ.0) GOTO 30 
        IF(TNEW-T) 29,31,30 
   31   IDGNEW = 1 
        GOTO 32 
   30   NNEW = N 
        TNEW = T 
        GOTO 27 
C 
C    Complete the scratchpad part of the new contiguity list 
C    for NCURR, and set the appropriate pointer to it from 
C    either JADDR or NADDR. 
C 
   29   L(LPTR) = NNEW 
        LPTR = LPTR+1 
   32   L(LPTR) = N 
        LPTR = LPTR +1 
        IF(LL.EQ.LOLD) GOTO 33 
        LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        GOTO 32 
   33   L(LCURR+1) = LPTR-LCURR-2 
        IF(NCURR) 34,707,35 
  707   STOP 724 
   34   JADDR(JCURR) = LCURR 
        GOTO 36 
   35   NADDR(NCURR) = LCURR 
C 
C    Add NCURR to the temporary scratchpad contiguity list 
C    for NMOD. 
C 
   36   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NEAR to the temporary scratchpad 
C    contiguity list for NMOD, then the main inner loop is 
C    completed. 
C 
        IF(NCURR.NE.NEAR) GOTO 19 
C 
C    Pick up the beginning and size of the old contiguity 
C    list for NMOD, and garbage if necessary, checking 
C    that enough space can be recovered for the extended 
C    list. 
C 
        LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+2.LE.LEFF) GOTO 37 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+2.LE.LEFF) GOTO 37 
        IFLAG = 6 
        RETURN 
C 
C    Find top of old list.  Set up back-pointer and size 
C    for extended list, reset address pointer to it, and 
C    kill old list.  Copy old list across and add marker. 
C 
   37   LHI = LLO-1+LSIZE 
        L(LPTR) = NMOD 
        L(LPTR+1) = LSIZE+1+LTOP-LEFF 
        NADDR(NMOD) = LPTR 
        L(LLO-2) = 0 
        LPTR = LPTR+2 
        DO 38  LL = LLO,LHI 
        L(LPTR) = L(LL) 
   38   LPTR = LPTR+1 
        L(LPTR) = 0 
        LPTR = LPTR+1 
C 
C    Copy the scratchpad list down from its temporary 
C    position at the top of the heap. 
C 
        LEFF = LEFF+1 
        DO 39  LL = LEFF,LTOP 
        L(LPTR) = L(LL) 
   39   LPTR = LPTR+1 
C 
C    Increment the count of points accepted in the 
C    scratchpad tessellation; this completes the outer 
C    main loop of the scratchpad calculation. 
C 
        NINX = NINX+1 
    3   CONTINUE 
C 
C    The scratchpad calculation is complete; the information 
C    it provides now has to be combined with that in the 
C    main data base. 
C 
C    The first stage is to merge main and scratchpad 
C    lists for objects contiguous to NINDEX.  In order 
C    to deal properly with degeneracy, a record has to 
C    be maintained of the objects preceding and 
C    following the current one; the scan through the 
C    contiguity list is designed to facilitate this. 
C 
        LLX = NADDR(NINDEX)+1+LSIZEX 
        NMOD = L(LLX-1) 
        NNEXT = L(LLX) 
        DO 40  LKX = 1,LSIZEX 
        LLX = NADDR(NINDEX)+1+LKX 
        NLAST = NMOD 
        NMOD = NNEXT 
        NNEXT = L(LLX) 
C 
C    Is NMOD a constraint or a point? 
C 
        IF(NMOD) 41,710,42 
  710   STOP 725 
C 
C    Constraint.  Pick up beginning and size of contiguity 
C    list, and garbage if necessary, checking that enough 
C    space can be recovered for merged list. 
C 
   41   JMOD = -NMOD 
        LLO = JADDR(JMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = JADDR(JMOD)+2 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        IFLAG = 6 
        RETURN 
C 
C    Point.  Pick up beginning and size of contiguity 
C    list, and garbage if necessary, checking that enough 
C    space can be recovered for merged list. 
C 
   42   LLO = NADDR(NMOD)+2 
        LSIZE = L(LLO-1) 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLO = NADDR(NMOD)+2 
        IF(LPTR+LSIZE+1.LE.LTOP) GOTO 43 
        IFLAG = 6 
        RETURN 
C 
C    Find top of old list.  Save value of LPTR as 
C    address of new list, and set back pointer.  Kill old 
C    list. 
C 
   43   LHI = LLO-1+LSIZE 
        L(LLO-2) = 0 
        LMOD = LPTR 
        L(LPTR) = NMOD 
        LPTR = LPTR+2 
C 
C    Find NINDEX in the main part of the extended list, 
C    copying across as far as that point, but not 
C    including it.  Then space one past NINDEX. 
C 
        DO 44  LL = LLO,LHI 
        N = L(LL) 
        IF(N.EQ.NINDEX) GOTO 45 
        L(LPTR) = N 
   44   LPTR = LPTR+1 
        STOP 726 
   45   LL = LL+1 
C 
C    Find the marker, and space on one to the beginning 
C    of the scratchpad list. 
C 
        DO 46  LMARK = LL,LHI 
        IF(L(LMARK).EQ.0) GOTO 47 
   46   CONTINUE 
        STOP 727 
   47   LMARK = LMARK+1 
C 
C    Find NNEXT in the scratchpad list. 
C 
        DO 48  LL1 = LMARK,LHI 
        IF(L(LL1).EQ.NNEXT) GOTO 49 
   48   CONTINUE 
        STOP 730 
C 
C    Copy the scratchpad list into the merged list, starting 
C    after (or, in the case of degeneracy, at) NNEXT and 
C    stopping at NLAST.  A special case arises if NINDEX 
C    occurs at position LLO. 
C 
   49   IF(LL.NE.(LLO+1).AND.NNEXT.EQ.L(LPTR-1)) LPTR = LPTR-1 
   50   N = L(LL1) 
        L(LPTR) = N 
        LPTR = LPTR+1 
        LL1 = LL1+1 
        IF(LL1.GT.LHI) LL1 = LMARK 
        IF(N.NE.NLAST) GOTO 50 
C 
C    Copy the rest of the main list into the merged list, 
C    starting after (or, in the case of degeneracy, at) 
C    NLAST.  Deal with special case. 
C 
        LHI = LMARK-2 
        IF(LL.GT.LHI) LL = LLO 
        IF(L(LL).EQ.NLAST) LPTR = LPTR-1 
        IF(LL.EQ.LLO) GOTO 51 
        DO 52  LL1 = LL,LHI 
        L(LPTR) = L(LL1) 
   52   LPTR = LPTR+1 
        IF(L(LPTR-1).EQ.L(LMOD+2)) LPTR = LPTR-1 
C 
C    Set up the length of the merged list and the address 
C    pointer to it. 
C 
   51   L(LMOD+1) = LPTR-LMOD-2 
        IF(NMOD) 53,714,54 
  714   STOP 731 
   53   JADDR(JMOD) = LMOD 
        GOTO 40 
   54   NADDR(NMOD) = LMOD 
   40   CONTINUE 
C 
C    Merging is now complete.  The next step is to kill 
C    the scratchpad lists for constraints not contiguous 
C    to NINDEX.  This is done, exceptionally, in situ, by 
C    making the marker and scratchpad section look like 
C    a dead data item, and resetting the length to that for 
C    the main section.  Nothing need be done to those 
C    constraints contiguous to NINDEX, as they have 
C    already been dealt with in the merging operation. 
C 
        DO 55  J = 1,JCNS 
        LLO = JADDR(J)+2 
        LHI = LLO-1+L(LLO-1) 
        DO 56  LMARK = LLO,LHI 
        IF(L(LMARK).EQ.0) GOTO 57 
   56   CONTINUE 
        GOTO 55 
   57   L(LMARK+1) = LHI-LMARK-1 
        L(LLO-1) = LMARK-LLO 
   55   CONTINUE 
C 
C    Copy the contiguity list for NINDEX, garbaging if 
C    necessary and checking that enough space can be 
C    recovered for this list. 
C 
        LLOX = NADDR(NINDEX)+2 
        LSIZEX = L(LLOX-1) 
        IF(LPTR+LSIZEX+1.LE.LTOP) GOTO 58 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        LLOX = NADDR(NINDEX)+2 
        IF(LPTR+LSIZEX+1.LE.LTOP) GOTO 58 
        IFLAG = 6 
        RETURN 
   58   LHIX = LLOX-1+LSIZEX 
        LL = LPTR+2 
        DO 59  LLX = LLOX,LHIX 
        L(LL) = L(LLX) 
   59   LL = LL+1 
        L(LPTR) = NINDEX 
        L(LPTR+1) = LSIZEX 
C 
C    Kill the original version of the contiguity list 
C    for NINDEX, chain NINDEX back into the beginning 
C    of the free space in NADDR, decrement NPTSIN, and 
C    return with IFLAG set to zero. 
C 
        L(LLOX-2) = 0 
        NADDR(NINDEX) = -NFREE 
        NFREE = NINDEX 
        NPTSIN = NPTSIN-1 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TRYPT 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TRYPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE, 
     1    NSTART,NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT, 
     2    IFLAG,IGARB,XPT,YPT,NINDEX) 
C 
C    Finds the contiguities which the trial point (XPT,YPT) 
C    would have if it were added to the tessellation, but 
C    does not actualy add it. 
C 
C    The subroutine constructs, at the exit value of LPTR, 
C    a data item (which, being above LPTR, is unprotected) 
C    whose back pointer is unset and which gives the 
C    contiguity list which the trial point (XPT,YPT) 
C    would have if it were added to the tessellation. 
C    Existing contiguity lists are not modified, although 
C    data items may be relocated as a result of garbage 
C    collection.  As with ADDPT, a check is made that the 
C    trial point is adequately within the window.  Then 
C    the data structure is interrogated to determine 
C    whether there are any accepted points; if not, the 
C    contiguity list for the trial point is obtained simply 
C    by reversing the boundary list.  If there are some 
C    accepted points the (a) nearest neighbour to the trial 
C    point is found and a check is made that it is not too 
C    close to the trial point.  The contiguity list for the 
C    trial point is constructed in a temporary position at 
C    the top of the heap - it is not necessary to take 
C    account of degeneracies, since other contiguity lists 
C    are not modified.  On succesful completion of the 
C    construction, the list is copied down to its final 
C    position at LPTR, and return takes place with IFLAG 
C    set to zero.  Nonzero values of IFLAG indicate failure 
C    to construct the contiguity list for the trial point, 
C    as follows. 
C 
C       4  Trial point outside window 
C       5  Trial point duplicates accepted point 
C       6  Heap overflow 
C       8  NSTART not the index of an accepted point 
C 
C    Note that heap overflow is not here a catastrophic error. 
C    "Outside" and "duplicate" are each interpreted in terms 
C    of the tolerance values EPSCN and EPSPT in the argument 
C    list.  If there are no accepted points, return with 
C    IFLAG set to 8 cannot, of course, occur. NINDEX is 
C    the index of the (a) nearest neighbour to the trial 
C    point if there are any accepted points. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS), 
     1    NADDR(NPTS),L(LTOP) 
C 
C    Call subroutine SEEPT to check that the trial point 
C    is far enough inside the window; return if not, with 
C    NINDEX giving the index (flagged with a minus sign) 
C    of an unsatisfied constraint. 
C 
        CALL SEEPT(CN,JCNS,JADDR,EPSCN,IFLAG,XPT,YPT, 
     1    NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    If any points are accepted, jump forward to the main 
C    part of the routine 
C 
    1   IF(NPTSIN.GT.0) GOTO 2 
C 
C    There are no accepted points.  Check that there is 
C    enough space available in the heap for the contiguity 
C    list of the trial point. 
C 
        IF(LPTR+L(2)+1.LE.LTOP) GOTO 3 
        IFLAG = 6 
        RETURN 
C 
C    Construct the contiguity list for the trial point, and 
C    its length, as part of a data item at LPTR, leaving 
C    the back pointer of that data item unset. 
C 
    3   L(LPTR+1) = L(2) 
        LTRIAL = LPTR+1+L(2) 
        LHI = 2+L(2) 
        DO 4  LL = 3,LHI 
        L(LTRIAL) = L(LL) 
    4   LTRIAL = LTRIAL-1 
C 
C    Return with IFLAG set to zero. 
C 
        IFLAG = 0 
        RETURN 
C 
C    The main part of the routine deals with the case where 
C    there are some accepted points.  Call subroutine LOCPT 
C    to find the (a) nearest neighbour of the trial point. 
C    If it is too close, return.  The index of the nearest 
C    neighbour is set in NINDEX. 
C 
    2   CALL LOCPT(PT,NPTS,NADDR,NSTART,L,LTOP,EPSPT, 
     1    IFLAG,XPT,YPT,NINDEX,TSTVAL) 
        IF(IFLAG.EQ.0) GOTO 5 
        RETURN 
C 
C    Initialise NNEW, the first item to be placed on the 
C    contiguity list for the trial point; it is the 
C    neighbour of the trial point clockwise from its 
C    nearest neighbour. 
C 
    5   NCURR = NINDEX 
        LLO = NADDR(NCURR)+2 
        LHI = LLO-1+L(LLO-1) 
        XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
        NNEW = 0 
        DO 6  LL = LLO,LHI 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.LE.0) GOTO 6 
        IF(NNEW.EQ.0) GOTO 7 
        IF(TNEW-T) 6,51,7 
   51   IF(LL.EQ.LHI.AND.NNEW.EQ.L(LLO)) NNEW = N 
        GOTO 6 
    7   NNEW = N 
        TNEW = T 
    6   CONTINUE 
C 
C    Initialise LEFF to LTOP; the contiguity list for 
C    the trial point is constructed in a temporary position 
C    working down from LTOP. 
C 
        LEFF = LTOP 
C 
C    Enter the main loop.  First update NOLD and NCURR. 
C 
    8   NOLD = NCURR 
        NCURR = NNEW 
C 
C    Check that there is enough space to add another point 
C    to the temporary list for the trial point, and garbage 
C    if necessary, checking that enough space can be 
C    recovered. 
C 
        IF(LPTR+2.LE.LEFF) GOTO 9 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IGARB = IGARB+1 
        IF(LPTR+2.LE.LEFF) GOTO 9 
        IFLAG = 6 
        RETURN 
C 
C     Is NCURR a constraint or a point? 
C 
    9   IF(NCURR) 10,701,11 
  701   STOP 732 
C 
C    Constraint: (XCURR,YCURR) is the reflexion of the 
C    trial point in it.  Pick up beginning of contiguity 
C    list. 
C 
   10   JCURR = -NCURR 
        AJ = CN(1,JCURR) 
        BJ = CN(2,JCURR) 
        DJ = 2.0*(AJ*XPT+BJ*YPT+CN(3,JCURR))/(AJ*AJ+BJ*BJ) 
        XCURR = XPT-AJ*DJ 
        YCURR = YPT-BJ*DJ 
        LLO = JADDR(JCURR)+2 
        GOTO 12 
C 
C    Point, with coordinates (XCURR,YCURR): pick up 
C    beginning of contiguity list. 
C 
   11   XCURR = PT(1,NCURR) 
        YCURR = PT(2,NCURR) 
        LLO = NADDR(NCURR)+2 
C 
C    Find the end of the contiguity list for NCURR and then, 
C    using symmetry of the contiguity relation, find NOLD 
C    in it. 
C 
   12   LHI = LLO-1+L(LLO-1) 
        DO 13  LOLD = LLO,LHI 
        IF(L(LOLD).EQ.NOLD) GOTO 14 
   13   CONTINUE 
        STOP 733 
C 
C    Update NNEW, the item to be added to the contiguity 
C    list for the trial point after NCURR; it is the 
C    neighbour of the trial point clockwise from NCURR. 
C 
   14   NNEW = 0 
        LL = LOLD 
   15   LL = LL+1 
        IF(LL.GT.LHI) LL = LLO 
        N = L(LL) 
        CALL TT(CN,JCNS,PT,NPTS,XPT,YPT,XCURR,YCURR,N,T,IWIPE) 
        IF(IWIPE.GT.0) GOTO 16 
        IF(NNEW) 17,15,17 
   16   IF(NNEW.EQ.0) GOTO 18 
        IF(TNEW.LE.T) GO TO 17 
   18   NNEW = N 
        TNEW = T 
        GOTO 15 
C 
C    Add NCURR to the contiguity list for the trial point. 
C 
   17   L(LEFF) = NCURR 
        LEFF = LEFF-1 
C 
C    If we have just added NINDEX to the contiguity list 
C    for the trial point, then the main loop is completed. 
C 
        IF(NCURR.NE.NINDEX) GOTO 8 
C 
C    Move the contiguity list for the trial point down to 
C    LPTR, as part of a data item at this location, and 
C    set up its size. 
C 
        L(LPTR+1) = LTOP-LEFF 
        LTRIAL = LPTR+2 
        LEFF = LEFF+1 
        DO 19  LL = LEFF,LTOP 
        L(LTRIAL) = L(LL) 
   19   LTRIAL = LTRIAL+1 
C 
C    The construction of the contiguity list for the trial 
C    point is complete.  Set IFLAG to zero and return. 
C 
        IFLAG = 0 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE TILE4M 4.2 DATED 13 JANUARY 1981 
C 
        SUBROUTINE TILE4M(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, 
     1    NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB) 
C 
C    A master routine for constructing a tessellation within an 
C    arbitrary window.  This routine gives access to the full TILE4 
C    database, and in new applications should be used in preference 
C    to subroutine TILE, which is supplied for compatibility with the 
C    routine of that name in TILE3.  Subroutine TILE4M can be followed 
C    by calls to ADDPT,SUBPT to manipulate the database, but this is 
C    not possible after a call to TILE.  Note that it is legitimate to 
C    call TILE4M with a temporary small value of NPTS in such cases - 
C    see below for details. 
C 
C    The user call to TILE4M must supply all the arguments.  The 
C    values of JCNS,NPTS,LTOP,EPSCN,EPSPT may be given as constants in 
C    the call, but all other arguments must be array or simple variable 
C    names, as appropriate.  The arguments are as follows. 
C 
C      CN      REAL array of dimension (3,JCNS).  User supplies 
C              constraints according to the convention 
C 
C                      CN(1,J)*x + CN(2,J)*y + CN(3,J) .LT. 0.0 
C 
C              for J = 1,...,JCNS. 
C 
C      JCNS    INTEGER variable in which user supplies number of 
C              constraints and also second dimension of CN and 
C              dimension of JADDR. 
C 
C      JADDR   INTEGER array of dimension JCNS supplied by user as 
C              uninitialised space for part of the TILE4 database. 
C              Contains heap addresses of constraint data items 
C              on return. 
C 
C      PT      REAL array of dimension (2,NPTS).  User supplies 
C              point coordinates according to the convention 
C 
C                      (PT(1,N),PT(2,N))  is  (x,y) 
C 
C              for N = 1,...,NPTS. 
C 
C      NPTS    INTEGER variable in which user supplies number of 
C              points and also second dimension of PT and 
C              dimension of NADDR.  In applications where a call to 
C              TILE4M is to be followed by subsequent manipulation 
C              of the database, it may be desirable to set up PT(2, ) 
C              and NADDR( ) to a larger true dimension in the 
C              calling program than the dimension actually passed 
C              to TILE4M.  This is quite legitimate provided that 
C              the call to TILE4M is followed by a call to 
C              subroutine EXTEND to initialise the remainder of 
C              NADDR before calls to ADDPT,SUBPT are made with 
C              the larger value of NPTS in the calling program. 
C 
C      NADDR   INTEGER array of dimension NPTS supplied by user 
C              as uninitialised space for part of TILE4 database. 
C              Contains heap addresses of point data items on 
C              return. 
C 
C      NFREE   INTEGER variable not initialised by user.  Returns a 
C              pointer to the first free location in PT,NADDR, or 
C              returns NPTS+1 if there are none. 
C 
C      NSTART  INTEGER variable holding index of starting-point for 
C              nearest-neighbour walk.  User initialisation controls 
C              updating as follows.  If NSTART is positive or zero on 
C              entry, then it is updated occasionally to lie near the 
C              centroid of the accepted points.  If it is negative on 
C              entry, then each nearest-neighbour walk commences from 
C              the previous accepted point.  Use 0,-1 as initial values 
C              by convention.  Return value is the index of the next 
C              starting-point which would have been used. 
C              The 0 option should be used when the order of the points 
C              in the array PT bears little relation to their actual 
C              position in the plane.  The -1 option should be used 
C              when the points in nearby locations in the array PT 
C              are likely to have nearby positions in the plane. 
C 
C      NPTSIN  INTEGER variable not initialised by user.  Returns the 
C              number of currently accepted points. 
C 
C      L       INTEGER (or INTEGER*2) array of dimension LTOP. 
C              Supplied by user as uninitialised space for part of the 
C              TILE4 database.  Holds a heap of data items on return. 
C 
C      LTOP    INTEGER variable in which user supplies dimension 
C              of L.  A value of 9*NPTS, or 1000 if that is greater, 
C              is recommended, although a larger value may be 
C              needed if the database is to be manipulated 
C              subsequently by ADDPT,SUBPT. 
C 
C      LPTR    INTEGER variable not initialised by user.  Returns a 
C              pointer to the first free location in L. 
C 
C      LBASE   INTEGER variable not initialised by user.  Returns a 
C              pointer to the base of the working space in L. 
C 
C      EPSCN   REAL variable in which user supplies tolerance 
C              value for "inside" window.  Negative sign is treated 
C              as zero. 
C 
C      EPSPT   REAL variable in which user supplies tolerance 
C              value for "identical" points.  Negative sign is treated 
C              as zero. 
C 
C      IFLAG   INTEGER variable not initialised by user.  Returns 
C              value zero on completely successful exit from TILE4M. 
C              Nonzero values indicate complete or partial failure 
C              of TILE4M as follows. 
C 
C              1  Dud constraint - CN(1,J) and CN(2,J) both zero 
C                 for some value of J.  TILE4M fails. 
C 
C              2  Unbounded window - constraints inadequate. 
C                 TILE4M fails. 
C 
C              3  Empty window - constraints inconsistent. 
C                 TILE4M fails. 
C 
C              4  Not all points "inside" window.  TILE4M succeeds 
C                 by rejecting such points. 
C 
C              5  Some points "identical".  TILE4M succeeds by 
C                 rejecting points which duplicate an already 
C                 accepted point.  If both errors 4 and 5 occur, 
C                 error 4 is reported. 
C 
C              6  Heap overflow - LTOP too small.  TILE4M fails. 
C 
C      IGARB   INTEGER variable not initialised by user.  Returns 
C              number of garbage collections necessarily carried out 
C              in TILE4M, that is, excluding the final "cosmetic" 
C              garbage collection carried out immediately before 
C              control is returned to the calling program. 
C 
C    On return from TILE4M, the TILE4 database is held as explained 
C    above.  Some of the values returned are of no direct interest 
C    to the user, and are returned simply so that they can be passed 
C    on to other routines.  The user resets items in the database 
C    at his own risk, and can easily produce chaos by attempting to 
C    do so other than by calling routines supplied as part of the 
C    TILE4 package.  Users intending to access the database directly, 
C    rather than via the supplied interrogation routines, will need 
C    to understand its structure, which is as follows.  Each 
C    effective constraint, and each accepted point, has a contiguity 
C    list of neighbouring points and constraints held as part of a 
C    data item in the heap L.  The data item for constraint J is 
C    at LL = JADDR(J), and L(LL) is set to -J.  The data item for 
C    point N is at LL = NADDR(N), and L(LL) is set to N.  Ineffective 
C    constraints have JADDR(J) set to zero.  Rejected points have 
C    NADDR(N) set to a negative value.  These negative values when 
C    stripped of the sign form a chain entered at NFREE and terminated 
C    by (-)(NPTS+1).  A data item in L consists of the back pointer 
C    already described, then a location holding the length of the 
C    succeeding list, then the contiguity list itself, in as many 
C    locations (always at least 3) as are required, with points 
C    referenced by their index and constraints by their index with 
C    a minus sign flagging them.  The order of objects within 
C    the lists is anticlockwise cyclic.  Storage in L is administered 
C    as a heap, that is to say, data items are not modified 
C    in situ but are replaced by new data items at the base of the 
C    free space.  Abandoned data items are flagged as dead by 
C    having their back-pointers set to zero, and subroutine GARBAJ 
C    returns the space occupied by such items to the top 
C    of the heap as free space when required.  At the very bottom 
C    of the heap is a special data item protected against 
C    removal by GARBAJ.  This has the following form.  L(1) 
C    is unset.  L(2) is the length of the succeeding list. 
C    Succeeding entries in L, as many as are needed, hold a 
C    list of effective constraints (each as usual flagged 
C    with a minus sign) in clockwise cyclic order round the 
C    window boundary.  LBASE is the first location above this data 
C    item. 
C 
C2      INTEGER*2 L 
        DIMENSION CN(3,JCNS),JADDR(JCNS),PT(2,NPTS),NADDR(NPTS),L(LTOP) 
C 
C    Hard error for nonpositive dimensions. 
C 
        IF(JCNS.LE.0.OR.NPTS.LE.0.OR.LTOP.LE.0) STOP 734 
C 
C    Set up the window and check its validity. 
C 
        CALL WINLD(CN,JCNS,JADDR,L,LTOP,LBASE,IFLAG) 
        IF(IFLAG.EQ.0) GOTO 1 
        RETURN 
C 
C    Save information from NSTART, then complete initialisation. 
C 
    1   IRESET = 0 
        IF(NSTART.LT.0) IRESET = -1 
        CALL CLEAR(NPTS,NADDR,NFREE,NSTART,NPTSIN,LPTR,LBASE,IGARB) 
        IFSAVE = 0 
        NRESET = 2 
C 
C    Scan the points one by one, attempting to insert each. 
C 
        DO 2  NPT = 1,NPTS 
        XPT = PT(1,NPT) 
        YPT = PT(2,NPT) 
        CALL ADDPT(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, 
     1    L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX) 
        IF(IFLAG.EQ.7) STOP 735 
        IF(IFLAG.EQ.8) STOP 736 
        IF(IFLAG.NE.6) GOTO 3 
        RETURN 
C 
C    Nothing disastrous has happened.  But something mildly nasty 
C    may have done. 
C 
    3   IF(IFLAG.EQ.0) GOTO 4 
        IF(IFLAG.EQ.4) IFSAVE = 4 
        IF(IFLAG.EQ.5.AND.IFSAVE.EQ.0) IFSAVE = 5 
        NSAVE = -NADDR(NPTS) 
        NADDR(NPTS) = -NPT 
        NFREE = -NADDR(NPT) 
        NADDR(NPT) = -NSAVE 
        GOTO 2 
C 
C    The point has been accepted.  Deal with NSTART. 
C 
    4   IF(IRESET.EQ.-1) GOTO 5 
        IF(NPTSIN.LT.NRESET) GOTO 2 
        NRESET = 2*NRESET 
        CALL MIDPT(PT,NPTS,NADDR,NSTART,NPTSIN,L,LTOP,IFLAG, 
     1    XBAR,YBAR,NINDEX) 
        IF(IFLAG.EQ.4) STOP 737 
        IF(IFLAG.EQ.8) STOP 740 
    5   NSTART = NINDEX 
    2   CONTINUE 
C 
C    Final garbage collection, set IFLAG, and return. 
C 
        CALL GARBAJ(JCNS,JADDR,NPTS,NADDR,L,LTOP,LPTR,LBASE) 
        IFLAG = IFSAVE 
        RETURN 
        END 
C 
C*********************************************************************** 
C 
C 
C    FORTRAN TILE 4 CREATED 1 SEPTEMBER 1977 
C    SUBROUTINE EXTEND 4.3 DATED 13 JANUARY 1981 
C 
        SUBROUTINE EXTEND(NADDR,NPTS,NSUB) 
C 
C    Extends the free space chain in NADDR from NSUB to NPTS. 
C 
        DIMENSION NADDR(NPTS) 
        IF(NSUB.GT.NPTS) STOP 741 
        IF(NSUB.EQ.NPTS) RETURN 
        NS1 = NSUB+1 
        DO 1  N = NS1,NPTS 
    1   NADDR(N) = -(N+1) 
        RETURN 
        END 
C 
C*********************************************************************** 
