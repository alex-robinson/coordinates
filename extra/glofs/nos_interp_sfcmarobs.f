C*****************************************************************************
C  PROGRAM Namme:   nos_interp_sfcmarobs.f (previous name  INTERPNN)
C
c  Location: /nwprod/sorc/nos_interp_sfcmarobs.fd 
c
C  Purpose: To interpolate meteorological data to grid points for lake
C   circulation models, wave models, and storm surge models
C
C  Algorithm: This program uses the natural neighbor technique (Sambridge et
C  al., 1995) to interpolate observations of wind, temperature, dew point,
C  and cloud cover to a grid covering one of the Great Lakes.  Observations
C  are obtained from a file created by the 'metedit' program.  All
C  observations have been adjusted to a common anemometer height and
C  observations from overland stations have been adjusted to be more
C  representative of overwater conditions.  The observations are
C  interpolated to a grid covering the lake using a natural neighbor
C  technique. Interpolated values are generated for each grid square in the
C  lake (defined by non-zero depth in the bathymetric data array).
C  Interpolated fields are stored in XDR (machine independent) format for
C  use in circulation models, wind wave prediction models, and storm surge
C  models. If there is not at least one observation of a parameter available
C  for a particular time step, the field for that time step for that
C  parameter is filled in by linear interpolation between the last available
C  field and next available field.
C
C  REFERENCE:  Sambridge, M., Braun, J., and H. McQueen, 1995.  Geophysical
C   parameterization and interpolation of irregular data using natural
C   neighbors.  Geophys. J. Int., 122, 837-857.
C
C   Input files:
C    Standard input (FORTRAN free format):  
C     Record 1 - File name for bathymetric (interpolation) grid (see
C        subroutine RGRID for a description of the file format)
C     Record 2 - File name  for 'metedit' file from which to read the
C        meteorological data records.
C     Record 3 - File name for output files.  Six files are created,
C        ---.wu = wind u-component (XDR format)
C        ---.wv = wind v-component (XDR format)
C        ---.at = air temperature (XDR format)
C        ---.dp = dew point (XDR format)
C        ---.cl = cloud cover 0=NONE 1=COMPLETE (XDR format)
C     Record 4 - Year, Julian day, and hour to start interpolation
C     Record 5 - Time interval (hours) between gridded fields
C     Record 6 - Total number of time intervals (number of fields
C        created in XDR files is this number plus one)
C
C   Output: Interpolated fields for wind u-component, wind v-component,
C     air temperature, cloud cover, and dew point are created in XDR format.
C     Also a printed summary of input parameters, various interpolation
C     actions, and the min, max, and mean of all interpolated fields.
C
C  Subroutines :  
C 
C        FILL - Fills gaps in interpolated fields by linear interpolation
C        INTPW - Interpolotes wind vector
C        INTPT - Interpolates scalar fields
C        NNINTERP - Sets up variables required for Sambridge and
C          Braun natural neighbor subroutine
C        RGRID - Reads the bathymetric data file
C        XDIST - Returns x-distance in grid coordinates for a
C          given latitude and longitude
C        YDIST - Returns y-distance in grid coordinates for a
C          given latitude and longitude
C        WXDR2 - Write a 2-d field to XDR file 
C        RXDR2 - Read a 2-d field from XDR file 
C
C  In addition to the subroutines included in this program, the
C   Sambridge and Braun routines in files 'nnq.f', 'stack.c',
C   'stackpair.c', and 'volume.c' must be compiled and linked to
C   this program, as well as the XDR Fortran interface routines.
C
C   On the HP C160 system (PA RISK 8000) the compiler command line is:
C    
C    f77 -o interpnn -O +DA2.0 +DS2.0 +Odataprefetch +Onomoveflops +E1 \
C       interpnn.f nnq.o stack.o stackpair.o volume.o hpxdrc.o -lm 
C
C   On other HP systems the compiler command lines are:
C    
C    F77 -O +E1 nnq.f stack.o stackpair.o volume.o -lm
C    F77 -o interpnn.x -O +E1 interpnn.f nnq.o stack.o stackpair.o volume.o hpxdrc.o -lm
C
C  History :  
C     Written by D.J. Schwab, NOAA Great Lakes Environmental Research
C      Laboratory, Ann Arbor, MI.  March, 1999.
c
c  Revisions:
c     6/6/2005   J. Kelley  Addition code to prevent wind speed and wind comps
c                           calculations when u and v comps have values of
c                           -9999. (missing value flag)      
c
c     2/2005-7/2005 G. Mott (CO-OPS) Added several additions for netcdf output,
c                    labelled as such 
c
c     1/24/2013  J.Kelley   Minor changes to the doc block made in
c                           preparation for migrating GLOFS to WCOSS
c  -----------------------------------------------------------------------------
C

C  INCLUDE STATEMENTS
      include 'rgrid.h'
      include 'xdist.h'
      include 'ydist.h'

C  PARAMETERS FOR MAX GRID DIMENSIONS, MAX OBS, MAX STATIONS, MAX HOURS
      PARAMETER(IDIM=131,JDIM=251,MAXOBS=250000,MAXHRS=8785)
      
cc julian and gregorian *8 variables

      real time
      real gregtime
      integer fc

      real*8 julian
      real*8 jfirst, jlast,jout,jbase
      real*8 ymet,monmet,daymet,hrmet
      real*8 fyra, fmon, fday, fhr
      real*8 YR,MN,DY,HR
      real*8 daynum
      real*8 bd(4),daysince
      real timewrite
      character*19 validtime
      character chyear*4,chmon*2,chdy*2,chhr*2

      integer lat,lon,i,j
      integer basedate(4),IYR,IMN,IDY,IHR

      LOGICAL DBERR
      COMMON /INT1/ D(IDIM,JDIM),SINTP(IDIM,JDIM),
     1 UINTP(IDIM,JDIM),VINTP(IDIM,JDIM)
      COMMON /INT2/ S(MAXOBS),T(MAXOBS),U(MAXOBS),V(MAXOBS),
     1 X(MAXOBS),Y(MAXOBS),AT(MAXOBS),DP(MAXOBS),CC(MAXOBS)
      COMMON /INT3/ DTI,IM,JM,DS,IDT,NTOT,NPTS,NRMAX,NDX
      COMMON /INT4/ FGRID,FINPUT,FOUTP,FNAME

      COMMON /GPARM/ RPARM(23),IPARM(54)

      CHARACTER*80 FGRID,FINPUT,FOUTP,FNAME,LATLONGRID

      INTEGER IERRWS(MAXHRS),IERRAT(MAXHRS),IERRDP(MAXHRS),
     1  IERRCC(MAXHRS)

C READ CONTROL FILE

C  GET GRID FILE NAME
      READ(5,'(A)') FGRID

C  GET NAME OF INPUT MAROBS FILE
      READ(5,'(A)') FINPUT

C  GET START DATE AND TIME
C      READ(5,*) IYRS,IJDS,IHRS
      READ(5,*) IYRS,IMONS,IDAYS,IHRS
C
      IMNS=0

C  GET TIME INTERVAL FOR GRIDDED FIELDS (HOURS)
      READ(5,*) DTI

C  GET NUMBER OF FIELDS TO INTERPOLATE
      READ(5,*) NTOT

C GET LAT LON GRID SIZE
      READ(5,*) lon, lat

      NHRS=IFIX(DTI*NTOT+0.5)
      IDT=IFIX(DTI+0.5)
C
C  CALCULATE JULIAN DAY OR DAY OF YEAR NUMBER FOR REQUESTED START DATE AND TIME
      CALL DBJDAT(IYRS,IMONS,IDAYS,IJDS,DBERR)
C
C  PRINT INTERPOLATION PARAMETERS

      WRITE(6,'('' ------MET OBS INTERPOLATION PROCEDURE-------'')')
      WRITE(6,'('' GRID FILE: '',A)') FGRID
      WRITE(6,'('' MET DATA FILE: '',A)') FINPUT
      WRITE(6,'('' LAT_LON GRID FILE: '',A)') LATLONGRID
C      WRITE(6,'('' START TIME (YR JD HR): '',I5,I4,I3)') IYRS,IJDS,IHRS
      WRITE(6,'('' START TIME (YR MON DY JD HR):'',
     1 I4,1X,I2,1X,I2,1X,I3,1X,I3)') IYRS,IMONS,IDAYS,IJDS,IHRS
      WRITE(6,'('' TIME INTERVAL FOR INTERPOLATED FIELDS (HR):''
     1   ,I3)') IDT
      WRITE(6,'('' TOTAL TIME (HR):'',I7)') NHRS

C  GET BATHYMETRIC GRID
      OPEN(1,FILE=FGRID,STATUS='OLD')
      CALL RGRID(1,D,IDIM,JDIM)
      IM=IPARM(1)
      JM=IPARM(2)
      DS=RPARM(3)
      CLOSE(1)

C  COUNT NUMBER OF GRID POINTS
      NS=0
      DO 5 I=1,IM
      DO 5 J=1,JM
       IF(D(I,J).GT.0) NS=NS+1
   5  CONTINUE

cc If not equal will give run time error

       if (lon .eq. IM) then
        
       write(*,*) 'IM,JM correct size for code to proceed'

      else
      write(*,*) 'EXIT IM JM different from IM,
     *JM read in from bath file'
      end if

C  GET OBSERVATIONS FROM SURFACE MARINE OBSERVATION FILE CREATED BY METEDIT
      OPEN(2,FILE=FINPUT,STATUS='OLD')
      I=1
  10  READ(2,'(20X,F10.2,2F8.0,6F6.1)',END=15) T(I),X(I),Y(I),
     1   U(I),V(I),S(I),AT(I),DP(I),CC(I)

      write(6,8888)T(I),X(I),Y(I),
     1   U(I),V(I),S(I),AT(I),DP(I),CC(I)


8888  format(20X,F10.2,2F8.0,6F6.1)
c
ccc  g.l. 10/19/1999
c   CONVERT CLOUD COVER FROM PERCENT (0-100) TO FRACTION (0-1)
      if (CC(I).ne.-99.9) CC(I)=CC(I)/100.
ccc
      I=I+1
      GOTO 10
  15  NPTS=I-1


       YR=DBLE(IYRS)
       MN=DBLE(IMONS)
       DY=DBLE(IDAYS)
       HR=DBLE(IHRS)

       jfirst= julian(YR,MN,DY,HR)
       jbase= julian(bd(1),bd(2),bd(3),bd(4))
       
       daysince=jfirst-jbase

C  OPEN OUTPUT FILES FOR WU, WV, AT, DP, CL

       OPEN(51,FILE='uwind.txt')
       OPEN(52,FILE='vwind.txt')
       OPEN(53,FILE='atemp.txt')
       OPEN(54,FILE='dtemp.txt')
       OPEN(55,FILE='cldcv.txt')


C  INTERPOLATION LOOP
      DO 500 IT=0,NTOT
      TIME=IT*DTI
      
 
c gregs addition to keep track of time
       jout = jfirst + IT*DTI/24.

cc compute calender date from julian day decimal
      call gregorian(jout,ymet,monmet,daymet,hrmet)

cc compute day of year (fake julian day)
      call daynumber(ymet,monmet,daymet,daynum,dberr)

c      write(*,*) 'model times',jout,NINT(ymet),NINT(monmet),
c     *                 NINT(daymet),NINT(hrmet),NINT(daynum)

c count for output file info
       fc=IT+1

c end gregs additions

C  WIND SPEED AND WIND COMPONENTS

   
      CALL INTPW(TIME,IERR)


      IERRWS(IT+1)=IERR
      DO 300 I=1,IM
      DO 300 J=1,JM
      IF(D(I,J).GT.0.) THEN
c
c   begin Kelleys additions
c   PREVENT CALCULATIONS WHEN INTERPOLATION YIELDS U AND V WIND COMPS 
c   FIELD CONTAINING MISSING VALUE FLAGS (-9999.)
       IF(UINTP(I,J).NE.-9999..AND.VINTP(I,J).NE.-9999.)THEN
        DENOM=SQRT(UINTP(I,J)**2+VINTP(I,J)**2)
        UINTP(I,J)=SINTP(I,J)*UINTP(I,J)/DENOM
        VINTP(I,J)=SINTP(I,J)*VINTP(I,J)/DENOM
       ENDIF
c   end Kelleys additions    
c 
      ENDIF
300   CONTINUE

      write(51,105)fc,NINT(ymet),NINT(monmet),NINT(daymet),
     1                              NINT(hrmet),NINT(daynum)
      write(51, "(10F10.2)")((UINTP(I,J), I=1, IM), J=1, JM)

      write(52,105)fc,NINT( ymet),NINT(monmet),NINT(daymet),
     2                            NINT(hrmet),NINT(daynum)
       write(52, "(10F10.2)")((VINTP(I,J), I=1, IM), J=1, JM)



C  AIR TEMPERATURE
       CALL INTPT(TIME,AT,'AT',IERR)
       IERRAT(IT+1)=IERR

      write(53,105)fc,NINT( ymet),NINT(monmet),NINT(daymet),
     1                            NINT(hrmet),NINT(daynum)
      write(53, "(10F10.2)")((SINTP(I,J), I=1, IM), J=1, JM)

C  DEW POINT
       CALL INTPT(TIME,DP,'DP',IERR)
       IERRDP(IT+1)=IERR

      write(54,105)fc, NINT( ymet),NINT(monmet),NINT(daymet),
     1                            NINT(hrmet),NINT(daynum)
      write(54, "(10F10.2)")((SINTP(I,J), I=1, IM), J=1, JM)

C  CLOUD COVER
       CALL INTPT(TIME,CC,'CC',IERR)
       IERRCC(IT+1)=IERR

      write(55,105)fc,NINT( ymet),NINT(monmet),NINT(daymet),
     1                            NINT(hrmet),NINT(daynum)
       write(55, "(10F10.2)")((SINTP(I,J), I=1, IM), J=1, JM)

        timewrite=daysince+time/24.

c convert integer to character

      write(chyear,71) NINT(ymet)
      write(chmon,72)  NINT(monmet)
      write(chdy,72)   NINT(daymet)
      write(chhr,72)   NINT(hrmet)

  71   format(I4)
  72   format(I2.2)

      validtime= chyear//'/'//chmon//'/'//chdy//':'//chhr
     *              //':00'//':00'

500   CONTINUE


c format statement for ascii date time one liner for each record

  105 format(i3.3,1x,i4,1x,i2.2,1x,i2.2,1x,i2.2,1x,'00',1x,
     *i3.3,' LAKE YYYY MM DD HH min DAYN')

C  FILL IN MISSING FIELDS BY LINEAR INTERPOLATION

C  CALCULATE STATS

      STOP
      END


C********************************************************************
      SUBROUTINE INTPW(TIME,IERR)

C  INTERPOLATE GRIDDED WIND FIELD FROM POINT OBSERVATIONS

C  PARAMETERS FOR MAX GRID DIMENSIONS, MAX OBS, MAX STATIONS, MAX HOURS

      PARAMETER(IDIM=131,JDIM=251,MAXOBS=250000,MAXHRS=8785)

      COMMON /INT1/ D(IDIM,JDIM),SINTP(IDIM,JDIM),
     1 UINTP(IDIM,JDIM),VINTP(IDIM,JDIM)
      COMMON /INT2/ S(MAXOBS),T(MAXOBS),U(MAXOBS),V(MAXOBS),
     1 X(MAXOBS),Y(MAXOBS),AT(MAXOBS),DP(MAXOBS),CC(MAXOBS)
      COMMON /INT3/ DTI,IM,JM,DS,IDT,NTOT,NPTS,NRMAX,NDX
      COMMON /INT4/ FGRID,FINPUT,FOUTP,FNAME

      COMMON /GPARM/ RPARM(23),IPARM(54)  

      CHARACTER*80 FGRID,FINPUT,FOUTP,FNAME

      REAL XX(1000),YY(1000),SS(1000),UU(1000),VV(1000)



C  DETERMINE TIME CUTOFF

      DTMAXS=DTI/2.
      
100   NSTART=1
      NSTOP=1
      DO 210 J=NSTART,NPTS
      IF(T(J).GE.TIME-DTMAXS) GO TO 220
210   CONTINUE
220   CONTINUE
      NSTART=J
      DO 230 J=NSTOP,NPTS
      IF(T(J).GT.TIME+DTMAXS) GO TO 240
230   CONTINUE
240   CONTINUE
      NSTOP=J-1

      NP=0
      DO 450 K=NSTART,NSTOP
      IF(S(K).EQ.-99.9) GO TO 450
      NP=NP+1
      XX(NP)=X(K)/DS
      YY(NP)=Y(K)/DS
      SS(NP)=S(K)
      UU(NP)=U(K)
      VV(NP)=V(K)
450   CONTINUE

      IF(NP.GT.0) THEN


       CALL NNINTERP(NP,XX,YY,SS,D,IDIM,IM,JM,SINTP)
       CALL NNINTERP(NP,XX,YY,UU,D,IDIM,IM,JM,UINTP)
       CALL NNINTERP(NP,XX,YY,VV,D,IDIM,IM,JM,VINTP)


       IERR=0
      ELSE
       IERR=1
       DO 470 I=1,IM
       DO 470 J=1,JM
        SINTP(I,J)=-9999.
        UINTP(I,J)=-9999.
        VINTP(I,J)=-9999.
470    CONTINUE
       PRINT *,' ***NO OBSERVATIONS FOR WS AT HOUR: ',TIME
      ENDIF     
      RETURN
      END


C********************************************************************
      SUBROUTINE INTPT(TIME,TEMP,ID,IERR)

C  INTERPOLATE GRIDDED TEMPERATURE FIELD (AIR TEMP OR DEW POINT OR CLOUD COVER)
C    FROM POINT OBSERVATIONS

C  PARAMETERS FOR MAX GRID DIMENSIONS, MAX OBS, MAX STATIONS, MAX HOURS

      PARAMETER(IDIM=131,JDIM=251,MAXOBS=250000,MAXHRS=8785)

      COMMON /INT1/ D(IDIM,JDIM),SINTP(IDIM,JDIM),
     1 UINTP(IDIM,JDIM),VINTP(IDIM,JDIM)
      COMMON /INT2/ S(MAXOBS),T(MAXOBS),U(MAXOBS),V(MAXOBS),
     1 X(MAXOBS),Y(MAXOBS),AT(MAXOBS),DP(MAXOBS),CC(MAXOBS)
      COMMON /INT3/ DTI,IM,JM,DS,IDT,NTOT,NPTS,NRMAX,NDX
      COMMON /INT4/ FGRID,FINPUT,FOUTP,FNAME

      COMMON /GPARM/ RPARM(23),IPARM(54)  

      REAL TEMP(MAXOBS)
      CHARACTER*80 FGRID,FINPUT,FOUTP,FNAME
      CHARACTER*2 ID

      REAL XX(1000),YY(1000),SS(1000)

C  DETERMINE TIME CUTOFF

      DTMAXS=DTI/2.
      
100   NSTART=1
      NSTOP=1
      DO 210 J=NSTART,NPTS
      IF(T(J).GE.TIME-DTMAXS) GO TO 220
210   CONTINUE
220   CONTINUE
      NSTART=J
      DO 230 J=NSTOP,NPTS
      IF(T(J).GT.TIME+DTMAXS) GO TO 240
230   CONTINUE
240   CONTINUE
      NSTOP=J-1

      NP=0
      DO 450 K=NSTART,NSTOP
c      IF(TEMP(K).EQ.999.9) GO TO 450
      IF(TEMP(K).EQ.-99.9) GO TO 450
      NP=NP+1
      XX(NP)=X(K)/DS
      YY(NP)=Y(K)/DS
      SS(NP)=TEMP(K)
      
450   CONTINUE

      IF(NP.GT.0) THEN
       CALL NNINTERP(NP,XX,YY,SS,D,IDIM,IM,JM,SINTP)
       IERR=0
      ELSE
       IERR=1
       DO 470 I=1,IM
       DO 470 J=1,JM
        SINTP(I,J)=-9999.
470    CONTINUE
       PRINT *,' ***NO OBSERVATIONS FOR '//ID//
     1   ' AT HOUR: ',TIME
      ENDIF     
      RETURN
      END

c-----------------------------------------------------------------------------
c  Subroutine to use Sambridge and Braun Natural Neighbors interpolation
c
c   Input to this subroutine is a list of the sample points and a
c   description of the interpolation grid.  The coordinates of the sample
c   points are slightly perturbed to avoid co-linear points. Four
c   artificial sample points with the mean sample value are added
c   to form a large bounding rectangle.    DJS - 8/98
c-----------------------------------------------------------------------------
c
c	Explanation of parameters:
c	Parameter	 Meaning			Used by 
c	np_max		:Max. number of data nodes.	All
c	nd_max		:Max. number of dimensions.     All
c	nt_max		:Max. number of triangles.      All
c	nh_max		:Max. number of triangles.      nn2d_setup 
c			 on convex hull.		(can be set to 1 if
c			                                extension to outside of
c							convex hull is not used)
c	nnpn_max	:Max. number of nearest		All (used by recurisve 
c			 neighbours per node.		     routines) 
c	nmax		:Max. number of nearest 	setup 
c			 neighbours per node times 
c			 maximum number of nodes.
c       nv_max          :Max. number of triangles       Only used by incremental
c                        visible from any point         Delaunay routine delaun.
c                        outside convex hull.		in nn2d_setup           
c	eps		:Tolerance used by delaun routine
c			 (see subroutine header)
c
c					M. Sambridge, RSES,  Oct. 1995.
c					(Last revised April 1996)
c
c-----------------------------------------------------------------------------
c                              
      subroutine nninterp(np,x,y,z,d,idim,im,jm,zgrid)
      
      real x(*),y(*),z(*),d(idim,*),zgrid(idim,*)

      parameter		(np_max=1000,nd_max=2,nt_max=2*np_max)
      parameter		(nh_max=1000,nnpn_max =50,nmax=3*nt_max+np_max)
      parameter		(nv_max=1000)
c
c						arrays and variables 
c						used by nn 2-D routines
c
      real*8		points(2,np_max)
      real*8		centres(3,nt_max)
      real*8		fvals(np_max)
      real*8		f,df(2),ddf(3),v
      real*8		xd,yd
      real*8		eps
      integer		vertices(3,nt_max)
      integer		neighbour(3,nt_max)
      integer		hulltriangles(nh_max)
      integer           int_method
      integer           dmode,nmode
      logical		out
      logical           nnext
      logical		rescale 
      logical		clockwise
      logical		consistent
c					
c						work arrays used by nn 2-D
c
      integer           vis_tlist(nv_max)
      integer           vis_elist(nv_max)
      integer           add_tlist(nv_max)
      integer           nnn(np_max+1)
      integer           nnlist(nmax)
      integer           ntrilist(nmax)
c                                               NN work arrays
      real*8            work_d1(2,nnpn_max)
      real*8            work_d2(2,nnpn_max)
      real*8            work_d3(2,nnpn_max)
      real*8            work_d4(2,nnpn_max)
      real*8            work_d5(2,nnpn_max)
      real*8            work_d6(2,nnpn_max)
      real*8            work_d7(2,nnpn_max)
      real              work_r1(nnpn_max,2)
      real              work_r2(nnpn_max)
      integer           work_i1(nnpn_max)
      integer           work_i2(2,nnpn_max)
      integer           work_i3(nnpn_max)
      logical*1         lt_work(nt_max)
      logical*1         ln_work(np_max)
      logical*1         nnwrite
c
      data              eps/1.E-10/
c
       common /nnswitches/lud,nnwrite
c
c					Set debug mode for nn routines
      nnwrite = .false.
      lud=0
c
      rescale = .false.
      clockwise=.false.
      nnext=.true.

c  Interpolation method - options are:
c   0=NN-Watson
c   1=NN-recursive
c   2=NN-recursive+derivatives
c   3=Linear interpolation
c   4=NN-closed formula
c   5=NN-closed formula + df(2)
c   6=NN-closed formula + df(2) + ddf(3)

      int_method=0
c
c  set up delaunay calculation mode
c    (0=qhull,-1=delaun+x-sort,-2=delaun,>0=read from file)
     
      dmode = -2
c
c  set up nn mode (0=calc delaun+setup)
c
      nmode = 0

      xmin=0
      xmax=im
      ymin=0
      ymax=jm
      favg=0
      do i=1,np
       points(1,i)=x(i)
       points(2,i)=y(i)
       fvals(i)=z(i)
       favg=favg+fvals(i)
      enddo
      favg=favg/np
c
c  perturb nodes - djs
c
c      delx=1.e-6*(xmax-xmin)
c      dely=1.e-6*(ymax-ymin)
c  max perturbation = 0.05 * grid size
      delx=0.1
      dely=0.1


      do i=1,np
c       write(6,'(i5,3e15.5)') i,points(1,i),points(2,i),fvals(i)
c       r=rand()  !f77
c       r=rrand() !lf95
       call random(r)
       points(1,i)=points(1,i)+delx*(0.5-r)
c       r=rand()
c      r=rrand()
       call random(r)
       points(2,i)=points(2,i)+dely*(0.5-r)
c       write(6,'(i5,3e15.5)') i,points(1,i),points(2,i),fvals(i)
      enddo


c
c  add distant points - djs
c
      points(1,np+1)=xmin-1.e3*(xmax-xmin)
      points(2,np+1)=ymin-1.e3*(ymax-ymin)
      points(1,np+2)=xmin-1.e3*(xmax-xmin)
      points(2,np+2)=ymax+1.e3*(ymax-ymin)
      points(1,np+3)=xmax+1.e3*(xmax-xmin)
      points(2,np+3)=ymax+1.e3*(ymax-ymin)
      points(1,np+4)=xmax+1.e3*(xmax-xmin)
      points(2,np+4)=ymin-1.e3*(ymax-ymin)
      fvals(np+1)=favg
      fvals(np+2)=favg
      fvals(np+3)=favg
      fvals(np+4)=favg
      np1=np+4


c					Perform setup of nn interpolation
c					(calculate Delaunay triangles,
c					 circumcentres, build neighbour
c					 and hulltriangles)
c
c					2-D setup
      call nn2dsetup
     &     (np1,nt_max,nh_max,np_max,nnpn_max,nmax,
     &     points,dmode,nmode,clockwise,fvals,nt,vertices,
     &     centres,neighbour,nh,hulltriangles,nohalt_hull,
     &     loc,nnn,nnlist,ntrilist,eps,nv_max,vis_tlist,
     &     vis_elist,add_tlist,lt_work,ln_work)

c                                       check consistency of
c                                       neighbour matrix
c					(for debug purposes)

c      write(6,*) 'NT    =',nt


      call check_neighbour(neighbour,nd,vertices,nt,consistent)

c      if(consistent)write(*,*)' neighbour matrix consistent'
c
c                                       write out vertices and neighbours
c					(for debug purposes)
c
      k = 0
c      open(21,file='list',status='unknown')
c      do 11 i=1,nt
c         write(21,*)i,' :',(vertices(j,i),j=1,3),
c     &              ' n:',(neighbour(j,i),j=1,3)
c 11   continue

c					Perform natural neighbour interpolation
      do j=1,jm
       do i=1,im
        zgrid(i,j)=0.
        if(d(i,j).ne.0) then
         xd=i-0.5
         yd=j-0.5
c					Perform interpolation at (xd,yd)

         call nn2D
     &     (xd,yd,points,vertices,neighbour,nnext,int_method,
     &      centres,hulltriangles,nh,loc,fvals,nnpn_max,work_r1,
     &      work_r2,work_d1,work_i1,work_i2,work_i3,work_d2,
     &      work_d3,work_d4,work_d5,work_d6,work_d7,lt_work,
     &      ln_work,out,f,df,ddf,v)

c					write out results
         zgrid(i,j)=f
         
	 if(out)write(*,*) "Point is outside convex hull",xd,yd
        endif
       enddo
      enddo
      return
      end
      
C ----------------------------------------------------------------------
      SUBROUTINE DBJDAT(IYEARS,IMON,IDAY,IJULS,DBERR)
C
C DBJDAT (IYEARS,IMON,IDAY,IJULS,DBERR)
C
C  Purpose: This subroutine converts a Gregorian calendar date to the
C  corresponding Julian day number 'IJUL'.  The Julian day number 'IJUL'
C  is computed from the given day 'IDAY', month 'IMON', and year 'IYEAR'.
C  without using tables.  The procedure is valid for any valid Gregorian
C  calendar date.  Leap year is defined to be any year divisible by 4
C  except centenary years not divisible by 400.  The routine is based on
C  Algorithm 199 in The Collected Algorithms of the ACM, as presented by
C  R. G. Tantzen in 1963.
C
C Arguments:
C On input:
C
C   IYEARS - Year for which Julian day is to be calculated.
C
C   IMON - Month for which Julian day is to be calculated.
C
C   IDAY - Day for which Julian day is to be calculated.
C
C ON OUTPUT:
C
C   IJULS - JULIAN DAY CALCULATED FROM IYEAR, IMON AND IDAY.
C
C   DBERR - LOGICAL variable returned with the value of .FALSE. if no errors
C     are detected and the value of .TRUE. if an error is detected.
C     Calling program must declare this variable to be of type LOGICAL and
C     is responsible for checking this variable to see if an error has occured.
C
C SUBROUTINES CALLED:
C
C   NONE
C
C I/O UNITS USED:
C
C   6 - USED FOR ERROR REPORT.
C
C History: Written by Edward W. Lynn, GLERL, April 1985.
C
      LOGICAL DBERR
C
C  SET DBERR TO FALSE INITIALLY
      DBERR = .FALSE.
C
C  CONFIRM THAT MONTH IS LEGAL.
      IF(IYEARS .LT. 0) THEN
       WRITE(6,'(A,I6)') ' DBJDAT: ILLEGAL INPUT PARAMETER YEAR ',IYEARS
       DBERR = .TRUE.
      END IF
      IF(DBERR) RETURN
C
      ID = IDAY
      IM = IMON
      IY = IYEARS
C
      IDF = 1
      IMF = 1 + 9
      IYF = IYEARS - 1
C
      IF(IM .GT. 2) THEN
       IM = IM - 3
      ELSE
       IM = IM + 9
       IY = IY-1
      END IF
C
      IJULF = (146097*(IYF/100))/4 + (1461*(IYF-100*(IYF/100)))/4 +
     1       (153*IMF+2)/5 + IDF + 1721119
      IJULL = (146097*(IY/100))/4 + (1461*(IY-100*(IY/100)))/4 +
     1       (153*IM+2)/5 + ID + 1721119
      IJULS = IJULL - IJULF + 1
      RETURN
      END

C Subroutine for generating random numbers having a uniform
c distribution, by the mixed multiplicative congruential method.
c Applied Numerical methods, Carnahan et al., pg 549
      subroutine random(z)
      save
      data i/1/
      integer a,x
      if (i.eq.0) go to 1
         i=0
         m=2**20
         fm=m
         x=566387
         a=2**10+3

1     x=mod(a*x,m)
      fx=x
      z=fx/fm
      return
      end



