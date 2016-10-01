      SUBROUTINE RGRID(LUN, D, IDIM, JDIM)

C PURPOSE:
C                    TO READ A BATHYMETRIC GRID DATA FILE
C                    AND RETURN GRID PARAMETERS AND DEPTHS.
C ARGUMENTS:
C ON INPUT:
C                    LUN - LOGICAL UNIT NUMBER OF BATHYMETRIC DATA FILE
C                    IDIM - FIRST DIMENSION OF ARRAY D IN
C                           DIMENSION STATEMENT OF CALLING PROGRAM
C                    JDIM - SECOND DIMENSION OF ARRAY D IN
C                           DIMENSION STATEMENT OF CALLING PROGRAM
C
C                      FORMAT OF BATHYMETRIC DATA FILE
C       ------------------------------------------------------------------
C                 FIELD                            FORTRAN FORMAT  COLUMNS
C       RECORD 1: LAKE NAME                            40A1          1-40
C       RECORD 2: FIRST (I) DIMENSION OF DEPTH ARRAY   I5            1-5
C                 SECOND (J) DIMENSION OF DEPTH ARRAY  I5            6-10
C                 BASE LATITUDE                        F12.7        11-22
C                 BASE LONGITUDE                       F12.7        23-34
C                 GRID SIZE                            F5.0         35-39
C                 MAXIMUM DEPTH                        F5.0         40-44
C                 MINIMUM DEPTH                        F5.0         45-49
C                 BASE ROTATION (CCW FROM E-W)         F6.2         50-55
C                 I DISPLACEMENT                       I5           56-60
C                 J DISPLACEMENT                       I5           61-65
C                 ROTATION FROM BASE (CCW)             F6.2         66-71
C                 POWER OF TEN TO CONVERT DEPTHS TO 
C                                       METERS         I3           72-74
C                 MAP PROJECTION INDICATOR
C                  0=APPROXIMATE POLYCONIC
C                  1=LAMBERT CONFORMAL CONIC           I2           75-76
C       RECORDS 3-6 (FOR APPROXIMATE POLYCONIC PROJECTION):
C                 MAP PROJECTION COORDINATE CONVERSION
C                  COEFFICIENTS                       4E15.6         1-60
C       RECORD 3 (FOR LAMBERT CONFORMAL CONIC PROJECTION):
C                 CENTRAL MERIDIAN OF PROJECTION       E15.6         1-15
C                 SOUTHERNMOST STANDARD PARALLEL       E15.6        16-30
C                 NORTHERNMOST STANDARD PARALLEL       E15.6        31-45
C       FOLLOWING RECORDS:
C                 DEPTHS IN ASCENDING I, ASCENDING J
C                  SEQUENCE, 19 TO A RECORD           19F4.0         1-76
C       ------------------------------------------------------------------
C
C ON OUTPUT:
C                     D - DEPTH ARRAY. ZERO FOR LAND, AVERAGE DEPTH
C                         OF GRID BOX IN METERS FOR WATER.
C                     RPARM - ARRAY CONTAINING REAL-VALUED BATHYMETRIC
C                            GRID PARAMETERS AS FOLLOWS:
C                        1.  BASE LATITUDE
C                        2.  BASE LONGITUDE
C                        3.  GRID SIZE (M)
C                        4.  MAXIMUM DEPTH (M)
C                        5.  MINIMUM DEPTH (M)
C                        6.  BASE ROTATION (COUNTERCLOCKWISE FROM E-W)
C                        7.  ROTATION FROM BASE (COUNTERCLOCKWISE)
C
C                      FOR IPARM(45)=0 (APPROXIMATE POLYCONIC PROJECTION):
C                        8-11. GEOGRAPHIC-TO-MAP COORDINATE CONVERSION
C                             COEFFICIENTS FOR X
C                        12-15.  GEOGRAPHIC-TO-MAP COORDINATE CONVERSION
C                                COEFFICIENTS FOR Y
C                        16-19.  MAP-TO-GEOGRAPHIC COORDINATE CONVERSION
C                                COEFFICIENTS FOR LONGITUDE
C                        20-23.  MAP-TO-GEOGRAPHIC COORDINATE CONVERSION
C                                COEFFICIENTS FOR LATITUDE
C
C                      FOR IPARM(45)=1 (LAMBERT CONFORMAL CONIC PROJECTION):
C                        8.  CENTRAL MERIDIAN OF PROJECTION (RADIANS)
C                        9.  SOUTHERNMOST STANDARD PARALLEL (RADIANS)
C                        10. NORTHERNMOST STANDARD PARALLEL (RADIANS)
C                        11. LOGARITHMIC COEFFICIENT FOR TRANSFORMATIONS
C                        12. DISTANCE SCALING FACTOR FOR TRANSFORMATIONS
C                        13. X DISPLACEMENT OF BATHYMETRIC GRID ORIGIN
C                             FROM MAP PROJECTION ORIGIN
C                        14. Y DISPLACEMENT OF BATHYMETRIC GRID ORIGIN
C                             FROM MAP PROJECTION ORIGIN
C
C                    IPARM - ARRAY CONTAINING INTEGER-VALUED BATHYMETRIC
C                            GRID PARAMETERS AS FOLLOWS:
C                        1.  NUMBER OF GRID BOXES IN X DIRECTION
C                        2.  NUMBER OF GRID BOXES IN Y DIRECION
C                        3.  I DISPLACEMENT - THE NUMBER OF NEW GRID
C                            SQUARES IN THE X-DIRECTION FROM THE NEW
C                            GRID ORIGIN TO THE OLD GRID ORIGIN
C                            (USED ONLY FOR IPARM(45)=0)
C                        4.  J DISPLACEMENT - THE NUMBER OF NEW GRID
C                            SQUARES IN THE Y-DIRECTION FROM THE NEW
C                            GRID ORIGIN TO THE OLD GRID ORIGIN
C                            (USED ONLY FOR IPARM(45)=0)
C                        5-44.  LAKE NAME (40A1)
C                        45. MAP PROJECTION USED FOR BATHYMETRIC GRID:
C                            0=APPROXIMATE POLYCONIC (GREAT LAKES GRIDS)
C                            1=LAMBERT CONFORMAL CONIC
C                        46-54. RESERVED FOR FUTURE USE
C
C                     NOTE: IF GRID IS TOO LARGE FOR DIMENSIONS OF D,
C                          THE IPARM ARRAY IS SET TO ZERO
C
C COMMON BLOCK:
C                    /GPARM/RPARM(23),IPARM(54)
C
      DIMENSION D(IDIM,JDIM)
      COMMON /GPARM/ RPARM(23), IPARM(54)
      DATA DTR /0.017453293/
      READ (LUN,'(40A1,I1)') (IPARM(I),I=5,44)
      READ (LUN,'(2I5, 2F12.7, 3F5.0, F6.2, 2I5, F6.2, I3, I2)')
     1  IPARM(1), IPARM(2), (RPARM(I),I=1,6), IPARM(3), IPARM(4),
     2  RPARM(7), IDEXP, IPARM(45)
      IM = IPARM(1)
      JM = IPARM(2)
      IF (IPARM(1) .GT. IDIM .OR. IPARM(2) .GT. JDIM) GO TO 10

C  READ CONSTANTS FOR APPROXIMATE POLYCONIC MAP PROJECTION

      IF(IPARM(45).EQ.0) READ (LUN,'(4E15.6)') (RPARM(I),I=8,23)

C  READ CONSTANTS FOR LAMBERT CONFORMAL CONIC MAP PROJECTION

      IF(IPARM(45).EQ.1) READ (LUN,'(3E15.6)') (RPARM(I),I=8,10)

C  READ DEPTHS

      READ (LUN,'(19F4.0, 4X)') ((D(I,J),I=1,IM),J=1,JM)
C
C  ADJUST DEPTHS
C
      DFAC=10.**IDEXP
      RPARM(4)=RPARM(4)*DFAC
      RPARM(5)=RPARM(5)*DFAC
      DO 5 I=1,IM
      DO 5 J=1,JM
5     D(I,J)=D(I,J)*DFAC
C
C  FOR LAMBERT PROJECTION COMPUTE REQUIRED CONSTANTS
C
      IF(IPARM(45).EQ.1) THEN
       A45=ATAN(1.)
       RPARM(8)=RPARM(8)*DTR
       RPARM(9)=RPARM(9)*DTR
       RPARM(10)=RPARM(10)*DTR
       ALON0=RPARM(8)
       A1=RPARM(9)
       A2=RPARM(10)
       RPARM(11)=(LOG(COS(A1))-LOG(COS(A2)))/
     1    (LOG(TAN(A45-A1/2.))-LOG(TAN(A45-A2/2.)))

C  SET SCALE FACTOR FOR N-S DISTANCE FROM A1 TO A2

       AEXP=RPARM(11)
       Y1=(TAN(A45-A1/2.))**AEXP
       Y2=(TAN(A45-A2/2.))**AEXP
       RPARM(12)=6378140.*(A2-A1)/(Y1-Y2)
       RPARM(13)=0.
       RPARM(14)=0.
       DX=XDIST(RPARM(1),RPARM(2))
       DY=YDIST(RPARM(1),RPARM(2))
       RPARM(13)=DX
       RPARM(14)=DY
      ENDIF
      
      RETURN

C  COME HERE ON ERROR

   10 DO 20 I = 1, 54
   20 IPARM(I) = 0
      WRITE (6,50)
   50 FORMAT (' BATHYMETRIC GRID TOO LARGE - INCREASE DIMENSIONS OF',
     1       ' NDEPTH AND DEPTH IN MAIN PROGRAM')
   70 RETURN
      END
