!-----------------------------------------------------------------------
      FUNCTION JULIAN(yr,month,day,hour)
      real*8 jday,yr,month,day,hour
      real*8 Y, m ,julian

      Y = yr
      if (yr.lt.100. .and. yr.gt.50.) then
        Y = yr+1900.
      elseif (yr.lt.100. .and. yr.le.50.) then
        Y = yr+2000.
      endif
            
      if (month.le. 2.) then
        Y = Y-1
        m = month +12.
      else
        Y = Y
        m = month
      endif
      
      JULIAN = aint(365.25*Y) + aint(30.6001*(m+1))                        
     &       + day + hour/24. + 1720981.50
!      write(*,'(''Julian ='',f12.4,4f8.2)') julian, yr,month,day,hour
      END FUNCTION JULIAN
      
      SUBROUTINE JULIANSUB(jday,yr,month,day,hour)
      real*8 jday,yr,month,day,hour,julian      
      jday = JULIAN(yr,month,day,hour)
      END
            
      SUBROUTINE GREGORIAN(jday,yr,month,day,hour)
      real*8 jday,yr,month,day,hour,dayoweek,week
       
      a = aint(jday+.5)
      b = a+1537
      c = aint((b-122.1)/365.25)
      d = aint(365.25*c)
      e = aint((b-d)/30.6001)
      day = b-d-aint(30.6001*e)+mod(jday+0.5d00,1.0d00)
      month = e-1-12*aint(e/14)
      yr = c-4715-aint((7+month)/10)
      hour = mod(jday+0.5,1.0d00)*24.
      day = aint(day)
      
      dayoweek = mod(aint(jday+0.5d00),7.d00)
      week = aint((jday-2444244.5)/7)
!  2000 leap year crashed.  
      if(month.eq.2. .and. day .eq. 31. ) then
       yr=yr-1.
       day = 29.
      endif
      
      END

      SUBROUTINE daynumber(IYEARS,IMON,IDAY,IJULS,DBERR)
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
      real*8 IYEARS
      real*8 IMON
      real*8 IDAY
      real*8 IJULS
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


