
module loess 

    use coord_constants
    use interp1D 
    use index 

    implicit none 

    ! Define the default missing value 
    real(sp), parameter :: LOESS_MISSING_VALUE = -9999.0 

    interface loess_smooth 
        module procedure loess_smooth_dble
        module procedure loess_smooth_float
    end interface 

    private
    public :: loess_smooth
!     public :: loess_smooth_dble

contains 

    function loess_smooth_dble(x,y,L,missing_value) result(yfit)

        implicit none 

        real(dp) :: x(:), y(:)
        real(dp) :: yfit(size(x)) 
        real(dp) :: L, missing  
        real(dp), optional :: missing_value 

        ! Assign missing value
        missing = dble(LOESS_MISSING_VALUE)
        if (present(missing_value)) missing = missing_value 

        ! Calculate loess smoothing (pass reals because `lowess` calculates in reals)
        yfit = loess_with_interp(real(x),real(y),real(L),real(missing))

        return 

    end function loess_smooth_dble


    function loess_smooth_float(x,y,L,missing_value) result(yfit)

        implicit none 

        real(sp) :: x(:), y(:)
        real(sp) :: yfit(size(x)) 
        real(sp) :: L, missing  
        real(sp), optional :: missing_value 

        ! Assign missing value
        missing = LOESS_MISSING_VALUE
        if (present(missing_value)) missing = missing_value 

        ! Calculate loess smoothing 
        yfit = loess_with_interp(x,y,L,missing)

        return 

    end function loess_smooth_float

    function loess_with_interp(x,y,L,missing_value) result(yfit)
        ! Internal interface to the lowess calculations
        ! This function ensures that unevenly spaced, randomly
        ! distributed x-values with missing y-values are handled properly.
        ! All input values are linearly interpolated to a common x-scale
        ! with a constant dx value. The smoothing window half-width, L, is
        ! given in units of x, and is used to generated the lowess smoothing
        ! parameter f - the fraction of points used to construct the local 
        ! polynomial.
        ! Missing values are omitted from the calculations.
        ! The output vector yfit has the same size as the input data. 

        implicit none 

        real(sp) :: x(:), y(:) 
        real(sp) :: yfit(size(x))
        real(sp) :: L, missing_value 

        real(sp), allocatable :: xtmp(:), ytmp(:), ys(:), rw(:), res(:)
        integer, allocatable :: ind(:)
        integer  :: np, np1, k
        real(sp) :: dx, x0, x1, frac

        ! Determine how many points we need 
        np  = size(x) 
        dx  = minval(x(2:np)-x(1:(np-1)))
        x0  = minval(x)
        x1  = maxval(x)
        np1 = ceiling( (x1-x0) / dx ) + 1
        
        ! Allocate temporary vectors
        allocate(xtmp(np1),ytmp(np1))

        ! Fill the x-axis 
        do k = 1, np1 
            xtmp(k) = x0 + (k-1)*dx 
        end do  

        ! Determine which indices are not missing values, to use for interpolation
        call which(y .ne. missing_value,ind)

        if (size(ind) .le. int(L*2.0+1.0)) then 
!             write(*,*) "loess_smooth:: error: not enough points available for smoothing."
!             write(*,*) "Missing value: ", missing_value 
!             write(*,*) "Missing points: ", count(y.eq.missing_value), " / ", size(y) 
!             write(*,*) "Smoothing window half-width L = ", L 
!             stop 

            yfit = missing_value 
        
        else 
            ! Perform lowess smoothing 

            ! Fill the y-axis by linear interpolation 
            ytmp = real(interp_linear(dble(x(ind)),dble(y(ind)),xout=dble(xtmp)))

            ! Determine the frac parameter for loess smoothing
            ! ie, the fraction of total points to use, where L is the half-width
            ! of the desired smoothing window in the units of x. 
            frac = ((L*2.0+1.0)/dx) / size(ytmp) 

            write(*,*) "L, frac = ", L, frac 

            ! Calculate loess curve 
            allocate(ys(np1),rw(np1),res(np1))
            call lowess(xtmp,ytmp,n=np1,f=frac,nsteps=2,delta=0.0,ys=ys,rw=rw,res=res) 

            ! Fit the smooth values back to the original time series 
            yfit = real(interp_linear(dble(xtmp),dble(ys),xout=dble(x)))

        end if 
        
        return

    end function loess_with_interp 

!                                                                       
!              COMPUTER PROGRAMS FOR LOCALLY WEIGHTED REGRESSION        
!                                                                       
!             This package consists  of  two  FORTRAN  programs  for    
!        smoothing    scatterplots   by   robust   locally   weighted   
!        regression, or lowess.   The  principal  routine  is  LOWESS   
!        which   computes   the  smoothed  values  using  the  method   
!        described in The Elements of Graphing Data, by William S.      
!        Cleveland    (Wadsworth,    555 Morego   Street,   Monterey,   
!        California 93940).                                             
!                                                                       
!             LOWESS calls a support routine, LOWEST, the code for      
!        which is included. LOWESS also calls a routine  SORT,  which   
!        the user must provide.                                         
!                                                                       
!             To reduce the computations, LOWESS  requires  that  the   
!        arrays  X  and  Y,  which  are  the  horizontal and vertical   
!        coordinates, respectively, of the scatterplot, be such  that   
!        X  is  sorted  from  smallest  to  largest.   The  user must   
!        therefore use another sort routine which will sort X  and  Y   
!        according  to X.                                               
!             To summarize the scatterplot, YS,  the  fitted  values,   
!        should  be  plotted  against X.   No  graphics  routines are   
!        available in the package and must be supplied by the user.     
!                                                                       
!                                                                       
!        Calling sequence                                               
!                                                                       
!        CALL LOWESS(X,Y,N,F,NSTEPS,DELTA,YS,RW,RES)                    
!                                                                       
!        Purpose                                                        
!                                                                       
!        LOWESS computes the smooth of a scatterplot of Y  against  X   
!        using  robust  locally  weighted regression.  Fitted values,   
!        YS, are computed at each of the  values  of  the  horizontal   
!        axis in X.                                                     
!                                                                       
!        Argument description                                           
!                                                                       
!              X = Input; abscissas of the points on the                
!                  scatterplot; the values in X must be ordered         
!                  from smallest to largest.                            
!              Y = Input; ordinates of the points on the                
!                  scatterplot.                                         
!              N = Input; dimension of X,Y,YS,RW, and RES.              
!              F = Input; specifies the amount of smoothing; F is       
!                  the fraction of points used to compute each          
!                  fitted value; as F increases the smoothed values     
!                  become smoother; choosing F in the range .2 to       
!                  .8 usually results in a good fit; if you have no     
!                  idea which value to use, try F = .5.                 
!         NSTEPS = Input; the number of iterations in the robust        
!                  fit; if NSTEPS = 0, the nonrobust fit is             
!                  returned; setting NSTEPS equal to 2 should serve     
!                  most purposes.                                       
!          DELTA = input; nonnegative parameter which may be used       
!                  to save computations; if N is less than 100, set     
!                  DELTA equal to 0.0; if N is greater than 100 you     
!                  should find out how DELTA works by reading the       
!                  additional instructions section.                     
!             YS = Output; fitted values; YS(I) is the fitted value     
!                  at X(I); to summarize the scatterplot, YS(I)         
!                  should be plotted against X(I).                      
!             RW = Output; robustness weights; RW(I) is the weight      
!                  given to the point (X(I),Y(I)); if NSTEPS = 0,       
!                  RW is not used.                                      
!            RES = Output; residuals; RES(I) = Y(I)-YS(I).              
!                                                                       
!                                                                       
!        Other programs called                                          
!                                                                       
!               LOWEST                                                  
!               SORT                                                    
!                                                                       
!        Additional instructions                                        
!                                                                       
!        DELTA can be used to save computations.   Very  roughly  the   
!        algorithm  is  this:   on the initial fit and on each of the   
!        NSTEPS iterations locally weighted regression fitted  values   
!        are computed at points in X which are spaced, roughly, DELTA   
!        apart; then the fitted values at the  remaining  points  are   
!        computed  using  linear  interpolation.   The  first locally   
!        weighted regression (l.w.r.) computation is carried  out  at   
!        X(1)  and  the  last  is  carried  out at X(N).  Suppose the   
!        l.w.r. computation is carried out at  X(I).   If  X(I+1)  is   
!        greater  than  or  equal  to  X(I)+DELTA,  the  next  l.w.r.   
!        computation is carried out at X(I+1).   If  X(I+1)  is  less   
!        than X(I)+DELTA, the next l.w.r.  computation is carried out   
!        at the largest X(J) which is greater than or equal  to  X(I)   
!        but  is not greater than X(I)+DELTA.  Then the fitted values   
!        for X(K) between X(I)  and  X(J),  if  there  are  any,  are   
!        computed  by  linear  interpolation  of the fitted values at   
!        X(I) and X(J).  If N is less than 100 then DELTA can be  set   
!        to  0.0  since  the  computation time will not be too great.   
!        For larger N it is typically not necessary to carry out  the   
!        l.w.r.  computation for all points, so that much computation   
!        time can be saved by taking DELTA to be  greater  than  0.0.   
!        If  DELTA =  Range  (X)/k  then,  if  the  values  in X were   
!        uniformly  scattered  over  the  range,  the   full   l.w.r.   
!        computation  would be carried out at approximately k points.   
!        Taking k to be 50 often works well.                            
!                                                                       
!        Method                                                         
!                                                                       
!        The fitted values are computed by using the nearest neighbor   
!        routine  and  robust locally weighted regression of degree 1   
!        with the tricube weight function.  A few additional features   
!        have  been  added.  Suppose r is FN truncated to an integer.   
!        Let  h  be  the  distance  to  the  r-th  nearest   neighbor   
!        from X(I).   All  points within h of X(I) are used.  Thus if   
!        the r-th nearest neighbor is exactly the  same  distance  as   
!        other  points,  more  than r points can possibly be used for   
!        the smooth at  X(I).   There  are  two  cases  where  robust   
!        locally  weighted regression of degree 0 is actually used at   
!        X(I).  One case occurs when  h  is  0.0.   The  second  case   
!        occurs  when  the  weighted  standard error of the X(I) with   
!        respect to the weights w(j) is  less  than  .001  times  the   
!        range  of the X(I), where w(j) is the weight assigned to the   
!        j-th point of X (the tricube  weight  times  the  robustness   
!        weight)  divided by the sum of all of the weights.  Finally,   
!        if the w(j) are all zero for the smooth at X(I), the  fitted   
!        value is taken to be Y(I).                                     
!                                                                       
                                                                        
      subroutine lowess(x, y, n, f, nsteps, delta, ys, rw, res) 
      integer n 
      integer nsteps 
      real x(n), y(n), f, delta, ys(n), rw(n) 
      real res(n) 
      integer nright, min0, max0, i, j, ifix 
      integer iter, last, m1, m2, ns, nleft 
      real abs, cut, cmad, r, d1, d2 
      real c1, c9, alpha, denom, float 
      logical ok 
      if (n .ge. 2) goto 1 
         ys(1) = y(1) 
         return 
! at least two, at most n points                                        
    1 ns = max0(min0(ifix(f*float(n)), n), 2) 
      iter = 1 
         goto  3 
    2    iter = iter+1 
    3    if (iter .gt. nsteps+1) goto  22 
! robustness iterations                                                 
         nleft = 1 
         nright = ns 
! index of prev estimated point                                         
         last = 0 
! index of current point                                                
         i = 1 
    4       if (nright .ge. n) goto  5 
! move nleft, nright to right if radius decreases                       
               d1 = x(i)-x(nleft) 
! if d1<=d2 with x(nright+1)==x(nright), lowest fixes                   
               d2 = x(nright+1)-x(i) 
               if (d1 .le. d2) goto  5 
! radius will not decrease by move right                                
               nleft = nleft+1 
               nright = nright+1 
               goto  4 
! fitted value at x(i)                                                  
    5       call lowest(x, y, n, x(i), ys(i), nleft, nright, res, iter  &
     &     .gt. 1, rw, ok)                                              
            if (.not. ok) ys(i) = y(i) 
! all weights zero - copy over value (all rw==0)                        
            if (last .ge. i-1) goto 9 
               denom = x(i)-x(last) 
! skipped points -- interpolate                                         
! non-zero - proof?                                                     
               j = last+1 
                  goto  7 
    6             j = j+1 
    7             if (j .ge. i) goto  8 
                  alpha = (x(j)-x(last))/denom 
                  ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last) 
                  goto  6 
    8          continue 
! last point actually estimated                                         
    9       last = i 
! x coord of close points                                               
            cut = x(last)+delta 
            i = last+1 
               goto  11 
   10          i = i+1 
   11          if (i .gt. n) goto  13 
! find close points                                                     
               if (x(i) .gt. cut) goto  13 
! i one beyond last pt within cut                                       
               if (x(i) .ne. x(last)) goto 12 
                  ys(i) = ys(last) 
! exact match in x                                                      
                  last = i 
   12          continue 
               goto  10 
! back 1 point so interpolation within delta, but always go forward     
   13       i = max0(last+1, i-1) 
   14       if (last .lt. n) goto  4 
! residuals                                                             
         do  15 i = 1, n 
            res(i) = y(i)-ys(i) 
   15       continue 
         if (iter .gt. nsteps) goto  22 
! compute robustness weights except last time                           
         do  16 i = 1, n 
            rw(i) = abs(res(i)) 
   16       continue 
         call sort(n,rw) 
         m1 = n/2+1 
         m2 = n-m1+1 
! 6 median abs resid                                                    
         cmad = 3.0*(rw(m1)+rw(m2)) 
         c9 = .999*cmad 
         c1 = .001*cmad 
         do  21 i = 1, n 
            r = abs(res(i)) 
            if (r .gt. c1) goto 17 
               rw(i) = 1. 
! near 0, avoid underflow                                               
               goto  20 
   17          if (r .le. c9) goto 18 
                  rw(i) = 0. 
! near 1, avoid underflow                                               
                  goto  19 
   18             rw(i) = (1.0-(r/cmad)**2)**2 
   19       continue 
   20       continue 
   21       continue 
         goto  2 
   22 return 
      end subroutine lowess   

      subroutine lowest(x, y, n, xs, ys, nleft, nright, w, userw        &
     &, rw, ok)                                                         
      integer n 
      integer nleft, nright 
      real x(n), y(n), xs, ys, w(n), rw(n) 
      logical userw, ok 
      integer nrt, j 
      real abs, a, b, c, h, r 
      real h1, sqrt, h9, amax1, range 
      range = x(n)-x(1) 
      h = amax1(xs-x(nleft), x(nright)-xs) 
      h9 = .999*h 
      h1 = .001*h 
! sum of weights                                                        
      a = 0.0 
      j = nleft 
         goto  2 
    1    j = j+1 
    2    if (j .gt. n) goto  7 
! compute weights (pick up all ties on right)                           
         w(j) = 0. 
         r = abs(x(j)-xs) 
         if (r .gt. h9) goto 5 
            if (r .le. h1) goto 3 
               w(j) = (1.0-(r/h)**3)**3 
! small enough for non-zero weight                                      
               goto  4 
    3          w(j) = 1. 
    4       if (userw) w(j) = rw(j)*w(j) 
            a = a+w(j) 
            goto  6 
    5       if (x(j) .gt. xs) goto  7 
! get out at first zero wt on right                                     
    6    continue 
         goto  1 
! rightmost pt (may be greater than nright because of ties)             
    7 nrt = j-1 
      if (a .gt. 0.0) goto 8 
         ok = .false. 
         goto  16 
    8    ok = .true. 
! weighted least squares                                                
         do  9 j = nleft, nrt 
! make sum of w(j) == 1                                                 
            w(j) = w(j)/a 
    9       continue 
         if (h .le. 0.) goto 14 
            a = 0.0 
! use linear fit                                                        
            do  10 j = nleft, nrt 
! weighted center of x values                                           
               a = a+w(j)*x(j) 
   10          continue 
            b = xs-a 
            c = 0.0 
            do  11 j = nleft, nrt 
               c = c+w(j)*(x(j)-a)**2 
   11          continue 
            if (sqrt(c) .le. .001*range) goto 13 
               b = b/c 
! points are spread out enough to compute slope                         
               do  12 j = nleft, nrt 
                  w(j) = w(j)*(b*(x(j)-a)+1.0) 
   12             continue 
   13       continue 
   14    ys = 0.0 
         do  15 j = nleft, nrt 
            ys = ys+w(j)*y(j) 
   15       continue 
   16 return 
      end subroutine lowest                                           
                                                                        
      subroutine sort(n,arr) 
      INTEGER n,M,NSTACK 
      REAL arr(n) 
      PARAMETER (M=7,NSTACK=50) 
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK) 
      REAL a,temp 
      jstack=0 
      l=1 
      ir=n 
    1 if(ir-l.lt.M)then 
        do 12 j=l+1,ir 
          a=arr(j) 
          do 11 i=j-1,1,-1 
            if(arr(i).le.a)goto 2 
            arr(i+1)=arr(i) 
   11     continue 
          i=0 
    2     arr(i+1)=a 
   12   continue 
        if(jstack.eq.0)return 
        ir=istack(jstack) 
        l=istack(jstack-1) 
        jstack=jstack-2 
      else 
        k=(l+ir)/2 
        temp=arr(k) 
        arr(k)=arr(l+1) 
        arr(l+1)=temp 
        if(arr(l+1).gt.arr(ir))then 
          temp=arr(l+1) 
          arr(l+1)=arr(ir) 
          arr(ir)=temp 
        endif 
        if(arr(l).gt.arr(ir))then 
          temp=arr(l) 
          arr(l)=arr(ir) 
          arr(ir)=temp 
        endif 
        if(arr(l+1).gt.arr(l))then 
          temp=arr(l+1) 
          arr(l+1)=arr(l) 
          arr(l)=temp 
        endif 
        i=l+1 
        j=ir 
        a=arr(l) 
    3   continue 
          i=i+1 
        if(arr(i).lt.a)goto 3 
    4   continue 
          j=j-1 
        if(arr(j).gt.a)goto 4 
        if(j.lt.i)goto 5 
        temp=arr(i) 
        arr(i)=arr(j) 
        arr(j)=temp 
        goto 3 
    5   arr(l)=arr(j) 
        arr(j)=a 
        jstack=jstack+2 
        if(jstack.gt.NSTACK) stop 'NSTACK too small in sort' 
        if(ir-i+1.ge.j-l)then 
          istack(jstack)=ir 
          istack(jstack-1)=i 
          ir=j-1 
        else 
          istack(jstack)=j-1 
          istack(jstack-1)=l 
          l=i 
        endif 
      endif 
      goto 1 
      end subroutine sort   

end module loess

