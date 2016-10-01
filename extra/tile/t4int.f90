

SUBROUTINE NNBRHG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART, &
          NPTSIN,L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT, &
          NINDEX,AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE,         &
          VAL,Z0,ZX0,ZY0,ZH,ZXH,ZYH) 

    !   Initial data values: PT, NPTS, VAL(NPTS),GRAD(2,NPTS)
    !   CN = constraints??
    !   Working arrays: SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX), 
    !   Working dimension: KNBMAX
    !   Calculates the C0 natural neighbour interpolant and its gradient, 
    !   and also the C1 natural neighbour interpolant and its gradient with 
    !   all data site gradients forced to zero, at the point (XPT,YPT). 
    !   This routine thus produces the same effect as a call to NNBR1G 
    !   with zeroes entered into GRAD, but the array GRAD does not need 
    !   to be passed to it and the calculation is more efficient.  The 
    !   C1 interpolant value is returned as ZH, its gradient as (ZXH,ZYH). 
    !   Other details are as for NNBR1G. 

        integer :: JCNS, NPTS, LTOP, KNBMAX
        real(4) :: CN(3,JCNS), JADDR(JCNS)
        real(4) :: PT(2,NPTS), VAL(NPTS) 
        integer :: NADDR(NPTS), L(LTOP)
        real(4) :: SBAREA(KNBMAX),DELSBA(2,KNBMAX),PTOFF(2,KNBMAX)
        integer :: NFREE, NSTART, NPTSIN, NINDEX, LBASE, IGARB, KNB
        real(4) :: EPSCN, EPSPT, AREA

        integer :: IFLAG, ICASE, LLO, LHI, LPTR, K, N, LL 
        real(4) :: Z0, ZX0, ZY0, ZH, ZXH, ZYH 
        real(4) :: S0, S1, S2, SM1, T0, TM1, S0X, S0Y, S1X, S1Y, S2X, S2Y
        real(4) :: SM1X, SM1Y, T0X, T0Y, TM1X, TM1Y
        real(4) :: EDGEX, EDGEY 
        real(4) :: UPT, VPT, DSQ, DPT, S, SX, SY, XPT, YPT, ZPT 
        real(4) :: ZHDEN 

        
    !   Call TRYSBG to calculate subareas and their gradients. 

       CALL TRYSBG(CN,JCNS,JADDR,PT,NPTS,NADDR,NFREE,NSTART,NPTSIN, &
          L,LTOP,LPTR,LBASE,EPSCN,EPSPT,IFLAG,IGARB,XPT,YPT,NINDEX, &
          AREA,SBAREA,DELSBA,PTOFF,KNBMAX,KNB,ICASE) 
!
!   Check IFLAG and return unless it is 0 or 5. 
!
        if (IFLAG .ne. 0 .and. IFLAG .ne. 5) then 
            write(*,*) "NNBRHG:: IFLAG = ", IFLAG 
            return 
        end if 

        if (IFLAG .eq. 5) then 
            !   Deal with the case where we have hit an accepted point. 

            Z0 = VAL(NINDEX) 
            ZX0 = 0.0 
            ZY0 = 0.0 
            ZH = Z0 
            ZXH = 0.0 
            ZYH = 0.0

        else if (IFLAG .eq. 0) then 
            !   Deal with the general case.  Initialise to zero all the 
            !   accumulators for the non-edge case. 

            S0 = 0.0 
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

            !   If needed, initialise accumulators for edge correction. 

            IF(ICASE .ne. 0) then  
                EDGEX = 0.0 
                EDGEY = 0.0 
            end if 

            !   Pick up the contiguity list, and enter a loop through the 
            !   neighbours of the trial point.  Pick up values 
            !   needed to update the accumulators. 

            LLO = LPTR+2 
            LHI = LLO-1+L(LLO-1) 
            K = 0 
            
            DO LL = LLO,LHI 
                K = K+1 
                N = L(LL) 
                IF(N.gt.0) exit
            end do

            UPT = -PTOFF(1,K) 
            VPT = -PTOFF(2,K) 
            DSQ = UPT**2+VPT**2 
            DPT = SQRT(DSQ) 
            S = SBAREA(K) 
            SX = DELSBA(1,K) 
            SY = DELSBA(2,K) 
            ZPT = VAL(N) 

            !   Accumulate values. 

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

            !   Accumulate edge corrector values if needed. 

            IF(ICASE.ne.0) then
                EDGEX = EDGEX+S*UPT 
                EDGEY = EDGEY+S*VPT 
            end if 

            !   Calculate values for return, making edge corrections if needed. 

            Z0 = T0/S0 
            ZX0 = (T0X-Z0*S0X)/S0 
            ZY0 = (T0Y-Z0*S0Y)/S0 
            IF(ICASE.ne.0) then 
                S2X = S2X+2.0*EDGEX 
                S2Y = S2Y+2.0*EDGEY 
            end if
            ZHDEN = S1*S0+S2*SM1 
            ZH = (S1*T0+S2*TM1)/ZHDEN 
            ZXH = S1X*T0+S1*T0X+S2X*TM1+S2*TM1X 
            ZXH = (ZXH-ZH*(S1X*S0+S1*S0X+S2X*SM1+S2*SM1X))/ZHDEN 
            ZYH = S1Y*T0+S1*T0Y+S2Y*TM1+S2*TM1Y 
            ZYH = (ZYH-ZH*(S1Y*S0+S1*S0Y+S2Y*SM1+S2*SM1Y))/ZHDEN 

        end if 

        RETURN 
        
    END subroutine NNBRHG