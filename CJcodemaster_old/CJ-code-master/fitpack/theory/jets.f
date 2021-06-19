      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adzint

C===========================================================================
C GroupName: Adzint
C Description: adaptive integration
C ListOfFiles: adzint adzcal adzspl intusz sglint totalz
C===========================================================================
C                                  Authors: Wu-Ki Tung and John C. Collins
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adzint.f,v 1.1 97/12/21 21:19:04 wkt Exp $
C #Log:	Adzint.f,v $
c Revision 1.1  97/12/21  21:19:04  wkt
c Initial revision
c 

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ (X, FX)
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ------------------------

C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
CError code is IER: 0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C                4 :  Error on Limits : B < A ; zero result returned.
C                5 :  Range of integration DX zero or close to roundoff
C                     returns DX * F(A+DX/2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL, Sml / 1E-20, 1E-12 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        AdzInt = 0D0
        Ier = 4
        If (DDX .Lt. 0D0)
     >     Print*, 'A=', A ,'   B=',B, 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        AdzInt = F(A + DDX/2) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZSPL(F,I,IER)
   40         CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      FUNCTION SMPSN2(FN, A, B, NX, ERR)
c modified version of smpsnf (jcp 11/27/01)
C
C                       Does integral of FN(X)*dx from A TO B by SIMPSON'S METHOD
C
C                       Input:          External function:      FN
C                                       Lower limit      :      A
C                                       Upper limit      :      B
C                                       Number of points :      Nx
C
C                       Uses (Nx-1) evenly spaced intervals.
C
C                       Output:         error estimate:         ERR
C                                       error code    :         IER
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / IOUNIT / NIN, NOUT, NWRT

      PARAMETER (MXPT = 10000)

      DIMENSION X(MXPT)
      external fn

      IF ((NX .LT. 0) .OR. (NX .GT. MXPT)) GOTO 99
C
      DX = (B - A) / FLOAT(NX-1)

c print warning if A > B; but routine SHOULD be ok either way.
c original routine smpsnf set integral to zero if this happened.
      IF (DX .LE. 0) THEN
        WRITE (NOUT, *) 'DX .LE. 0 in SMPSN2, DX =', DX
      ENDIF
C
      DO 10 I = 1, NX
      X(I) = (A*(NX-I) + B*(I-1)) / (NX-1)
   10 CONTINUE
C
      IF (NX .GT. 4) GOTO 50
C
c fast processing for very small NX -- give ultra-conservative error estimates --
      GOTO (20, 30, 40), NX-1
   20 SMPSN2 = (FN(X(1)) + FN(X(2))) * DX / 2.D0
      ERR = ABS(SMPSN2)
      RETURN
   30 SMPSN2 = (FN(X(1)) + 4.D0 * FN(X(2)) + FN(X(3))) * DX / 3.D0
      ERR = ABS(SMPSN2)
      RETURN
   40 SMPSN2 = (( FN(X(1)) + 4.D0 * FN(X(2)) +     FN(X(3))) / 3.D0
     > + (-FN(X(2)) + 8.D0 * FN(X(3)) + 5.D0 * FN(X(4))) / 12.D0 ) * DX
      ERR = ABS(SMPSN2)
      RETURN
C
   50 SE = FN(X(2))
      SO = 0
      NM1 = NX - 1
      DO 60 I = 4, NM1, 2
      IM1 = I - 1
      SE = SE + FN(X(I))
      SO = SO + FN(X(IM1))
   60 CONTINUE
      MS = MOD (NX, 2)
      IF (MS .EQ. 1) THEN
        SMPSN2 = (FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NX))) * DX/3.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO) + FN(X(NX))) * DX/2.D0
      ELSE
        SMPSN2 =(FN(X(1)) + 4.D0*SE + 2.D0*SO + FN(X(NM1))) * DX/3.D0
     > +(-FN(X(NM1-1)) + 8.D0*FN(X(NM1)) + 5.D0*FN(X(NX))) * DX/12.D0
        TRPZ = (FN(X(1)) + 2.D0*(SE + SO + FN(X(NM1))) + FN(X(NX)))
     >          * DX/2.D0
      ENDIF

      ERR = SMPSN2 - TRPZ
c ======================================================================
c print the points...
c	do i = 1, nx
c	   xx = x(i)
c	   ff = fn(xx)
c	   write(nout,666) i, xx, ff
666	   format(1x,'smpsn2:',i5,1x,e12.5,1x,e12.5)
c	enddo
c ======================================================================

      RETURN

   99 WRITE (NOUT, 999) NX
  999 FORMAT (/ 5X, 'NX = ', I6,
     >  'out of range in SIMPSON INTEGRATION SMPSN2')
      STOP
C                        ****************************
      END

      SUBROUTINE ADZCAL (F,I)
C                                                   -=-=- adzcal
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZWRK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)
        CALL SGLINT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGLINT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE TOTALZ
C                                                   -=-=- totalz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adz2nt
      SUBROUTINE ADZSPL (F, I, IER)
C                                                   -=-=- adzspl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZCAL (F, I)
      CALL ADZCAL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END
 
C
C                                                          =-=-= Charutl
      FUNCTION ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adz2nt
 
C===========================================================================
C GroupName: Adz2nt
C Description: second copy of adzint
C ListOfFiles: adz2nt adz2pl adz2al int2sz sgl2nt tot2lz
C=========================================================================== 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adz2nt.f,v 1.1 97/12/21 21:19:00 wkt Exp $
C #Log:	Adz2nt.f,v $
c Revision 1.1  97/12/21  21:19:00  wkt
c Initial revision
c 

C List of GLOBAL Symbols

C     FUNCTION   ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ2PL (F, I, IER)
C     SUBROUTINE ADZ2AL (F,I)
C     SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT2LZ
C     FUNCTION   INT2SZ (X, FX)
C
C     COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                   ------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA SMLL / 1E-20 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZ2NT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      NUMINT = 3
      DX = (B-A)/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZ2AL(F,I)
   10     CONTINUE
       CALL TOT2LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZ2PL(F,I,IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZ2NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

      SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sglint

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END
     
      SUBROUTINE ADZ2AL (F,I)
C                                                   -=-=- adz2al
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D1 = 1.0, D2 = 2.0, HUGE = 1.E15)
C                        Fill in details of interval I given endpoints
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB
 
      SAVE / ADZ2RK /

      DX =  V(I) - U(I)
      W  = (U(I) + V(I)) / 2.
     
      IF (I .EQ. 1 .AND. ICTA .GT. 0) THEN
C                                                                 Open LEFT end
        FW(I) = FA
        FA = F (U(I) + DX / 4.)

        CALL SGL2NT (ICTA, FA, FW(I), FV(I), DX, TEM, ER)
      ELSEIF (I .EQ. IB .AND. ICTB .GT. 0) THEN
C                                                                open RIGHT end
        FW(I) = FB
        FB = F (V(I) - DX / 4.)
        CALL SGL2NT (ICTB, FB, FW(I), FU(I), DX, TEM, ER)
      ELSE
C                                                                   Closed endS
        FW(I) = F(W)
        TEM = DX * (FU(I) + 4. * FW(I) + FV(I)) / 6.
C                                       Preliminary error Simpson - trapezoidal:
        ER  = DX * (FU(I) - 2. * FW(I) + FV(I)) / 12.
      ENDIF
 
      RESULT(I) = TEM         
      ERR   (I) = ABS (ER)
 
      RETURN
C                        ****************************
      END

      SUBROUTINE TOT2LZ
C                                                   -=-=- tot2lz
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      RES = 0.
      ERS = 0.
      DO 10  I = 1, NUMINT
          RES = RES + RESULT(I)
          ERS = ERS + ERR(I)
   10     CONTINUE
C                        ****************************
      END

C                                                          =-=-= Adz3nt
      SUBROUTINE ADZ2PL (F, I, IER)
C                                                   -=-=- adz2pl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                      Split interval I
C                                                   And update RESULT & ERR
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA TINY / 1.D-20 /
     
      IF (NUMINT .GE. MAXINT)  THEN
          IER = 1
          RETURN
          ENDIF
      NUMINT = NUMINT + 1
C                                                         New interval NUMINT
      IF (I .EQ. IB) IB = NUMINT
      U(NUMINT) = (U(I) + V(I)) / 2.
      V(NUMINT) = V(I)
 
      FU(NUMINT) = FW(I)
      FV(NUMINT) = FV(I)
C                                                             New interval I
       V(I) =  U(NUMINT)
      FV(I) = FU(NUMINT)
C                                                    Save old Result and Error
      OLDRES = RESULT(I)
      OLDERR = ERR(I)
     
      CALL ADZ2AL (F, I)
      CALL ADZ2AL (F, NUMINT)
C                                                               Update result
      DELRES = RESULT(I) + RESULT(NUMINT) - OLDRES
      RES = RES + DELRES
C                                  Good error estimate based on Simpson formula
      GODERR = ABS(DELRES) 
C                                                             Update new global 
      ERS = ERS + GODERR - OLDERR
C                                  Improve local error estimates proportionally
      SUMERR = ERR(I) + ERR(NUMINT)
      IF (SUMERR .GT. TINY) THEN
         FAC = GODERR / SUMERR 
      ELSE
         FAC = 1.
      ENDIF
      
      ERR(I)      = ERR(I) * FAC
      ERR(NUMINT) = ERR(NUMINT) * FAC
 
      RETURN
C                        ****************************
      END
 







      SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C                                                   -=-=- sgl2nt

C     Calculate end-interval using open-end algorithm based on function values
C     at three points at (1/4, 1/2, 1)DX from the indeterminant endpoint (0).

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)

      DATA HUGE / 1.E20 /
C                                                         Use quadratic formula
      TEM = DX * (4.*F1 + 3.*F2 + 2.*F3) / 9.
C                 Error est based on Diff between quadratic and linear integrals
      ER  = DX * (4.*F1 - 6.*F2 + 2.*F3) / 9.

C                          Invoke adaptive singular parametrization if IACT = 2
C                      Algorithm is based on the formula F(x) = AA + BB * x **CC
C                 where AA, BB & CC are determined from F(Dx/4), F(Dx/2) & F(Dx)

      IF (IACT .EQ. 2) THEN
          T1 = F2 - F1
          T2 = F3 - F2
          IF (T1*T2 .LE. 0.) GOTO 7
          T3  = T2 - T1
          IF (ABS(T3)*HUGE .LT. T1**2) GOTO 7
          CC  = LOG (T2/T1) / LOG(D2)
          IF (CC .LE. -0.8D0)  GOTO 7
          BB  = T1**2 / T3
          AA  = (F1*F3 - F2**2) / T3
C                                          Estimated integral based on A+Bx**C
          TMP = DX * (AA + BB* 4.**CC / (CC + 1.))
C                                       Error estimate based on the difference
          ER = TEM - TMP
C                                              Use the improved integral value
          TEM= TMP 
      ENDIF

    7 FINT = TEM
      ESTER= ER
      RETURN
C                        ****************************
      END
     
      FUNCTION ALPI (AMU)
C                                                   -=-=- alpi
C  These comments are enclosed in the lead subprogram to survive forsplit.

C ====================================================================
C GroupName: Alphas
C Description: Various callable functions for alpha_s and alpha_em
C ListOfFiles: alpi alpior alepi alpqcd g alphem
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Alphas.f,v 1.3 98/08/16 10:33:20 wkt Exp $ 
C #Log:	Alphas.f,v $
c Revision 1.3  98/08/16  10:33:20  wkt
c Warning in AlpQcd suppressed (cf. comments); numbers corrected.
c 
c Revision 1.2  98/08/11  21:39:31  wkt
c cross-line string argument corrected.
c 
c Revision 1.1  97/12/21  20:34:17  wkt
c Initial revision
c 

C               Returns effective g**2/(4pi**2) = alpha/pi.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
C                      Use the following as subroutine argument of type
C                      set by IMPLICIT statement:
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
 
      DATA IW1, IW2 / 2*0 /
C
      IF(.NOT.SET) CALL LAMCWZ
 
      NEFF = NFL(AMU)
      ALM  = ALAM(NEFF)
      ALPI = ALPQCD (NORDER, NEFF, AMU/ALM, IRT)
 
      IF (IRT .EQ. 1) THEN
         CALL QWARN (IW1, NWRT, 'AMU < ALAM in ALPI', 'MU', AMU,
     >              ALM, BIG, 1)
         WRITE (NWRT, '(A,I4,F15.3)') 'NEFF, LAMDA = ', NEFF, ALM
      ELSEIF (IRT .EQ. 2) THEN
         CALL QWARN (IW2, NWRT, 'ALPI > 3; Be aware!', 'ALPI', ALPI,
     >             D0, D1, 0)
         WRITE (NWRT, '(A,I4,2F15.3)') 'NF, LAM, MU= ', NEFF, ALM, AMU
      ENDIF
 
      RETURN
      END
C
C************
C
      FUNCTION NFL(AMU)
C                                                   -=-=- nfl
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C           NFL returns the number of 'light' flavors.
C
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      LOGICAL SET
C
      IF (.NOT. SET) CALL LAMCWZ
      NFL = NF - NHQ
      IF ((NFL .EQ. NF) .OR. (AMU .LE. AMN)) GOTO 20
      DO 10 I = NF - NHQ + 1, NF
         IF (AMU .GE. AMHAT(I)) THEN
            NFL = I
         ELSE
            GOTO 20
         ENDIF
10       CONTINUE
20    RETURN
      END
C
C***************************************************************
C
      SUBROUTINE LAMCWZ
C                                                   -=-=- lamcwz
C                       Set /CWZPRM/ from /QCDPAR/
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / QCDPAR / AL, NF, NORDER, SET
      LOGICAL SET
      CALL SETL1 (NF, AL)
      END
C
C***********************
C
      SUBROUTINE QWARN (IWRN, NWRT1, MSG, NMVAR, VARIAB,
C                                                   -=-=- qwarn
     >                  VMIN, VMAX, IACT)
 
C     Subroutine to handle warning messages.  Writes the (warning) message
C     and prints out the name and value of an offending variable to SYS$OUT
C     the first time, and to output file unit # NWRT1 in subsequent times.
C
C     The switch IACT decides whether the limits (VMIN, VMAX) are active or
C     not.
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
 
      CHARACTER*(*) MSG, NMVAR

      Save Iw
 
      IW = IWRN
      VR = VARIAB
 
      If (Iw .LT. 100) Then
         WRITE (NWRT1,'(I5, 3X,A/ 1X,A,'' = '',1PD16.7)') IW, MSG,
     >                  NMVAR, VR
      Else
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' !! Error # > 100 !! ; better check file unit #', NWRT1
      EndIf   

      IF  (IW .EQ. 0) THEN
         WRITE (NOUT, '(1X, A/1X, A,'' = '', 1PD16.7/A,I4)')
     >      MSG, NMVAR, VR,
     >      ' Complete set of warning messages on file unit #', NWRT1
         IF (IACT .EQ. 1) THEN
         WRITE (NOUT,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         WRITE (NWRT1,'(1X,A/2(1PD15.3))')'The limits are: ', VMIN,VMAX
         ENDIF
      ENDIF
 
      IWRN = IW + 1
 
      RETURN
C                         *************************
      END

      FUNCTION ALPQCD (IRDR, NF, RML, IRT)
C                                                   -=-=- alpqcd
 
C                                 Returns the QCD alpha/pi for RML = MU / LAMDA
C                                 using the standard perturbative formula for
C                                 NF flavors and to IRDR th order in 1/LOG(RML)
 
C                                 Return Code:  IRT
C                                                0:   O.K.
C                                                1:   Mu < Lamda; returns 99.
C                                                2:   Alpha > 10 ; be careful!
C                                                3:   IRDR out of range
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, BIG = 1.0D15)
      PARAMETER (CG = 3.0, TR = 0.5, CF = 4.0/3.0)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      DATA IW1, Iw2 / 0, 0 /
 
      IRT = 0
 
      IF (IRDR .LT. 1 .OR. IRDR .GT. 2) THEN
        WRITE(NOUT, *)
     >  'Order parameter out of range in ALPQCD; IRDR = ', IRDR
        IRT = 3
        STOP
      ENDIF
 
      B0 = (11.* CG  - 2.* NF) / 3.
      B1 = (34.* CG**2 - 10.* CG * NF - 6.* CF * NF) / 3.
      RM2 = RML ** 2

C           AlpQcd is used mainly as a mathematical function, for
C            inversion, as well as evaluation of the physical Alpi.
C           Warning should be deferred to the calling program.
 
      IF (RM2 .LE. 1.) THEN
         IRT = 1
C         CALL QWARN(IW1, NWRT,
C     >    'RM2 (=MU/LAMDA) < 1. not allowed in ALPQCD; Alpha->99.',
C     >    'RM2', RM2, D1, BIG, 1)
         ALPQCD = 99
         RETURN
      ENDIF
 
      ALN = LOG (RM2)
      AL = 4./ B0 / ALN
 
      IF (IRDR .GE. 2) AL = AL * (1.- B1 * LOG(ALN) / ALN / B0**2)
 
      IF (AL .GE. 3.) THEN
         IRT = 2
C         CALL QWARN(IW2, NWRT, 'ALPQCD > 3. in ALPQCD', 'ALPQCD', AL,
C     >              D0, D1, 1)
      ENDIF
 
      ALPQCD = AL
 
      RETURN
C                       *********************
      END
C
      SUBROUTINE SETL1  (NEF, VLAM)
C                                                   -=-=- setl1
C     Given LAMDA = VLAM for NEF flavors:
C                    (i) fills the array  ALAM (0:NF) with effective LAMDA's;
C                    (ii) fills the array AMHAT (NF) with threshold masses;
C                    (iii) count the # of "heavy quarks" (QMS > EFFLAM);
C                    (iv) fix the parameter AMN defined as MAX (ALAM),
C                         times safety factor;
C                    (v) set AL in / QCDPAR / equal to ALAM (NF);
C                    (vi) let SET = .TRUE.
C       Uses formula with NORDER (1 or 2) -- see /QCDPAR/
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      LOGICAL SET
 
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / QCDPAR / AL, NF, NORDER, SET
      COMMON / COMQMS / QMS(9)
      COMMON / IOUNIT / NIN, NOUT, NWRT
 
      IF (NEF .LT. 0 .OR. NEF .GT. NF) THEN
        WRITE(NOUT,*)'NEF out of range in SETL1, NEF, NF =',NEF,NF
        STOP
      ENDIF
C             Mass Thresholds are given by the Quark masses in the CWZ scheme
      AMHAT(0) = 0.
      DO 5 N = 1, NF
         AMHAT(N) = QMS(N)
    5    CONTINUE
      ALAM(NEF) = VLAM
      DO 10 N = NEF, 1, -1
         CALL TRNLAM(NORDER, N, -1, IR1)
   10    CONTINUE
      DO 20 N = NEF, NF-1
         CALL TRNLAM(NORDER, N, 1, IR1)
   20    CONTINUE
C=========================                Find first light quark:
      DO 30, N = NF, 1, -1
         IF ((ALAM(N) .GE. 0.7 * AMHAT(N))
     >       .OR. (ALAM(N-1) .GE. 0.7 * AMHAT(N)))THEN
            NHQ = NF - N
            GOTO 40
            ENDIF
   30    CONTINUE
      NHQ = NF
   40 CONTINUE
      DO 50, N = NF-NHQ, 1, -1
         AMHAT(N) = 0
         ALAM(N-1) = ALAM(N)
   50    CONTINUE
C========================               Find minimum mu
      AMN = ALAM(NF)
      DO 60, N = 0, NF-1
         IF (ALAM(N) .GT. AMN)  AMN = ALAM(N)
   60    CONTINUE
      AMN = AMN * 1.0001
      AL = ALAM(NF)
      SET = .TRUE.
      RETURN
C**************************************************************
      END

C                                                          =-=-= Setaux
      SUBROUTINE TRNLAM (IRDR, NF, IACT, IRT)
C                                                   -=-=- trnlam
 
C     This routine transforms LAMDA (N) to LAMDA (N+IACT) where IACT = 1/-1
C     The transformation is obtained by requiring the coupling constant to
C                be continuous at the scale Mu = Mass of the (N+1)th quark.
 
C                                         IRT is an return code.
C                                            (0 for OK)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 
      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / CWZPRM / ALAM(0:9), AMHAT(0:9), AMN, NHQ
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      EXTERNAL ZBRLAM
 
      DATA ALM0, BLM0, RERR / 0.01, 10.0, 0.0001 /
      DATA IR1, SML / 0, 1.E-5 /
 
      IRT = 0
 
      N = NF
      JRDR = IRDR
      JACT = IACT
      VLAM = ALAM(N)
 
      IF (JACT .GT. 0) THEN
         N1 = N + 1
         THMS = AMHAT(N1)
         ALM = LOG (THMS/VLAM)
         BLM = BLM0
      ELSE
         N1 = N -1
         THMS = AMHAT(N)
         ALM = ALM0
         THMS = MAX (THMS, SML)
         BLM = LOG (THMS/VLAM)
      ENDIF
C                          Fix up for light quark:
      IF (VLAM .GE. 0.7 * THMS) THEN
         IF (JACT . EQ. 1) THEN
            AMHAT(N1) = 0
         ELSE
            AMHAT(N) = 0
         ENDIF
         IRT = 4
         ALAM(N1) = VLAM
         RETURN
      ENDIF
 
C             QZBRNT is the root-finding function to solve ALPHA(N) = ALPHA(N1)
C             Since 1/Alpha is roughly linear in Log(Mu/Lamda), we use the
C             former in ZBRLAM and the latter as the function variable.
      IF (ALM .GE. BLM) THEN
         WRITE (NOUT, *) 'TRNLAM has ALM >= BLM: ', ALM, BLM
         WRITE (NOUT, *) 'I do not know how to continue'
         STOP
         ENDIF
      VMULM = THMS/VLAM
      ERR = RERR * LOG (VMULM)
      WLLN = QZBRNT (ZBRLAM, ALM, BLM, ERR, IR1)
      ALAM(N1) = THMS / EXP (WLLN)
 
      IF (IR1 .NE. 0) THEN
         WRITE (NOUT, *) 'QZBRNT failed to find VLAM in TRNLAM; ',
     >        'NF, VLAM =', NF, VLAM
         WRITE (NOUT, *) 'I do not know how to continue'
        STOP
      ENDIF
      RETURN
      END
C                             *************************
 
      FUNCTION ZBRLAM (WLLN)
C                                                   -=-=- zbrlam
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON / TRNCOM / VMULM, JRDR, N, N1
 
      WMULM = EXP (WLLN)
      TEM1 = 1./ ALPQCD(JRDR, N1, WMULM, I)
      TEM2 = 1./ ALPQCD(JRDR, N,  VMULM, I)
 
      ZBRLAM = TEM1 - TEM2
 
      END
 
 
C************************************************************
 
      FUNCTION QZBRNT(FUNC, X1, X2, TOLIN, IRT)
C                                                   -=-=- qzbrnt
 
C                          Return code  IRT = 1 : limits do not bracket a root;
C                                             2 : function call exceeds maximum
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /IOUNIT/ NIN, NOUT, NWRT
      PARAMETER (ITMAX = 1000, EPS = 3.E-12)

      external func
 
      TOL = ABS(TOLIN)
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF(FB*FA.GT.0.)  THEN
        WRITE (NOUT, *) 'Root must be bracketed for QZBRNT.'
        IRT = 1
      ENDIF
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          QZBRNT=B
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
      WRITE (NOUT, *) 'QZBRNT exceeding maximum iterations.'
      IRT = 2
      QZBRNT=B
      RETURN
C**************************************************
      END
C
      Subroutine Read1jet(kfile)

      Character kfile*(*)
      Parameter (MxNPJ=100,MxNYJ=9,MxJset=2)
      Integer NPJ(MxJset),NYJ(MxJset),Jet
      Double Precision Y1grid(MxNYJ,MxJset), PJ(MxNPJ,MxNYJ,MxJset), 
     >                 OnejKfac(MxNPJ,MxNYJ,MxJset)
      Common /OnejetArrays/ Y1grid,PJ,OnejKfac
      Common /OnejetNPY/ NPJ, NYJ
      Common /Onejetset/ Jset
C IO unit numbers
C
      INTEGER NIN,NOUT,NWRT
      COMMON /IOUNIT/ NIN,NOUT,NWRT
      Character Lin*80
      Integer IU, NextUN
      Integer I,J
      Double Precision Y1min, Y1max

      WRITE(NOUT,'(A,A)') ' reading ',kfile
	IU=NextUN()
      OPEN(IU,FILE=kfile,STATUS='OLD',ERR=99)

      Read(IU,'(A)') Lin
      Read(IU,*) NPJ(Jset), NYJ(Jset)
      If(NPJ(Jset).gt.MxNPJ .or. NYJ(Jset).gt.MxNYJ) goto 109
      Read(IU,'(A)') Lin
      Do J=1,NYJ(Jset)
         Do I=1,NPJ(Jset)
            Read(IU,*) PJ(I,J,Jset), Y1min, Y1max, OnejKfac(I,J,Jset)
         EndDo
         Y1grid(J,Jset) = Y1min
      EndDo

      close(IU)
      RETURN
99    WRITE(NOUT,'(A,A)')'Error opening data file for input:',Kfile
      STOP
 109  WRITE(NOUT,'(3A,2I9)')'NPJ or NYJ in ',Kfile,' is too large:'
     >   ,NPJ(Jset),NYJ(Jset)
      Stop

      End

      Double Precision Function OnejetKf(Et,y1)

      Double Precision Et, y1
      Parameter (MxNPJ=100,MxNYJ=9,MxJset=2)
      Integer NPJ(MxJset),NYJ(MxJset),Jset
      Double Precision Y1grid(MxNYJ,MxJset), PJ(MxNPJ,MxNYJ,MxJset),
     >                 OnejKfac(MxNPJ,MxNYJ,MxJset)
      Common /OnejetArrays/ Y1grid,PJ,OnejKfac
      Common /OnejetNPY/ NPJ, NYJ
      Common /Onejetset/ Jset

      Integer I,J,II,JJ
      Double Precision akf(3),Df(3),bkf

      Do 20 J=1,NYJ(Jset)
         If(y1 .eq. Y1grid(J,Jset)) goto 21
 20   Continue
 21   JJ=J

      Do 10 I=2,NPJ(Jset)-1
         If(Et .lt. PJ(I,JJ,Jset)) goto 11
 10   Continue
 11   II=I-1

      DO 30 I=1,3
         akf(I)= OnejKfac(II+I-1,JJ,Jset)
 30   Continue
      Call Polint(PJ(II,JJ,Jset),akf,3, Et, bkf, Df)
      OnejetKf = bkf
      Return
      End

      Subroutine JetdEt( Rs, Et, Smucoef, Result, y1mn, y1mx, Inlo)

      Implicit None

      Double Precision Rs, Et, Smucoef, Result, y1mn, y1mx
      Integer Inlo

      Double Precision RootS, Pt, y1, y2min, y2max
      Common /Kinematics/ RootS, Pt, y1, y2min, y2max
      DOUBLE PRECISION Mu
      COMMON /Scale/ Mu
c      Integer NLO
c      Common /Jetorder/ NLO

      Double Precision E2max, y1min, y1max, y1range, avgi
      Double Precision AERR, RERR, ERREST
      Integer IER, IACTA, IACTB

      Double Precision adzint
      Double Precision Dsigmay1
      External Dsigmay1
      Double Precision OnejetKf

      Data AERR, RERR , IACTA, IACTB/ 0D0, .001D0 ,0 ,0/
c  assumed in CDF Jets case (only calculate rapidity .1 ~ .7 and average over)
      Data y1min,y1max,y1range / .1D0, .7D0, .6D0 /

      RootS= Rs
      Pt = Et
      Mu = Smucoef * Pt
c  NLO set to be 0 for onejet. NLO K-fac is calculated in OnejetKf below.
c      NLO = Inlo
      y1min= y1mn
      y1max=y1mx
      y1range= y1max - y1min
      E2max= RootS - Pt
      y2max= log (( E2max + (E2max**2 - Pt**2)**0.5D0)/Pt)
      y2min= log (( E2max - (E2max**2 - Pt**2)**0.5D0)/Pt)
      avgi= adzint
     >  (Dsigmay1, y1min, y1max, AERR, RERR, ERREST, IER, IACTA, IACTB)

      Result= avgi / y1range
c      print*,result
      If (Inlo.gt.0) Result = Result * OnejetKf(Et,y1mn)
c      print*,result
      Return
      End


c ********************************************************************
      Double Precision Function Dsigmay1(y)

      Implicit None

      Double Precision y
      Double Precision AERR, RERR, ERREST
      Integer IER, IACTA, IACTB

      Double Precision adz2nt, Dsigmay2

      Double Precision RootS, Pt, y1, y2min, y2max
      Common /Kinematics/ RootS, Pt, y1, y2min, y2max

c ****************************************
      Double Precision DelYmax
      Common /DeltaYmax/ DelYmax
	double precision xmax,xt2,y2mx,y2mn
c ****************************************

      External Dsigmay2

      Data AERR, RERR, IACTA, IACTB / 0D0, .001D0, 0, 0 /

      y1=y

c include the kinematic limits directly...
	y2mn = y2min
	y2mx = y2max

	xMax = 0.98d0 		!taken from dsigmay2

	Xt2 = Pt / RootS
	y2mx = min(y2mx,log(xmax/xt2 - exp(y1)))
	y2mn = max(y2mn,-log(xmax/xt2 - exp(-y1)))

	y2mx = min(y2mx,y1 + delymax)
	y2mn = max(y2mn,y1 - delymax)

	if(y2mn .ge. y2mx) then
	   dsigmay1 = 0.d0
	   return
	endif
      Dsigmay1= adz2nt
     >  (Dsigmay2, y2mn, y2mx, AERR, RERR, ERREST, IER, IACTA, IACTB)

c      print*,dsigmay1,y1,y2mn,y2mx
      Return
      End

c ********************************************************************
      Double Precision Function Dsigmay2( Y )

      Implicit None
c
c   Jet cross section
c   d(Sigma)/[d(Pt) d(y1) d(y2)] * Jacobian(x1, x2, y1, y2)
c   
c

      Integer NumType, QQ, QQbar, QiQj, GQ, GG
      Integer NumPro, QiQj_QiQj, QQb_QiQj, QQ_QQ
      Integer QQb_QQb, GQ_GQ, QQb_GG, GG_QQb, GG_GG
      Double Precision Pi, GEV_NB
      Parameter (NumType=5, QQ=1, QQbar=2, QiQj=3, GQ=4, GG=5)
      Parameter (NumPro=8, QiQj_QiQj=1, QQb_QiQj=2, QQ_QQ=3 )
      Parameter (QQb_QQb=4, GQ_GQ=5, QQb_GG=6, GG_QQb=7, GG_GG=8)
      Parameter (Pi=3.14159265D0, GEV_NB=3.8937966D5)

      Double Precision Y

      Double Precision RootS, Pt, y1, y2min, y2max
      Common /Kinematics/ RootS, Pt, y1, y2min, y2max
      Double Precision Mu
      Common /Scale/ Mu
      Double Precision SumP
      Common /SumPn/ SumP(NumType)
      Double Precision DelYmax
      Common /DeltaYmax/ DelYmax

      Double Precision y2, Xt2, DeltaY, S, x1, x2
      Double Precision  Jac
      Double Precision sHat, tHat, uHat, sHat2, tHat2, uHat2,
     >    stHat, tuHat, usHat, stHat2, tuHat2, usHat2
      Double Precision Alpha, KinFac, Dsigma0
      Double Precision Process(NumPro)
      Double Precision oldMu, xMax, Kfac, Ystar

c      Integer final

c      Double Precision Alpqcd, RML, Lambda, Alamf
c      Integer Iorder, Nf, Ir
      Double Precision Alpi, JetKfac

      Integer I
c      Integer NLO
c      Common /Jetorder/ NLO

      Data xMax / .98D0 /
      Data oldMu /0D0/
      save oldMu, kinfac

c ********************************************
	integer nin, nout, nwrt
      COMMON /IOUNIT/ NIN,NOUT,NWRT
c ********************************************
      y2 = y

      Xt2 = Pt / RootS
      DeltaY = y1-y2

      S = RootS * RootS
      x1 = Xt2 * ( Exp(y1) + Exp(y2) )
      x2 = Xt2 * ( Exp(-y1) + Exp(-y2) )


      If (x1.ge.xMax .or. x2.ge.xMax .or. Abs(DeltaY).ge.DelYmax) then
         Dsigmay2 = 0.D0
         Return
      Endif

      Ystar = DeltaY / 2D0
c      Call SetMu(x1,x2,Ystar,MUUV,MUCO)

      sHat = x1 * x2 * S
      tHat = - sHat / ( 1.D0 + Exp(DeltaY) )
      uHat = -sHat -tHat

      stHat = sHat * tHat
      tuHat = tHat * uHat
      usHat = uHat * sHat
      sHat2 = sHat * sHat
      tHat2 = tHat * tHat
      uHat2 = uHat * uHat
      stHat2 = sHat2 + tHat2
      tuHat2 = tHat2 + uHat2
      usHat2 = uHat2 + sHat2

c      JacXY = Xt2 * Xt2 * Abs( Exp(DeltaY) - Exp(-DeltaY) )
c      JacPt_tHat =  2.D0*Pt * Abs( sHat/(tHat-uHat) )
c      To avoid DeltaY=0 which causes 0 * Infinity
c      combines JacXY and JacPt_tHat and gives a positive definite value.
c      Jac= JacXY * JacPt_tHat
c         = Xt2 * Xt2 * 2 * Pt * ( Exp(DeltaY) + Exp(-DeltaY) + 2 )

      Jac  = Xt2 * Xt2 * 2D0 * Pt
     >        * ( Exp(DeltaY) + Exp(-DeltaY) + 2D0 )

      Call SumPartons( x1, x2 )

c    consider t <-> u channel for final states

      Process(QiQj_QiQj)= Sump(QiQj) * 4.D0/9.D0 *
     >	  ( usHat2/tHat2 + stHat2/uHat2)
      Process(QQb_QiQj)= Sump(QQbar) * 4.D0/9.D0 * tuHat2/sHat2 * 8D0
c                                                         *2*4 for final states
      Process(QQ_QQ)= Sump(QQ) * ( 4.D0/9.D0 *
     >    ( usHat2/tHat2 + stHat2/uHat2 ) - 8.D0/27.D0* sHat2/tuHat )

      Process(QQb_QQb)= Sump(QQbar) * ( 4.D0/9.D0 *
     >    ( usHat2/tHat2 + 2D0 * tuHat2/sHat2 + stHat2/uHat2 )
     >    - 8.D0/27.D0 * ( uHat2/stHat + tHat2/usHat) )

      Process(GQ_GQ)= Sump(GQ) * 
     >    ( (usHat2/tHat2 + stHat2/uHat2 )
     >     - 4.D0/9.D0 * (usHat2/usHat + stHat2/stHat))

      Process(QQb_GG)= Sump(QQbar) * ( 32.D0/27.D0 * 
     >    tuHat2/tuHat - 8.D0/3.D0 * tuHat2 / sHat2 )

      Process(GG_QQb) =  Sump(GG) * ( 1.D0/6.D0 * 
     >    tuHat2/tuHat - 3.D0/8.D0 * tuHat2 / sHat2) * 10D0
c                                                          2*5 for final states
      Process(GG_GG)= Sump(GG) * 9.D0/2.D0 *
     >     ( 3.D0 - tuHat/sHat2 - usHat/tHat2 - stHat/uHat2 )
      Dsigma0= 0D0
      Do 10 I=1, NumPro
         Dsigma0= Dsigma0 + Process(I)
  10  Continue
      Dsigma0= Dsigma0/sHat2

      If (oldMu .ne. Mu) then
        Alpha= Alpi(Mu) * Pi
	KinFac= Pi*Alpha*Alpha * GEV_NB
        oldMu=Mu
c        print*,alpha, kinfac, mu
      Endif

      Dsigmay2= Dsigma0 * KinFac * Jac
c      print*,dsigmay2, y2, y1, dsigma0, kinfac, jac
c      If (NLO.gt.0) then
c         Kfac=JetKfac(x1,x2,Ystar)
c         If(Kfac.lt.0D0) Kfac=0D0
c         Dsigmay2 = Dsigmay2 * Kfac
c ******************************************************************
c	if(Kfac .eq. 0.d0) write(nout,16)x1,x2,ystar 	!*** strange that this can happen!
c16	format(1x,'dsigmay2: Kfac=0, x1,x2,ystar=',3e11.3)
c ******************************************************************
c      Endif

c      If(final().eq.1) then
c         Jac=Jac*JacVegas*KinFac/sHat2 * WGT
c         call Separate(Process,Jac)
c      Endif

c      If (Dsigmay2 .lt. 0D0) then
c 	 Print*, 'x1, x2, y1, y2=',x1, x2, y1, y2
c         Print*, 'Warning: Dsigmay2, Dsigam0='
c     >    , Dsigmay2, Dsigma0
c         Print*,'Rs,Pt,Mu,Kfac:',RootS,Pt,Mu,Kfac
c	 Print*,'sHat,Alpha:', sHat, Alpha

c         Dsigmay2 = 0D0
c      Endif

      Return
      End

c *******************************************************************

      Subroutine SumPartons( x1, x2 )

      Implicit None
      Integer NumType, QQ, QQbar, QiQj, GQ, GG, MaxParton
      Parameter (NumType=5, QQ=1, QQbar=2, QiQj=3, GQ=4, GG=5)
c      Parameter ( MaxParton=5 )

      Double Precision x1, x2

      Double Precision SumP(NumType)
      Common /SumPn/ SumP
      Double Precision Mu
      Common /Scale/ Mu
      Integer Iset, Hadron1, Hadron2, Ir
      Common /ForPDF/ Iset, Hadron1, Hadron2
      common/jetpdfs/ndrv,numfl,als
      integer ndrv,numfl
      Double Precision als,u,d,ub,db,sb,cb,bb,glue
c      Double Precision PDF
      Integer NFL
      Double Precision Parton(-6:6, 2)
      Integer I, J

      MaxParton = NFL(Mu)
c NFL(Mu) returns the number of `light' flavors at scale Mu - effective

c      Do 10 I= -MaxParton, MaxParton
c         Parton(I,1) = PDF(Iset, Hadron1, I, x1, Mu, Ir)
c         Parton(I,2) = PDF(Iset, Hadron2, I, x2, Mu, Ir)
c  10  Continue
c
c  interfaced to my routines
c  set up for p pbar only
c
      call fsupdf(ndrv,x1,als,u,d,ub,db,sb,cb,bb,glue)
      parton(-5,1)=bb/x1
      parton(-4,1)=cb/x1
      parton(-3,1)=sb/x1
      parton(-2,1)=db/x1
      parton(-1,1)=ub/x1
      parton(0,1)=glue/x1
      parton(1,1)=u/x1
      parton(2,1)=d/x1
      parton(3,1)=parton(-3,1)
      parton(4,1)=parton(-4,1)
      parton(5,1)=parton(-5,1)
c
c  pbar
c
      call fsupdf(ndrv,x2,als,u,d,ub,db,sb,cb,bb,glue)
      parton(-5,2)=bb/x2
      parton(-4,2)=cb/x2
      parton(-3,2)=sb/x2
      parton(-2,2)=d/x2
      parton(-1,2)=u/x2
      parton(0,2)=glue/x2
      parton(1,2)=ub/x2
      parton(2,2)=db/x2
      parton(3,2)=parton(-3,2)
      parton(4,2)=parton(-4,2)
      parton(5,2)=parton(-5,2)
c
      Do 11 I=1, NumType
         SumP(I) = 0.D0
  11  Continue

      Do 12 I=1, MaxParton
         SumP(QQ)= SumP(QQ) + Parton(I,1)*Parton(I,2)
     >               + Parton(-I,1)*Parton(-I,2)
         SumP(QQbar)= SumP(QQbar) + Parton(I,1)*Parton(-I,2)
     >               + Parton(-I,1)*Parton(I,2)
         SumP(GQ)= SumP(GQ) + Parton(0,1)*(Parton(I,2)+Parton(-I,2))
     >               + Parton(0,2)*(Parton(I,1)+Parton(-I,1))
         Do 21 J= 1, MaxParton
            If (I.ne.J) then
               SumP(QiQj)= SumP(QiQj)
     >          + (Parton(I,1)+Parton(-I,1))*(Parton(J,2)+Parton(-J,2))
            Endif
  21     Continue
  12  Continue
      SumP(GG) = Parton(0,1) * Parton(0,2)

      End

      Subroutine DiJetdEt
     >  ( Rs, Et, Smucoef, Result, y1mn, y1mx, y2mn, y2mx, Inlo)

      Implicit None

      Double Precision Rs, Et, Smucoef, Result
      Double Precision y1mn, y1mx, y2mn, y2mx
      Integer Inlo

      Double Precision RootS, Pt, y1, y2min, y2max
      Common /Kinematics/ RootS, Pt, y1, y2min, y2max
      DOUBLE PRECISION Mu
      COMMON /Scale/ Mu
      Integer NLO
      Common /Jetorder/ NLO

      Double Precision E2max, y1min, y1max, y1range, avgi
      Double Precision y2range
      Double Precision AERR, RERR, ERREST
      Integer IER, IACTA, IACTB

      Double Precision adzint, DijetKf
      Double Precision Dsigmay1
      External Dsigmay1

      Data AERR, RERR, IACTA, IACTB / 0D0, .001D0, 0, 0 /

      RootS= Rs
      Pt = Et
      Mu = Smucoef * Pt

c  NLO set to be 0 for dijet. NLO K-fac is calculated in DijetKf below.
      NLO = 0

      y1min= y1mn
      y1max=y1mx
      y1range= y1max - y1min
      E2max= RootS - Pt
c      y2max= log (( E2max + (E2max**2 - Pt**2)**0.5D0)/Pt)
c      y2min= log (( E2max - (E2max**2 - Pt**2)**0.5D0)/Pt)
      y2min=y2mn
      y2max=y2mx
      y2range= 2D0 * (y2max - y2min)

      avgi= adzint
     >  (Dsigmay1, y1min, y1max, AERR, RERR, ERREST, IER, IACTA, IACTB)

      y2min=-y2mx
      y2max=-y2mn
      avgi= avgi + adzint
     >  (Dsigmay1, y1min, y1max, AERR, RERR, ERREST, IER, IACTA, IACTB)

      Result= avgi / y1range /y2range
      If (Inlo.gt.0) Result = Result * DijetKf(Et,y2mn)

      Return
      End

      Subroutine ReadDijet(kfile)

      Character kfile*(*)
      Parameter (MxNPJ=100,MxNYJ=4)
      Integer I,J,NPJ,NYJ,IU
      Double Precision Y2min, Y2max
      Double Precision Y2grid(MxNYJ+1), PJ(MxNPJ), dijKfac(MxNPJ,MxNYJ)
      Common /DijetArrays/ Y2grid,PJ,dijKfac
      Common /DijetNPY/ NPJ, NYJ
C IO unit numbers
C
      INTEGER NIN,NOUT,NWRT
      COMMON /IOUNIT/ NIN,NOUT,NWRT
      Character Lin*80

      WRITE(NOUT,'(A,A)') ' reading ',kfile
      IU=NextUn()
      OPEN(IU,FILE=kfile,STATUS='OLD',ERR=99)

      Read(IU,'(A)') Lin
      Read(IU,*) NPJ, NYJ
      If(NPJ .gt. MxNPJ .or. NYJ.gt.MxNYJ) goto 109
      Read(IU,'(A)') Lin
      Do J=1,NYJ
         Do I=1,NPJ
            Read(IU,*) PJ(I), Y2min, Y2max, dijKfac(I,J)
         EndDo
         Y2grid(J) = Y2min
      EndDo
      Y2grid(NYJ+1) = Y2max

      close(IU)
      RETURN
99    WRITE(NOUT,'(A,A)')'Error opening data file for input:',Kfile
      STOP
 109  WRITE(NOUT,'(3A,2I9)')'NPJ or NYJ in ',Kfile,' is too large:'
     >   ,NPJ,NYJ
      Stop

      End

      Double Precision Function DijetKf(Et,y2)

      Double Precision Et, y2
      Parameter (MxNPJ=100,MxNYJ=4)
      Integer NPJ,NYJ
      Double Precision Y2grid(MxNYJ+1), PJ(MxNPJ), dijKfac(MxNPJ,MxNYJ)
      Common /DijetArrays/ Y2grid,PJ,dijKfac
      Common /DijetNPY/ NPJ, NYJ

      Integer I,J,II,JJ
      Double Precision akf(3),Df(3),bkf

      Do 10 I=2,NPJ-1
         If(Et .lt. PJ(I)) goto 11
 10   Continue
 11   II=I-1

      Do 20 J=2,NYJ+1
         If(y2 .lt. Y2grid(J)) goto 21
 20   Continue
 21   JJ=J-1

      DO 30 I=1,3
         akf(I)= dijKfac(II+I-1,JJ)
 30   Continue

      Call Polint(PJ(II),akf,3, Et, bkf, Df)
      DijetKf = bkf

      Return
      End

      FUNCTION ChiJet (Npt, ExDat, ExErr, dvJet)
C                                                   -=-=- chijet

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: ChiJet
C Description: routines to compute ChiSq for JET processes
C ListOfFiles: chijet chid0jet chicdfjet
C ====================================================================

C $Header: /group/cteqX/CVS/fitpack/theory/jets.f,v 1.3 2011/08/02 21:17:22 accardi Exp $ 
C $Log: jets.f,v $
C Revision 1.3  2011/08/02 21:17:22  accardi
C *** empty log message ***
C
C Revision 1.2  2011/03/15 18:24:03  jowens
C modified pdf call --> fsupdf
C
C Revision 1.1.1.4  2010/03/12 20:06:09  accardi
C NextUn function from jet.f put in util/io.f
C
C Revision 1.1.1.3  2010/01/26 22:53:07  jowens
C New files for the updated fitting and evolution packages
C
C Revision 1.1.1.2  2009/04/02 17:32:24  accardi
C New data sets, improvements in 'altpfit',  3D deuterium convolution, off-shell corrections, April_09 fitting subfolder
C
C Revision 1.1.1.1  2008/08/12 16:54:21  accardi
C Jeff Owen's global QCD fit package - May 2008
C
C Revision 1.1.1.1  2008/08/11 14:52:11  accardi
C Jeff Owen's global QCD fitting package - May 2008
C
C Revision 1.3  2002/01/25 15:06:59  wkt
C Add DIS scheme option
C
C Revision 1.2  2001/09/17 01:38:03  wkt
C The SysErr switch Lsw3 is changed to 10-based, rather than Liang's 4-based.
C
c Revision 65.1  99/05/19  10:54:20  wkt
c Synchronize to v65
c Cf. Log message in FcnEtc for v65 changes;
c Also, more changes made on Lprt use; Chizz4 now has one less argument.
c 
c Revision 6.9  99/02/13  16:28:22  wkt
c CDFjet module corrected for mistake (by CDF).  Use Lsw3=51 as switch to
c turn on convariance matrix.
c 
c Revision 6.8  98/12/24  11:50:37  wkt
c 1. Make use of Lsw2=51 to activate the use of ChiD0jet and ChiCdfJet
c routines to calculate correlated systematic errors.
c 2. The CDF module is addted for the first time.
c 
c Revision 6.7  98/10/03  10:46:49  wkt
c Normalization factor added by Liang; double precision function --> function
c for ease of handling by packgrp.r
c 
c Revision 6.0  98/03/08  14:03:17  wkt
c New start from Liang's "CTEQ4" package dated September 1997.
c 
C     --------------------------------------------
C                        CALCULATES ChiSq FOR Single Inclusive Jet Experiments
C
      Implicit Double Precision(A-H,O-Z)
 
      Parameter(MxDt = 5000, Mblk = 50, MxSys = 20)
 
      Character LiOut(Npt)*120, Kfile*40
      Character NmDt*10, FlPath*30, FlExplis*15
 
      COMMON
     >  /IoUnit/ Nin, Nout, Nwrt
     >  /LUnits/ IRmin, IWmin, IPunch, IUDta, IUSum, IUtmp, IUdtf
     >  /FitSwh/ Lfit, LErr, Lprt, Loutput
     >  /FitHrd/ AlfQ,Qalf,aMcbt(3),ISch,IorHrd,IorQcd,Lqcd(4),Iset
     >  /FcnHrd/ IB, ID, Ifg, N10, Irun
     >  /ThpJet/ ThJet(3), LJet(3)
     >  /DatSet/ Jsfn(Mblk),Jbsn(Mblk), Jtgt(Mblk), Jbem(Mblk)
     >  /ExpPar/ ExDis(4),ExVbp(4),ExDph(4),ExHhk(4),ExJet(4),ExXxx(4)
     >  /FlNams/ FlPath, FlExplis, NmDt(Mblk)
     >  /CorMtx/sys(MxSys,MxDt),Uncor(MxDt),Nptcor(Mblk),Ncor(Mblk)

      Common /ForPDF/ IPDF, NTAGT, NBEAM
      Common /DeltaYmax/ DelYmax
      Common /Onejetset/ Jset

      Double Precision ExDat(9,Npt),ExErr(Npt), FitT(Npt)

      Logical FirstFermi, FirstUA2, FirstCDFdijet, FirstF2
c   UA2 Jets case (only calculate rapidity .0 ~ .85 for y1)
c 4/28/01  new D0 Jets included (several intervals in rapidity y1)
c 8/21/02  Fermilab Run II kinematics implemented (FirstF2)
c          UA2's Inty1 is changed to be 11.
      Data y1mnCDF, y1mxCDF, y1mnUA2, y1mxUA2 /.1D0, .7D0, 0D0, .85D0/ 
      Data y1mnD0, y1mxD0 /.0D0, .5D0/ 
      Data RsFermi, RsFermi2 /1800D0, 1960D0/ 
      Data FirstFermi, FirstUA2, FirstCDFdijet, FirstF2 
     >      / 4*.true./
      save FirstFermi, FirstUA2, FirstCDFdijet, FirstF2

      If (Irun.eq.0) 
     >Print '(A, 2I4, A)', 'ChiJet ... IB, ID = ', IB, ID, '....'

C                                        Control of NLO calculations:  


c                      Translating parameters from the fitting package 
C                      to this application & general initiation 
C                                                              >>---->>

c Scale
      Smucoef = ThJet(1) + dvJet
      DelYmax = ThJet(2)

      K_f    = LJet(1)

      NBEAM = JBEM(IB)
      NTAGT = JTGT(IB)
      IPDF = Iset

c  Inty1=0, integrate over y2 only
c  Inty1=1/2, CDF/D0 --integrate over y1 and y2 
c  Inty1=3, Fermilab RunI Rs=1800 inclusive jets in multi-rapidity intervals
c  Inty1=4, Fermilab RunII Rs=1960 inclusive jets in multi- rapidity intervals
c  Inty1=11, UA2 (changed 4/28/01)
c  If Inty1=-1, CDF--integrate over y1range(.1~.7) and y2range
c  If Inty1=-2, D0--integrate over y1range(.0~.5) and y2range

      Inty1 = Jsfn(IB)

      CHI = 0.0
C                                                  Process data points:
C                        Starting point rotates (for NLO k-fac purpose)
      IoffSt = Mod (Irun-1, Npt) 

      If(Inty1.ge.-2 .and. Inty1.le.3) then
         Rs=RsFermi
         Jset=1
      Elseif(Inty1.eq.4) then
         Rs=RsFermi2
         Jset=2
      Endif

      if(IorHrd.gt.1) then
            if(FirstFermi .and. (Inty1.ge.0 .and. Inty1.le.3)) then
               Call TrmStr(FlPath,  Len)
               If(ISch.eq.0) then
                  Kfile=FlPath(1:Len)//'kf1jet.msb'
               Else
                  Kfile=FlPath(1:Len)//'kf1jet.dis'
               Endif
                    Call TrmStr(Kfile, Len)
               Call Read1Jet(Kfile(1:Len))
               FirstFermi = .false.
            elseif(FirstF2 .and. Inty1.eq.4) then
               Call TrmStr(FlPath,  Len)
               If(ISch.eq.0) then
                  Kfile=FlPath(1:Len)//'kf1jet2.msb'
               Else
                  Kfile=FlPath(1:Len)//'kf1jet2.dis'
               Endif
               Call TrmStr(Kfile, Len)
               Call Read1Jet(Kfile(1:Len))
               FirstF2 = .false.
c            elseif(FirstUA2 .and. Inty1.eq.11) then
c               Kfile='kfac.ua2'
c               Call Read1Jet(Kfile)
c               FirstUA2 = .false.
c               FirstFermi = .true.
            elseif(FirstCDFdijet .and. Inty1.eq.-1) then
               Kfile='kfdijet.cdf'
               Call ReadDijet(Kfile)
               FirstCDFdijet = .false.
            endif
      Endif
      DO 30 Icnt = 1, Npt
C                                          Must use +1 so 1 < Ipt < Npt
         Ipt = Mod (Icnt, Npt) + 1
C                                                             >>---->>
C                                                   Kinematic variables
         If ((Inty1.eq.-1) .or. (Inty1.eq.-2)) then
               y2mn=ExDat(1,Ipt)
               y2mx=ExDat(2,Ipt)
               Yj=y2mn
         Elseif (Inty1.eq.3) then
               y1mn=ExDat(1,Ipt)
               y1mx=ExDat(2,Ipt)
               Yj=y1mn
         Elseif (Inty1.eq.4) then
               y1mn=ExDat(1,Ipt)
               y1mx=ExDat(2,Ipt)
               Yj=y1mn
         Else
            Rs     = ExDat(1,Ipt)
            Yj     = ExDat(2,Ipt)
         Endif 
         Pt     = ExDat(3,Ipt)
C                                                   Experimental data
         ExpDd = ExDat(4,Ipt)
C                                            ---- Always calculate LO
         Call ChiZz1 ('ChiJet', IorHrd, K_f, Icnt, NLO, Nchk)
C                                                  <<--- calculate LO
         if (Inty1.eq.1) then 
            call JetdEt( Rs,Pt,Smucoef,Born,y1mnCDF,y1mxCDF,0 )
         elseif (Inty1.eq.2) then 
            call JetdEt( Rs,Pt,Smucoef,Born,y1mnD0,y1mxD0,0 )
         elseif (Inty1.eq.3 .or. Inty1.eq.4) then 
            call JetdEt( Rs,Pt,Smucoef,Born,y1mn,y1mx,0 )
         elseif (Inty1.eq.11) then 
            call JetdEt( Rs,Pt,Smucoef,Born,y1mnUA2,y1mxUA2,0 )
c         elseif (Inty1.eq.0) then
c            call JetdEtdy( Yj, Rs, Pt, Smucoef, Born, 0)
         elseif (Inty1.eq.-1) then
            call DiJetdEt
     >           (Rs,Pt,Smucoef,Born,y1mnCDF,y1mxCDF,y2mn,y2mx,0 )
c         elseif (Inty1.eq.-2) then
c            call DiJetdEt
c     >           (Rs,Pt,Smucoef,Born,y1mnD0,y1mxD0,y2mn,y2mx,0 )
         endif

         if(nlo.eq.3) then
            if (Inty1.eq.1) then 
               call JetdEt( Rs,Pt,Smucoef,OneLop,y1mnCDF,y1mxCDF,1 )
            elseif (Inty1.eq.2) then 
               call JetdEt( Rs,Pt,Smucoef,OneLop,y1mnD0,y1mxD0,1 )
            elseif (Inty1.eq.3 .or. Inty1.eq.4) then 
               call JetdEt( Rs,Pt,Smucoef,OneLop,y1mn,y1mx,1 )
            elseif (Inty1.eq.11) then 
               call JetdEt( Rs,Pt,Smucoef,OneLop,y1mnUA2,y1mxUA2,1 )
c            elseif (Inty1.eq.0) then
c               call JetdEtdy( Yj, Rs, Pt, Smucoef, OneLop, 1)
            elseif (Inty1.eq.-1) then
               call DiJetdEt
     >              (Rs,Pt,Smucoef,OneLop,y1mnCDF,y1mxCDF,y2mn,y2mx,1 )
c            elseif (Inty1.eq.-2) then
c               call DiJetdEt
c     >              (Rs,Pt,Smucoef,OneLop,y1mnD0,y1mxD0,y2mn,y2mx,1 )
            endif
         endif

C                        ----- Evaluate K-fac according to the switches
C                               Update K-fac and calculate Theory value

         Call ChiZz2
     $        ('ChiJet',NLO,Nchk,Born,OneLop, ExDat(7,Ipt), FitDt)
C                                Save for the CDF/D0 correlation matrix calculation
         FitT(Ipt)=FitDt
C                        ------  Calc Naive ChiSqr by ChiZz3
         Call ChiZz3 (FitDt,ExpDd,ExErr(Ipt),DL2, LiOut(Ipt), Pt, Yj)

C                                Accumulate Chi^2
         Chi = Chi + DL2

   30 Continue
C			                          Fifth digit of LErr controls SysErr of Jet; 
C                                       If CorSysErr, then over-write naive chi^2
      Ksys = Mod(LErr/10**4,10)
      If (Ksys.ge.1) Then
        If (ID .EQ. 510) Then
           Chi = ChiD0JetIa (Npt,ExDat,FitT)
        elseif (ID .Eq. 521) Then
           Chi = ChiD0jetIb (Npt,ExDat,FitT)
        ElseIf (ID .EQ. 504) Then
           Chi = ChiCdfJet (Npt,ExDat,FitT)
        ElseIf (Ncor(IB).ge.1) Then               ! must contain CorSysErr info
           Chi = Chicorsys(Npt,ExDat,FitT,Fac,Ksys)
        EndIf
      EndIf

      ChiJet = Chi
C                                                    Output to .dta file
      IF (Ifg.EQ.3 .OR. Ifg.EQ.5 .or. Lfit.eq.0) Then
         Call ChiZz4 ('   Yj     ', '    Pt    ', LiOut,ExDat,chi)
      EndIf

      Return
C                        ****************************
      End
C
      Subroutine ChiZz0 
C                                                   -=-=- chizz0
     > (Ifg0,Irun0,IB0,Id0,Npt0,Lprt0,Fac0,IOdta,IOsum,IOtmp)
C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: ChiZzz
C Description: A set of general subroutines to handle the logistics of ChiSq
C ListOfFiles: chizz0 chizz1 chizz2 chizz3 chizz4
C ====================================================================

C #Header: /Disk2a/wkt/h2a/4fit/RCS/ChiZzz.f,v 65.2 99/05/19 16:15:12 wkt Exp $ 
C #Log:	ChiZzz.f,v $
c Revision 65.2  99/05/19  16:15:12  wkt
c new kfac written to .tmp, rather than .sum
c 
c Revision 65.1  99/05/19  10:54:24  wkt
c Synchronize to v65
c Cf. Log message in FcnEtc for v65 changes;
c Also, more changes made on Lprt use; Chizz4 now has one less argument.
c 

C      Printout contral switch  	-- Lprt -- 
C 
C *        Temporary debug : Lprt = 				99
C 
C    ----------- old   --------------------------    new   -----------
C 
C * chizz2 (Kfac calculation)
C   > 3 : print kfac info					-->   = 12
C 
C * chizz3 :
C 
C   >=3 : print chisq for every 10 data points 		-->   = 13/14
C   >=4 : every point							= 14
C 
C ---------------------------------------------------------------------
C 
c Revision 6.2  98/12/23  23:25:27  wkt
c mislabel "NormExp" corrected
c 
c Revision 6.1  98/08/05  16:55:03  wkt
c Argument IUtmp added to ChiZz0; it is used in ChiZz4.
c 
c Revision 6.0  98/03/08  14:03:30  wkt
c New start from Liang's "CTEQ4" package dated September 1997.
c 
C ---------------------------------------------------------------------

c Revision 1.1  93/08/20  12:18:14  lai
c Initial revision
c 
c Revision 1.9  92/09/05  21:13:34  botts
c Htarget better
c 
c Revision 3.4  92/04/19  23:25:07  wkt
c Commons ThXxxx and ExXxx reorganized due to re-ordering of vbp Hhk..
c and other minor improvements.
c 
c Revision 3.5  92/04/19  23:11:11  wkt
c chizz4 call modified.
c 
c Revision 3.4  92/03/31  23:56:47  wkt
c ??
c 
c Revision 3.3  92/03/28  18:05:47  wkt
c chizz0 moved to fcn; N10 changed; Iset added
c 
c Revision 3.2  92/03/26  12:43:35  wkt
c Minor fine tuning
c 

C     This package of subroutines supports the calculation of Chi-squares for 
C     the various applications processes; allowing the corresponding applications
C     modules "ChiXxx" to be simple and transparent.

C                                         Initialization for the various entries
      Implicit Double Precision(A-H, O-Z)

      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lpnt,IUdta,IUsum,IUtmp

C                Make these common variables local and available for all entries
      Fac = Fac0
      Ifg = Ifg0
      Irun= Irun0
      IB = IB0
      ID  = Id0
      Npt = Npt0
      Lpnt= Lprt0
      IUdta = IOdta
      IUsum = IOsum
      IUtmp = IOtmp

      Return
C                                ===================
      End

C                                                   -=-=- chizz1
      Subroutine ChiZz1 (AppNam, IorHrd, K_f, Icnt, NLO, Nchk)

C     Input  variables: AppNam, IorHrd, K_f, Icnt

C     Output variables: NLO, Nchk

      Implicit Double Precision(A-H, O-Z)
      Character*(*) AppNam

      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lpnt,IUdta,IUsum,IUtmp

C                                      Initialize NLO and Check switches to 0
      NLO  = 0
      Nchk = 0

      If (IorHrd .eq. 2) Then
C                                           -----   Determine how to do NLO
        If     (K_f.Ge.1 .and. K_f.Le.10) Then
            KK = Mod(Icnt, Npt/K_f)
            If (KK .Eq. 1) Then
C                                       Do true calc at K_f evenly distributed 
C                                       points starting from the "entry point" 
               NLO = 3
               If (Icnt .Ne. 1)  Nchk = 1
            Else
C                        Use "ratio of ratio method" (see below) for other pts
               NLO = 2
            Endif
         Elseif (K_f.Eq.99) Then
C             For K_f = 99 calculate K-fac for all pts during the 1st run only
C                          use the calculated K-fac for the rest of fitting
              NLO = 1
              IF(IFG.EQ.1) NLO = 3
C                                      For K_f = 100 use true calc for all pts
        Elseif (K_f.Eq.100) Then

            NLO = 3
        Else
C                     For all other cases, use default K-fac method for all pts
            NLO = 1
            If (K_f .Ne. 0) Then
              Print *, 'Illegal value of K_f in ', AppNam,': K_f =',K_f
              Print *, 'Defaulting to K_f=0 (previous value) method.'
            Endif
        Endif

      Elseif (IorHrd .Ne. 1) Then
       Print *,'Illegal value of IorHrd in ', AppNam,'! IorHrd=',IorHrd

      Endif

      Return
C                                ================
      End

C                                                   -=-=- chizz2
      Subroutine ChiZz2 (AppNam, NLO,Nchk, Born,OneLp, OldKNew, FitDt)

C        Input  Variables: AppNam, NLO, Nchk, Born, OneLp, OldKNew

C        Output variables: OldKNew, FitDt

      Implicit Double Precision(A-H, O-Z)
      Character*(*) AppNam

      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lpnt,IUdta,IUsum,IUtmp

C                    ----- Evaluate and Update K-fac according to the switches
C                                              NLO = 0:  by definition K-fac=1
      If     (NLO .Eq. 0) Then
            FacK = 1.0
C                                              NLO = 1:  Use previous value 
      ElseIf (NLO .Eq. 1) Then

         FacK = OldKNew
C                                              NLO = 2: use ratio/ratio method 
      Elseif (NLO .Eq. 2) Then

C                               Shift Registers:          |Ref Run |Curr Run
C                                              -----------------------------
C                                              Ref Point  |  FK00     FK01
C                                              Cur Point  |  FK10     FK11

C                                            Use previous run value for FK10:
         FK10 = OldKNew
C                                  For now, Reference values of FK00 and FK01 
C                                             are set by Icnt = 1 calculation
C                                                       
C                                                     "Best estimate" of K-fac
         FK11 = FK01 * FK10 / FK00

C  ----               only useful if ratio is different from that of Icnt = 1:
C                                            Shift "current point" values into 
C                                    "reference point" registers for use later
C              FK00 = FK10
C              FK01 = FK11
C -----

         FacK = FK11

      Elseif (NLO .Eq. 3) Then
C                                        NLO = 3:  use true calculation results
         Onelop = Onelp  
C                                        K-fac from previous run for the record
         Fk0 = OldKNew
C                                                       K-fac from current run 
         Fk1 = Onelop / Born
C                          These two are kept as reference until next true calc.
         If (Ifg.Eq.1 .and. Nchk .Eq. 0) Then
          If (Lpnt .Eq. 12) Write (IUsum, '(2A/ A40, A30)') 
     >   ' K-fac info for ', AppNam, ' Entry Point:  DelK, FK0, FK1'
     >                             , 'Other Points: DelK, FK11, FK1'
         Endif
      If (Nchk.Eq.0 .and. Lpnt.Eq.12) Then
           DelK = (Fk1-Fk0)/(Fk1+Fk0) * 2
           Print '(A, F10.6, 2(1pE10.3))', 
     >    ' K-fac Calc: DelK, Fk0, Fk1 = ',     DelK, Fk0, Fk1
           Write (IUsum, '(F10.6, 2(1pE10.3))') DelK, Fk0, Fk1
        Endif

C                                 -------- This section performs a check on the 
C                                            accuracy of the ratio/ratio method.
C                               Has to be here - before (FK01 FK00) are updated
      If (Nchk .Eq. 1) Then
        FK10 = Fk0 
        FK11 = FK10 * FK01 / FK00
        DelK = (Fk1 - FK11)/ (Fk1 + FK11) * 2.
        DelK = Abs (DelK)
C                                                   DelKmx = Max value of error
        If (DelK .GT. DelKmx) Then
           IrunMx = Irun
           DelKmx = DelK
        Endif

        If (Lpnt .Eq. 12) Then
           Print '(25x, 2A, 5X, F10.6, 2(1pE10.3))', 
     >     AppNam, ' check : ' ,DelK, FK11, Fk1
           Write (IUsum, '(40x, F10.6, 2(1pE10.3))') DelK, FK11, Fk1
        Endif

C                                                IwarnK = # of pts > 0.01 error
        If (DelK .Ge. 0.02) Then 
           IwarnK = IwarnK + 1
           Print *,' K-fac method Warning: IwarnK, Irun =' 
     >         ,IwarnK, Irun
        Endif
C                                     Put Summary of K-fac check info to Output
        If (Ifg .eq. 3 .or. Ifg .eq. 5) Then
           Write (IUsum, '(3A, F10.6)') 
     >      ' Maximum K-fac % error in ', AppNam, ' = ', DelKmx
     >      ,'in run # ', IrunMx
           Print '(3A, F10.6, A, I5)',
     >      ' Maximum K-fac % error in ', AppNam, ' = ', DelKmx
     >      ,'in run # ', IrunMx
        Endif
C                                                             Reset to no Check
        Nchk = 0
      Endif
C                                 -------- End of check section

C                         Use new results to update "reference" shift registers
         FK00 = Fk0
         FK01 = Fk1

         FacK  = Fk1
C                                                NLO = none Of of above == error
      Else
         Print *, 'Illegal value for NLO in ',AppNam,': NLO=', NLO
         Print *, 'Defaulting to NLO=3 (previous K-fac) method.'
         FacK = OldKNew
      Endif
C                                                         ---- end of K-fac calc.
C                                                            Update K-fac except
      OldKNew = Fack

      FitDt = Born * FacK

      Return
C                         =======================
      End

C                                                   -=-=- chizz3
      Subroutine ChiZz3 (FitDt,ExpDd,ExErr, DL2, LiOut, V1, V2)
C
C     (i) computes the point-to-point \Delta_chi^2
C    (ii) saves LiOut which lists kinematic variables, Theory, expt, Errs, ... etc. for output;
C   (iii) To leave exptl value intact, Theory is divided by the normalization factor.
C                       
C
      Implicit Double Precision(A-H, O-Z)
      Character*(*) LiOut

      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lpnt,IUdta,IUsum,IUtmp
     >  /FitSwh/ Lfit, LErr, Lprt, Loutput
C                                                           Chi-Sq calculation
      ExpEr = ExErr
      ExpDt = ExpDd * Fac
      ExpEr = ExpEr * Fac
C                                                          compute ChiSq
      DEL = (ExpDt - FitDt) / ExpEr 
      DL2 = DEL **2
      if(del.le.0.0) then
         dljim = dl2
      else
         dljim = -dl2
      endif
C                                                   Save results for ChiZz4
      IF (Ifg.EQ.3 .OR. Ifg.EQ.5 .or. Lfit.eq.0) Then
C
      ExpJDt = ExpDd 
      ExpJer = ExErr
      FitJDt = FitDt / Fac
C                                                            Expt / Theory
      R1 = ExpDt/FitDt
      R2 = Abs(ExpEr/FitDt)

901   Format (5(1pE11.3), 3(0pF9.3))
902   Format (5(1pE11.3), 3(1pE9.1))
      If (dl2 .LE. 999.0) Then
        Write (LiOut, 901, Err=905) 
     >  V1,V2, ExpJDt,FitJDt,ExpJEr, R1,R2, dljim
      Else
        Write (LiOut, 902, Err=905) 
     >  V1,V2, ExpJDt,FitJDt,ExpJEr, R1,R2, dljim
        Print *, 'Large chi^2 warning: DelChi^2 =', dljim
      EndIf
C
      Endif

      Return
  905 Print '(A/A, 3I8/A/6G12.3)', 
     > 'Error writing to Output line in ChiZz3:'
     > , 'Irun, IB, ID =', Irun, IB, ID
     > , 'V1, V2, ExpJDt, FitJDt, ExpJEr, Del ='
     > ,  V1, V2, ExpJDt, FitJDt, ExpJEr, Del
      stop
C                             ==================
      End

C                                                   -=-=- chizz4
      Subroutine ChiZz4 (Var1Nm, Var2Nm, LiIn, ExDat, chi) 

C $Header: /group/cteqX/CVS/fitpack/theory/jets.f,v 1.3 2011/08/02 21:17:22 accardi Exp $
C $Log: jets.f,v $
C Revision 1.3  2011/08/02 21:17:22  accardi
C *** empty log message ***
C
C Revision 1.2  2011/03/15 18:24:03  jowens
C modified pdf call --> fsupdf
C
C Revision 1.1.1.4  2010/03/12 20:06:09  accardi
C NextUn function from jet.f put in util/io.f
C
C Revision 1.1.1.3  2010/01/26 22:53:07  jowens
C New files for the updated fitting and evolution packages
C
C Revision 1.1.1.2  2009/04/02 17:32:24  accardi
C New data sets, improvements in 'altpfit',  3D deuterium convolution, off-shell corrections, April_09 fitting subfolder
C
C Revision 1.1.1.1  2008/08/12 16:54:21  accardi
C Jeff Owen's global QCD fit package - May 2008
C
C Revision 1.1.1.1  2008/08/11 14:52:11  accardi
C Jeff Owen's global QCD fitting package - May 2008
C
C Revision 8.7  2003/07/23 03:46:53  wkt
C output .dta file only on demand.
C
C Revision 8.6  2003/07/23 03:09:44  wkt
C *** empty log message ***
C
      Implicit Double Precision(A-H, O-Z)
      
      Parameter (one = 1.0, Zero = 0.)
      Parameter (Maxpt = 500)
      
      Character LiIn*(*), Var1Nm*10, Var2Nm*10
      Logical Ldta, Ltmp, Lpds, Lplt, Lhip

      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lpnt,IUdta,IUsum,IUtmp
     >  /Loutpt/ Ldta, Ltmp, Lpds, Lplt
     
      Dimension LiIn(*), ExDat(9, *)

      Logical Lsys
      Character LiRval*132, LiCorSout*50, Line*120
      Common /CorSysOut/ Lsys, LiRval, LiCorSout(Maxpt)

C                                    >>   Write final results to .dta file
      chi2 = chi
      If (Ldta) then
C                                                                   Header:
      Write (IUdta,'(/A, I5, A, F6.3, A, I5, A, F10.2)') 
     >  ' DATA SET:', ID, ' ;  NORM Fac =', Fac, 
     >  ' ;  # of pts =', Npt, '  ; chi^2 =', chi2
      If (Lsys) Write (IUdta, '(A)') LiRval
      Write (IUdta, '(5x,4A)') Var1Nm, Var2Nm
     >  , '  Exp      Th./Norm    TOT ERR  EXP/FIT ERR/FIT  ChiSq'
     >  , '  Shift   ShiftedData  UnCorErr  ReducedChi2'

C                           Simplify the previous version (TopPlot format) to trivial case
      Do 11 Ipt = 1, Npt
        Line = LiIn(Ipt)
        Call TrmStr(Line,Len)
        Read (Line, *, Err=905) V1, V2, ExDt, ThDt, ExEr, R1, R2, Ch2
        If (Lsys) then
          Write (IUdta, '(A,3x,A)') Line(1:Len), LiCorSout(Ipt)
        Else
          Write(IUdta,'(A, 1pE11.3, 2E11.3, 0pF8.2)')
     >    Line(1:Len), zero, ExDt, ExEr, Ch2
        Endif
   11 Continue
      
      EndIf

      Lsys = .false.
      Return
      
901   Format (5(1pE11.3), 3(0pF10.3))
905   Print '(A/ I8, A)',
     >'Error reading variables from Line in ChiZz4; Ipt, Line ='
     >, Ipt, Line
      Stop
C                           *********************
      End
      Function ChiCDFjet(Mpt,ExDat,Fitjet)
C                                                   -=-=- chicdfjet
C----------------------------------------------------------------------
C    Corrected for mistake in initial version. From Huston 2/5/99

c    This routine is set up to copy what Liang did for the D0 jet
c     J. Huston 12/21/98
c
c    The unit in cross section is in nb/GeV. 
c    Therefore, the unit convertion constant "ConUnit" is set to be 1D00.
c
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Evaluate Chi^2 relative to CDF 0.1<|eta|<0.7data
C-
C-   Inputs  : THEORY R(33) - array of CS values (in pb) evaluated at the
C-             following ET values (in GeV):
C-
C-   Note: this CDF jet cross section has 33 bins
C-
C-   43.3,49.3,55.2,61.0,66.7,72.3,77.9,83.5,89.0,94.5,100.0,105.4,
C-   110.9,116.3,121.7,127.1,132.5,137.9,145.7,156.4,167.2,177.9,188.7,
C-   199.5,210.2,225.4,247.1,268.8,290.5,312.1,333.6,362.2,412.9
C-
C-   Outputs :
C-   Controls:
C-
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      Integer Mpt
      Double Precision ChiCDFjet, ExDat(9,Mpt), Fitjet(*)

      Integer Mblk, MaxjetDt
      Parameter (Mblk = 50) 
      Integer Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp
      Double Precision Fac
      Character NmDt*10, FlPath*30, FlExplis*15, Flmatx*50
      Common /FlNams/ FlPath, FlExplis, NmDt(Mblk)
      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp

      Integer IU, NextUn

      PARAMETER ( MaxjetDt = 50 )
      INTEGER i,j,k,n, Len
      Double Precision CovInV(MaxjetDt,Maxjetdt),c,e
      double precision fitjdt(MAXJETDT)
      Double Precision ConvUnit
      Double Precision chi2
      LOGICAL first
      DATA first /.true./
      save CovInV
C----------------------------------------------------------------------
      IF (FIRST) THEN
        FIRST=.FALSE.
        ConvUnit=1D00
        IU=NextUn()
        Call TrmStr(FlPath,  Len)
        Flmatx= FlPath(1:Len)//'cdf_01_07_covmtx.dat'
        OPEN(IU,file=Flmatx,status='old',err=999)
        READ(IU,*) n
        If(n.ne.Mpt) stop '# of jet data does not match Error Matrix.'
        DO k = 1,n*n
          READ (IU,*) i,j,c,e
          CovInv(i,j) = e * ConvUnit
        ENDDO
        CLOSE(IU)
      ENDIF

c     note that the CDF covariance matrix is stored in a different way
c     than that of D0
C      chi2= 0D0
C      DO j = 1,Mpt
C        Fitjet(j)=Fitjet(j)/Fac
C        DO k = 1,Mpt
C             chi2 = chi2 + (Fitjet(j) - ExDat(4,j)) * CovInv(j,k) *
C     &          (Fitjet(k) - ExDat(4,k))
C        ENDDO
C      ENDDO
C      ChiCDFjet=chi2

      do j = 1,mpt
         fitjdt(j) = fitjet(j) / fac
      end do


      chi2= 0.0D0
      DO j = 1,Mpt
         DO k = 1,Mpt
            chi2 = chi2 + (fitjdt(j) - ExDat(4,j))*CovInv(j,k)*
     &             (fitjdt(k) - ExDat(4,k))
         ENDDO
      ENDDO

      ChiCDFjet = chi2

      RETURN
      
 999  Print*,'Cannot open file:',Flmatx
      stop
C                      *******************************       
      END
C                                                          =-=-= ChiZzz
      Function ChiCorSys(Mpt,ExDat,FitT,FacNor,Ksys)
C----------------------------------------------------------------------
C      Translation between local variables and symbols in CTEQ6 paper, Appendix B
C     \sigma_i    :  ExDat(5,i)
C      D_i        :  ExDat(Idata,i) (Idata=3/4 for DIS/Jet)
C      T_i        :  FitT(i) / FacNor
C      u_i        :  Uncor(i)
C      \alpha_i^2 :  ErrUnc2(i)
C      \beta_ki   :  sys(k,i) * [FitT(i) | ExDat(Idata,i)]
C      A_k,k'     :  Akk(k,k')

      Implicit Double Precision (A-H, O-Z)

      Parameter (MxDt = 5000, Mblk = 50)            !unused ones removed
      Parameter (MxSys = 20, Maxpt = 500)

      Logical Lsys
      Character LiRval*132, LiCorSout*50

      Common 
     >  /FcnHrd/ IB, ID, Ifg, N10, Irun
     >  /FitSwh/ Lfit, LErr, Lprt, Loutput
     >  /CorMtx/sys(MxSys,MxDt),Uncor(MxDt),Nptcor(Mblk),Ncor(Mblk)
     >  /ExpSys/ Rik(MxSys, Mblk), Nsys(MxSys), LRik(MxSys, Mblk)
      Common /CorSysOut/ Lsys, LiRval, LiCorSout(Maxpt)

      Dimension Diff(Maxpt),ErrUnc2(Maxpt)
      Dimension Akk(MxSys,MxSys),Akkinv(MxSys,MxSys)
      Dimension Bk(MxSys), r(MxSys)
      Dimension ExDat(9,*),FitT(*)


C----------------------------------------------------------------------
      Ncorr=Ncor(IB)
      Itype=ID/100
      If(Itype.eq.5) then
         Idata=4
      Else
         Idata=3
      Endif

      Do i = 1,Mpt
         ErrUnc2(i) = Uncor(i+Nptcor(IB))**2 + ExDat(5,i)**2
         Diff(i) = (ExDat(Idata,i)-fitT(i)/FacNor)/ErrUnc2(i)
      Enddo

      Do k=1,Ncorr
      Do j=1,Ncorr
        If(j.eq.k) then
           Akk(j,k)=1D0
        Else
           Akk(j,k)=0D0
        Endif
        Do i=1, Mpt
           If     (Ksys.eq.1 .or. Ksys.eq.3) then
             Akk(j,k)=Akk(j,k)+   fitT(i)**2*
     >         sys(j,i+Nptcor(IB))*sys(k,i+Nptcor(IB))/ErrUnc2(i)
           ElseIf (Ksys.eq.2 .or. Ksys.eq.4) then
             Akk(j,k)=Akk(j,k)+ExDat(Idata,i)**2*
     >         sys(j,i+Nptcor(IB))*sys(k,i+Nptcor(IB))/ErrUnc2(i)
	     Else
           Stop 'ChiCorSys: Ksys must be 1-4 in LErr specification!'
           Endif
        Enddo
      Enddo
      Enddo

        Call invmatrix(Akk,Ncorr,MxSys,Akkinv)

      Do k=1,Ncorr
         Bk(k)=0D0
         Do i=1,Mpt
C                   Calculation of Akk already checked that Ksys = 1/2/3/4
           If (Ksys.eq.1 .or. Ksys.eq.3) then
              Bk(k)=Bk(k)+   fitT(i)*sys(k,i+Nptcor(IB))*Diff(i)
           Else
              Bk(k)=Bk(k)+ExDat(Idata,i)*sys(k,i+Nptcor(IB))*Diff(i)
           Endif
         Enddo
      Enddo
C                  Evaluate chi^2 due to uncorrelated errors: 1st term of Eq.(10)
      chi2= 0D0
      DO j = 1,Mpt
         chi2 = chi2 + (ExDat(Idata,j)-fitT(j)/FacNor)*Diff(j)
      ENDDO
C                  Now evaluate chi^2 contributions due to correlated errors 
      Do j=1,Ncorr
C      Evaluate the "r(k)" vector that gives the optimal systematic error shifts
         r(j) = 0D0
         Do k=1,Ncorr
            r(j) = r(j) + Akkinv(j,k)*Bk(k)
         Enddo
C                                             Save in /ExpSys/ for output to .new
         Rik(j,IB)= r(j)
C                              Add contributions to chi^2 due to correlated errors
         chi2=chi2-Bk(j)*r(j)
      Enddo
C                              Compute and store 
C                              (1) the optimal shifts of the data points;
C                              (2) the shifted data points (for plotting purpose);
C                              (3) point-to-point contribution to overal chi^2
C                                  as the results of correlated systematic errors
     
C                    Overall contribution to chi^2 due to non-vanishing r(k)
      r2 = 0
      Do k=1,Ncorr
            r2 = r2 + r(k)**2
      Enddo
C                                             Save for output to .dta
      Write (LiRval, '(A, 13F8.3)')
     >  'R^2, r(k) = ', r2, (r(k), k=1,Ncorr)                                                ! LiRval

      Chi22 = r2          ! Calculate Chi^2 according to the alternate formula, Eq.(7)
      Do i=1,Mpt
         Dd = 0           ! Shift of data point due to SysErr
         Do k = 1, Ncorr
            If(Mod(Lsw3,4).eq.1) then
              Dd = Dd + r(k) * ExDat(Idata,i)*sys(k,i+Nptcor(IB))        ! cf. Eq.(7)
            Else
              Dd = Dd + r(k) *    FitT(i)*sys(k,i+Nptcor(IB))
            Endif
         Enddo
         Ds = ExDat(Idata,i) - Dd              ! Shifted experimental value                                                    ! 'Shifted Data'
         Dchi  = Ds - FitT(i)/FacNor
         Dchi2 = Dchi **2 / ErrUnc2(i)                               !  cf. Eq.(7)
         Chi22 = Chi22 + Dchi2
C                               Write out signed chi^2; + if Ex < Th; - if Ex > Th
         if(Dchi .le. 0.0) then
           Chi2n = Dchi2
         else
           Chi2n =-Dchi2
         endif
         Write (LiCorSout(i), '(1pE11.3, 2E11.3, 0pF8.2)')
     >   Dd, Ds, Sqrt(ErrUnc2(i)), Chi2n                                                                              ! LiCorSout

      Enddo
C                        Consistency check between the two equivalent calculations
      D2 = (chi2 - chi22) / (chi2 + chi22)
      If (D2**2 .Ge. 1D-4) Print '(A, 1pE11.3)', 
     >'Warning: SysErr calc. in ChiDisCor is not consistent: D2 =', D2
     
C                              
      If (Ksys.LE.2) then
         ChiCorSys = chi2      ! Chi^2 calculated by Eq.(10) of CTEQ6, App.B
      Else
         ChiCorSys = chi22     ! Chi^2 calculated by Eq.(7) of CTEQ6, App.B
      Endif

      Lsys = .true.

      RETURN      
C                      *******************************       
      END      
      Function ChiD0jetIa (Mpt,ExDat,Fitjet)
C                                                   -=-=- chid0jet
C----------------------------------------------------------------------
c    Modified by H.L. Lai 9/2/98, 9/15/98
c
c    The unit in cross section is in nb/GeV. 
c    Therefore, the unit convertion constant "ConUnit" is set to be 1D12.
c
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Evaluate Chi^2 relative to D0 |eta|<0.5 data
C-
C-   Inputs  : THEORY R(24) - array of CS values (in pb) evaluated at the
C-             following ET values (in GeV):
C-
C-   64.60, 74.64  84.68, 94.71, 104.7, 114.8, 124.8, 134.8
C-  144.8, 154.8, 164.8, 174.8,  184.8, 194.8, 204.8, 214.8
C-  224.8, 239.4, 259.4, 279.5,  303.9, 333.9, 375.7, 461.1
C-
C-   Outputs :
C-   Controls:
C-
C-   Created  13-AUG-1998   Bob Hirosky
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      Integer Mpt
      Double Precision ExDat(9,Mpt),Fitjet(*),ChiD0jetIa

      Integer Mblk, MaxjetDt
      Parameter (Mblk = 50) 
      Integer Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp
      Double Precision Fac
      Character NmDt*10, FlPath*30, FlExplis*15, Flmatx*50
      Common /FlNams/ FlPath, FlExplis, NmDt(Mblk)
      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp

      Integer IU, NextUn

      PARAMETER ( MaxjetDt = 50 )
      INTEGER i,j,k,n, Len
      Double Precision CovInV(MaxjetDt,Maxjetdt),c,e
      Double Precision ConvUnit
      Double Precision chi2
      LOGICAL first
      DATA first /.true./
      save CovInV
C----------------------------------------------------------------------
      IF (FIRST) THEN
        FIRST=.FALSE.
        ConvUnit=1D12
        IU=NextUn()
        Call TrmStr(FlPath,  Len)
        Flmatx= FlPath(1:Len)//'d0_00_05_covmtx.dat'
        OPEN(IU,file=Flmatx,status='old',err=999)
        READ(IU,*) n
        If(n.ne.Mpt) stop '# of jet data does not match Error Matrix.'
        DO k = 1,n*n
          READ (IU,*) i,j,c,e
          CovInv(i,j) = e * ConvUnit
        ENDDO
        CLOSE(IU)
      ENDIF
C
      chi2= 0D0
      DO j = 1,Mpt
        Fitjet(j)=Fitjet(j)/Fac
        DO k = j,1,-1
          if (j.eq.k) then
             chi2= chi2 + (Fitjet(j) - ExDat(4,j))**2 * CovInv(j,j)
          else
             chi2 = chi2 + 2D0* (Fitjet(j) - ExDat(4,j)) * CovInv(j,k) *
     &          (Fitjet(k) - ExDat(4,k))
          endif
        ENDDO
      ENDDO
      ChiD0jetIa = Chi2
c
      RETURN
      
 999  Print*,'Cannot open file:',Flmatx
      stop
C                      *******************************       
      END
      
      Function ChiD0jetIb(Mpt,ExDat,Fitjet)
C                                                   -=-=- chid0jet
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Evaluate Chi^2 of D0 RunIb data
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
      Integer Mpt
      Double Precision ExDat(9,Mpt),Fitjet(*),ChiD0jetIb

      Integer Mblk, MaxDt
      Parameter (Mblk = 50) 
      Integer Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp
      Double Precision Fac, Tmp, Diag, Diaf
      Character NmDt*10, FlPath*30, FlExplis*15, Flmatx*50, Lin80*80
      Common /FlNams/ FlPath, FlExplis, NmDt(Mblk)
      Common /ChiZzz/ Fac, Ifg,Irun,IB,ID,Npt,Lprt,IUdta,IUsum,IUtmp

      Integer IU, NextUn, IU1

      PARAMETER ( MaxDt = 100 )
      INTEGER i,j,k,l,n, Len
      Double Precision CoVar(MaxDt,MaxDt), Beta(MaxDt), Gamma(MaxDt)
     >  , ErrMtx(MaxDt,MaxDt)
      Double Precision chi2, ConvUnit, c
      LOGICAL first
      DATA first /.true. /
      save ErrMtx
C----------------------------------------------------------------------
      IF (FIRST) THEN
        FIRST=.FALSE.
C                                           ConvUnit : from fb (data file) to nb (theory x-sec)
        ConvUnit=1D-12
        IU=NextUn()
        Call TrmStr(FlPath,  Len)
        Flmatx= FlPath(1:Len)//'D0CovIb.txt'
        OPEN(IU,file=Flmatx,status='old',err=999)
        IU1 = NextUn()
        Open(IU1, file='tmp.dat')
        Read (IU, '(A80)') Lin80
        Read (IU, *) n
        If(n.ne.Mpt) 
     >  stop '# of jet data does not match Error Matrix dimension.'
        DO k = 1,n*n
          Read (IU, '(A80)') Lin80
          READ (Lin80, *) i,j,c
C                      If (i.ne.90 .and. j.ne.90) Write(IU1, '(A)') Lin80
          tmp = ConvUnit *c     
          If (i .Eq. j) then
              diag = tmp - ExDat(5,i) **2
C                                                                                                                        Fractional error
              diaf = diag / ExDat(4,i)**2
              If (diaf .Ge. -2D-6) then
                beta(i) = Sqrt(Abs(diag))
                gamma(i)= Sqrt(Abs(diaf))
              Else 
                pause 'diagonal element non-positive'
              Endif
          Endif
          CoVar(i,j) =  Fitjet(i) * Fitjet(j)
     >                 * tmp /ExDat(4,i)  /ExDat(4,j)
        ENDDO
        CLOSE(IU)
        Write (IU1, '(I5, 2E15.5)') (i, beta(i), gamma(i), i=1,Mpt)
        
C                        Invert this matrix to get the Error matrix
        Print *, 'Inverting the D0 covariance matrix'

        Call InvMatrix(CoVar, N, MaxDt, ErrMtx)
C                CALL DLINRG (N, CoVar, MaxDt, ErrMtx, MaxDt)   ! ISML inversion program
        Print *, 'Done with inverting ... '
C                                             Check accuracy of inversion                                                                                                 Check accuracy of Inverse        
c      DO j = 1,Mpt
c        DO k = j,Mpt
c           E(j,k) = 0.
c           Do l = 1, Mpt
c             E(j,k) = E(j,k) + CoVar(j,l) * ErrMtx(l,k)
c           EndDo
c        ENDDO
c      ENDDO
C                                                                                                                        End of check
      ENDIF
C
      chi2= 0D0
      DO j = 1,Mpt
        Fitjet(j)=Fitjet(j)/Fac
        DO k = j,1,-1
          if (j.eq.k) then
             chi2= chi2 + (Fitjet(j) - ExDat(4,j))**2 * ErrMtx(j,j)
          else
             chi2 = chi2 + 2D0* (Fitjet(j) - ExDat(4,j)) * ErrMtx(j,k) *
     &          (Fitjet(k) - ExDat(4,k))
          endif
        ENDDO
      ENDDO
      ChiD0jetIb = Chi2
c
      RETURN
      
 999  Print*,'Cannot open file:',Flmatx
      stop
C                      *******************************       
      END
      
      Subroutine SetQCD
C                                                   -=-=- setqcd
C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: Setqcd
C Description: Set up the qcdpac of programs, initiate common blocks
C ListOfFiles: setqcd
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/2qcd/RCS/Setqcd.f,v 1.1 97/12/21 20:35:00 wkt Exp $ 
C #Log:	Setqcd.f,v $
c Revision 1.1  97/12/21  20:35:00  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      External DatQCD

      Dummy = 0.

      Return
C                        ****************************
      END

      BLOCK DATA DATQCD
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      COMMON /COMQCH/ VALQCH(9)
      COMMON /COMQMS/ VALQMS(9)
      COMMON /QCDPAR/ AL, NF, NORDER, SET
      COMMON /COMALP/ ALPHA
      LOGICAL SET
C
      DATA AL, NF, NORDER, SET / .226, 5, 2, .FALSE. /
      DATA VALQCH/ 0.66666667, -0.33333333,
     >  -0.33333333, 0.66666667,
     >  -0.33333333, 0.66666667,
     >  3*0./
      DATA VALQMS/  2*0., 0.2, 1.3, 4.5, 174., 3*0./
      DATA ALPHA/  7.29927E-3 /
 
C                       ******************************
      END

C                                                          =-=-= Alphas
