C **********************************************************************
      SUBROUTINE DEUTERONSF (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,Fd)
C
C  Compute 0.5* deuteron structure functions (Fd) per nucleon
C  at x (=v(1)) and Q^2 (=v(2)) by smearing nucleon structure function
C  (array F2Narr dimension nx defined over x array xarr) with nucleon
C  momentum distribution. Optionally smears the proton's PDFs at gamma=1
C  to calculate nuclear PDFs.
C
C    ibeam             = (i) 1=e 2=nu 3=nubar
C    istruc            = (i) 0=FL  1=F1  2=F2  3=x*F3
C    inuke             = (i) nuclear corrections: ABCDEF
C                            A: TMC in smearing functions:
C                                   gamma^2=1+4*xN^2*Lambda^2/Q^2
C                               0: Lambda^2=M^2
C                               1: Lambda^2=0  ("Bjorken limit", no TMC)
C                            B: Shadowing correction
C                               0 = no
C                               1 = Melnitchouk-Thomas
C                            C: Strength of off-shell corrections
C                               mKP: 1-9 = 0.3%-2.7% in steps of 0.3%
C                                    0 = 3%  (read "10" instead of "0")
C                               others: inactive flag
C                            D: Off-shell corrections
C                               0 =  on-shell nucleons
C                               1 =  KP - Kulagin-Petti fit, nucleon level
C                                    (OUTDATED --- use D=3 instead)
C                               2 =  MST - Melnitchouk-Thomas model,
C                                    parametrization at Deuterium level
C                               3 =  KP par - Kulagin-Petti fit,
C                                    parametrization at Deuterium level
C                              >4 =  mKP    - modified Kulagin-Petti model, c
C                                    parametrization at Deuterium level
C                            E: Deuteron wave-function
C                               0=Paris
C                               1=AV18
C                               2=CDBonn
C                               3=WJC1
C                               4=WJC2
C                            F: Deuteron correction model
C                               0=isospin average
C                               1=density model
C                               2= --- not supported anymore ----
C                               3=nuclear smearing WBA
C                               4=nuclear smearing WBAREL (WBA relativized)
C                               5=nuclear smearing AQV  (M/p0 norm)
C                               6=nuclear smearing AQVc (const. norm)
C    itm,iht,ndrv,v,xc = see 'strfn' routine
C    Fd                = deuteron structure function per nucleon
C
C  This version compues F2(xN) as needed iside the convolution
C  integral, instead of precomputing it and interpolating an array.
C  (for the latter, see deuteronSF_pre.f)
C
C **********************************************************************

      IMPLICIT NONE
      INTEGER	ibeam,inuke,istruc,itm,iht,ndrv
     &     ,idum0,idum1,idum2,idum3,idum4
      double precision v(4),xc(100),Fd
      INTEGER nx,ix, ismear
      PARAMETER (nx=100)
      REAL*8	x,Q2,varr(4),xarr(nx),FNarr(nx),dsf,dmc

*    *** Functions
      double precision gomez

*    *** common blocks

*     Proton, neutron fraction    
      double precision az,an
      common/target/az,an
      
      az=0.5D0                  ! p + n isoscalar, i.e., (p+n)/2
      an=0.5D0

      ! selects smearing model
      call split_nuke(inuke,idum0,idum1,idum2,idum3,idum4,ismear) 
                     
      if (ismear.le.1) then
*    *** p+n isoscalar 
         CALL strfn (ibeam,istruc,itm,iht,ndrv,v,xc,Fd)
         if (ismear.eq.1) then
*       *** Density model corrections (for backward compatibility)
            dmc=gomez(1d0,v(1))
            Fd = Fd/dmc
         end if
      else if ((ismear.ge.2).and.(ismear.le.9)) then
*       *** Nuclear smearing model
C       ...Compute nucleon str.fn.to be (interpolated and) smeared
*         varr(2) = v(2)	! store Q^2 value in new array
*	  DO ix=1,nx
*	     x = DFLOAT(ix)/100.D0 - 1d-4
*	     xarr(ix) = x	! x array for interpolation 
*	     varr(1) = x	! value of x at which SF to be calculated here
*	     CALL strfn (ibeam,istruc,itm,iht,ndrv,varr,xc,FNarr(ix))
*	  ENDDO
C        ...Compute deuteron str.fn. by smearing nucleon SF (F2 only for now)
         !if (istruc.ne.2) then
         !   write(*,*) 'ERROR (deuteronsf): F1 and F3 not implemented'
         !   stop
         !end if
         CALL SMEARF2 (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,Fd)
      else
         print*, 'ERROR(DEUTERONSF): inuke out of range =',inuke
      end if
	
      RETURN
      END


C **********************************************************************
      SUBROUTINE SMEARF2 (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,FFd)
C
C  Smear nucleon structure function (array FNarr dimension nx defined
C  over x array xarr) with nucleon light-cone momentum distribution.
C
C  xarr is free-nucleon x in [0,1]
C  Convolution set up for "deuteron" xD in [0,1], rather than
C    "nucleon" x in [0,2]
C    => need factor 2 conversion xN = 2 xD [change later].
C
C  If istruc=-1, smears the PDFs at gamma=1, results are returned
C  in common/nPDF/
C
C  NOTE: az,an must have been previously set in common/target/az,an
C
C **********************************************************************
      IMPLICIT NONE

      INTEGER	ibeam,inuke,istruc,itm,iht,ndrv,dRN
      double precision v(4),xc(100),FFd,weight
      
      INTEGER i,ibj,istrg,ioff,ismear,iwfn,ishad,ilam
      REAL*8  xD,x,yD,yDmin,yDmax,fy_diag,fy_off,FFN,F2N,err,Q2,gamma
      REAL*8  PHI_INT2D,PHI_INT3D
      REAL*8  OFF_MST,OFF_KP,OFF_mKP,DEL_SHAD,shad,offsh
      REAL*8  OFF_mKP_fit
      double precision S0,S,U,D,UB,DB,SB,CB,BB,GLUE

*    *** common blocks

*     Gaussian integration
*      double precision xi(16),WI(16),xx(17)
*      double precision xi(32),WI(32),xx(33)
      double precision xi(96),WI(96),xx(97)
      integer nterms
*      COMMON/GAUS16/XI,WI,NTERMS,XX
*      COMMON/GAUS32/XI,WI,NTERMS,XX
      COMMON/GAUS96/XI,WI,NTERMS,XX

*     Nuclear PDFS
      double precision UA,DA,UBA,DBA,SBA,CBA,BBA,GLUEA
      common/nPDF/UA,DA,UBA,DBA,SBA,CBA,BBA,GLUEA

*     Q^2 min and max
      double precision Q02,Q2MAX
      COMMON/Q2STUF/ Q02,Q2MAX

C    *** Constants

C    ... hc [MeV.fm]           - conversion factor
C    ... mN [MeV] mNGeV [GeV]  - nucleon mass 
C    ... MD [MeV]              - deuterium mass, including binding energy
      REAL*8  pi,hc,mN,mNGeV,MD 
      parameter(pi=3.141592653589793d0,hc=197.327D0,mN=938.91897D0
     $     ,mNGeV=mN/1d3,mD=2*mN-2.224575D0 )

c...Smearing model 

      ! selects Bj limit, off-shell and smearing model
      call split_nuke(inuke,ilam,ishad,istrg,ioff,iwfn,ismear)  

C...Value of x,gamma at which convolution to be made
      x = v(1)	
      xD = x/2                  ! deuteron Bjorken variable
      Q2 = v(2)

      if (ilam.ge.1) then
*    ... Bjorken limit
         gamma = 1d0
      else
*    ... finite Q^2
         gamma = DSQRT(1.D0 + 4*x**2 * mNGeV**2/Q2)
      end if
	
C...Convolution approximation (x used in convolution here is xD in [0,1])
      yDmax = 1d0 
      yDmin = xD	
      FFD = 0d0
      if(istruc.eq.-1) then
         UA = 0
         DA = 0
         UBA = 0
         DBA = 0
         SBA = 0
         CBA = 0
         BBA = 0
         GLUEA  = 0
      end if

      DO I=1,NTERMS
         yD=0.5*(yDmax-yDmin)*XI(I)+0.5*(yDmax+yDmin) ! y = yD in [x,1]
         v(1)=xD/yD		! computes FN at xN=xD/yD
         if ((ioff.eq.0).or.(ioff.ge.2)) then 
            ! 2D interpolate smearing fn. in y and gamma 
            if (istruc.eq.0) then 
               !!!! FL = f0 * FLN + f1 * F2N
               fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
               CALL strfn (ibeam,0,itm,iht,ndrv,v,xc,FFN)
               fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) 
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            else if (istruc.eq.1) then 
               !!!! (xF1) = (f0 * (xN*F1N) + (1/4)*f1 * F2N
               fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
               CALL strfn (ibeam,1,itm,iht,ndrv,v,xc,FFN)
               FFN = v(1)*FFN
               fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) / 4d0
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            else if (istruc.eq.2) then 
               !!! F2 = f2 * F2N
               fy_diag = PHI_INT2D (yD,gamma,2,ismear,iwfn)
               CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,FFN)
               fy_off  = 0d0
               F2N = 0d0
            else if (istruc.eq.3) then 
               !!!! xF3 = f0 * xF3N
               print*, 'ERROR(SMEARF2): F3 convolutin formula yet to be'
     &              //'checked analytically'
               stop
               !fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn)
               !CALL strfn (ibeam,3,itm,iht,ndrv,v,xc,FFN)
               !fy_off  = 0d0
               !F2N = 0d0
            else if (istruc.eq.-1) then
               !!! Nuclear PDFs as convolution of nucleon PDFs
               !!! Beware: it makes physical sense only at very large Q^2
               fy_diag = PHI_INT2D (yD,1d0,0,ismear,iwfn)  ! gamma=1
               !!! the following 3 lines are for compatibility with fitpack...
               S0=DLOG(Q02/XC(1)**2)
               S=DLOG(DLOG(Q2/XC(1)**2)/S0)
               call fsupdf(0,x,s,u,d,ub,db,sb,cb,bb,glue)
               !!! if working with stand-alone package, could use this 1-line:
               !!! call PDFsa(x,Q2,U,D,UB,DB,SB,CB,BB,GLUE) 
            else
               fy_diag = 0d0
               FFN = 0d0
               fy_off = 0d0
               F2N = 0d0
            end if
         else if (ioff.eq.1) then            ! 3D interpolate smearing fn
            print*, 'WARNING!'
            print*, 'KP off-shell smearing has become OBSOLETE!!! '
            print*, '(Needs to be rewritten as 2D interpolation '
            print*, 'using the p2-m2 averaged smearing functions)'
            print*, 'STOPPING HERE!'
            stop
*               fy = PHI_INT3D (yD,gamma,xD/yD,ismear,ioff)  
*                                                ! in y, gamma, xN=xB/y=xD/yD
         else 
            write(*,*) 'ERROR (smearf2): ioff out of range = ',ioff
            stop
         end if

         if (v(1).lt.0.995d0) then
            weight = .5*(yDmax-yDmin)*WI(I)
            if (istruc.ne.-1) then ! Structure function convolution
               FFd = FFd + weight * (fy_diag*FFN + fy_off*F2N)
            else                   ! PDF convolution
               UA    = UA    + weight * fy_diag*U
               DA    = DA    + weight * fy_diag*D
               UBA   = UBA   + weight * fy_diag*UB
               DBA   = DBA   + weight * fy_diag*DB
               SBA   = SBA   + weight * fy_diag*SB
               CBA   = CBA   + weight * fy_diag*CB
               BBA   = BBA   + weight * fy_diag*BB
               GLUEA = GLUEA + weight * fy_diag*GLUE
            end if
         end if
         
         !print*, '* fy_diag _off =',fy_diag, fy_off
         !print*, '* FFN, F2N =', FFN,F2N
         !print*, '* FFd =',FFd
         !stop

      ENDDO
      if (istruc.eq.1) then       ! convolution done for x*F1_A, returns F1_A
         FFd = FFd/x
      else if (istruc.eq.-1) then ! nPDFs are output in common/nPDF/  
         FFd = -9999999d0
      end if


*    *** Parametrized off-shell corrections (F2 only!!)
      offsh=0d0
      if (ioff.eq.2) then
         ! Off-shell corrections by Melnitchouk-Schreiber-Thomas
         ! (parametrization by W.Melnitchouk) 
         offsh = OFF_MST(x)
      else if (ioff.eq.3) then
         ! approximate Kulagin-Petti NPA(2007) 
         ! (parametrization by W.Melnitchouk) 
         offsh = OFF_KP(x)
      else if (ioff.ge.4) then
         ! modified Kulagin-Petti (by W.Melnitchouck) 
         if (istrg.eq.0) then
            dRN = 10             ! max nucleon swelling
         else
            dRN = istrg
         end if
         offsh = OFF_mKP_fit(x,iwfn,dRN)
      end if 
      FFd = FFd/(1.D0-offsh)

*    *** Parametrized deuterium shadowing
      shad = 0d0
      if (ishad.eq.1) then
         ! Melnitchouk and Thomas, PRD 47, 3783 (1993)
         shad = DEL_SHAD(x,Q2)         
      end if
      FFd = FFd/(1.d0-shad)


       v(1) = x

       RETURN
       END


************************************************
*     divides nuclear correction flag into its components
      subroutine split_nuke(inuke,ilam,ishad,istrg,ioff,iwfn,ismear) 
*     Author: A.Accardi
*     date: 18 Feb 2010
*
*     INPUT:
*
*       inuke   = (i) nuclear correction: ABCDEF
C                     see 'deuteronsf' routine
*
*     OUTPUT:
*
*       ilam    = (i) A
*       ishad   = (i) B
*       istrng  = (i) C
*       ioff    = (i) D
*       iwfn    = (i) E
*       ismear  = (i) F
*

      implicit none
      
      integer inuke,ilam,ishad,istrg,ioff,iwfn,ismear

      ilam   = inuke/100000              ! treatment of TMC scale Lambda 
      ishad  = (inuke-ilam*100000)/10000 ! shadowing corrections
      istrg  = (inuke-ilam*100000        ! strength of off-shell corrections
     &               -ishad*10000)/1000     
      ioff   = (inuke-ilam *100000       ! off-shell corrections
     &               -ishad*10000
     &               -istrg*1000)/100     
      iwfn   = (inuke-ilam *100000       ! nuclear wave function
     &               -ishad*10000
     &               -istrg*1000
     &               -ioff*100)/10     
      ismear = (inuke-ilam *100000       ! convolution formula
     &               -ishad*10000        ! or deuteron correction model 
     &               -istrg*1000
     &               -ioff*100
     &               -iwfn*10)     

      return
      end



************************************************
*     Density model correction gomez(f) = N/D
      function gomez(f,x)
      implicit double precision (a-h,o-z)
      data p1,p2,p3,p4,p5,p6,p7/1.0164d0,-.047762d0,-.13354d0,
     2.35303d0,.22719d0,-1.2906d0,5.6075d0/
c
c  Fit to the data for F2D/F2N where N=(p+n)/2
c  as extractred by Gomez et al PRD 49, 4348 (1994)
c  [Original reference by Frankfurt and Strikman]
c
      THEORY=P1+P2*X+P3*x**2+P4*x**3+P5*x**4
     2+P6*x**5
      deufac=p7*(1./x-1.)
      if(deufac.ge.20.)deufac=20.
      deucor=1.-exp(-deufac)
      theory=theory/deucor
      gomez=f/theory
      RETURN
      END


C ***********************************************************************
      FUNCTION OFF_MST (x)
C
C  Nucleon off-shell correction in MST model.
C
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_MST = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Parameters defined for x = nucleon scaling variable in [0,2].
C
C  Ref: Melnitchouk, Schreiber, Thomas, PRD 49, 1183 (1994);
C       PLB 335, 11 (1994)
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_MST,x
      REAL*8  aoff(0:5)
      DATA    aoff /-0.014D0, 3.D0, 20.D0, 1.067D0, 1.5D0, 18.D0/

      OFF_MST = 0.D0
      IF (x.LE.0.16D0) RETURN

      OFF_MST = aoff(0) * (1.D0 + aoff(1)*x**aoff(2))
     &     * (1.D0 - (aoff(3) - x**aoff(4))**aoff(5))

      RETURN
      END



C ***********************************************************************
      FUNCTION OFF_KP (x)
C
C  Analytic parametrization of nucleon off-shell correction from
C    Kulagin-Petti fit
C
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_KO = del^off F2d / F2d
C                   = off-shell correction w.r.t. F2d
C
C  Parameters defined for x = nucleon scaling variable in [0,2]. 
C
C  Ref: Kulagin, Petti, NPA 765, 126 (2006).
C
C  July 2010.
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_KP,x
      REAL*8  aoff(0:4),boff(0:13),coff(0:5)  
      DATA    aoff /-0.010675D0, 0.27315D0, -0.96047D0, 1.2396D0,
     &     -0.71897D0/
      DATA    boff /-0.11185D-1, 0.29485D0, -1.3024D0, 4.1337D0,
     &     -17.047D0, 66.378D0, -184.54D0, 293.49D0,
     &     -64.581D0, -785.17D0, 1779.8D0, -1908.6D0,
     &     1066.9D0, -250.54D0/
      DATA    coff /3.9913D0, 4.3786D0, 2.4799D0, 2.5043D0,
     &     -3.9996D0, 0.018932D0/

       OFF_KP = 0.D0

C...13th order polynomial fit valid for x < 0.91
c       OFF_KP = boff(0) + boff(1)*x + boff(2)*x**2 + boff(3)*x**3
c     &        + boff(4)*x**4 + boff(5)*x**5 + boff(6)*x**6
c     &        + boff(7)*x**7 + boff(8)*x**8 + boff(9)*x**9
c     &        + boff(10)*x**10 + boff(11)*x**11 + boff(12)*x**12
c     &        + boff(13)*x**13
c
c       IF (x.GT.0.85D0) RETURN
C...4th order polynomial fit valid for x < 0.86
c       OFF_KP = aoff(0) + aoff(1)*x + aoff(2)*x**2 + aoff(3)*x**3
c     &        + aoff(4)*x**4

C...6-parameter fit  (courtesy of Simona Malace)
      OFF_KP = coff(0) + coff(1)*x + coff(2)*x**coff(3)
     &     + coff(4)*DEXP(x) + coff(5)/DLOG(x)

      RETURN
      END



C ***********************************************************************
        FUNCTION OFF_mKP (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]).
C  Fit by O. Hen (5/24/11 - 6/6/11).
C  Fortran code by W. Melnitchouk (5/27/11).
C
C  wfn: 1 (AV18)
C       2 (CD-Bonn)
C       3 (WJC-1)
C       4 (WJC-2)
C
C  dRN: 0 (0.0%)        [% change in nucleon radius in the deuteron]
C       1 (0.3%)
C       2 (0.6%)
C       3 (0.9%)
C       4 (1.2%)
C       5 (1.5%)
C       6 (1.8%)
C       7 (2.1%)
C       8 (2.4%)
C       9 (2.7%)
C       10 (3.0%)
C       
C ***********************************************************************
        IMPLICIT NONE
        REAL*8  OFF_mKP, x
        INTEGER wfn, dRN
        REAL*8  p(0:9), q(0:3)

        OFF_mKP = 0.D0
        IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
        IF (dRN.LT.0 .OR. dRN.GT.10) RETURN

! .......................................................................
        IF (wfn.EQ.1) THEN              ! AV18
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00196328D0
            p(1) = -0.0406043D0
            p(2) = 0.0854386D0
            p(3) = 0.598314D0
            p(4) = -5.18768D0
            p(5) = 19.3408D0
            p(6) = -40.8099D0
            p(7) = 50.0192D0
            p(8) = -33.2016D0
            p(9) = 9.25861D0
            q(0) = 1.70778D0
            q(1) = -3.8263D0
            q(2) = 2.04641D0
            q(3) = 0.132957D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00337272D0
            p(1) = -0.0589391D0
            p(2) = 0.204111D0
            p(3) = -0.609015D0
            p(4) = 1.90664D0
            p(5) = -5.05949D0
            p(6) = 9.39415D0
            p(7) = -10.967D0
            p(8) = 7.17457D0
            p(9) = -2.01069D0
            q(0) = -3.19753D0
            q(1) = 9.92432D0
            q(2) = -10.1881D0
            q(3) = 3.44211D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00479856D0
            p(1) = -0.0789437D0
            p(2) = 0.365415D0
            p(3) = -2.28553D0
            p(4) = 11.7184D0
            p(5) = -38.5546D0
            p(6) = 77.7781D0
            p(7) = -93.3614D0
            p(8) = 61.2492D0
            p(9) = -16.9552D0
            q(0) = -2.58892D0
            q(1) = 5.16641D0
            q(2) = -1.68202D0
            q(3) = -1.00944D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00624526D0
            p(1) = -0.101004D0
            p(2) = 0.578828D0
            p(3) = -4.53405D0
            p(4) = 24.8387D0
            p(5) = -83.1084D0
            p(6) = 168.235D0
            p(7) = -201.705D0
            p(8) = 131.896D0
            p(9) = -36.336D0
            q(0) = 8.43383D0
            q(1) = -34.09D0
            q(2) = 44.9658D0
            q(3) = -19.5377D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00771879D0
            p(1) = -0.125631D0
            p(2) = 0.85681D0
            p(3) = -7.48852D0
            p(4) = 42.0315D0
            p(5) = -141.24D0
            p(6) = 285.726D0
            p(7) = -341.749D0
            p(8) = 222.735D0
            p(9) = -61.1069D0
            q(0) = 38.4052D0
            q(1) = -135.585D0
            q(2) = 159.821D0
            q(3) = -63.0094D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.00922543D0
            p(1) = -0.153372D0
            p(2) = 1.21284D0
            p(3) = -11.2957D0
            p(4) = 64.1424D0
            p(5) = -215.764D0
            p(6) = 435.837D0
            p(7) = -520.011D0
            p(8) = 337.89D0
            p(9) = -92.3594D0
            q(0) = 102.447D0
            q(1) = -348.309D0
            q(2) = 395.805D0
            q(3) = -150.485D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0107738D0
            p(1) = -0.184965D0
            p(2) = 1.66493D0
            p(3) = -16.1495D0
            p(4) = 92.2796D0
            p(5) = -310.337D0
            p(6) = 625.768D0
            p(7) = -744.838D0
            p(8) = 482.61D0
            p(9) = -131.475D0
            q(0) = 228.085D0
            q(1) = -761.2D0
            q(2) = 848.715D0
            q(3) = -316.366D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0123752D0
            p(1) = -0.22138D0
            p(2) = 2.23652D0
            p(3) = -22.3008D0
            p(4) = 127.866D0
            p(5) = -429.628D0
            p(6) = 864.673D0
            p(7) = -1026.79D0
            p(8) = 663.511D0
            p(9) = -180.187D0
            q(0) = 467.305D0
            q(1) = -1541.83D0
            q(2) = 1698.6D0
            q(3) = -625.146D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.014041D0
            p(1) = -0.2636D0
            p(2) = 2.95198D0
            p(3) = -30.0151D0
            p(4) = 172.431D0
            p(5) = -578.727D0
            p(6) = 1162.64D0
            p(7) = -1377.62D0
            p(8) = 888.008D0
            p(9) = -240.454D0
            q(0) = 923.217D0
            q(1) = -3021.86D0
            q(2) = 3301.08D0
            q(3) = -1203.91D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0157866D0
            p(1) = -0.312963D0
            p(2) = 3.84398D0
            p(3) = -39.6446D0
            p(4) = 227.989D0
            p(5) = -764.276D0
            p(6) = 1532.74D0
            p(7) = -1812.46D0
            p(8) = 1165.61D0
            p(9) = -314.77D0
            q(0) = 1814.91D0
            q(1) = -5904.72D0
            q(2) = 6408.9D0
            q(3) = -2321.16D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0176352D0
            p(1) = -0.371422D0
            p(2) = 4.9594D0
            p(3) = -51.687D0
            p(4) = 297.343D0
            p(5) = -995.421D0
            p(6) = 1992.79D0
            p(7) = -2351.74D0
            p(8) = 1509.02D0
            p(9) = -406.45D0
            q(0) = 3653.87D0
            q(1) = -11829.9D0
            q(2) = 12773.5D0
            q(3) = -4600.46D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.2) THEN         ! CD-Bonn
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00162711D0
            p(1) = -0.0327193D0
            p(2) = 0.0595688D0
            p(3) = 0.492763D0
            p(4) = -3.90178D0
            p(5) = 14.0632D0
            p(6) = -29.1224D0
            p(7) = 35.2469D0
            p(8) = -23.1699D0
            p(9) = 6.40719D0
            q(0) = -1.42342D0
            q(1) = 5.78443D0
            q(2) = -7.67182D0
            q(3) = 3.35475D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00303191D0
            p(1) = -0.0496808D0
            p(2) = 0.154688D0
            p(3) = -0.497176D0
            p(4) = 1.91967D0
            p(5) = -5.86469D0
            p(6) = 11.635D0
            p(7) = -13.9418D0
            p(8) = 9.17025D0
            p(9) = -2.55081D0
            q(0) = -0.596239D0
            q(1) = 1.40214D0
            q(2) = -0.882721D0
            q(3) = 0.05642D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00444406D0
            p(1) = -0.0675313D0
            p(2) = 0.273155D0
            p(3) = -1.74564D0
            p(4) = 9.23849D0
            p(5) = -30.8037D0
            p(6) = 62.4092D0
            p(7) = -74.9276D0
            p(8) = 49.0612D0
            p(9) = -13.5357D0
            q(0) = 4.81023D0
            q(1) = -18.1571D0
            q(2) = 22.6966D0
            q(3) = -9.44371D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00586546D0
            p(1) = -0.0864343D0
            p(2) = 0.418993D0
            p(3) = -3.29639D0
            p(4) = 18.3067D0
            p(5) = -61.5928D0
            p(6) = 124.867D0
            p(7) = -149.658D0
            p(8) = 97.7382D0
            p(9) = -26.8755D0
            q(0) = 17.3269D0
            q(1) = -61.129D0
            q(2) = 72.0038D0
            q(3) = -28.3804D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00729865D0
            p(1) = -0.106607D0
            p(2) = 0.597462D0
            p(3) = -5.20583D0
            p(4) = 29.4458D0
            p(5) = -99.2915D0
            p(6) = 201.094D0
            p(7) = -240.554D0
            p(8) = 156.726D0
            p(9) = -42.9721D0
            q(0) = 40.7797D0
            q(1) = -139.933D0
            q(2) = 160.481D0
            q(3) = -61.6055D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.008746D0
            p(1) = -0.128256D0
            p(2) = 0.813724D0
            p(3) = -7.53053D0
            p(4) = 42.9835D0
            p(5) = -144.996D0
            p(6) = 293.277D0
            p(7) = -350.176D0
            p(8) = 227.651D0
            p(9) = -62.2584D0
            q(0) = 81.0211D0
            q(1) = -273.524D0
            q(2) = 308.599D0
            q(3) = -116.49D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0102108D0
            p(1) = -0.151667D0
            p(2) = 1.0747D0
            p(3) = -10.3452D0
            p(4) = 59.3469D0
            p(5) = -200.119D0
            p(6) = 404.201D0
            p(7) = -481.759D0
            p(8) = 312.554D0
            p(9) = -85.2724D0
            q(0) = 147.141D0
            q(1) = -491.286D0
            q(2) = 548.022D0
            q(3) = -204.411D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0116974D0
            p(1) = -0.177204D0
            p(2) = 1.38924D0
            p(3) = -13.744D0
            p(4) = 79.0702D0
            p(5) = -266.411D0
            p(6) = 537.302D0
            p(7) = -639.272D0
            p(8) = 413.922D0
            p(9) = -112.667D0
            q(0) = 253.515D0
            q(1) = -839.606D0
            q(2) = 928.653D0
            q(3) = -343.264D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.0132096D0
            p(1) = -0.205202D0
            p(2) = 1.76563D0
            p(3) = -17.8183D0
            p(4) = 102.683D0
            p(5) = -345.645D0
            p(6) = 696.108D0
            p(7) = -826.844D0
            p(8) = 534.376D0
            p(9) = -145.138D0
            q(0) = 423.458D0
            q(1) = -1393.59D0
            q(2) = 1531.15D0
            q(3) = -561.925D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0147536D0
            p(1) = -0.236188D0
            p(2) = 2.21657D0
            p(3) = -22.7029D0
            p(4) = 130.947D0
            p(5) = -440.312D0
            p(6) = 885.492D0
            p(7) = -1050.08D0
            p(8) = 677.422D0
            p(9) = -183.603D0
            q(0) = 695.863D0
            q(1) = -2278.36D0
            q(2) = 2489.67D0
            q(3) = -908.341D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0163345D0
            p(1) = -0.270612D0
            p(2) = 2.75336D0
            p(3) = -28.5225D0
            p(4) = 164.588D0
            p(5) = -552.843D0
            p(6) = 1110.3D0
            p(7) = -1314.67D0
            p(8) = 846.662D0
            p(9) = -229.017D0
            q(0) = 1138.1D0
            q(1) = -3710.33D0
            q(2) = 4035.93D0
            q(3) = -1465.21D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.3) THEN         ! WJC-1
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00274105D0
            p(1) = -0.060672D0
            p(2) = 0.212256D0
            p(3) = 0.0219547D0
            p(4) = -3.19661D0
            p(5) = 14.5887D0
            p(6) = -33.5328D0
            p(7) = 43.4176D0
            p(8) = -30.1071D0
            p(9) = 8.74272D0
            q(0) = 7.79726D0
            q(1) = -23.5139D0
            q(2) = 23.1821D0
            q(3) = -7.38804D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00411184D0
            p(1) = -0.0770866D0
            p(2) = 0.273399D0
            p(3) = -0.512489D0
            p(4) = -0.00242581D0
            p(5) = 3.08803D0
            p(6) = -8.72114D0
            p(7) = 11.7869D0
            p(8) = -8.10611D0
            p(9) = 2.26714D0
            q(0) = -2.53615D0
            q(1) = 8.56149D0
            q(2) = -9.65752D0
            q(3) = 3.63462D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00549862D0
            p(1) = -0.0951343D0
            p(2) = 0.376461D0
            p(3) = -1.51089D0
            p(4) = 5.8943D0
            p(5) = -17.5102D0
            p(6) = 34.3868D0
            p(7) = -41.5302D0
            p(8) = 27.87D0
            p(9) = -7.98734D0
            q(0) = -12.1657D0
            q(1) = 37.7354D0
            q(2) = -38.6311D0
            q(3) = 12.9756D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.006907D0
            p(1) = -0.1153D0
            p(2) = 0.533373D0
            p(3) = -3.10283D0
            p(4) = 15.2396D0
            p(5) = -49.6885D0
            p(6) = 100.722D0
            p(7) = -122.297D0
            p(8) = 81.4779D0
            p(9) = -22.9919D0
            q(0) = -18.911D0
            q(1) = 56.7604D0
            q(2) = -55.7012D0
            q(3) = 17.6606D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00834428D0
            p(1) = -0.138215D0
            p(2) = 0.759563D0
            p(3) = -5.45509D0
            p(4) = 28.9898D0
            p(5) = -96.6165D0
            p(6) = 196.556D0
            p(7) = -237.821D0
            p(8) = 157.337D0
            p(9) = -43.9692D0
            q(0) = -18.5907D0
            q(1) = 51.8677D0
            q(2) = -45.7426D0
            q(3) = 12.1465D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.0098194D0
            p(1) = -0.16465D0
            p(2) = 1.07396D0
            p(3) = -8.77272D0
            p(4) = 48.3221D0
            p(5) = -162.2D0
            p(6) = 329.623D0
            p(7) = -397.11D0
            p(8) = 261.143D0
            p(9) = -72.4259D0
            q(0) = -3.2391D0
            q(1) = -2.99953D0
            q(2) = 19.6733D0
            q(3) = -13.9123D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0113437D0
            p(1) = -0.195584D0
            p(2) = 1.50055D0
            p(3) = -13.3152D0
            p(4) = 74.7236D0
            p(5) = -251.364D0
            p(6) = 509.657D0
            p(7) = -611.486D0
            p(8) = 400.039D0
            p(9) = -110.248D0
            q(0) = 42.4961D0
            q(1) = -157.808D0
            q(2) = 194.778D0
            q(3) = -80.1444D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.012932D0
            p(1) = -0.232291D0
            p(2) = 2.07042D0
            p(3) = -19.4176D0
            p(4) = 110.11D0
            p(5) = -370.435D0
            p(6) = 749.127D0
            p(7) = -895.409D0
            p(8) = 583.123D0
            p(9) = -159.831D0
            q(0) = 148.975D0
            q(1) = -510.974D0
            q(2) = 585.939D0
            q(3) = -224.882D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.014603D0
            p(1) = -0.276387D0
            p(2) = 2.82299D0
            p(3) = -27.5043D0
            p(4) = 156.902D0
            p(5) = -527.402D0
            p(6) = 1063.77D0
            p(7) = -1267.1D0
            p(8) = 821.851D0
            p(9) = -224.187D0
            q(0) = 378.701D0
            q(1) = -1264.44D0
            q(2) = 1410.67D0
            q(3) = -526.229D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0163792D0
            p(1) = -0.329812D0
            p(2) = 3.80591D0
            p(3) = -38.0918D0
            p(4) = 218.065D0
            p(5) = -732.077D0
            p(6) = 1472.94D0
            p(7) = -1749.06D0
            p(8) = 1130.38D0
            p(9) = -307.041D0
            q(0) = 868.106D0
            q(1) = -2857.75D0
            q(2) = 3141.02D0
            q(3) = -1153.18D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0182922D0
            p(1) = -0.395284D0
            p(2) = 5.08523D0
            p(3) = -51.8887D0
            p(4) = 297.626D0
            p(5) = -997.705D0
            p(6) = 2002.64D0
            p(7) = -2371.26D0
            p(8) = 1527.47D0
            p(9) = -413.31D0
            q(0) = 1941.33D0
            q(1) = -6332.67D0
            q(2) = 6892.95D0
            q(3) = -2504.2D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.4) THEN         ! WJC-2
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00214184D0
            p(1) = -0.0450598D0
            p(2) = 0.110368D0
            p(3) = 0.507098D0
            p(4) = -4.94836D0
            p(5) = 18.9418D0
            p(6) = -40.4816D0
            p(7) = 50.0558D0
            p(8) = -33.4725D0
            p(9) = 9.40103D0
            q(0) = 3.4267D0
            q(1) = -9.35273D0
            q(2) = 7.94702D0
            q(3) = -1.95589D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00354093D0
            p(1) = -0.0629637D0
            p(2) = 0.216152D0
            p(3) = -0.549377D0
            p(4) = 1.27716D0
            p(5) = -2.60462D0
            p(6) = 4.14477D0
            p(7) = -4.53293D0
            p(8) = 2.93507D0
            p(9) = -0.844147D0
            q(0) = -3.57622D0
            q(1) = 11.3618D0
            q(2) = -11.9905D0
            q(3) = 4.19102D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00495676D0
            p(1) = -0.0825665D0
            p(2) = 0.365339D0
            p(3) = -2.08408D0
            p(4) = 10.2762D0
            p(5) = -33.4462D0
            p(6) = 67.3774D0
            p(7) = -81.0641D0
            p(8) = 53.405D0
            p(9) = -14.8687D0
            q(0) = -6.54542D0
            q(1) = 18.3482D0
            q(2) = -16.3313D0
            q(3) = 4.42111D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00639435D0
            p(1) = -0.104301D0
            p(2) = 0.568512D0
            p(3) = -4.21151D0
            p(4) = 22.7055D0
            p(5) = -75.76D0
            p(6) = 153.525D0
            p(7) = -184.556D0
            p(8) = 121.11D0
            p(9) = -33.5121D0
            q(0) = -1.32741D0
            q(1) = -1.99123D0
            q(2) = 9.77254D0
            q(3) = -6.67373D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00785974D0
            p(1) = -0.128688D0
            p(2) = 0.838473D0
            p(3) = -7.07025D0
            p(4) = 39.3602D0
            p(5) = -132.183D0
            p(6) = 267.806D0
            p(7) = -321.086D0
            p(8) = 209.89D0
            p(9) = -57.7897D0
            q(0) = 19.4618D0
            q(1) = -73.7158D0
            q(2) = 92.4639D0
            q(3) = -38.5673D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.00936044D0
            p(1) = -0.156377D0
            p(2) = 1.19109D0
            p(3) = -10.8319D0
            p(4) = 61.2231D0
            p(5) = -205.969D0
            p(6) = 416.647D0
            p(7) = -498.125D0
            p(8) = 324.456D0
            p(9) = -88.9448D0
            q(0) = 69.1367D0
            q(1) = -240.043D0
            q(2) = 278.518D0
            q(3) = -108.141D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0109059D0
            p(1) = -0.188179D0
            p(2) = 1.64619D0
            p(3) = -15.71D0
            p(4) = 89.5156D0
            p(5) = -301.152D0
            p(6) = 608.003D0
            p(7) = -724.9D0
            p(8) = 470.617D0
            p(9) = -128.508D0
            q(0) = 172.305D0
            q(1) = -580.597D0
            q(2) = 653.83D0
            q(3) = -246.288D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0125082D0
            p(1) = -0.225133D0
            p(2) = 2.22902D0
            p(3) = -21.976D0
            p(4) = 125.784D0
            p(5) = -422.823D0
            p(6) = 851.883D0
            p(7) = -1012.99D0
            p(8) = 655.643D0
            p(9) = -178.388D0
            q(0) = 376.045D0
            q(1) = -1247.3D0
            q(2) = 1381.83D0
            q(3) = -511.62D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.0141801D0
            p(1) = -0.26835D0
            p(2) = 2.96709D0
            p(3) = -29.9298D0
            p(4) = 171.756D0
            p(5) = -576.732D0
            p(6) = 1159.68D0
            p(7) = -1375.67D0
            p(8) = 887.914D0
            p(9) = -240.797D0
            q(0) = 774.973D0
            q(1) = -2544.77D0
            q(2) = 2789.44D0
            q(3) = -1021.08D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0159409D0
            p(1) = -0.319476D0
            p(2) = 3.90024D0
            p(3) = -39.9977D0
            p(4) = 229.853D0
            p(5) = -770.824D0
            p(6) = 1546.97D0
            p(7) = -1830.9D0
            p(8) = 1178.66D0
            p(9) = -318.677D0
            q(0) = 1572.33D0
            q(1) = -5126.02D0
            q(2) = 5575.96D0
            q(3) = -2024.29D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0178144D0
            p(1) = -0.380548D0
            p(2) = 5.07793D0
            p(3) = -52.7094D0
            p(4) = 303.083D0
            p(5) = -1014.97D0
            p(6) = 2033.09D0
            p(7) = -2400.96D0
            p(8) = 1541.84D0
            p(9) = -415.674D0
            q(0) = 3248.08D0
            q(1) = -10530.3D0
            q(2) = 11386.6D0
            q(3) = -4107.33D0
          ENDIF
        ENDIF

        IF (x.LE.0.9D0) THEN
          OFF_mKP = p(0) + p(1)*x + p(2)*x**2 + p(3)*x**3 + p(4)*x**4
     &            + p(5)*x**5 + p(6)*x**6 + p(7)*x**7 + p(8)*x**8
     &            + p(9)*x**9
        ELSE IF (x.GT.0.9D0) THEN
          OFF_mKP = q(0) + q(1)*x + q(2)*x**2 + q(3)*x**3
        ENDIF

        RETURN
        END




C ***********************************************************************
        FUNCTION DEL_SHAD (x,Q2)
C
C  Nuclear shadowing correction to F2d, including VMD, Pomeron and
C    meson (antishadowing) exchange.
C
C  Defined such that F2d = F2d(noshad) + \delta(shad) F2d, with the
C    output DEL_SHAD = \delta(shad) F2d / F2d.
C
C  Parameters defined for x = nucleon(!) scaling variable in [0,2].
C          
C  Ref: Melnitchouk and Thomas, PRD 47, 3783 (1993)
C
C  Relative correction, based on Eric Christy's fit, implemented
C    (absolute correction removed)
C    [WM: July 5, 2012]
C
C ***********************************************************************
        IMPLICIT NONE
        REAL*8  DEL_SHAD,x,Q2
        REAL*8  DEL_V,DEL_P,DEL_M, p1,p2,p3,p4
        REAL*8  A(0:4),B(0:4),C(0:4)
        DATA    A /-0.038D0, -0.04D0, 1.8D0, -3.4D0, 0.9D0/
        DATA    B /-0.003D0, -0.13D0, 5.D0,  -2.2D0, 0.4D0/
        DATA    C / 0.002D0,  0.03D0, 6.D0,   0.D0,  0.D0 /

        DEL_SHAD = 0.D0
        IF (x.GE.0.3D0) RETURN          ! param. valid for x <~ 0.3
   
C...VMD (Q2 in GeV^2)
	p1 = -0.084D0*Q2**(-0.71D0)
	p2 = 2.18D0 + 3.11D0*DEXP(-0.689D0*Q2)
	DEL_V = p1 * x**0.032D0 * (0.3D0-x)**1.375D0 * DEXP(-p2*x)

C...Pomeron-exchange
        p1 = -0.0213D0
        p2 = -0.0883D0
        p3 = 1.45D0 
        p4 = 18.63D0
	DEL_P = p1 * x**p2 * (0.3D0-x)**p3 * DEXP(-p4*x)
      
C...Meson-exchange (antishadowing)
	p1 = 0.00308D0
        p2 = 0.103D0
        p3 = -0.108D0
        p4 = 7.15D0
	DEL_M = p1 * x**p2 * (0.3D0-x)**p3 * DEXP(-p4*x)

        DEL_SHAD = DEL_V + DEL_P + DEL_M
      
        RETURN
        END



**********************************************

      subroutine smearfile(outfile,outlen,ismear,iwfn,ioff)

      implicit none

      character outfile*150

      integer len,oldlen,outlen,ismear,iwfn,ioff


      len = 4
      outfile = 'phi.'
      oldlen = len+1

c$$$      if (ismear.eq.2) then
c$$$         len = oldlen + 1
c$$$         outfile(oldlen:len) = 'km'
c$$$         oldlen = len + 1
c$$$      else 
      if (ismear.eq.3) then
         len = oldlen + 2
         outfile(oldlen:len) = 'wba'
         oldlen = len + 1
      else if (ismear.eq.4) then
         len = oldlen + 5
         outfile(oldlen:len) = 'wbarel'
         oldlen = len + 1
      else if (ismear.eq.5) then
         len = oldlen + 5
         outfile(oldlen:len) = 'aqv'
         oldlen = len + 1
      else if (ismear.eq.6) then
         len = oldlen + 5
         outfile(oldlen:len) = 'aqvc'
         oldlen = len + 1
      else if ((ismear.ne.0).and.(ismear.ne.1)) then
         len = oldlen + 7
         outfile(oldlen:len) = '_smear??'
         oldlen = len + 1
      end if

      if (iwfn.eq.0) then
         len = oldlen + 5
         outfile(oldlen:len) = '_PARIS'
         oldlen = len + 1
      else if (iwfn.eq.1) then
         len = oldlen + 4
         outfile(oldlen:len) = '_AV18'
         oldlen = len + 1
      else if (iwfn.eq.2) then
         len = oldlen + 6
         outfile(oldlen:len) = '_CDBONN'
         oldlen = len + 1
      else if (iwfn.eq.3) then
         len = oldlen + 4
         outfile(oldlen:len) = '_WJC1'
         oldlen = len + 1
      else if (iwfn.eq.4) then
         len = oldlen + 4
         outfile(oldlen:len) = '_WJC2'
         oldlen = len + 1
      else 
         len = oldlen + 5
         outfile(oldlen:len) = '_wfn??'
         oldlen = len + 1
      end if

      if (ioff.eq.1) then
         len = oldlen + 4
         outfile(oldlen:len) = '_KPoff'
         oldlen = len + 1
      else if (ioff.ne.0) then
         len = oldlen + 7
         outfile(oldlen:len) = '_offsh??'
         oldlen = len + 1
      end if

      outlen = len

      return
      end



C ***********************************************************************
        FUNCTION off_mKP_fit (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]), constrained to vanish at x=0
C
C  Defined such that F2d = F2d(conv) + del^off F2d 
C       with off_mKP_fit = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Fit by Eric Christy (6/14/12).
C  Fortran code by W. Melnitchouk (6/17/12).
C
C  wfn: 1 (AV18)
C	2 (CD-Bonn)
C	3 (WJC-1)
C	4 (WJC-2)
C
C  dRN: 0 (0.0%)	[% change in nucleon radius in the deuteron]
C	1 (0.3%)
C	2 (0.6%)
C	3 (0.9%)
C	4 (1.2%)
C	5 (1.5%)
C	6 (1.8%)
C	7 (2.1%)
C	8 (2.4%)
C	9 (2.7%)
C      10 (3.0%)
C
C ***********************************************************************
	IMPLICIT NONE
	REAL*8	off_mKP_fit, x
	INTEGER	wfn, dRN
	REAL*8  p(0:9)

	off_mKP_fit = 0.D0
	IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
	IF (dRN.LT.0 .OR. dRN.GT.10) RETURN
! .......................................................................
	IF (wfn.EQ.1) THEN		! AV18
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02628D0
	    p(1) = 1.20029D0
	    p(2) = 7.49503D0
	    p(3) = 2.01901D0
	    p(4) = 0.00789D0
	    p(5) = 0.46739D0
	    p(6) = 0.73242D0
	    p(7) = 0.00328D0
	    p(8) = 0.87228D0
	    p(9) = 0.06400D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = 0.03638D0
	    p(1) = 0.38307D0
	    p(2) = 8.01156D0
	    p(3) = 2.30992D0
	    p(4) = 0.09027D0
	    p(5) = 0.69521D0
	    p(6) = 0.75973D0
	    p(7) = -0.05098D0
	    p(8) = 1.18963D0
	    p(9) = -0.19192D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02260D0
	    p(1) = 1.45377D0
	    p(2) = 0.50628D0
	    p(3) = 13.92200D0
	    p(4) = 0.03558D0
	    p(5) = 0.75147D0
	    p(6) = 0.86335D0
	    p(7) = -0.01383D0
	    p(8) = 1.04749D0
	    p(9) = 0.42099D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06410D0
	    p(1) = 1.18883D0
	    p(2) = 6.96799D0
	    p(3) = 8.87113D0
	    p(4) = 0.02603D0
	    p(5) = 0.70504D0
	    p(6) = 1.44139D0
	    p(7) = 0.00004D0
	    p(8) = -1.14305D0
	    p(9) = 0.73785D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06237D0
	    p(1) = 2.03192D0
	    p(2) = 4.01755D0
	    p(3) = 6.83741D0
	    p(4) = 0.04701D0
	    p(5) = -0.00457D0
	    p(6) = 1.30967D0
	    p(7) = -0.00996D0
	    p(8) = -0.42418D0
	    p(9) = 0.27524D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.06759D0
	    p(1) = 1.95103D0
	    p(2) = 3.54215D0
	    p(3) = 11.77533D0
	    p(4) = 0.09269D0
	    p(5) = 0.56534D0
	    p(6) = 0.98398D0
	    p(7) = -0.03031D0
	    p(8) = 3.26913D0
	    p(9) = -0.45923D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.07007D0
	    p(1) = 2.30938D0
	    p(2) = 4.94226D0
	    p(3) = 8.95701D0
	    p(4) = 0.06933D0
	    p(5) = 0.07145D0
	    p(6) = 1.94887D0
	    p(7) = -0.01210D0
	    p(8) = 5.92311D0
	    p(9) = 0.14312D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11965D0
	    p(1) = 2.06149D0
	    p(2) = 5.38881D0
	    p(3) = 12.08265D0
	    p(4) = 0.19668D0
	    p(5) = 0.61820D0
	    p(6) = 0.80489D0
	    p(7) = -0.08735D0
	    p(8) = 3.74802D0
	    p(9) = -0.70773D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.14735D0
	    p(1) = 2.27109D0
	    p(2) = 8.23092D0
	    p(3) = 7.31581D0
	    p(4) = 0.11953D0
	    p(5) = 0.67459D0
	    p(6) = 1.59118D0
	    p(7) = -0.02700D0
	    p(8) = 4.52840D0
	    p(9) = -1.77765D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.27194D0
	    p(1) = 2.01340D0
	    p(2) = 10.71380D0
	    p(3) = 8.84886D0
	    p(4) = 0.09345D0
	    p(5) = 0.49802D0
	    p(6) = 1.28523D0
	    p(7) = -0.00474D0
	    p(8) = 0.58703D0
	    p(9) = 0.88354D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.69848D0
	    p(1) = 1.48173D0
	    p(2) = 17.44991D0
	    p(3) = 12.73730D0
	    p(4) = 0.13118D0
	    p(5) = 0.34598D0
	    p(6) = 1.65884D0
	    p(7) = -0.02215D0
	    p(8) = 1.21306D0
	    p(9) = 0.96399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.2) THEN		! CD-Bonn
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02820D0
	    p(1) = 0.85879D0
	    p(2) = 9.48856D0
	    p(3) = 2.18885D0
	    p(4) = 0.00070D0
	    p(5) = -5.61817D0
	    p(6) = 14.80512D0
	    p(7) = 0.00348D0
	    p(8) = -1.30292D0
	    p(9) = -0.73075D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.02996D0
	    p(1) = 0.35717D0
	    p(2) = 6.53843D0
	    p(3) = 3.88389D0
	    p(4) = 0.00758D0
	    p(5) = -16.50399D0
	    p(6) = 77.60083D0
	    p(7) = 0.00320D0
	    p(8) = 0.42334D0
	    p(9) = 0.28545D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03261D0
	    p(1) = 0.91185D0
	    p(2) = 8.49348D0
	    p(3) = 10.19681D0
	    p(4) = 0.01598D0
	    p(5) = 0.83748D0
	    p(6) = 1.55960D0
	    p(7) = 0.00085D0
	    p(8) = -0.63447D0
	    p(9) = 0.65632D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.03034D0
	    p(1) = 1.58677D0
	    p(2) = 3.21753D0
	    p(3) = 11.66572D0
	    p(4) = 0.04999D0
	    p(5) = 0.56688D0
	    p(6) = 0.94941D0
	    p(7) = -0.01453D0
	    p(8) = -0.89157D0
	    p(9) = 0.27160D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.04831D0
	    p(1) = 1.75241D0
	    p(2) = 4.74662D0
	    p(3) = 8.29052D0
	    p(4) = 0.04730D0
	    p(5) = 0.33550D0
	    p(6) = 1.18790D0
	    p(7) = -0.00678D0
	    p(8) = -0.42800D0
	    p(9) = 0.36573D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.09019D0
	    p(1) = 1.22091D0
	    p(2) = 1.30114D0
	    p(3) = 17.58701D0
	    p(4) = 0.08312D0
	    p(5) = 0.66902D0
	    p(6) = 0.60767D0
	    p(7) = -0.02035D0
	    p(8) = 0.95978D0
	    p(9) = 1.11322D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.17060D0
	    p(1) = 1.42115D0
	    p(2) = 7.24672D0
	    p(3) = 5.80680D0
	    p(4) = 0.09200D0
	    p(5) = 0.43367D0
	    p(6) = 1.56378D0
	    p(7) = -0.02338D0
	    p(8) = 0.44968D0
	    p(9) = 0.29678D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11026D0
	    p(1) = 1.85213D0
	    p(2) = 6.74413D0
	    p(3) = 7.74362D0
	    p(4) = 0.08467D0
	    p(5) = 0.24708D0
	    p(6) = 1.12274D0
	    p(7) = -0.01505D0
	    p(8) = 0.44209D0
	    p(9) = 0.36126D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.15291D0
	    p(1) = 1.83333D0
	    p(2) = 7.76495D0
	    p(3) = 7.04783D0
	    p(4) = 0.09206D0
	    p(5) = 0.08655D0
	    p(6) = 1.27460D0
	    p(7) = -0.01659D0
	    p(8) = 0.45536D0
	    p(9) = 0.29407D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.24143D0
	    p(1) = 1.50401D0
	    p(2) = 9.33393D0
	    p(3) = 11.62779D0
	    p(4) = 0.09454D0
	    p(5) = 0.36361D0
	    p(6) = 0.82058D0
	    p(7) = -0.00802D0
	    p(8) = 0.34851D0
	    p(9) = 0.50844D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.22196D0
	    p(1) = 1.87228D0
	    p(2) = 10.18898D0
	    p(3) = 9.21038D0
	    p(4) = 0.11850D0
	    p(5) = 0.34360D0
	    p(6) = 1.28278D0
	    p(7) = -0.01754D0
	    p(8) = 0.54540D0
	    p(9) = 0.53457D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.3) THEN		! WJC-1
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.02322D0
	    p(1) = 0.11213D0
	    p(2) = 3.71079D0
	    p(3) = 5.51496D0
	    p(4) = 0.00877D0
	    p(5) = 0.84639D0
	    p(6) = 0.66227D0
	    p(7) = -0.00621D0
	    p(8) = -0.39896D0
	    p(9) = 0.32012D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.00058D0
	    p(1) = 2.33827D0
	    p(2) = 2.35664D0
	    p(3) = 36.75823D0
	    p(4) = -0.00752D0
	    p(5) = 0.05286D0
	    p(6) = 1.27262D0
	    p(7) = 0.01269D0
	    p(8) = 1.72720D0
	    p(9) = 0.20652D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03373D0
	    p(1) = 0.93858D0
	    p(2) = 0.15704D0
	    p(3) = 10.71630D0
	    p(4) = -0.00235D0
	    p(5) = -0.11937D0
	    p(6) = 0.74925D0
	    p(7) = 0.00452D0
	    p(8) = 2.96830D0
	    p(9) = -2.89070D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.08982D0
	    p(1) = 0.73060D0
	    p(2) = -0.16543D0
	    p(3) = 12.37035D0
	    p(4) = 0.04407D0
	    p(5) = 0.47361D0
	    p(6) = 0.74570D0
	    p(7) = -0.00933D0
	    p(8) = 0.53186D0
	    p(9) = 0.26943D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.11990D0
	    p(1) = 1.19824D0
	    p(2) = 3.06386D0
	    p(3) = 8.55017D0
	    p(4) = 0.05815D0
	    p(5) = 0.06123D0
	    p(6) = 1.45024D0
	    p(7) = -0.01414D0
	    p(8) = 0.48172D0
	    p(9) = 0.25171D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.15292D0
	    p(1) = 1.01991D0
	    p(2) = 1.20661D0
	    p(3) = 13.31860D0
	    p(4) = 0.02571D0
	    p(5) = -1.56438D0
	    p(6) = 2.69042D0
	    p(7) = -0.00000D0
	    p(8) = 0.29759D0
	    p(9) = 0.97967D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.35935D0
	    p(1) = 0.44637D0
	    p(2) = -0.25510D0
	    p(3) = 16.70057D0
	    p(4) = 0.10634D0
	    p(5) = 0.61659D0
	    p(6) = 0.58524D0
	    p(7) = -0.03335D0
	    p(8) = 0.93904D0
	    p(9) = 0.89819D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.97384D0
	    p(1) = -0.24934D0
	    p(2) = -0.61349D0
	    p(3) = 18.43254D0
	    p(4) = 0.18772D0
	    p(5) = 0.49599D0
	    p(6) = 0.61366D0
	    p(7) = -0.08116D0
	    p(8) = 0.87175D0
	    p(9) = 0.24026D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.21641D0
	    p(1) = 1.74710D0
	    p(2) = 5.19387D0
	    p(3) = 10.61285D0
	    p(4) = 0.06655D0
	    p(5) = 0.01300D0
	    p(6) = 0.94503D0
	    p(7) = -0.00642D0
	    p(8) = 0.48859D0
	    p(9) = 0.16331D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.32283D0
	    p(1) = 1.71708D0
	    p(2) = 7.51556D0
	    p(3) = 9.68202D0
	    p(4) = 0.09871D0
	    p(5) = 0.18788D0
	    p(6) = 0.80490D0
	    p(7) = -0.01673D0
	    p(8) = 0.48879D0
	    p(9) = 0.21016D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.32064D0
	    p(1) = 2.07514D0
	    p(2) = 9.34847D0
	    p(3) = 8.17225D0
	    p(4) = 0.10772D0
	    p(5) = 0.50272D0
	    p(6) = 1.30663D0
	    p(7) = -0.01215D0
	    p(8) = 0.59432D0
	    p(9) = 0.65399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.4) THEN		! WJC-2
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.03490D0
	    p(1) = 0.78902D0
	    p(2) = -0.25256D0
	    p(3) = 7.98679D0
	    p(4) = 0.00913D0
	    p(5) = 0.74835D0
	    p(6) = 0.60145D0
	    p(7) = -0.00464D0
	    p(8) = 0.41358D0
	    p(9) = 0.22524D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.01119D0
	    p(1) = 0.50514D0
	    p(2) = 19.35710D0
	    p(3) = 3.32395D0
	    p(4) = 0.00670D0
	    p(5) = 1.38279D0
	    p(6) = 1.24216D0
	    p(7) = 0.00049D0
	    p(8) = 0.38623D0
	    p(9) = 0.23497D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02653D0
	    p(1) = 1.27315D0
	    p(2) = -0.53410D0
	    p(3) = 14.08029D0
	    p(4) = 0.01474D0
	    p(5) = 1.82129D0
	    p(6) = 1.99455D0
	    p(7) = -0.00090D0
	    p(8) = 3.96583D0
	    p(9) = 4.61316D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06301D0
	    p(1) = 1.10373D0
	    p(2) = -0.26356D0
	    p(3) = 15.04038D0
	    p(4) = 0.02428D0
	    p(5) = -0.15349D0
	    p(6) = 3.03168D0
	    p(7) = 0.00127D0
	    p(8) = -0.73818D0
	    p(9) = 0.07474D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06150D0
	    p(1) = 2.15792D0
	    p(2) = 2.18241D0
	    p(3) = 9.84713D0
	    p(4) = 0.03608D0
	    p(5) = -0.13604D0
	    p(6) = 1.12241D0
	    p(7) = -0.00695D0
	    p(8) = -0.35646D0
	    p(9) = 0.31793D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.07179D0
	    p(1) = 1.97917D0
	    p(2) = 3.47662D0
	    p(3) = 10.00224D0
	    p(4) = 0.04587D0
	    p(5) = 0.06416D0
	    p(6) = 1.10677D0
	    p(7) = -0.00391D0
	    p(8) = -0.42677D0
	    p(9) = 0.26619D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.09883D0
	    p(1) = 1.96788D0
	    p(2) = 5.19182D0
	    p(3) = 8.82173D0
	    p(4) = 0.06468D0
	    p(5) = 0.11297D0
	    p(6) = 1.63850D0
	    p(7) = -0.00872D0
	    p(8) = 0.52753D0
	    p(9) = 0.41794D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.14258D0
	    p(1) = 2.00822D0
	    p(2) = 6.23508D0
	    p(3) = 7.81846D0
	    p(4) = 0.07064D0
	    p(5) = -0.05869D0
	    p(6) = 1.24848D0
	    p(7) = -0.01160D0
	    p(8) = 0.48932D0
	    p(9) = 0.22001D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.16184D0
	    p(1) = 2.16963D0
	    p(2) = 7.62378D0
	    p(3) = 7.33369D0
	    p(4) = 0.09197D0
	    p(5) = 0.15692D0
	    p(6) = 1.80734D0
	    p(7) = -0.01561D0
	    p(8) = 0.53224D0
	    p(9) = 0.39357D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.20205D0
	    p(1) = 2.28733D0
	    p(2) = 9.10375D0
	    p(3) = 7.24877D0
	    p(4) = 0.08325D0
	    p(5) = 0.36941D0
	    p(6) = 2.39131D0
	    p(7) = -0.00057D0
	    p(8) = 0.41640D0
	    p(9) = 0.90531D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.95664D0
	    p(1) = 1.11409D0
	    p(2) = 19.00631D0
	    p(3) = 15.97282D0
	    p(4) = 0.15616D0
	    p(5) = 0.40229D0
	    p(6) = 0.85878D0
	    p(7) = -0.03123D0
	    p(8) = 6.75437D0
	    p(9) = -3.83159D0
	  ENDIF

	ENDIF

	off_mKP_fit = -( p(0) * x**p(3) * DEXP(p(1) * x**p(2))
     &                + p(4) * x*DEXP(((x-p(5))/p(6))**2)
     &                + x**0.5D0 * p(7) * DEXP(((x-p(9))/p(8))**2) )

	RETURN
	END




