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
C                               4 =  mKP    - modified Kulagin-Petti model,
C                                    parametrization at Deuterium level
C                               5 =  fmKP    - parton level mKP modeling
C                                    (work by L.Brady, summer 2013)
C                            E: Deuteron wave-function
C                               0=Paris
C                               1=AV18
C                               2=CDBonn
C                               3=WJC1
C                               4=WJC2
C                               5=AV18_KP
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
      double precision foff_diag,foff_off,FFNoff,F2Noff,FFdoff !fmKP variables
      double precision S0,S,U,D,UB,DB,SB,CB,BB,GLUE
      double precision uoff,doff,uboff,dboff,sboff,cboff,bboff,glueoff !fmKP variables

*    *** common blocks

*     Gaussian integration
*      double precision xi(16),WI(16),xx(17)
C       double precision xi(32),WI(32),xx(33)
C      double precision xi(64),WI(64),xx(65)
       double precision xi(96),WI(96),xx(97)
      integer nterms
*      COMMON/GAUS16/XI,WI,NTERMS,XX
C       　COMMON/GAUS32/XI,WI,NTERMS,XX
C      COMMON/GAUS64/XI,WI,NTERMS,XX
      COMMON/GAUS96/XI,WI,NTERMS,XX

*     Nuclear PDFS
      double precision UAo,DAo,UBAo,DBAo,SBAo,CBAo,BBAo,GLUEAo !fmKP variables
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
C..      xD = x/2                  ! deuteron Bjorken variable 
      xD=x*(mN/mD)        !Changed by Ishara   (May, 2020)
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
      FFDoff=0d0
      if(istruc.eq.-1) then
         UA = 0
         DA = 0
         UBA = 0
         DBA = 0
         SBA = 0
         CBA = 0
         BBA = 0
         GLUEA  = 0
         if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
            UAo = 0
            DAo = 0
            UBAo = 0
            DBAo = 0
            SBAo = 0
            CBAo = 0
            BBAo = 0
            GLUEAo  = 0
         endif
      end if

      DO I=1,NTERMS
         yD=0.5*(yDmax-yDmin)*XI(I)+0.5*(yDmax+yDmin) ! y = yD in [x,1]
         v(1)=xD/yD		! computes FN at xN=xD/yD
         ! 2D interpolate smearing fn. in y and gamma
         if (istruc.eq.0) then 
            !!!! FL = f0 * FLN + f1 * F2N
            fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
            CALL strfn (ibeam,0,itm,iht,ndrv,v,xc,FFN)
            fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) 
            CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = PHI_INT2D (yD,gamma,10,ismear,iwfn) 
               foff_off = PHI_INT2D (yD,gamma,11,ismear,iwfn) 
               call set_offshell_on (.true.)
               CALL strfn(ibeam,0,itm,iht,ndrv,v,xc,FFNoff)
               CALL strfn(ibeam,2,itm,iht,ndrv,v,xc,F2Noff)
               call set_offshell_on (.false.)
            endif
         else if (istruc.eq.1) then 
            !!!! (xF1) = (f0 * (xN*F1N) + (1/4)*f1 * F2N
            fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn) 
            CALL strfn (ibeam,1,itm,iht,ndrv,v,xc,FFN)
            FFN = v(1)*FFN
            fy_off = PHI_INT2D (yD,gamma,1,ismear,iwfn) / 4d0
            CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,F2N)
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = PHI_INT2D (yD,gamma,10,ismear,iwfn) 
               foff_off = PHI_INT2D (yD,gamma,11,ismear,iwfn) 
               call set_offshell_on (.true.)
               CALL strfn(ibeam,1,itm,iht,ndrv,v,xc,FFNoff)
               FFNoff = v(1)*FFNoff
               CALL strfn(ibeam,2,itm,iht,ndrv,v,xc,F2Noff)
               call set_offshell_on (.false.)
            endif
         else if (istruc.eq.2) then 
            !!! F2 = f2 * F2N
            fy_diag = PHI_INT2D (yD,gamma,2,ismear,iwfn)
            CALL strfn (ibeam,2,itm,iht,ndrv,v,xc,FFN)
            fy_off  = 0d0
            F2N = 0d0
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = PHI_INT2D (yD,gamma,12,ismear,iwfn) 
               foff_off = 0d0
               call set_offshell_on (.true.)
               CALL strfn(ibeam,2,itm,iht,ndrv,v,xc,FFNoff)
               F2Noff = 0d0
               call set_offshell_on (.false.)
            endif
         else if (istruc.eq.3) then 
            fy_diag = PHI_INT2D (yD,gamma,0,ismear,iwfn)
            CALL strfn (ibeam,3,itm,iht,ndrv,v,xc,FFN)
            fy_off  = 0d0
            F2N = 0d0
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = PHI_INT2D (yD,gamma,13,ismear,iwfn) 
               foff_off = 0d0
               call set_offshell_on (.true.)
               CALL strfn(ibeam,3,itm,iht,ndrv,v,xc,FFNoff)
               F2Noff = 0d0
               call set_offshell_on (.false.)
            endif
         else if (istruc.eq.-1) then
            !!! Nuclear PDFs as convolution of nucleon PDFs
            !!! Beware: it makes physical sense only at very large Q^2
            fy_diag = PHI_INT2D (yD,1d0,0,ismear,iwfn) ! gamma=1
            !!! the following 3 lines are for compatibility with fitpack...
            S0=DLOG(Q02/XC(1)**2)
            S=DLOG(DLOG(Q2/XC(1)**2)/S0)
            call fsupdf(0,x,s,u,d,ub,db,sb,cb,bb,glue)
            !!! if working with stand-alone package, could use this 1-line:
            !!! call PDFsa(x,Q2,U,D,UB,DB,SB,CB,BB,GLUE) 
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = PHI_INT2D (yD,1d0,10,ismear,iwfn)
               call set_offshell_on (.true.)
               call fsupdf(0,x,s,uoff,doff,uboff,dboff,
     &              sboff,cboff,bboff,glueoff)
               call set_offshell_on (.false.)
            endif
         else
            fy_diag = 0d0
            FFN = 0d0
            fy_off = 0d0
            F2N = 0d0
            if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
               ! PDF-level off shell models e.g., (fm)KP, free par
               foff_diag = 0d0
               foff_off = 0d0
               FFNoff = 0d0
               F2Noff = 0d0
            endif
         end if

         if (v(1).lt.0.995d0) then
            weight = .5*(yDmax-yDmin)*WI(I)
            if (istruc.ne.-1) then ! Structure function convolution
               FFd = FFd + weight * (fy_diag*FFN + fy_off*F2N)
               if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
                  ! PDF-level off shell models e.g., (fm)KP, free par
                  FFdoff = FFdoff 
     &                 + weight*(foff_diag*FFNoff+foff_off*F2Noff)
               endif
            else                   ! PDF convolution
               UA    = UA    + weight * fy_diag*U
               DA    = DA    + weight * fy_diag*D
               UBA   = UBA   + weight * fy_diag*UB
               DBA   = DBA   + weight * fy_diag*DB
               SBA   = SBA   + weight * fy_diag*SB
               CBA   = CBA   + weight * fy_diag*CB
               BBA   = BBA   + weight * fy_diag*BB
               GLUEA = GLUEA + weight * fy_diag*GLUE
               if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
                  ! off-shell PDF convolution
                  UAo    = UAo    + weight * foff_diag*Uoff
                  DAo    = DAo    + weight * foff_diag*Doff
                  UBAo   = UBAo   + weight * foff_diag*UBoff
                  DBAo   = DBAo   + weight * foff_diag*DBoff
                  SBAo   = SBAo   + weight * foff_diag*SBoff
                  CBAo   = CBAo   + weight * foff_diag*CBoff
                  BBAo   = BBAo   + weight * foff_diag*BBoff
                  GLUEAo = GLUEAo + weight * foff_diag*GLUEoff
               endif
            end if
         end if
         
         !print*, '* fy_diag _off =',fy_diag, fy_off
         !print*, '* FFN, F2N =', FFN,F2N
         !print*, '* FFd =',FFd
         !stop

      ENDDO
      if (istruc.eq.1) then       ! convolution done for x*F1_A, returns F1_A
         FFd = FFd/x
         if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
            ! PDF-level off shell models e.g., (fm)KP, free par
            FFdoff = FFdoff/x
         endif
      else if (istruc.eq.-1) then ! nPDFs are output in common/nPDF/  
         FFd = -9999999d0
         if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
            FFdoff = 0d0
         endif
      end if

*    *** Off-shell corrections

      offsh=0d0

*    ... Parametrized off-shell corrections (F2 only!!)
      if (ioff.eq.2) then
         ! Off-shell corrections by Melnitchouk-Schreiber-Thomas
         ! (parametrization by W.Melnitchouk) 
         offsh = OFF_MST(x)

      else if (ioff.eq.3) then
         ! approximate Kulagin-Petti NPA(2007) 
         ! (parametrization by W.Melnitchouk) 
c         offsh = OFF_KP(x)

         ! free floating parametrization at F2 level
         offsh = OFF_KP(x,xc,ndrv)

      else if (ioff.eq.4) then
         ! modified Kulagin-Petti (by W.Melnitchouck)
         ! parametrized and applied at F2 level for speed
         if (istrg.eq.0) then
            dRN = 10             ! max nucleon swelling
         else
            dRN = istrg
         end if
         offsh = OFF_mKP_fit(x,iwfn,dRN)

*    ... PDF-level models (can be used with any observable)
      else if (ioff.eq.1.or.(ioff.ge.5.and.ioff.le.9)) then
         ! ioff = 5,6  fmKP precalculated dq(x) - fixed strength levels
         ! ioff = 7    fmKP dynamically calculated dq(x) 
         ! ioff = 8    KP-like (parametrized dq(x)) dynamically calculated
         ! ioff = 9    KP (fully paraetrized)
         ! ioff = 1    freepar : generic parametrization   
         offsh = 0d0 ! offshell corrections handled at PDF level
         FFd=FFd+FFdoff ! add the offshell part to the onshell part
         if (istruc.eq.-1) then ! offshell correction to the pdfs
            UA=UA+UAo
            DA=DA+DAo
            UBA=UBA+UBAo
            DBA=DBA+DBAo
            SBA= SBA + SBAo
            CBA= CBA + CBAo
            BBA= BBA + BBAo
            GLUEA = GLUEA + GLUEAo
            write(*,*) x,UAo/UA,DAo/DA,UBAo/UBA,DBAo/DBA,GLUEAo/GLUEA
         endif

*    ... No offshell corrections
      else if (ioff.eq.0)then
         offsh=0.d0

      else
         write(*,*) 'ERROR (smearf2-2): ioff out of range = ',ioff
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

      character outfile*150, uflag*100

      integer length,oldlen,outlen,ismear,iwfn,ioff
      integer itoy, read_appl, ipv, iny
*    *** common blocks 
      COMMON/commandline/ itoy, read_appl, ipv, iny

      length = 4
      outfile = 'phi.'
      oldlen = length+1

c$$$      if (ismear.eq.2) then
c$$$         length = oldlen + 1
c$$$         outfile(oldlen:length) = 'km'
c$$$         oldlen = length + 1
c$$$      else 
      if (ismear.eq.3) then
         length = oldlen + 2
         outfile(oldlen:length) = 'wba'
         oldlen = length + 1
      else if (ismear.eq.4) then
         length = oldlen + 5
         outfile(oldlen:length) = 'wbarel'
         oldlen = length + 1
      else if (ismear.eq.5) then
         length = oldlen + 5
         outfile(oldlen:length) = 'aqv'
         oldlen = length + 1
      else if (ismear.eq.6) then
         length = oldlen + 5
         outfile(oldlen:length) = 'aqvc'
         oldlen = length + 1
      else if ((ismear.ne.0).and.(ismear.ne.1)) then
         length = oldlen + 7
         outfile(oldlen:length) = '_smear??'
         oldlen = length + 1
      end if

      if (iwfn.eq.0) then
         length = oldlen + 5
         outfile(oldlen:length) = '_paris'
         oldlen = length + 1
      else if (iwfn.eq.1) then
         length = oldlen + 4
         outfile(oldlen:length) = '_AV18'
         oldlen = length + 1
      else if (iwfn.eq.2) then
         length = oldlen + 6
         outfile(oldlen:length) = '_CDBONN'
         oldlen = length + 1
      else if (iwfn.eq.3) then
         length = oldlen + 4
         outfile(oldlen:length) = '_WJC1'
         oldlen = length + 1
      else if (iwfn.eq.4) then
         length = oldlen + 4
         outfile(oldlen:length) = '_WJC2'
         oldlen = length + 1
      else if (iwfn.eq.5) then ! shujie: add AV18 calcualted from KP code
         length = oldlen + 9
         outfile(oldlen:length) = '_AV18_CALC'
         oldlen = length + 1
      else if (iwfn.eq.6) then ! shujie: add AV18 calcualted from KP code
         length = oldlen + 4
         outfile(oldlen:length) = '_test'
         oldlen = length + 1
      else 
         length = oldlen + 5
         outfile(oldlen:length) = '_wfn??'
         oldlen = length + 1
      end if

      if (ioff.eq.1) then
         length = oldlen + 4
         outfile(oldlen:length) = '_KPoff'
         oldlen = length + 1
      else if (ioff.ne.0) then
         length = oldlen + 7
         outfile(oldlen:length) = '_offsh??'
         oldlen = length + 1
      end if

      if (ipv.ne.1200.D0) then
        write(uflag, '(I8)') ipv
        outfile = trim(outfile)//'_'//trim(adjustl(uflag))
        length  = len(trim(outfile))
      end if
      if (iny.ne.1000.D0) then
        write(uflag, '(I8)') iny
        outfile = trim(outfile)//'_y'//trim(adjustl(uflag))
        length  = len(trim(outfile))
      end if
      outlen = length

      return
      end




*******************************************************
      subroutine set_nuke(jnuke)
*     stores the nuclear model code in a common block

      implicit none
      integer jnuke

      integer inuke
      common/nuke/inuke
      
      inuke = jnuke

      return 
      end

*******************************************************
      subroutine get_nuke(jnuke)
*     retrieves the nuclear model code from a common block

      implicit none
      integer jnuke

      integer inuke
      common/nuke/inuke
      
      jnuke = inuke

      return 
      end

         
