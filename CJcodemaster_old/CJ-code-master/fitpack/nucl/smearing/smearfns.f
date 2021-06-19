C ***********************************************************************
      subroutine smearfn_1D (iwfn,ibar,rel,y,gamma,phi0,phi1)
C
C     Nucleon fractional momentum distribution function in deuteron in
C     weak binding approximation (WBA) or AQV formalism, as a function of 
C     the fraction (y) of deuteron's momentum carried by nucleon, defined as
C        y = p.q/p_D.q
C     (sometimes labeled as y_D) so that in Bjorken limit y -> y0 in [0,1].
C     
C     Kulagin, Petti, NPA 765, 126 (2006), Eq. (43);
C     Kahn, WM, Kulagin, PRC 79, 035205 (2009).
C     
C     Relativistic kinematics implemented (W.M.): May 2010.
C     Relativistic convolution implemented (cf. AQV): Sep. 2010
C     Modified for CTEQ-JLab package (A.A): Mar 2011
C  
C     Note: All momenta in MeV.
*
*     INPUT:
*
*     iwfn   = (i) what wave function (see 'DN' for a legend)
*     ibar   = (i) Baryon number normalization: 0=M/p0  1=const
*                  Note: inly for AQV 
*     rel    = (i) relativ. corrections 0=WBA  1=WBAREL  2=AQV
*     y      = (dp) nucleon fractional momentum
*     gamma  = (dp) gamma factor (NOTE: not gamma^2)
*
*     OUTPUT:
* 
*     phi0(4) = (dp) smearing function -- integral of S(p))
*     phi1(4) = (dp) "1st moment" -- integral of p^2/m_N^2-1)*S(p) 
*
C ***********************************************************************
      IMPLICIT NONE
      INTEGER iwfn,ibar,ipT,npt,i,j,isimp
      REAL*8  phi0(3),phi1(3),y,gamma
      REAL*8  FYINT(3),FYOFF(3),ymax,pT,pTmin,pTmax,pTlim,pTint
     %     ,fyp0(3),fyp1(3)
      REAL*8  pi,hc,mN,MD,epsD
      COMMON  /con/ pi,hc
      COMMON  /mas/ mN,MD,epsD
      integer rel  


      do i=1,3         
         PHI0(i) = 0.D0
         PHI1(i) = 0.D0
      end do

C...Constants
      pi = 4*DATAN(1.D0)
      hc = 197.327D0 		! conversion factor (MeV.fm)
      mN = 938.91897D0          ! nucleon mass
      epsD = -2.224575D0	! deuteron binding energy
      MD = 2*mN + epsD          ! deuteron mass
      
C...Kinematic limits
      if (rel.eq.0) then        ! WBA (non relativistic)
         ymax = (mN/MD)*(1.D0 + gamma**2/2 + epsD/mN)
      else                      ! rel. kin => ymax = infty
         ymax = 1.D0		! use 1.0 as approx. to infty
      end if
      IF (y.LE.0.D0 .OR. y.GE.ymax) RETURN 

      pTmin = 0.D0
      pTmax = 2000.D0           ! INTEGRATION LIMITS
      npT = 2000
      pTint = pTmax/npt        ! integration step
      
      ipT = 0
C...Integrate over nucleon transverse momentum squared [d(pT^2) = 2 pT dpT]
      DO j = 0,npT
         pT = j*pTint + pTmin

C...Deuteron-nucleon-nucleon vertex
         CALL DN (iwfn,ibar,rel,y,gamma,pT,FYINT,FYOFF)

*C...Nucleon momentum distribution in deuteron [d(pT^2) => d(pT) 2 pT]
*         fyp0(2) = FYINT(2) * ( 2*pT )
*         fyp1(2) = FYOFF(2) * ( 2*pT )

C...Simpson's rule integration [d(pT^2) => d(pT) 2 pT]
         isimp = 2 + 2*(ipT - ipt/2*2)
         do i=1,3
            PHI0(i) = PHI0(i) + isimp*2*pT*fyint(i)
            PHI1(i) = PHI1(i) + isimp*2*pT*fyoff(i)
         end do
         ipT = ipT + 1
      ENDDO

      do i = 1,3
         PHI0(i) = (pTint/3) * PHI0(i)
         PHI1(i) = (pTint/3) * PHI1(i)
      end do

      RETURN
      END


C *********************************************************************
      SUBROUTINE DN (wfn,ibar,rel,y,gamma,pT,FYINT,FYOFF)
C
C     Integrand of WBA and AQV nucleon momentum distributions
C     (for nonrelativistic and relativistic kinematics).
C     
C     Includes also the p^2 weigted (e.g., for the KP off-shell model)
C     
C     Expressions valid only in the deuteron rest frame.
*
*
*     INPUT:
*
*     wfn      = (i) what wave function, see below for explanation
*     ibar     = (i) Baryon number normalization: 0=M/p0  1=const
*                    Note: inly for AQV 
*     rel      = (i) relativ. corrections 0=WBA  1=WBAREL  2=AQV
*     y        = (dp) nucleon fractional momentum
*     gamma    = (dp) gamma factor (NOTE: not gamma^2)
*     pT [MeV] = (dp) nucleon transverse momentum
*
*     OUTPUT:
* 
*     FYINT    = (dp) Integrand for the smearing function
*     FYOFF    = (dp) (p^2-mN^2)/mN^2 * FYINT
*
C *********************************************************************
      IMPLICIT NONE
      INTEGER wfn,ibar,imom,i
      REAL*8 y,gamma,pT,FYINT(4),FYOFF(3)
      REAL*8 ymax,y0,pz,pv,pv_max,Ep,p0,p2,eps
     &     ,Jac,flux,rhoN,rhoNT,corr,bnorm
      REAL*8 U,W,VS,VT,rho,CC,prefact
      REAL*8 pi,hc,mN,MD,epsD
      integer rel
      COMMON /con/ pi,hc
      COMMON /mas/ mN,MD,epsD

      do i = 1, 3
         FYINT(i) = 0.D0
         FYOFF(i) = 0.D0
      end do

C...Kinematics...........................................................
      IF (rel.eq.0) THEN	! non-relativistic

         ymax = (1.D0 + gamma**2/2 + epsD/mN - pT**2/(2*mN**2))*(mN/MD)
         IF (ymax.LE.0.D0 .OR. y.GE.ymax) RETURN
         pz = mN * ( gamma - DSQRT(2*MD/mN*(ymax-y)) )

      ELSE                      ! nonrelativistic (sqrt in pz > 0)

         IF (gamma.GT.1.D0) THEN
	    pz = ( DSQRT( (1.D0-y)**2*MD**2
     &           + (gamma**2-1.D0)*(pT**2 + mN**2) )
     &		 - MD*(1.D0-y)*gamma )
     &           / (gamma**2-1.D0)
         ELSE
	    y0 = y
	    pz = (pT**2 + mN**2 - (1.D0-y0)**2*MD**2) / (2*(1.D0-y0)*MD)
         ENDIF
         
      ENDIF

C...Total nucleon 3-momentum
	pv = DSQRT(pT**2 + pz**2)

C...On-shell spectator nucleon energy
	Ep = DSQRT(mN**2 + pv**2)

C...Off-shell interacting nucleon energy
	p0 = MD - Ep
!	IF (p0.LE.0.D0) RETURN		!!! for MST only !!!

C...Nucleon virtuality
	p2 = p0**2 - pv**2

C...Check that sqrt in y0 is non-negative
	IF (p2+pT**2.LT.-y**2*MD**2/(gamma**2-1.D0)) RETURN

C...Cut-off integration at pv=pv_max
	pv_max = 1200.D0
	IF (pv.GT.pv_max) RETURN

C...Separation energy
	IF (rel.eq.0) THEN
	  eps = epsD - pv**2/(2*mN)
	ELSE
	  eps = MD - Ep - mN
	ENDIF

C...light-cone variable 
C...solution of y = y0/2 (1 + gamma - (gamma-1)/y0^2 (p2+pT2)/M^2
	y0 = y				!!! corrected 4/21/11
     &	   * ( 1.D0
     &	     + DSQRT( 1.D0 + (gamma**2-1.D0)*(p2+pT**2)/MD**2/y**2 ) )
     &     / (1.D0+gamma) 

C...Jacobian for delta functions 
	Jac = gamma


C...Deuteron wavefunction................................................
C...NR:  1=Paris, 3=CDBonn, 4=AV18
C...Rel: 12=WJC-1, 13=WJC-2
C...[high-precision ones are: CDBonn, AV18, WJC-1, WJC-2]

      U = 0.D0                  ! S-wave
      W = 0.D0                  ! D-wave
      VS = 0.D0                 ! Singlet P-wave
      VT = 0.D0                 ! Triplet P-wave

      ! momenta in MeV

      !relativistic wave functions
      if (wfn.EQ.12) then
         imom = 0                       ! wfns normalized to ~105%
         CALL WJC (1,pv/hc,U,W,VS,VT)   ! -- 5% in V' term
      else IF (wfn.EQ.13) then 
         imom = 0                       ! wfns normalized to ~102%
         CALL WJC (2,pv/hc,U,W,VS,VT)  ! -- 2% in V' term 

      !Nonrelativistic wave functions (VT=VS=0)
      else IF (wfn.EQ.1) then 
         imom = 0
         CALL PARIS (1,pv/hc,U,W)
      else IF (wfn.EQ.3) then 
         imom = 0
         CALL CDBONN (pv/hc,U,W)
      else IF (wfn.EQ.4) then 
         imom = 1
         CALL AV18 (pv/hc,rho)
      end if

      !Output wavefunctions in fm^3/2 => MeV^-3/2
      if (imom.eq.0) then   
         U = U / hc**1.5D0
         W = W / hc**1.5D0
         VS = VS / hc**1.5D0
         VT = VT / hc**1.5D0
         CC = (U**2 + W**2 + VS**2 + VT**2)
      else                   !Wave fns in terms of momentum distributions
         CC = rho / hc**3.D0
      end if


C...Nucleon distributions in y=p.q/pA.q 
C...(building blocks for smearing functions)

C    ...flux factors and finite Q^2 corrections
      IF (rel.eq.2) THEN        ! AQV
C       ...Baryon number normalization factor
         if (ibar.eq.0) then    
            bnorm = mN/p0
         else                   ! taken care by normalizing for gamma=1
            bnorm = 1d0         ! in the calling routine
         end if
         flux = y*MD/mN * bnorm 
         corr = (gamma**2-1.D0)/(y*MD/mN)**2
     &	     * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2d0/mN**2 )
         ! i.e., corr = *gamma^2-1)/y^2 * (2p^2+3pT^2)/(2mN^2)
      else if (rel.eq.1) then   ! WBAREL
         if (ibar.eq.0) then
            flux = (1.D0 + gamma*pz/mN)
         else
            flux = (1.D0 + (eps+gamma*pz)/mN)
         end if
         corr = (gamma**2-1.D0)/(y*MD/mN)**2
     &	     * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2d0/mN**2 )
         ! same as AQV
      ELSE                      ! WBA
         if (ibar.eq.0) then
            flux = (1.D0 + gamma*pz/mN)
         else
            flux = (1.D0 + (eps+gamma*pz)/mN)
         end if
         flux = (1.D0 + gamma*pz/mN)
         corr = (gamma**2-1.D0)/(y*MD/mN)**2
     &        * (1.D0 + 2*eps/mN + (pv**2-3*pz**2)/(2*mN**2))
         ! same as WBAREL, but expanded to O(pvec^2)
      ENDIF

C...WBA at finite Q^2
      prefact =  Ep*mN/MD / (2*Jac*(1.D0-y0)) ! Jacobian etc.
ccc     &     * flux                    ! flux factor
     &     * CC                      ! wave function

      FYINT(1) = prefact        ! diagonal for 2xF1, xFL
     &	       * flux
      FYINT(2) = prefact        ! offdiagonal for 2xF1, xFL
     &	       * flux
     &         * (gamma**2-1d0)/(y*MD/mN)**2 * (pT/mN)**2
      FYINT(3) = prefact        ! F2 smearing function (diagaonal for F2)
     &	       * flux
     &         * (1.D0 + corr)/gamma**2  
      FYINT(4) = prefact	! xF3 smearing function (diagonal for xF3)
     &	       * (1.D0 + pz/mN/gamma)

C...Off-shell weighted part
      do i = 1, 3
         FYOFF(i) = FYINT(i) * (p2-mN**2)/mN**2
      end do
         
      RETURN
      END
