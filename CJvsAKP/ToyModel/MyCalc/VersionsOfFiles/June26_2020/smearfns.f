C ***********************************************************************
      subroutine smearfn_1D (iwfn,ibar,rel,yD,gamma,phi0,phi1)
C
C     Nucleon fractional momentum distribution function in deuteron in
C     weak binding approximation (WBA) or AQV formalism, as a function of 
C     the fraction (yD) of deuteron's momentum carried by nucleon, defined as
C        yD = p.q/p_D.q
C     so that in Bjorken limit yD -> y0D is in [0,1].
C     
C     Kulagin, Petti, NPA 765, 126 (2006), Eq. (43);
C     Kahn, WM, Kulagin, PRC 79, 035205 (2009).
C     Accardi, Overleaf notes (2020): www.overleaf.com/project/5ec054f4a8a9cc0001e9abc8
C     
C     Relativistic kinematics implemented (W.M.): May 2010.
C     Relativistic convolution implemented (cf. AQV): Sep. 2010
C     Modified for CTEQ-JLab package (A.A): Mar 2011
C     Revised after AKP-CJ benchmark (A.A): Jun 2020 
C     
C     Note: All momenta in MeV.
*
*     INPUT:
*
*     iwfn   = (i) what wave function (see 'DN' for a legend)
*     ibar   = (i) Baryon number normalization: 0=M/p0  1=const
*                  Note: only for AQV 
*     rel    = (i) relativ. corrections 0=WBA  1=WBAREL  2=AQV
*     yD     = (dp) nucleon fractional momentum
*     gamma  = (dp) gamma factor (NOTE: not gamma^2)
*
*     OUTPUT:
* 
*     phi0(4) = (dp) smearing functions -- integrals of S(p))
*     phi1(4) = (dp) "1st moments" -- integrals of (p^2/m_N^2-1)*S(p) 
*
C ***********************************************************************
      IMPLICIT NONE
      INTEGER iwfn,ibar,rel,ipT,npt,i,j,isimp
      REAL*8  phi0(4),phi1(4),yD,gamma
      REAL*8  FYINT(4),FYOFF(4),yDmax,pT,pTmin,pTmax,pvcut,pTlim,pTint
     %     ,fyp0(3),fyp1(3)
      REAL*8  pi,hc,mN,MD,epsD
      COMMON  /con/ pi,hc
      COMMON  /mas/ mN,MD,epsD
      ! user-definable pv=|\vect p| cut [MeV]
      integer*8 ipv
      common/flag/ ipv


      do i=1,4         
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

      pvcut = 1.D0 * ipv        ! pT<pv_cut)

      pTmin = 0.D0
      npT = 2000
      
      if (rel.eq.0) then        ! WBA (non relativistic)
         yDmax = (mN/MD)*(1.D0 + gamma**2/2 + epsD/mN) ! not exact
         ptmax = dsqrt(2*mN**2*(-MD/mN*yD+1+gamma**2/2+epsD/mN)) ! adaptive grids for pT
         if (ptmax.gt.pvcut) ptmax = pvcut
      else                      ! rel. kin => yDmax = infty
         yDmax = 1.D0    ! use 1.0 as approx. to infty
         pTmax = pvcut
      end if
      IF (yD.LE.0.D0 .OR. yD.GE.yDmax) RETURN 

      pTint = pTmax/npt        ! integration step
      
      ipT = 0
C...Integrate over nucleon transverse momentum squared [d(pT^2) = 2 pT dpT]
      DO j = 0,npT
         pT = j*pTint + pTmin

C...Deuteron-nucleon-nucleon vertex
         CALL DN (iwfn,ibar,rel,yD,gamma,pT,FYINT,FYOFF)

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

      do i = 1,4
         PHI0(i) = (pTint/3) * PHI0(i)
         PHI1(i) = (pTint/3) * PHI1(i)
      end do

      RETURN
      END


C *********************************************************************
      SUBROUTINE DN (wfn,ibar,rel,yD,gamma,pT,FYINT,FYOFF)
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
*     yD       = (dp) nucleon fractional momentum
*                     NOTE: yD=p.q/PD.q with PD the unrescaled Deuteron momentum
*     gamma    = (dp) gamma factor (NOTE: not gamma^2)
*     pT [MeV] = (dp) nucleon transverse momentum (unrescaled)
*
*     OUTPUT:
* 
*     FYINT    = (dp) Integrand for the smearing function
*     FYOFF    = (dp) (p^2-mN^2)/mN^2 * FYINT
*
C *********************************************************************
      IMPLICIT NONE
      INTEGER wfn,ibar,imom,i
      REAL*8 yD,gamma,pT,FYINT(4),FYOFF(4)
      REAL*8 yDmax,pz,pv,pv_cut,Es,p0,p2,eps
     &     ,Jac,flux,rhoN,rhoNT,corr,bnorm
      REAL*8 U,W,VS,VT,rho,CC,prefact
      REAL*8 pi,hc,mN,MD,epsD
      integer rel
      COMMON /con/ pi,hc
      COMMON /mas/ mN,MD,epsD
      integer*8 ipv
      common/flag/ ipv
      

      do i = 1, 4
         FYINT(i) = 0.D0
         FYOFF(i) = 0.D0
      end do

C...Kinematics...........................................................

      IF (rel.eq.0) THEN	! non-relativistic (sqrt in pz > 0)

         yDmax = (1.D0 + gamma**2/2 + epsD/mN - pT**2/(2*mN**2))*(mN/MD)
         pz = mN * ( gamma - DSQRT(2*MD/mN*(yDmax-yD)) )

      ELSE                      ! relativistic 

         yDmax = 1d10            ! proxy for yDmax = \infty
         IF (gamma.GT.1.D0) THEN
	    pz = ( DSQRT( (1.D0-yD)**2*MD**2
     &           + (gamma**2-1.D0)*(pT**2 + mN**2) )
     &		 - MD*(1.D0-yD)*gamma )
     &           / (gamma**2-1.D0)
         ELSE
	    pz = (pT**2 + mN**2 - (1.D0-yD)**2*MD**2) / (2*(1.D0-yD)*MD)
         ENDIF
         
      ENDIF

      IF (yDmax.LE.0.D0 .OR. yD.GE.yDmax) RETURN

      
C...Total nucleon 3-momentum
	pv = DSQRT(pT**2 + pz**2)
C...Cut-off integration at pv=pv_cut
C...Added command line flag to vary pv_cut values--SL, 05.2020
	pv_cut = ipv !1200.D0 
	IF (pv.GT.pv_cut) RETURN

C...On-shell spectator nucleon energy
        if (rel.eq.0) then
           Es = mN + pv**2/(2*mN)
        else 
           Es = DSQRT(mN**2 + pv**2)
        end if

C...Separation energy
        eps = MD - Es - mN  ! valid for both rel and non-rel options
        
C...Off-shell interacting nucleon energy
	p0 = MD - Es
!	IF (p0.LE.0.D0) RETURN		!!! for MST only !!!

C...Nucleon virtuality
	p2 = p0**2 - pv**2
        
!C...light-cone variable --- NO LONGER NEEDED [AA: 25 June 2020]
!C...solution of y = y0/2 (1 + gamma - (gamma-1)/y0^2 (p2+pT2)/M^2
!	IF (p2+pT**2.LT.-y**2*MD**2/(gamma**2-1.D0)) RETURN ! checks sign of sqrt argument below
!	y0 = y				!!! corrected 4/21/11
!     &	   * ( 1.D0
!     &	     + DSQRT( 1.D0 + (gamma**2-1.D0)*(p2+pT**2)/MD**2/y**2 ) )
!     &     / (1.D0+gamma) 


C...Deuteron wavefunction................................................

!...  Calculates (u^2+v^2+w^2) = |psi|^2 / (4\pi)
!...(--> the 1/(4\pi) term is included in the Jacobian below)  
!...NR:  1=Paris, 3=CDBonn, 4=AV18
!...Rel: 12=WJC-1, 13=WJC-2
!...[high-precision ones are: CDBonn, AV18, WJC-1, WJC-2]
        
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
         CALL WJC (2,pv/hc,U,W,VS,VT)   ! -- 2% in V' term 

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
      else IF (wfn.EQ.5) then 
         imom = 1
         CALL AV18_CALC (pv/hc,rho)
      else IF (wfn.EQ.6) then 
         imom = 0
         CALL TESTWF (pv/hc,U,W)
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


C...Nucleon distributions in yD=p.q/pA.q 
C...(building blocks for smearing functions)

C    ...flux factors and finite Q^2 corrections
      IF (rel.eq.2) THEN        ! AQV
C       ...Baryon number normalization factor
         if (ibar.eq.0) then    
            bnorm = mN/p0 ! -> for the future: onsider mN/(p0+gamma*pz) = 1/y, instead
         else                   ! taken care by normalizing for gamma=1
            bnorm = 1d0         ! in the calling routine
         end if
         flux = yD*MD/mN * bnorm 
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &	     * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2d0/mN**2 )
         ! i.e., corr = (gamma^2-1)/y^2 * (2p^2+3pT^2)/(2mN^2)
      else if (rel.eq.1) then   ! WBAREL
         if (ibar.ne.0) then
            write(*,*) 'ERROR(SMEARFNS): ibar out of bounds for rel=1:'
     &           , ibar
            stop
         end if
         ! bnorm=mN/p0 already included in flux factor
         flux = (1.D0 + gamma*pz/mN)    ! neglects O(eps_D/mN) compared to AQV
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &          * ( (1.D0 + eps/mN)**2 + (pv**2-3*pz**2)/2d0/mN**2 )
         ! corr is same as AQV 
      ELSE                      ! WBA
         if (ibar.ne.0) then
            write(*,*) 'ERROR(SMEARFNS): ibar out of bounds for rel=0:'
     &           , ibar
            stop
         end if
         flux = (1.D0 + gamma*pz/mN)  ! same as in WBAREL
         corr = (gamma**2-1.D0)/(yD*MD/mN)**2
     &        * (1.D0 + 2*eps/mN + (pv**2-3*pz**2)/(2*mN**2))
         ! corr is same as WBAREL, but expanded to O(eps/mN)=O(pvec^2/mN^2)
      ENDIF

C...  "Jacobian" .....................................
!...  d3p -> dy*dpT2 & normalization of |u^2+v^2| read in from file
      prefact =  gamma*Es / (4*((1-yD)+(gamma**2-1)*Es/MD))


C...  FINAL CALCULATION


C...On-shell
      FYINT(1) = prefact        ! diagonal for 2xF1, xFL
     &         * flux * CC
      FYINT(2) = prefact        ! offdiagonal for 2xF1, xFL
     &	       * flux * CC
     &         * (gamma**2-1d0)/(yD*MD/mN)**2 * (pT/mN)**2
      FYINT(3) = prefact        ! F2 smearing function (diagonal for F2)
     &         * flux * CC
     &         * (1.D0 + corr)/gamma**2  
      FYINT(4) = prefact	! xF3 smearing function (diagonal for xF3)
     &         * (1.D0 + pz/mN/gamma) * CC
               ! ^^^^^^^^^^^^^^^^^^^^ CHECK implementation of bnorm here!! [AA 25 June 2020]

C...Off-shell weighted part
      do i = 1, 4
         FYOFF(i) = FYINT(i) * (p2-mN**2)/mN**2
      end do
         
      RETURN
      END
