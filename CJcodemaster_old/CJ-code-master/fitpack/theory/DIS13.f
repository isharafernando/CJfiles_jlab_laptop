      SUBROUTINE strfn(ibeam,istruc,itm,iht,ndrv,v,xc,sfn)
*     (28 Nov 08) A.Accardi
*     Calculates DIS structure functions or cross section including 
*     Target Mass Corrections and Higher-Twist terms
*     It is a wrapper that calls "theory2" with variables depending on TMC, 
*     and applies TMC kinematic factors, and HT corrections on top of that.
*
*     NOTE: the target's number of protons and neutrons (Z,N) must be 
*           set previously to calling strfn in 'common/target/az,an' 
*           where az=Z/(Z+N) and an=N/(Z+N). This can be accomplished by 
*           calling 'targ(az,an)' without using the common block in the 
*           calling routine.
*
*     NOTE: theory2 is called with its argument itarg=2, so that proper 
*           weighting of protons and neutrons is performed. This trick 
*           was needed for backwards compatibility.  
*
*     INPUT:
*
*       ibeam     = (i)  beam type: 1=e 2=nu 3=nubar 4=e with gamma* Z included
*                                   5 = e reduced cross section for charm
*                                   6 = e interference terms for A_pv
*       istruc    = (i)  structure fn: 0=FL 1=F1 2=F2 3=xF3  
*       itm       = (i)  TMC  0=no  1=approx GP  2=CF
*       iht       = (i)  HT   0=no  
*                             1= multiplicative: F2 = LT*(1+C(xB)/Q^2)
*                             2= additive: F2 = LT + C(xB)/Q^2
*                             11-12 = SAME BUT C=C(xi)
*                        NOTE NOTE: additive only experimental: correct
*                                   HT parametrization not yet explored
*       ndrv      = (i)  needed to determine the proper index in the pdf array
*                        when derivatives are calculated during the fitting
*       v(4)      = (dp) variables - for DIS v(1)=xB  v(2)=Q^2;  
*                        v(3) is unused for str.fns. but in the fitting
*                        code is = W^2 or y depending if the data is for
*                        str.fns. or cross sections; v(4) unused
*       xc(100)   = (dp) current fit parameters
*
*     OUTPUT:
*
*       sfn        = (dp) The required structure function:
*                         FL, F1, F2, x*F3
*
*     EXAMPLE calling sequence for elm proton F2. The first example requires
*             the common block 'common/target/az,an' with double precision 
*             az,an in the calling routine; the second does not. 
*
*             az = 1.          ! fraction of protons 
*             an = 0.          ! fraction of neutrons
*             call strfn(1,2,NDRV,V,XC,sfn)
*
*             OR
*
*             call targ(1.,0.) 
*             call strfn(1,2,NDRV,V,XC,sfn)
*
*     TO DO LIST (2 Sep 10)
*
*       - HT for F1, F3, FL  (for the moment they are == HT for F2)
*       - cross sections
*
      implicit none
      
      integer ibeam,itarg,istruc,ndrv,itm,iht,ist,i
      double precision v(4),xc(100),sfn,add

      double precision vp(4),x,q2,xtmc,amu,r,fLcor,f1cor,f2cor,f3cor,xht
     &     ,ans,ansf2,ansf1,ht_f1,ht_f2,ht_f3,ht_fL,C(0:7),HT(0:7)

*     Higher-twist parameters
      integer nht(4),nhtn(4)
      common/ht/nht,nhtn
*     FL Higher-twist parameters      
      integer nhtFL(4),nhtFLn(4)
      common/htFL/nhtFL,nhtFLn

*     Proton and neutron numbers in the target
      double precision az,an
      common/target/az,an

*      integer itype(35),nflag(35),inorm(35),itgt(35),icorr(35),itmc(35)
*     &     ,iht(35)
*      common/flags/itype,nflag,inorm,itgt,icorr,itmc,iht

*    *** Sets up TMC
*    ... [AA 080519] inclusion of DIS TMC in col.fact. Modified the 
*    ... meaning of VP(3), instead of VP(3)=V(3)=W^2 as previously
*    ... done. Not a problem yet for "THEORY2" subroutine because
*    ... W^2 is not used in the str.fn. computations. 
*    ... [AA 120518] keeps VP(3)=V(3)=W^2,y because it is used in DIS
*    ... cross section computations. Instad, uses VP(4) to perform TMCs
*    ... since V(4)=data type is needed only in global fits by THEORY10.f 
*    ... to select the appropriate calculation for a given data point, but
*    ... is not needed by the structure function subroutines. In fact, the 
*    ... DIS kinematics is completely specified by x,Q^2,y = V(1),V(2),V(3)

      x=v(1)
      q2=v(2)
      amu=.8836/q2
      r=sqrt(1.+4.*x**2*amu)
      xtmc=2.*x/(1.+r)

      vp(2) = q2
      if(itm.eq.0)then
*       *** No TMC: 
*       ... vp(1)=xB  vp(4)=xB
         vp(1) = x
         vp(4) = x
         fLcor=1.
         f1cor=1.
         f2cor=1.
         f3cor=1.
      else if(itm.eq.1)then
*       *** approximate GP TMC:
*       ... VP(1) = xi_Nachtmann, and VP(4) = xi_Nachtmann
         vp(1)=xtmc
         vp(4)=xtmc
*       ... f1cor & fLcor GP not yet implemented !!!!
c         fLcor= -1d0
c         f1cor = -1d0
c         f2cor=(x/xtmc)**2/r**3*(1.+6.*amu*x*xtmc/r*(1.-xtmc)**2)
c         f3cor=x/xtmc/r**2*(1.-amu*x*xtmc/r*(1.-xtmc)*dlog(xtmc))
c jfo
c jfo     GP TMC now calculated in routine GPtmc
c jfo     All coefficients taken care of there, so they are set equal
c jfo     to one here
c jfo
         fLcor=1d0
         f1cor=1d0
         f2cor=1d0
         f3cor=1d0
      else if (itm.eq.2) then
*       *** TMC in col.fact.
*       ... VP(1) = xi_Nachtmann, and VP(4) = xB
         vp(1)=xtmc
         vp(4)=x
         fLcor=x/xtmc
         f1cor=1.
         f2cor=x/(xtmc*r**2)
         f3cor=x/(xtmc*r)
      else 
*       *** xi-scaling (Kretzer-Reno) TMC in col.fact.
*       ... VP(1) = xi_Nachtmann, and VP(4) = xi_Nachtmann
         vp(1)=xtmc
         vp(4)=xtmc
         fLcor=x/xtmc
         f1cor=1.
         f2cor=x/(xtmc*r**2)
         f3cor=x/(xtmc*r)
      endif

*    *** Higher-twist term
*    ... iht  =  0      no HT
*    ... iht  =  1-10   F2(TMC) = F2(TMC) * (1+ht_f2(xB))
*    ... iht  = 11-20   F2(TMC) = F2(TMC) * (1+ht_f2(xi))
*                       with same ht_f2 as iht-10
      if (iht.le.10) then
         xht = x
      else 
         xht = xtmc
      end if

*    ... computes C(xB)/Q^2 for the F2 structure function
*    ... NOTE: C=C(ist) with  ist=istruc*2 + {0 for p, 1 for n}
*              ist=  0,1  for FL(p,n)   
*              ist = 2,3  for F1(p,n)
*              ist = 4,5  for F2(p,n)
*              ist = 6,7  for F3(p,n)

      if (iht.eq.0) then
         C(4) = 0d0
         C(5) = 0d0
         C(0) = 0d0
         C(1) = 0d0
      else if ((iht.eq.1).or.(iht.eq.11)) then
*    ... multiplicative parametrization
         C(4) = ( xc(nht(4))                       ! F2 proton
     &        + xc(nht(1)) * xht**xc(nht(2)) 
     &        * (1d0 + xc(nht(3))*xht) )
     &        / Q2
         C(5) = ( xc(nhtn(4))                      ! F2 neutron
     &        + xc(nhtn(1)) * xht**xc(nhtn(2)) 
     &        * (1d0 + xc(nhtn(3))*xht) )
     &        / Q2
         C(0) = ( xc(nhtFL(4))                       ! FL proton
     &        + xc(nhtFL(1)) * xht**xc(nhtFL(2)) 
     &        * (1d0 + xc(nhtFL(3))*xht) )
     &        / Q2
         C(1) = ( xc(nhtFLn(4))                      ! FL neutron
     &        + xc(nhtFLn(1)) * xht**xc(nhtFLn(2)) 
     &        * (1d0 + xc(nhtFLn(3))*xht) )
     &        / Q2
      else if ((iht.eq.2).or.(iht.eq.12)) then
*    ... additive parametrization --- Still experimental (June 2010)
         C(4) = xc(nht(1))                      ! F2 proton
     &        * xht**xc(nht(2)) 
     &        * (1d0 + xc(nht(3))*xht) 
     &        * (1d0-xht)**xc(nht(4))
     &        / Q2
         C(5) = xc(nhtn(1))                     ! F2 neutron
     &        * xht**xc(nhtn(2)) 
     &        * (1d0 + xc(nhtn(3))*xht) 
     &        * (1d0-xht)**xc(nhtn(4))
     &        / Q2
         C(0) = xc(nhtFL(1))                      ! FL proton
     &        * xht**xc(nhtFL(2)) 
     &        * (1d0 + xc(nhtFL(3))*xht) 
     &        * (1d0-xht)**xc(nhtFL(4))
     &        / Q2
         C(1) = xc(nhtFLn(1))                     ! FL neutron
     &        * xht**xc(nhtFLn(2)) 
     &        * (1d0 + xc(nhtFLn(3))*xht) 
     &        * (1d0-xht)**xc(nhtFLn(4))
     &        / Q2
      else
         print*, 'ERROR (strfn): iht out of range:', iht  
      end if
*    ... C(xB) for FL, F1 and F3 == C(xB) for F2, for the moment 
      C(2) = C(4)       ! F1 proton
      C(3) = C(5)       ! F1 neutron
      C(6) = C(4)       ! F3 proton
      C(7) = C(5)       ! F3 neutron
*    ... defines HT parameter to pass to THEORY2, and additive term
      if ((iht.eq.2).or.(iht.eq.12)) then                    
         do i=0,7
            HT(i) = 1d0         ! additive corrections
         end do
         ist = istruc*2
         add = az*C(ist)+an*C(ist+1)
      else if (iht.eq.0.or.iht.eq.1.or.iht.eq.11) then   
         do i=0,7
            HT(i) = 1d0+C(i)    ! mult. or no correction
         end do
         add = 0d0
      else 
         print*, 'ERROR(strfn): iht out of range:', iht
         stop
      end if

*    *** computes str.fn. with TMC+HT, if any

*    ... NOTE: theory2 with itarg=2 as below always corrects 
*    ... for proton/neutron ratio (the common vars az=Z/A and 
*    ... an=1-Z/A must always be set prior to any call, also for a proton)
*    ... The multiplicative HT correction is computed by theory2, 
*    ... the additive one below.


      if(itm.eq.1)then
         call GPtmc(ibeam,2,istruc,HT,ndrv,vp,x,r,xc,ans)
      else if(itm.eq.0)then
         if(istruc.eq.0)then
            call theory2(ibeam,2,1,HT,ndrv,vp,xc,ansf1)
            call theory2(ibeam,2,2,HT,ndrv,vp,xc,ansf2)
            ans=r**2*ansf2-2.*vp(1)*ansf1  
            !call theory2(ibeam,2,0,HT,ndrv,vp,xc,ans) ! <--- this is the "normal" one; above uses kineamtic coefficients, but if tmc=0 no TMC inside F2 and F1
         else
            call theory2(ibeam,2,istruc,HT,ndrv,vp,xc,ans)
         endif
      else if(itm.eq.2.or.itm.eq.3)then
         call theory2(ibeam,2,istruc,HT,ndrv,vp,xc,ans)
      else
         print*,'Incorrect value for itm'
         stop
      endif
*    ...  Final answer
      if (istruc.eq.0) then
         sfn = ans * fLcor + add
      else if (istruc.eq.1) then
         sfn = ans * f1cor + add
      else if (istruc.eq.2) then
         sfn = ans * f2cor + add
      else if (istruc.eq.3) then
         sfn = ans * f3cor + add
      end if

      return
      end

      subroutine GPtmc(ibeam,itarg,istruc,HT,ndrv,vp,x,r,xc,ftmc)
      implicit real*8 (a-h,o-z)
      dimension vp(4), xc(100), HT(0:7)
      data xm, xm2/0.94d0, .8836d0/
c jfo October 22, 2013
c jfo calculates approximate GP target mass corrections 
c jfo following arXiv:0709.1775[hep-ph]
c jfo retains subleading terms in FL
c jfo vp(1)=vp(4)=xi
      q2=vp(2)
      xi=vp(1)
      amu=xm2/q2
      call theory2(ibeam,2,istruc,HT,ndrv,vp,xc,ans)
      if(istruc.eq.0.or.istruc.eq.1)then
         call theory2(ibeam,2,2,HT,ndrv,vp,xc,ansf2)
      endif
      if(istruc.eq.0)then
c         ftmc=(x/xi)**2/r*ans*(1.+4.*amu*x*xi/r*(1.-xi)**2*ansf2/ans)
         ftmc=(x/xi)**2/r*(ans+4.*amu*x*xi/r*(1.-xi)**2*ansf2)
c         if ((ftmc/ans-1d0).gt.0.001) write(6,*) '* dFL =',x,Q2,ftmc/ans       
      else if(istruc.eq.1)then
         ftmc=x/xi/r*ans*(1.+amu*x/r*(1.-xi)**2*ansf2/ans)
c         ftmc=x/xi/r*(ans+amu*x/r*(1.-xi)**2*ansf2)
      else if(istruc.eq.2)then
         ftmc=(x/xi/r)**2/r*ans*(1.+6.*amu*x*xi/r*(1.-xi)**2)
      else if(istruc.eq.3)then
         ftmc=(x/xi/r)**2*ans*(1.-amu*x*xi/r*(1.-xi)*dlog(xi))
      else
         print*,'Incorrect value for istruc',istruc
         stop
      endif
      return
      end

      subroutine targ(dz,dn)
*     Selects the target's atomic and neutron numbers
**     Programmer: A.Accardi
*     Created: 18 Mar 2010
*     
*     INPUT
*
*       dz  = (dp) Fraction of protons  Z/(Z+N)
*       dn  = (dp) Fraction of neutrons N/(Z+N)
*
      implicit none
      
      double precision dz,dn

      double precision az,an
      common/target/az,an

      az = dz
      an = dn

      return
      end
      
      SUBROUTINE THEORY2(ibeam,itarg,istruc,HT,NDRV,V,XC,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(4),XC(100),FTEMP(2),HT(0:7)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
c      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/FLAVOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/Q2STUF/Q02,Q2MAX 
      common/flags/itype(49),nflag(49),inorm(49),itgt(49),icorr(49)
     2     ,itmc(49),iht(49)
      COMMON/CONSTANTS/PI,PI2
      common/q2scale/q2
      common/sacot/ihq,a0u,a0d,pzcoef,s2tw
*     Pion mass in GeV
      parameter(xmpi=0d0,xmpi2=xmpi**2,xmn=0.939d0,xmn2=xmn**2)
*      parameter(xmpi=0.139d0,xmpi2=xmpi**2,xmn=0.939d0,xmn2=xmn**2)

C  CALCULATES DEEP INELASTIC STRUCTURE FUNCTIONS
C  FETCH IS USED TO UNFOLD THE PROPER QUARK DISTRIBUTIONS
C  FROM THOSE WHICH ARE EVOLVED BY INTQCD.
C  FOR THE HIGHER ORDER CALCULATION THE CONVOLUTION WITH
C  THE COEFFICIENT FUNCTIONS IS HANDLED HERE. THE CONVENTIONS
C  FOR THE GLUON COEFFICIENT FUNCTION ARE EXPLAINED BELOW.
c
c    istruc: 0=FL  1=F1  2=F2  3=F3
c    ibeam : 1=e  2=nu  3=nubar 4=electron NC with Z included
c            5 = e for reduced cross section for charm
c            6 = e interference term for PVDIS asymmetry
c
c
c  [A.Accardi] Updated 5/19/08 to include TMC in collinear factorization 
c  and approximate Georgi-Politzer TMC in a unified way.
c
c  [A.Accardi] Updated 26 feb 2009 to compute FL
c
c  updated 12/05 to make the package more streamlined and flexible and 
c  to allow the calculation of the nutev cross section data
c

*    *** input vars, and derived quantities
*    ... NOTE: x here can be x=xB or x=xi_Nachtmann depending on 
*    ... the calling routine 'strfns'
      X=V(1)
      Q2=V(2)
      S0=DLOG(Q02/XC(1)**2)
      S=DLOG(DLOG(Q2/XC(1)**2)/S0)
      xmc2=xmc**2
      xmb2=xmb**2
*    *** tree level str.fn. 
      CALL FETCH(ibeam,itarg,istruc,HT,NDRV,X,S,FTEMP)
      if (istruc.eq.0) then
*       ... FL at LO
         theory=0d0
      else
*       ... F1,F2,F3 at LO
         THEORY=FTEMP(1)
      end if
      IF(IORD.EQ.0)then
         if(istruc.eq.1)theory=theory/2./x
         return
      endif

*    *** 1 loop contribution
c
c  flags for 2xF1, F2, xF3 or FL
c
*    ... flags for F2
      ff1=0.
      ff2=1.
      ff3=0.
      if(istruc.eq.0) then
*       ... FL defined such that FL = F2 - 2xF1 (standard definition)
         ff1 = -1.
         ff2 =  0.
      else if(istruc.eq.1)then
*       ... 2xF1
         ff1 = 1.
      else if(istruc.eq.3)then
*       ... xF3
         ff3=1.
      endif
      CF=4./3.
      AL1=DLOG(1.-X)
      al=alpha_s(iord+1,q2,xc(1),neff)/(4.*pi)
c
c
c     (26 Mar 10) integration limits modified to respect pion threshold
c     rather than xB<1, which is obtained setting the pion mass xmpi=0.
c     The hadron mass parameters are defined in the subroutine 
c     definitions section.
c     NOTE: keep xmpi=0 for the moment, unless for testing purposes.
c           Inclusion of non-zero hadron masses should be reconsidered
c           carefully, to include kaons, heavy flavors, heavy quarks and so on
c
      ymax = 1. / (1.+xmpi2/Q2)  
      ymin = v(4) / (1.-v(4)*2*xmpi*xmn/Q2)
      if (ymin.ge.ymax) then
         THEORY=0d0
      else         
         ! flavor starts at 4 at Q0 and switches to 5 at the b threshold
         ! flavor=3 or 6 are never encountered in the current kinematic range
         flavor=neff
         ! for electron beams, uses fac=sum over quark charges squared.
         ! otherwise fac=flavor
c
c Factor should be 2*(number of active doublets) for neutrinos or antineutrinos
c W(+)g-->u dbar, c sbar, but not t bbar would give FAC=4
c We would never get to FAC=5 (or 3 for that matter). Below charm threshold, 
c FAC=2, above t threshold FAC=6
c
c         FAC=FLAVOR
         FAC=4.
         IF(ibeam.eq.1)then
            fac=10./9.+(flavor-4.)/9.  ! for gamma str.fns
         else if(ibeam.eq.4.or.ibeam.eq.6)then
            fac=1.
         else if(ibeam.eq.5)then
            fac=4./9.  ! charm reduced cross section
         endif
         FX=THEORY
         THEORY=FX+FX*AL*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
         DO 23 I=1,NTERMS
            Y=0.5*(ymax-ymin)*XI(I)+0.5*(ymax+ymin)
            XY=X/Y
            AL1=DLOG(1.-Y)
            CALL FETCH(ibeam,itarg,istruc,HT,NDRV,XY,S,FTEMP)
            C22=CF*(ff2*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*DLOG(Y)
     2           -2.*(1.+Y)*AL1)-ff3*2.*(1.+Y)-ff1*4.*y)
            C23=CF*ff2*(-3.+4.*AL1)/(1.-Y)
c            CG2=2.*FAC*(ff2*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)
c     2         *DLOG(1./Y-1.)) - ff1*4.*Y*(1.-Y))

c
c jfo S-ACOT scheme implemented for charged lepton neutral current case
c jfo 5/29/13
c jfo Controlled by flag IHQ = 0 (massless) = 1 (S-ACOT)
c
            TCG2=2.*(ff2*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)
     2         *DLOG(1./Y-1.)) - ff1*4.*Y*(1.-Y))
            CG2=FAC*TCG2
            if(IHQ.eq.0.or.ibeam.eq.2.or.ibeam.eq.3)then
               CG2=FAC*TCG2
            else 
               if(ibeam.eq.1)then
                  CG2=6./9.*TCG2+4./9.*CG2M(xmc2,y,Q2,ff1,ff2)
     2               +1./9.*CG2M(xmb2,y,Q2,ff1,ff2)
               else if(ibeam.eq.4.or.ibeam.eq.6)then
                  CG2=(a0u+2.*a0d)*TCG2+a0u*CG2M(xmc2,y,Q2,ff1,ff2)
     2               +a0d*CG2M(xmb2,y,Q2,ff1,ff2)
               else if(ibeam.eq.5)then
                  CG2=FAC*CG2M(xmc2,y,Q2,ff1,ff2)
               endif
            endif
C
C     THE ABOVE GLUON COEFFICIENT FUNCTION CORRESPONDS TO THE
C     CONVENTIONS OF FLORATOS,HERROD,WADA,ETC. THE FOLLOWING
C     EXPRESSION CORRESPONDS TO THE CONVENTION OF ALTARELLI,
C     ELLIS,AND MARTINELLI.
C     CG2=2.*FAC*(6.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*DLOG(1./Y-1.))
C
            THEORY = THEORY + .5*(ymax-ymin)*WI(I)*AL*(C22*FTEMP(1)
     &           +C23*(FTEMP(1)-FX))
            if(istruc.lt.3)THEORY=THEORY
     &           +.5*(ymax-ymin)*WI(I)*AL*CG2*FTEMP(2)
 23      end do
      end if
      if(istruc.eq.1)theory=theory/2./x
c      if(istruc.eq.3)print*,x,q2,theory
      RETURN
      END 
c
      function CG2M(xm2,x,q2,ff1,ff2)
      implicit real*8 (a-h,o-z)
      sh=q2*(1./x-1.)
      if(sh.le.4.*xm2)then
         cg2m=0.d0
c         return
      else
         del=sqrt(sh*(sh-4.*xm2))
c      endif
         xm=sqrt(xm2)
         al=2.*log((sqrt(sh)+sqrt(sh-4.*xm2))/(2.*xm))
         xmbar=x*(1.+xm2/q2)
         qps=q2+sh
         g1=al*(q2**2+sh**2)/qps**2-(sh-q2)**2*del/sh/qps**2
     2     +al*4.*xm**2*(sh-2.*xm2)/qps**2-4.*xm2*del/qps**2
         g2=-al*8.*xm2*q2/qps**2+4.*q2*del/qps**2
         cg2m=2.*(ff2*g1+(ff2-ff1)*g2)
c     &       -2.*log(q2/xm2)*(xmbar**2+(1.-xmbar)**2)
      endif
      sub=0.d0
      if(q2.gt.xm2)then
c        sub=log(q2/xm2)*(xmbar**2+(1.d0-xmbar)**2)
        sub=ff2*log(q2/xm2)*(x**2+(1.d0-x)**2)
      endif
      cg2m=cg2m-2.*sub 
      return
      end
c
      SUBROUTINE FETCH(ibeam,itarg,istruc,HT,NDRV,X,S,FTEMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FTEMP(2),HT(0:7)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/GRPTHY/FLAVOR
      common/target/az,an
      common/q2scale/q2
      common/sacot/ihq,a0u,a0d,pzcoef,s2tw

*    *** weights for isospin and HT corrections
*    ... NOTE: ist = 0,1  for FL(p,n)
*              ist = 2,3  for F1(p,n)
*              ist = 4,5  for F2(p,n)
*              ist = 6,7  for F3(p,n)
      ist = istruc*2
      aaz = az * HT(ist)
      aan = an * HT(ist+1)

C
C  UNFOLDS QUARK AND GLUON DISTRIBUTIONS FOR USE IN THEORY2.
C
      call fsupdf(ndrv,x,s,u,d,ub,db,sb,cb,bb,glue)
c
c  s=sbar, c=cbar, b=bbar
c
      sq=sb
      cq=cb
      bq=bb
c
c  correct for neutron/proton ratio unless it is a proton target
c
      if(itarg.eq.2)then
         tu=u
         tub=ub
         td=d
         tdb=db
         u=aaz*tu+aan*td
         ub=aaz*tub+aan*tdb
         d=aaz*td+aan*tu
         db=aaz*tdb+aan*tub
      endif
c
c  change sign of antiquarks if for xf3
c
      iqbar=1
      if(istruc.eq.3) iqbar=-1
c
c  set up q and qbar coefficients
c
      if(ibeam.eq.1.)then
         fu=4./9.
         fub=4./9.
         fd=1./9.
         fdb=1./9.
         ftemp(1)=fu*(u+cq)+fub*(ub+cb)+fd*(d+sq+bq)+fdb*(db+sb+bb)
      else if(ibeam.eq.5)then
         fu=4./9.
         fub=4./9.
         ftemp(1)=fu*cq+fub*cb
      else if(ibeam.eq.2)then
         fu=0.
         fub=2.*iqbar
         fd=2.
         fdb=0.
         ftemp(1)=fu*(u+cq)+fub*(ub+cb)+fd*(d+sq)+fdb*(db+sb)
      else if(ibeam.eq.3)then
         fu=2.
         fub=0.
         fd=0.
         fdb=2.*iqbar
         ftemp(1)=fu*(u+cq)+fub*(ub+cb)+fd*(d+sq)+fdb*(db+sb)
      else if(ibeam.eq.4)then
         s2tw=0.2315
         den=4.*s2tw*(1-s2tw)
         xmz=91.187
         pz=q2/(q2+xmz**2)/den
         pzcoef=pz
         ve=-0.5+2.*s2tw    
         ae=-0.5
         vu=0.5-4./3.*s2tw
         au=0.5
         vd=-0.5+2./3.*s2tw
         ad=-0.5

         if(istruc.eq.2.or.istruc.eq.0.or.istruc.eq.1)then
c
c  Note: this includes the Z exchange, gamma* exchange and interference terms
c
            a0u=4./9.-2.*(2./3.)*vu*ve*pz+(ve**2+ae**2)
     &           *(vu**2+au**2)*pz**2
            a0d=1./9.-2.*(-1./3.)*vd*ve*pz+(ve**2+ae**2)
     &           *(vd**2+ad**2)*pz**2
            ftemp(1)=a0u*(u+cq+ub+cb)+a0d*(d+sq+db+sb+bq+bb)
c            ftemp(1)=a0u*(cq+cb)
         else if(istruc.eq.3)then
            b0u=-2.*2./3.*au*ae*pz+4.*au*vu*ve*ae*pz**2
            b0d=2.*1./3.*ad*ae*pz+4.*ad*vd*ve*ae*pz**2
            ftemp(1)=b0u*(u+cq-ub-cb)+b0d*(d+sq+bq-db-sb-bb)
c            ftemp(1)=0.d0   
         endif
      else if(ibeam.eq.6)then
         s2tw=0.2315
c         s2tw=.2223
         den=4.*s2tw*(1-s2tw)
         xmz=91.187
         pz=q2/(q2+xmz**2)/den
         pzcoef=pz
c         ve=-0.5+2.*s2tw    
c         ae=-0.5
c         vu=0.5-4./3.*s2tw
c         au=0.5
c         vd=-0.5+2./3.*s2tw
c         ad=-0.5
          c1u=-.5+4./3.*s2tw
          c1d=.5-2./3.*s2tw
          c2u=-.5+2.*s2tw
          c2d=.5-2.*s2tw

c
c  Note: this includes only the gamma Z intereference terms
c  Note: remove pz and include as an overall factor if needed (PVDIS)
c
         if(istruc.lt.3)then
            a0u=c1u*2./3
            a0d=-c1d/3.
            ftemp(1)=a0u*(u+cq+ub+cb)+a0d*(d+sq+db+sb+bq+bb)
c            ftemp(1)=a0u*(cq+cb)
         else if(istruc.eq.3)then
            b0u=c2u*2./3.
            b0d=-1./3.*c2d
            ftemp(1)=b0u*(u+cq-ub-cb)+b0d*(d+sq+bq-db-sb-bb)
c            ftemp(1)=0.d0   
         endif

      endif
      ftemp(2)=glue*(aaz+aan)
      if(ihq.eq.0.and.(ibeam.eq.4.or.ibeam.eq.6))then
         ftemp(2)=(2.*a0u+2.*a0d+(flavor-4.)*a0d) * glue*(aaz+aan)
      else if(ihq.eq.1.and.(ibeam.eq.4.or.ibeam.eq.6))then

c
c  Note: The a0u and a0d terms are handled in the heavy quark routines
c
         ftemp(2)=glue*(aaz+aan)
      endif
      return
      end
