      FUNCTION THEORY(MODE,NPT,V,XC)
      IMPLICIT REAL*8 (A-H,O-Z)

      integer nsf,ntg

      
      character*200,fastnlo_table,fln
      REAL*4 PAR(100),V4(4),ANSWER
      real*8 nuke_cteq
      DIMENSION V(4),XC(100),GF(11,60,2160)
     2,VP(4),FTEMP(2),NFLAG(49),INORM(49),itype(49),ITGT(49),icorr(49)
     3,itmc(49),iht(49), tv(4),tmp(1008), sumd0(36)
      dimension xsect1(1188,3),xsect2(3240,3),xsect3(3960,3),
     2xsect4(2592,3)
      dimension temp1(33,3),temp2(90,3),temp3(110,3),temp4(77,3)
      COMMON/GINT/GF
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/GFUNC/CALC(11,60,60)
      COMMON/Q2STUF/ Q02,Q2MAX
      COMMON/GAUSS4/XI(4),WI(4),NTERMS,XX(5)
      COMMON/GRPTHY/FLAVOR
      COMMON/FLAGS/ITYPE,NFLAG,INORM,ITGT,icorr,itmc,iht
      COMMON/MINIMB/UNCRT(100),ERMIN(50),PWATE(100),
     *     IFREP(50),NPARAM,NVAR,IREJ
      COMMON/CONSTANTS/PI,PI2
      common/target/az,an
      common/e866xf/ixfx
      common/jetpdfs/ndrv,nflav,als
      common/jetpdfs2/amu0,alambda5
      common/stepsize/delta
      common/threshold/sb
      common/d0runIINP/anp(110)
      common/numjet/n201,n202,n203,n204
      common/nd0z/n141
      common/sacot/ihq,a0u,a0d,pzcoef,s2tw
 
c-jfo
c-jfo  data xmc,xmb,ncb/1.3d0,4.5d0,12/
c-jfo  these values are now set in ALTPFIT10
c-jfo
*     Nucleon mass [GeV]
      double precision mn,mn2
      data mN,mN2/0.939,0.882/      
*     Initialization flags
      data itime/0/
      data ijet/0/
      data id0z/0/
      data ifastnlo/0/
      data iztest/0/
      save ijet,itime,ifastnlo,id0z,sumd0,npt1
C  THIS ROUTINE CONTROLS THE CALCULATION OF THE VARIOUS
C  OBSERVABLES WHICH ARE TO BE FITTED. IT MAKES USE OF A NUMBER OF
C  ASSOCIATED ROUTINES AND CAN EASILY BE EXPANDED FOR THE CALCULATION 
C  OF ADDITIONAL QUANTITIES.
c
c  [A.Accardi] Updated 5/19/08 to include TMC in collinear factorization 
c  and approximate Georgi-Politzer TMC in a unified way
c
c  [A.Accardi] 9/9/08 inclusion of Higher-Twist corrections
c
c  Approximate Georgi-Politzer target mass corrections added 5/5/08
c
c  modified 12/05 to calculate nutev cross section data
c  corrected DIS bug for calculating f_2^d/f_2^p
c  corrected bug in gluon coefficient for DIS with nf=5
c  generalized 'theory2' and 'fetch' to make the routines more streamlined 
c  and versatile
C
C  THEORY2 CALCULATES THE DIS OBSERVABLES
C  DYANTH CALCULATES THE DRELL-YAN AND J/PSI OBSERVABLES
C
C  DECIDE WHETHER OR NOT TO CALL INTQCD 
C
C  MODE=1....NO DERIVATIVES NEEDED
C      =0....CALCULATE DERIVATIVES
C
C  NPT=1....FILL GF 
C     >1....PROCEED 
C
C  MODIFIED 11/87 SO THAT EVOLUTION IS NOT PERFORMED WHEN 
C  NORMALIZATION PARAMETERS ARE VARIED...SAVES TIME.
C
*  NOTES:
*
*  (26 Feb 08) HT only for F2  (F1==F3==F2 for now)
*              GP only for F2, F3
*              CF, NV only for F1, F2, F3
*              Nuke smearing only for F2 (There is none for F1, F3)
*
c
c  normalization for observables with zero normalization error
c
      xc(100)=1.d0
      IF(NPT-1)1,1,2 ! NPT = 0,1 --> fillgf (perform DGLAP); >1 proceeds
      !-------- DGLAP evolution performed -----
    1 IF(MODE-1) 3,4,3   
    4 NDRV=0
      NOLD=1
      GO TO 7
    3 CONTINUE
      NDRV=NDRV+1
      IF(IFREP(NDRV)-35)7,7,5  ! If 'PDF' parameters (NDRV<=35) is being 
                               ! changed performs QCD evolution, and offshell
                               ! initialization, else proceeds
    7 CONTINUE
      S0=DLOG(Q02/XC(1)**2)
      SMAX=DLOG(Q2MAX/XC(1)**2)
      SMAX=DLOG(SMAX/S0)
c
c  Have to calculate delta before it can be used
c  Moved calcualtion from ALTPAR to here
c
      sb=dlog(xmb**2/xc(1)**2)
      sb=dlog(sb/s0)
      delta=sb/ncb
      NMAX=SMAX/DELTA+3
      if(nmax.gt.60) nmax=60
c      if(itime.eq.0)then
c         print*,s0,sb,delta,ncb,smax,nmax
c         itime=1
c      endif

C
C  FILL ARRAY GF
C
      CALL INTQCD(XC)
      DO J=1,60
         DO K=1,60
            KK=K+NDRV*60
            DO IL=1,11
               GF(IL,J,KK)=CALC(IL,J,K)
            end do
         end do
      end do
      call setoffshell(ndrv,xc) 
      GO TO 5 ! Skips to theory calculations section
C
C  Proceeds when there is no need to re-evolve the PDFs
    2 IF(MODE-1) 8,9,8
    9 NDRV=0 ! when MODE=1 sets ndrv=0 to use central PDFs
      GO TO 5
    8 IF(NOLD-NPT) 10,11,10 ! increments ndrv until needed
   11 NDRV=NDRV+1
      GO TO 5
   10 NOLD=NPT
      NDRV=1

C
C  Theory calculations start here
C

 5    CONTINUE
      index=V(4)
      it=itype(index)
      nfl=nflag(index)
      itg=itgt(index)
      itm=itmc(index)
      ihtfl=iht(index)
*    ... old (pre-080519) vp assignments
      vp(1)=v(1)
      vp(2)=v(2)
      vp(3)=v(3)
      vp(4)=v(4)
      r=1
*      vp(1) = 0.0699596366
*      vp(2) = 100.
*      vp(3) = vp(1)
*      itm=0
*      ihtfl=0
*      az = 1.
*      an = 0.
*      call strfn(1,2,itm,ihtfl,NDRV,VP,XC,F222)
*      print*,'*strfn -- z,Q2, F2 =',vp(1),vp(2),vp(3),f222
*      call theory2(1,2,2,ndrv,vp,xc,ans)
*      print*,'*thry2 -- z,Q2, F2 =',vp(1),vp(2),vp(3),ans
*      stop

****  Jeff's PVDIS scheme
*$$$      IF(IT.LT.80.and.it.ne.60.and.it.ne.61)then
      IF(IT.LT.60.or.IT.eq.78.or.it.eq.79)then
*    *** DIS structure functions
         nsf = nfl/10
         ntg = nfl-nfl/10*10
         
*       ... mapping of flags for the 'strfn' routine
         if (nsf.eq.0) then
            if(ntg.eq.4.or.ntg.eq.8)then
               nsf = 3
            else 
               nsf = 2             ! F2 structure function
            endif
         else if (nsf.eq.2) then
            nsf = 0             ! FL structure function
         end if                 ! all others are already OK

         IF(ntg.eq.9)then
C       *** D/P RATIO CALCULATED HERE for charged lepton DIS
C          ... first compute F_2(D)
	    call deuteronSF(1,nsf,itg,itm,ihtfl,ndrv,vp,xc,F2deut)
C          ... CALCULATE PROTON SF HERE
            az=1.
            an=0.
            CALL strfn(1,nsf,itm,ihtfl,NDRV,VP,XC,F2P)
            theory=f2DEUT/f2p   ! modified to F2D/F2P  10/21/05
            IN=Inorm(index)
            theory=theory/xc(in) ! NOTE theory is normalized to data
c
c  Normalization parameter added to F2D/F2P
c
            RETURN
         else if(ntg.eq.6)then
c
c        F2n/F2D for Bonus data calculated here
c        Have to multiply our F2deut by 2 to martch that of Bonus
c
C          ... first compute F_2(D)
	    call deuteronSF(1,nsf,itg,itm,ihtfl,ndrv,vp,xc,F2deut)
            F2deut=F2deut*2.
C          ... CALCULATE NEUTRON SF HERE
            az=0.
            an=1.
            CALL strfn(1,nsf,itm,ihtfl,NDRV,VP,XC,F2N)
            theory=f2n/f2deut   
            IN=Inorm(index)
            theory=theory/xc(in)
            RETURN            
         else if(ntg.eq.7)then
c       *** D/(p+n) ratio F#D/F#N for EMC-like plots calculated here
C          ... first compute F_2(D) with smearing
	    call deuteronSF(1,nsf,itg,itm,ihtfl,ndrv,vp,xc,F2deut)
C          ... then F_2(p+n), i.e. F_2(D) with only isospin correction
C          ... and call it f2n for lazyness
	    call deuteronSF(1,nsf,0,itm,ihtfl,ndrv,vp,xc,F2n) 
C          ... then take the ratio
            theory=f2deut/f2n ! this is D/(p+n) ratio
c          ... normalization correction
            IN=Inorm(index)
            theory=theory/xc(in)
            RETURN            
         endif
C
C  DIS STRUCTURE FUNCTIONS CALCULATED HERE
C
c  the logic here must be updated as new data sets are added 
c
c  note: s=sbar, c=cbar, and b=bar are assumed for this version
c
c        higher-twist only for F2 (09 sep 2008)
c

         if(ntg.eq.0)then
*          ... FF(e+n)
            az=0.
            an=1.
            call strfn(1,nsf,itm,ihtfl,ndrv,vp,xc,ans)
         else if(ntg.eq.1)then
*          ... FF(e+p)
            az=1.
            an=0.
c SLAC_R
c [By Jeff, ~ summer 2013 ????]
            if(it.ge.21.and.it.lt.30)then
c
c  GP Target Mass Corrections not yet implemented for FL
c  For now, turn off HT (same for F2, FL) and TMC
c
               if(it.eq.21)then
                  call strfn(1,1,itm,ihtfl,ndrv,vp,xc,ans1)
                  call strfn(1,2,itm,ihtfl,ndrv,vp,xc,ans2)
                  rho2=1.+4.*0.8836*vp(1)**2/vp(2)
c                  ans2=ans2*rho2
c                  ans=ans0/(ans2-ans0)             
                  ans=-1. +ans2/(2.*vp(1)*ans1)*rho2
                  theory=ans
               else if(it.eq.22.or.it.eq.23.or.it.eq.24)then
                  call strfn(1,0,itm,ihtfl,ndrv,vp,xc,ans)
                  theory=ans
               endif
               return
            else
               call strfn(1,nsf,itm,ihtfl,ndrv,vp,xc,ans)       
            endif
         else if(ntg.eq.2)then
*          ... F2(e+D)
	    call deuteronSF(1,nsf,itg,itm,ihtfl,ndrv,vp,xc,ans)
         else if(ntg.eq.5)then
*          ... F2(nu) == [F2(nu)+F2(nubar)]/2
*          ... Only for D or Fe targets, only for F2
*          ... --- legacy code by Jeff Owens---
            if (nsf.ne.2) then
               print*, 'ERROR(theory): ntg=5 works only for F2'
               print*,nfl,nsf,ntg,it
               stop 
            end if
            if (itg.ge.0) then 
               iitg = itg
            else
               iitg = 0
            end if
	    call deuteronSF(2,nsf,iitg,itm,ihtfl,ndrv,vp,xc,ans1)
	    call deuteronSF(3,nsf,iitg,itm,ihtfl,ndrv,vp,xc,ans2)
            if(itg.eq.-3) then
*             ... Kulagin-Petti nuclear corrections for Fe targets
               ans1=ans1 * nuke_cteq(v(1),v(2),2,0,2,3,2)
               ans2=ans2 * nuke_cteq(v(1),v(2),2,0,-2,3,2)
            endif
            ans=(ans1+ans2)/2.
            if(itg.eq.-1)then
*             ... EMC/NMC parametrization of Fe corrections
               cor = emcnmc(1d0,v(1),v(2))
               ans=ans/cor
            else if(itg.eq.-2)then
*             ... EMC/NMC + density model correction
*             ... for backward compatibility
               cor  = emcnmc(1d0,v(1),v(2))
               dmc  = gomez(1d0,v(1))
               ans=ans/(cor*dmc)
            endif
         else if(ntg.eq.8)then
*          ... F3(nu) == [F3(nu)+F3(nubar)]/2 
*          ... Only for D or Fe targets, only F3
*          ... No deuterium corrections yet beside isospin averaging
*          ... --- legacy code ---
            if (nsf.ne.3) then
               print*, 'ERROR(theory): ntg=8 works only for F3'
               stop 
            end if
            az=0.5
            an=0.5
            call strfn(2,3,itm,ihtfl,ndrv,vp,xc,ans1)
            call strfn(3,3,itm,ihtfl,ndrv,vp,xc,ans2)
            if(itg.eq.-3)then
*             ... Kulagin-Petti nuclear corrections for Fe targets
               ans1=ans1 *nuke_cteq(v(1),v(2),3,0,2,3,2)
               ans2=ans2 *nuke_cteq(v(1),v(2),3,0,-2,3,2)
            end if
            ans=(ans1+ans2)/2.
            if(itg.eq.-1)then
*             ... EMC/NMC parametrization of Fe corrections
               cor = emcnmc(1d0,v(1),v(2))
               ans=ans/cor
            else if(itg.eq.-2)then
*             ... EMC/NMC + density model correction
*             ... for backward compatibility
               cor  = emcnmc(1d0,v(1),v(2))
               dmc  = gomez(1d0,v(1))
               ans=ans/(cor*dmc)
            endif 
c
c  WA-25 neutrino structure functions
c
         else if(ntg.eq.3)then
            if(it.eq.42)then
               az=0.5
               an=0.5
               call strfn(2,2,itm,ihtfl,ndrv,vp,xc,ans)
            else if(it.eq.44)then
               az=1.
               an=0.
               call strfn(2,2,itm,ihtfl,ndrv,vp,xc,ans)
            else if(it.eq.46)then
               az=0.
               an=1.
               call strfn(2,2,itm,ihtfl,ndrv,vp,xc,ans)
            endif
         else if(ntg.eq.4)then
            print*,'* xF3'
            if(it.eq.43)then
               az=0.5
               an=0.5
               call strfn(2,3,itm,ihtfl,ndrv,vp,xc,ans1)
               call strfn(3,3,itm,ihtfl,ndrv,vp,xc,ans2)
               ans=(ans1+ans2)/2.
               ! print*,'* ',vp(1),vp(2),ans
            else if(it.eq.45)then
               az=1.
               an=0.
               call strfn(2,3,itm,ihtfl,ndrv,vp,xc,ans1)
               az=0.
               an=1.
               call strfn(3,3,itm,ihtfl,ndrv,vp,xc,ans2)
               ans=(ans1+ans2)/2.
            else if(it.eq.47)then
               az=0.
               an=1.
               call strfn(2,3,itm,ihtfl,ndrv,vp,xc,ans1)
               az=1.
               an=0.
               call strfn(3,3,itm,ihtfl,ndrv,vp,xc,ans2)
               ans=(ans1+ans2)/2.
            endif
         endif


****  the commented out section is Jeff's PVDIS code from Dec 2017 *****
         
c$$$      else if (it.eq.60.or.it.eq.61)then
c$$$c
c$$$c  Parity violating asymmetry
c$$$c  It=60 is for PVDIS (deuterium)
c$$$c  It=61 is for the EIC simulated data (proton)
c$$$c
c$$$         if(it.eq.60)then !<--- keep for backward compatibility with Jeff
c$$$            !az=0.5   ! [AA] NOTE: deuteronSF already sets az=an=0.5
c$$$            !an=0.5
c$$$            call deuteronSF(6,1,itg,itm,ihtfl,ndrv,vp,xc,F1int)
c$$$            call deuteronSF(6,3,itg,itm,ihtfl,ndrv,vp,xc,F3int)
c$$$            call deuteronSF(1,1,itg,itm,ihtfl,ndrv,vp,xc,F1)
c$$$         else if(it.eq.61) then ! proton APV_red
c$$$            az=1.
c$$$            an=0.
c$$$            call strfn(6,1,itm,iht,ndrv,vp,xc,F1int)
c$$$            call strfn(6,3,itm,iht,ndrv,vp,xc,F3int)
c$$$            call strfn(1,1,itm,iht,ndrv,vp,xc,F1)
c$$$         else
c$$$            print*,'Wrong value of itype for A_pv'
c$$$            stop
c$$$      endif
c$$$c
c$$$c program returns xF3
c$$$c
c$$$         x=vp(1)
c$$$         y=vp(3)
c$$$         F3int=F3int/x
c$$$         yp=1+(1-y)**2
c$$$         ym=1-(1-y)**2
c$$$         yf=ym/2./yp
c$$$         apv=-(F1int+yf*F3int)/F1
c$$$         if(it.eq.60)apv=apv*pzcoef
c$$$         ans=apv
c$$$         print*,x,vp(2),y,F1int/F1,yf,F3int/F1

***** End of Jeff's code  *********************


***** Beginning of Alberto's PVDIS code, Feb 2017  *************

      else if (it.ge.60.and.it.le.77)then
         nsf = nfl/10
         ntg = nfl-nfl/10*10
         
*    *** Parity violating asymmetry
*    ... it = 60-79
*    ... ntg = 0,1,2  [n,p,d targets]
*    ... nsf = 4  "reduced" A_PV     
*    ...      5  full A_PV
*     NOTE: this scheme decouples itgt = {code for a given data set}
*           from nfl = (nsf)(ntg) = observable & target
*           as opposed to Jeff's scheme, in which itgt controls everything
*     Translation Jeff --> Alberto
*         it(Jeff)=60 ---> it(AA)=whatever, ntg=2, nsf=5
*         it(Jeff)=61 ---> it(AA)=whatever, ntg=1, nsf=4
*     [AA Feb 2017]
         
         if(ntg.eq.0.or.ntg.eq.1) then ! neutron/proton str.fns. for APV
            az=real(ntg)
            an=1.-real(ntg)
            call strfn(6,1,itm,iht,ndrv,vp,xc,F1int)
            call strfn(6,3,itm,iht,ndrv,vp,xc,F3int)
            call strfn(1,1,itm,iht,ndrv,vp,xc,F1)
         else if (ntg.eq.2) then ! Deuteron str. fns. for APV
            call deuteronSF(6,1,itg,itm,ihtfl,ndrv,vp,xc,F1int)
            call deuteronSF(6,3,itg,itm,ihtfl,ndrv,vp,xc,F3int)
            call deuteronSF(1,1,itg,itm,ihtfl,ndrv,vp,xc,F1)
         else
            print*,'[theory15] Wrong NTG value '//
     &           '[0,1,2 = n,p,d targets] =', NTG
            stop
         endif
c
         x=vp(1)
         y=vp(3)
         F3int=F3int/x  ! This is b/c 'strfn' routines return x*F3
         yp=1+(1-y)**2
         ym=1-(1-y)**2
         yf=ym/2./yp
         apv=-(F1int+yf*F3int)/F1
         !print*, '* x,Q2,y =', vp(1),vp(2),vp(3)
         !print*, '* yp,ym =',yp,ym
         !print*, '* F1int,F3int,F1 =',F1int,F3int,F1
         if (nsf.eq.4) then 
            ans=apv             ! "reduced" APV
         else if (nsf.eq.5) then
            ans=apv*pzcoef      ! full APV
         else
            print*,'[theory15] Wrong NSF value for PVDIS '//
     &           '[4,5 = reduced, full APV] =', NSF
            stop
         end if

         !print*,'* ',x,vp(2),y,F1int/F1,yf,F3int/F1

C     ****************
C     *  HERA DIS cross sections
C     *  it=80-84 is for HERA I combination
C     *  it=93-99 is for HERA I+II combination
C     ****************
      else if ( it.eq.80.or.it.eq.81.or.it.eq.84
     &        .or.(it.ge.93.and.it.le.97) ) then
c     it=80 is HERA I electron neutral current
c     it=81 is HERA I positron neutral current
c     it=93 is HERA I+II electron neutral current
c     it=94-97 are HERA I+II four positron neutral current sets
         az = 1.                ! proton target
         an = 0.
         x=vp(1)
         Q2=vp(2)
         y=vp(3)
         is=1  ! NC electrons
         if(it.eq.81.or.(it.ge.93.and.it.le.97))is=-1  ! NC positrons
         yp=1.+(1.-y)**2+2*x*x*y*y*mN2/Q2                  
         ym=1.-(1.-y)**2
         
         if(it.eq.80.or.it.eq.81.or.(it.ge.93.and.it.le.96))then ! NC sigma
            call strfn(4,2,itm,ihtfl,ndrv,vp,xc,ansf2)
            call strfn(4,3,itm,ihtfl,ndrv,vp,xc,ansxf3)
            call strfn(4,0,itm,ihtfl,ndrv,vp,xc,ans0)

            !write(6,85),'* ',x,Q2,y,ansf2,ansxf3,ans0
 85         FORMAT(A,e9.4,F10.3,F7.3,2X,3E14.6)
            
            ans=ansf2+is*ym/yp*ansxf3-y**2/yp*ans0
         else if(it.eq.84)then
c
c           ibeam=5 is for charm production
c           This is the HERA-I combined reduced cross section for charm
c
            call strfn(5,2,itm,ihtfl,ndrv,vp,xc,ansf2)
            call strfn(5,0,itm,ihtfl,ndrv,vp,xc,ans0)
            ans=ansf2-y**2/yp*ans0
c            print*,'* ',x,Q2,ans,ansf2,ans0
         endif
         
      else if(it.eq.82.or.it.eq.83.or.it.eq.98.or.it.eq.99)then
         az = 1.                ! proton target
         an = 0.
         x=vp(1)
         Q2=vp(2)
         y=vp(3)
         is=1
         if(it.eq.83.or.it.eq.99)is=-1
         yp=1.+(1.-y)**2+2*x*x*y*y*mN2/Q2
         ym=1.-(1.-y)**2
c
c  it=82, 98 is HERA I (HERA I+II) electron charged current
c  which is like antineutrino scattering -->ibeam=3 
c  it=83,99 HERA I (HERA I+II) is positron charged current
c  which is like neutrino scattering -->ibeam=2
c
         ibeam=3    !  CC electrons  
         if(it.eq.83.or.it.eq.99)ibeam=2   ! CC positrons
         call strfn(ibeam,2,itm,ihtfl,ndrv,vp,xc,ansf2)
         call strfn(ibeam,3,itm,ihtfl,ndrv,vp,xc,ansxf3)
         call strfn(ibeam,0,itm,ihtfl,ndrv,vp,xc,ans0)
c         ans0=ansf2-2.*x*ansf1
c         print*,'* ',x,y,vp(2),ansf2,ansf2,ans0
c         print*,'* ',x,y,vp(2),ansf2,ans0
         ans=0.25*yp*ansf2+is*0.25*ym*ansxf3-0.25*y**2*ans0

c     ************
c     * DIS cross section in fixed target experiments
c     *   it = 85-92
c     *   it = odd --> electron-proton DIS sigma
c     *   it = even --> electron-Deuteron DIS sigma
c     * There are 4 pairs of it values for 4 different DIS experiments
c     ************
      else if(it.eq.85.or.it.eq.87.or.it.eq.89.or.it.eq.91)then
c     electron-proton neutral current cross-sections
         az = 1.                ! proton target
         an = 0.
         x=vp(1)
         Q2=vp(2)
         y=vp(3)
         is=1
         yp=1.+(1.-y)**2+2*x*x*y*y*mN2/Q2
c         yp=1.+(1.-y)**2
         ym=1.-(1.-y)**2
c         call strfn(4,2,itm,ihtfl,ndrv,vp,xc,ansf2)
c         call strfn(4,3,itm,ihtfl,ndrv,vp,xc,ansxf3)
c         call strfn(4,0,itm,ihtfl,ndrv,vp,xc,ans0)
         call strfn(1,2,itm,ihtfl,ndrv,vp,xc,ansf2)
         call strfn(1,0,itm,ihtfl,ndrv,vp,xc,ans0)
c         if (it.eq.91) then ! Hermes
c            ansxf3=0D0
c         endif
c         ans=ansf2+is*ym/yp*ansxf3-y**2/yp*ans0
         ans=ansf2-y**2/yp*ans0
      else if(it.eq.86.or.it.eq.88.or.it.eq.90.or.it.eq.92) then
c     electron-proton neutral current cross-sections
         x=vp(1)
         Q2=vp(2)
         y=vp(3)
         is=1
         yp=1.+(1.-y)**2+2*x*x*y*y*mN2/Q2
c         yp=1.+(1.-y)**2
         ym=1.-(1.-y)**2
c         call DeuteronSF(4,2,itg,itm,ihtfl,ndrv,vp,xc,ansf2)
c         call DeuteronSF(4,3,itg,itm,ihtfl,ndrv,vp,xc,ansxf3)
c         call DeuteronSF(4,0,itg,itm,ihtfl,ndrv,vp,xc,ans0)
         call DeuteronSF(1,2,itg,itm,ihtfl,ndrv,vp,xc,ansf2)
         call DeuteronSF(1,0,itg,itm,ihtfl,ndrv,vp,xc,ans0)
c         if (it.eq.92) then ! Hermes
c            ansxf3=0D0
c         endif
c         ans=ansf2+is*ym/yp*ansxf3-y**2/yp*ans0
         ans=ansf2-y**2/yp*ans0


c     NEUTRINO CROSS SECTIONS intead of structure functions
c        (no TMC, nor HT; Kulagin-Petti corrections only for Fe)
c
c     The nfl=10,11 case was written by Jeff to provide cross section
c     computations for the NuTeV collaboration. It has not yet been 
c     modified to include TMC or HT corrections.
c$$$         else if(nfl.eq.10)then
c$$$            x=v(1)
c$$$            q2=v(2)
c$$$            y=v(3)
c$$$            if(itg.eq.-3)then
c$$$               az=0.5
c$$$               an=0.5
c$$$               call strfn(2,1,itm,ihtfl,ndrv,v,xc,ans1)
c$$$               call strfn(2,2,itm,ihtfl,ndrv,v,xc,ans2)
c$$$               call strfn(2,3,itm,ihtfl,ndrv,v,xc,ans3)
c$$$               ans1=ans1*nuke_cteq(x,q2,1,0,2,3,2)
c$$$               ans2=ans2*nuke_cteq(x,q2,2,0,2,3,2)
c$$$               ans3=ans3*nuke_cteq(x,q2,3,0,2,3,2)
c$$$            else if (itg.eq.-1.or.itg.eq.-2) then
c$$$               az=0.465
c$$$               an=0.535
c$$$               call strfn(2,1,itm,ihtfl,ndrv,v,xc,ans1)
c$$$               call strfn(2,2,itm,ihtfl,ndrv,v,xc,ans2)
c$$$               call strfn(2,3,itm,ihtfl,ndrv,v,xc,ans3)
c$$$            endif
c$$$            ans=(1.-y-(.938*x*y)**2/q2)*ans2+y**2*x*ans1
c$$$     2          +y*(1-y/2.)*ans3
c$$$            ans=ans*1.5816/(1.+q2/80.22**2)**2
c$$$         else if(nfl.eq.11)then
c$$$            x=v(1)
c$$$            q2=v(2)
c$$$            y=v(3)
c$$$            if(itg.eq.-3)then
c$$$               az=0.5
c$$$               an=0.5
c$$$               call strfn(3,1,itm,ihtfl,ndrv,v,xc,ans1)
c$$$               call strfn(3,2,itm,ihtfl,ndrv,v,xc,ans2)
c$$$               call strfn(3,3,itm,ihtfl,ndrv,v,xc,ans3nubar)
c$$$               call strfn(2,3,itm,ihtfl,ndrv,v,xc,ans3nu)
c$$$               ans1=ans1*nuke_cteq(x,q2,1,0,-2,3,2)
c$$$               ans2=ans2*nuke_cteq(x,q2,2,0,-2,3,2)
c$$$               ans3=(ans3nu+ans3nubar)*nuke_cteq(x,q2,3,0,-2,3,2)
c$$$     2                         -ans3nu*nuke_cteq(x,q2,3,0,2,3,2)
c$$$            else if (itg.eq.-1.or.itg.eq.-2) then
c$$$               az=0.465
c$$$               an=0.535
c$$$               call strfn(3,1,itm,ihtfl,ndrv,v,xc,ans1)
c$$$               call strfn(3,2,itm,ihtfl,ndrv,v,xc,ans2)
c$$$               call strfn(3,3,itm,ihtfl,ndrv,v,xc,ans3)
c$$$            endif
c$$$            ans=(1.-y-(.938*x*y)**2/q2)*ans2+y**2*x*ans1
c$$$     2          -y*(1-y/2.)*ans3
c$$$            ans=ans*1.5816/(1.+q2/80.22**2)**2
         


      else if(it.gt.100.and.it.le.200)then
C  AS CURRENTLY SET UP ITYPE < 100 CORRESPONDS TO DEEP
C  INELASTIC SCATTERING WHILE DRELL-YAN AND W PRODUCTION HAVE ITYPE = 100
C  OR GREATER. THIS CONVENTION MUST BE REMEMBERED IF ADDITIONS OR
C  ALTERATIONS ARE MADE.
         az=1.
         an=0. 
         ixfx=0
         IF(nfl.eq.1.or.nfl.eq.2)THEN
            if(it.eq.106)then ! E605 p+Cu DY 
               az=0.456
               an=0.544
            else if(it.eq.108)then ! E886 p+p DY
               az=1.
               an=0.
               ixfx=1
            else if(it.eq.110)then ! E866 p+d DY
               az=0.5
               an=0.5 
               ixfx=1
            else if(it.eq.111)then ! ????
               az=1.
               an=0.
               ixfx=1
            else if(it.eq.140)then
c
c  Z dsig/dy
c
               az=1.
               an=0.
               ixfx=2
            else if(it.eq.141)then
c
c  D0 Z normalized rapidity distribution
c  calculate total of all bins first
c  rescale calculated data
c  total of 28 data points
c
               if(npt.eq.n141)then
                  if(ndrv.ne.0)then
                     if(ifrep(ndrv).gt.35)goto 1410
                  endif
                  sumd0(ndrv+1)=0.
                  az=1.
                  an=0.
                  ixfx=2
                  do j=1,28
                     y=.05+.1*(j-1)
                     tv(1)=v(1)
                     tv(2)=v(2)
                     tv(3)=y
                     tv(4)=v(4)
                     call dyanth(ndrv,tv,xc,tmp(j+28*ndrv))
                     sumd0(ndrv+1)=sumd0(ndrv+1)+tmp(j+28*ndrv)
                  enddo
               endif
 1410          continue
               k=npt-n141+1
               ic=0
               if(ndrv.ne.0)then
                  if(ifrep(ndrv).gt.35)then
                     ic=0
                  else
                     ic=ndrv
                  endif
               endif
               theory=tmp(k+28*ic)*4.991/sumd0(ic+1)
               return
c
c  pp Z rapidity distribution
c
            else if(it.eq.142)then
               az=1.
               an=0.
               ixfx=3
            endif
            if(iztest.eq.0.and.it.eq.140)then
               do jz=1,15
                  y=.1+.2*(jz-1)
                  tv(1)=v(1)
                  tv(2)=v(2)
                  tv(3)=y
                  tv(4)=v(4)
                  call dyanth(ndrv,tv,xc,tstz)
                  !print*,'* ',y,tstz
               enddo
            iztest=1
            endif     
            CALL DYANTH(NDRV,V,XC,THEORY)
            IF(IT.EQ.106.or.it.eq.108.or.it.eq.110.or.it.eq.111) THEN
               THEORY=THEORY*(V(2)/SQRT(V(1)))**3
               in=inorm(index)
               theory=theory/xc(in)
            ENDIF
c
c  Correct ATLAS Z rapidity distribution for the effects of the 
c  lepton pt cut (observable type 142)
c 
            if(it.eq.140)then
               in =inorm(index)
               theory=theory/xc(in)
            else if(it.eq.142)then
               in=inorm(index)
c               theory=theory*(0.776+0.011*v(3))/xc(in)  
               theory=theory*(0.79303-0.044921*v(3)
     2                +0.014762*V(3)**2)/xc(in)  
            endif
            if(it.eq.151.or.it.eq.152)then
               S=V(1)
               Y=V(3)
               RS=SQRT(S)
               ALS=DLOG(DLOG(80.**2/XC(1)**2)/dlog(q02/xc(1)**2))
               CALL WASYM(RS,Y,ALS,NDRV,THEORY,it)
               in=inorm(index)
               theory=theory/xc(in)
             endif   
         ELSE IF(IT.EQ.125)THEN         
            CALL DYANTH(NDRV,V,XC,THEORYP)
            az=0.
            an=1.
            CALL DYANTH(NDRV,V,XC,THEORYN)
            THEORY=(THEORYP-THEORYN)/(THEORYP+THEORYN)
         else if(it.eq.133.or.it.eq.136)then
            ixfx=1

            az=1. ! p+p cross section
            an=0. 
            nfl = nflag(index) 
            nflag(index)=0 ! sets no nuclear correction for p+p 
            call dyanth(ndrv,v,xc,theoryp)
            nflag(index)=nfl ! restores the user selected nuclear correction

            az=0.5 ! p+d cross section
            an=0.5
            call dyanth(ndrv,v,xc,theoryd)
c
c  E866 data are for sig(pd)/2sig(pp) where sig(pd)~sig(pp)+sig(pn)
c  Here we must have az+an=1. so the correct theory is sig(pd)/sig(pp)
c  since the theory calculatio is on a per nucleon basis.
c
            theory=theoryd/theoryp
c            s=v(1)
c            q=v(2)
c            xf=v(3)
c            rt=sqrt(xf**2+4.*q**2/s)
c            x2=(-xf+rt)/2.
c            print*,x2,q,xf,theoryp,theoryd,theory
         ELSE IF(IT.GE.126.and.IT.LE.132.or.it.eq.134.or.it.eq.135
     >   .or.it.eq.153)THEN
            S=V(1)
            Y=V(3)
            RS=SQRT(S)
            ALS=DLOG(DLOG(80.**2/XC(1)**2)/dlog(q02/xc(1)**2))
            CALL WASYM(RS,Y,ALS,NDRV,THEORY,it)
         ENDIF
         return

      else if(it.gt.200.and.it.le.300)then
c
c  Needed for fastNLO
c
         amu0=sqrt(q02)
         alambda5=xc(1)
         if (ifastnlo.eq.0) then
            call GETENV('cteqx',fln)
            call trmstr(fln,ifastnlo)
            fln=fln(1:ifastnlo)//'/theory/'
            ifastnlo = ifastnlo + 8
         end if

c          if(it.eq.201)then
c         if(it.eq.201.or.it.eq.202)then
* JET observables for Run I D0 or CDF
*
c            ymin=v(1)
c            ymax=v(2)
c            pt=v(3)
c
c  hardwire these for now - may pass from input file later on
c
c            rs=1800.
c            smucoef=.5
c
c            q2=(smucoef*pt)**2
c            s0=dlog(q02/xc(1)**2)
c            als=dlog(dlog(q2/xc(1)**2)/s0)
c            call jetdet(rs,pt,smucoef,ans,ymin,ymax,iord)
         if(it.eq.201)then
c
c  CDF Run I cone algorithm as calculated by fastNLO
c
            if(npt.eq.n201)then
c               ijet=1
c               smucoef=1.0
               if(ndrv.ne.0)then
                  if(ifrep(ndrv).gt.35) goto 251
               endif
               smucoef=1.0      
               fastnlo_table = fln(1:ifastnlo)//'/fnt1001rsep.tab'
               call ft1001cc(fastnlo_table,smucoef,smucoef,0,temp1)
               do jn=1,33
                  do jm=1,3
                     xsect1(jn+33*ndrv,jm)=temp1(jn,jm)
                  enddo
               enddo
            endif
 251        continue
            kpt=npt-n201+1
            ic=0
            if(ndrv.ne.0)then
               if(ifrep(ndrv).gt.35)then
                  ic=0
               else 
                  ic=ndrv
               endif
            endif
            if(iord.eq.0)then
               ans=xsect1(kpt+33*ic,1)
            else
               ans=(xsect1(kpt+33*ic,1)+xsect1(kpt+33*ic,2))
            endif
         else if(it.eq.202)then
c
c  D0 Run I cone algorithm as calculated by fastNLO
c
            if(npt.eq.n202)then
c               ijet=1
c               smucoef=1.0       
               if(ndrv.ne.0)then
                  if(ifrep(ndrv).gt.35) goto 252
               endif
               smucoef=1.0       
               fastnlo_table = fln(1:ifastnlo)//'/fnt1002rsep.tab'
               call ft1002cc(fastnlo_table,smucoef,smucoef,0,temp2)
               do jn=1,90
                  do jm=1,3
                     xsect2(jn+90*ndrv,jm)=temp2(jn,jm)
                  enddo
               enddo
            endif
 252        continue
            kpt=npt-n202+1
            ic=0
            if(ndrv.ne.0.and.ifrep(ndrv).gt.35)then
               ic=0
            else 
               ic=ndrv
            endif
            if(iord.eq.0)then
               ans=xsect2(kpt+90*ic,1)
            else
               ans=(xsect2(kpt+90*ic,1)+xsect2(kpt+90*ic,2))
            endif
c 
c Convert to nb/GeV
c
            ans=ans/1000000.
         else if(it.eq.203)then
c
c  D0 Run II cone algorithm as calculated by fastNLO
c
c  Nonperturbative corrections included
c
            if(npt.eq.n203)then
c               ijet=1
c               smucoef=1.0
               if(ndrv.ne.0)then
                  if(ifrep(ndrv).gt.35) goto 253
               endif
               smucoef=1.0       
               fastnlo_table = fln(1:ifastnlo)//'/fnt2009midp.tab'
               call ft2009cc(fastnlo_table,smucoef,smucoef,0,temp3)
               do jn=1,110
                  do jm=1,3
                     xsect3(jn+110*ndrv,jm)=temp3(jn,jm)
                  enddo
               enddo
            endif
 253        continue
            kpt=npt-n203+1
            ic=0
            if(ndrv.ne.0.and.ifrep(ndrv).gt.35)then
               ic=0
            else 
               ic=ndrv
            endif
            if(iord.eq.0)then
               ans=xsect3(kpt+110*ic,1)
     &              *(1.d0+anp(kpt)/100.)
            else
               ans=(xsect3(kpt+110*ic,1)+xsect3(kpt+110*ic,2))
     &              *(1.d0+anp(kpt)/100.)
            endif
c 
c Convert to nb/GeV
c
            ans=ans/1000.
         else if(it.eq.204)then         
            if(npt.eq.n204)then
c               ijet=1
c               smucoef=1.0
               if(ndrv.ne.0)then
                  if(ifrep(ndrv).gt.35) goto 254
               endif
               smucoef=1.0       
               fastnlo_table = fln(1:ifastnlo)//'/fnt2007midp.tab'
               call fx9999cc(fastnlo_table,smucoef,smucoef,0,temp4)
c
c The fastNLO scenario for the CDF Run II data has an extra pt point 
c for each rapidity bin. Skip these points when filling the array xsect4
c
               jns=0
               do jn=1,77
               jns=jns+1
                  if(jn.eq.1.or.jn.eq.18.or.jn.eq.35.or.jn.eq.51
     &              .or.jn.eq.66)then
                     jns=jns-1
                     goto 354
                  endif
                  do jm=1,3
                     xsect4(jns+72*ndrv,jm)=temp4(jn,jm)
                  enddo
 354              continue
               enddo
            endif
 254        continue
            kpt=npt-n204+1
            ic=0
            if(ndrv.ne.0.and.ifrep(ndrv).gt.35)then
               ic=0
            else 
               ic=ndrv
            endif
            if(iord.eq.0)then
               ans=xsect4(kpt+72*ic,1)
            else
               ans=xsect4(kpt+72*ic,1)+xsect4(kpt+72*ic,2)
            endif
c
c convert to nb/GeV
c
            ans=ans/1000.
         endif         
      else if(it.gt.300.and.it.lt.400)then
c
c  gamma + jet
c
         rs=v(1)
         pt=v(2)
         nregion=v(3)
c
c  hardwire scale for now
c
         scale=1.0
         q2=(scale*pt)**2
         s0=dlog(q02/xc(1)**2)
         als=dlog(dlog(q2/xc(1)**2)/s0)
         s=rs**2
         call d0gamjet(s,pt,nregion,ndrv,als,ans,scale)
         if(ans.eq.0.d0)then
            print*,'* ',pt,nregion
         endif
      endif   
      IN=INORM(index)
      THEORY=ANS/XC(IN)
      RETURN
      END 



*****************************************************
*     structure function routines have been moved
*     to 'DIS10.f'
*****************************************************



      SUBROUTINE DYANTH(NDRV,V,XC,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(4),XC(100)
      COMMON/Q2STUF/Q02,Q2MAX 
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GRPTHY/FLAVOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      common/flags/itype(49),nflag(49),inorm(49),itgt(49),icorr(49)
     2     ,itmc(49),iht(49)
      index=V(4)
      IFL=nflag(index)
      INU=itgt(index)
      S=V(1)
      Q2=V(2)**2
      IF(IFL.EQ.5) Q2=Q2/2.
      ALS0=DLOG(Q02/XC(1)**2) 
      ALS=DLOG(DLOG(Q2/XC(1)**2)/ALS0)
      Y=V(3)
      CALL HODY(IFL,INU,NDRV,S,Q2,Y,ALS,THEORY)
      RETURN
      END 


      SUBROUTINE HODY(IFL,INU,NDRV,S,QS,Y,ALS,THEORY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Q10(-5:5),Q1(-5:5,32),Q20(-5:5),Q2(-5:5,32),QT(-5:5)
     2,HQ1Q(32),HQQ2(32),HQG2(32),HG1Q(32)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/CONSTANTS/PI,PI2
      COMMON/CURPAR/XC(100)
      COMMON/DYTEST/XSECLO,GLUCOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      common/e866xf/ixfx
C
C  HIGHER ORDER DRELL-YAN SUBROUTINE ADDED APRIL 29, 1991
C  CALCULATES S**1.5 DSIGMA/DQ/DY
C
c  modifications added to allow calculation of s**1.5 dsigma/dq dxf
c  for E866  3/5/03
c
      TAU=QS/S
      if(ixfx.eq.1)then
         xf=y
         rt=sqrt(xf**2+4.*tau)
         x10=(xf+rt)/2.
         x20=(-xf+rt)/2.
         xffac=1./(x10+x20)
         qu=4./9.
         qd=1./9.
      else
         EX=EXP(Y)
         X10=SQRT(TAU)*EX
         X20=X10/EX**2
         xffac=1.
         qu=4./9.
         qd=1./9.
      endif
c
c  Z production
c
      if(ixfx.eq.2.or.ixfx.eq.3)then
         sin2thw=0.2312
         qu=.25+(0.5-4./3.*sin2thw)**2
         qd=.25+(-0.5+2./3.*sin2thw)**2
      endif 
      al=alpha_s(iord+1,qs,xc(1),neff)
      flavor=neff
      ACF=AL*4./3./(2.*PI)
      AT=3.*ACF/4.
      if(ixfx.eq.2.or.ixfx.eq.3)then
         GF=1.1664e-05
         XMZ=91.1
         BR=0.034
         fac=BR*pi*1.414/3.*GF*XMZ**2/s*389.e06
      else
         FAC=8.*PI/137.**2/9./SQRT(TAU)*389.E03
      endif
      CALL HOQUARK(1,IFL,INU,X10,ALS,NDRV,Q10)
      CALL HOQUARK(2,IFL,INU,X20,ALS,NDRV,Q20)
      HQQ=(qu*Q10(1)*Q20(-1)+qd*(Q10(2)*Q20(-2)+Q10(3)*Q20(-3))
     2    +qu*Q10(-1)*Q20(1)+qd*(Q10(-2)*Q20(2)+Q10(-3)*Q20(3))
     3    +2.*qu*Q10(4)*Q20(4)+2.*qd*q10(5)*q20(5))
      THEORY=FAC*HQQ*(1.+
     2       ACF*(-8.+PI2+DLOG(X10*X20/(1.-X10)/(1.-X20))**2))
      theory=theory*xffac
      XSECLO=FAC*HQQ*xffac
      tmp0=(qs/s)**1.5*xseclo
      if(iord.eq.0)then
         theory=xseclo
         return
      endif
      QQ1=0.
      QG1=0.
      DO 100 I=1,NTERMS
      X2=.5*(1.-X20)*XI(I)+.5*(1.+X20)
      Z=X20/X2
      CALL HOQUARK(2,IFL,INU,X2,ALS,NDRV,QT)
      DO 101 IQ=-5,5
  101 Q2(IQ,I)=QT(IQ)
      HQQ2(I)=(qu*Q10(1)*Q2(-1,I)+qd*(Q10(2)*Q2(-2,I)+Q10(3)*Q2(-3,I))
     2    +qu*Q10(-1)*Q2(1,I)+qd*(Q10(-2)*Q2(2,I)+Q10(-3)*Q2(3,I))
     3    +2.*qu*Q10(4)*Q2(4,I)+2.*qd*q10(5)*q2(5,i))
C        1         2         3         4         5         6         7 *
      if(ixfx.eq.1)then
         TEMP=((1.+Z*Z)*DLOG((x10+x20)*(1.-X10)/X10/x20/(X2+X10))
     2        *HQQ2(I)-2.*DLOG((1.-X10)/X10/X20)
     3         *HQQ)/(X2-X20)+2.*DLOG(X2-X20)/(X2-X20)*(Z*HQQ2(I)-HQQ)
     4         +HQQ2(I)*(X2-X20)/X2**2*(1.+DLOG(X2-X20))
      else
         TEMP=((1.+Z*Z)*DLOG(2.*(1.-X10)/X10/(X2+X20))*HQQ2(I)
     2        -2.*DLOG((1.-X10)/X10/X20)
     3         *HQQ)/(X2-X20)+2.*DLOG(X2-X20)/(X2-X20)*(Z*HQQ2(I)-HQQ)
     4         +HQQ2(I)*(X2-X20)/X2**2*(1.+DLOG(X2-X20))
      endif
      QQ1=QQ1+.5*(1.-X20)*WI(I)*FAC*ACF*TEMP*xffac
      HQG2(I)=(qu*Q10(1)+qd*(Q10(2)+Q10(3))+qu*Q10(-1)
     2        +qd*(Q10(-2)+Q10(-3))
     2        +2.*qu*Q10(4)+2.*qd*q10(5))*Q2(0,I)
      if(ixfx.eq.1)then
         TEMP=((X20**2+(X2-X20)**2)/(2.*X2**3)*
     2         DLOG((X2-X20)*(1.-X10)*(x10+x20)/X10/x20/(X2+X10))
     3         +X20*(X2-X20)/X2**3)*HQG2(I)
      else
         TEMP=((X20**2+(X2-X20)**2)/(2.*X2**3)*
     2         DLOG(2.*(X2-X20)*(1.-X10)/X10/(X2+X20))
     3         +X20*(X2-X20)/X2**3)*HQG2(I)
      endif
      QG1=QG1+.5*(1.-X20)*WI(I)*FAC*AT*TEMP*xffac
  100 CONTINUE
      QQ2=0.
      QG2=0.
      DO 200 I=1,NTERMS
      X1=.5*(1.-X10)*XI(I)+.5*(1.+X10)
      Z=X10/X1
      CALL HOQUARK(1,IFL,INU,X1,ALS,NDRV,QT)
      DO 201 IQ=-5,5
  201 Q1(IQ,I)=QT(IQ)
      HQ1Q(I)=(qu*Q1(1,I)*Q20(-1)+qd*(Q1(2,I)*Q20(-2)+Q1(3,I)*Q20(-3))
     2    +qu*Q1(-1,I)*Q20(1)+qd*(Q1(-2,I)*Q20(2)+Q1(-3,I)*Q20(3))
     3    +2.*qu*Q1(4,I)*Q20(4)+2.*qd*q1(5,i)*q20(5))
      if(ixfx.eq.1)then
         TEMP=((1.+Z*Z)*DLOG((x10+x20)*(1.-X20)/x10/X20/(X1+X20))
     2        *HQ1Q(I)-2.*DLOG((1.-X20)/X10/X20)
     3        *HQQ)/(X1-X10)+2.*DLOG(X1-X10)/(X1-X10)*(Z*HQ1Q(I)-HQQ)
     4        +HQ1Q(I)*(X1-X10)/X1**2*(1.+DLOG(X1-X10))
      else
         TEMP=((1.+Z*Z)*DLOG(2.*(1.-X20)/X20/(X1+X10))*HQ1Q(I)
     2         -2.*DLOG((1.-X20)/X10/X20)
     3         *HQQ)/(X1-X10)+2.*DLOG(X1-X10)/(X1-X10)*(Z*HQ1Q(I)-HQQ)
     4         +HQ1Q(I)*(X1-X10)/X1**2*(1.+DLOG(X1-X10))
      endif
      QQ2=QQ2+.5*(1.-X10)*WI(I)*FAC*ACF*TEMP*xffac
      HG1Q(I)=Q1(0,I)*(qu*Q20(1)+qd*(Q20(2)+Q20(3))+qu*Q20(-1)
     2       +qd*(Q20(-2)+Q20(-3))+2.*qu*Q20(4)+2.*qd*q20(5))
      if(ixfx.eq.1)then
         TEMP=((X10**2+(X1-X10)**2)/(2.*X1**3)*
     2         DLOG((x10+x20)*(X1-X10)*(1.-X20)/x10/X20/(X1+X20))
     3         +X10*(X1-X10)/X1**3)*HG1Q(I)
      else
         TEMP=((X10**2+(X1-X10)**2)/(2.*X1**3)*
     2         DLOG(2.*(X1-X10)*(1.-X20)/X20/(X1+X10))
     3         +X10*(X1-X10)/X1**3)*HG1Q(I)
      endif
      QG2=QG2+.5*(1.-X10)*WI(I)*FAC*AT*TEMP*xffac
  200 CONTINUE
C        1         2         3         4         5         6         7 *
      QQ3=0.
      QG3=0.
      DO 300 I=1,NTERMS
      X2=.5*(1.-X20)*XI(I)+.5*(1.+X20)
      TEMPQQ=0.
      TEMPQG=0.
      TEMPGQ=0.
      DO 400 J=1,NTERMS
      X1=.5*(1.-X10)*XI(J)+.5*(1.+X10)
      HQ1Q2=(qu*Q1(1,J)*Q2(-1,I)+qd*(Q1(2,J)*Q2(-2,I)+Q1(3,J)*Q2(-3,I))
     2      +qu*Q1(-1,J)*Q2(1,I)+qd*(Q1(-2,J)*Q2(2,I)+Q1(-3,J)*Q2(3,I))
     3      +2.*qu*Q1(4,J)*Q2(4,I)+2.*qd*q1(5,j)*q2(5,i))
      HQ1G2=(qu*Q1(1,J)+qd*(Q1(2,J)+Q1(3,J))+qu*Q1(-1,J)
     2       +qd*(Q1(-2,J)+Q1(-3,J))
     2       +2.*qu*Q1(4,J)+2.*qd*q1(5,j))*Q2(0,I)
      HG1Q2=(qu*Q2(1,I)+qd*(Q2(2,I)+Q2(3,I))+qu*Q2(-1,I)
     2       +qd*(Q2(-2,I)+Q2(-3,I))
     2      +2.*qu*Q2(4,I)+2.*qd*q2(5,i))*Q1(0,J)
      TEMP=HQ1Q2*HA(X1,X2,X10,X20)+(GA(X1,X2,X10,X20)*HQ1Q2
     2    +GA(X10,X20,X10,X20)*HQQ-GA(X1,X20,X10,X20)*HQ1Q(J)
     3    -GA(X10,X2,X10,X20)*HQQ2(I))/(X1-X10)/(X2-X20)
      TEMPQQ=TEMPQQ+.5*(1.-X10)*WI(J)*TEMP*FAC*ACF*2.
      TEMP=HQ1G2*HC(X1,X2,X10,X20)+(GC(X1,X2,X10,X20)*HQ1G2
     2+   -GC(X10,X2,X10,X20)*HQG2(I))/(X1-X10)
      TEMPQG=TEMPQG+.5*(1.-X10)*WI(J)*TEMP*FAC*AT
      TEMP=HG1Q2*HC(X2,X1,X20,X10)+(GC(X2,X1,X20,X10)*HG1Q2
     2+   -GC(X20,X1,X20,X10)*HG1Q(J))/(X2-X20)
      TEMPGQ=TEMPGQ+.5*(1.-X10)*WI(J)*TEMP*FAC*AT
  400 CONTINUE
      QQ3=QQ3+.5*(1.-X20)*WI(I)*TEMPQQ
      QG3=QG3+.5*(1.-X20)*WI(I)*(TEMPQG+TEMPGQ)
  300 CONTINUE
      THEORY=THEORY+QQ1+QG1+QQ2+QG2+QQ3+QG3
      GLUCOR=QG1+QG2+QG3
      RETURN
      END


      FUNCTION GA(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=(x1+x2)*(tau**2+(x1*x2)**2)
         ga=temp/2./(x1*x2)**2/(x1+x20)/(x2+x10)
      else
         TEMP=(X1*X2+TAU)*(TAU**2+(X1*X2)**2)
         GA=TEMP/(X1*X2)**2/(X1+X10)/(X2+X20)
      endif
      RETURN
      END
      FUNCTION HA(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         ha=-1./x1/x2/(x1+x2)
      else
         HA=-2.*TAU*(X1*X2+TAU)/(X1*X2)/(X10*X2+X20*X1)**2
      endif
      RETURN
      END


      FUNCTION GC(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
c
c  this term corresponds to my notation with 1/(x1-x10)_R
c  so 1<-->2 compared to KLMP here and in HC
c 
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=tau**2+(x1*x2-tau)**2
         gc=temp/2./x1**2/x2**3/(x1+x20)
      else
         TEMP=X10*(X1*X2+TAU)*(TAU**2+(TAU-X1*X2)**2)
         GC=TEMP/X1**2/X2**3/(X1+X10)/(X10*X2+X20*X1)
      endif
      RETURN
      END
      FUNCTION HC(X1,X2,X10,X20)
      implicit real*8 (a-h,o-z)
      common/e866xf/ixfx
      TAU=X10*X20
      if(ixfx.eq.1)then
         temp=x2*(x1+x20)*(x1-x10)+2.*tau*(x1+x2)
         hc=temp/2./(x1*x2)**2/(x1+x2)**2
      else
         TEMP=X20*X2*X1**2+TAU*(X10*X2+2.*X1*X20)
         HC=TEMP*TAU*(X1*X2+TAU)/(X1*X2)**2/(X10*X2+X20*X1)**3
      endif
      RETURN
      END


      SUBROUTINE HOQUARK(INDEX,IFL,INU,X,ALS,NDRV,QUARK)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION QUARK(-5:5)
      double precision xi(96),WI(96),xx(97)
      common/grpthy/flavor
      common/target/az,an
      common/e866xf/ixfx
      COMMON/GAUS96/XI,WI,NTERMS,XX
C
c  unfolds quark distributions
c  
c  index = 1  --> "projectile" PDFs
c        = 2  --> "target" PDFs
c  IFL   = INACTIVE
c  inu   = nuclear correction flag
c  x     = (dp) parton's momentum fraction
c  s     = (dp) DGLAP s variable
c  ndrv  = (i) PDF parameter being varied
c  quark(-5:5) = (dp) PDF array (-5 = bbar .... -1=ubar, 0=g, ....)

c<<<<<<< theory13.f
c      call fsupdf(ndrv,x,als,u,d,ub,db,sb,cb,bb,glue)
c      QUARK(0)=GLUE/X
c=======
      ! selects smearing model
      call split_nuke(inu,i0,i1,ir,ioff,iwfn,ism)
c
c Disable smearing for now
c
      !ism=0  ! remove to apply smearing
      if ((az.ne.an).or.(index.eq.1)) then
         ism=0                     ! only apply smearing to target deuterons
      end if

      if (ism.le.1) then           ! no smearing, no offshell
         call fsupdf(ndrv,x,als,u,d,ub,db,sb,cb,bb,glue)
      else if (ism.ge.2) then      ! smearing, allows for fmKP offshell
         xD = x/2                           ! deuteron Bjorken variable
         yDmax = 1d0 
         yDmin = xD	
         u = 0
         d = 0
         ub = 0
         db = 0
         sb = 0
         cb = 0
         bb = 0
         glue  = 0
         if (ioff.ge.5) then
            uo = 0
            do = 0
            ubo = 0
            dbo = 0
            sbo = 0
            cbo = 0
            bbo = 0
            glueo  = 0
         endif
         DO I=1,NTERMS
            yD=0.5*(yDmax-yDmin)*XI(I)+0.5*(yDmax+yDmin) ! y = yD in [xD,1]
            xoy=xD/yD
            fy_diag = PHI_INT2D (yD,1d0,2,ism,iwfn)
            call fsupdf(ndrv,xoy,als,up,dp,ubp,dbp,sbp,cbp,bbp,gluep)
            if (ioff.ge.5) then    ! (fm)KP off-shell model
               foff_diag = PHI_INT2D (yD,1d0,12,ism,iwfn)
               call set_offshell_on (.true.)
               call fsupdf(ndrv,xoy,als,uop,dop,ubop,dbop,
     &                     sbop,cbop,bbop,glueop)
               call set_offshell_on (.false.)
            endif

            if (xoy.lt.0.995d0) then
               weight = .5*(yDmax-yDmin)*WI(I)
               u    = u    + weight * fy_diag*up
               d    = d    + weight * fy_diag*dp
               ub   = ub   + weight * fy_diag*ubp
               db   = db   + weight * fy_diag*dbp
               sb   = sb   + weight * fy_diag*sbp
               cb   = cb   + weight * fy_diag*cbp
               bb   = bb   + weight * fy_diag*bbp
               glue = glue + weight * fy_diag*gluep
               if (ioff.ge.5) then !(fm)KP offshell correction
                  uo    = uo    + weight * foff_diag*uop
                  do    = do    + weight * foff_diag*dop
                  ubo   = ubo   + weight * foff_diag*ubop
                  dbo   = dbo   + weight * foff_diag*dbop
                  sbo   = sbo   + weight * foff_diag*sbop
                  cbo   = cbo   + weight * foff_diag*cbop
                  bbo   = bbo   + weight * foff_diag*bbop
                  glueo = glueo + weight * foff_diag*glueop
               endif
            end if
         ENDDO
         if (ioff.ge.5) then ! apply offshell corrections to pdfs
            u=u+uo
            d=d+do
            ub=ub+ubo
            db=db+dbo
            sb=sb+sbo
            cb=cb+cbo
            bb=bb+bbo
            glue=glue+glueo
            !write(*,*) '* ',x,uo/u,do/d,ubo/ub,dbo/db,glueo/glue
         end if 
      end if
      QUARK(0)=GLUE/X
      quark(5)=bb/x
      QUARK(4)=cb/X
      QUARK(3)=sb/X
      quark(-5)=quark(5)
      quark(-4)=quark(4)
      QUARK(-3)=QUARK(3)

      IF(INDEX.EQ.1) THEN ! projectile: always a proton
         QUARK(-1)=ub/X
         QUARK(-2)=db/X
         QUARK(1)=u/x
         QUARK(2)=d/x
      ELSE IF(INDEX.EQ.2) THEN  ! target
         TUB2=ub/x              ! isospin average 
         TDB2=db/x
         TU2=u/x
         TD2=d/x
         QUARK(-1)=AZ*TUB2+AN*TDB2
         QUARK(-2)=AZ*TDB2+AN*TUB2
         QUARK(1)=AZ*TU2+AN*TD2
         QUARK(2)=AZ*TD2+AN*TU2
c
c  ixfx=2 for p pbar Z production
c  interchange u, ubar and d, dbar
c
         if(ixfx.eq.2)then      ! exchanges u and ubar, d and dbar
            tmp=quark(-1)
            quark(-1)=quark(1)
            quark(1)=tmp
            tmp=quark(-2)
            quark(-2)=quark(2)
            quark(2)=tmp
         endif

      ENDIF

      RETURN
      END

      SUBROUTINE TEST_EV(XC)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XC(100),FV(32),XMSTV(10),CN31(10),DN(10),DN1(10)
     2,V(4)
      COMMON/Q2STUF/Q02,Q2MAX 
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/FLAVOR
      COMMON/NORMAL/INORM(49)
      DATA CN31/0.,4.282,6.047,7.209,8.094,8.820,9.440,9.983,
     210.468,10.907/
      DATA DN/0.,.427,.667,.837,.971,1.080,1.174,1.255,1.327,1.392/
      COMMON/CONSTANTS/PI,PI2
C
C  THIS ROUTINE TESTS THE ACCURACY OF THE EVOLUTION ROUTINE BY
C  COMPARING THE EVOLVED MOMENTS WITH THE MOMENTS OF THE EVOLVED
C  DISTRIBUTIONS. AT THIS TIME IT IS IMPLEMENTED ONLY FOR THE
C  NONSINGLET DISTRIBUTION X(UV+DV).
C  FOR TEST PURPOSES SOME SINGLET AND GLUON DISTRIBUTIONS ARE
C  PRINTED OUT AS WELL.
C
C      IF(INS) 200,200,100
  100 CONTINUE
      WRITE(6,1001)
 1001 FORMAT(///,' NON-SINGLET MOMENTS TEST') 
      WRITE(6,1002)
 1002 FORMAT(5X,' N',5X,'EVOLVED MOMENT',1X,'INTQCD MOMENT',
     24X,'PERCENT ERROR') 
      B0=11.-2.*FLAVOR/3.
      B1=102.-38.*FLAVOR/3.
      DO 10 N=1,10
      DN1(N)=CN31(N)*B0/B1
   10 CONTINUE
c      T=DLOG(Q02/XC(1)**2)
c      A0=ALPHA(T)
      a0=alpha_s(iord+1,q02,xc(1),neff)
      flavor=neff
      ALB0=1.+B1/B0*A0/4./PI
C      DO 1 I=1,NTERMS
C      V(1)=XX(I)
C      V(2)=Q02
C      V(4)=7
C      FV(I)=THEORY(1,I,V,XC)
C    1 CONTINUE
C      DO 2 N=2,10
C      XMSTV(N)=XMNT32(FV,N)
C    2 CONTINUE
      C=3./(BETA(XC(2),XC(3)+1.)+XC(4)*BETA(XC(2)+1.,XC(3)+1.)
     2+XC(5)*BETA(XC(2)+2.,XC(3)+1.))
      DO 2 N=1,10
      XMSTV(N)=C*(BETA(XC(2)+N-1.,XC(3)+1.)+XC(4)*BETA(XC(2)+N,XC(3)+1.)
     2+XC(5)*BETA(XC(2)+N+1.,XC(3)+1.))/XC(INORM(7))
    2 CONTINUE
c      T=DLOG(Q2MAX/XC(1)**2)
c      AL=ALPHA(T)
      al=alpha_s(iord+1,q2,xc(1),neff)
      flavor=neff
      ALB=1.+B1*AL/B0/4./PI
      DO 4 N=1,10
      IF(N.LE.2) THEN
      DO 3 I=1,NTERMS
      V(1)=XX(I)
C
C  MODIFIED N=1 MOMENT INTEGRATION
C
      IF(N.EQ.1) V(1)=XX(I)**2
      V(2)=Q2MAX
      V(4)=7
      FV(I)=THEORY(1,I,V,XC)
    3 CONTINUE
      ENDIF
      XMV=XMNT32(FV,N)
      S1=SS(N,1)
      S2=SS(N,2)
      CN=4./3.*(2.*S1**2-2.*S2+3.*S1-2.*S1/N/(N+1.)
     2+1./N+2./(N+1.)+2./N/N-9.)
      IF(IORD) 6,6,7
    6 CONTINUE
      ALB=1.
      ALB0=1.
      CN=0.
    7 CONTINUE
      TV=XMSTV(N)*(AL/A0)**DN(N)*(ALB/ALB0)**(DN1(N)-DN(N)) 
C      TV=TV*(1.+CN*AL/4./PI)/(1.+CN*A0/4./PI)
      TV=TV*(1.+CN*AL/4./PI)
      ERV=(TV-XMV)/TV*100.
    4 WRITE(6,5) N,TV,XMV,ERV
    5 FORMAT(5X,I2,5X,3E15.4) 
      RETURN
C  200 CONTINUE
C      WRITE(6,201)
C  201 FORMAT(///,' SINGLET AND GLUON TEST OUTPUT')
C      ITIME=1
C      Q2=Q02
C      DO 211 K=1,2
C      S=DLOG(DLOG(Q2/XC(1)**2)/DLOG(Q02/XC(1)**2))
C      WRITE(6,202) Q2
C  202 FORMAT(///,' Q2=',F8.2) 
C      WRITE(6,203)
C  203 FORMAT(/,5X,' X',8X,'SINGLET',9X,'GLUON',12X,'F2',
C     212X,'XF3',12X,'F2EN') 
C      DO 210 J=1,19 
C      X=.05*J
C      V(1)=X
C      V(2)=Q2
C      V(4)=6
C      F2=THEORY(1,ITIME,V,XC) 
C      V(4)=7
C      ITIME=2
C      XF3=THEORY(1,ITIME,V,XC)
C     V(4)=3
C     F2EN=THEORY(1,ITIME,V,XC)
C      F2EN=0.
C      CALL GINTERP(7,0,X,S,SING)
C      CALL GINTERP(8,0,X,S,GLUE)
C      WRITE(6,204) X,SING,GLUE,F2,XF3,F2EN 
C  204 FORMAT(F10.4,5E15.4)
C  210 CONTINUE
C      Q2=Q2MAX
C  211 CONTINUE
C      RETURN
      END 


      FUNCTION SS(N,I)
      implicit real*8 (a-h,o-z)
      SS=0.
      DO 1 J=1,N
      SS=SS+1./J**I 
    1 CONTINUE
      RETURN
      END 


      SUBROUTINE PLOTLL(X1,X2,Y1,Y2,X,Y,ER,XTH,TH,NEXP,NTH) 
      implicit real*8 (a-h,o-z)
      REAL LN(100)
      DIMENSION X(1),Y(1),ER(1),XTH(1),TH(1)
      DATA AB,AS,AX,AE,AP,AI/1H ,1H*,1HX,1H-,1H+,1HI/
C
C  THIS IS A RELATIVELY CRUDE PLOTTER WHICH WORKS ON A
C  LINE PRINTER. IT CAN BE MODERATELY USEFUL BUT IT
C  DOES RESULT IN A LOT OF OUTPUT IF YOU FIT A LOT OF
C  DIFFERENT X POINTS.
C
      B=100./(DLOG(X2)-DLOG(X1))
      A=1.-B*DLOG(X1)
      D=49./(DLOG(Y1)-DLOG(Y2))
      C=1.-D*DLOG(Y2)
      WRITE(6,9)
    9 FORMAT(1H1)
      DO 1 I=1,100
    1 LN(I)=AS
      DO 2 I=10,90,10
    2 LN(I)=AP
      WRITE(6,7) (LN(I),I=1,100) 
    7 FORMAT(1H ,15X,100A1)
      DO 100 I=1,50 
      LN(1)=AS
      LN(100)=AS
      DO 5 J=2,99
    5 LN(J)=AB
      YY=EXP((I-C)/D)
      DO 99 K=1,NEXP
      IX=A+B*DLOG(X(K))
      IYU=C+D*DLOG(Y(K)+ER(K))
      IY=C+D*DLOG(Y(K))
      IF(ER(K).GE.Y(K)) ER(K)=Y(K)-.001 
      IYL=C+D*DLOG(Y(K)-ER(K))
      IF(I.EQ.IYU.OR.I.EQ.IYL) LN(IX)=AE
      IF(I.EQ.IY) LN(IX)=AX
   99 CONTINUE
      DO 98 L=1,NTH 
      IX=A+B*DLOG(XTH(L))
      ITH=C+D*DLOG(TH(L))
      IF(I.EQ.ITH) LN(IX)=AS
   98 CONTINUE
  100 WRITE(6,3) YY,(LN(K),K=1,100)
    3 FORMAT(1H ,F10.4,5X,100A1)
      DO 6 I=1,100
      LN(I)=AS
      IP=I/10*10
      IF(I.EQ.IP) LN(I)=AP
    6 CONTINUE
      LN(100)=AS
      WRITE(6,7) (LN(K),K=1,100) 
      WRITE(6,8) X1,X2 
    8 FORMAT(1H ,10X,F10.3,90X,F10.3)
      RETURN
      END 



      SUBROUTINE WATE4
      implicit real*8 (a-h,o-z)
      COMMON/GAUSS4/XI(4),WI(4),NTERMS,XX(5)
      NTERMS=4
C
C  4 POINT GAUSSIAN QUADRATURE
C
      XI(1)=-.8611363116
      XI(2)=-.3399810436
      XI(3)= .3399810436
      XI(4)= .8611363116
      WI(1)=.3478548451
      WI(2)=.6521451549
      WI(3)=.6521451549
      WI(4)=.3478548451
      DO 1 J=1,4
    1 XX(J)=.5*(XI(J)+1.)
      XX(5)=1.
      RETURN
      END 


      subroutine wasym(rs,y,ALS,NDRV,answer,it)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension xa(24),argp(24),argm(24),qa(-5:5),qb(-5:5)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB

C      common/qs/q2
C      common/cteqpdf/ist
c
c  calculates the w lepton asymmetry. See Barger and Phillips 
c  "Collider Physics" pgs 254-256. Note: theta-hat on pgs 255, 256 
c  is - (theta-hat) on page 254!  
c
c  J.F. Owens    June 30, 1994
c              & April 1, 2009 (for it=131 CDF 2009 data)
c
c  rs = sqrt(s)  y=lepton rapidity
c
c  Returns CDF W asymmetry (it=131) or W lepton asymmetries for 
c  CDF or D0 (it=126-130)
c
c  1/8/13 Modified to return ATLAS W+ (it=151) or W- (it=152)
c  lepton rapidity distributions or the W lepton asymmetry (it=153)
c  K factors from MCFM
c
      Ak1ResCt4m(x)= 1.07949- .0041764*x+ .017066*x**2+ .016434*x**3
     >    + .00074073*x**4 - .0015393*x**5 - .00040555*x**6
      CDF1800(x)= 1.0904 + 0.02835 * x + 0.023048 * x**2 
     >- 0.0075403 * x**3 - 0.0022601 * x**4 + 0.0020092 * x**5
c      CDF1960(x) = 1.0795 + 0.014884 * x + 0.025591 * x**2 
c     >+ 0.00015113 * x**3 - 0.0029232 * x**4 + 0.00082595 * x**5
      D01960(x) = 1.1209 + 0.022313 * x + 0.018531 * x**2 
     >- 0.001236 * x**3 - 0.0019652 * x**4 + 0.00041654 * x**5
      Akf2525(x)=1.0962 + 0.0067064*x + 0.016858*x**2 + 0.0077185*x**3 
     >- 0.0016489*x**4 - 0.00052809*x**5
      Akf3535(x)=0.85241 - 0.0096781*x + 0.033227*x**2 
     >+0.010615*x**3-0.002623*x**4-0.00049916*x**5
      ATLASp(x) = 1.1073 + 0.0039 * x
      ATLASm(x) = 1.1030 + 0.0043 * x
c  cteq2m flag
      ist=1
c  ckm factors (using just the Cabibbo angle)
      c2tc=.95
      s2tc=.05
c  w mass
      xmw=80.
      q2=xmw**2
      s=rs**2
c  w width
      gam=2.141
c  Gf
      gf=1.166e-5
c normalization factor (cross section in picobarns)
      fac=1./3./8.*(gf*xmw**2/1.414)**2*xmw/gam/rs**2*389400000.
c
c ATLAS W lepton asymmetry - use it=151, 152
c Must now allow for ppbar (Tevatron) or pp (LHC)
c
c CDF W asymmetry is 131
c D0  W asymmetry is 132
      is=-1
      if(it.eq.151.or.it.eq.152.or.it.eq.153)is=1
      if(it.eq.131.or.it.eq.132)then
         r=xmw/rs
         xa0=r*exp(y)
         xb0=r*exp(-y)
         CALL HOQUARK(1,16,0,XA0,ALS,NDRV,QA)
         CALL HOQUARK(1,16,0,XB0,ALS,NDRV,QB)
c 
c  CDF W asymmetry
c
         sigp=c2tc*(qa(1)*qb(2)+qa(4)*qb(3)
     2   +(qa(-2)*qb(-1)+qa(-3)*qb(4)))
     3   +s2tc*(qa(1)*qb(3)+qa(4)*qb(2)
     4   +qa(-3)*qb(-1)+qa(-2)*qb(4))
         sigm=c2tc*(qa(2)*qb(1)+qa(3)*qb(4)
     2   +qa(-1)*qb(-2)+qa(4)*qb(-3))
     3   +s2tc*(qa(3)*qb(1)+qa(2)*qb(4)
     4   +qa(-1)*qb(-3)+qa(4)*qb(-2))
         asym=(sigp-sigm)/(sigp+sigm)
         answer=asym
         return
      endif
c
c  CDF uses pt(muon)>25 x=xmw/(2.*ptmin)
c  Old D0 used 20
c  New D0 electron data uses 25
c  ATLAS uses 25 for missing mass
c  D0 muon data also has an ET>35 set
c  
      x=xmw/50.
      if(it.eq.129) x=xmw/40.
      if(it.eq.135) x=xmw/70.
      w=x+sqrt(x**2-1.)
      xi=Dlog(w)
c
c  In lowest order, the pt cut translates into simple 
c  limits on the xa integral
c
      xamin=xmw/rs*dexp(y-xi)
      xamax=xmw/rs*dexp(y+xi)
      if(xamax.gt.1.d0)xamax=1.d0
c
c  6-pt gaussian quadrature routine -- supplied below
c
      call gq11(xamin,xamax,4,xa,argp,ansp)
      do 100 j=1,24
      xb=xmw**2/(s*xa(j))
c
c  my parton distribution calls
c  replace with yours
c
C      call dist(xa(j),qa)
C      call dist(xb,qb)
      CALL HOQUARK(1,16,0,XA(J),ALS,NDRV,QA)
      CALL HOQUARK(1,16,0,XB,ALS,NDRV,QB)
      yh=y-.5*Dlog(xa(j)/xb)
      st=1./cosh(yh)
      ct=tanh(yh)
      facp=(1.-ct)**2*st**2*fac
      facm=(1.+ct)**2*st**2*fac
c
c  s=sbar, c=cbar, b=bbar used here
c

      argp(j)=c2tc*((qa(1)*qb(-2*is)+qa(4)*qb(-3))*facp
     2+(qa(-2)*qb(1*is)+qa(-3)*qb(4))*facm)
     3+s2tc*((qa(1)*qb(-3)+qa(4)*qb(-2*is))*facp
     4+(qa(-3)*qb(1*is)+qa(-2)*qb(4))*facm)
      argm(j)=c2tc*((qa(2)*qb(-1*is)+qa(3)*qb(-4))*facm
     2+(qa(-1)*qb(2*is)+qa(-4)*qb(-3))*facp)
     3+s2tc*((qa(3)*qb(-1*is)+qa(2)*qb(-4))*facm
     4+(qa(-1)*qb(3)+qa(-4)*qb(2*is))*facp)
      argp(j)=argp(j)/xa(j)
      argm(j)=argm(j)/xa(j)
  100 continue
      call gq11(xamin,xamax,0,xa,argp,ansp)
      call gq11(xamin,xamax,0,xa,argm,ansm)
      if(it.eq.151)then
         if(iord.eq.0)then
            answer=ansp
         else
            answer=ansp*ATLASp(y)
         endif
         return
      else if(it.eq.152)then
         if(iord.eq.0)then
            answer=ansm
         else
            answer=ansm*ATLASm(y)
         endif
         return
      else if(it.eq.153)then
         if(iord.gt.0)then
            ansp=ansp*ATLASp(y)
            ansm=ansm*ATLASm(y)
         endif
      else if(it.eq.126)then
         if(iord.gt.0)then
            ansp=ansp*Ak1ResCt4m(y)
            ansm=ansm*Ak1ResCt4m(-y)
         endif
      else if(it.eq.127)then
         if(iiord.gt.0)then
            ansp=ansp*CDF1800(y)
            ansm=ansm*CDF1800(-y)
         endif
c      else if(it.eq.128.or.it.eq.130)then
c         ansp=ansp*CDF1960(y)
c         ansm=ansm*CDF1960(-y)
      else if(it.eq.129.or.it.eq.134)then
         if(iord.gt.0)then
            ansp=ansp*D01960(y)
            ansm=ansm*D01960(-y)
         endif
      else if(it.eq.128.or.it.eq.130.or.it.eq.134)then
         if(iord.gt.0)then
            ansp=ansp*Akf2525(y)
            ansm=ansm*Akf2525(-y)
         endif
      else if(it.eq.135)then
         if(iord.gt.0)then
            ansp=ansp*Akf3535(y)
            ansm=ansm*Akf3535(-y)
         endif
      endif
      answer=(ansp-ansm)/(ansp+ansm)
      return
      end


      SUBROUTINE GQ11(XMIN,XMAX,N,X,Y,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(24),Y(24),V(6),R(6) 
      DATA V/0.03376524,0.16939531,0.38069041,0.61930959,0.83060469,
     10.96623476/,R/0.08566225,0.18038079,0.23395697,0.23395697,
     20.18038079,0.08566225/
      save d,nsave
      IF(N) 250,250,150 
  150 K=0 
      NSAVE=N 
      D=(XMAX-XMIN)/N 
      XL=XMIN-D 
      DO 200 I=1,N
      XL=XL+D 
      DO 200 J=1,6
      K=K+1 
  200 X(K)=XL+D*V(J)
      RETURN
  250 W=0.
      K=0 
      DO 300 I=1,NSAVE
      DO 300 J=1,6
      K=K+1 
  300 W=W+Y(K)*R(J) 
      W=W*D 
      RETURN
      END 
      FUNCTION BETA(X1,X2)
      IMPLICIT REAL*8 (A-H,O-Z)
C  EULER BETA FUNCTION SUBROUTINE
      CALL GAMMA(X1,G1,IER)
      CALL GAMMA(X2,G2,IER)
      X3=X1+X2
      CALL GAMMA(X3,G3,IER)
      BETA=G1*G2/G3 
      RETURN
      END 
