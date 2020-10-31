************************************************************************
*     ALTPAR08-ALTPAR10                                                *
*                                                                      * 
*     PDF parametrization and evolution                                *
*                                                                      *
*     HISTORY:                                                         *
*                                                                      *
*     altpar08   : based on altpar02                                   *
*     (12 Nov 08)  modified the d-quark parametrization to             *
*                  allow d/u->finite as x->1                           *
*                  NOTE: requirs a new d-quark parameter in position   *
*                  #32, called 'a6dv'.                                 *
*                                                                      *
*     altpar10   : new xgrid, new treatment of b threshold,            *
*     (12 Jan 10)  new evolution scheme, new interpolation scheme,     *
*                  NNLO splitting functions added - jfo                *
*                                                                      *
*     (3/22/11)    added new flag (ipconv) to select parametrization   *
*                  1=CTEQ6X  2=MSTW                                    *
************************************************************************


      subroutine findht(jmax,par,ipname,pwate)
*     finds which parameters are for HT and stores their index in nht

      implicit real*8 (a-h,o-z)
      character*20 ipname(100)
      dimension par(100),pwate(100)

      logical stopflag

*     F2 Higher-twist parameters
      integer nht(4),nhtn(4),ihtFL(4)
      common/ht/nht,nhtn
*     FL Higher-twist parameters      
      integer nhtFL(4),nhtFLn(4)
      common/htFL/nhtFL,nhtFLn

      do i=1,4
         nhtn(i) = -1
         nhtFL(i) = -1
         nhtFLn(i) = -1
         ihtFL(i) = 0
      end do

      do 40 j=1,jmax
         if(ipname(j).eq.'ht1') then
            nht(1)=j
         else if(ipname(j).eq.'ht2') then
            nht(2)=j
         else if(ipname(j).eq.'ht3') then
            nht(3)=j
         else if(ipname(j).eq.'ht4') then
            nht(4)=j
         else if(ipname(j).eq.'ht1n') then
            nhtn(1)=j
         else if(ipname(j).eq.'ht2n') then
            nhtn(2)=j
         else if(ipname(j).eq.'ht3n') then
            nhtn(3)=j
         else if(ipname(j).eq.'ht4n') then
            nhtn(4)=j
         else if(ipname(j).eq.'htFL1') then
            nhtFL(1)=j
         else if(ipname(j).eq.'htFL2') then
            nhtFL(2)=j
         else if(ipname(j).eq.'htFL3') then
            nhtFL(3)=j
         else if(ipname(j).eq.'htFL4') then
            nhtFL(4)=j
         else if(ipname(j).eq.'htFL1n') then
            nhtFLn(1)=j
         else if(ipname(j).eq.'htFL2n') then
            nhtFLn(2)=j
         else if(ipname(j).eq.'htFL3n') then
            nhtFLn(3)=j
         else if(ipname(j).eq.'htFL4n') then
            nhtFLn(4)=j
         end if
 40   continue

*    ... checks that HT parameters for the proton are provided in input
      stopflag = .false.
      do i=1,4
         if (nht(i).eq.-1) then
            print*, 'ERROR(findht): HT paramter no.',i
     &           ,'absent in input file'
         end if
      end do
      if (stopflag) stop
 
*    ... sets HT(n) = HT(p) if no HT parameters present in the input list
*    ... or the parameter and its variation are both set to 0
      !print *, '** nhtn =', (nhtn(i),i=1,4)
      !print *, '** htn  =', (par(nhtn(i)),i=1,4)
      !print *, '** wate =', (pwate(nhtn(i)),i=1,4)
      do i=1,4
         if ( (nhtn(i).eq.-1)
     &        .or. (par(nhtn(i)).eq.0d0.and.pwate(nhtn(i)).eq.0d0) ) 
     &        then
            nhtn(i) = nht(i)
         end if
      end do
*     ... sets HT_FL(p)=HT(p) if no FL HT parameters are present in the input
*     ... list or ALL of the parameters and their variations are set to 0   
      do i=1,4
         if ( (nhtFL(i).eq.-1)
     &        .or. (par(nhtFL(i)).eq.0d0.and.pwate(nhtFL(i)).eq.0d0) ) 
     &        then
cjfo            nhtFL(i) = nht(i)
            ihtFL(i)=1 
         end if
      end do 
      itest=ihtFL(1)+ihtFL(2)+ihtFL(3)+ihtFL(4)
      if(itest.eq.4)then
         do i=1,4
            nhtFL(i)=nht(i)
         enddo
      endif
*     ... sets HT_FL(n)=HT_FL(p) if no FL n HT parameters are present in the
*     ... input list or the parameter and its variation are both set to 0   
      do i=1,4
         if ( (nhtFLn(i).eq.-1)
     &        .or. (par(nhtFLn(i)).eq.0d0.and.pwate(nhtFLn(i)).eq.0d0) ) 
     &        then
            nhtFLn(i) = nhtFL(i)
         end if
      end do 
      
      !print *, '** ht   =', (par(nht(i)),i=1,4)
      !print *, '** htn  =', (par(nhtn(i)),i=1,4)
      !print *, '** ht_FL   =', (par(nhtFL(i)),i=1,4)
      !print *, '** htn_FL  =', (par(nhtFLn(i)),i=1,4)

      return
      end


      SUBROUTINE INTQCD(XPAR) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FVL(60),FVL1(60),FVLSTO(60),PRD(60),PRD1(60),
     $ COR(60),COR1(60),FUDGE(60),PRDM(60)
      DIMENSION FS(60),FS1(60),FSSTO(60),PRDS(60),PRDS1(60),CORS(60), 
     $ CORS1(60),FUDGES(60),PRDMS(60)
      DIMENSION FG(60),FG1(60),FGSTO(60),PRDG(60),PRDG1(60),CORG(60), 
     $ CORG1(60),FUDGEG(60),PRDMG(60)
      DIMENSION XPAR(100),tmp(8,60),tmp2(8,60)
      COMMON/CONSTANTS/PI,PI2
      COMMON/ALQ2/T 
      COMMON/GRPTHY/FLAVOR
      COMMON/PARAM/PARA(100)
c      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB,NDEG
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/GFUNC/CALC(11,60,60)
      COMMON/GRID/NX,XGRID(60)
      common/threshold/sb
      common/stepsize/delta
      common/paramconv/ipconv
      save s0
c      data xmc,xmb,ncb/1.3d0,4.5d0,12/
c        1         2         3         4          5         6         7  *
      DATA NX,XGRID/60,.00001,.000015,.00002,.00004,.00006,.00008,
     2.0001,.00015,.0002,.0004,.0006,.0008,.001,.0015,.002,.004,.006,
     3.008,.01,.015,.02,.04,.06,.08,.1,.125,.15,.175,.2,.225,.25,.275,
     4.3,.325,.35,.375,.4,.425,.45,.475,.50,.525,.55,.575,.60,.625,
     5.65,.675,.70,.725,.75,.775,.80,.825,.85,.875,.90,.925,.95,.975/
C
C  UPDATED 7/27/05 TO MODIFY INITIALIZATION OF THE EVOLUTION TO REMOVE 
C  THE NEED FOR A BACKWARD STEP. ALSO REVAMPED FLAVOR THRESHOLD TREATMENT TO 
C  HANDLE MC<Q<MB AND Q>MB SEPARATELY. 
C
C  UPDATED 6/17/05 TO USE 3-LOOP SPLITTING FUNCTIONS
C  CHANGES TO ILOOP, IORD
C  ILOOP DETERMINES THE BETA FUNCTION FOR ALPHA_S (ILOOP=1 OR 2)
C  IORD DETERMINES THE ORDER OF THE SPLITTING FUNCTIONS (IORD=0, 1, OR 2)
C  IORD CONVENTION CHANGED (IORD-->IORD-1) TO AGREE WITH THE ALTPFIT 
C  CONVENTIONS.  1/15/10
C  CALC ENTRIES HAVE BEEN RENUMBERED TO MAKE THE NONSINGLET ARRAY STRUCTURE 
C  CLEANER
C
C  GFUNC IS THE ARRAY WHICH HOLDS THE EVOLVED DISTRIBUTIONS.
C  ALLOWANCE HAS BEEN MADE FOR EIGHT PARTON DISTRIBUTIONS.
C  THE SECOND DIMENSION IS THE NUMBER OF STEPS IN S
C  S=LN(LN(Q2/LAMBDA2)/LN(Q02/LAMBDA2)) 
C  THE THIRD DIMENSION IS THE NUMBER OF X POINTS POINTS IN THE GRID 
      DO 61 I=1,11
      DO 61 J=1,60
      DO 61 K=1,NX
   61 CALC(I,J,K)=0.D0
C  ALLOWANCE HAS BEEN MADE FOR UP TO 100 PARAMETERS
      DO 48 J=1,100             ! stores parameters in /param/ 
 48      PARA(J)=XPAR(J)        ! for use in evolution subroutines 
      PI=4.*DATAN(1.D0)
      PI2=PI**2
      S0=DLOG(XMC**2/PARA(1)**2)
c
c  set b threshold
c  moved to routine THEORY prior to the call to INTQCD
c
c      sb=dlog(dlog(XMB**2/para(1)**2)/s0)
c      DELTA=SB/NCB
c
c  renormalize uv, dv, g to satisfy sum rules
c
      para(2)=1.
      para(8)=1.
      para(25)=1.
      call pdfnrm
      DO 70 JFL=1,2
      IF(JFL.EQ.1)THEN
         SV=0.
         SST=S0
         N1=2
         N2=NCB
         NN=2
      ELSE
         SV=SB
         SST=S0*DEXP(SB)
         N1=NCB+3
         N2=NMAX-1
         NN=NCB+3
      ENDIF
C  CALCULATE THE SINGLET (FS) AND GLUON (FG) DISTRIBUTIONS 
      S=SV
C**  CALCULATE FS AND FG AT S=SV
      DO 1 I=1,NX
      X=XGRID(I)
      IF(JFL.EQ.1)THEN
         FS(I)=FCNFS(X)
         FG(I)=FCNFG(X)
         CALC(2,1,I)=FS(I)
         CALC(1,1,I)=FG(I)
c         tmp(1,i)=fs(i)
      ELSE
         FS(I)=CALC(2,NCB+1,I)
         FG(I)=CALC(1,NCB+1,I)
         CALC(2,NCB+2,I)=FS(I)
         CALC(1,NCB+2,I)=FG(I)
c         tmp2(1,i)=fs(i)
      ENDIF
    1 CONTINUE
C**  CALCULATE D(FS,FG)/DS AT S=SV
      T=SST
      DO 2 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,FS,FG,FS1(I))
      CALL FINTGG(X,FG,FS,FG1(I))
    2 CONTINUE
C**  CALCULATE PREDICTOR AT S=SV+DELTA 
      DO 3 I=1,NX
      PRDS(I)=FS(I)+DELTA*FS1(I)
      PRDG(I)=FG(I)+DELTA*FG1(I)
    3 CONTINUE
C**  CALCULATE D(PRD)/DS AT S=SV+DELTA 
      T=SST*DEXP(DELTA)
      DO 4 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,PRDS,PRDG,PRDS1(I)) 
      CALL FINTGG(X,PRDG,PRDS,PRDG1(I)) 
C**  CALCULATE CORRECTOR AT S=SV+DELTA 
      CORS(I)=FS(I)+0.5*DELTA*(FS1(I)+PRDS1(I))
c      if(jfl.eq.1)then
c         tmp(2,i)=fs1(i)
c         tmp(3,i)=prds1(i)
c         tmp(4,i)=cors(i)
c      else 
c         tmp2(2,i)=fs1(i)
c         tmp2(3,i)=prds1(i)
c         tmp2(4,i)=cors(i)
c      endif
      CORG(I)=FG(I)+0.5*DELTA*(FG1(I)+PRDG1(I))
      FUDGES(I)=PRDS(I)-CORS(I)
      FUDGEG(I)=PRDG(I)-CORG(I)
    4 CONTINUE
      DO 6 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,CORS,CORG,CORS1(I)) 
      CALL FINTGG(X,CORG,CORS,CORG1(I)) 
    6 CONTINUE
C**  RESHUFFLE THE DECK
      DO 7 I=1,NX
      FSSTO(I)=FS(I)
      FGSTO(I)=FG(I)
      FS(I)=CORS(I) 
      FG(I)=CORG(I) 
      FS1(I)=CORS1(I)
      FG1(I)=CORG1(I)
      CALC(2,NN,I)=FS(I)
      CALC(1,NN,I)=FG(I)
    7 CONTINUE
      DO 10 NTIMES=N1,N2
      NT1=NTIMES+1
      DO 11 I=1,NX
      PRDS(I)=FSSTO(I)+2.*DELTA*FS1(I)
      PRDG(I)=FGSTO(I)+2.*DELTA*FG1(I)
      PRDMS(I)=PRDS(I)-0.8*FUDGES(I)
      PRDMG(I)=PRDG(I)-0.8*FUDGEG(I)
   11 CONTINUE
      T=T*DEXP(DELTA)
      DO 12 I=1,NX
      X=XGRID(I)
      FSSTO(I)=FS(I)
      FGSTO(I)=FG(I)
      CALL FINTGS(X,PRDMS,PRDMG,PRDS1(I))
      CALL FINTGG(X,PRDMG,PRDMS,PRDG1(I))
      CORS(I)=FS(I)+0.5*DELTA*(FS1(I)+PRDS1(I))
      CORG(I)=FG(I)+0.5*DELTA*(FG1(I)+PRDG1(I))
      FUDGES(I)=PRDS(I)-CORS(I)
      FUDGEG(I)=PRDG(I)-CORG(I)
      FS(I)=CORS(I)+0.2*FUDGES(I)
      FG(I)=CORG(I)+0.2*FUDGEG(I)
      CALC(2,NT1,I)=FS(I)
      CALC(1,NT1,I)=FG(I)
   12 CONTINUE
      DO 14 I=1,NX
      X=XGRID(I)
      CALL FINTGS(X,FS,FG,FS1(I))
      CALL FINTGG(X,FG,FS,FG1(I))
   14 CONTINUE      
      S=NTIMES*DELTA
 10   CONTINUE
C**  START THE NONSINGLET CALCULATION
C**  CALCULATE FVL AT S=0.
      NVL=9
      DO 100 KIVL=1,NVL
      IVL=KIVL+2
C
C  IPM IS A PHASE THAT APPEARS IN THE
C  NEXT-TO-LEADING ORDER CALCULATION
C  IPM=+1 FOR THE 'PLUS' TYPE DISTRIBUTIONS DEFINED AS
C         (Q+QBAR)-SINGLET/FLAVORS
C  IPM=-1 FOR THE 'MINUS' TYPE DISTRIBUTIONS DEFINESD AS 
C          (Q-QBAR)-Q_NS^V/FLAVORS
C
      IPM=1
      IF(IVL.GE.7) IPM=-1
      S=SV
      DO 41 I=1,NX
      X=XGRID(I)
      IF(JFL.EQ.1)THEN
         FVL(I)=FCNVL(IVL,X)
         CALC(IVL,1,I)=FVL(I)
c         if(ivl.eq.7) tmp(5,i)=fvl(i)
      ELSE
         IF(IVL.EQ.7)THEN
            CALC(IVL,NCB+2,I)=CALC(IVL,NCB+1,I)
            FVL(I)=CALC(IVL,NCB+2,I)
         ELSE IF(IVL.GT.2.AND.IVL.LT.7)THEN
c            sing=calc(2,ncb+1,i)
c            singm=calc(7,ncb+1,i)
c            if(sing.lt.singm)sing=singm
c            CALC(IVL,NCB+2,I)=CALC(IVL,NCB+1,I)+sing/20.
            CALC(IVL,NCB+2,I)=CALC(IVL,NCB+1,I)+CALC(2,NCB+1,I)/20.
            FVL(I)=CALC(IVL,NCB+2,I)
         ELSE IF(IVL.GT.7)THEN
            CALC(IVL,NCB+2,I)=CALC(IVL,NCB+1,I)+CALC(7,NCB+1,I)/20.
            FVL(I)=CALC(IVL,NCB+2,I)
         ENDIF
c         if(ivl.eq.7) tmp2(5,i)=fvl(i)
      ENDIF
   41 CONTINUE
C**  CALCULATE D(FVL)/DS AT S=SV
      T=SST
      DO 42 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,FVL,FVL1(I))
   42 CONTINUE
C**  CALCULATE PREDICTOR AT S=SV+DELTA
      DO 43 I=1,NX
      PRD(I)=FVL(I)+DELTA*FVL1(I) 
   43 CONTINUE
C**  CALCULATE D(PRD)/DS AT S=SV+DELTA
      T=SST*DEXP(DELTA)
      DO 44 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,PRD,PRD1(I))
      COR(I)=FVL(I)+0.5*DELTA*(FVL1(I)+PRD1(I))
c      if(ivl.eq.7.and.jfl.eq.1)then
c         tmp(6,i)=fvl1(i)
c         tmp(7,i)=prd1(i)
c         tmp(8,i)=cor(i)
c      else if(ivl.eq.7.and.jfl.eq.2)then
c         tmp2(6,i)=fvl1(i)
c         tmp2(7,i)=prd1(i)
c         tmp2(8,i)=cor(i)
c      endif
      FUDGE(I)=PRD(I)-COR(I)
   44 CONTINUE
      DO 45 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,COR,COR1(I))
   45 CONTINUE
      DO 46 I=1,NX
      FVLSTO(I)=FVL(I)
      FVL(I)=COR(I) 
      FVL1(I)=COR1(I)
      CALC(IVL,NN,I)=FVL(I)
   46 CONTINUE
C**  INITIALIZATION COMPLETE...ITERATE UNTIL S=SMAX
      DO 47 NTIMES=N1,N2
      NT1=NTIMES+1
      DO 51 I=1,NX
      X=XGRID(I)
      PRD(I)=FVLSTO(I)+2.*DELTA*FVL1(I) 
      PRDM(I)=PRD(I)-0.8*FUDGE(I)
   51 CONTINUE
      T=T*DEXP(DELTA)
      DO 52 I=1,NX
      X=XGRID(I)
      FVLSTO(I)=FVL(I)
      CALL FINTGV(X,IPM,PRDM,PRD1(I))
      COR(I)=FVL(I)+0.5*DELTA*(FVL1(I)+PRD1(I))
      FUDGE(I)=PRD(I)-COR(I)
      FVL(I)=COR(I)+0.2*FUDGE(I)
   52 CONTINUE
      DO 53 I=1,NX
      X=XGRID(I)
      CALL FINTGV(X,IPM,FVL,FVL1(I))
      CALC(IVL,NT1,I)=FVL(I)
   53 CONTINUE
      S=NTIMES*DELTA
   47 CONTINUE
  100 CONTINUE
 70   CONTINUE
c      open(unit=4,file='cmp_sns.dat',status='unknown')
c      do 71 j=1,60
c      write(4,72),xgrid(j),(tmp(i,j),i=1,8)
c      write(4,72),xgrid(j),(tmp2(i,j),i=1,8)
c 71   continue
c 72   format(9e16.9)
      DO 49 J=1,100
   49 XPAR(J)=PARA(J)
c      open(unit=12,file='tst10_calc.dat',status='unknown')
c      do 5001 k=1,60
c      do 5001 j=1,60
c      do 5001 i=1,11
c      write(12,5002) calc
c 5001 continue
 5002 format(8e15.9)
      RETURN
      END 



      SUBROUTINE FSUPDF(NDRV,X,S,U,D,UB,DB,SB,CB,BB,GLUE)
      IMPLICIT REAL*8 (A-H,O-Z)
      logical onoff
c      double precision dfv,dfsq,dfg,dfsb,dfcb,dfbb ! for deuteron off-shell 
c                                                    calculations
c      common/input/ins,nvl,nmax,delta,s0,iord
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      common/threshold/sbbar
c
c  modified 1/8/10 for the new evolution scheme
c  note: q not equal to qbar not yet implemented, but the 
c  arrays are there if we want to go to NNLO at some future 
c  date -- jfo

c
      CALL GINTERP(1,NDRV,X,S,glue)
      CALL GINTERP(2,NDRV,X,S,sing)
      CALL GINTERP(3,NDRV,X,S,up)
      CALL GINTERP(4,NDRV,X,S,dp)
      CALL GINTERP(5,NDRV,X,S,sp)
      call ginterp(6,ndrv,x,s,cp)
      CALL GINTERP(7,NDRV,X,S,singm)
      CALL GINTERP(8,NDRV,X,S,um)
      CALL GINTERP(9,NDRV,X,S,dm)
      nf=5
      if(s.lt.sbbar)nf=4
      uv=um+singm/nf
      dv=dm+singm/nf
      upub=up+sing/nf
      dpdb=dp+sing/nf
      ub=.5*(upub-uv)
      db=.5*(dpdb-dv)       
      sb=.5*(sp+sing/nf)
      cb=.5*(cp+sing/nf)
      if(nf.eq.5)then
         bp=-(up+dp+sp+cp)
         bb=.5*(bp+sing/nf)
      else
         bb=0.
      endif
C     Now gets u,d out of uv,dv and checks for negative PDFs
      if(x.lt.0.80d0)then
         u=uv+ub
         d=dv+db
      else
         u=uv
         d=dv
      endif
      if(ub.lt.0.d0) ub=0.d0
      if(db.lt.0.d0) db=0.d0
      if(sb.lt.0.d0) sb=0.d0
      if(cb.lt.0.d0) cb=0.d0
      if(bb.lt.0.d0) bb=0.d0

C!!!!!!! The following code relates to the deuteron off-shell corrections
C        If you are just computing the structure functions, all the
C        following variables will be one.
C        Only if you are calculating an off-shell part of the structure
C        function will these not be one.
c         call get_delta_v(x,dfv) ! valence correction
c         call get_delta_sq(x,dfsq) ! ubar and dbar correction
c         call get_delta_g(x,dfg) ! gluon correction
c         call get_delta_sb(x,dfsb) ! strange correction - not implemented
c         call get_delta_cb(x,dfcb) ! charm correction - not implemented
c         call get_delta_bb(x,dfbb) ! bottom correction - not implemented
c
c Modified 3/24/15 to all the off_shell strength to be fitted
c Combines all 'get_delta' calls into one routine
c
      call get_offshell_on(onoff)
      if(onoff)then
         call get_delta(ndrv,x,dfv,dfsq,dfg,dfsb,dfcb,dfbb)
	 write(*,*) "Let's use off-shell at Quark level"
	 ! write(*,*) dfv
         u = (u-ub)*dfv + ub*dfsq
         d = (d-db)*dfv + db*dfsq
         ub = ub*dfsq
         db = db*dfsq
         sb = sb*dfsb
         cb = cb*dfcb
         bb = bb*dfbb
         glue = glue*dfg
      endif   
C        end of off-shell code
         !write(*,*) "***** This is printing because onoff condition is false but get_offshell_on is calling"

      RETURN
      END

      SUBROUTINE GINTERP(I,NDRV,X,S,ANS)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F1(60),F2(60)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      COMMON/MINIMB/UNCRT(100),ERMIN(50),PWATE(100),
     *     IFREP(50),NPARAM,NVAR,IREJ
      COMMON/GINT/GF(11,60,2160)
      common/threshold/sb
      common/stepsize/delta

C
C  THIS ROUTINE INTERPOLATES AS NEEDED IN THE ARRAY GF TO OBTAIN THE
C  EVOLVED DISTRIBUTIONS AT THE REQUIRED Q2 AND X VALUES. NOTE THAT ONLY
C  50 PARAMETERS CAN BE VARIED IN A GIVEN FIT. OF THESE, AT MOST 35 CAN 
C  BE PARAMETERS RELATED TO THE PARTON DISTRIBUTIONS. THIS IS GOVERNED BY THE
C  THIRD DIMENSION OF GF (2160 ABOVE). THIS IS 60+60*(NO. OF VARIED PARTON 
C  PARAMETERS).
C
c
c  modified 1/8/10 for new evolution scheme and increased array
c  dimensions -- jfo
c
      if(s.lt.sb)then         
         is=s/delta+1
         s1=(is-1)*delta
         s2=s1+delta
      else
c         is=(s-sb)/delta+14
         is=(s-sb)/delta+ncb+2
c         s1=(is-14)*delta+sb
         s1=(is-(ncb+2))*delta+sb
         s2=s1+delta
      endif
      IS1=IS+1
      NDRV2=NDRV                    ! ndrv2 = PDF direction in par space
      IF(IFREP(NDRV).GT.35) NDRV2=0 ! but for non-PDF directions, sets ndrv2=0
      DO 1 L=1,60
      KL=L+60*NDRV2 ! Note: uses the central grid for non-PDF directions
      F1(L)=GF(I,IS,KL)
      F2(L)=GF(I,IS1,KL)
    1 CONTINUE
      A1=GETFV(X,F1)
      A2=GETFV(X,F2)
      ANS=A1*(S-S2)/(S1-S2)+A2*(S-S1)/(S2-S1)
      RETURN
      END 




      subroutine pdfnrm
      implicit real*8 (a-h,o-z)
c      common/gaus16/xi(16),wi(16),nterms,xx(17)
      common/gaus32/xi(32),wi(32),nterms,xx(33)
      common/param/xp(100)
c
c  quantum number sum rule for uv normalization
c
      auv=0.
      do 100 l=1,nterms
         z=.5*(1.+xi(l))
         x=z**3
c        # integral of uv from x= 0-1 without overall coefficient
         auv=auv+.5*wi(l)*fuv(x)/x/xp(2)*3.*z**2
 100  continue
c     # renormalize uv to satisfy sum rules
      xp(2)=2./auv

c
c  quantum number sum rule for dv normalization
c
      adv=0.
      do 110 l=1,nterms
         z=.5*(1.+xi(l))
         x=z**3
c        # integral of dv from x= 0-1 without overall coefficient
         adv=adv+.5*wi(l)*fdv(x)/x/xp(8)*3.*z**2
 110  continue
c     # renormalize dv to satisfy sum rules
      xp(8)=1./adv

c
c  momentum sum rule for gluon normalization
c
      as=0.
      ag=0.
      do 200 l=1,nterms
      z=.5*(1.+xi(l))
      x=z**3
c     # integral of x(singlet)
      as=as+.5*wi(l)*fcnfs(x)*3.*z**2
c     # integral of xgluon without overall coefficient
      ag=ag+.5*wi(l)*fcnfg(x)/xp(25)*3.*z**2
 200  continue
c     # renormalize gluon distribution
      xp(25)=(1.-as)/ag
      return
      end


      FUNCTION FCNVL(IVL,X)
      IMPLICIT REAL*8 (A-H,O-Z) 
C  CALCULATES THE INPUT QUARK DISTRIBUTIONS AT Q0**2
C  THE SINGLET AND GLUON DISTRIBUTIONS ARE HANDLED IN
C  SEPARATE ROUTINES (FCNFS AND FCNFG)
      COMMON/PARAM/XP(100)
      COMMON/GRID/NX,XGRID(60)
      COMMON/GRPTHY/FLAVOR
      common/iparam/iparamc
C
c  IVL = 3  XU+
C  IVL = 4  XD+
C  IVL = 5  XS+
C  IVL = 6  XC+
C  IVL = 7  XQNSV
C  IVL = 8  XU-
C  IVL = 9  XD-
C  IVL = 10 XS-
C  IVL = 11 XC-
C
C  Revised 6/3/02 for five flavors
c  Parametrization changed to CTEQ6 style
c  Revised 2/27/03 to implement improved treatment of b threshold
c
C  REVISED 6/17/05 TO INLCUDE THE 3-LOOP SPLITTING FUNCTIONS
C  ARRAY ELEMENTS RENUMBERED
C
C      if(ivl.eq.1)then
c
c  xuv
c
C         fcnvl=fuv(x)
C         return
C      else if(ivl.eq.2)then
c
c  xdv
c
C         fcnvl=fdv(x)
C         return

      if(ivl.eq.3)then
c
c  x(uv+2.*ubar-.25*singlet)
c
         ubpdb=fubpdb(x)
         dboub=fdboub(x)
         fcnvl=fuv(x)+2.*ubpdb/(1.+dboub)-.25*fcnfs(x)
         return

      else if(ivl.eq.4)then
c
c  x(dv+2.*dbar-.25*singlet)
c
         ubpdb=fubpdb(x)
         dboub=fdboub(x)
         fcnvl=fdv(x)+2.*ubpdb/(1.+1./dboub)-.25*fcnfs(x)
         return

      else if(ivl.eq.5)then
c
c  x(s+sbar-.25*singlet)
c
c  iapramc<=7 used in CTEQ-JLab fits, and earlier fits by Jeff
c  iparamc=8 uses parameters 7 and 13 to fit for x-dependence in kappa
c  (kappa = xp(31)) -- jfo 11/26/16 
c  Extended also params 19, 30 -- aa 2/13/2017
c  iparamc=9 parametrizes s=s(x) directly instead of s=k(x)*(ybar+dbar)
c  -- aa 2/13/2017
         if(iabs(iparamc).le.7
     &        .or.iabs(iparamc).eq.10
     &        .or.iabs(iparamc).eq.11) then
            ! x(s+sbar) = kappa * (ubar+dbar)
            fcnvl=xp(31)*fubpdb(x)-.25*fcnfs(x)
            !print*, '* fcnvl - 10'
         else if(iabs(iparamc).eq.8)then
            ! x(s+sbar) = kappa(x) * (ubar+dbar)
            fcnvl=xp(31)*x**xp(7)*(1.-x)**xp(13)
     &           *(1d0+xp(19)*dsqrt(x)+xp(30)*x)*fubpdb(x)-.25*fcnfs(x)
         else if(iabs(iparamc).eq.9)then
            ! x(s+sbar) = 2*s(x) -- directly parametrized 
            fcnvl=xp(31)*x**xp(7)*(1.-x)**xp(13)
     &           *(1d0+xp(19)*dsqrt(x)+xp(30)*x) - .25*fcnfs(x)
         else
            print*, 'ERROR (altpar15): strange iparamc out of range:'
     *           , iparamc
         endif
         !print*, '* fcnvl - iparamc =', iparamc
         return
         
      else if(ivl.eq.6)then
c
c  xcplus = -1/4*singlet
c
         fcnvl=-.25*fcnfs(x)
         RETURN

      ELSE IF(IVL.EQ.7)THEN
C
C  XQNSV
C  X(U-UBAR + D-DBAR + S-SBAR + C-CBAR)
C  ASSUME S=SBAR AND C=CBAR FOR NOW
C
         FCNVL = XQNSV(X)
         RETURN

      ELSE IF(IVL.EQ.8)THEN
C
C  XU- = X(U-UBAR) - .25*XQNSV
C
         FCNVL=FUV(X)-.25*XQNSV(X)
         RETURN

      ELSE IF(IVL.EQ.9)THEN
C
C  XD- = X(D-DBAR)-.25*XQNSV
C
         FCNVL=FDV(X)-.25*XQNSV(X)
         RETURN

      ELSE IF(IVL.EQ.10)THEN
C
C  XS- = X(S-SBAR) - .25*XQNSV
C 
         FCNVL=-.25*XQNSV(X)
         RETURN

      ELSE IF(IVL.EQ.11)THEN
C
C  XC- = X(C-CBAR) - .25*XQNSV
C 
         FCNVL=-.25*XQNSV(X)
         RETURN

      ENDIF
      end


      function fuv(x)
*     Up valence parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
c
c     Added parametrization choices 5/17/11 jfo
c
      if(iabs(iparamc).eq.1)then
         fuv=xp(2)*x**xp(3)*(1.-x)**xp(4)*exp(xp(5)*x)
     $      *(1.+exp(xp(6))*x)**xp(7)
      else if(iabs(iparamc).eq.2)then
         fuv=xp(2)*x**xp(3)*(1.-x)**xp(4)*(1.+xp(5)*x+xp(6)*x**2)
      else if(iabs(iparamc).ge.3.and.iabs(iparamc).le.7) then 
         fuv=xp(2)*x**xp(3)*(1.-x)**xp(4)*(1.+xp(5)*sqrt(x)+xp(6)*x
     &       +xp(7)*x*sqrt(x))
c
c iparamc=8,9 frees up xp(7,13,19,30) for use in strange sea studies
c -- jfo 11/26/16  /  aa 2/13/17
c iparamc=10,11 frees up xp(7,13,19,30) for use in dbar/ubar studies
c -- aa April 2017
c
      else if(iabs(iparamc).ge.8.and.iabs(iparamc).le.11) then
         fuv=xp(2)*x**xp(3)*(1.-x)**xp(4)*(1.+xp(5)*sqrt(x)+xp(6)*x)
         !print*, '* fuv - 10'
      endif  
      !print*, '* fuv - iparamc =', iparamc
      return
      end


      function fdv(x)
*     Down valence -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
* (02 Jul 09) d-valence for d/u->finite as x->1. 
*     New functional form: d'=d+c*x^alpha*u
*     with 2 free parameter to allow any  lim_(x->1) d/u. 
*     NOTE: alpha=1.4 --> reproduces Bodek-Yang, PRL82(99)2467
*           alpha=2.3 --> reproduces Melnitchouk-Peng, PLB400(97)220
*           c=0 --> standard CTEQ6.1 d-quark
*     (A.Accardi 2-9 Jul 2009 - bug corrected 17 May 2010)
*     (AA Oct-Dec 2016) added "iparamc<0" option to remove absolute
*           values from the parametrization
c
c     Added parametrization choices 5/17/11 jfo
c

C     ### AA Dec 2016:
C     ### iparamc<0 --> don't use the absolute values for #32 and #33
C     ### (in the future may want to use squares insted of abs values)
      if (iparamc.gt.0) then    ! params as originally designed
         const = dabs(xp(32))
         alpha = dabs(xp(33))
      else  ! remove the absolute values
         const = xp(32)
         alpha = xp(33)
      end if
         
      if(iabs(iparamc).eq.1)then
         fdv=xp(8)*(x**xp(9)*(1.-x)**xp(10)*exp(xp(11)*x)
     &      *(1.+exp(xp(12))*x)**xp(13)
     &      + const*x**alpha*fuv(x))
      else if(iabs(iparamc).eq.2)then
         fdv=xp(8)*(x**xp(9)*(1.-x)**xp(10)*(1.+xp(11)*x+xp(12)*x**2)
     &      + const*x**alpha*fuv(x))
      else if(iabs(iparamc).ge.3.and.iabs(iparamc).le.7) then
         fdv=xp(8)*(x**xp(9)*(1.-x)**xp(10)*(1.+xp(11)*sqrt(x)+xp(12)*x
     &      + xp(13)*x*sqrt(x))
     &      + const*x**alpha*fuv(x))
c
c iparamc=8,9 frees up xp(7,13,19,30) for use in strange sea studies
c -- jfo 11/26/16  /  aa 2/13/17
c iparamc=10,11 frees up xp(7,13,19,30) for use in dbar/ubar studies
c -- aa April 2017
c
      else if(iabs(iparamc).ge.8.or.iabs(iparamc).le.11)then
         fdv=xp(8)*(x**xp(9)*(1.-x)**xp(10)*(1.+xp(11)*sqrt(x)+xp(12)*x)
     &        + const*x**alpha*fuv(x))
         !print*, '* fdv - 10'
      endif
      !print*, '* fdv - iparamc =', iparamc
      return
      end

      function fubpdb(x)
*     ubar+dbar -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
c
c     Added parametrization choices 5/17/11 jfo
c
      if(iabs(iparamc).eq.1)then
         fubpdb=xp(14)*x**xp(15)*(1.-x)**xp(16)*exp(xp(17)*x)*
     $         (1.+exp(xp(18))*x)**xp(19)
      else if(iabs(iparamc).eq.2)then
         fubpdb=xp(14)*x**xp(15)*(1.-x)**xp(16)*
     $         (1.+xp(17)*x+xp(18)*x**2)
      else if(iabs(iparamc).ge.3)then
         fubpdb=xp(14)*x**xp(15)*(1.-x)**xp(16)*
     $         (1.+xp(17)*sqrt(x)+xp(18)*x)
         !print*, '* fubpdb - 10'
      endif
      !print*, '* fubpdb - iparamc =', iparamc
      return
      end


      function fdboub(x)
c
c     Added parametrization choices 5/17/11 jfo
c
*     dbar/ubar -  parametrization at Q0
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
      if(iabs(iparamc).eq.1.or.iabs(iparamc).eq.6)then
         fdboub=xp(20)*x**(xp(21))*(1.-x)**xp(22)
     $   +(1.+xp(23)*x)*(1.-x)**xp(24)
c     $   +1.+xp(23)*x*(1.-x)**xp(24)

c
c  add piece from CTEQ parametrization (gfun.f)
c
         if(iabs(iparamc).eq.6)return
         bconst = 10.d0
         fac=fdboub
         if(fac.gt.bconst)then
            fdboub=fac
         else if(fdboub.lt.-bconst)then
            fdboub=0.
            return
         else
            tmp = 1.d0 + exp(-bconst*fac) - exp(-bconst)
            fdboub = fdboub + (1.d0/bconst)*log(tmp)
         endif
      else if (iabs(iparamc).ge.7.and.iabs(iparamc).le.9)then
         ! CJ15 parametrization
         fdboub=xp(20)*x**(xp(21))*(1.-x)**xp(22)
c     $   +(1.+xp(23)*x)*(1.-x)**xp(24)
     $   +1.+xp(23)*x*(1.-x)**xp(24)
      else if (iabs(iparamc).eq.10) then
c        ! extended CJ15 - borrows params 7,13,30 from other PDFs
c        ! [AA April 2017]
         fdboub = xp(20)*x**(xp(21))*(1.-x)**xp(22) ! primary bump
     $        +1. + xp(7)*x**xp(13) ! Asymptotic at x-->0,1
     $        +xp(23)*x**xp(30)*(1.-x)**xp(24) ! secondary bump
         print*, '* fdboub - xp(7,13,30)=',xp(7),xp(13),xp(30)
      else if (iabs(iparamc).eq.11) then
c        ! "regular" db/ub parametrization - borrows params 13,30
c        ! from other PDFs [AA April 2017]
         fdboub = xp(20)*x**(xp(21))*(1.-x)**xp(22)   ! bump
     $        * (1. + xp(23)*dsqrt(x) + xp(24)*x )    ! deformation
     $        + 1. + xp(7)*x**xp(13)          ! Asymptotic at x-->0,1
      else if(iabs(iparamc).gt.1.and.iabs(iparamc).le.5)then
         ! CJ12 and earlier (db-ub parametrizaed instead, see "fdbmub")
c         fdboub=1.+xp(20)*x**xp(21)*(1.-x)**xp(22)
c     $            +xp(23)*x*(1.-x)**xp(24)
         fdboub=(fubpdb(x)+fdbmub(x))/(fubpdb(x)-fdbmub(x))
c      else if(iabs(iparamc).eq.5)then
c         fdboub=1.+xp(20)**2*x**xp(21)*(1.-x)**xp(22)*exp(xp(23)*x)
c     2*(1.+exp(xp(24))*x)
c         fdboub=(1.+xp(20)**2*x**xp(21))*(1.-x)**xp(22)*exp(xp(23)*x)
c     2*(1.+sqrt(xp(24)**2)*x)
      else
         print*, 'ERROR(altpar_fdboub): iparamc out of range: ',iparamc
      endif
      !print*, '* fdboub - iparamc =',iparamc
      return
      end

      function fdbmub(x)
*
*     x(dbar-ubar) -  parametrization at Q0
*
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
c  only used for iparamc=2 or 3 or 4 or 5
      if(iabs(iparamc).eq.2.or.iabs(iparamc).eq.3)then
         fdbmub=xp(20)*x**xp(21)*(1.-x)**xp(22)*(1.+xp(23)*x
     2         +xp(24)*x**2)
      else if(iabs(iparamc).eq.4)then
c         xp(22)=xp(16)+2.5d0
         fdbmub=xp(20)*x**xp(21)*(1.-x)**xp(22)*(1.+xp(23)*sqrt(x)
     2         +xp(24)*x)
      else if(iabs(iparamc).eq.5)then
c         print*,'error with iparamc=',iparamc
c         stop
         fdbmub=xp(20)*x**xp(21)*(1.-x)**xp(22)*exp(xp(23)*x)
     2         *(1.+sqrt(xp(24)**2)*x)
      else if(iabs(iparamc).ge.6)then
         print*,'ERROR(fdbmub): iparamc out of bound ',iparamc
         stop
      endif
      return
      end
c

      function fcnfs(x)
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
c     
c  singlet at Q0
c
      if (iabs(iparamc).le.7
     &     .or.(iabs(iparamc).ge.10.and.iabs(iparamc).le.11)) then 
         fcnfs=fuv(x)+fdv(x)+fubpdb(x)*(2.+xp(31))
         !print*, '* fcnfs - 10'
      else if (iabs(iparamc).eq.8) then
         fcnfs=fuv(x)+fdv(x)+fubpdb(x)
     &        * ( 2. + xp(31)*x**xp(7)*(1.-x)**xp(13)
     &           *(1d0+xp(19)*dsqrt(x)+xp(30)*x) )
      else if (iabs(iparamc).eq.9) then
         fcnfs=fuv(x)+fdv(x)+2.*fubpdb(x)
     &        + 2*xp(31)*x**xp(7)*(1.-x)**xp(13)
     &          *(1d0+xp(19)*dsqrt(x)+xp(30)*x)
      else
         print*,'ERROR(fcnfs): iparamc out of bound ',iparamc
         stop
      endif
      !print*, '* fcnfs - iparamc =',iparamc
      return
      end

      function fcnfg(x)
      implicit real*8 (a-h,o-z)
      common/param/xp(100)
      common/iparam/iparamc
c
c  gluon at Q0
c
c
c     Added parametrization choices 5/17/11 jfo
c
      if(iabs(iparamc).eq.1)then
         fcnfg=xp(25)*x**xp(26)*(1.-x)**xp(27)*exp(xp(28)*x)*
     $   (1.+exp(xp(29))*x)**xp(30)
      else if(iabs(iparamc).eq.2)then
         fcnfg=xp(25)*x**xp(26)*(1.-x)**xp(27)*(1.+xp(28)*x+xp(29)*x**2)
      else if(iabs(iparamc).ge.3)then
         fcnfg=xp(25)*x**xp(26)*(1.-x)**xp(27)*(1.+xp(28)*sqrt(x)
     *        +xp(29)*x)
         !print*, '* fcnfg - 10'
      endif
      !print*, '* fcnfg - iparamc =',iparamc
c      fcnfg=xp(25)*(x**xp(26)*(1.-x)**xp(27)*(1.-xp(28)*sqrt(x)+xp(29)*x)
c     2+xp(30)*x**xp(26)*(1.-x)**32.)
      return
      end

      FUNCTION XQNSV(X)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  X(U-UBAR + D-DBAR + S-SBAR + C-CBAR)
C  ASSUME S=SBAR AND C=CBAR FOR NOW
C
      XQNSV=FUV(X)+FDV(X)
      RETURN
      END

      SUBROUTINE FINTGV(X,IPM,FVL,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES VALENCE TERM(S) - SEE NOTES
      DIMENSION FVL(60)
      COMMON/LAGUER/XL(8),WL(8),NTERML
c      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(100)
c      AL=ALPHA(T)/(2.*PI)
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(ILOOP,xmu2,xc(1),neff)/(2.d0*pi)
      flavor=neff
      FX=GETFV(X,FVL)
      FINTG=T*AL*FX*(2.d0+8.d0/3.d0*DLOG(1.d0-X))
c      IF(IORD.GT.1)THEN
      IF(IORD.GT.0)THEN
         C1=16.d0/9.d0
         C2=2.d0
         C3=2.d0/3.d0*FLAVOR
         C4=1.202056903D0
         C5=-5.d0/8.d0*C4
         A1=C2*2.d0*(67.d0/9.d0-PI2/3.d0)-C3*20.d0/9.d0
         C=C1*(3.d0/8.d0-PI2/2.d0+C4-8.d0*C5)
     2   +C2*(17.d0/12.d0+11.d0*PI2/9.d0
     2   -C4+8.d0*C5)-C3*(1.d0/6.d0+2.d0*PI2/9.d0)
         IF(IORD.EQ.2)THEN
            IF(IPM.EQ.1)THEN
               C=C+AL/8.d0*P2NSPC(X,NEFF)
            ELSE IF(IPM.EQ.-1)THEN
               C=C+AL/8.d0*P2NSMC(X,NEFF)
            ENDIF
         ENDIF
         FINTG=FINTG+T*AL**2*FX*(C+A1*DLOG(1.d0-X))
      ENDIF
      DO 1 I=1,NTERMS
      Z=0.5d0*(1.d0-X)*XI(I)+0.5d0*(1.d0+X)
      FZ=GETFV(X/Z,FVL)
      FINTG=FINTG+0.5d0*(1.-X)*T*AL*4.d0/3.d0*((1.d0+Z*Z)*FZ
     2  -2.d0*FX)/(1.d0-Z)*WI(I)
      IF (IORD.GT.0)THEN
         ALZ=DLOG(Z)
         AL1=DLOG(1.d0-Z)
         Z2=1.d0+Z*Z
         ZM=1.d0-Z
         ZP=1.d0+Z
         AZ=C1*(-2.d0*Z2*ALZ*AL1-3.d0*ALZ)+C2*Z2*(ALZ*(ALZ+11.d0/3.d0)
     2   +67.d0/9.d0-PI2/3.d0)-C3*2.d0/3.d0*Z2*(ALZ+5.d0/3.d0)
         PA=4.d0*ZM+2.d0*ZP*ALZ+2.d0*Z2/ZP*S2(Z)
         BZ=C1*(-ALZ*(ALZ*ZP/2.d0+2.d0*Z)-5.d0*ZM+IPM*PA)
     2   +C2*(2.d0*ZP*ALZ
     2   +40.d0/3.d0*ZM-IPM*PA)-C3*4.d0/3.d0*ZM
         FINTG=FINTG+0.5d0*(1.-X)*T*AL**2*((AZ*FZ-A1*FX)/ZM+BZ*FZ)*WI(I)
         IF(IORD.EQ.2)THEN
            IF(IVL.EQ.7)THEN
               PFF2=P2NSMA(Z,NEFF)+P2NSSA(Z,NEFF)
               PFF3=P2NSB(Z,NEFF)
            ELSE IF(IPM.EQ.1)THEN
               PFF2=P2NSPA(Z,NEFF)
               PFF3=P2NSB(Z,NEFF)
            ELSE IF(IPM.EQ.-1)THEN
               PFF2=P2NSMA(Z,NEFF)
               PFF3=P2NSB(Z,NEFF)
            ENDIF
            FINTG=FINTG+.5d0*(1.-X)*T*AL**3/8.d0*(FZ*PFF2+PFF3*(FZ-FX))
     2      *WI(I)
         ENDIF
      ENDIF
    1 CONTINUE
      RETURN
      END 

      FUNCTION S2(X)
      IMPLICIT REAL*8 (A-H,O-Z)
c      DIMENSION FN(51)
c      DATA FN/1.644934076,1.588625448,1.545799712,1.507899041,
c     X1.473125860,1.440633797,1.409928300,1.380685041,1.352675161,
c     X1.325728728,1.299714723,1.274529160,1.250087584,1.226320101,
c     X1.203167961,1.180581124,1.158516487,1.136936560,1.115808451,
c     X1.095103088,1.074794600,1.054859830,1.035277934,1.016030062,
c     X0.997099088,0.978469393,0.960126675,0.942057798,0.924250654,
c     X0.906694053,0.889377624,0.872291733,0.855427404,0.838776261,
c     X0.822330471,0.806082689,0.790026024,0.774153992,0.758460483,
c     X0.742939737,0.727586308,0.712395042,0.697361058,0.682479725,
c     X0.667746644,0.653157631,0.638708705,0.624396071,0.610216108,
c     X0.596165361,0.582240526/
      COMMON/CONSTANTS/PI,PI2
C
C  THESE ARE THE VALUES OF F(X)=LI2(1-X) FOR X BETWEEN O AND .50 IN STEPS
C  OF O.O1 TAKEN FROM ABRAMOWITZ AND STEGUN, PG. 1005, TABLE 27.7.
C
C  S2(X)=INTEGRAL OF LN((1-Z)/Z)/Z FROM X/(1+X) TO 1/(1+X)
C
C  REWRITE S2(X) IN TERMS OF F(X/(1+X))
C
C  USE A LINEAR INTERPOLATION TO OBTAIN F(Z), Z=X/(1.+X)
C  IN THE REGION OF X BELOW 0.8
C
C  NEAR X=1 SMALL ERRORS IN THE INTERPOLATION ARE AMPLIFIED 
C  SINCE THE ANSWER IS A SMALL DIFFERENCE BETWEEN MUCH LARGER 
C  NUMBERS. USE A TAYLOR SERIES FOR F(X) EXPANDED ABOUT THE 
C  NEAREST X VALUE.
C
      Z=X/(1.+X)
c      N=100.*Z+1
c      Z1=(N-1)/100.
c      Z2=N/100.
c      IF(X.LT.0.8) THEN
c         F=FN(N)*(Z-Z2)/(Z1-Z2)+FN(N+1)*(Z-Z1)/(Z2-Z1)
c      ELSE
c         DELT1=DABS(Z-Z1)
c         DELT2=DABS(Z-Z2)
c         ZT=Z1
c         IF(DELT2.LT.DELT1) THEN
c            ZT=Z2
c            N=N+1
c         ENDIF
c         OMZT=1.-ZT
c         F0=FN(N)
c         F1=DLOG(ZT)/OMZT
c         F2=(1./ZT+F1)/OMZT
c         F3=(1./ZT**2+2.*F2)/OMZT
c         F=F0+(Z-ZT)*F1+.5*(Z-ZT)**2*F2+(Z-ZT)**3/6.*F3
c      ENDIF
c      S2=PI**2/6.+.5*DLOG(X)**2-DLOG(1.+X)**2-2.*F
      S2=PI2/6.+.5*DLOG(X)**2-DLOG(1.+X)**2-2.*ddilog(1./(1.+x))
      RETURN
      END

      SUBROUTINE FINTGG(X,FG,FS,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES GLUON TERM - SEE NOTES
      DIMENSION FG(60),FS(60) 
c      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(100)
c      AL=ALPHA(T)/2./PI
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(ILOOP,xmu2,xc(1),neff)/(2.*pi)
      flavor=neff
      FX=GETFV(X,FG)
      FINTG=T*AL*FX*((33.-2.*FLAVOR)/6.+6.*DLOG(1.-X))
      CF=4./3.
      IF(IORD.GT.0)THEN
         TR=FLAVOR/2.
         CA=3.
         Z3=1.202056903D0
         AL1=DLOG(1.-X)
         PGG1=-CF*TR-4./3.*CA*TR*(1.+5./3.*AL1)+CA*CA/3.*
     2   (8.+9.*Z3+(67/3.-PI2)*AL1)
         IF(IORD.EQ.2)THEN
            PGG1=PGG1+AL/8.*P2GGC(X,NEFF)
         ENDIF
         FINTG=FINTG+T*AL*AL*FX*PGG1
      ENDIF
      DO 1 I=1,NTERMS
      Z=0.5*(1.-X)*XI(I)+0.5*(1.+X)
      FGZ=GETFV(X/Z,FG)
      FSZ=GETFV(X/Z,FS)
      FINTG=FINTG+.5*(1.-X)*T*AL*WI(I)*(CF*FSZ*(Z*Z-2.*Z+2.)/Z+
     26.*(FGZ*Z-FX)/(1.-Z)+6.*FGZ*(1.-Z+Z*Z-Z**3)/Z)
      IF(IORD.GT.0)THEN
         AL1=DLOG(1.-Z)
         ALZ=DLOG(Z)
C     *  1         2         3         4         5         6         7 *
         PGG2=CF*TR*(4./3.*(1./Z-12.+6.*Z+5.*Z*Z)-2.*(3.+5.*Z)*ALZ
     2   -2.*(1.+Z)*ALZ*ALZ)+4./3.*CA*TR*(-1./6.*(23./Z-29.+19.*Z
     3   -23.*Z*Z)-(1.+Z)*ALZ)+CA*CA*(-1./18.*(25.+109.*Z)
     4   -1./3.*(25.-11.*Z+44.*Z*Z)*ALZ+4.*(1.+Z)*ALZ*ALZ
     5   -PI2/3.*(1./Z-2.+Z-Z*Z)+(1./Z-2.+Z-Z*Z+
     6   1./(1.-Z))*ALZ*(ALZ-4.*AL1)-2.*(1./Z+2.+Z+Z*Z-1./(1.+Z))*S2(Z))
         PGG3=(-20./9.*CA*TR+CA*CA/3.*(67./3.-PI2))/(1.-Z)
         PGQ=4./3.*CF*TR*(2./3.*(5./Z-5.+4.*Z)+(2./Z-2.+Z)*AL1)
     2   +.5*CF*CF*(5.+7.*Z-(4.+7.*Z)*ALZ+2.*(6./Z-6.+5*Z)*AL1 
     3   +(2.-Z)*ALZ**2+2.*(2./Z-2+Z)*AL1**2)+CA*CF*(-1./9.*(9./Z+19.
     4   +37.*Z+44.*Z*Z)+1./3.*(36.+15.*Z+8.*Z*Z)*ALZ-.5*(2./Z+6.
     5   +3.*Z)*ALZ**2-1./3.*(22./Z-22.+17.*Z)*AL1+(2./Z-2.+Z)*(
     6   2.*ALZ*AL1+PI2/6.-AL1**2)+(2./Z+2.+Z)*S2(Z))
         PGQ=-PGQ
         IF(IORD.EQ.2)THEN
            PGQ=PGQ+AL/8.*P2GQA(Z,NEFF)
            PGG2=PGG2+AL/8.*P2GGA(Z,NEFF)
            PGG3=PGG3+AL/8.*p2GGB(Z,NEFF)
         ENDIF
         FINTG=FINTG+.5*(1.-X)*T*AL*AL*WI(I)*(FSZ*PGQ+FGZ*PGG2+PGG3*(
     2   FGZ-FX))
      ENDIF
    1 CONTINUE
      RETURN
      END 

      SUBROUTINE FINTGS(X,FS,FG,FINTG)
      IMPLICIT REAL*8 (A-H,O-Z)
C  INTEGRATES SINGLET TERM - SEE NOTES
      DIMENSION FS(60),FG(60) 
c      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      COMMON/GRPTHY/ FLAVOR
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/ALQ2/T 
      COMMON/CONSTANTS/PI,PI2
      common/param/xc(100)
      xmu2=xc(1)**2*dexp(t)
      al=alpha_s(ILOOP,xmu2,xc(1),neff)/(2.*pi)
      flavor=neff
      FX=GETFV(X,FS)
      FINTG=FX*T*AL*(2.d0+8.d0/3.d0*DLOG(1.-X))
      CF=4.d0/3.d0
      IF(IORD.GT.0)THEN
         TR=FLAVOR/2.d0
         CA=3.d0
         Z3=1.202056903D0
         AL1=DLOG(1.d0-X)
         PFF1=-2.d0/9.d0*CF*TR*(3.d0/4.d0+PI2+10.d0*AL1)+CF*CF*(3.d0-
     2   4.d0*PI2+48.d0*Z3)/8.d0+CA*CF*(17.d0/3.d0+44.d0/9.d0*PI2
     3   -24.d0*Z3+(536.d0/9.d0-8.d0*PI2/3.)*AL1)/8.d0
         IF(IORD.EQ.2)THEN
            PFF1=PFF1+AL/8.d0*P2NSPC(X,NEFF)
         ENDIF
         FINTG=FINTG+T*AL*AL*FX*PFF1
      ENDIF
      DO 1 I=1,NTERMS
      Z=0.5d0*(1.d0-X)*XI(I)+0.5d0*(1.d0+X)
      FSZ=GETFV(X/Z,FS)
      FGZ=GETFV(X/Z,FG)
      FINTG=FINTG+.5d0*(1.d0-X)*T*AL*WI(I)*(CF*((1.d0+Z*Z)*FSZ
     2-2.*FX)/(1.d0-Z)
     2+FLAVOR*((1.d0-Z)**2+Z*Z)*FGZ)
      IF(IORD.GT.0)THEN
         AL1=DLOG(1.d0-Z)
         ALZ=DLOG(Z)
         PFF2=2.d0/9.d0*CF*TR*(20.d0/Z-19.d0+65.d0*Z-56.d0*Z*Z
     2   +6.d0*(2.d0+8.d0*Z+4.d0*
     2   Z*Z-1.d0/(1.d0-Z))*ALZ-9.d0*(1.d0+Z)*ALZ**2)
     3   +CF*CF*(-1.d0+Z+(2.d0-3.d0/
     3   (1.d0-Z))*ALZ-.5d0*(1.+Z)*ALZ**2
     4   -2.d0*(1.d0+Z*Z)/(1.-Z)*ALZ*AL1+
     4   2.d0*(1.d0+Z*Z)/(1.d0+Z)*S2(Z))+CA*CF/8.d0
     5   *(4.d0/9.d0*(17.d0-151.d0*Z)
     5   +4.d0*(1.d0+Z*Z)/(1.d0-Z)*ALZ*(11.d0/3.d0+ALZ)
     6   +4.d0*(1.d0+Z)*PI2/3.d0-8.d0*
     6   (1.d0+Z*Z)/(1.d0+Z)*S2(Z))
c        1         2         3         4         5         6         7 *
         PFF3=(-20.d0/9.d0*CF*TR+(67.d0/9.d0-PI2/3.d0)*CA*CF)/(1.d0-Z)
         PQG=-CF*TR*(-14.d0+29.d0*Z-20.d0*Z*Z
     2   -(3.d0-4.d0*Z+8.d0*Z*Z)*ALZ-
     2   (1.d0-2.d0*Z+4.d0*Z*Z)*ALZ**2-8.d0*Z*(1.d0-Z)*AL1
     3   +2.d0*(1.d0-2.d0*Z+2.d0*Z*Z)*
     3   (2.d0*ALZ*AL1+PI2/3.d0-AL1**2))
     4   -CA*TR*(-2.d0/9.d0*(20.d0/Z-18.d0+225.d0*Z-
     4   218.d0*Z*Z)-2.d0/3.d0*(3.d0+24.d0*Z+44.d0*Z*Z)*ALZ
     5   +(3.d0+6.d0*Z+2.d0*Z*Z)*
     5   ALZ**2+8.d0*Z*(1.d0-Z)*AL1+(1.d0-2.d0*Z+2.d0*Z*Z)
     6   *(2.d0*AL1**2-PI2/3.d0)
     6   -2.d0*(1.d0+2.d0*Z+2.d0*Z*Z)*S2(Z))
         IF(IORD.EQ.2)THEN
            PQG=PQG+AL/8.d0*P2QGA(Z,NEFF)
            PFF2=PFF2+AL/8.d0*(P2PSA(Z,NEFF)+P2NSPA(Z,NEFF))
            PFF3=PFF3+AL/8.d0*P2NSB(Z,NEFF)
         ENDIF
         FINTG=FINTG+.5d0*(1.d0-X)*T*AL*AL*WI(I)*(FSZ*PFF2+PFF3*(FSZ-FX)
     2   +PQG*FGZ)
      ENDIF
    1 CONTINUE
      RETURN
      END 

      FUNCTION GETFV(X,FVL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FVL(60),xx(4),fx(4)
      COMMON/GRID/NX,XGRID(60)
      DO 1 I=1,NX 
      IF(X.LT.XGRID(I)) GO TO 2
    1 CONTINUE
c    2 I=I-1
    2 I=I-2
      IF(I.le.0) I=1
      if(i.gt.(nx-2))i=nx-2
      xx(1)=xgrid(i)
      xx(2)=xgrid(i+1)
      xx(3)=xgrid(i+2)
      fx(1)=fvl(i)
      fx(2)=fvl(i+1)
      fx(3)=fvl(i+2)
      if(i.eq.(nx-2))then
         xx(4)=1.
         fx(4)=0.
      else
         xx(4)=xgrid(i+3)
         fx(4)=fvl(i+3)
      endif
c      CALL POLINT(XX,FX,4,X,ANS,DY)
      CALL POLINT4(XX,FX,X,ANS)
      getfv=ans
      RETURN
      END 

      SUBROUTINE WATE32
      IMPLICIT REAL*8 (A-H,O-Z)
C  32 POINT GAUSSIAN QUADRATURE ROUTINE 
      DIMENSION X(16),W(16)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      NTERMS=32
      X(1)=0.048307665687738316235
      X(2)=0.144471961582796493485
      X(3)=0.239287362252137074545
      X(4)=0.331868602282127649780
      X(5)=0.421351276130635345364
      X(6)=0.506899908932229390024
      X(7)=0.587715757240762329041
      X(8)=0.663044266930215200975
      X(9)=0.732182118740289680387
      X(10)=0.794483795967942406963
      X(11)=0.849367613732569970134
      X(12)=0.896321155766052123965
      X(13)=0.934906075937739689171
      X(14)=0.964762255587506430774
      X(15)=0.985611511545268335400
      X(16)=0.997263861849481563545
      W(1)=0.096540088514727800567
      W(2)=0.095638720079274859419
      W(3)=0.093844399080804565639
      W(4)=0.091173878695763884713
      W(5)=0.087652093004403811143
      W(6)=0.083311924226946755222
      W(7)=0.078193895787070306472
      W(8)=0.072345794108848506225
      W(9)=0.065822222776361846838
      W(10)=0.058684093478535547145
      W(11)=0.050998059262376176196
      W(12)=0.042835898022226680657
      W(13)=0.034273862913021433103
      W(14)=0.025392065309262059456
      W(15)=0.016274394730905670605
      W(16)=0.007018610009470096600
      DO 1 I=1,16
      XI(I)=-X(17-I)
      WI(I)=W(17-I) 
      XI(I+16)=X(I) 
      WI(I+16)=W(I) 
    1 CONTINUE
      DO 2 I=1,32
    2 XX(I)=0.5*(XI(I)+1.)
      XX(33)=1.0
      RETURN
      END 

      SUBROUTINE WATE16
      IMPLICIT REAL*8 (A-H,O-Z)
C  16 POINT GAUSSIAN QUADRATURE ROUTINE 
      DIMENSION X(8),W(8)
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
      NTERMS=16
      X(1)=0.095012509837637440185D0
      X(2)=0.281603550779258913230D0
      X(3)=0.458016777657227386342D0
      X(4)=0.617876244402643748447D0
      X(5)=0.755404408355003033895D0
      X(6)=0.865631202387831743880D0
      X(7)=0.944575023073232576078D0
      X(8)=0.989400934991649932596D0
      W(1)=0.189450610455068496285D0
      W(2)=0.182603415044923588867D0
      W(3)=0.169156519395002538189D0
      W(4)=0.149595988816576732081D0
      W(5)=0.124628971255533872052D0
      W(6)=0.095158511682492784810D0
      W(7)=0.062253523938647892863D0
      W(8)=0.027152459411754094852D0
      DO 1 I=1,8
      XI(I)=-X(9-I)
      WI(I)=W(9-I) 
      XI(I+8)=X(I) 
      WI(I+8)=W(I) 
    1 CONTINUE
      DO 2 I=1,16
    2 XX(I)=0.5*(XI(I)+1.)
      XX(17)=1.0
      RETURN
      END 

      SUBROUTINE WATE8
      IMPLICIT REAL*8 (A-H,O-Z)
C  8 POINT GAUSSIAN QUADRATURE ROUTINE
      COMMON/GAUSS8/XI(8),WI(8),NTERMS,XX(9)
      NTERMS=8
      XI(4)=-0.183434642495650D0
      XI(3)=-0.525532409916329D0
      XI(2)=-0.796666477413627D0
      XI(1)=-0.960289856497536D0
      WI(4)=0.362683783378362D0
      WI(3)=0.313706645877887D0
      WI(2)=0.222381034453374D0
      WI(1)=0.101228536290376D0
      DO 1 I=5,8
      XI(I)=-XI(9-I)
    1 WI(I)=WI(9-I) 
      DO 2 I=1,8
    2 XX(I)=0.5*(XI(I)+1.D0)
      XX(9)=1.
      RETURN
      END 


      FUNCTION XMNT32(F,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C  CALULATES NTH MOMENT OF F USING 32 POINT GAUSSIAN
C  WARNING!! THIS IS NOT ENTIRELY ADEQUATE FOR N=2
C  AT HIGH Q**2 WHERE A LOW-X SPIKE DEVELOPS. A VARIABLE
C  TRANSFORMATION IS NEEDED IF ONE WANTS TO TEST THE
C  PROGRAM AT LOW X AND HIGH Q**2.
      DIMENSION F(1)
      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)
      XMNT32=0.0D0
      M=N-2
C
C  MODIFIED FOR N=1 TO OBTAIN GREATER ACCURACY
C
      IF(N.GT.1)THEN
      DO 1 I=1,NTERMS
      X=XX(I)
    1 XMNT32=XMNT32+0.5*X**M*F(I)*WI(I)
      ELSE
      DO 2 I=1,NTERMS
      X=XX(I)**2
      XMNT32=XMNT32+0.5*XX(I)**M*F(I)*WI(I)*2.
    2 CONTINUE
      ENDIF
      RETURN
      END 


      SUBROUTINE ADIMEN
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ASFMOM/XMSTRS(25),XMSTRG(25),GPSI(25),GPSIA(25),GAPSI(25),
     $ GAA(25),GPLUS(25),GMINUS(25),ALPHAN(25),BETAN(25),EPSILN(25)
      COMMON/GRPTHY/ FLAVOR
 2000 FORMAT(10X,'CONSTANTS OF ASYMPTOTIC FREEDOM FOR',F3.0,' FLAVORS',
     $/10X,'NOTATION FOLLOWS A. J. BURAS, NUCL. PHYS. B125,125(1977)')
C    $/10X,'NOTATION FOLLOWS GROSS AND WILCZEK')
C    $/10X,'NOTATION FOLLOWS FLORATOS, ROSS AND SACHRAJDA') 
      WRITE(6,2001) 
 2001 FORMAT(/3X,'N',5X,'G(PSI,PSI)',2X,'G(PSI,A)',4X,'G(A,PSI)',4X,
     $'G(A,A)',6X,'G(+)',8X,'G(-)',8X,'ALPHA',7X,'BETA',8X,'EPSILON'/)
      B=(33.-2.*FLAVOR)/3.
      S=0.0
      DO 1 N=2,25
      S=S+1./N
      GPSI(N)=8./(6.*B)*(1.-2./N/(N+1.)+4.*S)
      IF(IPOLZN) 60,60,70
   60 CONTINUE
      GPSIA(N)=FLAVOR/(2.*B)*(8./(N+2.)+16./N/(N+1.)/(N+2.))
C     GPSIA(N)=-0.5*GPSIA(N)
      GAPSI(N)=8./(6.*B)*(1./(N+1.)+2./N/(N-1.))
C     GAPSI(N)=-2.*GAPSI(N)
      GAA(N)=3./B*(1./3.-4./N/(N-1.)-4./(N+1.)/(N+2.)+4.*S)+2.*FLAVOR/3.
     $/B
      GO TO 80
   70 CONTINUE
      GPSIA(N)=FLAVOR/(2.*B)* 8.*(N-1.)/N/(N+1.)
      GAPSI(N)=8./(6.*B)*(N+2.)/N/(N+1.)
      GAA(N)=3./B*(1./3.-8./N/(N+1.)+4.*S)+2.*FLAVOR/3./B
   80 CONTINUE
C        1         2         3         4         5         6         7 *
      GPLUS(N)=0.5*(GPSI(N)+GAA(N)+DSQRT((GPSI(N)-GAA(N))**2+4.*GPSIA(N)
     $ *GAPSI(N)))
      GMINUS(N)=0.5*(GPSI(N)+GAA(N)-DSQRT((GPSI(N)-GAA(N))**2+
     $4.*GPSIA(N)*GAPSI(N)))
      ALPHAN(N)=GPSIA(N)*GAPSI(N)/(GPSIA(N)*GAPSI(N)+(GMINUS(N)-GPSI(N))
     $**2)
      BETAN(N)=0.5*ALPHAN(N)*(GPSI(N)-GMINUS(N))/GAPSI(N)
      EPSILN(N)=ALPHAN(N)*(1.-ALPHAN(N))/BETAN(N) 
      WRITE(6,1000) N,GPSI(N),GPSIA(N),GAPSI(N),GAA(N),GPLUS(N),
     $GMINUS(N),ALPHAN(N),BETAN(N),EPSILN(N)
 1000 FORMAT(1X,I3,3X,9(F10.7,2X))
    1 CONTINUE
      RETURN
      END 


      SUBROUTINE WATE96
  !*******************************************************************
  !*****              *****
  !***** THE X(I) AND W(I) ARE THE DIRECT OUTPUT FROM A PROGRAM  *****
  !***** USING NAG ROUTINE D01BCF TO CALCULATE THE        *****
  !***** GAUSS-LEGENDRE WEIGHTS FOR 96 POINT INTEGRATION.        *****
  !***** THEY AGREE TO TYPICALLY 14 DECIMAL PLACES WITH THE      *****
  !***** TABLE IN ABRAMOWITZ & STEGUN, PAGE 919.         *****
  !*****              *****
  !***** ---->   PETER HARRIMAN, APRIL 3RD 1990.         *****
  !*****              *****
  !*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(48),W(48)
      COMMON/GAUS96/XI(96),WI(96),nterms,XX(97)
      NTERMS=96
      X( 1)=   0.01627674484960183561
      X( 2)=   0.04881298513604856015
      X( 3)=   0.08129749546442434360
      X( 4)=   0.11369585011066471632
      X( 5)=   0.14597371465489567682
      X( 6)=   0.17809688236761733390
      X( 7)=   0.21003131046056591064
      X( 8)=   0.24174315616383866556
      X( 9)=   0.27319881259104774468
      X(10)=   0.30436494435449495954
      X(11)=   0.33520852289262397655
      X(12)=   0.36569686147231213885
      X(13)=   0.39579764982890709712
      X(14)=   0.42547898840729897474
      X(15)=   0.45470942216774136446
      X(16)=   0.48345797392059470382
      X(17)=   0.51169417715466604391
      X(18)=   0.53938810832435567233
      X(19)=   0.56651041856139533470
      X(20)=   0.59303236477757022282
      X(21)=   0.61892584012546672523
      X(22)=   0.64416340378496526886
      X(23)=   0.66871831004391424358
      X(24)=   0.69256453664216964528
      X(25)=   0.71567681234896561582
      X(26)=   0.73803064374439816819
      X(27)=   0.75960234117664555964
      X(28)=   0.78036904386743123629
      X(29)=   0.80030874413913884180
      X(30)=   0.81940031073792957139
      X(31)=   0.83762351122818502758
      X(32)=   0.85495903343459936363
      X(33)=   0.87138850590929436968
      X(34)=   0.88689451740241818933
      X(35)=   0.90146063531585023110
      X(36)=   0.91507142312089592706
      X(37)=   0.92771245672230655266
      X(38)=   0.93937033975275308073
      X(39)=   0.95003271778443564022
      X(40)=   0.95968829144874048809
      X(41)=   0.96832682846326217918
      X(42)=   0.97593917458513455843
      X(43)=   0.98251726356301274934
      X(44)=   0.98805412632962202890
      X(45)=   0.99254390032376081654
      X(46)=   0.99598184298720747465
      X(47)=   0.99836437586317963722
      X(48)=   0.99968950388322870559
      W( 1)=   0.03255061449236316962
      W( 2)=   0.03251611871386883307
      W( 3)=   0.03244716371406427668
      W( 4)=   0.03234382256857594104
      W( 5)=   0.03220620479403026124
      W( 6)=   0.03203445623199267876
      W( 7)=   0.03182875889441101874
      W( 8)=   0.03158933077072719007
      W( 9)=   0.03131642559686137819
      W(10)=   0.03101033258631386231
      W(11)=   0.03067137612366917839
      W(12)=   0.03029991542082762553
      W(13)=   0.02989634413632842385
      W(14)=   0.02946108995816795100
      W(15)=   0.02899461415055528410
      W(16)=   0.02849741106508543861
      W(17)=   0.02797000761684838950
      W(18)=   0.02741296272602931385
      W(19)=   0.02682686672559184485
      W(20)=   0.02621234073567250055
      W(21)=   0.02557003600534944960
      W(22)=   0.02490063322248370695
      W(23)=   0.02420484179236479915
      W(24)=   0.02348339908592633665
      W(25)=   0.02273706965832950717
      W(26)=   0.02196664443874448477
      W(27)=   0.02117293989219144572
      W(28)=   0.02035679715433347898
      W(29)=   0.01951908114014518992
      W(30)=   0.01866067962741165898
      W(31)=   0.01778250231604547316
      W(32)=   0.01688547986424539715
      W(33)=   0.01597056290256253144
      W(34)=   0.01503872102699521608
      W(35)=   0.01409094177231515264
      W(36)=   0.01312822956696188190
      W(37)=   0.01215160467108866759
      W(38)=   0.01116210209983888144
      W(39)=   0.01016077053500880978
      W(40)=   0.00914867123078384552
      W(41)=   0.00812687692569928101
      W(42)=   0.00709647079115442616
      W(43)=   0.00605854550423662775
      W(44)=   0.00501420274292825661
      W(45)=   0.00396455433844564804
      W(46)=   0.00291073181793626202
      W(47)=   0.00185396078894924657
      W(48)=   0.00079679206555731759
      DO I=1,48
         XI(I)=-X(49-I)
         WI(I)=W(49-I)
         XI(I+48)=X(I)
         WI(I+48)=W(I)
      END DO
      DO I=1,96
         XX(I)=0.5*(XI(I)+1.)
      END DO
      XX(97)=1.0
      EXPON=1.0
      DO I=1,96
         YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
         WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
         XI(I)=YI
         XX(I)=0.5*(1.+YI)
      END DO
      RETURN
      END SUBROUTINE WATE96


      SUBROUTINE WATE64
  !*******************************************************************
  !***** THE X(I) AND W(I) ARE THE DIRECT OUTPUT FROM A PROGRAM  *****
  !***** USING NAG ROUTINE D01BCF TO CALCULATE THE               *****
  !***** GAUSS-LEGENDRE WEIGHTS FOR 64 POINT INTEGRATION.        *****
  !***** https://keisan.casio.com/exec/system/1329114617         *****
  !***** Shujie Li, 052020                                       *****
  !*******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(32),W(32)
      COMMON/GAUS64/XI(64),WI(64),nterms,XX(65)
      NTERMS=64
      X(1)  =  0.0243502927
      X(2)  =  0.0729931218
      X(3)  =  0.1214628193
      X(4)  =  0.1696444204
      X(5)  =  0.2174236437
      X(6)  =  0.2646871622
      X(7)  =  0.311322872
      X(8)  =  0.3572201583
      X(9)  =  0.402270158
      X(10) =  0.4463660173
      X(11) =  0.4894031457
      X(12) =  0.531279464
      X(13) =  0.5718956462
      X(14) =  0.6111553552
      X(15) =  0.6489654713
      X(16) =  0.6852363131
      X(17) =  0.7198818502
      X(18) =  0.7528199073
      X(19) =  0.7839723589
      X(20) =  0.8132653151
      X(21) =  0.8406292963
      X(22) =  0.8659993982
      X(23) =  0.889315446
      X(24) =  0.9105221371
      X(25) =  0.9295691721
      X(26) =  0.9464113749
      X(27) =  0.9610087997
      X(28) =  0.9733268278
      X(29) =  0.9833362539
      X(30) =  0.9910133715
      X(31) =  0.9963401168
      X(32) =  0.9993050417
      W(1)  =  0.048690957
      W(2)  =  0.0485754674
      W(3)  =  0.0483447622
      W(4)  =  0.0479993886
      W(5)  =  0.0475401657
      W(6)  =  0.0469681828
      W(7)  =  0.0462847966
      W(8)  =  0.0454916279
      W(9)  =  0.0445905582
      W(10) =  0.0435837245
      W(11) =  0.0424735151
      W(12) =  0.0412625632
      W(13) =  0.0399537411
      W(14) =  0.0385501532
      W(15) =  0.0370551285
      W(16) =  0.0354722133
      W(17) =  0.0338051618
      W(18) =  0.0320579284
      W(19) =  0.0302346571
      W(20) =  0.0283396726
      W(21) =  0.0263774697
      W(22) =  0.0243527026
      W(23) =  0.0222701738
      W(24) =  0.0201348232
      W(25) =  0.0179517158
      W(26) =  0.0157260305
      W(27) =  0.0134630479
      W(28) =  0.0111681395
      W(29) =  0.0088467598
      W(30) =  0.006504458
      W(31) =  0.0041470333
      W(32) =  0.0017832807

      DO I=1,32
         XI(I)=-X(33-I)
         WI(I)=W(33-I)
         XI(I+32)=X(I)
         WI(I+32)=W(I)
      END DO
      DO I=1,64
         XX(I)=0.5*(XI(I)+1.)
      END DO
      XX(65)=1.0
C       EXPON=1.0
C       DO I=1,64
C          YI=2.*(0.5*(1.+XI(I)))**EXPON-1.
C          WI(I)=WI(I)/(1.+XI(I))*(1.+YI)*EXPON
C          XI(I)=YI
C          XX(I)=0.5*(1.+YI)
C       END DO
      RETURN
      END SUBROUTINE WATE64



C     =======================================
      SUBROUTINE WATE8L
      IMPLICIT REAL*8 (A-H,O-Z)
C  8 POINT LAGUERRE INTGRATOR (FOR EXPONENTIALS OR OTHER
C  STEEPLY FALLING FUNCTIONS
      COMMON/LAGUER/XL(8),WL(8),NTERML
      XL(1)=0.1702796
      XL(2)=0.9037018
      XL(3)=2.2510866
      XL(4)=4.2667002
      XL(5)=7.0459054
      XL(6)=10.758516
      XL(7)=15.740679
      XL(8)=22.863132
      WL(1)=0.4377234
      WL(2)=1.0338693
      WL(3)=1.6697098
      WL(4)=2.3769247
      WL(5)=3.2085409
      WL(6)=4.2685755
      WL(7)=5.8180834
      WL(8)=8.9062262
      NTERML=8
      RETURN
      END 
