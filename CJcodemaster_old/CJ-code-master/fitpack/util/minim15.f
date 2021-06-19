      SUBROUTINE MINIM(NX,NNPTS,DATA,IDRV,nofit) 
      implicit real*8(a-h,o-z)
      CHARACTER*10 IPNAME(100)
      character*4 fnot,froz,blnk,tol,err
c      REAL*8 MODALF 
c      DIMENSION ALPHA(50,50),MODALF(50,50),BETA(50),COV(50,50) 
      REAL*8 MODALF(333,333)  ! dimension = 333 to use with 'minv' subroutine 
      DIMENSION ALPHA(50,50),BETA(50),COV(50,50) 
      DIMENSION OLDPAR(100),DATA(1,1)
      COMMON/ERRORS/COV
      COMMON/AMAT/ALPHA,BETA
      COMMON/MINIMC/IPNAME,PLB(100),PUB(100),TOLRNC,IPAR(100),NSTEPS,
     1IPRNT,IERR,IBND,CHISQR 
      COMMON/MINIMB/UNCRT(100),ERMIN(50),PWATE(100),
     1     IFREP(50),NPARAM,NVAR,IREJ
      COMMON/CURPAR/CURPAR(100)
      DATA FROZ/'FIX '/,BLNK/'    '/,FNOT/' NOT'/ 
      dimension bevalpha(50,50) ! Bevington's alpha matrix (Hessian)
      common/amatrix/bevalpha          !'ALPHA' in 'amat' is the normalized
                                       ! error matrix
      common /Step/ NIStep
C
C  MINIMIZATION ROUTINE USING THE MARQUARDT ALGORITHM
C  FROM BEVINGTON, PAGE 235. SEE PROGRAM CURFIT.
C  VERSION CREATED BY D.W. DUKE CIRCA 1978.
C  COMMENTS AND SOME MODS ADDED BY J.F. OWENS, A.Accardi
C

      ! (aa) debugging messages to screen?
      !      0 = no 
      !      1 = prints only chi^2 and params at end of each minimization step
      !      2 = also breaks down the Marquardt algorithms into components
      !      3 = also prints beta and alpha to screen 
      idebug = 0  

C INPUT AND OUTPUT UNIT NUMBERS
      MI=5
      MO=4
      FLAMDA=0.001 

C COUNT NUMBER OF PARAMETERS VARIED
      NVAR=0 
      DO 340 I=1,NPARAM 
         IF(PWATE(I).EQ.0.0) GO TO 340 
         NVAR=NVAR+1 
         IFREP(NVAR)=I 
C IFREP IS THE SEQUENCE OF PARAMETER NUMBERS FOR THOSE BEING VARIED
 340  CONTINUE



      
***** 1st step: *****
***** calculates initial chi^2
      if(idebug.ge.2) then
         write(6,'(A,I3)') '-----------------------------'
         write(6,'(A,I3)') '  Initial chi^2 calculation'
         write(6,'(A,I3)') '-----------------------------'
      end if
       
C     CALCULATE CHI SQUARE
      if (nofit.eq.-1) then
         ! if calculating with errorbands (nofit=-1) then:
         ! (1) calculates the derivatives of the observables
         CALL DATSCN(idrv,CHISQR,NX,NNPTS,DATA)
         ! (2) calculates also the correlated shifts
         CALL DATSCN(1,CHISQR,NX,NNPTS,DATA)
      else
         if(nofit.eq.0) nvar=0  ! If doing calculations w/o uncertainties
         CALL DATSCN(1,CHISQR,NX,NNPTS,DATA)
      end if
      if(idebug.ge.2) write(6,*) '----- initial chi^2 =', chisqr       
c      write(6,*), ierr,iprnt
      IF(IERR.EQ.1) GO TO 1
      IF(IPRNT.EQ.0) GO TO 1
      ! outputs chi^2 and params to file
      WRITE(MO,1111)(IP,IPNAME(IP),CURPAR(IP),PWATE(IP),IP=1,NPARAM)
 1111 FORMAT(5X,I5,2X,A10,'  VALUE=',G14.5,'  WEIGHT=',G14.5)
      WRITE(MO,1305) TOLRNC,NSTEPS 
      WRITE(MO,1300) CHISQR 
      if (nofit.eq.1) WRITE(MO,1310)
      if (nofit.eq.0.or.nofit.eq.-1) return ! exits when only calculating


***** 2nd step: *****
***** applies Marquardt's algorithm (Bevington, 3rd ed., page 162)

 1    DO 100 NOSTEP=1,NSTEPS    ! max of nstep times
         if(idebug.ge.2) then
            write(6,'(A,I3)') '--------------------------'
            write(6,'(A,I3)') '  Minimization step #',nostep
            write(6,'(A,I3)') '--------------------------'
         end if
         IF(FLAMDA.GT.0.001) FLAMDA=0.001   ! restarst from lamda=0.001
         IF(FLAMDA.LT.1.E-10) FLAMDA=1.E-10 ! unless previous step was
                                            ! successful with that or less
C SAVE OLD PARAMETER VALUES
         OLDCHI=CHISQR 
         DO 5 I=1,NPARAM 
 5          OLDPAR(I)=CURPAR(I) 
C EVALUATE ALPHA AND BETA
         if (idebug.ge.2) write(6,*)' --- evaluating alpha, beta'
         CALL DATSCN(IDRV,DUM,NX,NNPTS,DATA)
         NTRIES=0
         
         if (idebug.ge.3) then
            write(6,*) ' --- beta ='
               write(6,'(6X,35ES12.4)') (beta(i),i=1,nvar)
            write(6,*) ' --- alpha ='
            do i = 1, nvar
               write(6,'(6X,35ES12.4)') (alpha(i,j),j=1,i)
            end do
         end if
         
         if (idebug.ge.2) write(6,*)' --- lambda search begins'
**** lambda search iterations
 10      continue
         DO 20 I=1,NVAR       ! calc normalized alfa' from Bevington (8.39) 
            DO 15 J=1,NVAR    ! (normalized for easier inversion)
 15            MODALF(I,J)=ALPHA(I,J)/ SQRT(ALPHA(I,I)*ALPHA(J,J)) 
 20         MODALF(I,I)=1.+FLAMDA
            
C INVERT MODALF 
         CALL MINV(MODALF,NVAR,DET)  ! modalf now = modalf^-1
         NTRIES=NTRIES+1 
         NISTEP=NOSTEP
C UPDATE PARAMETERS [Bevington's (8.41)]
         DO 25 I=1,NVAR
            K=IFREP(I) 
            DO 26 J=1,NVAR    ! NOTE: undoes the normalization in modalfa  
               CURPAR(K)=CURPAR(K)+MODALF(I,J)*BETA(J)
     1              / SQRT(ALPHA(I,I)*ALPHA(J,J))
 26         continue
 25      continue
* calculates new chi^2
 32      CALL DATSCN(1,CHISQR,NX,NNPTS,DATA)
         if(idebug.ge.2) write(6,'(A,I2,A,E10.3,A,E12.5)')
     &        ' -- step #',ntries, ' lamda=',flamda,' chi^2=',chisqr 
* if less than previous chi^2 go to covergence check...
         IF(CHISQR.LT.OLDCHI) GO TO 40 
* ...otherwise, tries a bigger lambda (max of 7 times)
         ! --- longer trial for testing --- IF(NTRIES.GT.10) GO TO 110 
         IF(NTRIES.GT.6) GO TO 110 
         FLAMDA=10.*FLAMDA 
         DO 35 I=1,NPARAM 
 35         CURPAR(I)=OLDPAR(I) 
         GO TO 10 
C CHECK FOR CONVERGENCE
 40      IF((OLDCHI-CHISQR).LT.TOLRNC) GO TO 110 ! if tol reached, exits
C CHECK NUMBER OF STEPS 
         IF(NOSTEP.EQ.NSTEPS) GO TO 110 ! if no tol, but too many steps, exits
C prints steps results to file (if iprnt=/0)
         IF(IPRNT.EQ.0) GO TO 100
         ERR=BLNK 
         IF(PWATE(1).EQ.0.0) ERR=FROZ
         IF(IPRNT.NE.0) WRITE(MO,1320)NISTEP,CHISQR
         DO 60 I=1,NPARAM 
            ERR=BLNK 
            IF(PWATE(I).EQ.0.0) ERR=FROZ 
            IF(IPRNT.NE.0) WRITE(MO,1330)IPAR(I),IPNAME(I),CURPAR(I),ERR 
 60      CONTINUE
         ! ... writes also to screen for debugging purposes
         if (idebug.ge.1) then
            IF(IPRNT.NE.0) WRITE(6,1320)NISTEP,CHISQR
            DO I=1,NPARAM 
               ERR=BLNK 
               IF(PWATE(I).EQ.0.0) ERR=FROZ 
               IF(IPRNT.NE.0)
     &              WRITE(6,1330) IPAR(I),IPNAME(I),CURPAR(I),ERR 
            END DO
            write(6,*)
         end if
         
***** end of step(2) *****

C At the succesful end of a good lambda search,
C steps back the last increment of alfa
 100     FLAMDA=FLAMDA/10.0

       
***** step 3 *****         
***** Minimumfound: Evaluates errors and prints final results to file
 110  continue
      if (idebug.ge.1) then
         write(6,*) '------------------------------'
         write(6,*) '       MINIMUM FOUND!         '
         write(6,*) '------------------------------'
      end if         
      IF(IERR.EQ.1) RETURN  ! ...unless ierr=1
C EVALUATE ERRORS AND PRINT OUT RESULTS OF FIT 
      if (idebug.ge.1) then
         write(6,*) '---   evaluating errors   ----'
      end if         
      CALL DATSCN(IDRV,DUM,NX,NNPTS,DATA) 
      DO 120 I=1,NVAR 
      DO 120 J=1,I 
      bevalpha(i,j) = alpha(i,j) ! stores the Hessian in 'bevalpha'
      bevalpha(j,i) = bevalpha(i,j) 
      MODALF(I,J)=ALPHA(I,J)/ SQRT(ALPHA(I,I)*ALPHA(J,J)) 
  120 MODALF(J,I)=MODALF(I,J) 
      CALL MINV(MODALF,NVAR,DET) 
      DO 130 I=1,NVAR 
      AMOD=MODALF(I,I)/ALPHA(I,I) 
      IMOD=1 
      IF(AMOD.LT.0.0) IMOD=-1 
      UNCRT(IFREP(I))=IMOD*SQRT(IMOD*AMOD)
 130  continue
      CHIPNT=CHISQR/(NNPTS-IREJ-NVAR) 
      IF(IPRNT.EQ.0) GO TO 170
      WRITE(MO,1320)NISTEP,CHISQR 
      WRITE(MO,1335)(IPAR(I),IPNAME(I),CURPAR(I),UNCRT(I),I=1,NPARAM) 
      WRITE(MO,1340) NISTEP
      ! ... writes also to screen for debugging purposes
      if (idebug.ge.1) then
         WRITE(6,1320)NISTEP,CHISQR 
         WRITE(6,1335)(IPAR(I),IPNAME(I),CURPAR(I),UNCRT(I),I=1,NPARAM) 
         WRITE(6,1340) NISTEP
      end if
      TOL=BLNK 
      IF(dabs(OLDCHI-CHISQR).GT.TOLRNC) TOL=FNOT 
      WRITE(MO,1350) TOL 
      IF(IREJ.NE.0) WRITE(MO,1360) IREJ 
  160 WRITE(MO,1345)CHIPNT 
      DO 180 I=1,NVAR 
      DO 180 J=1,NVAR 
  180 MODALF(I,J)=MODALF(I,J)/SQRT(ALPHA(I,I)*ALPHA(J,J)) 
      DO 182 I=1,NVAR 
      DO 181 J=1,NVAR 
      COV(I,J)=MODALF(I,J)
  181 ALPHA(I,J)=MODALF(I,J)/SQRT(MODALF(I,I)*MODALF(J,J))
  182 WRITE(MO,1380)(ALPHA(I,J),J=1,NVAR) 
      WRITE(MO,1370)
  170 CONTINUE
      RETURN
 1300 FORMAT(28H0THE INITIAL CHI SQUARED IS ,1PE12.5) 
 1305 FORMAT(1H0,17HTHE TOLERANCE IS ,F9.4/32H THE MAXIMUM NUMBER OF STEPS IS ,I3) 
 1310 FORMAT(33H0  STEP  CHI SQUARED    PARAMETER,6X,5HVALUE,9X,5HERROR) 
 1320 FORMAT(1H ,60(1H-)/3H   ,I3,3X,1PE12.5) 
 1330 FORMAT(22X,I2,2X,A8,2X,1PE11.4,7X,A4) 
 1335 FORMAT(22X,I2,2X,A8,2X,1PE11.4,3X,1PE11.4) 
 1340 FORMAT(35H0 MINIMIZATION TERMINATED WITH STEP,I3) 
 1345 FORMAT(39H0 CHI SQUARED PER DEGREE OF FREEDOM IS ,1PE10.3) 
 1350 FORMAT(22H  TOLERANCE CRITERION ,A4,10H SATISFIED)
 1360 FORMAT(12H0 THERE WERE,I3,40H DATA POINTS WITH ZERO OR NEGATIVE ERROR)
 1370 FORMAT(/1H0,16X,32H***   ***   ***  ***   ***   ***) 
 1380 FORMAT(1X,50F7.3) 
      END 


      SUBROUTINE DATSCN(MODE,CHISQR,NX,NNPTS,DATA)
      implicit real*8 (a-h,o-z)
      character*10 ititle(49)
      dimension errn(49),iserr(49),itype(49),nflag(49),inorm(49),
     2itgt(49),icorr(49)
      DIMENSION ALPHA(50,50),BETA(50),C(49,100,333) ! is C unused ???
      COMMON/MINIMB/UNCRT(100),ERMIN(50),PWATE(100),
     1     IFREP(50),NPARAM,NVAR,IREJ
      COMMON/CURPAR/CURPAR(100)
      COMMON/AMAT/ALPHA,BETA
      common/expts/ititle,errn,iserr
      common/flags/itype,nflag,inorm,itgt,icorr
      common/normalization/nfit 
      common/cdf/covinv(33,33),cdfchisqr,ncdf,ncdfl,ncov
      common/minimd/nset(49),ncor(49),correr(49,400,333),
     2ainv(49,333,333),bvec(49,333),bab(49),am(333,333)
      DIMENSION DATA(NX,NNPTS),V(10),DRV(100),cdfxsec(100),
     2cdfdrv(100,100)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
      common/obs/oprime(10001,50) ! 1st derivative of observable wrt params
      save cdfxsec
      NN=NX-2 
C CALCULATE CHI SQUARE
      IF(MODE.NE.1) GO TO 100 
      CHISQR=0.0 
      do j=1,35
      do l=1,333
         bvec(j,l)=0.
      enddo
      enddo
c      DO 10 NPT=1,NNPTS 
      npt=0
      do 10 j=1,nfit
      do 10 l=1,nset(j)
         npt=npt+1
C FILL ARRAY V FOR DATA POINT LABELLED BY NPT
         DO 5 I=1,NN 
            V(I)=DATA(I+2,NPT)
 5       continue
         IF(DATA(2,NPT).LE.0.d0) GO TO 10 ! skips data with zero/neg error
         index=v(4)      
         if(itype(index).eq.201.and.icorr(index).ne.0)goto 10 
         fval=theory(1,npt,v,curpar)
         chisqr = chisqr 
     &        + chi2(fval,DATA(1,NPT),DATA(2,NPT))
c
c Calculate B array for use in chi square offset when
c correlated errors are used
c see Eq.(9) of CTEQ6.1 paper (JHEP 0207 (2002) 012)
c
         if(icorr(index).eq.2)then
            do kl=1,ncor(index)
               bvec(index,kl)=bvec(index,kl)+(data(1,npt)-fval)
     2         *correr(index,l,kl)/data(2,npt)**2
            enddo
         endif
   10 CONTINUE 
c
c  Add CDF chi square if correlations are used
c
      chisqr0 =chisqr
      if(ncov.ne.0)then
         chisqr=chisqr0+cdfchi2(data,NX,NNPTS,curpar,cdfxsec)
      endif
c
c  add normalization chi square penalty
c
      do 11 i=1,nfit
      if(errn(i).eq.0.d0.or.inorm(i).eq.100)goto 11
      if(iord.eq.0)goto 11
      chisqr=chisqr+chi2(curpar(inorm(i)),1d0,errn(i))
 11   continue
c
c Calculate chi square offset and add to chi square if 
c correlated errors are used
c see Eq.(10) of CTEQ6.1 paper (JHEP 0207 (2002) 012)
c     
      do j=1,nfit
         bab(j)=0.
         if(icorr(j).eq.2.and.itype(j).ne.201)then
            do l1=1,ncor(j)
               do l2=1,ncor(j)
                  bab(j)=bab(j)+bvec(j,l2)*ainv(j,l2,l1)*bvec(j,l1)
               enddo
            enddo
         endif
         chisqr=chisqr-bab(j)
      enddo   
      
      if (l.eq.1) then
            write(*,*) fval
      endif
      RETURN ! returns if mode=1 (only calculates the chi^2)

 100  continue  ! if mode=/=1 calulates derivative of chi^2 and alpha matrix

      IREJ=0                    
C CALCULATE ALPHA AND BETA
      do i=1,35
      if(icorr(i).eq.2)then
         do j=1,nvar
         do k=1,ncor(i)
            c(i,j,k)=0.
         enddo
         enddo
      endif
      enddo
      ! initializes alpha, beta and drv
      DO 20 I=1,NVAR
         BETA(I)=0.0 
         DRV(I)=0.0 
         DO 20 J=1,NVAR
   20    ALPHA(I,J)=0.0
c      DO 50 NPT=1,NNPTS 
c
c  Redo indexing by summing over experiments
c
      npt=0
      do 50 j=1,nfit
      do 50 l=1,nset(j)
      npt=npt+1
C IGNORE DATA POINTS WITH NEGATIVE OR ZERO ERRORS
      IF(DATA(2,NPT).LE.0.0) GO TO 45 
      DO 41 II=1,NN 
   41 V(II)=DATA(II+2,NPT) 
      FVAL=THEORY(1,NPT,V,CURPAR) 
C IF MODE EQUALS 3 USER SUPPLIES DERIVATIVES IN ARRAY DRV
C OTHERWISE, PROGRAM CALCULATES DERIVATIVES USING FINITE DIFFERENCE
      IF(MODE.EQ.3) GO TO 125 
      DO 30 I=1,NVAR
      K=IFREP(I)
      SV=CURPAR(K)
c
c  Current version has normalization parameters for k = 50-99. k=100 is 
c  reserved as the value of the normalization parameter for ratios, i.e., 
c  for the case where there is none.
c
c  For the normalization parameters, the derivative is easily calculated 
c  by hand, which will speed up the program.
c
      if(k.ge.50)then
         index=v(4)
         in=inorm(index)
         if(k.eq.in.and.k.ne.100)then
            drv(k)=-fval/curpar(k)
         else
            drv(k)=0.d0
         endif
      else
         CURPAR(K)=SV+PWATE(K)
C CALCULATE DERIVATIVE numerically
         DRV(K)=(THEORY(MODE,NPT,V,CURPAR)-FVAL)/PWATE(K)
      endif
c      if(itype(index).eq.201.or.itype(index).eq.202)then
c      endif
      index=v(4)
      if(itype(index).eq.201.and.icorr(index).ne.0)then
         nc=npt-ncdfl+1
         cdfdrv(k,nc)=drv(k)
      endif
      oprime(npt,i) = drv(k) ! saves the derivative of the observable
                             ! k=curpar indx; i=free par indx; npt=data pt. no.
   30 CURPAR(K)=SV ! restores original (central) value of curpar(k)
c
c 
  125 continue
c
c  choose method to calculate alpha and beta
c
      if(itype(index).eq.201.and.ncov.ne.0)then
         goto 50
      else
         DO 40 I=1,NVAR
         K=IFREP(I)
         BETA(I)=BETA(I)+(DATA(1,NPT)-FVAL)/DATA(2,NPT)*DRV(K)/
     2   DATA(2,NPT)
         if(icorr(j).eq.2)then
            do m=1,ncor(j)
               c(j,i,m)=c(j,i,m)+drv(k)*correr(j,l,m)/data(2,npt)**2
            enddo
         endif 
         DO 40 IJ=1,I 
         KL=IFREP(IJ) 
         ALPHA(I,IJ)=ALPHA(I,IJ)+DRV(K)*DRV(KL)/DATA(2,NPT)**2
 40      ALPHA(IJ,I)=ALPHA(I,IJ)
      endif
      GO TO 50 
   45 IREJ=IREJ+1 
   50 CONTINUE
c
c  add pieces to beta and alpha for the correlated errors
c
      do j=1,nfit
      if(icorr(j).eq.2)then
         do i=1,nvar
            do k1=1,ncor(j)
            do k2=1,ncor(j)
               beta(i)=beta(i)-c(j,i,k1)*ainv(j,k1,k2)*bvec(j,k2)
            enddo
            enddo
         do ij=1,nvar
            do k1=1,ncor(j)
            do k2=1,ncor(j)
               alpha(i,ij)=alpha(i,ij)-c(j,i,k1)*ainv(j,k1,k2)
     2                     *c(j,ij,k2)
            enddo
            enddo
         enddo
         enddo
      endif 
      enddo
c
c  fake call to 'theory' to renormalize the valence and gluon 
c  normalization parameters
c
      v(1)=data(3,1)
      v(2)=data(4,1)
      v(3)=data(5,1)
      v(4)=data(6,1)
      dum=theory(1,1,v,curpar) 

c
c  calculate alpha and beta for the cdf points if required
c
      if(ncov.ne.0)then
      do ik=1,nvar
         k=ifrep(ik)
         do i=1,ncdf
            do ip=1,ncdf
               beta(ik)=beta(ik)+(data(1,ncdfl+i-1)-cdfxsec(i))
     2         *covinv(i,ip)*cdfdrv(k,ip)
            if(ik.eq.9)then
            endif
            enddo
         enddo
         do il=1,ik
            l=ifrep(il)
            do i=1,ncdf
               do ip=1,ncdf
                  alpha(ik,il)=alpha(ik,il)+
     2            cdfdrv(k,i)*covinv(i,ip)*cdfdrv(l,ip)
               enddo
            enddo
         enddo
         alpha(il,ik)=alpha(ik,il)
      enddo
      endif
c
c  calculate alpha and beta for the normalization parameter terms
c
      do 12 i=1,nvar
      do 13 k=1,nfit
      if(errn(k).eq.0.d0.or.inorm(k).eq.100)goto 13
      if(ifrep(i).eq.inorm(k))then
         beta(i)=beta(i)+(1.d0-curpar(ifrep(i)))/errn(k)**2
         alpha(i,i)=alpha(i,i)+1.d0/errn(k)**2
      endif
 13   continue
 12   continue
      RETURN 
      END 



      SUBROUTINE MINV(ARRAY,NORDER,DET) 
      implicit real*8 (a-h,o-z)
      DIMENSION ARRAY(333,333),IK(333),JK(333)
C
C  MATRIX INVERSION ROUTINE FROM BEVINGTON, PAGE 302.
C
      DET=1. 
      DO 100 K=1,NORDER 
      AMAX=0.D0 
   21 DO 30 I=K,NORDER 
      DO 30 J=K,NORDER 
c      IF( ABS(AMAX).GT. ABS(ARRAY(I,J))) GO TO 30 
      IF(DABS(AMAX).GT.DABS(ARRAY(I,J))) GO TO 30 
   24 AMAX=ARRAY(I,J)
      IK(K)=I 
      JK(K)=J 
   30 CONTINUE 
      IF(AMAX.NE.0.D0) GO TO 41 
   32 DET=0.
      GO TO 140 
   41 I=IK(K) 
      IF(I-K) 21,51,43 
   43 DO 50 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=ARRAY(I,J) 
   50 ARRAY(I,J)=-SAVE 
   51 J=JK(K) 
      IF(J-K) 21,61,53 
   53 DO 60 I=1,NORDER
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=ARRAY(I,J) 
   60 ARRAY(I,J)=-SAVE 
   61 DO 70 I=1,NORDER 
      IF(I.EQ.K) GO TO 70 
   63 ARRAY(I,K)=-ARRAY(I,K)/AMAX 
   70 CONTINUE 
      DO 80 I=1,NORDER 
      DO 80 J=1,NORDER 
      IF(I.EQ.K.OR.J.EQ.K) GO TO 80 
   75 ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J) 
   80 CONTINUE 
      DO 90 J=1,NORDER 
      IF(J.EQ.K) GO TO 90 
   83 ARRAY(K,J)=ARRAY(K,J)/AMAX 
   90 CONTINUE 
      ARRAY(K,K)=1./AMAX 
  100 DET=DET*AMAX 
      DO 130 L=1,NORDER 
      K=NORDER-L+1
      J=IK(K) 
      IF(J.LE.K) GO TO 111 
  105 DO 110 I=1,NORDER 
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=-ARRAY(I,J) 
  110 ARRAY(I,J)=SAVE 
  111 I=JK(K) 
      IF(I.LE.K) GO TO 130 
  113 DO 120 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=-ARRAY(I,J) 
  120 ARRAY(I,J)=SAVE 
  130 CONTINUE 
  140 RETURN 
      END 
