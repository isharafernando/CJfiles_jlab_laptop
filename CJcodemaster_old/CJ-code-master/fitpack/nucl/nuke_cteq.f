c*************************************************************
c  This is a collection of nuclear corrections to Fe targets 
c
c  Routines:
c
c     * nuke_cteq: Kulagin-Petti correction  
c     * emcslac    parametrization of SLAC data 
c     * F2DTOA     some parametrization
c     * emcnmc     parametrization of emc/nmc data
c
c************************************************



c************************************************
c   Nuclear correction to structure functions
c   according to the calculation by
c 
c   S.A. Kulagin and R. Petti, 
c       Nucl.Phys.A 765 (2006) 126-187 [hep-ph/0412425]  
c 
c   Please refer to the above paper for results 
c   based upon the present code. 
c  
c   The routine interpolates tables for the ratio: 
c      R_sf = (Nuclear SF) / (nucleon SF)  
c   and therefore the result should be multiplied 
c   by the corresponding value of the nucleon structure 
c   function for isoscalar target, including target mass 
c   corrections. Note the isovector correction is already 
c   taken into account. The target mass corrections used 
c   for the present calculation are based on a 1/Q^2 expansion 
c   of the Georgi-Politzer formula, as described in the paper 
c   cited above (NPA 765, 126-187).  
c
c   xb   = Bjorken x
c   q2   = Q**2 in GeV**2
c
c   nsf  = structure function
c        = 1 2*x*F1
c        = 2 F2
c        = 3 x*F3
c            Note: for antineutrinos (kint<0) the sum of xF3(nu)+xF3(nubar)
c                  is provided instead to avoid numerical problems due to
c                  negative SF at small x values
c
c   ityp = nucleon type
c        = 0  actual average in nucleus
c
c   kint = interaction
c        = 2 neutrino (-2 antineutrino) CC
c
c   kord = order of QCD calculation
c        =  3 NNLO
c
c  ftyp  = target nucleus 
c        = 2 Iron (A=56, Z=26)
c    
c************************************************

      real*8 function nuke_cteq(xb,q2,nsf,ityp,kint1,kord,ftyp)
      implicit none 
*  
      integer nsf,ityp,kint,kord,kint0,kint1,ftyp,Nsys,ifirst
      integer nxb,nq,np,ntarg,i,n,m,kk,tt,nxbb,nqb,kx,kq
      integer len
      parameter(nxb=99,nq=20,np=81,ntarg=3)
      data ifirst /0/

      real syst,Nsigmas,rvec
      real*8 xb,q2,x,qsq,a,b,aa,f0,fp,fm
      real*8 qq(nq+1),xx(nxb)
      real*8 x1,delx,delx1,ss,xd,xlog1,dels 
* 
      real*8 xmin,xmax,qsqmin,qsqmax
      data xmin,xmax,qsqmin,qsqmax/1.d-4,1.0d0,5.00d0,1d3/

      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 f_nuke(nxb,nq+1,0:np,ntarg),b_nuke(nxb,nq+1,0:np,ntarg) 
      real*8 c_nuke(nxb,nq+1,0:np,ntarg),d_nuke(nxb,nq+1,0:np,ntarg)

      real*8 f_dnuke(nxb,nq+1,0:np,ntarg),b_dnuke(nxb,nq+1,0:np,ntarg)
      real*8 c_dnuke(nxb,nq+1,0:np,ntarg),d_dnuke(nxb,nq+1,0:np,ntarg)

      character locfile*120

      integer nport
c I/O channel to read the data
      data nport/23/

      integer ityps,kints,kords
      data ityps,kints,kords /-999,-999,-999/
      save ityps,kints,kords,xx,qq,f_nuke
      save b_nuke,c_nuke,d_nuke,b_dnuke,c_dnuke,d_dnuke
* 
      nuke_cteq=0.d0
      if ((nsf.lt.1).or.(nsf.gt.3)) return 
      if (q2.gt.2.d8) return 
* 
*...Check input flags 
* 
      if (ifirst.eq.0) then 
        ifirst = 1
        if ((ftyp.le.0).or.(ftyp.gt.3)) then 
          write(*,*) 'Type of target out of range: ftyp = ',ftyp
          return
        end if 
      end if
*  
      kint = kint1
      if(IABS(kint1).EQ.4.OR.IABS(kint1).EQ.5) kint=2*kint1/IABS(kint1)
* 
      nxbb=nxb/2
      x1=0.3d0
      xlog1=dlog(x1)
      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
      DELX1=(1.d0-x1)**2/dble(nxbb+1)
* 
      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)
*  
*...Load tables and grid 
* 
      if (kords.eq.kord) goto 10 
      ityps=ityp
      kints=kint
      kords=kord
*  
*...X GRID 
      do kx=1,nxbb
        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
      end do
      do kx=nxbb+1,nxb-1
        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
      end do
      xx(nxb)=1.d0
* 
*...Q2 GRID 
      do kq=1,nq
        qq(kq)=0.04d0*dexp(dexp(dlog(dlog(qsqmin/0.04d0))
     +        +dels*dble(kq-1)))
      end do
*  
      print *,'***** Reading Fe tables *****'
      call GETENV('cteqx_nucl',locfile)
      call trmstr(locfile,len)
      locfile =  locfile(1:len)//'/nuke_Fe.dat'
      print *,'Opening the file ',locfile
      open(unit=nport,status='old',err=200,file=locfile)
* 
      tt = 2
      do n=1,nxb-1
        do m=1,nq
          read(nport,*) (f_nuke(n,m,i,tt),i=1,np)
        end do
      end do
      close(unit=nport)
* 
      do i=0,np
        do m=1,nq
          if (i.ne.0) then
            f_nuke(nxb,m,i,tt)=0d0
          else
            f_nuke(nxb,m,i,tt)=f_nuke(nxb-1,m,i,tt)
          end if
          do n=1,nxb
            fsp(n)=f_nuke(n,m,i,tt)
          end do
          call spline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            b_nuke(n,m,i,tt)=bs(n)
            c_nuke(n,m,i,tt)=cs(n)
            d_nuke(n,m,i,tt)=ds(n)
          end do
        end do
      end do
* 
  111 format (27f12.6)
* 
  10  continue
* 
      if((q2.lt.qsqmin).or.(q2.gt.qsqmax)) then
         print 99,q2,qsqmin,qsqmax
         return
      end if
      if((xb.lt.xmin).or.(xb.gt.xmax)) then
         print 98,xb,xmin,xmax
         return
      end if
  99  format(' NUKE_FAST WARNING:  Q^2 VALUE IS OUT OF RANGE   ',3g12.3)
  98  format(' NUKE_FAST WARNING:   X  VALUE IS OUT OF RANGE   ',3g12.3)
* 
*...Now actual interpolation
* 
      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)
* 
      if (x.gt.x1) then
        xd=(1d0-x1)**2-(1d0-x)**2
        n=int(xd/delx1)+nxbb
      else
        xd=dlog(x)-xlog1
        n=nxbb+int(xd/DELX)-1
      end if
      aa=x-xx(n)
* 
      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0
* 
      kint0=kint
      if(kint.eq.-2) kint0=6
      if(kint.eq.-3) kint0=7
      if(kint.eq.-4) kint0=8
      if(kint.eq.-5) kint0=9

      kk = 27*(nsf-1) + 9*(ityp+1) + kint0

      if ((f_nuke(n,m,kk,ftyp).ne.0d0).or.
     +    (f_nuke(n+1,m,kk,ftyp).ne.0d0)) then
        f0 = f_nuke(n,m,kk,ftyp) + aa*b_nuke(n,m,kk,ftyp) 
     +     + aa**2*c_nuke(n,m,kk,ftyp) + aa**3*d_nuke(n,m,kk,ftyp)
        fp = f_nuke(n,m+1,kk,ftyp) + aa*b_nuke(n,m+1,kk,ftyp)
     _     + aa**2*c_nuke(n,m+1,kk,ftyp) + aa**3*d_nuke(n,m+1,kk,ftyp)
        if (m.ge.2) then
          fm = f_nuke(n,m-1,kk,ftyp) + aa*b_nuke(n,m-1,kk,ftyp) 
     +       + aa**2*c_nuke(n,m-1,kk,ftyp) + aa**3*d_nuke(n,m-1,kk,ftyp)
          nuke_cteq = fm*b*(b-1d0)/2d0 +f0*(1d0-b**2) +fp*b*(b+1d0)/2d0
        else
          nuke_cteq = f0*(1d0-b) + fp*b
        end if
      else 
        nuke_cteq = 0.d0
      end if
*. 
      if (nuke_cteq.lt.0.d0) nuke_cteq = 0.d0
* 
      return 
 200  print *,'The Nuke set is unavailable (FILE: ',locfile,')'

      end 


**************************************************************************************
      SUBROUTINE SPLINE(N,X,Y,B,C,D)
* ---------------------------------------------------------------------
* CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
* INTERPOLATION SUBROUTINES ARE TAKEN FROM
* G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
* COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION X(N), Y(N), B(N), C(N), D(N)
*
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 210 K=2,NM1
         D(K)=X(K+1)-X(K)
         B(K)=2.0D0*(D(K-1)+D(K))
         C(K+1)=(Y(K+1)-Y(K))/D(K)
         C(K)=C(K+1)-C(K)
  210 CONTINUE
      B(1)=-D(1)
      B(N)=-D(N-1)
      C(1)=0.0D0
      C(N)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1)=C(1)*D(1)**2.0D0/(X(4)-X(1))
      C(N)=-C(N)*D(N-1)**2.0D0/(X(N)-X(N-3))
 215  CONTINUE
      DO 220 K=2,N
         T=D(K-1)/B(K-1)
         B(K)=B(K)-T*D(K-1)
         C(K)=C(K)-T*C(K-1)
 220  CONTINUE
      C(N)=C(N)/B(N)
      DO 230 IB=1,NM1
         K=N-IB
         C(K)=(C(K)-D(K)*C(K+1))/B(K)
 230  CONTINUE
      B(N)=(Y(N)-Y(NM1))/D(NM1)
     1     +D(NM1)*(C(NM1)+2.0D0*C(N))
      DO 240 K=1,NM1
         B(K)=(Y(K+1)-Y(K))/D(K)
     1        -D(K)*(C(K+1)+2.0D0*C(K))
         D(K)=(C(K+1)-C(K))/D(K)
         C(K)=3.0D0*C(K)
 240  CONTINUE
      C(N)=3.0D0*C(N)
      D(N)=D(N-1)
      RETURN
 250  CONTINUE
      B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0D0
      D(1)=0.0D0
      B(2)=B(1)
      C(2)=0.0D0
      D(2)=0.0D0
      RETURN
      END


*-----------------------------------------------------------

      function emcslac(f,x,q2)

c                     The function below "EMC(x, iver)" gives the
c	Fe/D correction for the "x" bin. Note that the
c	.07< x < 1 is what you must use this function for.
	
c                   For x=.045, the Fe/D is 0.95. It is the
c	"shadowing" region and the factor was obtained for
c	Ca/D, and correcting Ca to Fe using a prescription
c	by Strikman.


c two versions of SLAC fits for the F2(Fe)/F2(D2) ratio

      implicit double precision (a-z)

c version 2 is the fit to SLAC E-139 and E-140 and STEIN

	data d0, d1, d2, d3, d4, d5, d6, d7, d8
     $  /4.58558707D-01,
     $   1.62185596D+01,
     $  -1.79392859D+02,
     $   1.04313998D+03,
     $  -3.53408342D+03,
     $   7.17002801D+03,
     $  -8.56431003D+03,
     $   5.54039709D+03,
     $  -1.49318167D+03/

      x2 = x*x
      x4 = x2*x2
      x8 = x4*x4
      
      if(x.ge..07) then
         emc1 = d0 + d1*x + d2*x2 + d3*x*x2 + d4*x4
     $        +  d5*x4*x + d6*x4*x2 + d7*x4*x2*x + d8*x8
      else
c     
c the following form is a linear interpolation between the x = .07
c value of the above polynomial and the Strikman nuclear shadowing
c calculation - emc1 = .95 @ x = .045
c
         emc1 = .95 + (x - .045)*1.9544

      endif
      emcslac = f / emc1
      
      return
      end


      FUNCTION F2DTOA(F,X,Q2)
 

      Implicit Double Precision (A-H, O-Z)
 

C *** GIVEN F2 MEASURED ON FE OR C AT A GIVEN X IT WILL GIVE THE
C *** EQUIVALENT DEUTERIUM VALUE.  IT WILL ALSO CALCULATE THE POSSIBLE
C *** Q**2 DEPEndENCE OF THE EFFECT FOR IRON IF THE VALUE OF Q GIVEN TO
C *** FUNCTION IS GREATER THAN SQRT(5).
 

C     DATA Q20 / 5.0 /
 

      X2 = X * X
      X4 = X2 * X2
 

      DF2 = 1.18 - 4.02*X + 24.35*X2 - 61.5*X2*X + 72.61*X4 -32.35*X4*X
      F2D = DF2 * F
 

C      IF(Q2 .LE. Q20) GO TO 20
C      DQ2 = Q2 / 5.0
C      ALDQ2 = LOG(DQ2)
C      F2D = F2D + (0.077 - 0.244*X) * ALDQ2
C  20  Continue
 

      F2DTOA = F2D
      Return
C                        ****************************
      End

      function emcnmc(f,x,q2)

c Concerning the "heavy target correction" I used a parametrisation of
c our measurement of Ca/D2 and the Slac result for Fe/D2.
c As far as I can see it,this should be also valid for the iron data
c of CCFR (from our preliminay data on Sn/C we conclude that the
c amount of shadowing seems to saturate at values of A around 40).
c I include this parametrisation into this mail.


      implicit double precision(a-z)

      emc = 1.118 - 0.4199*x - 0.3597*exp(-22.88*x) +
     &        1.872*(x**11.27)

      emcnmc = f /  emc

      return
      end


      function emce665(f,x,q2)

      implicit double precision (a-z)

c Data from E-665 (Xenon) and NMC (Ca scaled) used for large x fit
c and E-665 for very small x.  The small x fit went to x = 10**-5

      x2 = x*x
      x4 = x2*x2

      if(x.ge..01) then
         rat = 0.73 + 7.0*x - 52.9*x2 + 173.8*x2*x - 263.9*x4 +
     $        150.*x4*x
      else
         rat = 0.71 + 6.8*x - 31.1*x2
      endif

      emce665 = f / rat

      return

      end




