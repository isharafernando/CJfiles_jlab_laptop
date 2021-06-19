****************************************************************************
*     Structure function calcs and output 
*
*     Programmers:  aa - A.Accardi (accard@jlab.org)
*
*     HISTORY:
*
*     sfn.f: :    -- the first one-- (aa)
*     (27 Apr 18) 
*
*****************************************************************************


**************************************************************************
      subroutine FFcalc(x,Q2,par,ibeam,istruc,itm,iht,inuke,ffp,ffn,ffd)
*     
*     Calculates p,n,d structure functions at once
*
*     INPUT:
*
*     x        = (dp) Bjorken x
*     Q2       = (dp) [GeV^2] photon virtuality squared 
*     par(100) = (dp) PDF parameters
*     ibeam    = (i)  beam type: 1=e 2=nu 3=nubar 
*     istruc   = (i)  structure fn: 0=FL 1=F1 2=F2 3=xF3  
*     itm      = (i)  TMC  0=no  1=approx GP  2=CF
*     iht      = (i)  HT   0=no  1=C=C(xB)  2=C=C(xi_Nachtmann)
*     inuke    = (i)  nuclear correction flag
*      
*     OUTPUT:
*
*     ffp    = proton structure function
*     ffn    = neutron structure function
*     ffd    = deuteron structure function
*

      implicit none

      double precision x,Q2,par(100),ffp,ffn,ffd
      integer ibeam,istruc,inuke,itm,iht

      double precision v(4)

      v(1) = x
      v(2) = Q2
      ! F2 proton
      call targ(1d0,0d0)        ! target proton, neutron fractions: proton
                                ! for free neutron use (0d0,1d0)  
                                ! for deuteron=p+n use (1d0,1d0)
                                ! for a nucleon=(p+n)/2 use (0.5d0,0.5d0)
                                ! No call to 'targ' needed with 'deuteronsf'
      call strfn(ibeam,istruc,itm,iht,0,v,par,ffp)
      ! F2 neutron
      call targ(0d0,1d0)        ! target Z,A numbers: neutron
      call strfn(ibeam,istruc,itm,iht,0,v,par,ffn)
      ! F2 Deuteron (automatically sets A,Z)
      call deuteronsf(ibeam,istruc,inuke,itm,iht,0,v,par,ffd)

      return
      end

      
**************************************************************************
      subroutine writetbl_sfn(tblfile,ibeam,istruc)
C     Tabulates and writes to file the CJ DIS structure functions in
C     the "tbl" format used by the CTEQ6 series, so that they can be
C     used in CTEQ functions and in the CJ grid interpolator. 
C
C     The resultant grid files allow also to use the CJ PDFs in the 
C     MCFM monte carlo.  All you have to do is to call the cteq functions 
C     in mcfm to read the obtained pds file.
C
C     NOTE: It sssumes that a CJ parameter files has already been read 
C     and evolved, e.g., 
C           call readpar(parfile,par,inuke,itmc,iht)
C           call QCDev      
C
C     INPUT:
C 
C       tblfile = (ch*50) pds file name to be written
C       ibeam    = (i)  beam type: 1=e 2=nu 3=nubar 
C       istruc   = (i)  structure fn: 0=FL 1=F1 2=F2 3=xF3  
C     
C     Programmers: Alberto Accardi
C     Date: Apr 2018
C

      implicit none

      character*50 tblfile
      integer ibeam,istruc
      
      integer iQ,ix,i,j,outn,NextUn,nu2
      real*8 xb,Q,Q2,ffp,ffn,ffd,ffplt,ffnlt,ffdlt,ffp0,ffn0,ffd0
      real*8 ord,nfl,lambda,qm1,qm2,qm3,qm4,qm5,qm6
      integer nfmx
      real*8 qini,qmax,xmin
      
*    *** CJ parameters
      double precision par(100)
      common/curpar/par
      integer ILOOP,IORD,NMAX,IVL,ncb
      double precision xmc,xmb
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      double precision FLAVOR
      COMMON/GRPTHY/FLAVOR

*    *** Theory correction flag
      integer inuke,itm,iht
      common/thcor/inuke,itm,iht

*    *** x and Q grids
      integer nQ,nx
      parameter(nx=96) 
      parameter(nQ=26)
      double precision fullpdf(9*nx*nQ)      
      double precision Qgr(nQ), xgr(nx)
      data Qgr/   1.30000D+00, 1.50159D+00, 1.75516D+00, 2.07810D+00,
     &            2.49495D+00, 3.04086D+00, 3.76715D+00, 4.50000D+00,
     &            4.75000D+00,
     &            6.23113D+00, 8.37423D+00, 1.15549D+01, 1.64076D+01,
     &            2.40380D+01, 3.64361D+01, 5.73145D+01, 9.38707D+01,
     &            1.60654D+02, 2.88438D+02, 5.45587D+02, 1.09231D+03,
     &            2.32646D+03, 5.30043D+03, 1.29956D+04, 3.45140D+04,
     &            1.00000D+05 /
      data xgr/   0.00000D+00, 1.00000D-06, 1.28121D-06, 1.64152D-06,
     &            2.10317D-06, 2.69463D-06, 3.45242D-06, 4.42329D-06,
     &            5.66715D-06, 7.26076D-06, 9.30241D-06, 1.19180D-05,
     &            1.52689D-05, 1.95617D-05, 2.50609D-05, 3.21053D-05,
     &            4.11287D-05, 5.26863D-05, 6.74889D-05, 8.64459D-05,
     &            1.10720D-04, 1.41800D-04, 1.81585D-04, 2.32503D-04,
     &            2.97652D-04, 3.80981D-04, 4.87518D-04, 6.26039D-04,
     &            8.00452D-04, 1.02297D-03, 1.30657D-03, 1.66759D-03,
     &            2.12729D-03, 2.71054D-03, 3.44865D-03, 4.37927D-03,
     &            5.54908D-03, 7.01192D-03, 8.83064D-03, 1.10763D-02,
     &            1.38266D-02, 1.71641D-02, 2.11717D-02, 2.59364D-02,
     &            3.15062D-02, 3.79623D-02, 4.53425D-02, 5.36750D-02,
     &            6.29705D-02, 7.32221D-02, 8.44039D-02, 9.64793D-02,
     &            1.09332D-01, 1.23067D-01, 1.37507D-01, 1.52639D-01,
     &            1.68416D-01, 1.84794D-01, 2.01731D-01, 2.19016D-01,
     &            2.36948D-01, 2.55242D-01, 2.73927D-01, 2.92954D-01,
     &            3.12340D-01, 3.32036D-01, 3.52019D-01, 3.72282D-01,
     &            3.92772D-01, 4.13533D-01, 4.34326D-01, 4.55495D-01,
     &            4.76836D-01, 4.98342D-01, 5.20006D-01, 5.41818D-01,
     &            5.63773D-01, 5.85861D-01, 6.08077D-01, 6.30459D-01,
     &            6.52800D-01, 6.75387D-01, 6.98063D-01, 7.20830D-01,
     &            7.43683D-01, 7.66623D-01, 7.89636D-01, 8.12791D-01,
     &            8.35940D-01, 8.59175D-01, 8.82485D-01, 9.05866D-01,
     &            9.29311D-01, 9.52817D-01, 9.76387D-01, 1.00000D+00 /


      ! CJ fit parameters
      ord=iord+1
      nfl=flavor
      qm1=0d0   
      qm2=0d0
      qm3=0d0    
      qm4=xmc
      qm5=xmb
      qm6=180d0
      lambda=par(1)

      ! No. of str.fns in output (minus 3 b/c of grid interpolator routine)
      ! 9 == (F2p,n,d for full,noHT,noHTnoTMC)
      nfmx=9-3

      qini=Qgr(1)
      qmax=Qgr(nQ)
      xmin=xgr(1)

      ! file headers
      outn=NextUn()
      open (outn,FILE=tblfile,status='new',err=200)
      
      write(outn,*) " CJ structure function tables"
      write(outn,*) "  Ordr, Nfl, lambda        Qmass"//
     &              " 1,  2,  3,         4,  5,  6"
      write(outn,22) ord,nfl,lambda,qm1,qm2,qm3,qm4,qm5,qm6
      write(outn,*) "   NX,  NQ,  Nsfns-3"
      write(outn,23) nx,nq,nfmx  ! 6 = 9-3 structure functions 
      write(outn,*) "QINI, QMAX, (QV(I), I =1, NQ)"
      write(outn,24) qini,qmax
      
      write(outn,24) qgr

      write(outn,*)"XMIN, (XV(I), I =1, NX)"
      write(outn,24)xmin

      write(outn,24) xgr
      write(outn,*)"Structure Function Table:"
      
      ! calculates grid
      do iQ=1,nQ
          Q=Qgr(iQ)
          q2=q*q
          do ix=1,nx
              xb=xgr(ix)

              ! structure function

              itm=1
              iht=1
              call FFcalc(xb,Q2,par,ibeam,istruc,itm,iht,inuke
     &           ,ffp,ffn,ffd)
              itm=1
              iht=0
              call FFcalc(xb,Q2,par,ibeam,istruc,itm,iht,inuke
     &           ,ffplt,ffnlt,ffdlt)
              itm=0
              iht=0
              call FFcalc(xb,Q2,par,ibeam,istruc,itm,iht,inuke
     &           ,ffp0,ffn0,ffd0)

              i=(0*nQ+(iQ-1))*nx+ix
              
              fullpdf(i)=ffp
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(1*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffn
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(2*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffd
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(3*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffplt
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(4*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffnlt
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(5*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffdlt
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(6*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffp0
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(7*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffn0
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(8*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ffd0
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
          enddo
          
      enddo

      ! writes to file
      write(outn,25) fullpdf
      
      close (outn)     
      
 22   format(F6.0,F6.0,F8.4,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)
 23   format(I6,I5,I5)
 24   format(6es12.5)
 25   format(5es12.5)
      return
      
 200  print*, 'WARNING(writetbl_sfn): tblfile already exists: "'
     &     //trim(tblfile)//'"'
      return
      end


**************************************************************************
      subroutine writedat_sfn(datfile,type,flav,ibeam)  
C     Tabulates and writes to file the CJ pdfs in LHAPDF6 "dat" format 
C     see http://lhapdf.hepforge.org/design.html
C
C     NOTE: It sssumes that a CJ parameter files has already been read 
C     and evolved, e.g., 
C           call readpar(parfile,par,inuke,itmc,iht)
C           call QCDev      
C
C     INPUT:
C 
C       datfile = (ch*50) pds file name to be written
C       type = (ch*50) should be 'central' or 'error'
C       flav = PDG flavor codes for output files
C       ibeam = 1 for electron NC scattering; 2,3 for CC scattering
C      
C     Programmers: Alberto Accardi
C     Date: May 2018
C     

      implicit none

      character*50 datfile(3),type
      character*1000 flav
      integer ibeam
      
      integer iQ,ix,i,j,idx,outp,outn,outd,NextUn,nu2,neff
     &     ,istruc,itm0,iht0
      real*8 xb,Q,Q2,f2p,f2n,f2d,f2plt,f2nlt,f2dlt,f2p0,f2n0,f2d0
      real*8 lambda

      double precision alphas_sa
      
*    *** CJ parameters
      double precision par(100)
      common/curpar/par
      integer ILOOP,IORD,NMAX,IVL,ncb
      double precision xmc,xmb
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      double precision FLAVOR
      COMMON/GRPTHY/FLAVOR

*     Heavy quark treatment
      integer ihq
      double precision a0u,a0d
      common/sacot/ihq,a0u,a0d
      
*    *** Theory correction flag
      integer inuke,itm,iht
      common/thcor/inuke,itm,iht

*    *** x and Q grids
      integer nQ,nx,nsf
      parameter(nx=95) 
      parameter(nQ=26)
      parameter(nsf=9)
      double precision fullSFp(nsf*nx*nQ)
     &     ,fullSFn(nsf*nx*nQ),fullSFd(nsf*nx*nQ)
      double precision Qgr(nQ), xgr(nx)
      data Qgr/   
     &            1.30000D+00, 1.50159D+00, 1.75516D+00, 2.07810D+00,
     &            2.49495D+00, 3.04086D+00, 3.76715D+00, 4.50000D+00,
     &            4.75000D+00,
     &            6.23113D+00, 8.37423D+00, 1.15549D+01, 1.64076D+01,
     &            2.40380D+01, 3.64361D+01, 5.73145D+01, 9.38707D+01,
     &            1.60654D+02, 2.88438D+02, 5.45587D+02, 1.09231D+03,
     &            2.32646D+03, 5.30043D+03, 1.29956D+04, 3.45140D+04,
     &            1.00000D+05 /
      data xgr/   1.00000D-06, 1.28121D-06, 1.64152D-06,
     &            2.10317D-06, 2.69463D-06, 3.45242D-06, 4.42329D-06,
     &            5.66715D-06, 7.26076D-06, 9.30241D-06, 1.19180D-05,
     &            1.52689D-05, 1.95617D-05, 2.50609D-05, 3.21053D-05,
     &            4.11287D-05, 5.26863D-05, 6.74889D-05, 8.64459D-05,
     &            1.10720D-04, 1.41800D-04, 1.81585D-04, 2.32503D-04,
     &            2.97652D-04, 3.80981D-04, 4.87518D-04, 6.26039D-04,
     &            8.00452D-04, 1.02297D-03, 1.30657D-03, 1.66759D-03,
     &            2.12729D-03, 2.71054D-03, 3.44865D-03, 4.37927D-03,
     &            5.54908D-03, 7.01192D-03, 8.83064D-03, 1.10763D-02,
     &            1.38266D-02, 1.71641D-02, 2.11717D-02, 2.59364D-02,
     &            3.15062D-02, 3.79623D-02, 4.53425D-02, 5.36750D-02,
     &            6.29705D-02, 7.32221D-02, 8.44039D-02, 9.64793D-02,
     &            1.09332D-01, 1.23067D-01, 1.37507D-01, 1.52639D-01,
     &            1.68416D-01, 1.84794D-01, 2.01731D-01, 2.19016D-01,
     &            2.36948D-01, 2.55242D-01, 2.73927D-01, 2.92954D-01,
     &            3.12340D-01, 3.32036D-01, 3.52019D-01, 3.72282D-01,
     &            3.92772D-01, 4.13533D-01, 4.34326D-01, 4.55495D-01,
     &            4.76836D-01, 4.98342D-01, 5.20006D-01, 5.41818D-01,
     &            5.63773D-01, 5.85861D-01, 6.08077D-01, 6.30459D-01,
     &            6.52800D-01, 6.75387D-01, 6.98063D-01, 7.20830D-01,
     &            7.43683D-01, 7.66623D-01, 7.89636D-01, 8.12791D-01,
     &            8.35940D-01, 8.59175D-01, 8.82485D-01, 9.05866D-01,
     &            9.29311D-01, 9.52817D-01, 9.76387D-01, 1.00000D+00 /

*    *** Info for .info file

      double precision xmin,xmax,Qmin,Qmax,MZ
      common/lhainfo1/xmin,xmax,Qmin,Qmax,MZ

      double precision alphas_MZ,alphas_L4,alphas_L5
      common/lhainfo2/alphas_MZ,alphas_L4,alphas_L5

      integer nalpha
      double precision alphas_Q(100),alphas_v(100)
      common/lhainfo3/alphas_Q,alphas_v,nalpha

      integer nfl,ord,hq
      common/lhainfo4/nfl,ord,hq

      double precision qm1,qm2,qm3,qm4,qm5,qm6
      common/lhainfo5/qm1,qm2,qm3,qm4,qm5,qm6

      
      lambda=par(1)
      
      i = 0
      do ix=1,nx
         xb=xgr(ix)
         do iQ=1,nQ
            Q=Qgr(iQ)
            Q2 = Q*Q
            
              ! structure function

              itm0=1
              iht0=1
              istruc = 2
              call FFcalc(xb,Q2,par,ibeam,istruc,itm0,iht0,inuke
     &           ,f2p,f2n,f2d)
              itm0=1
              iht0=0
              istruc = 2
              call FFcalc(xb,Q2,par,ibeam,istruc,itm0,iht0,inuke
     &           ,f2plt,f2nlt,f2dlt)
              itm0=0
              iht0=0
              istruc = 2
              call FFcalc(xb,Q2,par,ibeam,istruc,itm0,iht0,inuke
     &           ,f2p0,f2n0,f2d0)

              ! fills the lhapdf grid
              ! (if negative writes 0d0)

              ! proton
              fullsfp(i+1)  = dmax1(f2p  ,0d0) ! F2
              fullsfp(i+2)  = dmax1(f2plt,0d0)
              fullsfp(i+3)  = dmax1(f2p0 ,0d0)
              do idx = 4,9 ! FL and F3 = 0 for now
                 fullsfp(i+idx) = 0d0
              end do
              ! neutron
              fullsfn(i+1)  = dmax1(f2n  ,0d0) ! F2
              fullsfn(i+2)  = dmax1(f2nlt,0d0)
              fullsfn(i+3)  = dmax1(f2n0 ,0d0)
              do idx = 4,9 ! FL and F3 = 0 for now
                 fullsfn(i+idx) = 0d0
              end do
              ! deuteron
              fullsfd(i+1)  = dmax1(f2d  ,0d0) ! F2
              fullsfd(i+2)  = dmax1(f2dlt,0d0)
              fullsfd(i+3)  = dmax1(f2d0 ,0d0)
              do idx = 4,9 ! FL and F3 = 0 for now
                 fullsfd(i+idx) = 0d0
              end do

c$$$              if (Q2.lt.10d0.and.xb.gt.0.1.and.xb.lt.0.5) then
c$$$                 print*, '* ix,iQ,i =',ix,iq,i
c$$$                 print*, '* xb,Q =', xb, Q2
c$$$                 print*, '* F2   =',f2p,f2n,f2d
c$$$                 print*, '* F2lt =',f2plt,f2nlt,f2dlt
c$$$                 print*, '* F20  =',f2p0,f2n0,f2d0
c$$$                 print*, fullsfp(i+1:i+9)
c$$$                 print*, fullsfn(i+1:i+9)
c$$$                 print*, fullsfd(i+1:i+9)
c$$$              end if 
              
              ! sets the counter for next row
              i = i+9


          enddo
          
      enddo

*    ... collects information for .info file

      ord=iord
      nfl=flavor
      hq = ihq
      
      qm1=0d0   
      qm2=0d0
      qm3=0d0    
      qm4=xmc
      qm5=xmb
      qm6=180d0

      xmin = xgr(1)
      xmax = xgr(nx)
      Qmin = Qgr(1)
      Qmax = Qgr(nQ)
      MZ = 91.1876
      alphas_MZ = alphas_sa(MZ**2,iord+1,neff)
      alphas_L4 = alphas_sa(xmc**2,iord+1,neff)
      alphas_L5 = alphas_sa(xmb**2,iord+1,neff)
      nalpha=nq
      do iQ = 1, nq
         alphas_Q(iq) = Qgr(iq)
         alphas_v(iq) = alphas_sa(Qgr(iq)**2,iord+1,neff)
      end do

      ! writes to file   
       
      outp=NextUn()
      outn=NextUn()
      outd=NextUn()

      open (outp,FILE=datfile(1),status='new',err=200) ! proton
      write(outp,'(A)') "PdfType: "//trim(type)
      write(outp,'(A)') "Format: lhagrid1"
      write(outp,'(A)') "---"
      write(outp,31) xgr
      write(outp,31) Qgr
      write(outp,'(A)') trim(flav)
      write(outp,32) fullsfp
      write(outp,'(A)') "---"
      goto 201
 200  print*, 'WARNING(writedat_sfn): datfile already exists: "'
     &     //trim(datfile(1))//'"'
      
 201  open (outn,FILE=datfile(2),status='new',err=202) ! neutron
      write(outn,'(A)') "PdfType: "//trim(type)
      write(outn,'(A)') "Format: lhagrid1"
      write(outn,'(A)') "---"
      write(outn,31) xgr
      write(outn,31) Qgr
      write(outn,'(A)') trim(flav)
      write(outn,32) fullsfn
      write(outn,'(A)') "---"
      goto 203
 202  print*, 'WARNING(writedat_sfn): datfile already exists: "'
     &     //trim(datfile(2))//'"'


 203  open (outn,FILE=datfile(3),status='new',err=204) ! deuteron
      write(outd,'(A)') "PdfType: "//trim(type)
      write(outd,'(A)') "Format: lhagrid1"
      write(outd,'(A)') "---"
      write(outd,31) xgr
      write(outd,31) Qgr
      write(outd,'(A)') trim(flav)
      write(outd,32) fullsfd
      write(outd,'(A)') "---"
      goto 205
 204  print*, 'WARNING(writedat_sfn): datfile already exists: "'
     &     //trim(datfile(3))//'"'
      
 205  close (outp)
      close (outn)
      close (outd)

 31   format(es11.5,200es12.5)
 32   format(es11.5,8es12.5)

      return
      end
