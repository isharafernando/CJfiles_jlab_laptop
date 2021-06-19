****************************************************************************
*     PDF evolution and interpolation 
*        -- stand-alone package -- 
*
*     Programmers:  aa - A.Accardi (accard@jlab.org)
*                  jfo - J.Owens
*
*     HISTORY:
*
*     pdf10.f:    -- the first one-- (aa)
*     (11 Mar 10) Builds upon 'dist10.f' and 'fillgf10.f' by jfo. 
*                 Contains stand-alone versions of 'pdf' which interpolates 
*                 array 'gf' ('calc' in the fitting code)
*                 Still needs to be linked to 'altpar10' for the actual 
*                 evolution routines.
*
*****************************************************************************


      subroutine readpar(parfile,params,inuke,itmc,iht)
*     Reads PDf parameter from file
*     INPUT:  parfile   = (ch*50)  PDF parameters file 
*     OUTPUT: params    = (dp*100) the PDF parameters 
*             inuke     = nuclear corretion code
*             itmc      = TMC code
*             iht       = HT code

      implicit none
      character*50 parfile,dummy
      integer nu,NextUn,i,j,inuke,iht,itmc,ios,npdfpar
      double precision params(100),epsilon(50,50)
*     QCD evolution parameters
      double precision Q02,Q2MAX
      COMMON/Q2STUF/ Q02,Q2MAX
        ! Q02 must be set up here; Q2max will is set by 'QCDev'
      integer ILOOP,IORD,NMAX,IVL,ncb
      double precision xmc,xmb
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
        ! NMAX set by 'QCDev'
        ! IVL does not need to be set here
      double precision FLAVOR
      COMMON/GRPTHY/FLAVOR
      ! if using this with jet routines 'jet.f', needs also the following
      ! common/qcdpar/al,nf,norder,set ! needed by 'jets.f'
      ! and setting norder=iord+1
*     PDF params
      double precision par(100)
      common/param/par
      integer iparamc
      common/iparam/iparamc
      integer npar              ! stores the total number of parameters
      common/npar/npar
*     Other PDF parameters      ! Needed when using 'writepar'
      double precision uncrt(100)
      character*10 pname(100)
      common/oparam/pname,uncrt
      integer INS,NVL,iserr
      common/oparam1/INS,NVL,iserr


*     alpha matrix and paramter errors
      integer pmax
      parameter(pmax=50) ! must be = dimension of alpha in 'minim'
      double precision alfa(pmax,pmax),central(pmax),eps(pmax)
     &     ,cov(pmax,pmax)
      integer pos(pmax),nfree
*      common/PDFerr/alfa,cov,central,eps,pos,nfree
      common/PDFerr/alfa,cov,central,eps,pos,npdfpar
*     comment lines in .par file
      character*100 line(3)
      common/parfil/line


      nu=NextUn()
      open(unit=nu,file=parfile,status='old')
*    *** reads QCD parameters from .par file
      read(nu,'(A)') line(1)
      read(nu,'(A)') line(2)
      read(nu,'(A)') line(3)
      print*, line
      read(nu,*) dummy
      read(nu,*) dummy,INS
      read(nu,*) dummy,NVL
      read(nu,*) dummy,IORD
      ILOOP = iord + 1
      read(nu,*) dummy,Q02
c jfo      read(nu,*) dummy,ISERR
      read(nu,*) dummy,IPARAMC
      read(nu,*) dummy,FLAVOR
      read(nu,*) dummy,xmc
      read(nu,*) dummy,xmb
      read(nu,*) dummy,ncb
 
*    *** Reads theory corrections codes
      read(nu,*) dummy
      read(nu,*) dummy,inuke
      read(nu,*) dummy,itmc
      read(nu,*) dummy,iht

*    *** reads PDF params from file
      DO J=1,100
         PAR(J)=0.
         uncrt(j)=0.
      end do
      read(nu,*) dummy
      nfree = 0
      do 31 j=1,99       ! only up to 99 because pname(100) is reserved
                         ! for datasets that do not require a normalization
         read(nu,*,end=32,err=32) pname(j),par(j),uncrt(j)
         params(j) = par(j)     ! par(j) ends up in common/param/ for use 
                                ! by fitpack, such as in 'altpar.f'
                                ! params(j) returned to thr user who does
                                ! not need to know the previous details
         if (uncrt(j).gt.0d0) then ! fills in arrays for error evaluation
            nfree = nfree + 1       ! counts number of free parameters
            central(nfree) = par(j) ! parameter central value
!            eps(nfree) = uncrt(j)
            pos(nfree) = j          ! position of central par in 'par' array
         end if
 31   continue
 32   npar=j-2                  ! total number of parameters read in
*    ...looks for higher-twist parameters
      call findht(npar,par,pname,uncrt)

c jfo Count number of PDF parameters

      npdfpar=0
      do j=1,35
         if(uncrt(j).ne.0.)npdfpar=npdfpar+1
      enddo

*    *** reads in the alpha matrix, calculates errors
c jfo      do i=1,nfree
c jfo         read(nu,*,iostat=ios,end=33,err=33)(alfa(i,j),j=1,nfree) 
      do i=1,npdfpar
         read(nu,*,iostat=ios,end=33,err=33)(alfa(i,j),j=1,npdfpar) 
      end do
 33   if (ios.ne.0) then
         write(*,*) 'WARNING(readpdf): missing or wrong error matrix'
      end if
      
c jfo      call invmatrix(alfa,nfree,pmax,epsilon) ! inverts alfa 
c jfo      do i=1,nfree
      call invmatrix(alfa,npdfpar,pmax,epsilon) ! inverts alfa 
      do i=1,npdfpar
         eps(i) = dsqrt(epsilon(i,i)) 
         if ((eps(i)/uncrt(pos(i))-1d0).gt.1d-3) then
            print*, 'WARNING(readpar): i=',i,' eps=',eps(i)
     &           ,' uncrt=', uncrt(pos(i))
         end if
      end do
 1000 format('WARNING(readpar): i=',i4,'  eps=',E12.5,'  uncrt=',E12.5)

*    *** reads in the cov matrix
      read(nu,*,end=34,err=34) dummy
      do i=1,nfree
         read(nu,*,iostat=ios,end=34,err=34)(cov(i,j),j=1,nfree) 
      end do
 34   if (ios.ne.0) then
         write(*,*) 'WARNING(readpara): missing or wrong cov matrix'
      end if

      close(nu) 

      return 
      end



      subroutine writepar(parfile,line,pname,params,uncrt,alfa,nfree
     &     ,inuke,itmc,iht)
*     Writes PDF parameter file
*     INPUT:  parfile   = (ch*50)  file name
*             line      = (ch*100*3) comment lines
*             pname     = (ch*10*100) parameter name
*             params    = (dp*100) the PDF parameters 
*             uncrt     = (dp*100) the parameter error
*             alfa      = (dp*100,100) the 'alpha matrix' of Bevington
*             nfree     = (i) physical dimension of alfa (no. of free params.)
*             inuke     = nuclear correction code
*             itmc      = TMC code
*             iht       = HT code
*     All the rest is passed in common

      implicit none

      character*50 parfile
      character*100 line(3)

      integer inuke,itmc,iht,nfree
     &     ,nu,NextUn,len,i,j

      double precision params(100),uncrt(100),alfa(100,100)
      character*10 pname(100)


*     QCD evolution parameters
      double precision Q02,Q2MAX
      COMMON/Q2STUF/ Q02,Q2MAX
      integer ILOOP,IORD,NMAX,IVL,ncb
      double precision xmc,xmb
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      double precision FLAVOR
      COMMON/GRPTHY/FLAVOR
      integer INS,NVL,iserr
      common/oparam1/INS,NVL,iserr

*     PDF params
      integer iparamc
      common/iparam/iparamc
      integer npar              
      common/npar/npar


      nu = NextUn()
      call trmstr(parfile,len)
      open(unit=nu,file=parfile(1:len),status='new',err=200)

      call trmstr(line(1),len)
      write(nu,'(A)') line(1)(1:len)
      call trmstr(line(2),len)
      write(nu,'(A)') line(2)(1:len)
      call trmstr(line(3),len)
      write(nu,'(A)') line(3)(1:len)
      
      write(nu,*)
      write(nu,*) '# QCD setup'
      write(nu,*) ' INS	               ',INS
      write(nu,*) ' NVL                ',NVL
      write(nu,*) ' IORD               ',IORD
      write(nu,*) ' Q02               ',Q02
      write(nu,*) ' ISERR              ',ISERR
      write(nu,*) ' IPARAMC            ',IPARAMC
      write(nu,*) ' FLAVOR            ',FLAVOR
      write(nu,*) ' xmc               ',xmc
      write(nu,*) ' xmb               ',xmb
      write(nu,*) ' ncb                ',ncb
      
      write(nu,*) 
      write(nu,*) '# Theory corrections'
      write(nu,*) ' inuke              ',inuke   
      write(nu,*) ' itmc               ',itmc 
      write(nu,*) ' iht                ',iht   
      
      write(nu,*) 
      write(nu,*) '# fit parameters'
      write(nu,1111) (pname(j),params(j),uncrt(j),j=1,npar)
      write(nu,*) 'END 0. 0.'
 1111 FORMAT(1X,A10,2G18.10)


      if (alfa(1,1).ne.0d0) then ! writes alfa only if alfa(1,1)=/=0
         write(nu,*) 
         write(nu,*) '# alpha matrix'
         do I=1,nfree
            WRITE(nu,1390)(alfa(I,J),J=1,nfree) 
 1390       FORMAT(1X,100E18.10) 
            ! 14 significant figures seem a bit exxaggerated, but 
            ! are needed for precise inversion of the matrix
            ! (10 are sufficient in the test examples, but one never knows)
         end do
      end if

      goto 210
 200  print*, 'WARNING(writepar): parfile already exists: "'
     &     //parfile(1:len)//'"'
 210  close(nu)

      return
      end



      function alphas_sa(Q2,iloop,neff)
*     returns alpha_s at scale Q2, with Lambda_QCD as read in from 
*     .par or .gf file.
*    
*     INPUT:  Q2 [GeV^2]  = (dp) scale
*             iloop       = (i)  0=constant, 1=one loop, 2=two loops
*     OUTPUT: alpha_cteqx = (dp) alpha_s
*             neff        = (i)  number of active flavors at scale Q2

      implicit none

      double precision alphas_sa,ALFAS5

      integer iloop,neff
      double precision alpha_cteqx,Q2

      double precision s0,alambda5
      common/scales/s0,alambda5

      alphas_sa = ALFAS5(Q2,alambda5,iloop,neff)

      return
      end


      subroutine writegf(gffile)
*     Saves the 'gf' array created by the QCD evolution 'QCDev'
      implicit real*8 (a-h,o-z)  
      character gffile*50 
      common/param/par(100)
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/GFUNC/CALC(11,60,60) 

      nu = NextUn()
      open(nu,file=gffile,status='new',err=200)
      write(nu,11) xmc,xmb,par(1),ncb
 11   format(3f10.4,i5)
      WRITE(nu,10) (((calc(i,j,k),i=1,11),j=1,60),k=1,60)
 10   FORMAT(8E15.9)
      goto 210
 200  call trmstr(gffile,len)
      print*, 'WARNING(fillgf): gffile already exists: "'
     &     //gffile(1:len)//'"'
 210  close(nu)

      return 
      end


      subroutine readgf(gffile)
*     Reads a 'gf' array containing the evolved PDF grid 
      implicit real*8 (a-h,o-z)  
      character gffile*50 
      COMMON/GFUNC/CALC(11,60,60) 
*    ... needed for PDF interpolation
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      common/threshold/sb
      common/stepsize/delta
      common/scales/s0,alambda5
*    ... needed for F2 computations
      COMMON/Q2STUF/Q02,Q2MAX

      nu = NextUn()
      open(nu,file=gffile,status='old')
      read(nu,11) xmc,xmb,alambda5,ncb
 11   format(3f10.4,i5)
      read(nu,10) (((calc(i,j,k),i=1,11),j=1,60),k=1,60)
*      READ(nu,10) calc
 10   FORMAT(8E15.9)
      close(nu)
      call fillgf
      q02 = xmb**2
      s0=dlog(xmc**2/alambda5**2)
      sb=dlog(dlog(xmb**2/alambda5**2)/s0)
      delta=sb/ncb

      return 
      end


      subroutine writepdf(pdffile)
*     writes a .pdf sample file of PDFs, w/o errorbands
*     (but with PDFerror columns for back-compatibility

      implicit none
      
      character pdffile*50 

      integer nu,NextUn,jq2,jx,j,len,ipt,nu2

      double precision q2,x
     &     ,u0,d0,ub0,db0,sb0,cb0,bb0,glue0

      integer nx,nQ2
      double precision pdf0(-5:2,50),xgrid(50),q2grid(30)
      DATA NX,XGRID/30,1d-4,2d-4,4d-4,6d-4,8d-4,1d-3,2d-3,4d-3
     2     ,8d-3,1.6d-2,3.2d-2,6.4d-2,1d-1,1.5d-1,2d-1,2.5d-1,3d-1
     3     ,3.5d-1,4d-1,4.5d-1,5d-1,5.5d-1,6d-1,6.5d-1,7d-1,7.5d-1
     $     ,8d-1,8.5d-1,9d-1,9.5d-1,20*0d0/
      DATA NQ2,Q2GRID/29
     &     ,1.69d0,2d0,2.5d0,3.2d0,4.0d0,5.0d0,6.3d0,7.9d0,10d0
     &     ,13d0,17d0,20d0,25d0,32d0,40d0,50d0,63d0,79d0,100d0
     &     ,130d0,170d0,200d0,250d0,320d0,400d0,500d0,630d0,790d0,1000d0
     &     ,0d0 /

*     Q2, x values for PDF output
      double precision q2pdf(10),xpdf(10)
      integer nq2pdf,nxpdf
      data  q2pdf/ 1.69,2.,10.,25.,64.,100.,1328.,3*0 /
      data nq2pdf/ 7 /
      data   xpdf/ 0.1, 0.3, 0.5, 0.7, 0.85,5*0 /
      data  nxpdf/ 5 /


      nu = NextUn()
      call trmstr(pdffile,len)
      open(unit=nu,file=pdffile(1:len),status='new',err=200)
      
*    ... PDF at given Q2 vs. x
      DO jq2=1,nq2pdf
         Q2 = Q2pdf(jq2)
         WRITE(nu,5103) Q2 
 5103    FORMAT(//,' Q2=',F10.2)
         WRITE(nu,5104)
 5104    FORMAT(/,3X,' X',10X,'xd',18X,'xu',18X,'xg',18X,'xub',
     2        17X,'xdb',17X,'xs',18X,'xc',18X,'xb')
         do j=1,nx
            x = xgrid(j)
            call pdfsa(x,Q2,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
            pdf0(-5,j)=bb0
            pdf0(-4,j)=cb0
            pdf0(-3,j)=sb0
            pdf0(-2,j)=db0
            pdf0(-1,j)=ub0
            pdf0(-0,j)=glue0
            pdf0(1,j)=u0
            pdf0(2,j)=d0
            WRITE(nu,5031) xgrid(j)
     &           ,(pdf0(ipt,j),0d0,ipt=2,-5,-1)
 5031       FORMAT(E10.3,16E10.3) 
         end do
      end do
c      
*    ... PDF at given x vs. Q2
      WRITE(nu,*)
      WRITE(nu,*) '-----------------------------------------------'
      DO jx=1,nxpdf
         x = xpdf(jx)
         WRITE(nu,5203) x 
 5203    FORMAT(//,'  x=',E15.3)
         WRITE(nu,5204)
 5204    FORMAT(/,4X,'Q2',9X,'xd',18X,'xu',18X,'xg',18X,'xub',
     2        17X,'xdb',17X,'xs',18X,'xc',18X,'xb')
         do j=1,nq2
            Q2 = q2grid(j)
            call pdfsa(x,Q2,u0,d0,ub0,db0,sb0,cb0,bb0,glue0)
            pdf0(-5,j)=bb0
            pdf0(-4,j)=cb0
            pdf0(-3,j)=sb0
            pdf0(-2,j)=db0
            pdf0(-1,j)=ub0
            pdf0(-0,j)=glue0
            pdf0(1,j)=u0
            pdf0(2,j)=d0
            WRITE(nu,5031) q2grid(j)
     &           ,(pdf0(ipt,j),0d0,ipt=2,-5,-1)
         end do
      end do                    ! on to the next xB value...
      
      goto 210
 200  print*, 'WARNING(writepdf): pdffile already exists: "'
     &     //pdffile(1:len)//'"'
 210  close(nu)
      
      return
      end






      subroutine writetbl(tblfile)
C     Tabulates and writes to file the CJ pdfs in the "tbl" format 
C     used by the CTEQ6 series, so that they can be used in CTEQ functions. 
C     It follows the correct CTEQ6 formatting and everything
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
C
C     Programmers: Lucas Brady, Alberto Accardi
C     Date: Jul 2011
C

      implicit none

      character*50 tblfile

      integer iQ,ix,i,j,len,outn,NextUn,nu2
      real*8 x,Q,Q2,u0,d0,ub0,db0,cb0,sb0,bb0,glue
      real*8 ord,nfl,lambda,qm1,qm2,qm3,qm4,qm5,qm6
      integer nfmx
      real*8 qini,qmax,xmin
      
*    *** CJ parameters
      double precision par(100)
      common/param/par
      integer ILOOP,IORD,NMAX,IVL,ncb
      double precision xmc,xmb
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      double precision FLAVOR
      COMMON/GRPTHY/FLAVOR

*    *** x and Q grids
      integer nQ,nx
      parameter(nx=96) 
      parameter(nQ=25)
      double precision fullpdf(8*nx*nQ)      
      double precision Qgr(nQ), xgr(nx)
      data Qgr/   1.30000D+00, 1.50159D+00, 1.75516D+00, 2.07810D+00,
     &            2.49495D+00, 3.04086D+00, 3.76715D+00, 4.75000D+00,
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


      ord=iord+1
      nfl=flavor
      qm1=0d0   
      qm2=0d0
      qm3=0.2d0    
      qm4=xmc
      qm5=xmb
      qm6=180d0
      nfmx=5
      qini=Qgr(1)
      qmax=Qgr(nQ)
      xmin=xgr(1)

      lambda=par(1)
      
      do iQ=1,nQ
          Q=Qgr(iQ)
          do ix=1,nx
              x=xgr(ix)
              ! PDFs
              q2=q*q
              call PDFsa(x,Q2,U0,D0,UB0,DB0,SB0,CB0,BB0,GLUE)
              i=(0*nQ+(iQ-1))*nx+ix
c              if(bb0.eq.0.d0.and.x.ne.0.d0.and.x.ne.1.d0)then
c                 print*,x,q2,bb0/x
c              endif
              fullpdf(i)=bb0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(1*nQ+(iQ-1))*nx+ix
              fullpdf(i)=cb0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(2*nQ+(iQ-1))*nx+ix
              fullpdf(i)=sb0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(3*nQ+(iQ-1))*nx+ix
              fullpdf(i)=db0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(4*nQ+(iQ-1))*nx+ix
              fullpdf(i)=ub0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(5*nQ+(iQ-1))*nx+ix
              fullpdf(i)=glue/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(6*nQ+(iQ-1))*nx+ix
              fullpdf(i)=u0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
              i=(7*nQ+(iQ-1))*nx+ix
              fullpdf(i)=d0/x
              if ((ix.eq.1).or.(ix.eq.nx)) fullpdf(i)=0d0
          enddo
          
      enddo

      !writes to file
       
      outn=NextUn()
      open (outn,FILE=tblfile,status='new',err=200)
      
      write(outn,*) " Parton Distribution Table : CJ"
      write(outn,*) "  Ordr, Nfl, lambda        Qmass"//
     &              " 1,  2,  3,         4,  5,  6"
      write(outn,22) ord,nfl,lambda,qm1,qm2,qm3,qm4,qm5,qm6
      write(outn,*) "   NX,  NQ,  NfMx"
      write(outn,23) nx,nq,nfmx
      write(outn,*) "QINI, QMAX, (QV(I), I =1, NQ)"
      write(outn,24) qini,qmax
      
      write(outn,24) qgr

      write(outn,*)"XMIN, (XV(I), I =1, NX)"
      write(outn,24)xmin

      write(outn,24) xgr
      write(outn,*)"Parton Distribution Table:"

      write(outn,25) fullpdf      
      goto 210
 200  call trmstr(tblfile,len)
      print*, 'WARNING(writetbl): tblfile already exists: "'
     &     //tblfile(1:len)//'"'
 210  close (outn)
      
      
 22   format(F6.0,F6.0,F8.4,F9.3,F9.3,F9.3,F9.3,F9.3,F9.3)
 23   format(I6,I5,I5)
 24   format(6E14.4)
 25   format(5E14.6)
      return
      end




      subroutine QCDev
*     Initializes the QCD evolution by setting appropriate constants 
*     and computing needed variables. Then pereforms QCD evolution itself.
*     It also calls WATE16 and WATE32 which are needed in the
*     subroutines called by 'intqcd' (Note: no harm is done if the WATExx 
*     series of routine were already, or will be called again somewhere else.)

      implicit real*8 (a-h,o-z)
*    ... thresholds and other parameters for evolution
      common/threshold/sb
      common/stepsize/delta
      common/scales/s0,alambda5
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,xmc,xmb,ncb
      COMMON/Q2STUF/ Q02,Q2MAX
*    ... 'par' must be initialized using 'readpar'
      common/param/par(100)
*    ... matrices of evolved PDFs
      COMMON/GFUNC/CALC(11,60,60) 
      COMMON/GINT/GF(11,60,2160)

*    *** Gaussian integrations needed, e.g., by 'intqcd' and 'pdfnrm'
      CALL WATE16
      CALL WATE32
      CALL WATE96

*    *** initializes GF array
      DO 10 I=1,11 
      DO 10 J=1,60
      DO 10 K=1,60
         calc(I,J,K)=0.
 10   continue

*    *** DGLAP evolution of PDFs
*    ... sets parameters for evolution 


      alambda5=par(1)
      S0=DLOG(Q02/alambda5**2) 
      Q2max = 1d10               ! arbitrary value; covers up to LHC
      SMAX=DLOG(DLOG(q2max/alambda5**2)/S0) 
      SB=DLOG(DLOG(XMB**2/alambda5**2)/S0)
      DELTA=SB/NCB
      NMAX=SMAX/DELTA+3 
      IF(NMAX.GT.60) NMAX=60
      CALL INTQCD(PAR)
      call fillgf 

      return 
      end

   
      subroutine fillgf
*     Fills the first part of the GF matrix with CALC
*     (CALC is used by the QCD evolution, GF by the 'pdf' routine 
*     in the fitting package)

      implicit real*8 (a-h,o-z)
*    ... matrices of evolved PDFs
      COMMON/GFUNC/CALC(11,60,60) 
      COMMON/GINT/GF(11,60,2160)

*    ... fills the gf matrix
      DO I=1,11
         DO J=1,60
            DO K=1,60
               GF(I,J,K)=CALC(I,J,K)
            end do
            DO K=60,2160
               GF(I,J,K)=0d0
            end do
         end do
      end do
 
      return
      end


      SUBROUTINE PDFsa(X,Q2,U0,D0,UB0,DB0,SB0,CB0,BB0,GLUE)
C     End-user routine for PDF interpolation. Assumes that the calc and gf
C     matrices are obtained by QCD evolution of PDF's at Q0 with 'QCDev',
C     or by reading them from an external file with 'readgf'.
C     The starting scale s0, and Lambda_QCD with 5 flavors need to be 
C     precomputed in the common /scales/
      implicit none
      integer iflagb
      double precision x,Q2,u0,d0,ub0,db0,sb0,cb0,bb0,glue
      double precision s
      double precision s0,alambda5
      common/scales/s0,alambda5
      common/test/iflagb
      s=dlog(dlog(q2/alambda5**2)/s0)
      if(s.lt.0d0) then ! this freezes PDF at Q0 when called at smaller Q
         s=0d0
      endif
      if(q2.eq.1000d0)then
         iflagb=1
      else
         iflagb=0
      endif
      call FSUPDF(0,X,S,U0,D0,UB0,DB0,SB0,CB0,BB0,GLUE)
      return
      end
