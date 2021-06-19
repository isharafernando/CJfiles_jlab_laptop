****************************************************************
*     Calculation of the CJ PDF error sets
*
*     produce 2*NFREE .par files and .gf files corresponidng
*     to the "eigenPDFs" needed for PDF error calculations
*     (Diagonalizes the error matrix, and shifts the parameters according
*     to the eigenvectors normalized to the eigenvalue).
*
*     Programmers: aa = A.Accardi (accardi@jlab.org)
*
*     HISTORY:
*
*     PDFerr      _the first one_
*     (28 Jul 11)
*
****************************************************************
      program PDFerrors

      implicit none

      character*50 file,parfile,eigenfile,errfile,gffile,tblfile,pdffile
     &     ,datfile(3),datdir,datdirp,datdirn,datdird
     &     ,infofile,infofilep,infofilen,infofiled
     &     ,type,tbldir
     &     ,pardir,pdfdir
      character*1000 desc,auth,ref,flav
      integer part,partp,partn,partd
     &     ,index,indexp,indexn,indexd
     &     ,ver,n,nu,NextUn
      character*2 cpdf 
      character*4 dpdf
      character*100 uflag,sfn
 
      integer len,len1,ipar,i,j,npdf,k
      integer itest,nfg,ierrPDF,istruc,ibeam
      integer*4 today(3),now(3)

      ! Command line flags and user choices 
      integer jpdf,jpar,jeigen,jtbl,jgf,jlhapdf,jtblsfn,jlhasfn,nfree
      data jpdf,jpar,jeigen,jtbl,jgf,jlhapdf,jtblsfn,jlhasfn,nfree/ 
     &     0   ,0   ,0     ,0   ,0  ,0      ,0      ,0      ,0    /
      double precision deltachi
      data deltachi/1.645/  ! nominal Gaussian 90% c.l. 
      logical errorPDFflag
      data errorPDFflag/.true./
      
      integer np
      parameter (np=50)
      double precision Amat(np,np),Bvec(np),epsilon(np,np),eval(np)
     &     ,E(np),norm,norm2,par0(100),indx(np),d,delta,tmp,fac 

      double precision del(np,np),det,myeps(np),nrm(np)
       
*     alpha matrix and paramter errors
      integer pmax
      parameter(pmax=50) ! must be = dimension of alpha in 'minim'
      double precision alfa(pmax,pmax),central(pmax),eps(pmax)
      integer pos(pmax),npdfpar
      common/PDFerr/alfa,central,eps,pos,npdfpar
      double precision cov(pmax,pmax)
      common/PDFcov/cov
c$$$*     alpha matrix and parameter errors
c$$$      integer pmax
c$$$      parameter(pmax=50) ! must be = dimension of alpha in 'minim'
c$$$      double precision alfa(pmax,pmax),alfainv(pmax,pmax)
c$$$     &     ,central(pmax),eps(pmax)
c$$$      integer pos(pmax),npdfpar
c$$$      double precision cov(pmax,pmax)
c$$$      common/PDFerr/alfa,cov,central,eps,pos,npdfpar

*     parameters for QCD evolution
      double precision par(100)
      common/curpar/par

*     Other PDF parameters 
      character*10 pname(100)  ! needed if writing a par file with 'writepar'
      double precision uncrt(100)
      common/oparam/pname,uncrt

*     Total number of parameters
      integer npar              
      common/npar/npar

*     Theory correction flag
      integer inuke,itmc,iht
      common/thcor/inuke,itmc,iht
      
*     comment lines in .par file
      character*100 line(3)
      common/parfil/line


      if (iargc().eq.0) then    ! HELP
         print*, 'USAGE: PDFerr <file> [options]'
         print*, 'where file = .par file from fitpack (no suffix)'
         print*, 'Options:'
         print*, '-nfree nn         : no. of eigendirections to invert' 
         print*, '-deltachi r       : deltachi for PDF errors '
     &        //'(default r=1.645)'
         print*, '-noerrorsets      : only writes central PDFs'
         print*, '-igf              : writes native CJ format grids'
         print*, '-itbl             : writes CTEQ6 format grids'
         print*, '-itbl_sfn ibm ist : writes str.fns. instead '
     &        //'(ibm=1 ist=2 for NC F2)'
         print*, '-ilhapdf          : writes LHAPDF6 grids '
     $        //'(headers and .info are hard coded)'
         print*, '-ilhapdf_sfn ibm  : writes str.fns. instead '
     &        //'(ibm=1 for NC)'
         print*, '-ipar             : writes .par parameter file'
         print*, '-ipdf             : writes sample .pdf file'
         print*, '-ieigen           : writes eigenvectors in .eigen '
     &        //'file'        
         stop
      end if

      call getarg(1,file)        ! gets file name (no suffix, please)
      call trmstr(file,len)
      parfile = file(1:len)//'.par'

      nfg=1                     ! gets optional flags
      do
         nfg=nfg+1
         call getarg(nfg,uflag)
         if (uflag.eq.'') then
            EXIT
         else if (uflag.eq.'-nfree') then
            nfg=nfg+1
            call getarg(nfg,uflag)
            read (uflag,*) nfree
         else if (uflag.eq.'-deltachi') then
            nfg=nfg+1
            call getarg(nfg,uflag)
            read (uflag,*) deltachi
         else if (uflag.eq.'-noerrorsets') then
            errorPDFflag = .false.
         else if (uflag.eq.'-igf') then
            jgf = 1
         else if (uflag.eq.'-itbl') then
            jtbl = 1
         else if (uflag.eq.'-ilhapdf') then
            jlhapdf = 1
         else if (uflag.eq.'-ipar') then
            jpar = 1
         else if (uflag.eq.'-ipdf') then
            jpdf = 1
         else if (uflag.eq.'-ieigen') then
            jeigen = 1
         else if (uflag.eq.'-itbl_sfn') then
            jtbl=1
            jtblsfn=1
            nfg=nfg+1
            call getarg(nfg,uflag)
            read (uflag,*) ibeam
            nfg=nfg+1
            call getarg(nfg,uflag)
            read (uflag,*) istruc
            if (istruc.eq.0) sfn='Fl'
            if (istruc.eq.1) sfn='F1'
            if (istruc.eq.2) sfn='F2'
            if (istruc.eq.3) sfn='F3'
            if (ibeam.eq.1) sfn=trim(sfn)//'NC'
            if (ibeam.eq.2) sfn=trim(sfn)//'nu'
            if (ibeam.eq.3) sfn=trim(sfn)//'nubar'
         else if (uflag.eq.'-ilhapdf_sfn') then
            jlhapdf=1
            jlhasfn=1
            nfg=nfg+1
            call getarg(nfg,uflag)
            read (uflag,*) ibeam
         end if
      end do

      ! reads the central value pars, and inverts the Hessian
      write(*,*) 'Working on '//file
      call readpar(parfile,par,inuke,itmc,iht)      
      ! default choice of minor to invert
      if (nfree.eq.0) then
         !nfree = 18          ! CJ15 PDF        params
         !nfree = 21          ! CJ15 PDF+off    params
         !nfree = 24          ! CJ15 PDF+off+HT params
         nfree = npdfpar      ! all parameters (including normalizations) 
      end if
      
c     NOTE: alfa is the Hessian matrix from the fit (from /PDFerr/)
c     tredn2 reduces alpha to a tridiagonal matrix
c     tqli finds the eigenvalues
c
      call tredn2(alfa,nfree,pmax,eval,E) ! NOTE: alfa is destroyed 
      call tqli(eval,E,nfree,pmax,alfa)
      ! on exit, eval contains eigenvalues, and the Kth column of epsilon 
      ! the Kth eigenvector, i.e., alfa(n,K) = Kth-eigenvec(n)
      
c jfo Eigenvectors are normalized to unity
c jfo Change normalization to be equivalent to choosing the scale factors 
c jfo s_k=1/dsqrt(eval(k)) as introduced in the Tung Hessian paper. 
c jfo 
      do k=1,nfree                ! cycles thorugh the eigenvectors
         norm = 1./dsqrt(eval(k))
         do i=1,nfree
            epsilon(i,k) = norm*alfa(i,k)
         end do
      end do                    ! Now epsilon(*,k) contains the normalized 
                                ! k-th eigenvector
c jfo      print*,(eval(i),i=1,nfree)
c jfo      do k=1,nfree
c jfo         print*,k,(epsilon(i,k),i=1,nfree)
c jfo      enddo


      ! Output .par file
      if (jpar.eq.1) then
         pardir = 'par_'//trim(file)//'/'
         call system('mkdir '//pardir)
         parfile = trim(pardir)//file(1:len)//'_00.par'
         call trmstr(line(1),len1)
         tmp = alfa(1,1)
         alfa(1,1) = 0d0        ! to suppress output of alfa
         call writepar(parfile,line,pname,par,uncrt,alfa,nfree
     &        ,inuke,itmc,iht)
         alfa(1,1) = tmp
      end if

      ! Write eigenvalues and eigenvectors to file
      if (jeigen.eq.1) then 
         eigenfile = file(1:len)//'.eigen'
         nu = NextUn()
         open(unit=nu,file=eigenfile,status='new',err=1400)
         write(nu,*) 'Eigenvalues:'
         write(nu,1390) (eval(n),n=1,nfree)
         write(nu,*) 'Eigenvectors:'
         do n = 1,nfree
            write(nu,1390) (alfa(n,ipar),ipar=1,nfree)
         end do
         write(nu,*) 'Normalization:'
         do ipar=1,nfree
            nrm(ipar)=0d0
            do n=1,nfree
               nrm(ipar)=nrm(ipar)+alfa(n,ipar)**2
            end do
            nrm(ipar)=dsqrt(nrm(ipar))
         end do
         write(nu,1390) (nrm(n),n=1,nfree)      
 1390    format(1X,100E18.10)
         goto 1500
 1400    print*, 'WARNING(PDFerr): eigenfile already exists: "'
     &        //trim(eigenfile)//'"'
 1500    close(nu)
      end if

      call QCDev
      
      ! Native CJ .gf grids
      if (jgf.eq.1) then 
         gffile  = file(1:len)//'_00.gf'
         call writegf(gffile)
      end if

      ! CTEQ6 format .tbl grids
      if (jtbl.eq.1) then
         if (jtblsfn.eq.0) then
            ! PDF grids
            tbldir = 'tbl_'//trim(file)//'/'
            call system('mkdir '//tbldir)
            tblfile  = trim(tbldir)//file(1:len)
            call writetbl(trim(tblfile)//'_00.tbl')
         else
            ! str.fns. grids
            tbldir = 'tbl'//'_'//trim(file)//'-'//trim(sfn)//'/'
            call system('mkdir '//tbldir)
            tblfile  = trim(tbldir)//trim(file)//'-'//trim(sfn)
            call writetbl_sfn(trim(tblfile)//'_00.tbl',ibeam,istruc)
         end if
      end if
      
      ! LHAPDF grids and info
      if (jlhapdf.eq.1) then 

         if (jlhasfn.eq.0) then
            ! PDF grids
!           # Official CJ indexes:
!           # 12000-12200 : CJ12 min mid max
!           # 12300-12400 : CJ15 lo nlo
!           # 12500-12549 : CJ15 proton NC str.fns
!           # 12550-12599 : CJ15 neutron NC str.fns
!           # 12600-12649 : CJ15 deuteron NC str.fns
c            desc = "CJ12 series global PDFs with nuclear and "
c     &        //"finite-Q^2 corrections. "
c     &        //"These grids correspond to a tolerance factor T=1; users "
c     &        //"are welcome to scale these with a T factor appropriate "
c     &        //"for their intended applications, see Eqs.(6)-(7) in "
c     &        //"Owens et al., PRD87(2013)094012 [arXiv:1212.1702]." 
c            auth = "Owens J.F, Accardi A, and Melnitchouk W"
c            ref = "Phys.Rev. D87 (2013) 9, 094012 [arXiv:1212.1702]"
c     &        //" - e-mail accardi@jlab.org"
c            index = 12000
c            ver = 2
            desc = "CJ15 series global PDFs with nuclear and "
     &           //"finite-Q^2 corrections. The 48 error grids"
     &           //"are calculated with Delta\chi=1.645; if "
     &           //"experimental errors were Gaussian and all "
     &           //"experiments compatible, this would correspond"
     &           //"to a 90% confidence level. Errors on observables "
     &           //"can be calculated using Eqs.(6)-(7) with T=1 from "
     &           //"Owens et al., PRD87(2013)094012 [arXiv:1212.1702]."
            auth= "Accardi A, Brady L T, Melnitchouk W, Owens J F"
     &           //", Sato N"
            ref = "arXiv:1602.03154v2 - e-mail accardi@jlab.org"
            index = 12400
            ver = 2
            part = 2212         ! proton PDFs
            flav = '-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21'

            datdir = trim(file)//'/'
            call system('mkdir '//datdir)
            datfile(1)  = trim(datdir)//file(1:len)//'_0000.dat'
            type = 'central'
            call writedat(datfile(1),type,flav)
            infofile = trim(datdir)//file(1:len)//'.info'
            call writeinfo(infofile,desc,part,flav,auth,ref,index
     &           ,ver,nfree)

         else 
            ! Structure functions
            desc = "CJ15 DIS structure functions. The 48 error grids"
     &           //"are calculated with Delta\chi=1.645; if "
     &           //"experimental errors were Gaussian and all "
     &           //"experiments compatible, this would correspond"
     &           //"to a 90% confidence level. Errors on observables "
     &           //"can be calculated using Eqs.(6)-(7) with T=1 from "
     &           //"Owens et al., PRD87(2013)094012 [arXiv:1212.1702]."
     &           //"Flavor indexes (NC structure functions): "
     &           //"901-903=F2(full),F2(LT),F2(massles);"
     &           //"904-906=same for FL; "
     &           //"907-909=same for F3"
            auth= "Accardi A, Brady L T, Melnitchouk W, Owens J F"
     &           //", Sato N"
            ref = "arXiv:1602.03154v2 - e-mail accardi@jlab.org"
            indexp = 12550
            indexn = 12600
            indexd = 12650
            ver = 1
            partp = 2212         ! proton SFs
            partn = 2112         ! neutron SFs
            partd = 1000010020  ! deuteron SFs
            flav = '901, 902, 903, 904, 905, 906, 907, 908, 909'

            datdirp = trim(file)//'-FpNC/' 
            datdirn = trim(file)//'-FnNC/' 
            datdird = trim(file)//'-FdNC/' 
            call system('mkdir '//datdirp)
            call system('mkdir '//datdirn)
            call system('mkdir '//datdird)
            datfile(1)  = trim(datdirp)//file(1:len)
     &           //'-FpNC_0000.dat'
            datfile(2)  = trim(datdirn)//file(1:len)
     &           //'-FnNC_0000.dat'
            datfile(3)  = trim(datdird)//file(1:len)
     &           //'-FdNC_0000.dat'
            type = 'central'
            call writedat_sfn(datfile,type,flav,ibeam)   
            infofilep = trim(datdirp)//file(1:len)//'-FpNC.info'            
            infofilen = trim(datdirn)//file(1:len)//'-FnNC.info'            
            infofiled = trim(datdird)//file(1:len)//'-FdNC.info'            
            call writeinfo(infofilep,desc,partp,flav,auth,ref,indexp
     &           ,ver,nfree)            
            call writeinfo(infofilen,desc,partn,flav,auth,ref,indexn
     &           ,ver,nfree)            
            call writeinfo(infofiled,desc,partd,flav,auth,ref,indexd
     &           ,ver,nfree)            

         end if

      end if

      ! Sample PDF file
      if (jpdf.eq.1) then
         pdfdir = 'pdf_'//trim(file)//'/'
         call system('mkdir '//pdfdir)
         pdffile = trim(pdfdir)//file(1:len)//'_00.pdf'
         call writepdf(pdffile)
      end if

      
*    *** Calculates the PDF error sets

      itest=1
      if(itest.eq.2) stop ""

      ierrPDF = (1-jpdf)*(1-jpar)*(1-jtbl)*(1-jgf)*(1-jlhapdf)      

      if ((errorPDFflag).and.(ierrPDF.eq.0)) then

         npdf = 0
         do i = 1, npar
            par0(i) = par(i)
         end do

         do ipar = 1, nfree
            print*
            print*, 'Working on parameter #',ipar,' of',nfree
            !print*, 'eigenvalue =',eval(ipar)
            !print*, 'eigenvector ='
            !print*, (epsilon(n,ipar),n=1,nfree)

            ! creates and writes new par file
            do i = 1, npar
               par(i) = par0(i)
            end do
            do i=1,-1,-2        ! positive and negative shifts
               npdf = npdf + 1

               do j=1,nfree     ! shifted parameter set with tolerance
c
c jfo For a chi square tolerance T, one should shift the parameters 
c jfo by an amount sqrt(T)/2 times the error
c jfo T=50 is used here   AA: why divided by 2?
c     jfo
c               fac=3.5355  ! for T=50
c               fac=10.0    ! for T=100
c               fac=1.0     ! for Dchi=1 (T=1) == 68% Gauss.equiv. conf level
c               fac=1.645d0   ! for 90% Gauss.equiv. conf level
c
c  fac=1 corresponds to delta chi = 1 provided that one divides the
c  error on an observable by 2
c  delta sigma = 1/2 sqrt(sum(sigma(a+)-sigma(a-)**2)
c  where the sum runs over the varied parameters denoted by 'a'
c  and a+ (a-) are the parameters in the odd(even) numbered PDF error sets
c
                  fac = deltachi ! deltachi is user-decided, see flags
                  par(pos(j)) = par0(pos(j)) + i*fac*epsilon(j,ipar)
               end do

               call QCDev       ! performs QCD evolution 
                                ! (and normalizes valence and glue)
            
               ! CJ parameters file
               if (jpar.eq.1) then    
                  write(cpdf,'(I2.2)') npdf
                  errfile = trim(pardir)//file(1:len)//'_'//cpdf//'.par'
                  line(1) = line(1)(1:len1)//' - pdf set #'//cpdf
                  call idate(today) ! today(1)=day, (2)=month, (3)=year
                  call itime(now) ! now(1)=hour, (2)=minute, (3)=second
                  write(line(3),1000) today(2), today(1), today(3), now
 1000             format ( '# Created on ', i2.2, '/', i2.2, '/'
     &                 , i4.4, ' at ',i2.2, ':', i2.2, ':', i2.2 )
                  tmp = alfa(1,1)
                  alfa(1,1) = 0d0 ! to suppress output of alfa
                  call writepar(errfile,line,pname,par,uncrt,alfa,nfree
     &                 ,inuke,itmc,iht)
                  alfa(1,1) = tmp
               end if

               ! CJ native PDF grids
               if (jgf.eq.1) then
                  write(cpdf,'(I2.2)') npdf
                  gffile  = file(1:len)//'_'//cpdf//'.gf'
                  call writegf(gffile)
               end if
               
               ! CTEQ6-series grids
               if (jtbl.eq.1) then
                  write(cpdf,'(I2.2)') npdf
                  if (jtblsfn.eq.0) then
                     ! PDF grids
                     call writetbl(trim(tblfile)//'_'//cpdf//'.tbl')
                  else
                     ! str.fns. grids
                     call writetbl_sfn(trim(tblfile)//'_'//cpdf//'.tbl'
     &                    ,istruc,ibeam)
                  end if
               end if
               
               ! LHAPDF6 format grids
               if (jlhapdf.eq.1) then
                  write(dpdf,'(I4.4)') npdf
                  if (jlhasfn.eq.0) then
                     ! PDF grids
                     datfile(1) = trim(datdir)//file(1:len)//'_'
     &                    //dpdf//'.dat'
                     type = 'error'
                     call writedat(datfile(1),type,flav)
                  else
                     ! str. fns. grid
                     datfile(1)  = trim(datdirp)//file(1:len)
     &                    //'-FpNC_'//dpdf//'.dat'
                     datfile(2)  = trim(datdirn)//file(1:len)
     &                    //'-FnNC_'//dpdf//'.dat'
                     datfile(3)  = trim(datdird)//file(1:len)
     &                    //'-FdNC_'//dpdf//'.dat'
                     call writedat_sfn(datfile,type,flav,ibeam)
                  end if
               end if
                  
               ! sample .pdf files for quality control
               if (jpdf.eq.1) then
                  write(cpdf,'(I2.2)') npdf
                  pdffile = trim(pdfdir)//file(1:len)//'_'//cpdf//'.pdf'
                  call writepdf(pdffile)
               end if

            end do

         end do

      else

         print*
         print*, 'WARNING: no error PDF file format selected'
         print*, '         --> no output written to file'
         print*
         
      end if

      stop "END of PDFerr"
      end
