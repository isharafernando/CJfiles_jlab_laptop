*I.   Truncated moments with CTEQX structure functions
      program truncmom
*     programmer: Alberto Accardi 
*     v. 1:   18 Mar 2010
* 
*  A. COMMENTARY
*
*     Computes truncated moments of F2 vs Q^2
*
*     USAGE: !under construction!
*     
*      See "parameters file" for options and parameters
* 
*     COMPILATION (see makefile):  
*
*       make truncmom
*
*     INPUT FILES:
*
*      -  truncmom.dat     -- parameters file
*      -  +++.par          -- cteqX PDF parameter file 
*
*     HYSTORY:
*     
*     * truncomom  : the _first_ one                             
*       (18 Mar 10)  
*                   
*     TO DO LIST: 
*     -----------
*
*     KNOWN BUGS: 
*     -----------
*
*  B. DECLARATIONS
      implicit none

      logical lflag

      character outfile*150,tmpfile*150,string*100

      integer ndat,outlen,out,i,ireg,nreg,iQ2,len,nout(2),comlen
     $     ,foldlen 

      double precision v(4),par(100),sfn,xB,xBmax,xBmin,delQ2,sfn_p
     $     ,sfn_D 

*    *** functions      

      integer NextUn

      double precision xW

*    *** test zone declarations 

*    *** common blocks

*     File name routine
      character foldname*80,comment*80,parfile*80,linlog*3
      double precision Q2min,Q2max
      integer nQ2,istruc,ibeam,itm,iht,inuke
      common /filnam/ Q2min,Q2max
     &     ,nQ2,istruc,ibeam,itm,iht,inuke
     &     ,foldname,comment,parfile,linlog
      
      
*     Gaussian Integration
      double precision XI,WI,XX
      integer NTERMS
      COMMON/GAUS16/XI(16),WI(16),NTERMS,XX(17)
c      COMMON/GAUS32/XI(32),WI(32),NTERMS,XX(33)

*    *** Parameters and initial values

      integer nQ2max,nregmax
      parameter(nQ2max=100,nregmax=50)
      double precision Q2(0:nQ2max),xBlim(0:nQ2max,nregmax)
     $     ,mom_p(0:nQ2max,nregmax),mom_D(0:nQ2max,nregmax)
     $     ,mom_Dp(0:nQ2max,nregmax)

      integer nWbar
      parameter(Nwbar=20)
      double precision W2bar(Nwbar)
*      data W2bar/1.3,1.9,2.5,3.1,3.9,4.5,6.0,9.0,15.,11*0d0/
      data W2bar/1.3,1.9,2.5,3.1,3.9,4.5,6.0,9.0,15.,30.,10*0d0/


*    ***

* 
*  D. ACTION
*

*    ********************
*    *  Initialization  *
*    ********************

*    *** output file
      
      ndat = NextUn()
      open(ndat,file='truncmom.dat',status='old',form='formatted')

*    ... output file name
      read(ndat,*) string, foldname
      call trmstr(foldname,foldlen)
      print*, 'folder:  '''//foldname(1:foldlen)//''''
      read(ndat,*) string, comment
      call trmstr(comment,comlen)
      print*, 'comment:  '''//comment(1:max(1,comlen))//''''

*    *** main parameters

      read(ndat,*) string
*    ... Q2 range
      read(ndat,*) string
      read(ndat,*) Q2min, Q2max, nQ2,linlog
      if (linlog.eq.'log') then
         Q2max = dlog(Q2max)/dlog(10d0)
         Q2min = dlog(Q2min)/dlog(10d0)
      end if
      delQ2 = (Q2max-Q2min)/nQ2

      print*, ' Q2min, Q2max, nQ2 =', Q2min, Q2max, nQ2
*    ... PDF parametr file
      read(ndat,*),string,parfile
*    ... Str.Fn. switches
      call trmstr(parfile,len)
      print*, 'parfile =',parfile(1:len)
      read(ndat,*),string,ibeam
      read(ndat,*),string,istruc
      print*, 'istruc, ibeam =',istruc,ibeam
      print*, '*1'
      call filename(outfile,outlen)
      print*, '*2'
      tmpfile = outfile(1:outlen)
      inquire(file=tmpfile,exist=lflag)
      if (lflag) then
         print*, 'ERROR: outfile exists: '''//outfile(1:outlen)//''''
         stop
      else
         write(*,'(/,A,/)') 'Output in: '''//outfile(1:outlen)//''''
      end if
      close (ndat)

*    *** Reads in parameter file, and store params in 'par(100)'
      call readpar(parfile,par,inuke,itm,iht)
      print*, '*3'
*    ... evolves the PDFs
      print*, 'Performing QCD evolution...'
      call QCDev
      print*, '*4'
      print*, '... done!'
      print*

*    ***************
*    *  MAIN LOOP  *
*    ***************

*    *** counts W regions
      do i = 1, nWbar
         if (W2bar(i).gt.0d0) nreg=i-1
      end do

*    *** computes truncated moments
      do iQ2 = 0, nQ2
         if (linlog.eq.'lin') then
            Q2(iQ2) = Q2min + iQ2*delQ2
         else
            Q2(iQ2) = 10d0**(Q2min + iQ2*delQ2)
         end if
         do ireg = 1, nreg
            xBmin = xW(Q2(iQ2),W2bar(ireg+1))
            xBmax = xW(Q2(iQ2),W2bar(ireg))
            mom_p(iQ2,ireg) = 0d0
            mom_D(iQ2,ireg) = 0d0
            DO I=1,NTERMS
               xB=0.5*(xBmax-xBmin)*XI(I)+0.5*(xBmax+xBmin)
               v(1)=xB
               v(2)=Q2(iQ2)
               if (xBmax.lt.0.985) then
                  call targ(1d0,0d0) ! proton target
                  call strfn(ibeam,istruc,itm,iht,0,v,par,sfn_p)
                  call deuteronsf(ibeam,istruc,inuke,itm,iht,0,v,par
     $                 ,sfn_D)
               else
                  sfn_p = 0d0
                  sfn_D = 0d0
               end if
               mom_p(iQ2,ireg) = mom_p(iQ2,ireg)
     $              + .5*(xBmax-xBmin)*WI(I)*sfn_p
               mom_D(iQ2,ireg) = mom_D(iQ2,ireg)
     $              + .5*(xBmax-xBmin)*WI(I)*sfn_D
            end do
            xBlim(iQ2,ireg)   = xBmax
            xBlim(iQ2,ireg+1) = xBmin
            if (mom_p(iQ2,ireg).ne.0d0) then
               mom_Dp(iQ2,ireg) = mom_D(iQ2,ireg)/mom_p(iQ2,ireg)
            else
               mom_Dp(iQ2,ireg) = 0d0
            end if
         end do
      end do

*    *** Writes output
      
      nout(1) = 6               ! screen
      nout(2) = NextUn()        ! file
      open(nout(2),file=outfile(1:outlen),status='new',form='formatted')

      do i = 1, 2
         write(nout(i),*) 'Truncated moments calculation'
         write(nout(i),*) '-----------------------------'
         write(nout(i),*) 
         write(nout(i),*) ' PROTON'
         write(nout(i),*) '          W ---->'
         write(nout(i),100) '    Q2',(W2bar(ireg),ireg=1,nreg+1)
         do iQ2 = 0, nQ2
            write(nout(i),110) Q2(iQ2),(mom_p(iQ2,ireg),ireg=1,nreg)
         end do
         write(nout(i),*) 
         write(nout(i),*) ' D/p ratio'
         write(nout(i),*) '          W ---->'
         write(nout(i),100) '    Q2',(W2bar(ireg),ireg=1,nreg+1)
         do iQ2 = 0, nQ2
            write(nout(i),110) Q2(iQ2),(mom_Dp(iQ2,ireg),ireg=1,nreg)
         end do
         write(nout(i),*) 
         write(nout(i),*) ' xB limits'
         write(nout(i),*) '          W ---->'
         write(nout(i),100) '    Q2',(W2bar(ireg),ireg=1,nreg+1)
         do iQ2 = 0, nQ2
            write(nout(i),110) Q2(iQ2),(xBlim(iQ2,ireg),ireg=1,nreg+1)
         end do
      end do
 100  format(A6,50(F8.1,3X))
 110  format(F6.1,50E11.3)

      write(*,'(/,A,/)') 'Output in: '''//outfile(1:outlen)//''''

      call exit
      end


**********************************************

      function xW(Q2,W2)
*     calculates x as a fn. of Q^2 at fixed W
*     IN:  Q2  [GeV^2] = (dp) photon virtuality squared
*          W2  [Gev^2] = (dp) invariant mass
*     OUT: xW          = (dp) Bjorken x at given W

      implicit none
      double precision xW,Q2,W2,M2
      parameter (M2=0.939**2)

      xW = Q2/(Q2+W2-M2)

      return 
      end


**********************************************

      subroutine filename(outfile,outlen)

      implicit none

      character outfile*150,cdum*80

      integer i,max,ismear,ioff,ibj
      parameter (max=9)
      character car*(max),car1*(max)

      integer ll,len,oldlen,outlen,min,min1
      double precision mmin,mmax

*    *** common blocks

*     File name routine
      character foldname*80,comment*80,parfile*80,linlog*3
      double precision Q2min,Q2max
      integer nQ2,istruc,ibeam,itm,iht,inuke
      common /filnam/ Q2min,Q2max
     &     ,nQ2,istruc,ibeam,itm,iht,inuke
     &     ,foldname,comment,parfile,linlog

      character*80  nstruc(0:9),nbeam(0:9),ntm(0:9),nht(0:9)
      character*80  nnuke1(0:9),nnuke2(0:9),nnuke3(0:9)
      data nstruc / 'FL','F1','F2','F3',6*'  '    /
      data nbeam  / ' ','e','nu','nub',6*' '      /
      data ntm    / '00','GP','CF','xi',6*' '     /
      data nht    / 'LT','HT'          ,8*' '     /
      data nnuke1 / 'free','dmc','wba','lc',6*' ' / 
      data nnuke2 / ' ','_KP','_MST',7*''         /
      data nnuke3 / ' ','_bj',8*''                /

      outfile = 'no name'
      oldlen = 1
      len = 1

      print*, '*11'
      call addname(foldname,outfile,oldlen)
      print*, '*12'
      oldlen=oldlen-1
      cdum='/mom'
      call addname(cdum,outfile,oldlen)

      call addname(nstruc(istruc),outfile,oldlen)
      call addname(nbeam(ibeam),outfile,oldlen)
      call addname(ntm(itm),outfile,oldlen)
      call addname(nht(iht),outfile,oldlen)
      print*, '*13'

      !call split_nuke(inuke,ibj,ioff,ismear)
      !print*, '*14'
      !call addname(nnuke1(ismear),outfile,oldlen)
      !call addname(nnuke2(ioff)  ,outfile,oldlen)
      !call addname(nnuke3(ibj)   ,outfile,oldlen)

      if (linlog.eq.'log') then
         mmin = 10d0**Q2min
         mmax = 10d0**Q2max
      else
         mmin = Q2min
         mmax = Q2max
      end if
      !call inttochar(nint(mmin),car,min)
      !call inttochar(nint(mmax),car1,min1)
      write(car,'(I4.4)') nint(mmin)
      write(car1,'(I4.4)') nint(mmax)
      len = oldlen+24-min-min1
      outfile(oldlen:len)='QQ_'//trim(car)
     &     //'_'//trim(car1)//'.'
      oldlen = len + 1

      cdum='log'
      if (linlog.eq.'log') call addname(cdum,outfile,oldlen)

      call addname(parfile,outfile,oldlen)
      if (outfile(oldlen-5:oldlen-2).eq.'.par') oldlen=oldlen-4

      call addname(comment,outfile,oldlen)
 
      outlen = oldlen-2

      return
      end

*-------------

      subroutine addname(string,outfile,oldlen)

      implicit none

      character string*80,outfile*150
      integer oldlen,len,ll

      call trmstr(string,ll)
      len = oldlen+ll
      if (ll.gt.0) then
         outfile(oldlen:len)=string(1:ll)//'.'
         oldlen=len+1
      else
         oldlen = len
      end if

      return
      end

