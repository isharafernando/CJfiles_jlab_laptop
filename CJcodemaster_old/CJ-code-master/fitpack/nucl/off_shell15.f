C  This File contains code necessary for the flavor modified off shell
C  corrections.  These off shell corrections follow the same general
C  model as the modified Kulagin Petti corrections but modified even
C  further to allow different flavors of quarks.
C
C  Created 7/10/13 by L.T.Brady [lucas_brady@hmc.edu]
C  Currently supports valence quarks, sea quarks (u and d), and gluons
C  Implemented in such a way that it will propogate correctly through
C  any possible observable in DIS13.f
C
C  Usage:
C        Access this option in the input file by using an inuke of the
C        form:   inuke = FED5BA   (i.e. include a 5 in the C space)
C        
C        Before accessing anything else, call set_offshell with the
C        correct parameters and .true. in the last spot
C
C IMPORTANT: whenever you are done with calculating offshell corrections
C            call set_offshell again with .false. in the last flag
C            this turns off the off shell corrections and reverts
C            DIS10 to its original state.  Things will be bad otherwise



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_qder(ndrv,x,iflavor,qder)
C
C     A specific combination of things including the quark derivative
C     this combination is needed for the offshell correction model
C     specifically it calculates
C     (1/q)*\pder{q}{x}*h(x)
C     iflavor is an integer flag indicating which flavor to use:
C           1     =     valence
C           2     =     sea
C           3     =     gluons
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit None
      integer ndrv,iflavor
      double precision x,qder
      double precision h, delt
      double precision fdv, fuv, fubpdb, fcnfg
      double precision a0,a1,a2,a3,a4
      double precision b,c,a0d,a1d,a2d,a3d,a4d
      
      call get_h(ndrv,x,iflavor,h)
      delt = 1.0d-4
c
c  Don't allow x+delt to be greater than 1.d0
c
      if(1.0-x.le.delt)delt=-delt

      if (iflavor.eq.1) then        ! valence
         qder=(fdv(x+delt)+fuv(x+delt))/(fdv(x)+fuv(x))-1.
         qder=qder/delt
c
c  fsupdf returns xq, xg, etc the -1./x corrects the derivative to 
c  1/q dq/dx etc
c
         qder=h*(qder-1./x)
      else if (iflavor.eq.2) then   ! sea
         qder=fubpdb(x+delt)/fubpdb(x)-1.
         qder=qder/delt
         qder=h*(qder-1./x)
      else if (iflavor.eq.3) then   ! gluons
         qder=fcnfg(x+delt)/fcnfg(x)-1.
         qder=qder/delt
         qder=h*(qder-1./x)
      else
            write(*,*) "Error: Offshell model: iflavor out of range"
      endif

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_h (ndrv,x,iflavor,h)
C
C     Computes h(x) which is a set of constants needed for the fmKP offshell
C     corrections.  This function relies on sbar, the average spectator mass squared
C     outputs the result as h
C     iflavor is an integer flag:
C           1     =     valence
C           2     =     sea
C           3     =     gluons
C     Note, h(x) here is equivalent to x(1-x)h(x) from Kulagin-Petti
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
   
      integer ndrv,iflavor
      double precision x, h
      double precision s, sv, ssq, sg, lambda
      
      double precision Mp
      parameter(Mp=0.938272D0)

      ! fmKP PDF-dependent parameters
      double precision cv(0:35), csq(0:35), cg(0:35), xlambda(0:35)
      common/C_array/cv, csq, cg, xlambda
      
      
      call get_offshell_s(sv,ssq,sg,lambda)
      if (iflavor.eq.1) then
            s=sv
      else if (iflavor.eq.2) then
            s=ssq
      else if (iflavor.eq.3) then
            s=sg
      else
            write(*,*) "Error(get_h): iflavor out of range", iflavor
      endif
      
      lambda = xlambda(ndrv)  ! needs the "current" xlambda stored in C_array
      h=x*(1-x)*((lambda*s)/Mp**2 + (1 - lambda)*(1 - x))/
     &  (-(s/Mp**2) + (1 - x)**2)

      !print*, '* ndrv,lambda,iflv,s =', ndrv,lambda,iflavor,s

     
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_delta (ndrv,x,dfv,dfsq,dfg,dfsb,dfcb,dfbb)
C
C     determines the offshell correction to the valence quark
C     input parameters are a specific x, s, c, and lambda
C     outputs in dfv.  This was made in Mathematica
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Implicit None
      double precision x,dfv,dfsq,dfg,dfsb,dfcb,dfbb
      double precision sv, ssq, sg, lambda, old_lambda
      double precision C,x0
      double precision logis
C     Determine whether we need this offshell correction      
      logical onoff
      integer ndrv
      data old_lambda/-1.0d9/

      double precision OFF_KP
      
      ! fmKP PDF-dependent parameters
      double precision cv(0:35), csq(0:35), cg(0:35), xlambda(0:35)
      common/C_array/cv, csq, cg, xlambda

      ! KP parameters
      double precision x1(0:35),x_1
      common/KP_param/x1,x_1

      ! current fit parameters --- NOT TRUE!!!!
      double precision xc(100)  ! used in the fixed par position 
      common/curpar/xc           ! implementation only; not needed w/ findoff
      
      integer ioff_flag
      common/offshell_flag/ioff_flag

      double precision offpar(14)
      
      !save old_lambda
      !save c_v,c_sq,c_g

c
c modified get_delta routines --jfo 3/24/15
c allows for variable off_shell strength to be fitted
c combines all 'get_delta' type calls into one routine
c
c [aa 9/25/15] added KP-like model (parametrized PDF offshell mod)
c For offshell models that do not modify the PDFs, sets all delta_q = 0      


      ! OCS (formerly fmKP) model
      if (ioff_flag.ge.5.and.ioff_flag.le.7) then
      
c      call get_offshell_on(onoff)
      
c      call get_offshell_c(c_v,c_sq,c_g)
           
c      if (onoff) then 
c         call set_offshell_par
c         call get_offshell_s(sv,ssq,sg,lambda)
c         if(lambda.ne.old_lambda)then
c            old_lambda=lambda
c            call init_offshell
c         endif

         call get_logistic(x, logis)

         call get_qder(ndrv,x,1,dfv)
         dfv= dfv + cv(ndrv)
         dfv = dfv*logis
         call get_qder(ndrv,x,2,dfsq)
         dfsq= dfsq + csq(ndrv)
         dfsq = dfsq*logis
         call get_qder(ndrv,x,3,dfg)
         dfg= dfg + cg(ndrv)
         dfg = dfg*logis

         dfsb=dfsq
         dfcb=0.
         dfbb=0.

      ! KP-like model (parametrized PDF off-shell coefficient)
      else if (ioff_flag.eq.8) then 
         
         call get_off(offpar)
         C=offpar(1)
         x0=offpar(2)
         dfv = C*(x-x0)*(x-x1(ndrv))*(1.+x0-x)

         dfsq = dfv
         dfg = dfv
         dfsb = dfv
         dfcb = 0.
         dfbb = 0.
         
      ! no offshell mod at PDF level needed
      else 

         dfv = 0.d0
         dfsq=0.d0
         dfg=0.d0
         dfsb=0.d0
         dfcb=0.d0
         dfbb=0.d0

      endif


      return
      end


cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_v (x,dfv)
cC
cC     determines the offshell correction to the valence quark
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfv.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfv
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
cC     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,1,dfv)
c          dfv= dfv + c_v
c          call get_logistic(x, logis)
c          dfv = dfv*logis
c      else
c          dfv = 1.D0
c      endif
c      
c      return
c      end



cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_sq (x,dfsq)
cC
cC     determines the offshell correction to the sea quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfsq.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfsq
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
cC     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,2,dfsq)
c          dfsq= dfsq + c_sq
c          call get_logistic(x, logis)
c          dfsq = dfsq*logis
c      else
c          dfsq = 1.D0
c      endif
c
c     return
c     end



cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     subroutine get_delta_g (x,dfg)
cC
cC     determines the offshell correction to the gluons
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.  This was made in Mathematica
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfg
c      double precision sv, ssq, sg, lambda
c      double precision c_v, c_sq, c_g
c      double precision logis
C     Determine whether we need this offshell correction      
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      call get_offshell_c(c_v,c_sq,c_g)
c      call get_offshell_s(sv,ssq,sg,lambda)
c      
c      if (onoff) then 
c          call get_qder(x,3,dfg)
c          dfg= dfg + c_g
c          call get_logistic(x, logis)
c          dfg = dfg*logis
c      else
c          dfg = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_cb (x,dfcb)
cC
cC     determines the offshell correction to the charm quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.
cC     Not yet implemented fully - assumes offshell charm contribution is zero
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfcb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          dfcb = 0.d0
c      else
c          dfcb = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_sb (x,dfsb)
cC
cC     determines the offshell correction to the strange quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfsb.
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfsb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          call get_delta_sq(x,dfsb)
c      else
c          dfsb = 1.D0
c      endif
c
c      return
c      end
      
      
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      subroutine get_delta_bb (x,dfbb)
cC
cC     determines the offshell correction to the bottom quarks
cC     input parameters are a specific x, s, c, and lambda
cC     outputs in dfg.
cC     Not yet implemented fully - assumes offshell bottom contribution is zero
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      Implicit None
c      double precision x,dfbb
c      logical onoff
c      call get_offshell_on(onoff)
c      
c      if (onoff) then 
c          dfbb = 0.d0
c      else
c          dfbb = 1.D0
c      endif
c
c      return
c      end
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_offshell_on (onoff)
C     returns the current value of offshell_on
C     onoff - Logical - whether or not the offshell corrections
C                          are on or off.  If true, then DIS10 will 
C                          compute the off shell corrections instead of
C                          the normal DIS results
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      Logical onoff
      logical offshell_on
C     If the offshell common block has not bee initialized offshell_on=false
      data offshell_on /.false./
      common /offshell/ offshell_on
      
      onoff = offshell_on
      
      return
      end
      
      
      
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$      subroutine get_offshell_c (cv, csq, cg)
c$$$C     returns the current value of offshell normalization constants
c$$$C     All the outputs are double precision
c$$$CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c$$$      implicit none
c$$$      double precision cv, csq, cg
c$$$      double precision c_v, c_sq, c_g
c$$$      common /C_block/ c_v, c_sq, c_g
c$$$      
c$$$      cv=c_v
c$$$      csq=c_sq
c$$$      cg=c_g
c$$$      
c$$$      return
c$$$      end
c$$$      
c$$$


     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_offshell_s (s_v, s_sq, s_g,l)
C     returns the current value of offshell s and lambda constants
C     All the outputs are double precision
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      double precision s_v,s_sq,s_g,l
      double precision sv, ssq, sg, lambda
      common /S_block/ sv, ssq, sg, lambda
      
      s_v=sv
      s_sq=ssq
      s_g=sg
      l=lambda
      
      return
      end
    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mk_cv (ndrv,cv)
C     
C     makes the normalization constant for the valence offshell correction
C     This requires an integral which is currently computed using 
C     Simpson's rule.  Returns the normalization constant as cv
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      implicit none
c      double precision x, cv, fdv, fuv, qder, dx, coeff, N
c      integer ix,nx
      implicit real*8 (a-h,o-z)
      integer ndrv
      common/gaus32/xi(32),wi(32),nterms,xx(33)
c      nx=25 ! No. of points for Simpson's rule integration
c      dx=1.D0/nx
      
      cv=0.D0
      den=0.D0
      
c      do ix=1,nx-1 ! Simpson's rule integration
c            if (ix.Eq.1.or.ix.Eq.nx-1) then 
c                  coeff=1.
c            else if (mod(ix,2).Eq.0) then
c                  coeff=4.
c            else
c                  coeff=2.
c            endif
            
c            x=ix*dx
c            call get_qder(x,1,qder)
c            qder = qder*(1./x)*(fdv(x)+fuv(x))
c            cv = cv+coeff*qder
c            N = N + coeff*(1./x)*(fdv(x)+fuv(x))
c      enddo
      do l=1,nterms
         z=0.5*(1.+xi(l))
         x=z**3
         call get_qder(ndrv,x,1,qder)
         tmp=(fuv(x)+fdv(x))/x*3.*z**2         
         qder=qder*tmp
         cv=cv+0.5*wi(l)*qder
         den=den+0.5*wi(l)*tmp
       enddo      
c      cv=-cv/N
      cv=-cv/den      
      return
      end   
      
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine mk_cg_csq (ndrv,cv, cg, csq)
C     
C     makes the normalization constant for the gluon and sea offshell
C     corrections.  This requires an integral which is currently 
C     computed using Simpson's rule.  Returns the normalization 
C     constants as cg and csq.  Takes in cv as an input
C     At the moment, this assumes that cg=csq
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit real*8 (a-h,o-z)
      integer ndrv
      common/gaus32/xi(32),wi(32),nterms,xx(33)

      cg=0.D0
      csq=0.D0
      anum=0.D0
      denom=0.D0
      Delta=1.D0
      do l=1,nterms
         z=0.5*(1.+xi(l))
         x=z**3
         fac=0.5*wi(l)*3.*z**2
         call get_qder(ndrv,x,1,qder)
         c1=fdv(x)+fuv(x)
         anum=anum+fac*(cv+qder)*c1
         call get_qder(ndrv,x,2,qder) ! could rename qder --> qder2, use below
         c2=fubpdb(x)
         anum=anum+fac*2.*qder*c2
         call get_qder(ndrv,x,3,qder)
         c3=fcnfg(x)
         anum=anum+fac*qder*c3
         call get_qder(ndrv,x,2,qder) ! this call could be saved by returning
                                      ! qder2 instead of qder above  
         c4=(fcnvl(5,x)+.25*fcnfs(x))
         anum=anum+fac*qder*c4
         denom=denom+fac*Delta*(2.*c2+c3+c4)
      enddo
      cg=-anum/denom
      csq=Delta*cg  
      return
      end   
    
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_logistic (x, ans)
C     
C     returns a logistic curve.  This curve is designed to keep large x
C     behavior unchanged but suppress small x behavior to zero.  At the
C     moment, the function turns at x=0.15.  It is 92% by x=0.2 and is
C     7% by x=0.1.  If we want to change this function in the future,
C     the parameters are tunable.
C     A = Lower Asymptote, currently 0
C     m = mid-point of the downturn, currently x=0.15
C     B = strength of the downturn, currently 50
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      implicit none
      double precision x, ans
      double precision A, m, B
      
      A = 0D0
      m = 0.15D0
      B = 50D0
      
      ans = A + (1D0-A)/(1D0+exp(-B*(x-m)))

      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      subroutine findoff(jmax,par,ipname,pwate)
*     finds which parameters are for off-shell corrections (if any)
*     and stores their index in noff. If no parameter noff(i)=0.

      implicit none
      integer jmax
      character*10 ipname(100)
      double precision par(100),pwate(100)

      integer i, j
      character*2 num

*     Off-shell model parameter indexes
      integer maxpar
      parameter(maxpar=14)
      integer noff(maxpar) ! max of 14 off-shell parameters 
      common/offpar/noff

      do i=1,maxpar
         noff(i)=0
         write(num,'(I0)') i
         do j=1,jmax
            if(ipname(j).eq.'off'//trim(num)) then
               noff(i)=j
              end if
         end do
      end do

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
      subroutine get_off(offpar)
*     returns the current value of the offshell parameters.
*     If a parameter is not present in the original input file
*     it returns -1D99 for that parameter. The calling routine shoudl take
*     care to test for missing parameters by checking against this value 

      implicit none

      double precision offpar(14) ! Note: size should be = maxpar below

      integer i

*     Off-shell model parameter indexes
      integer maxpar
      parameter(maxpar=14)
      integer noff(maxpar) ! max of 14 off-shell parameters 
      common/offpar/noff
*     Current fit parameters  ---- NOT TRUE!!! 
      double precision xc(100)
      common/curpar/xc 

      do i=1,maxpar
         if (noff(i).gt.0) then
            offpar(i) = xc(noff(i))
         else
            offpar(i) = -1D99
         end if
      end do

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine set_offshell_on (onoff)
C
C     This function sets on and off the off-shell corrections to the PDFs
C     as needed by the user. 
c
C        onoff - Logical - sets whether or not the offshell corrections
C                          are on or off.  If true, then DIS10 will 
C                          compute the off shell corrections instead of
C                          the normal DIS results
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC NOTE CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     onoff should be set to false using this routine as soon as you are
C     done calculating the offshell corrections.  It WILL interfere in
C     the rest of the code otherwise.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Implicit None
      logical onoff

      logical offshell_on
      common /offshell/ offshell_on

      offshell_on = onoff

      ! Sets the needed offshell parameters 
      ! (Needed for the free-floating option, to make sure the current dRN
      ! is used - negligible overhead otherwise)
      call set_offshell_par   ! really needed? For ioff=7,8 the intializations 
                              ! are done elsewhere; maybe needed for ioff=4,5,6 ?
         
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine set_offshell_par
C
C     Sets the offshell model constants
C
C     *** fmKP model: sets the values of sbar for each parton,
C     lambda, and the normalization
C
C     *** 
C     
C
C     It retrieves the nuclear model previously stored in a common block
C     and uses the following switches. (NOTE: it assumes the inuke code
C     was stored using 'set_nuke' in the calling routines.)
C
C        wfn - Integer - Denotes the wavefunction being used
C                    0 = Paris
C                    1 = AV18
C                    3 = CDBonn   
C                    4 = WJC-1  
C                    5 = WJC-2   
C
C        ioff - offshell model  4 = mKP, 5 = fmKP, 6 = "negative" fmKP
C                               7 = free floating fmKP (iDRN ignored)
C 
C        idRN - Integer - Denotes the size of the Nuclear Radius decrease
C
C                      mkP          fmKP      neg fmKP
C
C                    1 = 0.3%     1 = 0.1%    1 = -0.1% 
C                    2 = 0.6%     2 = 0.2%    2 = -0.2% 
C                    3 = 0.9%     3 = 0.3%    3 = -0.3% 
C                    4 = 1.2%     4 = 0.4%    4 = -0.4% 
C                    5 = 1.5%     5 = 0.5%    5 = -0.5% 
C                    6 = 1.8%     6 = 0.6%    6 = -0.6% 
C                    7 = 2.1%     7 = 0.7%    7 = -0.7% 
C                    8 = 2.4%     8 = 0.8%    8 = -0.8% 
C                    9 = 2.7%     9 = 0.9%    9 = -0.9% 
C                    0 = 3.0%     0 = 0%      0 =  0% 
c
      Implicit None
      integer inuke,ilam,ishad,idRN,ioff,iwfn,wfn,itgt
      double precision dp2_array(5)
      double precision dRN_array_mKP(10),dRN_array_fmKP(10)
      double precision offpar(14)
      double precision c_uv_array(5,10),c_dv_array(5,10)
      double precision c_g_array(5,10),c_s_array(5,10),c_v_array(5,10)
      double precision sv, ssq, sg, lambda
      common /S_block/ sv, ssq, sg, lambda
      data dp2_array/-3.60D-2, -4.26D-2, -6.21D-2, -4.86D-2, -4.29D-2/
      data dRN_array_mKP/0.3D-2, 0.6D-2, 0.9D-2, 1.2D-2, 1.5D-2, 1.8D-2,
     &                   2.1D-2, 2.4D-2, 2.7D-2, 3.0D-2/
      data dRN_array_fmKP/0.1D-2, 0.2D-2, 0.3D-2, 0.4D-2, 0.5D-2, 0.6D-2,
     &                    0.7D-2, 0.8D-2, 0.9D-2, 0D-2/

      double precision old_lambda
      data old_lambda/-2000d0/
      save old_lambda
      
      ! Retrieves and unfolds the nuclear model code
      call get_nuke(inuke)
      call split_nuke(inuke,ilam,ishad,idRN,ioff,wfn,itgt)

      ! if fmKP off-shell model requested, sets its parameters
      if (ioff.ge.4.and.ioff.le.7) then

         if (wfn.Eq.2) iwfn=1
         if (wfn.Eq.1) iwfn=2
         if (wfn.Eq.3) iwfn=3
         if (wfn.Eq.4) iwfn=4
         if (wfn.Eq.0) iwfn=5
         
         if (idRN.eq.0) then
            idRN = 10
         endif
         
         ssq = 5.5              ! 2.907 commented GRV, using GJR
         sg  = 8.0              ! 1.0
         sv =  2.2              ! 2.4

         if (ioff.eq.4) then    ! mKP off-shell parameter grid 
            lambda = -2*dRN_array_mKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.5) then ! finer grid for fmKP model
            lambda = -2*dRN_array_fmKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.6) then ! "negative" fmKP model (plus sign below)
            lambda = 2*dRN_array_fmKP(idRN)/dp2_array(iwfn)
         else if (ioff.eq.7) then ! fmKP with running offshell param
            call get_off(offpar) 
            lambda = -2*offpar(1)/dp2_array(iwfn)
         end if

         if (isnan(offpar(1))) then
            print*, 'ERROR(set_offshell_par): offpar(1)=NaN'
            stop
         end if

      end if

      return 
      end




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine init_offshell_fmKP(ndrv)

C     Initializes the fmKP off-shell model. 
C     (Does no harm if called for other models)

      implicit none

      integer ndrv

      double precision c_v, c_sq, c_g
      common /C_block/ c_v, c_sq, c_g

      double precision s_v, s_sq, s_g, lambda
      common/S_block/s_v, s_sq, s_g, lambda                        

      double precision cv(0:35), csq(0:35), cg(0:35), xlambda(0:35)
      common/C_array/cv, csq, cg, xlambda                          
                                                                   
      ! for debugging purposes only
      integer count
      data count/0/
      save count


      ! Switches off offshell effects as a precaution, so the right
      ! sum rules can be applied
      !   (The user is always required to switch them on only when needed, 
      !   and off as done with off-shell PDF calculations as to not mess up
      !   possible subsequent on-shell PDF evaluations.)
      call set_offshell_on(.false.) 

      ! Sets the needed (PDF-independent) fmKP coefficients 
      call set_offshell_par 
      xlambda(ndrv)=lambda
      
      ! calculuates PDF-dependent coefficients, store in common block
      call mk_cv(ndrv,c_v)
      call mk_cg_csq(ndrv,c_v,c_g,c_sq) 
      cv(ndrv)=c_v
      csq(ndrv)=c_sq
      cg(ndrv)=c_g

      write(6,'(A,I4,X,4ES13.5)'),'* fmKP: '
     &     ,ndrv,xlambda(ndrv),cv(ndrv),csq(ndrv),cg(ndrv)

      
      count = count+1   ! counter is for dubugging purposes only
C      print*, '* init_offshell call number', count, c_v,c_g,c_sq
      return
      end


*******************************************************
      subroutine init_offshell_KP(xc)
*     Initializes the KP offshell parametrization
*       x1 : determined by valence normalization
*       x_1 = ???
*     [jfo: June 2015]

      implicit real*8 (a-h,o-z)
      dimension xc(100) 

      double precision offpar(14)

      common/gaus32/xi(32),wi(32),nterms,xx(33)
      common/KP_param/x1(0:35),x_1

      
      !x0=xc(34)      
      call get_off(offpar)
      x0=offpar(2)   ! Note the sum rule for x1 is independent of the
                     ! normalization C of the offshell function delta_q(x)
      
      anum=0.
      den=0.
      do l=1,nterms
         z=0.5*(1.+xi(l))
         x=z**3
         tmp=(fuv(x)+fdv(x))/x*3.*z**2
         tmp_den=tmp*(x-x0)*(1.+x0-x)
         tmp_anum=tmp_den*x
*         print*,z,x,tmp,tmp_anum,tmp_den
         den=den+0.5*wi(l)*tmp_den
         anum=anum+0.5*wi(l)*tmp_anum
      enddo
      x_1=anum/den
*      print*,anum,den,x_1
      return
      end


***************************************************************
      subroutine setoffshell(ndrv,xc)
*     Initializes various constants for offshell corrections
*
*       ndrv    = (i) free parameter index for derivatives (0=central value)
*       xc(100) = (dp) fit parameters 
*
*     (aa, 08 Aug 2015)

      implicit none

      integer ndrv
      double precision xc(100)

      double precision offpar(14)  ! needed for debugging only

      integer ioff_flag
      common/offshell_flag/ioff_flag

c$$$      double precision c_v, c_sq, c_g
c$$$      common/C_block/c_v, c_sq, c_g
c$$$
c$$$      double precision cv(0:35), csq(0:35), cg(0:35), xlambda(0:35)
c$$$      common/C_array/cv, csq, cg, xlambda
c$$$
c$$$      double precision s_v, s_sq, s_g, lambda
c$$$      common/S_block/s_v, s_sq, s_g, lambda 

      double precision x1(0:35),x_1
      common/KP_param/x1,x_1

!     Only calculates parameters that are PDF-dependent;
!     for models with PDF-independent parameters this call has no effect
      
      ! fmKP model 
      if(ioff_flag.ge.5.and.ioff_flag.le.7)then  ! not ioff=4 b/c 4 is mKP at F2 level
         call init_offshell_fmKP(ndrv) ! initializes the offshell mods for PDFs
c$$$         cv(ndrv)=c_v
c$$$         csq(ndrv)=c_sq
c$$$         cg(ndrv)=c_g
c$$$         xlambda(ndrv)=lambda
c$$$         write(6,'(A,I4,X,4ES13.5)'),'* fmKP: '
c$$$     &        ,ndrv,lambda,cv(ndrv),csq(ndrv),cg(ndrv)
     
c     MMHT-like (Parametrizes F2D/F2N-1)
c             
c     ! Normalizations and offshell strength are now stored for each value
c     ! of ndrv. They may be accessed through the common block C_array
c     ! jfo 3/25/15 
      else if(ioff_flag.eq.3)then
         call init_offshell_KP(xc) !Applies valence quark constraint
         x1(ndrv)=x_1
         write(6,'(A,I4,X,4ES13.5)'),'* MMHT-like: ', ndrv,x1(ndrv),x_1

!     KP-like at PDF level
!     ! aa 9/24/15
      else if(ioff_flag.eq.8)then
         call init_offshell_KP(xc) !Applies valence quark constraint
         x1(ndrv)=x_1
         call get_off(offpar)
         write(6,'(A,I4,X,4ES13.5)'),'* KP-like:'
     &        ,ndrv,offpar(1),offpar(2),x_1
         
      end if 

      return 
      end


C ***********************************************************************
      FUNCTION OFF_KP(x,xc,ndrv)
C
C  Analytic parametrization of nucleon off-shell correction from
C    Kulagin-Petti fit
C
C  Defined such that F2d = F2d(conv) + del^off F2d
C       with OFF_KO = del^off F2d / F2d
C                   = off-shell correction w.r.t. F2d
C
C  Parameters defined for x = nucleon scaling variable in [0,2]. 
C
C  Ref: Kulagin, Petti, NPA 765, 126 (2006).
C
C  July 2010.
C
C  Modified to allow fitting the KP parameters - 5/13/15 jfo
C ***********************************************************************
      IMPLICIT NONE
      REAL*8  OFF_KP,x,xc(100)
      REAL*8  aoff(0:4),boff(0:13),coff(0:5)
      integer ndrv  

      double precision x1(0:35),x_1
      common/KP_param/x1,x_1

      DATA    aoff /-0.010675D0, 0.27315D0, -0.96047D0, 1.2396D0,
     &     -0.71897D0/
      DATA    boff /-0.11185D-1, 0.29485D0, -1.3024D0, 4.1337D0,
     &     -17.047D0, 66.378D0, -184.54D0, 293.49D0,
     &     -64.581D0, -785.17D0, 1779.8D0, -1908.6D0,
     &     1066.9D0, -250.54D0/
      DATA    coff /3.9913D0, 4.3786D0, 2.4799D0, 2.5043D0,
     &     -3.9996D0, 0.018932D0/

       OFF_KP = 0.D0

C...13th order polynomial fit valid for x < 0.91
c       OFF_KP = boff(0) + boff(1)*x + boff(2)*x**2 + boff(3)*x**3
c     &        + boff(4)*x**4 + boff(5)*x**5 + boff(6)*x**6
c     &        + boff(7)*x**7 + boff(8)*x**8 + boff(9)*x**9
c     &        + boff(10)*x**10 + boff(11)*x**11 + boff(12)*x**12
c     &        + boff(13)*x**13
c
c       IF (x.GT.0.85D0) RETURN
C...4th order polynomial fit valid for x < 0.86
c       OFF_KP = aoff(0) + aoff(1)*x + aoff(2)*x**2 + aoff(3)*x**3
c     &        + aoff(4)*x**4

C...6-parameter fit  (courtesy of Simona Malace)
c      OFF_KP = coff(0) + coff(1)*x + coff(2)*x**coff(3)
c     &     + coff(4)*DEXP(x) + coff(5)/DLOG(x)

      OFF_KP=xc(35)*(x-xc(34))*(x-x1(ndrv))*(1.+xc(34)-x)

      RETURN
      END



C ***********************************************************************
        FUNCTION OFF_mKP (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]).
C  Fit by O. Hen (5/24/11 - 6/6/11).
C  Fortran code by W. Melnitchouk (5/27/11).
C
C  wfn: 1 (AV18)
C       2 (CD-Bonn)
C       3 (WJC-1)
C       4 (WJC-2)
C
C  dRN: 0 (0.0%)        [% change in nucleon radius in the deuteron]
C       1 (0.3%)
C       2 (0.6%)
C       3 (0.9%)
C       4 (1.2%)
C       5 (1.5%)
C       6 (1.8%)
C       7 (2.1%)
C       8 (2.4%)
C       9 (2.7%)
C       10 (3.0%)
C       
C ***********************************************************************
        IMPLICIT NONE
        REAL*8  OFF_mKP, x
        INTEGER wfn, dRN
        REAL*8  p(0:9), q(0:3)

        OFF_mKP = 0.D0
        IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
        IF (dRN.LT.0 .OR. dRN.GT.10) RETURN

! .......................................................................
        IF (wfn.EQ.1) THEN              ! AV18
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00196328D0
            p(1) = -0.0406043D0
            p(2) = 0.0854386D0
            p(3) = 0.598314D0
            p(4) = -5.18768D0
            p(5) = 19.3408D0
            p(6) = -40.8099D0
            p(7) = 50.0192D0
            p(8) = -33.2016D0
            p(9) = 9.25861D0
            q(0) = 1.70778D0
            q(1) = -3.8263D0
            q(2) = 2.04641D0
            q(3) = 0.132957D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00337272D0
            p(1) = -0.0589391D0
            p(2) = 0.204111D0
            p(3) = -0.609015D0
            p(4) = 1.90664D0
            p(5) = -5.05949D0
            p(6) = 9.39415D0
            p(7) = -10.967D0
            p(8) = 7.17457D0
            p(9) = -2.01069D0
            q(0) = -3.19753D0
            q(1) = 9.92432D0
            q(2) = -10.1881D0
            q(3) = 3.44211D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00479856D0
            p(1) = -0.0789437D0
            p(2) = 0.365415D0
            p(3) = -2.28553D0
            p(4) = 11.7184D0
            p(5) = -38.5546D0
            p(6) = 77.7781D0
            p(7) = -93.3614D0
            p(8) = 61.2492D0
            p(9) = -16.9552D0
            q(0) = -2.58892D0
            q(1) = 5.16641D0
            q(2) = -1.68202D0
            q(3) = -1.00944D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00624526D0
            p(1) = -0.101004D0
            p(2) = 0.578828D0
            p(3) = -4.53405D0
            p(4) = 24.8387D0
            p(5) = -83.1084D0
            p(6) = 168.235D0
            p(7) = -201.705D0
            p(8) = 131.896D0
            p(9) = -36.336D0
            q(0) = 8.43383D0
            q(1) = -34.09D0
            q(2) = 44.9658D0
            q(3) = -19.5377D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00771879D0
            p(1) = -0.125631D0
            p(2) = 0.85681D0
            p(3) = -7.48852D0
            p(4) = 42.0315D0
            p(5) = -141.24D0
            p(6) = 285.726D0
            p(7) = -341.749D0
            p(8) = 222.735D0
            p(9) = -61.1069D0
            q(0) = 38.4052D0
            q(1) = -135.585D0
            q(2) = 159.821D0
            q(3) = -63.0094D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.00922543D0
            p(1) = -0.153372D0
            p(2) = 1.21284D0
            p(3) = -11.2957D0
            p(4) = 64.1424D0
            p(5) = -215.764D0
            p(6) = 435.837D0
            p(7) = -520.011D0
            p(8) = 337.89D0
            p(9) = -92.3594D0
            q(0) = 102.447D0
            q(1) = -348.309D0
            q(2) = 395.805D0
            q(3) = -150.485D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0107738D0
            p(1) = -0.184965D0
            p(2) = 1.66493D0
            p(3) = -16.1495D0
            p(4) = 92.2796D0
            p(5) = -310.337D0
            p(6) = 625.768D0
            p(7) = -744.838D0
            p(8) = 482.61D0
            p(9) = -131.475D0
            q(0) = 228.085D0
            q(1) = -761.2D0
            q(2) = 848.715D0
            q(3) = -316.366D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0123752D0
            p(1) = -0.22138D0
            p(2) = 2.23652D0
            p(3) = -22.3008D0
            p(4) = 127.866D0
            p(5) = -429.628D0
            p(6) = 864.673D0
            p(7) = -1026.79D0
            p(8) = 663.511D0
            p(9) = -180.187D0
            q(0) = 467.305D0
            q(1) = -1541.83D0
            q(2) = 1698.6D0
            q(3) = -625.146D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.014041D0
            p(1) = -0.2636D0
            p(2) = 2.95198D0
            p(3) = -30.0151D0
            p(4) = 172.431D0
            p(5) = -578.727D0
            p(6) = 1162.64D0
            p(7) = -1377.62D0
            p(8) = 888.008D0
            p(9) = -240.454D0
            q(0) = 923.217D0
            q(1) = -3021.86D0
            q(2) = 3301.08D0
            q(3) = -1203.91D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0157866D0
            p(1) = -0.312963D0
            p(2) = 3.84398D0
            p(3) = -39.6446D0
            p(4) = 227.989D0
            p(5) = -764.276D0
            p(6) = 1532.74D0
            p(7) = -1812.46D0
            p(8) = 1165.61D0
            p(9) = -314.77D0
            q(0) = 1814.91D0
            q(1) = -5904.72D0
            q(2) = 6408.9D0
            q(3) = -2321.16D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0176352D0
            p(1) = -0.371422D0
            p(2) = 4.9594D0
            p(3) = -51.687D0
            p(4) = 297.343D0
            p(5) = -995.421D0
            p(6) = 1992.79D0
            p(7) = -2351.74D0
            p(8) = 1509.02D0
            p(9) = -406.45D0
            q(0) = 3653.87D0
            q(1) = -11829.9D0
            q(2) = 12773.5D0
            q(3) = -4600.46D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.2) THEN         ! CD-Bonn
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00162711D0
            p(1) = -0.0327193D0
            p(2) = 0.0595688D0
            p(3) = 0.492763D0
            p(4) = -3.90178D0
            p(5) = 14.0632D0
            p(6) = -29.1224D0
            p(7) = 35.2469D0
            p(8) = -23.1699D0
            p(9) = 6.40719D0
            q(0) = -1.42342D0
            q(1) = 5.78443D0
            q(2) = -7.67182D0
            q(3) = 3.35475D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00303191D0
            p(1) = -0.0496808D0
            p(2) = 0.154688D0
            p(3) = -0.497176D0
            p(4) = 1.91967D0
            p(5) = -5.86469D0
            p(6) = 11.635D0
            p(7) = -13.9418D0
            p(8) = 9.17025D0
            p(9) = -2.55081D0
            q(0) = -0.596239D0
            q(1) = 1.40214D0
            q(2) = -0.882721D0
            q(3) = 0.05642D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00444406D0
            p(1) = -0.0675313D0
            p(2) = 0.273155D0
            p(3) = -1.74564D0
            p(4) = 9.23849D0
            p(5) = -30.8037D0
            p(6) = 62.4092D0
            p(7) = -74.9276D0
            p(8) = 49.0612D0
            p(9) = -13.5357D0
            q(0) = 4.81023D0
            q(1) = -18.1571D0
            q(2) = 22.6966D0
            q(3) = -9.44371D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00586546D0
            p(1) = -0.0864343D0
            p(2) = 0.418993D0
            p(3) = -3.29639D0
            p(4) = 18.3067D0
            p(5) = -61.5928D0
            p(6) = 124.867D0
            p(7) = -149.658D0
            p(8) = 97.7382D0
            p(9) = -26.8755D0
            q(0) = 17.3269D0
            q(1) = -61.129D0
            q(2) = 72.0038D0
            q(3) = -28.3804D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00729865D0
            p(1) = -0.106607D0
            p(2) = 0.597462D0
            p(3) = -5.20583D0
            p(4) = 29.4458D0
            p(5) = -99.2915D0
            p(6) = 201.094D0
            p(7) = -240.554D0
            p(8) = 156.726D0
            p(9) = -42.9721D0
            q(0) = 40.7797D0
            q(1) = -139.933D0
            q(2) = 160.481D0
            q(3) = -61.6055D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.008746D0
            p(1) = -0.128256D0
            p(2) = 0.813724D0
            p(3) = -7.53053D0
            p(4) = 42.9835D0
            p(5) = -144.996D0
            p(6) = 293.277D0
            p(7) = -350.176D0
            p(8) = 227.651D0
            p(9) = -62.2584D0
            q(0) = 81.0211D0
            q(1) = -273.524D0
            q(2) = 308.599D0
            q(3) = -116.49D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0102108D0
            p(1) = -0.151667D0
            p(2) = 1.0747D0
            p(3) = -10.3452D0
            p(4) = 59.3469D0
            p(5) = -200.119D0
            p(6) = 404.201D0
            p(7) = -481.759D0
            p(8) = 312.554D0
            p(9) = -85.2724D0
            q(0) = 147.141D0
            q(1) = -491.286D0
            q(2) = 548.022D0
            q(3) = -204.411D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0116974D0
            p(1) = -0.177204D0
            p(2) = 1.38924D0
            p(3) = -13.744D0
            p(4) = 79.0702D0
            p(5) = -266.411D0
            p(6) = 537.302D0
            p(7) = -639.272D0
            p(8) = 413.922D0
            p(9) = -112.667D0
            q(0) = 253.515D0
            q(1) = -839.606D0
            q(2) = 928.653D0
            q(3) = -343.264D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.0132096D0
            p(1) = -0.205202D0
            p(2) = 1.76563D0
            p(3) = -17.8183D0
            p(4) = 102.683D0
            p(5) = -345.645D0
            p(6) = 696.108D0
            p(7) = -826.844D0
            p(8) = 534.376D0
            p(9) = -145.138D0
            q(0) = 423.458D0
            q(1) = -1393.59D0
            q(2) = 1531.15D0
            q(3) = -561.925D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0147536D0
            p(1) = -0.236188D0
            p(2) = 2.21657D0
            p(3) = -22.7029D0
            p(4) = 130.947D0
            p(5) = -440.312D0
            p(6) = 885.492D0
            p(7) = -1050.08D0
            p(8) = 677.422D0
            p(9) = -183.603D0
            q(0) = 695.863D0
            q(1) = -2278.36D0
            q(2) = 2489.67D0
            q(3) = -908.341D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0163345D0
            p(1) = -0.270612D0
            p(2) = 2.75336D0
            p(3) = -28.5225D0
            p(4) = 164.588D0
            p(5) = -552.843D0
            p(6) = 1110.3D0
            p(7) = -1314.67D0
            p(8) = 846.662D0
            p(9) = -229.017D0
            q(0) = 1138.1D0
            q(1) = -3710.33D0
            q(2) = 4035.93D0
            q(3) = -1465.21D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.3) THEN         ! WJC-1
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00274105D0
            p(1) = -0.060672D0
            p(2) = 0.212256D0
            p(3) = 0.0219547D0
            p(4) = -3.19661D0
            p(5) = 14.5887D0
            p(6) = -33.5328D0
            p(7) = 43.4176D0
            p(8) = -30.1071D0
            p(9) = 8.74272D0
            q(0) = 7.79726D0
            q(1) = -23.5139D0
            q(2) = 23.1821D0
            q(3) = -7.38804D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00411184D0
            p(1) = -0.0770866D0
            p(2) = 0.273399D0
            p(3) = -0.512489D0
            p(4) = -0.00242581D0
            p(5) = 3.08803D0
            p(6) = -8.72114D0
            p(7) = 11.7869D0
            p(8) = -8.10611D0
            p(9) = 2.26714D0
            q(0) = -2.53615D0
            q(1) = 8.56149D0
            q(2) = -9.65752D0
            q(3) = 3.63462D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00549862D0
            p(1) = -0.0951343D0
            p(2) = 0.376461D0
            p(3) = -1.51089D0
            p(4) = 5.8943D0
            p(5) = -17.5102D0
            p(6) = 34.3868D0
            p(7) = -41.5302D0
            p(8) = 27.87D0
            p(9) = -7.98734D0
            q(0) = -12.1657D0
            q(1) = 37.7354D0
            q(2) = -38.6311D0
            q(3) = 12.9756D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.006907D0
            p(1) = -0.1153D0
            p(2) = 0.533373D0
            p(3) = -3.10283D0
            p(4) = 15.2396D0
            p(5) = -49.6885D0
            p(6) = 100.722D0
            p(7) = -122.297D0
            p(8) = 81.4779D0
            p(9) = -22.9919D0
            q(0) = -18.911D0
            q(1) = 56.7604D0
            q(2) = -55.7012D0
            q(3) = 17.6606D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00834428D0
            p(1) = -0.138215D0
            p(2) = 0.759563D0
            p(3) = -5.45509D0
            p(4) = 28.9898D0
            p(5) = -96.6165D0
            p(6) = 196.556D0
            p(7) = -237.821D0
            p(8) = 157.337D0
            p(9) = -43.9692D0
            q(0) = -18.5907D0
            q(1) = 51.8677D0
            q(2) = -45.7426D0
            q(3) = 12.1465D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.0098194D0
            p(1) = -0.16465D0
            p(2) = 1.07396D0
            p(3) = -8.77272D0
            p(4) = 48.3221D0
            p(5) = -162.2D0
            p(6) = 329.623D0
            p(7) = -397.11D0
            p(8) = 261.143D0
            p(9) = -72.4259D0
            q(0) = -3.2391D0
            q(1) = -2.99953D0
            q(2) = 19.6733D0
            q(3) = -13.9123D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0113437D0
            p(1) = -0.195584D0
            p(2) = 1.50055D0
            p(3) = -13.3152D0
            p(4) = 74.7236D0
            p(5) = -251.364D0
            p(6) = 509.657D0
            p(7) = -611.486D0
            p(8) = 400.039D0
            p(9) = -110.248D0
            q(0) = 42.4961D0
            q(1) = -157.808D0
            q(2) = 194.778D0
            q(3) = -80.1444D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.012932D0
            p(1) = -0.232291D0
            p(2) = 2.07042D0
            p(3) = -19.4176D0
            p(4) = 110.11D0
            p(5) = -370.435D0
            p(6) = 749.127D0
            p(7) = -895.409D0
            p(8) = 583.123D0
            p(9) = -159.831D0
            q(0) = 148.975D0
            q(1) = -510.974D0
            q(2) = 585.939D0
            q(3) = -224.882D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.014603D0
            p(1) = -0.276387D0
            p(2) = 2.82299D0
            p(3) = -27.5043D0
            p(4) = 156.902D0
            p(5) = -527.402D0
            p(6) = 1063.77D0
            p(7) = -1267.1D0
            p(8) = 821.851D0
            p(9) = -224.187D0
            q(0) = 378.701D0
            q(1) = -1264.44D0
            q(2) = 1410.67D0
            q(3) = -526.229D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0163792D0
            p(1) = -0.329812D0
            p(2) = 3.80591D0
            p(3) = -38.0918D0
            p(4) = 218.065D0
            p(5) = -732.077D0
            p(6) = 1472.94D0
            p(7) = -1749.06D0
            p(8) = 1130.38D0
            p(9) = -307.041D0
            q(0) = 868.106D0
            q(1) = -2857.75D0
            q(2) = 3141.02D0
            q(3) = -1153.18D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0182922D0
            p(1) = -0.395284D0
            p(2) = 5.08523D0
            p(3) = -51.8887D0
            p(4) = 297.626D0
            p(5) = -997.705D0
            p(6) = 2002.64D0
            p(7) = -2371.26D0
            p(8) = 1527.47D0
            p(9) = -413.31D0
            q(0) = 1941.33D0
            q(1) = -6332.67D0
            q(2) = 6892.95D0
            q(3) = -2504.2D0
          ENDIF
! .......................................................................
        ELSE IF (wfn.EQ.4) THEN         ! WJC-2
          IF (dRN.EQ.0) THEN            ! 0.0%
            p(0) = 0.00214184D0
            p(1) = -0.0450598D0
            p(2) = 0.110368D0
            p(3) = 0.507098D0
            p(4) = -4.94836D0
            p(5) = 18.9418D0
            p(6) = -40.4816D0
            p(7) = 50.0558D0
            p(8) = -33.4725D0
            p(9) = 9.40103D0
            q(0) = 3.4267D0
            q(1) = -9.35273D0
            q(2) = 7.94702D0
            q(3) = -1.95589D0
          ELSE IF (dRN.EQ.1) THEN       ! 0.3%
            p(0) = 0.00354093D0
            p(1) = -0.0629637D0
            p(2) = 0.216152D0
            p(3) = -0.549377D0
            p(4) = 1.27716D0
            p(5) = -2.60462D0
            p(6) = 4.14477D0
            p(7) = -4.53293D0
            p(8) = 2.93507D0
            p(9) = -0.844147D0
            q(0) = -3.57622D0
            q(1) = 11.3618D0
            q(2) = -11.9905D0
            q(3) = 4.19102D0
          ELSE IF (dRN.EQ.2) THEN       ! 0.6%
            p(0) = 0.00495676D0
            p(1) = -0.0825665D0
            p(2) = 0.365339D0
            p(3) = -2.08408D0
            p(4) = 10.2762D0
            p(5) = -33.4462D0
            p(6) = 67.3774D0
            p(7) = -81.0641D0
            p(8) = 53.405D0
            p(9) = -14.8687D0
            q(0) = -6.54542D0
            q(1) = 18.3482D0
            q(2) = -16.3313D0
            q(3) = 4.42111D0
          ELSE IF (dRN.EQ.3) THEN       ! 0.9%
            p(0) = 0.00639435D0
            p(1) = -0.104301D0
            p(2) = 0.568512D0
            p(3) = -4.21151D0
            p(4) = 22.7055D0
            p(5) = -75.76D0
            p(6) = 153.525D0
            p(7) = -184.556D0
            p(8) = 121.11D0
            p(9) = -33.5121D0
            q(0) = -1.32741D0
            q(1) = -1.99123D0
            q(2) = 9.77254D0
            q(3) = -6.67373D0
          ELSE IF (dRN.EQ.4) THEN       ! 1.2%
            p(0) = 0.00785974D0
            p(1) = -0.128688D0
            p(2) = 0.838473D0
            p(3) = -7.07025D0
            p(4) = 39.3602D0
            p(5) = -132.183D0
            p(6) = 267.806D0
            p(7) = -321.086D0
            p(8) = 209.89D0
            p(9) = -57.7897D0
            q(0) = 19.4618D0
            q(1) = -73.7158D0
            q(2) = 92.4639D0
            q(3) = -38.5673D0
          ELSE IF (dRN.EQ.5) THEN       ! 1.5%
            p(0) = 0.00936044D0
            p(1) = -0.156377D0
            p(2) = 1.19109D0
            p(3) = -10.8319D0
            p(4) = 61.2231D0
            p(5) = -205.969D0
            p(6) = 416.647D0
            p(7) = -498.125D0
            p(8) = 324.456D0
            p(9) = -88.9448D0
            q(0) = 69.1367D0
            q(1) = -240.043D0
            q(2) = 278.518D0
            q(3) = -108.141D0
          ELSE IF (dRN.EQ.6) THEN       ! 1.8%
            p(0) = 0.0109059D0
            p(1) = -0.188179D0
            p(2) = 1.64619D0
            p(3) = -15.71D0
            p(4) = 89.5156D0
            p(5) = -301.152D0
            p(6) = 608.003D0
            p(7) = -724.9D0
            p(8) = 470.617D0
            p(9) = -128.508D0
            q(0) = 172.305D0
            q(1) = -580.597D0
            q(2) = 653.83D0
            q(3) = -246.288D0
          ELSE IF (dRN.EQ.7) THEN       ! 2.1%
            p(0) = 0.0125082D0
            p(1) = -0.225133D0
            p(2) = 2.22902D0
            p(3) = -21.976D0
            p(4) = 125.784D0
            p(5) = -422.823D0
            p(6) = 851.883D0
            p(7) = -1012.99D0
            p(8) = 655.643D0
            p(9) = -178.388D0
            q(0) = 376.045D0
            q(1) = -1247.3D0
            q(2) = 1381.83D0
            q(3) = -511.62D0
          ELSE IF (dRN.EQ.8) THEN       ! 2.4%
            p(0) = 0.0141801D0
            p(1) = -0.26835D0
            p(2) = 2.96709D0
            p(3) = -29.9298D0
            p(4) = 171.756D0
            p(5) = -576.732D0
            p(6) = 1159.68D0
            p(7) = -1375.67D0
            p(8) = 887.914D0
            p(9) = -240.797D0
            q(0) = 774.973D0
            q(1) = -2544.77D0
            q(2) = 2789.44D0
            q(3) = -1021.08D0
          ELSE IF (dRN.EQ.9) THEN       ! 2.7%
            p(0) = 0.0159409D0
            p(1) = -0.319476D0
            p(2) = 3.90024D0
            p(3) = -39.9977D0
            p(4) = 229.853D0
            p(5) = -770.824D0
            p(6) = 1546.97D0
            p(7) = -1830.9D0
            p(8) = 1178.66D0
            p(9) = -318.677D0
            q(0) = 1572.33D0
            q(1) = -5126.02D0
            q(2) = 5575.96D0
            q(3) = -2024.29D0
          ELSE IF (dRN.EQ.10) THEN      ! 3.0%
            p(0) = 0.0178144D0
            p(1) = -0.380548D0
            p(2) = 5.07793D0
            p(3) = -52.7094D0
            p(4) = 303.083D0
            p(5) = -1014.97D0
            p(6) = 2033.09D0
            p(7) = -2400.96D0
            p(8) = 1541.84D0
            p(9) = -415.674D0
            q(0) = 3248.08D0
            q(1) = -10530.3D0
            q(2) = 11386.6D0
            q(3) = -4107.33D0
          ENDIF
        ENDIF

        IF (x.LE.0.9D0) THEN
          OFF_mKP = p(0) + p(1)*x + p(2)*x**2 + p(3)*x**3 + p(4)*x**4
     &            + p(5)*x**5 + p(6)*x**6 + p(7)*x**7 + p(8)*x**8
     &            + p(9)*x**9
        ELSE IF (x.GT.0.9D0) THEN
          OFF_mKP = q(0) + q(1)*x + q(2)*x**2 + q(3)*x**3
        ENDIF

        RETURN
        END

C ***********************************************************************
        FUNCTION off_mKP_fit (x,wfn,dRN)
C
C  Polynomial fit to calculated off-shell correction in mKP model
C    (see CJ11, arXiv:1102.3686 [hep-ph]), constrained to vanish at x=0
C
C  Defined such that F2d = F2d(conv) + del^off F2d 
C       with off_mKP_fit = del^off F2d / F2d
C                    = off-shell correction w.r.t. F2d
C
C  Fit by Eric Christy (6/14/12).
C  Fortran code by W. Melnitchouk (6/17/12).
C
C  wfn: 1 (AV18)
C	2 (CD-Bonn)
C	3 (WJC-1)
C	4 (WJC-2)
C
C  dRN: 0 (0.0%)	[% change in nucleon radius in the deuteron]
C	1 (0.3%)
C	2 (0.6%)
C	3 (0.9%)
C	4 (1.2%)
C	5 (1.5%)
C	6 (1.8%)
C	7 (2.1%)
C	8 (2.4%)
C	9 (2.7%)
C      10 (3.0%)
C
C ***********************************************************************
	IMPLICIT NONE
	REAL*8	off_mKP_fit, x
	INTEGER	wfn, dRN
	REAL*8  p(0:9)

	off_mKP_fit = 0.D0
	IF (wfn.LT.1 .OR. wfn.GT.4) RETURN
	IF (dRN.LT.0 .OR. dRN.GT.10) RETURN
! .......................................................................
	IF (wfn.EQ.1) THEN		! AV18
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02628D0
	    p(1) = 1.20029D0
	    p(2) = 7.49503D0
	    p(3) = 2.01901D0
	    p(4) = 0.00789D0
	    p(5) = 0.46739D0
	    p(6) = 0.73242D0
	    p(7) = 0.00328D0
	    p(8) = 0.87228D0
	    p(9) = 0.06400D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = 0.03638D0
	    p(1) = 0.38307D0
	    p(2) = 8.01156D0
	    p(3) = 2.30992D0
	    p(4) = 0.09027D0
	    p(5) = 0.69521D0
	    p(6) = 0.75973D0
	    p(7) = -0.05098D0
	    p(8) = 1.18963D0
	    p(9) = -0.19192D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02260D0
	    p(1) = 1.45377D0
	    p(2) = 0.50628D0
	    p(3) = 13.92200D0
	    p(4) = 0.03558D0
	    p(5) = 0.75147D0
	    p(6) = 0.86335D0
	    p(7) = -0.01383D0
	    p(8) = 1.04749D0
	    p(9) = 0.42099D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06410D0
	    p(1) = 1.18883D0
	    p(2) = 6.96799D0
	    p(3) = 8.87113D0
	    p(4) = 0.02603D0
	    p(5) = 0.70504D0
	    p(6) = 1.44139D0
	    p(7) = 0.00004D0
	    p(8) = -1.14305D0
	    p(9) = 0.73785D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06237D0
	    p(1) = 2.03192D0
	    p(2) = 4.01755D0
	    p(3) = 6.83741D0
	    p(4) = 0.04701D0
	    p(5) = -0.00457D0
	    p(6) = 1.30967D0
	    p(7) = -0.00996D0
	    p(8) = -0.42418D0
	    p(9) = 0.27524D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.06759D0
	    p(1) = 1.95103D0
	    p(2) = 3.54215D0
	    p(3) = 11.77533D0
	    p(4) = 0.09269D0
	    p(5) = 0.56534D0
	    p(6) = 0.98398D0
	    p(7) = -0.03031D0
	    p(8) = 3.26913D0
	    p(9) = -0.45923D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.07007D0
	    p(1) = 2.30938D0
	    p(2) = 4.94226D0
	    p(3) = 8.95701D0
	    p(4) = 0.06933D0
	    p(5) = 0.07145D0
	    p(6) = 1.94887D0
	    p(7) = -0.01210D0
	    p(8) = 5.92311D0
	    p(9) = 0.14312D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11965D0
	    p(1) = 2.06149D0
	    p(2) = 5.38881D0
	    p(3) = 12.08265D0
	    p(4) = 0.19668D0
	    p(5) = 0.61820D0
	    p(6) = 0.80489D0
	    p(7) = -0.08735D0
	    p(8) = 3.74802D0
	    p(9) = -0.70773D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.14735D0
	    p(1) = 2.27109D0
	    p(2) = 8.23092D0
	    p(3) = 7.31581D0
	    p(4) = 0.11953D0
	    p(5) = 0.67459D0
	    p(6) = 1.59118D0
	    p(7) = -0.02700D0
	    p(8) = 4.52840D0
	    p(9) = -1.77765D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.27194D0
	    p(1) = 2.01340D0
	    p(2) = 10.71380D0
	    p(3) = 8.84886D0
	    p(4) = 0.09345D0
	    p(5) = 0.49802D0
	    p(6) = 1.28523D0
	    p(7) = -0.00474D0
	    p(8) = 0.58703D0
	    p(9) = 0.88354D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.69848D0
	    p(1) = 1.48173D0
	    p(2) = 17.44991D0
	    p(3) = 12.73730D0
	    p(4) = 0.13118D0
	    p(5) = 0.34598D0
	    p(6) = 1.65884D0
	    p(7) = -0.02215D0
	    p(8) = 1.21306D0
	    p(9) = 0.96399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.2) THEN		! CD-Bonn
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = -0.02820D0
	    p(1) = 0.85879D0
	    p(2) = 9.48856D0
	    p(3) = 2.18885D0
	    p(4) = 0.00070D0
	    p(5) = -5.61817D0
	    p(6) = 14.80512D0
	    p(7) = 0.00348D0
	    p(8) = -1.30292D0
	    p(9) = -0.73075D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.02996D0
	    p(1) = 0.35717D0
	    p(2) = 6.53843D0
	    p(3) = 3.88389D0
	    p(4) = 0.00758D0
	    p(5) = -16.50399D0
	    p(6) = 77.60083D0
	    p(7) = 0.00320D0
	    p(8) = 0.42334D0
	    p(9) = 0.28545D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03261D0
	    p(1) = 0.91185D0
	    p(2) = 8.49348D0
	    p(3) = 10.19681D0
	    p(4) = 0.01598D0
	    p(5) = 0.83748D0
	    p(6) = 1.55960D0
	    p(7) = 0.00085D0
	    p(8) = -0.63447D0
	    p(9) = 0.65632D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.03034D0
	    p(1) = 1.58677D0
	    p(2) = 3.21753D0
	    p(3) = 11.66572D0
	    p(4) = 0.04999D0
	    p(5) = 0.56688D0
	    p(6) = 0.94941D0
	    p(7) = -0.01453D0
	    p(8) = -0.89157D0
	    p(9) = 0.27160D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.04831D0
	    p(1) = 1.75241D0
	    p(2) = 4.74662D0
	    p(3) = 8.29052D0
	    p(4) = 0.04730D0
	    p(5) = 0.33550D0
	    p(6) = 1.18790D0
	    p(7) = -0.00678D0
	    p(8) = -0.42800D0
	    p(9) = 0.36573D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.09019D0
	    p(1) = 1.22091D0
	    p(2) = 1.30114D0
	    p(3) = 17.58701D0
	    p(4) = 0.08312D0
	    p(5) = 0.66902D0
	    p(6) = 0.60767D0
	    p(7) = -0.02035D0
	    p(8) = 0.95978D0
	    p(9) = 1.11322D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.17060D0
	    p(1) = 1.42115D0
	    p(2) = 7.24672D0
	    p(3) = 5.80680D0
	    p(4) = 0.09200D0
	    p(5) = 0.43367D0
	    p(6) = 1.56378D0
	    p(7) = -0.02338D0
	    p(8) = 0.44968D0
	    p(9) = 0.29678D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.11026D0
	    p(1) = 1.85213D0
	    p(2) = 6.74413D0
	    p(3) = 7.74362D0
	    p(4) = 0.08467D0
	    p(5) = 0.24708D0
	    p(6) = 1.12274D0
	    p(7) = -0.01505D0
	    p(8) = 0.44209D0
	    p(9) = 0.36126D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.15291D0
	    p(1) = 1.83333D0
	    p(2) = 7.76495D0
	    p(3) = 7.04783D0
	    p(4) = 0.09206D0
	    p(5) = 0.08655D0
	    p(6) = 1.27460D0
	    p(7) = -0.01659D0
	    p(8) = 0.45536D0
	    p(9) = 0.29407D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.24143D0
	    p(1) = 1.50401D0
	    p(2) = 9.33393D0
	    p(3) = 11.62779D0
	    p(4) = 0.09454D0
	    p(5) = 0.36361D0
	    p(6) = 0.82058D0
	    p(7) = -0.00802D0
	    p(8) = 0.34851D0
	    p(9) = 0.50844D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.22196D0
	    p(1) = 1.87228D0
	    p(2) = 10.18898D0
	    p(3) = 9.21038D0
	    p(4) = 0.11850D0
	    p(5) = 0.34360D0
	    p(6) = 1.28278D0
	    p(7) = -0.01754D0
	    p(8) = 0.54540D0
	    p(9) = 0.53457D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.3) THEN		! WJC-1
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.02322D0
	    p(1) = 0.11213D0
	    p(2) = 3.71079D0
	    p(3) = 5.51496D0
	    p(4) = 0.00877D0
	    p(5) = 0.84639D0
	    p(6) = 0.66227D0
	    p(7) = -0.00621D0
	    p(8) = -0.39896D0
	    p(9) = 0.32012D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.00058D0
	    p(1) = 2.33827D0
	    p(2) = 2.35664D0
	    p(3) = 36.75823D0
	    p(4) = -0.00752D0
	    p(5) = 0.05286D0
	    p(6) = 1.27262D0
	    p(7) = 0.01269D0
	    p(8) = 1.72720D0
	    p(9) = 0.20652D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.03373D0
	    p(1) = 0.93858D0
	    p(2) = 0.15704D0
	    p(3) = 10.71630D0
	    p(4) = -0.00235D0
	    p(5) = -0.11937D0
	    p(6) = 0.74925D0
	    p(7) = 0.00452D0
	    p(8) = 2.96830D0
	    p(9) = -2.89070D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.08982D0
	    p(1) = 0.73060D0
	    p(2) = -0.16543D0
	    p(3) = 12.37035D0
	    p(4) = 0.04407D0
	    p(5) = 0.47361D0
	    p(6) = 0.74570D0
	    p(7) = -0.00933D0
	    p(8) = 0.53186D0
	    p(9) = 0.26943D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.11990D0
	    p(1) = 1.19824D0
	    p(2) = 3.06386D0
	    p(3) = 8.55017D0
	    p(4) = 0.05815D0
	    p(5) = 0.06123D0
	    p(6) = 1.45024D0
	    p(7) = -0.01414D0
	    p(8) = 0.48172D0
	    p(9) = 0.25171D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.15292D0
	    p(1) = 1.01991D0
	    p(2) = 1.20661D0
	    p(3) = 13.31860D0
	    p(4) = 0.02571D0
	    p(5) = -1.56438D0
	    p(6) = 2.69042D0
	    p(7) = -0.00000D0
	    p(8) = 0.29759D0
	    p(9) = 0.97967D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.35935D0
	    p(1) = 0.44637D0
	    p(2) = -0.25510D0
	    p(3) = 16.70057D0
	    p(4) = 0.10634D0
	    p(5) = 0.61659D0
	    p(6) = 0.58524D0
	    p(7) = -0.03335D0
	    p(8) = 0.93904D0
	    p(9) = 0.89819D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.97384D0
	    p(1) = -0.24934D0
	    p(2) = -0.61349D0
	    p(3) = 18.43254D0
	    p(4) = 0.18772D0
	    p(5) = 0.49599D0
	    p(6) = 0.61366D0
	    p(7) = -0.08116D0
	    p(8) = 0.87175D0
	    p(9) = 0.24026D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.21641D0
	    p(1) = 1.74710D0
	    p(2) = 5.19387D0
	    p(3) = 10.61285D0
	    p(4) = 0.06655D0
	    p(5) = 0.01300D0
	    p(6) = 0.94503D0
	    p(7) = -0.00642D0
	    p(8) = 0.48859D0
	    p(9) = 0.16331D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.32283D0
	    p(1) = 1.71708D0
	    p(2) = 7.51556D0
	    p(3) = 9.68202D0
	    p(4) = 0.09871D0
	    p(5) = 0.18788D0
	    p(6) = 0.80490D0
	    p(7) = -0.01673D0
	    p(8) = 0.48879D0
	    p(9) = 0.21016D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.32064D0
	    p(1) = 2.07514D0
	    p(2) = 9.34847D0
	    p(3) = 8.17225D0
	    p(4) = 0.10772D0
	    p(5) = 0.50272D0
	    p(6) = 1.30663D0
	    p(7) = -0.01215D0
	    p(8) = 0.59432D0
	    p(9) = 0.65399D0
	  ENDIF
! .......................................................................
	ELSE IF (wfn.EQ.4) THEN		! WJC-2
	  IF (dRN.EQ.0) THEN		! 0.0%
	    p(0) = 0.03490D0
	    p(1) = 0.78902D0
	    p(2) = -0.25256D0
	    p(3) = 7.98679D0
	    p(4) = 0.00913D0
	    p(5) = 0.74835D0
	    p(6) = 0.60145D0
	    p(7) = -0.00464D0
	    p(8) = 0.41358D0
	    p(9) = 0.22524D0
	  ELSE IF (dRN.EQ.1) THEN	! 0.3%
	    p(0) = -0.01119D0
	    p(1) = 0.50514D0
	    p(2) = 19.35710D0
	    p(3) = 3.32395D0
	    p(4) = 0.00670D0
	    p(5) = 1.38279D0
	    p(6) = 1.24216D0
	    p(7) = 0.00049D0
	    p(8) = 0.38623D0
	    p(9) = 0.23497D0
	  ELSE IF (dRN.EQ.2) THEN	! 0.6%
	    p(0) = 0.02653D0
	    p(1) = 1.27315D0
	    p(2) = -0.53410D0
	    p(3) = 14.08029D0
	    p(4) = 0.01474D0
	    p(5) = 1.82129D0
	    p(6) = 1.99455D0
	    p(7) = -0.00090D0
	    p(8) = 3.96583D0
	    p(9) = 4.61316D0
	  ELSE IF (dRN.EQ.3) THEN	! 0.9%
	    p(0) = 0.06301D0
	    p(1) = 1.10373D0
	    p(2) = -0.26356D0
	    p(3) = 15.04038D0
	    p(4) = 0.02428D0
	    p(5) = -0.15349D0
	    p(6) = 3.03168D0
	    p(7) = 0.00127D0
	    p(8) = -0.73818D0
	    p(9) = 0.07474D0
	  ELSE IF (dRN.EQ.4) THEN	! 1.2%
	    p(0) = 0.06150D0
	    p(1) = 2.15792D0
	    p(2) = 2.18241D0
	    p(3) = 9.84713D0
	    p(4) = 0.03608D0
	    p(5) = -0.13604D0
	    p(6) = 1.12241D0
	    p(7) = -0.00695D0
	    p(8) = -0.35646D0
	    p(9) = 0.31793D0
	  ELSE IF (dRN.EQ.5) THEN	! 1.5%
	    p(0) = 0.07179D0
	    p(1) = 1.97917D0
	    p(2) = 3.47662D0
	    p(3) = 10.00224D0
	    p(4) = 0.04587D0
	    p(5) = 0.06416D0
	    p(6) = 1.10677D0
	    p(7) = -0.00391D0
	    p(8) = -0.42677D0
	    p(9) = 0.26619D0
	  ELSE IF (dRN.EQ.6) THEN	! 1.8%
	    p(0) = 0.09883D0
	    p(1) = 1.96788D0
	    p(2) = 5.19182D0
	    p(3) = 8.82173D0
	    p(4) = 0.06468D0
	    p(5) = 0.11297D0
	    p(6) = 1.63850D0
	    p(7) = -0.00872D0
	    p(8) = 0.52753D0
	    p(9) = 0.41794D0
	  ELSE IF (dRN.EQ.7) THEN	! 2.1%
	    p(0) = 0.14258D0
	    p(1) = 2.00822D0
	    p(2) = 6.23508D0
	    p(3) = 7.81846D0
	    p(4) = 0.07064D0
	    p(5) = -0.05869D0
	    p(6) = 1.24848D0
	    p(7) = -0.01160D0
	    p(8) = 0.48932D0
	    p(9) = 0.22001D0
	  ELSE IF (dRN.EQ.8) THEN	! 2.4%
	    p(0) = 0.16184D0
	    p(1) = 2.16963D0
	    p(2) = 7.62378D0
	    p(3) = 7.33369D0
	    p(4) = 0.09197D0
	    p(5) = 0.15692D0
	    p(6) = 1.80734D0
	    p(7) = -0.01561D0
	    p(8) = 0.53224D0
	    p(9) = 0.39357D0
	  ELSE IF (dRN.EQ.9) THEN	! 2.7%
	    p(0) = 0.20205D0
	    p(1) = 2.28733D0
	    p(2) = 9.10375D0
	    p(3) = 7.24877D0
	    p(4) = 0.08325D0
	    p(5) = 0.36941D0
	    p(6) = 2.39131D0
	    p(7) = -0.00057D0
	    p(8) = 0.41640D0
	    p(9) = 0.90531D0
	  ELSE IF (dRN.EQ.10) THEN	! 3.0%
	    p(0) = 0.95664D0
	    p(1) = 1.11409D0
	    p(2) = 19.00631D0
	    p(3) = 15.97282D0
	    p(4) = 0.15616D0
	    p(5) = 0.40229D0
	    p(6) = 0.85878D0
	    p(7) = -0.03123D0
	    p(8) = 6.75437D0
	    p(9) = -3.83159D0
	  ENDIF

	ENDIF

	off_mKP_fit = -( p(0) * x**p(3) * DEXP(p(1) * x**p(2))
     &                + p(4) * x*DEXP(((x-p(5))/p(6))**2)
     &                + x**0.5D0 * p(7) * DEXP(((x-p(9))/p(8))**2) )

	RETURN
	END


