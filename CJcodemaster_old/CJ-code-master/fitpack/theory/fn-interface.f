******************************************************************
* M. Wobisch - July 26, 2005           fn-interface.f
*
*   fastNLO user interface to PDF and alpha_s code
*   --> to be edited by user 
*       to interface own PDF/alphas code
*
* included:    
*      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
*      SUBROUTINE FNPDF(X,MUF,XPDF)
*
*  in the default version the PDF interface FNPDF is set up
*  to access LHAPDF - the alpha_s routine FNALPHAS calls
*  alphas-demo.f which is an iterative solution of the 2-loop RGE
*******************************************************************


      DOUBLE PRECISION FUNCTION FNALPHAS(MUR)
*-----------------------------------------------------------------
* MW 06/29/2005  - alphas interface to the fastNLO code
*
* alpha_s computation
*   input:   MUR     renormalization scale in GeV
*   output:  value of (alpha_s/2Pi) at scale MUR
*
*  !!!! again: this function must return  alpha_s/(2Pi)  !!!!!
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nloop, neff
      INTEGER ILOOP,IORD,NMAX,IVL,NCB
      DOUBLE PRECISION xmc,xmb
      DOUBLE PRECISION MUR, ALPSMZ, PI
      DOUBLE PRECISION ALPS_IT, alam5, alpha_s, amu0 
c      common/alambda/alam5
      common/jetpdfs2/amu0,alam5
      COMMON/INPUT/ILOOP,IORD,NMAX,IVL,XMC,XMB,NCB
c      DOUBLE PRECISION MY_FAVOURITE-ALPHAS
      PARAMETER (PI=3.1415927d0)

c === example: exact, iterative solution of the 2-loop RGE 
      nloop=2
      alpsmz=0.118              ! set here the value of alpha_s(Mz)
c      alpsmz=0.1185              ! for H1-2000 MSbar
c      alpsmz=0.1205             ! for MRST2004
c      FNALPHAS = ALPS_IT(MUR,ALPSMZ,NLOOP)/2d0/PI
c      alam5=0.226
      FNALPHAS=alpha_s(iord+1,MUR**2,alam5,neff)/2./pi
c - here you can call your own alpha_s code
c           -> remember to divide by 2Pi
c
c     FNALPHAS = MY_FAVOURITE-ALPHAS(MUR)
c

      RETURN
      END

C *****************************************************************

      SUBROUTINE FNPDF(X,MUF,XPDF)
*-----------------------------------------------------------------
* MW 06/29/2005 
*
* PDF interface to the fastNLO usercode
*
*   input   X       parton momentum fraction 
*           MUF     factorization scale in GeV
*
*   output  XPDF(-6:6) array of PDF momentum densities i.e. x*pdf!
*                      using the LHAPDF numbering convention:
*        tbar, bbar, cbar, sbar, ubar, dbar, g, d, u, s, c, b, t
*         -6 ,  -5 ,  -4 ,  -3 ,  -2 ,  -1 , 0, 1, 2, 3, 4, 5, 6
*
*-----------------------------------------------------------------
      IMPLICIT NONE
      integer ndrv,numfl
      DOUBLE PRECISION X, MUF, XPDF(-6:6), ctq6pdf, als
      DOUBLE PRECISION u,d,ub,db,sb,cb,bb,glue,al,al0,amu0,alambda5
      common/jetpdfs/ndrv,numfl,als
      common/jetpdfs2/amu0,alambda5

c ======= example: interface to LHAPDF ============================
c remember to initialize the LHAPDF set first (in the main routine)
c                                       -> see the fastNLO example
c             call InitPDFset("cteq61.LHgrid")
c             call InitPDF(0)
c      call evolvePDF(X,MUF,XPDF)
c


c === here you can call your own PDF code
c  -> remember that these are *momentum densities* i.e. x * PDF
c     and the scale is in GeV
c
c     call MY-FAVORITE-PDFS(....)
c
      al0=2.*log(amu0/alambda5)
      al=2.*log(muf/alambda5)
      als=log(al/al0)
      call fsupdf(ndrv,x,als,u,d,ub,db,sb,cb,bb,glue)
c
c  Note: pdf returns momentum densities - no additional factor of x needed
c
      xpdf(6)=0.d0
      xpdf(5)=bb
      xpdf(4)=cb
      xpdf(3)=sb
      xpdf(2)=u
      xpdf(1)=d
      xpdf(0)=glue
      xpdf(-1)=db
      xpdf(-2)=ub
      xpdf(-3)=xpdf(3)
      xpdf(-4)=xpdf(4)
      xpdf(-5)=xpdf(5)
      xpdf(-6)=xpdf(6)

c      call setctq6(200)
c      xpdf(6)=0.d0
c      xpdf(5)=x*ctq6pdf(5,x,MUF)
c      xpdf(4)=x*ctq6pdf(4,x,MUF)
c      xpdf(3)=x*ctq6pdf(3,x,MUF)
c      xpdf(2)=x*ctq6pdf(1,x,MUF)
c      xpdf(1)=x*ctq6pdf(2,x,MUF)
c      xpdf(0)=x*ctq6pdf(0,x,MUF)
c      xpdf(-1)=x*ctq6pdf(-2,x,MUF)
c      xpdf(-2)=x*ctq6pdf(-1,x,MUF)
c      xpdf(-3)=xpdf(3)
c      xpdf(-4)=xpdf(4)
c      xpdf(-5)=xpdf(5)
c      xpdf(-6)=xpdf(6)



      RETURN
      END
