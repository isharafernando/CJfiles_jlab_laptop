*******************************************************************
* These are old routine to perform simplified nuclear corrections *
* to the deuteron F2 structure function                           *
* Superseded in 2008 by 'deuteronSF.f'                            *
*******************************************************************


      function wally(f,x)
c
c     Fit to the Wally's offshell curve for F2D/F2N where N=(p+n)/2 
c     plotted in 
c     http://hallcweb.jlab.org/elogs/CTEQ+Large+x+Working+Group/81
c
c     INPUT:
c
c     f  = (dp) Deuterium F2 (exp data)
c     x  = (dp) Bjorken's xB
c  
c     OUTPUT
c
c     wally = (dp) Corrected F2(D) = [F2(p)+F2(n)]/2
c
c     Divides deuterium data by the function 'theory' or multiply the 
c     theoretical expression by the function 'theory'
c     the function 'wally' is set up to convert the DATA
c

      implicit double precision (a-h,o-z)



*     *** Smearing function in Weak Binding Approx - Q^2=2 GeV^
      data p1,p2,p3,p4,p5,p6,p7/1.0010d0,4.6625d-02,-.18022d0,
     &     -.21123d0,.41024d0,0.00d0,5.9409d0/

*     *** Smearing function in Weak Binding Approx - Q^2=20 GeV^2
*      data p1,p2,p3,p4,p5,p6,p7/.98322d0,.10631d0,-8.6856d-2,
*     &     -.57174d0,.86329d0,0.00d0,7.4018d0/

*     *** Smearing function in Bj limit
*      data p1,p2,p3,p4,p5,p6,p7/1.0187d0,-.14497d0,.04005d0,
*     &     .14832d0,.05776d0,-.15067d0,7.8420d0/


      THEORY=P1+P2*X+P3*x**2+P4*x**3+P5*x**4
     2+P6*x**5
      deufac=p7*(1.d0/x-1.d0)
      if(deufac.ge.20.d0)deufac=20.d0
      deucor=1.d0-exp(-deufac)
      theory=theory/deucor
      wally=f/theory
      RETURN
      END


      function gomez(f,x)
      implicit double precision (a-h,o-z)
      data p1,p2,p3,p4,p5,p6,p7/1.0164d0,-.047762d0,-.13354d0,
     2.35303d0,.22719d0,-1.2906d0,5.6075d0/
c
c  Fit to the data for F2D/F2N where N=(p+n)/2
c  as extractred by Gomez et al PRD 49, 4348 (1994)
c
      THEORY=P1+P2*X+P3*x**2+P4*x**3+P5*x**4
     2+P6*x**5
      deufac=p7*(1./x-1.)
      if(deufac.ge.20.)deufac=20.
      deucor=1.-exp(-deufac)
      theory=theory/deucor
      gomez=f/theory
      RETURN
      END


      function deucor(f,x,q2)
*
*  An old routine used by Jeff in exploratory fits of deuterium corrections.
*  Can probably by deleted. As of 10/22/08, it is not needed for compilation
*  of 'altpfit02'
*

      implicit real*8 (a-h,o-z)
      dimension par(7)
      data par/.98878,.20046,-.68491,.85263,-.42238,-.94523
     2,4.3473/
      THEORY=PAR(1)+PAR(2)*X+PAR(3)*x**2+PAR(4)*x**3+PAR(5)*x**4
     2+PAR(6)*x**5
      deufac=par(7)*(1./x-1.)
      if(deufac.ge.20.)deufac=20.
      deu=1.-exp(-deufac)
      theory=theory/deu
      deucor=f/theory
      return
      end      
