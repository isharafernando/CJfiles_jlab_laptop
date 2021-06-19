C **********************************************************************
	SUBROUTINE DEUTERONSF (ibeam,istruc,inuke,itm,iht,ndrv,v,xc,Fd)
C
C  Compute deuteron F2 structure function (F2d) per nucleon at x (=v(1)) and
C  Q^2 (=v(2)) by smearing nucleon structure function (array F2Narr
C  dimension nx defined over x array xarr) with nucleon light-cone
C  momentum distribution.
C
C    ibeam     = (i) 1=e 2=nu 3=nubar
C    istruc    = (i) 1=F1  2=F2  3=F3
C    inuke     = (i) nuclear corrections: 0=isospin average
C                                         1=density model
C                                         3=nuclear smearing
C
C  this version precomputes F2(xN) in an array, which is interpolated 
C  inside the convolution integral.
C
C  Prepared for CTEQx study (Dec. 2008).
C **********************************************************************
	IMPLICIT NONE
	INTEGER	ibeam,inuke,istruc,itm,iht,ndrv
	double precision v(4),xc(100),Fd
	INTEGER nx,ix
	PARAMETER (nx=100)
	REAL*8	x,Q2,varr(4),xarr(nx),FNarr(nx),dsf,dmc

*      *** Functions
	double precision gomez

*      *** common blocks

*       Proton, neutron fraction    
	double precision az,an
	common/target/az,an

	az=0.5D0		! p + n isoscalar
	an=0.5D0
        
	if (inuke.le.1) then
*         *** p+n isoscalar 
	   CALL strfn (ibeam,istruc,itm,iht,ndrv,v,xc,Fd)
	   if (inuke.eq.1) then
*            *** Density model corrections (for backward compatibility)
	      dmc=gomez(1d0,v(1))
	      Fd = Fd/dmc
	   end if
	else if (inuke.eq.2) then
*         *** Nuclear smearing model
C          ...Compute nucleon str.fn.to be (interpolated and) smeared
	   varr(2) = v(2)	! store Q^2 value in new array
	   DO ix=1,nx
	      x = DFLOAT(ix)/100.D0 - 1d-4
	      xarr(ix) = x	! x array for interpolation 
	      varr(1) = x	! value of x at which SF to be calculated here
	      CALL strfn (ibeam,istruc,itm,iht,ndrv,varr,xc,FNarr(ix))
	   ENDDO
*	   do ix = 95,100
*           print *, '* ix,x,FN=',ix,xarr(ix),FNarr(ix)
*	   end do
C         ...Compute deuteron str.fn. by smearing nucleon SF (F2 only for now)
	   CALL SMEARF2 (v,xarr,FNarr,Fd)
*            print *, '*xB, F2d=',v(1),fd
	else
	   print*, 'ERROR(DEUTERONSF): inuke out of range =',inuke
	end if
	
	RETURN
	END


C **********************************************************************
	SUBROUTINE SMEARF2 (v,xarr,FNarr,F2d)
C
C  Smear nucleon structure function (array FNarr dimension nx defined
C  over x array xarr) with nucleon light-cone momentum distribution.
C
C  xarr is free-nucleon x in [0,1]
C  Convolution set up for "deuteron" xD in [0,1], rather than
C    "nucleon" x in [0,2]
C    => need factor 2 conversion xN = 2 xD [change later].
C
C **********************************************************************
	IMPLICIT NONE
	INTEGER ix,nx,iy,ny
	PARAMETER (nx=100)
	REAL*8  x,FNarr(nx),xarr(nx),F2d, xD
	REAL*8	v(4)
	REAL*8  PHI_INT2D, yD,yint,ymax,fy, F2N,err,sy
	REAL*8  Q2,pi,hc,mN,MD, gamma,mNGeV

C...Constants
	pi = 4*DATAN(1.D0)
	hc = 197.327D0  ! conversion factor (MeV.fm)
	mN = 938.91897D0
	mNGeV = mN/1.D3
	MD = 2*mN - 2.224575D0  ! including binding energy


C...Value of x at which convolution to be made
	x = v(1)	! check whether it is x or 2*x here!!
	xD = x/2	! deuteron Bjorken variable
	Q2 = v(2)

	gamma = DSQRT(1.D0 + 4*x**2 * mNGeV**2/Q2)

C...Convolution approximation (x used in convolution here is xD in [0,1])
	ny = 1000
	yint = (1.D0-xD)/DFLOAT(ny)
	iy = 0
	DO yD=xD+yint,1.D0-yint,yint	! y = yD in [x,1]

	  fy = PHI_INT2D (yD,gamma)	! interpolate smearing function
					! in y and gamma

	  CALL pinterp (xarr,FNarr,nx,xD/yD,F2N,err,0)
				! interpolates FN at xN=xD/yD
				! use linear interpolation (to start with)
	  
	  sy = fy * F2N

          IF (iy.EQ.ny) THEN
	    F2d = F2d + sy
          ELSE IF (iy/2*2.NE.iy) THEN
	    F2d = F2d + 4*sy
	  ELSE IF (iy/2*2.EQ.iy) THEN
	    F2d = F2d + 2*sy
	  ENDIF

*	  print*, '* xN, F2N=', xD/yD,F2N
*	  print*, '* yD, fy=', yD,fy,gamma

	  iy = iy + 1

	ENDDO
	F2d = (yint/3.D0) * F2d

	RETURN
	END


************************************************
*     Density model correction gomez(f) = N/D
      function gomez(f,x)
      implicit double precision (a-h,o-z)
      data p1,p2,p3,p4,p5,p6,p7/1.0164d0,-.047762d0,-.13354d0,
     2.35303d0,.22719d0,-1.2906d0,5.6075d0/
c
c  Fit to the data for F2D/F2N where N=(p+n)/2
c  as extractred by Gomez et al PRD 49, 4348 (1994)
c  [Original reference by Frankfurt and Strikman]
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
