C ***********************************************************************
	PROGRAM smeargrids
C
C  Calculate smearing functions in WBA and WBAREL formalisms
C
C  adapted from W.Melnitchouk's "PHI_norm.f", May 13, 2010
C
C  Programmer: A.Accardi
C  Date: 14 March 2011
C
C  Compilation:
C  f77 smeargrids.f smearfns.f av18.f wjc.f paris.f cdbonn.f ../../util/pinterp.f
C ***********************************************************************
	IMPLICIT NONE
        INTEGER iy,ny,ig,iwfn,ibar,isimp,i
        PARAMETER (ny=1000)
        REAL*8	y,yint,yplot(0:ny),fy0(3,0:ny),fy1(3,0:ny)
     &      ,norm_overall,norm0,norm1
	REAL*8	PHI_INT2D,PHI_MST,PHI_WBA,PHI_LC
     &      ,gamma,gammamin,gammamax,gammadel
	REAL*8	mN,MD,epsD,phi0(4),phi1(4)
	integer rel,igamma,ngamma
	character model*7,pot*7,cdum*100,outfil*30
	integer*4 today(3), now(3)

C...Constants
        !mN = 938.91897D0        ! nucleon mass
        !epsD = -2.224575D0      ! deuteron binding energy
        !MD = 2*mN + epsD        ! deuteron mass


C...Open data file
	open (12,file='smeargrids.dat',form='formatted',status='old') 
	read (12,*), cdum

c...Cycles until there are valid lines in the input file
c...(and exits becuase of intended 'end of file' error)
 10	read(12,*) iwfn,ibar,rel,outfil

c...    Attention: overwrites previous files
	OPEN (11,FILE=outfil,FORM='FORMATTED',STATUS='unknown')

	if (rel.eq.0) then
	   if (ibar.eq.0) then 
	      model = 'WBA    '
	   else
	      model = 'WBAc   '
	   end if
	else if (rel.eq.1) then
	   if (ibar.eq.0) then 
	      model = 'WBAREL '
	   else
      	      model = 'WBARELc'
	   end if
	else
	   if (ibar.eq.0) then 
	      model = 'AQV    '
	   else
	      model = 'AQVc   '
	   end if
	end if

	if (iwfn.eq.12) then 
	   pot = 'WJC-1  '
	else if (iwfn.eq.13) then 
	   pot = 'WJC-2  '
	else if (iwfn.eq.1) then 
	   pot = 'Paris  '
	else if (iwfn.eq.3) then 
	   pot = 'CD-Bonn'
	else if (iwfn.eq.4) then 
	   pot = 'AV18   '
	end if

	write(11,*) 'Smearing function S(y,gamma) in 1D '//model
     &       //' approximation, finite-Q^2'
	write(11,*) 'Uses '//pot//' potential'
	write(11,*) 'Ref: Kahn et al.arXiv:0809.4308 '
     &       //'+ notes by Melnitchouk and Accardi Apr-May2010'
	call idate(today)	! today(1)=day, (2)=month, (3)=year
	call itime(now)		! now(1)=hour, (2)=minute, (3)=second
	write (11,1000)  today(2), today(1), today(3), now
 1000	format ( ' Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &      i2.2, ':', i2.2, ':', i2.2 )
	write (11,*)
	write (11,*) 'gamma   y       f0         f0_off  '
     &      //'     f1         f1_off       f2         f2_off'
	write (11,*)
	
C...Loop in gamma = sqrt[1 + 4 M^2 x^2/Q^2]
	ig = 1
	norm_overall = 1.D0

	gammamin = 1.0d0
	gammamax = 2.2d0
	gammadel = 0.1d0
	ngamma   = (gammamax-gammamin)/gammadel 
	DO igamma = 0, ngamma
	   gamma = igamma*gammadel + gammamin
	   print *, 'ig,gamma=',ig,sngl(gamma)

C...Integrate over nucleon light-cone momentum fraction
          norm0 = 0.D0
          norm1 = 0.D0
          yint = 1/DFLOAT(ny)
	  DO iy=1,ny-1
	    y = DFLOAT(iy)/DFLOAT(ny) 
	    yplot(iy) = y
	    call smearfn_1D(iwfn,ibar,rel,y,gamma,phi0,phi1)
	    do i = 1,3
	       fy0(i,iy) = phi0(i)
	       fy1(i,iy) = phi1(i)
	    end do
	    isimp = 2 + 2*(iy - iy/2*2)
	    norm0 = norm0 + isimp*fy0(3,iy)
	    norm1 = norm1 + isimp*fy1(3,iy)
	    IF (iy.EQ.iy/10*10) PRINT *, 'gamma,y,f(y)='
     &           ,sngl(gamma),sngl(y),sngl(fy0(3,iy)),sngl(fy1(3,iy))

	  ENDDO
	  yplot(ny) = 1d0
	  do i = 1,3
	     fy0(i,ny) = 0.D0
	     fy1(i,ny) = 0.D0
	  end do
C...Check normalization
          norm0 = (yint/3) * norm0
          norm1 = (yint/3) * norm1
          WRITE (*,*) 'f0(y) normalization =',norm0/norm_overall
     &  ,norm_overall
          WRITE (*,*) 'f1(y) normalization =',norm1/norm_overall
	  IF (ig.EQ.1) norm_overall = norm0
	  write(*,*) 'norm_overall', norm_overall

C...Write calculated distribution function to file
	  DO iy=0,ny
	    WRITE (11,300) gamma, yplot(iy),
     &	       (fy0(i,iy)/norm_overall, fy1(i,iy)/norm_overall,i=1,3)
	  ENDDO

 234	  ig = ig + 1
C...gamma
	ENDDO

	close(11)

C...goes to next line in input file
	goto 10

	CLOSE (12)
 100    FORMAT (1001(F15.12,1X))
 200    FORMAT (2(F15.12,1X))
 201    FORMAT (F15.12)
 300    FORMAT (F3.1,1X,F7.4,1X,6(E11.4,1X))
	STOP
	END
