C **********************************************************************
	FUNCTION PHI_INT3D (yD,gamma,xN,ismear,ioff)
C
C  Interpolate finite-Q^2 nucleon light-cone momentum distribution
C  function in deuteron, including nucleon off-shell effects,
C  in yD, gamma = Sqrt{1 + 4 xB^2 M^2/Q^2} and xN = xB/y = xD / yD.
C
C    yD       = (dp) Deuteron's near light-cone variable yD = y/2
C    gamma    = (dp) gamma = 4*xB^2*m_N^2/Q^2
C    xN       = (dp) nucleon's invariant x        
C    ismear   = (i)  smearing model:   2=nuclear smearing KMK
C                                      3=nuclear smearing LC
C  Written: W. Melnitchouk, A.Accardi (Mar. 2009)
C **********************************************************************
	IMPLICIT NONE
        INTEGER ix,nx, igam,ngam, iy,ny,ismear,ioff,IU,NextUn
        PARAMETER (nx=10, ngam=13, ny=1000)
	REAL*8  PHI_INT3D,xN,yD,gamma,
     &		X_ARRAY(nx),G_ARRAY(ngam),Y_ARRAY(ny),
     &		F_ARRAY(nx,ngam,ny),RINTERP3D
	INTEGER	i,NA(3)
	REAL*8	ARG(3),A(nx+ngam+ny),fint
	character*200 fil,fln,dummy
	integer len1
	integer ismear_old, ioff_old 
	data ismear_old, ioff_old / 9999, 9999 /
*        LOGICAL read /.FALSE./
*	save na,a,f_array,read
	save na,a,f_array,ismear_old,ioff_old

C...Check if data has already been read
*        IF (read) GOTO 300
        IF ((ismear_old.eq.ismear).and.(ioff_old.eq.ioff)) GOTO 300

C...Read f(yD,gamma,xD) array from file
C...Store values of x, gamma, y & f(y,gamma) in X_ARRAY, G_ARRAY, 
C...Y_ARRAY & F_ARRAY
C...note the order of the variables in the data file (x, gamma, y, f)

*        WRITE (*,*) '*** reading f(y,gamma,x) array ...........'

	call GETENV('cteqx_nucl',fln)
	call trmstr(fln,len1)
	fln = fln(1:len1)//'/'
	len1 = len1+1
	
	IU=NextUn()
	if ((ismear.eq.2).and.(ioff.eq.1)) then	    ! KP off-sh. with smearing
	   fil = fln(1:len1)//'phi.km_paris_KPoff'  ! by Kahn-Melnit.-Kulagin
	   OPEN(IU,FILE=fil,FORM='FORMATTED',STATUS='OLD')
	   DO ix=1,nx
	      DO igam=1,ngam
		 DO iy=1,ny
		    READ (IU,*) X_ARRAY(ix), G_ARRAY(igam)
     &			  , Y_ARRAY(iy),F_ARRAY(ix,igam,iy)
		 ENDDO
	      ENDDO
	   ENDDO
	else
	   write(*,*) 'ERROR (PHI_INT3D): ismear,ioff out of range ='
     &         ,ismear,ioff
	end if
        CLOSE (IU)
	call TRMSTR(fil,len1)
	write(6,*),'Nuclear smaring '''//fil(1:len1)//''' read in' 
*        read = .TRUE.
	ismear_old = ismear
	ioff_old = ioff

C...Dimensions of grid
	NA(1) = nx
	NA(2) = ngam
	NA(3) = ny

C...Construct 1-d vector (length nx+ngam+ny) with coordinates of grid
*	print*, '*1'

	DO i=1,nx+ngam+ny
	  IF (i.LE.nx)		          A(i) = X_ARRAY(i)
	  IF (i.GT.nx .AND. i.LE.nx+ngam) A(i) = G_ARRAY(i-nx)
	  IF (i.GT.nx+ngam)		  A(i) = Y_ARRAY(i-nx-ngam)
	ENDDO

*	print*, '*2'

C...Values at which grid to be interpolated
 300	ARG(1) = xN
	ARG(2) = gamma
	ARG(3) = yD

*	print*, '*3 xN,gamma,yD =',xN,gamma,yD
*	print*, '*3 na =',na

	PHI_INT3D = FINT (3,ARG,NA,A,F_ARRAY)

*	print*, '*4 phi3d=',phi_int3d

        RETURN
        END


C **********************************************************************
          FUNCTION FINT (NARG,ARG,NENT,ENT,TABLE)
C
C  Multidimensional interpolation (not just 3) code from
C    Alex Sibirtsev (March, 2009).
C  Uses CERN library routine FINT, based on Newton method.
C
C  NARG: number of dimensions
C  ARG:  vector (dim NARG) at which function to be evaluated
C  NENT: vector (dim NARG) with length of each dimension
C  ENT:  1-d vector with coordinates of grid 
C  TABLE: values of function on grid
C
C CERN PROGLIB# E104    FINT            .VERSION KERNFOR  4.02  820723
C ORIG. 09/08/65 CHL.
C
C   INTERPOLATION ROUTINE. AUTHOR C. LETERTRE.
C   MODIFIED BY B. SCHORR, 1.07.1982.
C   Converted to double precision: A.Accardi, 26 Mar 2009
C **********************************************************************
	  implicit double precision (a-h,o-z)

	  integer narg,NENT(NARG),INDEX(32)
	  double precision fint,ARG(NARG),WEIGHT(32)
	  double precision ENT(9),TABLE(9)

*	  print*, '*00 nent =',nent

          FINT  =  0.
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  RETURN
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1  
          INDEX(1)  =  1
          WEIGHT(1) =  1. 
*	  print*, '*01  nent =',nent
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
*	     print*, '*001  nent =',nent
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
*	     print*, '*002  nent =',nent
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0.)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
*	     print*, '*003 locc =', locc,loca,locb
*	     print*, '*003 nent =',nent
*	     print*, '*003 x,ent(locc)=',x,ent(locc)
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
*	     print*, '*004'
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
*	     print*, '*005'
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
*	     print*, '*006'
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT  
  22            CONTINUE
             GOTO 90
*	     print*, '*007'
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
*	  print*, '*02'
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             FINT  =  FINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
*	  print*, '*03'
          RETURN
          END
