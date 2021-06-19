C **********************************************************************
      FUNCTION PHI_INT2D (y,gamma,ifunc,ismear,iwfn)
C
C  Interpolate finite-Q^2 nucleon light-cone momentum distribution
C  function in deuteron, in y and gamma = Sqrt{1 + 4 x^2 M^2/Q^2}.
C
C  Ref: Kahn, Melnitchouk, Kulagin, arXiv:0809.4308 [nucl-th].
C
C    y        = (dp) near light-cone variable
C    gamma    = (dp) gamma = 4*xB^2*m_N^2/Q^2
C    ifunc    = (i)  which smearing function to interpolate
C                    (0,1,2) = smearing functions
C                    (10,11,12) = (p^2/m^2-1) weighted
C    ismear   = (i)  smearing model:   2=nuclear smearing KMK
C                                      3=nuclear smearing LC
C    iwfn     = (i)  wave-function to be used
C
C  Originally was using IMSL routine
C  > use imsl
C  > f77 -o objectfile file.f
C        -R/site/vni/lib/lib.solaris
C        -L/site/vni/lib/lib.solaris -limsl -lsocket -lnsl
C  Now it uses RINTERP2D
C
C  Written: W. Melnitchouk, CSSM, Oct. 24, 2008
C  New interpolating routine RINTERP2D implemented, Dec. 10, 2008
C  Adapted: A.Accardi, Mar'09
C  New smearing functions for F1, FL, F3 added: Jul'11
C **********************************************************************
      IMPLICIT NONE
      INTEGER ifunc,ismear,iwfn
      integer igam,ngam, iy,ny,IFAIL,len1,i,IU,NextUn,nuklen
      
      PARAMETER (ngam=13, ny=1001)
      character*200 fil,fln,dummy
      character nukfil*150
      double precision  PHI_INT2D,y,gamma
     &     ,G_ARRAY(ngam),Y_ARRAY(ny)
     &     ,F0_ARRAY(ngam,ny),F1_ARRAY(ngam,ny),F2_ARRAY(ngam,ny)
     &     ,G0_ARRAY(ngam,ny),G1_ARRAY(ngam,ny),G2_ARRAY(ngam,ny)
     &     ,DQD2VL,RINTERP2D
      integer ismear_old,iwfn_old 
      data ismear_old,iwfn_old / 9999,9999 /
      save ismear_old,iwfn_old,g_array,y_array
     &     ,f0_array,g0_array,f1_array,g1_array,f2_array,g2_array
      
C...Check if data has already been read
      IF ((ismear.eq.ismear_old).and.(iwfn.eq.iwfn_old)) GOTO 300
      
C...Read f(y,gamma) array from file
C...Store values of gamma, y & f(y,gamma) in G_ARRAY, Y_ARRAY & F_ARRAY
C...note the order of the variables in the data file (gamma, y, f)
      
      call GETENV('cteqx',fln)
c      write(6,*),'data directory path = ',fln
      call trmstr(fln,len1)
      fln = fln(1:len1)//'/nucl/'
      len1 = len1+6
      
      IU=NextUn()
      if (ismear.ge.2) then	! nuclear smearing
         call smearfile(nukfil,nuklen,ismear,iwfn,0)
         fil = fln(1:len1)//nukfil(1:nuklen)
         OPEN(IU,FILE=fil,FORM='FORMATTED',STATUS='OLD')
         do i=1,5
            read(IU,*) dummy
         end do
         DO igam=1,ngam
            DO iy=1,ny
               READ (IU,*) G_ARRAY(igam), Y_ARRAY(iy)
     &              ,F0_ARRAY(igam,iy),G0_ARRAY(igam,iy)
     &              ,F1_ARRAY(igam,iy),G1_ARRAY(igam,iy)
     &              ,F2_ARRAY(igam,iy),G2_ARRAY(igam,iy)
            ENDDO
         ENDDO
      else
         write(*,*) 'ERROR (PHI_INT2D): ismear out of range =',ismear
      end if
      call TRMSTR(fil,len1)
      write(6,*),'Nuclear smearing '''//fil(1:len1)//''' read in' 
      CLOSE (IU)
      ismear_old = ismear
      iwfn_old = iwfn

C...Interpolate at required y value using IMSL routine DQDVAL
C...Inputs are:  DQD2VL (X, Y, NXDATA, XDATA, NYDATA, YDATA, FDATA,
C...                   LDF, CHECK)
c 300    PHI_INT2D = DQD2VL (gamma,y, ngam,G_ARRAY, ny,Y_ARRAY, F_ARRAY,
c     &                      ngam,.FALSE.)

 300  continue

      if (ifunc.eq.0) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,F0_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else if (ifunc.eq.1) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,F1_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else if (ifunc.eq.2) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,F2_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else if (ifunc.eq.10) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,G0_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else if (ifunc.eq.11) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,G1_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else if (ifunc.eq.12) then 
         PHI_INT2D = RINTERP2D (G_ARRAY,Y_ARRAY,G2_ARRAY,gamma,y,
     &                         ngam,ny,0.1D0,0.001D0)
      else 
         write(*,*) 'ERROR(phi_int2d): ifunc out of ramge =',ifunc
         stop
      end if
      
      RETURN
      END


C **********************************************************************
C     2-Dimensional Linear Interpolation Routine:
C     (modified from version sent by Peter Monaghan, Dec. 2008)
C
C       Calculates F(X0,Y0) given array F(X,Y) by two successive
C       interpolations, first in X and then in Y.
C
C       Assumes uniform spacing of array in X and Y.
C
C	(NX,NY) = number of points in (X,Y) arrays.
C	(DELX,DELY) = size of (evenly spaced) intervals in (X,Y) arrays.
C
C **********************************************************************
      DOUBLE PRECISION FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY,DELX,DELY)
      IMPLICIT NONE
      INTEGER	NX,NY
      INTEGER	I_EXTRAP,I,J
      REAL*8	X(NX),Y(NY),F(NX,NY)
      REAL*8	X0,Y0,DELX,DELY
      REAL*8	RINTX1,RINTX2
      
      I_EXTRAP = 0
      I = DINT((X0-X(1))/DELX) + 1
      J = DINT((Y0-Y(1))/DELY) + 1
      IF((I+1.GT.NX).OR.(I.LT.1).OR.(J+1.GT.NY).OR.(J.LT.1))THEN
         I_EXTRAP = I_EXTRAP + 1
         RINTERP2D = 0.D0
         RETURN
      ENDIF

      RINTX1 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J)-F(I,J))
     &     + F(I,J)
      RINTX2 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J+1)-F(I,J+1))
     &     + F(I,J+1)
      RINTERP2D = ((Y0-Y(J))/(Y(J+1)-Y(J)))*(RINTX2-RINTX1) + RINTX1
      
      RETURN
      END
