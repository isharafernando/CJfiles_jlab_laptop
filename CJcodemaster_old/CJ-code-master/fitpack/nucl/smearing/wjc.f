C **********************************************************************
      SUBROUTINE WJC (iwjc,q,u,w,vt,vs)
C
C  Deuteron wavefunction from WJC-1 (Gross et al.) NN potential model.
C  Gross & Stadler, Phys. Rev. C78, 014005 (2008).
C
C  Note: q in 1/fm, wfns in fm^3/2.
C
*       iwjc = 1 --> WJC1
*            = 2 --> WJC2
*
C  Wave function normalization
C    \int dq q^2 (u^2+w^2+vt^2+vs^2) + V' term = 1
C    so that wave functions themselves normalized to ~102% (WJC-1)
C    and ~105% (WJC-2)
C    => renormalize to 1 for structure functions
C
C **********************************************************************
        IMPLICIT NONE
        INTEGER id,nq,iwjc
        PARAMETER (nq=60)
        REAL*8  q,u,w,vt,vs,du,dw,dvt,dvs,
     &		qgrid(nq),ugrid(nq),wgrid(nq),vtgrid(nq),vsgrid(nq),
     &		DQDVAL
        integer iwjc_old/0/
	save iwjc_old,qgrid,ugrid,wgrid,vtgrid,vsgrid
        REAL*8  pi,hcM,hcG

        pi = 4*DATAN(1.D0)
        hcM = 197.327D0		! GeV.fm conversion factor
        hcG = 0.197327D0	! MeV.fm conversion factor

        IF (iwjc.eq.iwjc_old) GO TO 999     ! Data already read
C...Read data from file
	if (iwjc.eq.1) then 
	   OPEN (10, FORM='FORMATTED',
     &          FILE='wjc-1.dat',
     &          STATUS='OLD')
	else if (iwjc.eq.2) then 
	   OPEN (10, FORM='FORMATTED',
     &          FILE='wjc-2.dat',
     &          STATUS='OLD')
	else
	   print*, 'ERROR(wjc): iwjc out of range:', iwjc
	   stop
	end if

C...Momentum space [qgrid in MeV, uqgrid in GeV^-3/2]
        DO id=1,nq
          READ (10,*) qgrid(id),
     &		      ugrid(id), wgrid(id), vtgrid(id), vsgrid(id)
          qgrid(id) = qgrid(id) / hcM		! MeV => 1/fm
          ugrid(id)  = ugrid(id)  * hcG**1.5D0	! GeV^-3/2 => fm^3/2
          wgrid(id)  = wgrid(id)  * hcG**1.5D0
          vtgrid(id) = vtgrid(id) * hcG**1.5D0
          vsgrid(id) = vsgrid(id) * hcG**1.5D0
        ENDDO
	if (iwjc.eq.1) then
	   PRINT *, '... WJC-1 model read...'
	else if (iwjc.eq.2) then
	   PRINT *, '... WJC-2 model read...'
	end if
	iwjc_old = iwjc

C...Evaluate wavefunction
 999	continue
        
	!u  = DQDVAL (q,nq,qgrid,ugrid,.FALSE.)
	!w  = DQDVAL (q,nq,qgrid,wgrid,.FALSE.)
	!vt = DQDVAL (q,nq,qgrid,vtgrid,.FALSE.)
	!vs = DQDVAL (q,nq,qgrid,vsgrid,.FALSE.)

	call pinterp(qgrid,ugrid,nq,q,u,du,2)
	call pinterp(qgrid,wgrid,nq,q,w,dw,2)
	call pinterp(qgrid,vtgrid,nq,q,vt,dvt,2)
	call pinterp(qgrid,vsgrid,nq,q,vs,dvs,2)

        RETURN
        END
