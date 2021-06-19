      function chi2(th,dt,er)
*     Calculates the chi^2 given
*       th    = (dp) the theoretical value
*       dt    = (dp) the data point
*       er    = (dp) the experimental error

      implicit none
      double precision chi2,th,dt,er

      chi2 = ((dt-th)/er)**2

      return
      end


      function residual(th,dt,er)
*     Calculates the residual given
*       th    = (dp) the theoretical value
*       dt    = (dp) the data point
*       er    = (dp) the experimental error

      implicit none
      double precision residual,th,dt,er

      residual = (dt-th)*th/er**2

      return
      end


      function cdfchi2(data,NX,NNPTS,curpar,cdfxsec)
*     Calculates chi^2 for cdf data with correlations

      implicit none
      integer nx,nnpts
      double precision cdfchi2
     &     ,data(nx,nnpts),CURPAR(*),V(10),cdfxsec(100)

      integer i,j,k,npt

      double precision theory

      double precision covinv(33,33),cdfchisqr
      integer ncdf,ncdfl,ncov
      common/cdf/covinv,cdfchisqr,ncdf,ncdfl,ncov

      do 14 i=1,ncdf
         npt=ncdfl+i-1
         do 15 j=1,nx-2
            v(j)=data(j+2,npt)
 15      end do
         cdfxsec(i)=theory(1,npt,v,curpar)
 14   continue
      cdfchisqr=0.d0
      DO j = 1,ncdf
         DO k = 1,ncdf
            cdfchisqr=cdfchisqr+(cdfxsec(j) - data(1,ncdfl+j-1))
     2           *CovInv(j,k)*(cdfxsec(k) - data(1,ncdfl+k-1))
         ENDDO
*         print*,'* ncdf,ncdfl,nx,cdfxsec,data=',ncdf,ncdfl,nx
*     &        ,cdfxsec(j),data(1,ncdfl+j-1)
      ENDDO

      cdfchi2=cdfchisqr

      return
      end


      function chi2wt(itype,v)
*     Determines weights for chi^2 calculation
*     Can be changed by the user at compilation time to suit any 
*     weighting needs. 
*
*     INPUT:
*
*       itype = (dp) the process type 
*                    (0,99=DIS 100-199=DY, 200-299=jets, 300-399=gamma+jet)
*       v(4)  = (dp) the kinematics of the measurement
*                   (process dependent; e.g., 
*                    DIS: v(1)=xB v(2)=Q2, v(3)=W^2 v(4)=process type)
*     OUTPUT:
*
*       chi2wt = (dp) the weight to be applied to the standard chi^2

      implicit none
      double precision chi2wt,v(4),xB
      integer itype

      chi2wt = 1d0
     
*     *** DIS weigths
       if (itype.lt.100) then
         xB = v(1)
         if (xB.gt.0.7) chi2wt = 1d0
      end if

*     *** W asymmetry weigths
      if (itype.eq.131) then
*       ... CDF data
         chi2wt = 1d0
      end if

      return 
      end
