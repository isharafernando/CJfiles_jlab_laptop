**************************************************************
*
* II  Polynomial function interpolation of given order
      subroutine pinterp(xa,ya,n,x,y,dy,order)
*     programmer: Alberto Accardi
*     date: 2/05/01
*
*  A. COMMENTARY
*
*     Performs an interpolation using a polynomial function
*     interpolation at a given order: given an x, it uses "order" points 
*     to its left and "order" to its right to perform the interpolation
*
*     xa(*) = (DP) array with tabulated abscissae (any dimension)
*     ya(*) = (DP) array with tabulated function  (any dimension)
*     n     = (I)  number of tabulated points  
*     x     = (DP) abscissa at which compute the interpolated function
*     y     = (DP) value of the function at x
*     dy    = (DP) estimated error (usually larger than real error)
*     order = (I)  order of the interpolation (see intro)  
*                  If order = 0 performs a linear interpolation
*                  between the nearest neighbours lattice point
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      integer n, order

      double precision xa(*), ya(*), x, y, dy, tempx(n)
     :     , x1(2*order), y1(2*order), xmax, xmin, ymax, ymin  

      integer i, nlow, nmin

*
*  C. ACTION
*

      do i = 1, n
         tempx(i) = xa(i)
      end do
      call hunt(tempx,n,x,nlow)

      if (order.ge.1) then
         if (nlow.lt.order) then
            nmin = 0
         else if (nlow.le.n-order) then
            nmin = nlow-order
         else
            nmin = n-2*order
         end if
         do i = 1, 2*order
            x1(i) = xa(nmin+i) 
            y1(i) = ya(nmin+i) 
         end do
         call polintnum(x1,y1,2*order,x,y,dy)
      else
         ymax = ya(nlow+1)
         ymin = ya(nlow)
         xmax = xa(nlow+1)
         xmin = xa(nlow)
         y = ymin + (ymax-ymin)/(xmax-xmin) * (x-xmin)
      end if

      return
      end


************************************************************************
*
* III search in an ordered table 
      SUBROUTINE hunt(xx,n,x,jlo)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given an array xx(1:n) and given a value x, returns a value j
*     suchthat x is between xx(j) and xx(j+1). xx(1:n) must be monotonic,
*     either decreasing or increasing. j=0 or j=n is returned to
*     indicate that x is out of range.
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER jlo,n

      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

*
*  C. ACTION
*

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo 
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END


************************************************************************
*
* IV  Polynomial interpolation and extrapolation
      SUBROUTINE polintnum(xa,ya,n,x,y,dy)
*     from "NUMERICAL RECIPES IN F77"
*
*  A. COMMENTARY
*
*     Given arrays xa and ya of length n, and given a value x, this
*     routine returns a value y and an error estimate dy. If P(x) is the
*     polynomial of degree N-1 such that P(xa_i) = ya_i, i=1,...,n
*     then the returned value y = P(x).
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables
      
      INTEGER n,NMAX
      double precision dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

*
*  C. ACTION
*

      if (n.gt.nmax) then
         print*, 'ERROR(polintnum): order larger than max', n,'>', nmax 
         stop
      end if
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             print*, 'failure in polintnum'
             stop
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      
