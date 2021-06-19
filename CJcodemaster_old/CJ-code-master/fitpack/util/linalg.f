*****************************************************************
*
*     Linear algebra package from Numerical Recipes 
*
*     CONTENTS:
*
*      I. Linear systmes and matri inversion
*     II. Eigensystems (eigenvalues and eigenvectors)
*
*****************************************************************


*****************************************************************
*
* I. Linear systems and matrix inversion
*
*****************************************************************

*****************************************************************
      subroutine invmatrix(A,N,NP,Y)

*     Given an N x N matrix A, with physical dimension NP, this routine 
*     finds its inverse Y (with physical dimensions NP) such that A.Y=1. 
*     The inverse is found column by column by repeated applications of 
*     the LU decomposition method for solving sets of linear equations. 
*     A and N are input, Y is output. 

      implicit none
      
      integer N,NP,INDX(NP),i,j
      double precision A(NP,NP),B(NP,NP),Y(NP,NP),d

      do i=1,n                  ! sets up the identity matrix
         do j=1,n
            y(i,j) = 0d0
            b(i,j) = a(i,j)     ! A is preserved, B used (and destroyed) below
         end do
         y(i,i) = 1d0
      end do
      
      call ludcmp(b,n,np,indx,d) ! decomposes the matrix just once
      
      do j=1,n                  ! finds inverse column by column
         call lubksb(b,n,np,indx,y(1,j))
         ! note that fortran stores 2-dimensional matrices by column, 
         ! so that y(1,j) is the address of the jth column of y.
      end do

      return
      end


*****************************************************************
      subroutine ludcmp(A,N,NP,INDX,D)

*     given a N x N matrix A, with physical dimension NP, this routine 
*     replaces it by the LU decomposition of a rowwise permutation of itself.
*     A and N are input, A is output, arranged as in Eq.(2.3.14); INDX is
*     an output vector which records teh row permutation effected by the 
*     partial pivoting; D is output as +-1 depending on whether the number of 
*     row interchanges was even or odd, respectively. This routine is used in 
*     combination with LUBKSB to solve linear equations or invert a matrix.
*
*     REF: Numerical recipes, Cambridge U. press, 1986

      implicit none
      
      integer nmax
      double precision tiny
      parameter (nmax=100,tiny=1.0d-20)
      
      integer N,NP,INDX(NP),i,j,k,imax
      double precision A(NP,NP),D,VV(NMAX),aamax,sum,dum

      D=1d0

      do 12 i=1,n
         aamax = 0d0
         do 11 J=1,n
            if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
 11      end do
         if (aamax.eq.0d0) then
            print*, 'Error(ludcmp): singular matrix!'
            stop
         end if
         vv(i) = 1d0/aamax
 12   end do

      do 19 j=1,n
         do 14 i=1,j-1
            sum = a(i,j)
            do k=1,i-1
               sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j) = sum
 14      end do
         aamax = 0d0
         do 16 i=j,n
            sum = a(i,j)
            do k=1,j-1
               sum = sum - a(i,k)*a(k,j)
            end do
            a(i,j) = sum
            dum = vv(i)*dabs(sum)
            if (dum.ge.aamax) then
               imax=i
               aamax = dum
            end if
 16      end do
         if (j.ne.imax) then 
            do 17 k=1,n
               dum = a(imax,k)
               a(imax,k) = a(j,k)
               a(j,k) = dum
 17         end do
            d = -d 
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (j.ne.n) then
            if (a(j,j).eq.0d0) A(J,J)=tiny
            dum = 1d0/a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            end do
         end if
 19   end do

      if (A(n,n).eq.0d0) a(n,n)=tiny
      
      return
      end


*****************************************************************
      subroutine lubksb(A,N,NP,INDX,B)

*     Solves the set of linear equations A.X=B. Here A is input, not as the
*     matrix A but rather its LU decomposition determined by calling LUDCMP.
*     INDX is input as the permutation vector returned by LUDCMP. B is input
*     as the right-hand side vector B, and returns with the solution vector X.
*     A, N, NP, INDX are not modified by this routine, and can be left in 
*     place for wuccessive calls with different right-hand side vector B. 
*     This routine takes into account teh possibility that B begin with 
*     many zero elements, so it is efficient for use in matrix inversion.


      implicit none
      
      integer N,NP,INDX(NP),i,j,ii,ll
      double precision A(NP,NP),B(NP),sum

      ii = 0
      do 12 i=1,n
         ll = indx(i)
         sum = b(ll) 
         b(ll)=b(i)
         if (ii.ne.0) then
            do j=ii,i-1
               sum = sum - a(i,j)*b(j)
            end do
         else if (sum.ne.0d0) then
            ii = i
         end if
         b(i) = sum
 12   end do

      do 14 i=n,1,-1
         sum = b(i)
            do j=i+1,n
               sum = sum - a(i,j)*b(j)
            end do
         b(i) = sum / a(i,i)
 14   end do

      return
      end


      SUBROUTINE INVMAT(ARRAY,NORDER,DET) 
      implicit real*8 (a-h,o-z)
      DIMENSION ARRAY(50,50),IK(50),JK(50)
C
C  MATRIX INVERSION ROUTINE FROM BEVINGTON, PAGE 302.
C  It takes 3*Norder times longer than 'invmatrix'
C
      DET=1. 
      DO 100 K=1,NORDER 
      AMAX=0.D0 
   21 DO 30 I=K,NORDER 
      DO 30 J=K,NORDER 
c      IF( ABS(AMAX).GT. ABS(ARRAY(I,J))) GO TO 30 
      IF(DABS(AMAX).GT.DABS(ARRAY(I,J))) GO TO 30 
   24 AMAX=ARRAY(I,J)
      IK(K)=I 
      JK(K)=J 
   30 CONTINUE 
      IF(AMAX.NE.0.D0) GO TO 41 
   32 DET=0.
      GO TO 140 
   41 I=IK(K) 
      IF(I-K) 21,51,43 
   43 DO 50 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=ARRAY(I,J) 
   50 ARRAY(I,J)=-SAVE 
   51 J=JK(K) 
      IF(J-K) 21,61,53 
   53 DO 60 I=1,NORDER
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=ARRAY(I,J) 
   60 ARRAY(I,J)=-SAVE 
   61 DO 70 I=1,NORDER 
      IF(I.EQ.K) GO TO 70 
   63 ARRAY(I,K)=-ARRAY(I,K)/AMAX 
   70 CONTINUE 
      DO 80 I=1,NORDER 
      DO 80 J=1,NORDER 
      IF(I.EQ.K.OR.J.EQ.K) GO TO 80 
   75 ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J) 
   80 CONTINUE 
      DO 90 J=1,NORDER 
      IF(J.EQ.K) GO TO 90 
   83 ARRAY(K,J)=ARRAY(K,J)/AMAX 
   90 CONTINUE 
      ARRAY(K,K)=1./AMAX 
  100 DET=DET*AMAX 
      DO 130 L=1,NORDER 
      K=NORDER-L+1
      J=IK(K) 
      IF(J.LE.K) GO TO 111 
  105 DO 110 I=1,NORDER 
      SAVE=ARRAY(I,K) 
      ARRAY(I,K)=-ARRAY(I,J) 
  110 ARRAY(I,J)=SAVE 
  111 I=JK(K) 
      IF(I.LE.K) GO TO 130 
  113 DO 120 J=1,NORDER 
      SAVE=ARRAY(K,J) 
      ARRAY(K,J)=-ARRAY(I,J) 
  120 ARRAY(I,J)=SAVE 
  130 CONTINUE 
  140 RETURN 
      END 




*****************************************************************
*
* II. Eigensystems
*
*****************************************************************


*****************************************************************
      subroutine tqli(D,E,N,NP,Z)
*     QL algorithm with implicit shifts to determine the eigenvalues 
*     and eigenvectors of a real, symmetric, tridiagonal NxN matrix of 
*     physical dimension NP, possibly previously reduced by TREDN2. D is 
*     a vector of length NP, On input its first N elements are the diagonal 
*     elements of the tridiagonal matrix. On output it returns the eigenvalues.
*      The vector E, of length NP, inputs the subdiagonal elements in E(2:N) 
*     with E(1) arbitrary. On output E is destroyed. When finding only 
*     eigenvalues several lines can be omitted, as indicated in the comments. 
*     If the eigenvectors of the tridiagonal matrix are desired, the Z matrix 
*     (NxN stored in a a NPxNP array) is input as the identity matrix. If the 
*     eigenvectors of a matrix reduced by TREDN2 are required, then Z is input 
*     as the matrix output by TREDN2. In either case, the Kth column of Z 
*     returns the normalized eigenvector corresponding to D(K).
 
      implicit none
      
      integer N,NP,i,j,k,l,m,iter

      double precision  D(NP),E(NP),Z(NP,NP),dd,g,r,s,c,p,b,f

      if (N.gt.1) then
         do i=2,N
            E(i-1) = E(i)       ! convenient to renumber elements of E
         end do
         E(N) = 0d0
         do 15 l=1,N
            iter = 0
 1          do m=l,N-1          ! look for a small subdiagonal
               dd = dabs(D(m))+dabs(D(m+1)) ! element to split the matrix
               if (dabs(E(m))+dd.eq.dd) goto 2
            end do
            m=N
 2          if (m.ne.l) then
               if (iter.eq.30) then
                  print*, 'ERROR(tqli): too many iterations'
                  stop
               end if
               iter = iter + 1
               g = (D(l+1)-D(l))/(2d0*E(L)) ! Form shift
               r = dsqrt(g**2+1d0) 
               g = D(m)-D(l)+E(l)/(g+sign(r,g)) ! This is d_m-k_s
               s = 1d0
               c = 1d0
               p = 0d0
               do 14 i=m-1,l,-1 ! A plane rotation as in the original QL, 
                  f = s*E(i)    ! followed by Givens rotation to restore 
                  b = c*E(i)    ! tridiagonal form
                  if (dabs(f).ge.dabs(g)) then
                     c = g/f
                     r = dsqrt(c**2+1d0)
                     E(i+1) = f*r
                     s = 1d0/r
                     c = c*s
                  else
                     s = f/g
                     r = dsqrt(s**2+1d0)
                     E(i+1) = g*r
                     c = 1d0/r
                     s = s*c
                  end if
                  g = D(i+1)-p
                  r = (D(i)-g)*s+2d0*c*b
                  p = s*r
                  D(i+1) = g+p
                  g = c*r-b
                  ! omit lines from here...
                  do 13 k=1,N
                     f = Z(k,i+1)
                     Z(k,i+1) = s*Z(k,i)+c*f
                     Z(k,i) = c*Z(k,i)-s*f
 13               end do
                  ! ...to here when finding only eigenvalues
 14            end do
               D(l) = D(l)-p
               E(l) = g
               E(m) = 0d0
               goto 1
            end if
 15      end do
      end if

      return
      end


*****************************************************************
      subroutine tredn2(A,N,NP,D,E)
*     Housholder reduction of a real, symmetric, NxN matrix A stored 
*     in a NPxNP physical array. On output, A is replaced by the orthogonal 
*     matrix Q effecting teh transformation. D returns the diagonal elements 
*     of the tridiagonal matrix, and E the off-diagonal elements, with 
*     arbitrary E(1)=0. The lines indicated in comments can be omitted for 
*     speed if only eigenvalues are to be found.

      implicit none

      integer N,NP,i,j,k,l

      double precision A(NP,NP),D(NP),E(NP),f,g,h,hh,scale

      if (N.gt.1) then
         do 18 i=N,2,-1
            l = i-1
            h = 0d0
            scale = 0d0         ! skip transformation
            if (l.gt.1) then
               do k=1,l
                  scale = scale + dabs(a(i,k))
               end do
               if (scale.eq.0d0) then
                  E(i) = A(i,l)
               else
                  do k=1,l      
                     A(i,k) = A(i,k)/scale ! use scales a's for transformation
                     h = h+A(I,K)**2       ! form sigma in H
                  end do
                  f = A(i,l)
                  g = - sign(sqrt(h),f)
                  E(i) = scale*g
                  h = h-f*g
                  A(i,l) = f-g
                  f = 0d0
                  do 15 j=1,l
                     ! omit following line if finding only eigenvalues
                     A(j,i) = A(i,j)/h  !store u/H in i^th column of A
                     g = 0d0            ! Form an element of A.u in G
                     do k=1,j
                        g = g+A(j,k)*A(i,k)
                     end do
                     if (l.gt.j) then
                        do k=j+1,l
                           g = g+A(k,j)*A(i,k)
                        end do
                     end if
                     E(j) = g/h ! Form element of p in temporarily unused E
                     f = f+E(j)*A(i,j)
 15               end do
                  hh = f/(h+h)  ! Form K, Eq.(11.2.11)
                  do 17 j=1,l
                     f = A(i,j)
                     g = E(j)-hh*f
                     E(j) = g
                     do k=1,j
                        A(j,k) = A(j,k)-f*E(k)-g*A(i,k)
                     end do
 17               end do
               end if
            else
               E(i) = A(i,l)
            end if
            D(i) = h
 18      end do
      end if

      ! Omit following line  if finding only eigenvalues
      D(1) = 0d0
      E(1) = 0d0
      do 23 i=1,N               ! begin accumulation of transformation matrices
         ! Delete lines from here...
         l = i-1
         if(D(i).ne.0d0) then   ! This block skipped when i=1
            do 21 j=1,l
               g = 0d0
               do k=1,l         ! Use u and u/H stored in A to form P.Q
                  g = g+A(i,k)*A(k,j)
               end do
               do k=1,l
                  A(K,j) = A(k,j)-g*A(k,i)
               end do
 21         end do
         end if
         ! ...to here when finding only eigenvalues
         D(i) = A(i,i)
         ! Also delete from here...
         A(i,i) = 1d0           ! reset row and column of A to identity 
         if (l.ge.1) then
            do j=1,l
               A(i,j) = 0d0
               A(j,i) = 0d0
            end do
         end if
         ! ...to here when finding only eigenvalues
 23   end do

      return 
      end



            
         
