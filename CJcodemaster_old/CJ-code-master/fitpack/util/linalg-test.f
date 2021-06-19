      program linalg_test

      implicit none

      integer np
      parameter (np=10)

      integer N,INDX(NP),i,j,imax
      double precision A(NP,NP),B(NP),Y(NP,NP),D(NP),E(NP)

*     matrix A (stored by columns, so it looks like A transposed)
      data n/2/
      data a/ 1.,2.,3. ,7*0d0
     &       ,2.,1.,1. ,7*0d0
     &       ,3.,1.,2. ,7*0d0
     &       ,70*0d0/


      print*
      print*, 'A-matrix'
      print*
      do i=1,n
         print 200, (a(i,j),j=1,n)
 200     format(10G12.4)
      end do

      print*
      print*, 'Y-matrix, i.e., A^-1'
      print*
      call invmatrix(A,N,NP,Y)  ! A is preserved, Y=A^-1
      do i=1,n
         print 200, y(1,i),y(2,i), y(3,i)
      end do


      print*
      print*, 'Eigenvalues and vectors of A'
      print*
      call tredn2(A,N,NP,D,E)   ! on exit, A is orthonormal transformation 
                                ! to bring original A to tridiagonal form
      call tqli(D,E,N,NP,A)     ! on exit, D contains eigenvalues, and the 
                                ! columns of A the corresponding eigenvector
      print*, '    e-value       e-vector ----> '
      do i=1,n
         print 100, i,D(i),(A(j,i),j=1,n)
 100     format('#',I1,' ',g12.4,'   ',10g12.4)
      end do
      print*

      stop
      end


      include 'linalg.f'
