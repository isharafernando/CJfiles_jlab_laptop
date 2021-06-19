C ******************************************************************************
        SUBROUTINE PARIS (mode,q,u0,u2)
C
C  Deuteron wavefunction from Paris NN potential model, parameterised in
C  Lacombe, Loiseau et al. Phys.Lett. 101B (1981) 139.
C
C  mode=1: momentum space (q im 1/fm, u0,u2 in fm**3/2)
C	   normalization \int dq q^2 (u0^2+u2^2) = 1
C
C  mode=2: coordinate space (q=r in fm, u0,u2 in fm**-1/2)
C ******************************************************************************
        IMPLICIT NONE	!UNDEFINED (A-Z)

        INTEGER mode
        REAL*8  q,u0,u2
        INTEGER j,n
        PARAMETER (n=13)
        REAL*8  C(n),D(n),d1,d2,d3,m(n),alpha,m0
        REAL*8  pi,c1
        LOGICAL initcoef /.FALSE./
        save initcoef,alpha,m0,m,C,D

C...Cj,Dj in fm**(-1/2)
        DATA    C /0.88688076D0, -0.34717093D0, -0.30502380D1,
     &             0.56207766D2, -0.74957334D3,  0.53365279D4,
     &            -0.22706863D5,  0.60434469D5, -0.10292058D6,
     &             0.11223357D6, -0.75925226D5,  0.29059715D5,
     &             0.00000000D0/
        DATA    D /0.23135193D-1, -0.85604572D0, 0.56068193D1,
     &            -0.69462922D2,  0.41631118D3, -0.12546621D4,
     &             0.12387830D4,  0.33739172D4, -0.13041151D5,
     &             0.19512524D5,  0.00000000D0,  0.00000000D0,
     &             0.00000000D0/


C...Value of pi
        pi = 4*DATAN(1.D0)


C...Initialise coefficients.....................................................
        IF (.NOT.initcoef) THEN
          initcoef = .TRUE.
C...alpha, m0 in 1/fm
          alpha = 0.23162461D0
          m0 = 1.D0
C...Set masses (units of 1/fm)
          DO j=1,n
            m(j) = alpha + DFLOAT(j-1)*m0
          ENDDO

C...Set last Cj coefficient
          DO j=1,n-1
            C(n) = C(n) - C(j)
          ENDDO

C...Set last 3 Dj coefficients
          d1 = 0.D0
          d2 = 0.D0
          d3 = 0.D0
          DO j=1,n-3
            d1 = d1 + D(j)/m(j)**2
            d2 = d2 + D(j)
            d3 = d3 + D(j)*m(j)**2
          ENDDO
          D(n-2) = m(n-2)**2 / (m(n)**2 - m(n-2)**2)
     &                       / (m(n-1)**2 - m(n-2)**2)
     &           * ( - m(n-1)**2 * m(n)**2 * d1
     &               + (m(n-1)**2 + m(n)**2) * d2
     &               - d3 )
          D(n-1) = m(n-1)**2 / (m(n-2)**2 - m(n-1)**2)
     &                       / (m(n)**2 - m(n-1)**2)
     &           * ( - m(n)**2 * m(n-2)**2 * d1
     &               + (m(n)**2 + m(n-2)**2) * d2
     &               - d3 )
          D(n)   = m(n)**2 / (m(n-1)**2 - m(n)**2)
     &                     / (m(n-2)**2 - m(n)**2)
     &           * ( - m(n-2)**2 * m(n-1)**2 * d1
     &               + (m(n-2)**2 + m(n-1)**2) * d2
     &               - d3 )

C...Check constraints on Cj,Dj
          c1 = 0.D0
          d1 = 0.D0
          d2 = 0.D0
          d3 = 0.D0
          DO j=1,n
            c1 = c1 + c(j)
            d1 = d1 + d(j)/m(j)**2
            d2 = d2 + d(j)
            d3 = d3 + d(j)*m(j)**2
          ENDDO
          print *, 'ZERO?? c1,d1,d2,d3=',c1,d1,d2,d3

        ENDIF



        u0 = 0.D0
        u2 = 0.D0
C...Momentum space wavefunction
        IF (mode.EQ.1) THEN
          DO j=1,n
            u0 = u0 + DSQRT(2.D0/pi) * C(j) / (q**2 + m(j)**2)
C...Note: u2 should be positive with the wfn defined as u0 - u2 S12...
            u2 = u2 - DSQRT(2.D0/pi) * D(j) / (q**2 + m(j)**2)
          ENDDO
C...Coordinate space wavefunction
        ELSE IF (mode.EQ.2) THEN
          u0 = 0.D0
          u2 = 0.D0
          IF (q.GT.0.D0) THEN
            DO j=1,n
              u0 = u0 + C(j) * DEXP(-m(j)*q)
              u2 = u2 + D(j) * DEXP(-m(j)*q)
     &                  * (1.D0 + 3.D0/(m(j)*q) + 3.D0/(m(j)**2*q**2))
            ENDDO
          ENDIF
        ENDIF

        RETURN
        END
