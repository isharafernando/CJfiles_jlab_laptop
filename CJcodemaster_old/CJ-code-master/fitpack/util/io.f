
      FUNCTION NEXTUN()
C                                                   -=-=- nextun
C                                    Returns an unallocated FORTRAN i/o unit.
      LOGICAL EX
C
      DO 1  N = 10, 98
          INQUIRE (UNIT=N, OPENED=EX)
          IF (.NOT. EX) then
             nextun = n
             RETURN
          end if
   1  CONTINUE

      Stop 'NextUnit number not found!  Stopped in NextUn().'

C               *************************
      END
