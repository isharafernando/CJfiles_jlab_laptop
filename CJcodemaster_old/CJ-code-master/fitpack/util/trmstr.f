       SUBROUTINE TRMSTR(STRING,ILEN)
C                                                   -=-=- trmstr

C      Removes leading spaces and returns true length of a character string

       CHARACTER STRING*(*),SPACE*1

       DATA SPACE/' '/

       ILEN=0

       IF (STRING.EQ.SPACE) RETURN
C                                          Remove leading spaces
1      IF (STRING(1:1).NE.SPACE) GOTO 2
          STRING=STRING(2:)
          GOTO 1
2      CONTINUE
C                                          Count up trailing spaces
       DO 3 I=LEN(STRING),1,-1
          IF (STRING(I:I).NE.SPACE) THEN
             ILEN=I
             RETURN
          END IF
3      CONTINUE
C               *************************
       END

