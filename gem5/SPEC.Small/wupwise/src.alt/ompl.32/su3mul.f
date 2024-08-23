C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE SU3MUL(U,TRANSU,X,RESULT)
C     ====================================

C --------------------------------------------------------------------
C
C  MULTIPLIES TWO SU(3) MATRICES: RESULT = U(TRANSU) * X
C
C --------------------------------------------------------------------

      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      CHARACTER*1  TRANSU

      COMPLEX*16   U(3,*),X(*),RESULT(*)

C     %------------%
C     | PARAMETERS |
C     %------------%

      COMPLEX*16   ZERO,ONE

      PARAMETER   (ZERO=(0.0D+0, 0.0D+0),ONE=(1.0D+0, 0.0D+0))

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL       ZGEMM

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      CALL ZGEMM(TRANSU, 'NO TRANSPOSE',3,4,3,
     .           ONE,U,3,X,3,ZERO,RESULT,3)

C     %---------------%
C     | END OF SU3MUL |
C     %---------------%

      RETURN
      END
