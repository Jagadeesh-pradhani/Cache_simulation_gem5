C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE GAMMUL(MU,MODE,X,RESULT)
C     ===================================

      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER      MU

      INTEGER      MODE

      COMPLEX*16   RESULT(*),X(*)

C     %------------%
C     | PARAMETERS |
C     %------------%

      COMPLEX*16   ZERO,ONE

      PARAMETER   (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))

      COMPLEX*16   TWO,I

      PARAMETER   (TWO=(2.0D+0,0.0D+0),I=(0.0D+0,1.0D+0))

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL ZCOPY,ZAXPY,ZSCAL

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      IF (MODE.EQ.0) THEN

      IF (MU.EQ.1) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3,  I, X(10), 1, RESULT( 1), 1)
       CALL ZAXPY( 3,  I, X( 7), 1, RESULT( 4), 1)
       CALL ZAXPY( 3, -I, X( 4), 1, RESULT( 7), 1)
       CALL ZAXPY( 3, -I, X( 1), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.2) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3,  ONE, X(10), 1, RESULT( 1), 1)
       CALL ZAXPY( 3, -ONE, X( 7), 1, RESULT( 4), 1)
       CALL ZAXPY( 3, -ONE, X( 4), 1, RESULT( 7), 1)
       CALL ZAXPY( 3,  ONE, X( 1), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.3) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3,  I, X( 7), 1, RESULT( 1), 1)
       CALL ZAXPY( 3, -I, X(10), 1, RESULT( 4), 1)
       CALL ZAXPY( 3, -I, X( 1), 1, RESULT( 7), 1)
       CALL ZAXPY( 3,  I, X( 4), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.4) THEN
       CALL ZCOPY( 6, X,    1, RESULT,    1)
       CALL ZCOPY( 6, ZERO, 0, RESULT(7), 1)
       CALL ZSCAL( 6, TWO, RESULT, 1)

      END IF

      ELSE IF (MODE.EQ.1) THEN

      IF (MU.EQ.1) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3, -I, X(10), 1, RESULT( 1), 1)
       CALL ZAXPY( 3, -I, X( 7), 1, RESULT( 4), 1)
       CALL ZAXPY( 3,  I, X( 4), 1, RESULT( 7), 1)
       CALL ZAXPY( 3,  I, X( 1), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.2) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3, -ONE, X(10), 1, RESULT( 1), 1)
       CALL ZAXPY( 3,  ONE, X( 7), 1, RESULT( 4), 1)
       CALL ZAXPY( 3,  ONE, X( 4), 1, RESULT( 7), 1)
       CALL ZAXPY( 3, -ONE, X( 1), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.3) THEN
       CALL ZCOPY(12, X, 1, RESULT, 1)
       CALL ZAXPY( 3, -I, X( 7), 1, RESULT( 1), 1)
       CALL ZAXPY( 3,  I, X(10), 1, RESULT( 4), 1)
       CALL ZAXPY( 3,  I, X( 1), 1, RESULT( 7), 1)
       CALL ZAXPY( 3, -I, X( 4), 1, RESULT(10), 1)

      ELSE IF (MU.EQ.4) THEN
       CALL ZCOPY( 6, ZERO, 0, RESULT   , 1)
       CALL ZCOPY( 6, X(7), 1, RESULT(7), 1)
       CALL ZSCAL( 6,  TWO,    RESULT(7), 1)

      END IF

      END IF

C     %---------------%
C     | END OF GAMMUL |
C     %---------------%

      RETURN
      END
