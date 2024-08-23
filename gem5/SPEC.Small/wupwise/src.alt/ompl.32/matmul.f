C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE MATMUL(U,N1,N2,N3,N4,KAPPA,X,AUX,RESULT)
C     ===================================================

C --------------------------------------------------------------------
C
C     MATMULT MULTIPLIES MATRIX ME(U) = I-KAPPA^2 * DE(U) BY VECTOR X
C     AND GIVES BACK VECTOR RESULT.
C
C --------------------------------------------------------------------


      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER         N1,N2,N3,N4

      COMPLEX*16      KAPPA

      COMPLEX*16      X(12,N1/2,N2,N3,N4)

      COMPLEX*16      AUX(12,N1/2,N2,N3,N4)

      COMPLEX*16      RESULT(12,N1/2,N2,N3,N4)

      COMPLEX*16      U(3,3,4,N1,N2,N3,N4)

C     %------------------------%
C     | LOCAL SCALARS & ARRAYS |
C     %------------------------%

      INTEGER      SIZE

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL  MULDEO,MULDOE
      EXTERNAL  P_ZAXPY,P_ZCOPY

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      SIZE=6*N1*N2*N3*N4

      CALL MULDOE(U,N1,N2,N3,N4,X,RESULT)
      CALL MULDEO(U,N1,N2,N3,N4,RESULT,AUX)

      CALL P_ZCOPY(SIZE,X,1,RESULT,1)
      CALL P_ZAXPY(SIZE,-KAPPA**2,AUX,1,RESULT,1)

C     %----------------%
C     | END OF  MATMUL |
C     %----------------%

      RETURN
      END


