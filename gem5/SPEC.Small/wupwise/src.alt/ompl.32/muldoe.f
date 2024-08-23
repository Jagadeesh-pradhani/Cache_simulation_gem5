C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE MULDOE(U,N1,N2,N3,N4,X,RESULT)
C     =========================================

C --------------------------------------------------------------------
C
C     MULDOE MULTIPLIES MATRIX DOE BY VECTOR X
C     AND GIVES BACK VECTOR RESULT.
C
C --------------------------------------------------------------------


      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER         N1,N2,N3,N4

      COMPLEX*16      X(12,N1/2,N2,N3,N4)

      COMPLEX*16      RESULT(12,N1/2,N2,N3,N4)

      COMPLEX*16      U(3,3,4,N1,N2,N3,N4)

C     %------------%
C     | PARAMETERS |
C     %------------%

      COMPLEX*16      ONE
      PARAMETER      (ONE=(1.0D+0,0.0D+0))

C     %------------------------%
C     | LOCAL SCALARS & ARRAYS |
C     %------------------------%

      INTEGER      I,J,K,L,JKL
      INTEGER      IP,JP,KP,LP
      INTEGER      IM,JM,KM,LM

      COMPLEX*16   AUX1(12),AUX2(12),AUX3(12)

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL  GAMMUL,SU3MUL
      EXTERNAL  ZAXPY,ZCOPY

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

C     %---------------------%
C     | LOOP OVER ODD SITES |
C     %---------------------%

C$OMP PARALLEL
C$OMP+        PRIVATE (AUX1, AUX2, AUX3),
C$OMP+        PRIVATE (I, IM, IP, J, JM, JP, K, KM, KP, L, LM, LP),
C$OMP+        SHARED (N1, N2, N3, N4, RESULT, U, X)

C$OMP DO
      DO 100 JKL = 0, N2 * N3 * N4 - 1

       L = MOD (JKL / (N2 * N3), N4) + 1
       LP=MOD(L,N4)+1

        K = MOD (JKL / N2, N3) + 1
        KP=MOD(K,N3)+1

         J = MOD (JKL, N2) + 1
         JP=MOD(J,N2)+1

         DO 100 I=(MOD(J+K+L,2)+1),N1,2

           IP=MOD(I,N1)+1

           CALL GAMMUL(1,0,X(1,(IP+1)/2,J,K,L),AUX1)
           CALL SU3MUL(U(1,1,1,I,J,K,L),'N',AUX1,AUX3)

           CALL GAMMUL(2,0,X(1,(I+1)/2,JP,K,L),AUX1)
           CALL SU3MUL(U(1,1,2,I,J,K,L),'N',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL GAMMUL(3,0,X(1,(I+1)/2,J,KP,L),AUX1)
           CALL SU3MUL(U(1,1,3,I,J,K,L),'N',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL GAMMUL(4,0,X(1,(I+1)/2,J,K,LP),AUX1)
           CALL SU3MUL(U(1,1,4,I,J,K,L),'N',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL ZCOPY(12,AUX3,1,RESULT(1,(I+1)/2,J,K,L),1)

 100  CONTINUE
C$OMP END DO


C$OMP DO
      DO 200 JKL = 0, N2 * N3 * N4 - 1

       L = MOD (JKL / (N2 * N3), N4) + 1
       LM=L-1
       IF(LM.EQ.0) LM=N4

        K = MOD (JKL / N2, N3) + 1
        KM=K-1
        IF(KM.EQ.0) KM=N3

         J = MOD (JKL, N2) + 1
         JM=J-1
         IF(JM.EQ.0) JM=N2

         DO 200 I=(MOD(J+K+L,2)+1),N1,2

           IM=I-1
           IF(IM.EQ.0) IM=N1

           CALL GAMMUL(1,1,X(1,(IM+1)/2,J,K,L),AUX1)
           CALL SU3MUL(U(1,1,1,IM,J,K,L),'C',AUX1,AUX3)

           CALL GAMMUL(2,1,X(1,(I+1)/2,JM,K,L),AUX1)
           CALL SU3MUL(U(1,1,2,I,JM,K,L),'C',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL GAMMUL(3,1,X(1,(I+1)/2,J,KM,L),AUX1)
           CALL SU3MUL(U(1,1,3,I,J,KM,L),'C',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL GAMMUL(4,1,X(1,(I+1)/2,J,K,LM),AUX1)
           CALL SU3MUL(U(1,1,4,I,J,K,LM),'C',AUX1,AUX2)
           CALL ZAXPY(12,ONE,AUX2,1,AUX3,1)

           CALL ZAXPY(12,ONE,AUX3,1,RESULT(1,(I+1)/2,J,K,L),1)

 200  CONTINUE
C$OMP END DO

C$OMP END PARALLEL

C     %---------------%
C     | END OF MULDOE |
C     %---------------%

      RETURN
      END
