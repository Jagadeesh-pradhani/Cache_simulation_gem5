C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE UINITH(N1,N2,N3,N4,SEED,U)
C     =====================================

C --------------------------------------------------------------------
C
C INITIALIZES HOT START
C
C --------------------------------------------------------------------


      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER           N1,N2,N3,N4

      INTEGER           SEED(4)

      COMPLEX*16        U(3,3,4,N1,N2,N3,N4)

C     %------------%
C     | PARAMETERS |
C     %------------%

      COMPLEX*16        ZERO,ONE

      PARAMETER        (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))

C     %---------------%
C     | LOCAL SCALARS |
C     %---------------%

      INTEGER           I,J,K,L

      INTEGER           MU

      INTEGER           LENGTH

      COMPLEX*16        ALPHA, BETA

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL  RNDCNF

      EXTERNAL  ZAXPY,ZSCAL,ZCOPY

C     %--------------------%
C     | EXTERNAL FUNCTIONS |
C     %--------------------%

      EXTERNAL     DZNRM2
      REAL*8       DZNRM2

      EXTERNAL     ZDOTC
      COMPLEX*16   ZDOTC

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      LENGTH=72*N1*N2*N3*N4

      CALL RNDCNF(U,LENGTH,SEED)

C     %---------------------------------------------------%
C     | UNITARIZES SU(3)-MATRIZES USING SCHMIDT-ALGORITHM |
C     %---------------------------------------------------%
!$OMP PARALLEL DO PRIVATE(L,K,J,I,MU,BETA,ALPHA)
      DO 100 L=1,N4
        DO 100 K=1,N3
          DO 100 J=1,N2
            DO 100 I=1,N1

              DO 100 MU=1,4

                BETA = DZNRM2(3,U(1,1,MU,I,J,K,L),1)

                CALL ZSCAL(3,ONE/BETA,U(1,1,MU,I,J,K,L),1)

                ALPHA = ZDOTC(3,U(1,1,MU,I,J,K,L),1,
     .                          U(1,2,MU,I,J,K,L),1)

                CALL ZAXPY(3,-ALPHA,U(1,1,MU,I,J,K,L),1,
     .                              U(1,2,MU,I,J,K,L),1)

                BETA = DZNRM2(3,U(1,2,MU,I,J,K,L),1)

                CALL ZSCAL(3,ONE/BETA,U(1,2,MU,I,J,K,L),1)

                CALL ZCOPY(3,ZERO,0,U(1,3,MU,I,J,K,L),1)

                U(1,3,MU,I,J,K,L)=
     .          DCONJG(U(2,1,MU,I,J,K,L)*U(3,2,MU,I,J,K,L)
     .                -U(3,1,MU,I,J,K,L)*U(2,2,MU,I,J,K,L))

                U(2,3,MU,I,J,K,L)=
     .          DCONJG(U(3,1,MU,I,J,K,L)*U(1,2,MU,I,J,K,L)
     .                -U(1,1,MU,I,J,K,L)*U(3,2,MU,I,J,K,L))

                U(3,3,MU,I,J,K,L)=
     .          DCONJG(U(1,1,MU,I,J,K,L)*U(2,2,MU,I,J,K,L)
     .                -U(2,1,MU,I,J,K,L)*U(1,2,MU,I,J,K,L))

 100  CONTINUE

C     %---------------%
C     | END OF UINITH |
C     %---------------%

      RETURN
      END
