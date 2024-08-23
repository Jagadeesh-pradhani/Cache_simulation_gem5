C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE PHINIT(N1,N2,N3,N4,MODE,IX,IY,IZ,IT,IC,ID,SEED,PHI)
C     ==============================================================

C --------------------------------------------------------------------
C
C SET UP SOURCE VECTOR
C
C --------------------------------------------------------------------


      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER           N1,N2,N3,N4

      INTEGER           MODE

      INTEGER           IX,IY,IZ,IT,IC,ID

      INTEGER           SEED(4)

      COMPLEX*16        PHI(12,N1,N2,N3,N4)

C     %------------%
C     | PARAMETERS |
C     %------------%

      COMPLEX*16        ZERO,ONE

      PARAMETER        (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))

C     %---------------%
C     | LOCAL SCALARS |
C     %---------------%

      INTEGER           SIZE

C     %----------------------%
C     | EXTERNAL SUBROUTINES |
C     %----------------------%

      EXTERNAL  RNDPHI
      EXTERNAL  ZCOPY

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      SIZE=12*N1*N2*N3*N4

      IF(MODE.EQ.1) THEN

C       %---------------------%
C       | SET UP POINT SOURCE |
C       %---------------------%

        CALL P_ZCOPY(SIZE,ZERO,0,PHI,1)

        PHI(3*(ID-1)+IC,IX,IY,IZ,IT)=ONE

      ELSE IF(MODE.EQ.2) THEN

C       %----------------------%
C       | SET UP RANDOM SOURCE |
C       %----------------------%

        CALL RNDPHI(PHI,2*SIZE,SEED)

      ENDIF

C     %---------------%
C     | END OF PHINIT |
C     %---------------%

      RETURN
      END


