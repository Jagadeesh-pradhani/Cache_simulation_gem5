C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE RNDCNF(U,LENGTH,SEED)
C     ================================

C --------------------------------------------------------------------
C
C SET UP HOT GAUGE FIELD
C
C --------------------------------------------------------------------

      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER           LENGTH

      INTEGER           SEED(4)

      REAL*8            U(*)

C     %------------------------%
C     | LOCAL SCALARS & ARRAYS |
C     %------------------------%

      INTEGER           I

C     %--------------------%
C     | EXTERNAL FUNCTIONS |
C     %--------------------%

      EXTERNAL         DLARND
      REAL*8           DLARND

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

C     %-------------------------------------%
C     | DLARND RETURNS A RANDOM REAL NUMBER |
C     | FROM A UNIFORM (-1,1) DISTRIBUTION. |
C     %-------------------------------------%
!$OMP PARALLEL DO
      DO I=1,LENGTH
        U(I) = 0.0
      ENDDO
      DO 100 I=1,LENGTH
        U(I) = DLARND(2,SEED)
 100  CONTINUE

C     %---------------%
C     | END OF RNDCNF |
C     %---------------%

      RETURN
      END


