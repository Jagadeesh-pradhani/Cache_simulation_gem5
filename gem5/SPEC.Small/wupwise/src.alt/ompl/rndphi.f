C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE RNDPHI(PHI,LENGTH,SEED)
C     ==================================

C --------------------------------------------------------------------
C
C SET UP RANDOM SOURCE VECTOR
C
C --------------------------------------------------------------------


      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER           LENGTH

      INTEGER           SEED(4)

      REAL*8            PHI(*)

C     %------------------------%
C     | LOCAL SCALARS & ARRAYS |
C     %------------------------%

      INTEGER           I

C     %--------------------%
C     | EXTERNAL FUNCTIONS |
C     %--------------------%

      EXTERNAL    DLARAN
      REAL*8      DLARAN

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

C     %-------------------------------------%
C     | DLARAN RETURNS A RANDOM REAL NUMBER |
C     | FROM A UNIFORM (0,1) DISTRIBUTION.  |
C     %-------------------------------------%

!$OMP PARALLEL DO
      DO I=1,LENGTH
        PHI(I) = 0
      ENDDO
      DO 100 I=1,LENGTH
        PHI(I) = DLARAN(SEED)
 100  CONTINUE

C     %---------------%
C     | END OF RNDPHI |
C     %---------------%

      RETURN
      END
