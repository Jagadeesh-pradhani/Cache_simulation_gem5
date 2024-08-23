C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012



      SUBROUTINE INIT(SOURCE,IX,IY,IZ,IT,IC,ID,SEED,NITER,KAPPA)
C     ==========================================================

C --------------------------------------------------------------------
C
C INITIALIZATION AND INPUT FOR LINSOLV
C
C --------------------------------------------------------------------

      IMPLICIT NONE

C     %-----------%
C     | ARGUMENTS |
C     %-----------%

      INTEGER           SOURCE

      INTEGER           IX,IY,IZ,IT,IC,ID

      INTEGER           NITER

      INTEGER           SEED(4)

      COMPLEX*16        KAPPA

C     %-----------------------%
C     | EXECUTABLE STATEMENTS |
C     %-----------------------%

      OPEN(10,FILE='wupwise.in',STATUS='OLD')

        READ(10,*)     KAPPA
        READ(10,*)     NITER
        READ(10,*)     SOURCE
        READ(10,*)     IX,IY,IZ,IT
        READ(10,*)     IC,ID
        READ(10,*)     SEED(1),SEED(2),SEED(3),SEED(4)

      CLOSE(10, STATUS='KEEP')

      IF(SOURCE.EQ.1) THEN
         WRITE(*,*) ' POINT SOURCE'
         WRITE(*,*) ' SPACE INDEX = ',IX,IY,IZ,IT
         WRITE(*,*) ' COLOR INDEX = ',IC
         WRITE(*,*) ' DIRAC INDEX = ',ID
      ELSE IF(SOURCE.EQ.2) THEN
         WRITE(*,*) ' RANDOM SOURCE'
      ENDIF

C     WRITE(*,*) ' SEED  : ',SEED   ! KJS

      WRITE(*,991) '  SEED : ',SEED
991   FORMAT(A9,4I8)

      WRITE(*,*) ' KAPPA : ',KAPPA
      WRITE(*,*) ' NITER : ',NITER

C     %-------------%
C     | END OF INIT |
C     %-------------%

      RETURN
      END

