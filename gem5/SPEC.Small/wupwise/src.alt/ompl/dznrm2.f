*     DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
      REAL*8 FUNCTION DZNRM2( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      COMPLEX*16                        X( * )
*     ..
*
*  DZNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DZNRM2 := sqrt( conjg( x' )*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to ZLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*     Modified for SPEC CPU2000: Replaced
*	DOUBLE PRECICION
*     by
*	REAL*8
*     since for benchmark fairness reasons, SPEC wants explicit control
*     over the length of variables, even though the Fortran standards
*     do not define REAL*8.
*     All common Fortran compilers, however, recognize this keyword.
*
*     Reinhold Weicker, Sept. 24, 1999
*
*
*     .. Parameters ..
*     DOUBLE PRECISION      ONE         , ZERO
      REAL*8	           ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
*     DOUBLE PRECISION      NORM, SCALE, SSQ, TEMP
      REAL*8	            NORM, SCALE, SSQ, TEMP
*     .. Intrinsic Functions ..
C     INTRINSIC             ABS, DIMAG, DBLE, SQRT
C     DBLE not allowed in strict Fortran77
      INTRINSIC             ABS, DIMAG, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DZNRM2 = NORM
      RETURN
*
*     End of DZNRM2.
*
      END
*=======================================================================
*
* The following is a parallel (OpenMP) version written by Paul Kinney
* (paul.kinney@sun.com).
*
*     DOUBLE PRECISION FUNCTION P_DZNRM2 (N, X, INCX)
      REAL*8 FUNCTION P_DZNRM2 (N, X, INCX)
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      COMPLEX*16                        X( * )
*
      REAL*8               ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*
      INTEGER               IX
*
      REAL*8                LSCALE, LSSQ, NORM, SCALE, SSQ, TEMP
*
      INTRINSIC             ABS, DIMAG, SQRT
*
      NORM  = ZERO
      IF ((N .GE. 1) .AND. (INCX .GE. 1)) THEN
         SCALE = ZERO
         SSQ   = ONE
C$OMP    PARALLEL
C$OMP+     PRIVATE (IX, LSCALE, LSSQ, TEMP) ,
C$OMP+     SHARED (SCALE, SSQ, X)
         LSCALE = ZERO
         LSSQ   = ONE
*
C$OMP    DO
         DO IX = 1, 1 + (N - 1) *INCX, INCX
            IF (DBLE (X(IX)) .NE. ZERO) THEN
               TEMP = ABS (DBLE (X(IX)))
               IF (LSCALE.LT.TEMP) THEN
                  LSSQ   = ONE   + LSSQ* (LSCALE / TEMP) ** 2
                  LSCALE = TEMP
               ELSE
                  LSSQ   = LSSQ   + (TEMP / LSCALE) ** 2
               END IF
            END IF
            IF (DIMAG (X(IX)) .NE. ZERO) THEN
               TEMP = ABS (DIMAG (X(IX)))
               IF (LSCALE.LT.TEMP) THEN
                  LSSQ   = ONE   + LSSQ* (LSCALE / TEMP) ** 2
                  LSCALE = TEMP
               ELSE
                  LSSQ   = LSSQ   + (TEMP / LSCALE) ** 2
               END IF
            END IF
         END DO
C$OMP    END DO
C$OMP    CRITICAL
         IF (SCALE .LT. LSCALE) THEN
            SSQ = ((SCALE / LSCALE) ** 2) * SSQ + LSSQ
            SCALE = LSCALE
         ELSE
            SSQ = SSQ + ((LSCALE / SCALE) ** 2) * LSSQ
         END IF
C$OMP    END CRITICAL
C$OMP    END PARALLEL
         NORM  = SCALE * SQRT (SSQ)
      END IF
*
      P_DZNRM2 = NORM
      RETURN
*
      END
