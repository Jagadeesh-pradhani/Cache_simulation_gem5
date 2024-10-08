*     DOUBLE PRECISION FUNCTION DLARND( IDIST, ISEED )
      REAL*8 FUNCTION DLARND( IDIST, ISEED )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
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
*     .. Scalar Arguments ..
      INTEGER            IDIST
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARND returns a random real number from a uniform or normal
*  distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  uniform (0,1)
*          = 2:  uniform (-1,1)
*          = 3:  normal (0,1)
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine DLARAN to generate a random
*  real number from a uniform (0,1) distribution. The Box-Muller method
*  is used to transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. Parameters ..
*     DOUBLE PRECISION   ONE, TWO
      REAL*8		 ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
*     DOUBLE PRECISION   TWOPI
      REAL*8		 TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663D+0 )
*     ..
*     .. Local Scalars ..
*     DOUBLE PRECISION   T1, T2
      REAL*8 		 T1, T2
*     ..
*     .. External Functions ..
*     DOUBLE PRECISION   DLARAN
      REAL*8	         DLARAN
      EXTERNAL           DLARAN
*     ..
*     .. Intrinsic Functions ..
C     INTRINSIC          COS, LOG, SQRT
C     Not allowed in strict Fortran77. RW
*     ..
*     .. Executable Statements ..
*
*     Generate a real random number from a uniform (0,1) distribution
*
      T1 = DLARAN( ISEED )
*
      IF( IDIST.EQ.1 ) THEN
*
*        uniform (0,1)
*
         DLARND = T1
      ELSE IF( IDIST.EQ.2 ) THEN
*
*        uniform (-1,1)
*
         DLARND = TWO*T1 - ONE
      ELSE IF( IDIST.EQ.3 ) THEN
*
*        normal (0,1)
*
         T2 = DLARAN( ISEED )
         DLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 )
      END IF
      RETURN
*
*     End of DLARND
*
      END
