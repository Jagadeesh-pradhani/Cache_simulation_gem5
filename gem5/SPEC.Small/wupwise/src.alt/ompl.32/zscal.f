      subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex*16 za,zx(*)
      integer i,incx,ix,n
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
C=======================================================================
C
C The following is a parallel (OpenMP) version written by Paul Kinney
C (paul.kinney@sun.com).
C
      SUBROUTINE  P_ZSCAL (N, ZA, ZX, INCX)
C
      COMPLEX*16 ZA,ZX(*)
      INTEGER BX, I, INCX, N
C
      IF( N.LE.0 .OR. INCX.LE.0 )RETURN
C
C$OMP PARALLEL
C$OMP+  PRIVATE (BX, I),
C$OMP+  SHARED (INCX, N, ZA, ZX)
      IF (INCX .EQ. 1) THEN
C
C        code for increment equal to 1
C
C$OMP   DO
        DO I = 1, N
          ZX(I) = ZA * ZX(I)
        END DO
C$OMP   END DO
      ELSE
C
C        code for increment not equal to 1
C
        IF (INCX .LT. 0) THEN
          BX = 1 - N * INCX
        ELSE
          BX = 1 - INCX
        END IF
C$OMP   DO
        DO I = 1, N
          ZX(BX+I*INCX) = ZA * ZX(BX+I*INCX)
        END DO
C$OMP   END DO
      END IF
C$OMP END PARALLEL

      RETURN
      END
