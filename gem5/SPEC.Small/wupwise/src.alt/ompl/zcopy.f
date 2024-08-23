      subroutine  zcopy(n,zx,incx,zy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 4/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex*16 zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
C=======================================================================
C
C The following is a parallel (OpenMP) version written by Paul Kinney
C (paul.kinney@sun.com).
C
      SUBROUTINE P_ZCOPY (N, ZX, INCX, ZY, INCY)
C
      COMPLEX*16 ZX(*), ZY(*)
      INTEGER BX, BY, I, INCX, INCY, N
C
      IF (N .LE. 0) RETURN
C
C$OMP PARALLEL
C$OMP+  PRIVATE (BX, BY, I),
C$OMP+  SHARED (INCX, INCY, N, ZX, ZY)
      IF ((INCX .EQ. 1) .AND. (INCY .EQ. 1)) THEN
C
C        code for increment not equal to 1
C
C$OMP   DO
        DO I = 1, N
          ZY(I) = ZX(I)
        END DO
C$OMP   END DO
      ELSE
C
C        code for increment equal to 1
C
        IF (INCX .LT. 0) THEN
          BX = 1 - N * INCX
        ELSE
          BX = 1 - INCX
        END IF
        IF (INCY .LT. 0) THEN
          BY = 1 - N * INCY
        ELSE
          BY = 1 - INCY
        END IF
C$OMP   DO
        DO I = 1, N
          ZY(BY+I*INCY) = ZX(BX+I*INCX)
        END DO
C$OMP   END DO
      END IF
C$OMP END PARALLEL

      RETURN
      END
