
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine l2norm ( ldx, ldy, ldz, 
     >                    nx0, ny0, nz0,
     >                    ist, iend, 
     >                    jst, jend,
     >                    v, sum )
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   to compute the l2-norm of vector v.
c---------------------------------------------------------------------

      implicit none


c---------------------------------------------------------------------
c  input parameters
c---------------------------------------------------------------------
      integer ldx, ldy, ldz
      integer nx0, ny0, nz0
      integer ist, iend
      integer jst, jend
c---------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
c---------------------------------------------------------------------
      double precision  v(5,ldx/2*2+1,ldy/2*2+1,*), sum(5)

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      double precision  sum1,sum2,sum3,sum4,sum5
      integer i, j, k, m


      sum1 = 0.0d+00
      sum2 = 0.0d+00
      sum3 = 0.0d+00
      sum4 = 0.0d+00
      sum5 = 0.0d+00

!$omp parallel do default(shared) private(i,j,k)
!$omp& reduction(+:sum1,sum2,sum3,sum4,sum5)
      do k = 2, nz0-1
         do j = jst, jend
            do i = ist, iend
               sum1 = sum1 + v(1,i,j,k) * v(1,i,j,k)
               sum2 = sum2 + v(2,i,j,k) * v(2,i,j,k)
               sum3 = sum3 + v(3,i,j,k) * v(3,i,j,k)
               sum4 = sum4 + v(4,i,j,k) * v(4,i,j,k)
               sum5 = sum5 + v(5,i,j,k) * v(5,i,j,k)
            end do
         end do
      end do
!$omp end parallel do

      sum(1) = sum1
      sum(2) = sum2
      sum(3) = sum3
      sum(4) = sum4
      sum(5) = sum5

      do m = 1, 5
         sum(m) = sqrt ( sum(m) / ( (nx0-2)*(ny0-2)*(nz0-2) ) )
      end do

      return
      end
