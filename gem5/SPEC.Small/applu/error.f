c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine error

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   compute the solution error
c
c---------------------------------------------------------------------

      implicit none

      include 'applu.incl'

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, m
      double precision  tmp(5)
      double precision  u000ijk(5)
      double precision  errnm1,errnm2,errnm3,errnm4,errnm5


      errnm1 = 0.0d+00
      errnm2 = 0.0d+00
      errnm3 = 0.0d+00
      errnm4 = 0.0d+00
      errnm5 = 0.0d+00

!$omp parallel do default(shared) private(i,j,k,m,tmp,u000ijk)
!$omp& reduction(+:errnm1,errnm2,errnm3,errnm4,errnm5)
      do k = 2, nz-1
         do j = jst, jend
            do i = ist, iend
               call exact( i, j, k, u000ijk )
               do m = 1, 5
                  tmp(m) = ( u000ijk(m) - u(m,i,j,k) )
               end do
               errnm1 = errnm1 + tmp(1) * tmp(1)
               errnm2 = errnm2 + tmp(2) * tmp(2)
               errnm3 = errnm3 + tmp(3) * tmp(3)
               errnm4 = errnm4 + tmp(4) * tmp(4)
               errnm5 = errnm5 + tmp(5) * tmp(5)
            end do
         end do
      end do
!$omp end parallel do

      errnm(1) = errnm1
      errnm(2) = errnm2
      errnm(3) = errnm3
      errnm(4) = errnm4
      errnm(5) = errnm5

      do m = 1, 5
         errnm(m) = sqrt ( errnm(m) / ( (nx0-2)*(ny0-2)*(nz0-2) ) )
      end do

c        write (*,1002) ( errnm(m), m = 1, 5 )

 1002 format (1x/1x,'RMS-norm of error in soln. to ',
     > 'first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'second pde = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of error in soln. to ',
     > 'fifth pde  = ',1pe12.5)

      return
      end
