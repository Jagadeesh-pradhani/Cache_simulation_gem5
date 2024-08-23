c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine ssor

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   to perform pseudo-time stepping SSOR iterations
c   for five nonlinear pde's.
c---------------------------------------------------------------------

      implicit none


      include 'applu.incl'

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, l, m, mt, n, nt
      integer ipx, npx, i0, i1
      integer ipy, npy, j0, j1
      integer istep
      double precision  tmp, tmp2
      double precision  delunm(5), tv(5,isiz1/2*2+1,isiz2)
      double precision  tvu(5,isiz1/2*2+1,isiz2)


c---------------------------------------------------------------------
c  Useful OpenMP routines
c---------------------------------------------------------------------
!$      external omp_get_thread_num
!$      integer omp_get_thread_num
!$      external omp_get_num_threads
!$      integer omp_get_num_threads


c---------------------------------------------------------------------
c   begin pseudo-time stepping iterations
c---------------------------------------------------------------------
      tmp = 1.0d+00 / ( omega * ( 2.0d+00 - omega ) )

c---------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
c---------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(M,N,I,J)
!$OMP DO
      do j=jst,jend
         do i=ist,iend
            do m=1,5
               do n=1,5
                  a(m,n,i,j) = 0.d0
                  b(m,n,i,j) = 0.d0
                  c(m,n,i,j) = 0.d0
                  d(m,n,i,j) = 0.d0
               enddo
            enddo
         enddo
      enddo
!$OMP END DO NOWAIT
!$OMP DO
      do j=jend,jst,-1
         do i=iend,ist,-1
            do m=1,5
               do n=1,5
                  au(m,n,i,j) = 0.d0
                  bu(m,n,i,j) = 0.d0
                  cu(m,n,i,j) = 0.d0
                  du(m,n,i,j) = 0.d0
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
      call rhs

c---------------------------------------------------------------------
c   compute the L2 norms of newton iteration residuals
c---------------------------------------------------------------------
      call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >             ist, iend, jst, jend,
     >             rsd, rsdnm )


      if ( ipr .eq. 1 ) then
         write (*,*) '          Initial residual norms'
         write (*,*)
         write (*,1007) ( rsdnm(m), m = 1, 5 )
        write (*,'(/a)') 'Iteration RMS-residual of 5th PDE'
      end if



c---------------------------------------------------------------------
c   the timestep loop
c---------------------------------------------------------------------
      do istep = 1, itmax


         if ( ( mod ( istep, inorm ) .eq. 0 ) .and.
     >          ipr .eq. 1 ) then
             write ( *, 1001 ) istep
         end if
         if (mod ( istep, 20) .eq. 0 .or.
     >         istep .eq. itmax .or.
     >         istep .eq. 1) then
            write( *, 200) istep
 200        format(' Time step ', i4)
         endif

c---------------------------------------------------------------------
c   perform SSOR iteration
c---------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(M,I,J,K,tmp2)
         tmp2 = dt
!$omp do
         do k = 2, nz - 1
            do j = jst, jend
               do i = ist, iend
                  do m = 1, 5
                     rsd(m,i,j,k) = tmp2 * rsd(m,i,j,k)
                  end do
               end do
            end do
         end do
!$omp end do
!$OMP END PARALLEL


!$omp parallel
!$omp&  default (shared)
!$omp&  private (i0, i1, ipx, ipy, j0, j1, k, l, mt, nt, npx, npy)
!$omp&  shared (nx, ny, nz, omega)

c-----------------------------------------------------------------------
c Partition X-Y plane amongst processors in two dimensions.
c-----------------------------------------------------------------------
c
c		without OpenMP
c
c	Vishal:  Need to change mt to 0
         nt = 1
         mt = 0
c
c		with OpenMP
c
!$       nt = omp_get_num_threads ()
!$       mt = omp_get_thread_num ()
         npx = int (sqrt (dble (nt)))
         do while ((npx .gt. 1) .and. (mod (nt, npx) .ne. 0))
           npx = npx - 1
         end do
         npy = nt / npx
         ipx = mod (mt, npx)
         ipy = mt / npx
         i0 = (ipx + 0) * (nx - 2) / npx + 2
         i1 = (ipx + 1) * (nx - 2) / npx + 1
         j0 = (ipy + 0) * (ny - 2) / npy + 2
         j1 = (ipy + 1) * (ny - 2) / npy + 1

         DO l = 2, npx + npy + nz - 3
          k = l - ipx - ipy
          if ((1 .lt. k) .and. (k .lt. nz)) then
c---------------------------------------------------------------------
c   form the lower triangular part of the jacobian matrix
c---------------------------------------------------------------------
            call jacld (i0, i1, j0, j1, k)

c---------------------------------------------------------------------
c   perform the lower triangular solution
c---------------------------------------------------------------------
            call blts( isiz1, isiz2, isiz3,
     >                 nx, ny, nz, i0, i1, j0, j1, k,
     >                 omega,
     >                 rsd, tv,
     >                 a, b, c, d )
           end if
!$omp barrier
          end do

         DO l = npx + npy + nz - 3, 2, -1
          k = l - ipx - ipy
          if ((1 .lt. k) .and. (k .lt. nz)) then

c---------------------------------------------------------------------
c   form the strictly upper triangular part of the jacobian matrix
c---------------------------------------------------------------------
            call jacu(i0, i1, j0, j1, k)

c---------------------------------------------------------------------
c   perform the upper triangular solution
c---------------------------------------------------------------------
            call buts( isiz1, isiz2, isiz3,
     >                 nx, ny, nz, i0, i1, j0, j1, k,
     >                 omega,
     >                 rsd, tvu,
     >                 du, au, bu, cu )
           end if
!$omp barrier
          END DO
!$omp end parallel

c---------------------------------------------------------------------
c   update the variables
c---------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(M,I,J,K,tmp2)
         tmp2 = tmp
!$omp do
         do k = 2, nz-1
            do j = jst, jend
               do i = ist, iend
                  do m = 1, 5
                     u( m, i, j, k ) = u( m, i, j, k )
     >                    + tmp2 * rsd( m, i, j, k )
                  end do
               end do
            end do
         end do
!$omp end do nowait
!$OMP END PARALLEL

c---------------------------------------------------------------------
c   compute the max-norms of newton iteration corrections
c---------------------------------------------------------------------
         if ( mod ( istep, inorm ) .eq. 0 ) then
            call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >                   ist, iend, jst, jend,
     >                   rsd, delunm )
            if ( ipr .eq. 1 ) then
                write (*,1006) ( delunm(m), m = 1, 5 )
            else if ( ipr .eq. 2 ) then
                write (*,'(i5,f15.6)') istep,delunm(5)
            end if
         end if

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
         call rhs

c---------------------------------------------------------------------
c   compute the max-norms of newton iteration residuals
c---------------------------------------------------------------------
         if ( ( mod ( istep, inorm ) .eq. 0 ) .or.
     >        ( istep .eq. itmax ) ) then
            call l2norm( isiz1, isiz2, isiz3, nx0, ny0, nz0,
     >                   ist, iend, jst, jend,
     >                   rsd, rsdnm )
            if ( ipr .eq. 1 ) then
                write (*,1007) ( rsdnm(m), m = 1, 5 )
            end if
         end if

c---------------------------------------------------------------------
c   check the newton-iteration residuals against the tolerance levels
c---------------------------------------------------------------------
         if ( ( rsdnm(1) .lt. tolrsd(1) ) .and.
     >        ( rsdnm(2) .lt. tolrsd(2) ) .and.
     >        ( rsdnm(3) .lt. tolrsd(3) ) .and.
     >        ( rsdnm(4) .lt. tolrsd(4) ) .and.
     >        ( rsdnm(5) .lt. tolrsd(5) ) ) then
            if (ipr .eq. 1 ) then
               write (*,1004) istep
            end if
            return
         end if

      end do




      return

 1001 format (1x/5x,'pseudo-time SSOR iteration no.=',i4/)
 1004 format (1x/1x,'convergence was achieved after ',i4,
     >   ' pseudo-time steps' )
 1006 format (1x/1x,'RMS-norm of SSOR-iteration correction ',
     > 'for first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for second pde = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of SSOR-iteration correction ',
     > 'for fifth pde  = ',1pe12.5)
 1007 format (1x/1x,'RMS-norm of steady-state residual for ',
     > 'first pde  = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'second pde = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'third pde  = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'fourth pde = ',1pe12.5/,
     > 1x,'RMS-norm of steady-state residual for ',
     > 'fifth pde  = ',1pe12.5)

      end
