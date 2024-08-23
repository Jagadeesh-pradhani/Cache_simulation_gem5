
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine read_input

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      implicit none

      include 'applu.incl'

      integer  fstatus


c---------------------------------------------------------------------
c    if input file does not exist, it uses defaults
c       ipr = 1 for detailed progress output
c       inorm = how often the norm is printed (once every inorm iterations)
c       itmax = number of pseudo time steps
c       dt = time step
c       omega 1 over-relaxation factor for SSOR
c       tolrsd = steady state residual tolerance levels
c       nx, ny, nz = number of grid points in x, y, z directions
c---------------------------------------------------------------------

         write(*, 1000)

chomp         open (unit=3,file='inputlu.data',status='old',
chomp     >         access='sequential',form='formatted', iostat=fstatus)
chomp         if (fstatus .eq. 0) then

chomp            write(*, *) ' Reading from input file inputlu.data'
            read (5,*)
            read (5,*)
            read (5,*) ipr, inorm
            read (5,*)
            read (5,*)
            read (5,*) itmax
            read (5,*)
            read (5,*)
            read (5,*) dt
            read (5,*)
            read (5,*)
            read (5,*) omega
            read (5,*)
            read (5,*)
            read (5,*) tolrsd(1),tolrsd(2),tolrsd(3),tolrsd(4),tolrsd(5)
            read (5,*)
            read (5,*)
            read (5,*) nx0, ny0, nz0
chomp            close(3)
chomp         else
chomp            ipr = ipr_default
chomp            inorm = inorm_default
chomp            itmax = itmax_default
chomp            dt = dt_default
chomp            omega = omega_default
chomp            tolrsd(1) = tolrsd1_def
chomp            tolrsd(2) = tolrsd2_def
chomp            tolrsd(3) = tolrsd3_def
chomp            tolrsd(4) = tolrsd4_def
chomp            tolrsd(5) = tolrsd5_def
chomp            nx0 = isiz1
chomp            ny0 = isiz2
chomp            nz0 = isiz3
chomp         endif

c---------------------------------------------------------------------
c   check problem size
c---------------------------------------------------------------------

         if ( ( nx0 .lt. 4 ) .or.
     >        ( ny0 .lt. 4 ) .or.
     >        ( nz0 .lt. 4 ) ) then

            write (*,2001)
 2001       format (5x,'PROBLEM SIZE IS TOO SMALL - ',
     >           /5x,'SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5')
            stop

         end if

         if ( ( nx0 .gt. isiz1 ) .or.
     >        ( ny0 .gt. isiz2 ) .or.
     >        ( nz0 .gt. isiz3 ) ) then

            write (*,2002)
 2002       format (5x,'PROBLEM SIZE IS TOO LARGE - ',
     >           /5x,'NX, NY AND NZ SHOULD BE EQUAL TO ',
     >           /5x,'ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY')
            stop

         end if


         write(*, 1001) nx0, ny0, nz0
         write(*, 1002) itmax


 1000 format(//,' Programming Baseline for NPB (PBN-O)',
     >          ' - LU Benchmark', /)
 1001    format(' Size: ', i3, 'x', i3, 'x', i3)
 1002    format(' Iterations: ', i3)
 1003    format(' Number of processes: ', i5, /)
         


      return
      end


