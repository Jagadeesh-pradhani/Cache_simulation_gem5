
      subroutine print_results(name, class, n1, n2, n3, niter, 
     >               verified, npbversion)
      
      implicit none
      character*2 name
      character*1 class
      integer n1, n2, n3, niter, j
      character size*13
      logical verified
      character*(*) npbversion 

         write (*, 2) name 
 2       format(//, ' ', A2, ' Benchmark Completed.')

         write (*, 3) Class
 3       format(' Class           = ', 12x, a12)

c   If this is not a grid-based problem (EP, FT, CG), then
c   we only print n1, which contains some measure of the
c   problem size. In that case, n2 and n3 are both zero.
c   Otherwise, we print the grid size n1xn2xn3

         if ((n2 .eq. 0) .and. (n3 .eq. 0)) then
            if (name(1:2) .eq. 'EP') then
               write(size, '(f12.0)' ) 2.d0**n1
               do j =13,1,-1
                  if (size(j:j) .eq. '.') size(j:j) = ' '
               end do
               write (*,42) size
 42            format(' Size            = ',12x, a14)
            else
               write (*,44) n1
 44            format(' Size            = ',12x, i12)
            endif
         else
            write (*, 4) n1,n2,n3
 4          format(' Size            =  ',12x, i3,'x',i3,'x',i3)
         endif

         write (*, 5) niter
 5       format(' Iterations      = ', 12x, i12)
         
CHOMP         if (verified) then 
CHOMP            write(*,12) '  SUCCESSFUL'
CHOMP         else
CHOMP            write(*,12) 'UNSUCCESSFUL'
CHOMP         endif
CHOMP 12      format(' Verification    = ', 12x, a)
CHOMP
CHOMP         write(*,13) npbversion
CHOMP 13      format(' Version         = ', 12x, a12)


CHOMP         write (*,130)
CHOMP 130     format(//' Please send all errors/feedbacks to:'//
CHOMP     >            ' PBN Working Team '/
CHOMP     >            ' pbn@nas.nasa.gov'//)
c 130     format(//' Please send the results of this run to:'//
c     >            ' NPB Development Team '/
c     >            ' Internet: npb@nas.nasa.gov'/
c     >            ' '/
c     >            ' If email is not available, send this to:'//
c     >            ' MS T27A-1'/
c     >            ' NASA Ames Research Center'/
c     >            ' Moffett Field, CA  94035-1000'//
c     >            ' Fax: 415-604-3957'//)


         return
         end

