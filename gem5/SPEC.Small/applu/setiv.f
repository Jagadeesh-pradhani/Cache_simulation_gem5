
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine setiv

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   set the initial values of independent variables based on tri-linear
c   interpolation of boundary values in the computational space.
c
c---------------------------------------------------------------------

      implicit none

      include 'applu.incl'

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, m
      double precision  xi, eta, zeta
      double precision  pxi, peta, pzeta
      double precision  ue_1jk(5),ue_nx0jk(5),ue_i1k(5),
     >        ue_iny0k(5),ue_ij1(5),ue_ijnz(5)


!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(PZETA,PETA,PXI,M,UE_IJNZ,
!$OMP& UE_IJ1,UE_INY0K,UE_I1K,UE_NX0JK,UE_1JK,XI,I,ETA,J,ZETA,K)
!$OMP& SHARED(NX0,NY0,NZ)
      do k = 2, nz - 1
         zeta = ( dble (k-1) ) / (nz-1)
         do j = 2, ny - 1
            eta = ( dble (j-1) ) / (ny0-1)
            do i = 2, nx - 1
               xi = ( dble (i-1) ) / (nx0-1)
               call exact (1,j,k,ue_1jk)
               call exact (nx0,j,k,ue_nx0jk)
               call exact (i,1,k,ue_i1k)
               call exact (i,ny0,k,ue_iny0k)
               call exact (i,j,1,ue_ij1)
               call exact (i,j,nz,ue_ijnz)
               do m = 1, 5
                  pxi =   ( 1.0d+00 - xi ) * ue_1jk(m)
     >                              + xi   * ue_nx0jk(m)
                  peta =  ( 1.0d+00 - eta ) * ue_i1k(m)
     >                              + eta   * ue_iny0k(m)
                  pzeta = ( 1.0d+00 - zeta ) * ue_ij1(m)
     >                              + zeta   * ue_ijnz(m)

                  u( m, i, j, k ) = pxi + peta + pzeta
     >                 - pxi * peta - peta * pzeta - pzeta * pxi
     >                 + pxi * peta * pzeta

               end do
            end do
         end do
      end do
!$OMP END PARALLEL DO

      return
      end
