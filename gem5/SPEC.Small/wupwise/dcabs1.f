C
C A SPEC-modified DCABS1 which is intended to avoid problems with
C compilers that object to equivalences in the original.
C
C      double precision function dcabs1(z)
C      double complex z,zz
C      double precision t(2)
C      equivalence (zz,t(1))
C      zz = z
C      dcabs1 = dabs(t(1)) + dabs(t(2))
C      return
C      end

      REAL*8 FUNCTION DCABS1 (QVAR)
      REAL*8 QVAR(2)
      DCABS1 = DABS(QVAR(1)) + DABS(QVAR(2)) 
      RETURN
      END
