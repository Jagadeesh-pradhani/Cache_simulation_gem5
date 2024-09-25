c
c       This benchmark tests cache/memory speeds in keeping up
c       with vector saxpy memory references.  The parameter that
c       changes is the stride in the saxpy loop.  Times, MFLOPS
c       and speed ratios (larger is worse) with unit stride are
c       printed out.  The timing for this package uses Cray second()
c       function and a sample 4.3BSD implementation is also provided.
c       If the system clock is of poor quality this will show up
c       as clumps of the MFLOP and ratio data as plotted against the
c       stride.  If this occurs simply increase NREP and rerun the
c       benchmark.
c               Send Results to:
c               Dr. Mark K. Seager (seager@llnl.gov)
c               LLNL
c               PO BOX 808, L-60
c               Livermore, CA 94550
c               925-423-3141
c
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER(M=1024, N=16, M1=M+1,NREP=1024000)
      COMMON /FOOBAR/ X(M1*N), Y(M1*N)
      REAL TIC, TOC, TIR, TOR, TTC, TTR, ETIME, TARRAY(2), TMC, TMR
      DATA PI/3.141592653589793/, IZERO /0/
C
      DO 10 J = 1, N*M1
         T = 2. * PI * (J-1) / N
         X(J) = SIN(T)
         Y(J) = COS(T)
 10   CONTINUE
      T = 1.0/T
      TIC = ETIME( TARRAY )
C$$$      TIC = second( tarray )
      CALL WOT( 64, T, X, 1, Y, NREP, DOT )
      TOC = ETIME( TARRAY )
C$$$      TOC = SECOND( TARRAY )
      TTC = TOC - TIC
      TMC = 2.0*1.0E-6*FLOAT(64)*FLOAT(NREP)/TTC
      WRITE(6,1000)
      WRITE(6,1010) 64, TTC, TMC, 1.0
 1000 FORMAT(1X,'#LENGTH',13X,'TIME',11X,'MFLOPS',12X,'RATIO')
 1010 FORMAT(2X,I6,3(1X,1PE16.7))
      DO 20 LEN = 8, N*M, 8
         IREP = MAX( 2, 64*NREP/LEN )
         TIR = ETIME( TARRAY )
C$$$         TIR = SECOND( TARRAY )
         CALL WOT( LEN, T, X, INCX, Y, IREP, DOT )
         TOR = ETIME( TARRAY )
C$$$         TOR = SECOND( TARRAY )
         TTR = TOR - TIR
         TMR = 2.0*1.0E-6*FLOAT(LEN)*FLOAT(IREP)/TTR
         WRITE(6,1010) LEN, TTR, TMR,
     $        (64.0*FLOAT(NREP)/(FLOAT(LEN*IREP)))*TTR/TTC
 20   continue
      CALL EXIT( IZERO )
      END
      SUBROUTINE WOT( N, t, x, incx, y, irep, dot )
C
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(N), Y(N)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 J = 1, IREP
CCDIR$ VECTOR
         DO 100 I = 1, n
            DOT = DOT + X(I)*Y(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE FOOI( N, t, x, incx, y, irep, dot )
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(*), Y(*)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 j = 1, irep
CCDIR$ VECTOR
         DO 100 i = 1, N*incx,incx
            DOT = DOT + X(I)*Y(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
