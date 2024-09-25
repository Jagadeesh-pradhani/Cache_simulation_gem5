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
c               LLNL/Integrated Computing & Communications Department
c               PO BOX 808, L-554
c               Livermore, CA 94551-9900
c               925-423-3141-Voice
c               925-423-8911-Fax
c
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER(M=1024, N=256, M1=M+1,NREP=102400)
      COMMON /ARRAYS/  X(M1*N), Y(M1*N)
      REAL TIC, TOC, TIR, TOR, TTC, TTR, ETIME, TARRAY(2), TMC, TMR
      DATA PI/3.141592653589793/
C
      DO 10 J = 1, N*M1
         T = 2. * PI * (J-1) / N
         X(J) = SIN(T)
         Y(J) = COS(T)
 10   CONTINUE
      T = 0.5
      TIC = ETIME( TARRAY )
C$$$      tic = second( tarray )
      CALL WOTV1( N, T, X, 1, Y, NREP )
      TOC = ETIME( TARRAY )
C$$$      toc = second( tarray )
      TTC = TOC - TIC
      TMC = 2.0*1.0E-6*FLOAT(N)*FLOAT(NREP)/TTC
      WRITE(6,1000)
      WRITE(6,1010) 1, TTC, TMC, 1.0
      DO 20 INCX = 1, M
         TIR = ETIME( TARRAY )
C$$$         TIR = SECOND( TARRAY )
         CALL WOTVI( N, T, X, INCX, Y, NREP )
         TOR = ETIME( TARRAY )
C$$$         TOR = SECOND( TARRAY )
         TTR = TOR - TIR
         TMR = 2.0*1.0E-6*FLOAT(N)*FLOAT(NREP)/TTR
         WRITE(6,1010) INCX, TTR, TMR, TTR/TTC
 20   CONTINUE
      CALL EXIT( 0 )
 1000 FORMAT(1X,'STRIDE',13X,'TIME',11X,'MFLOPS',12X,'RATIO')
 1010 FORMAT(1X,I6,3(1X,1PE16.7))
      END
      SUBROUTINE WOTV1( N, T, X, INCX, Y, IREP )
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(*), Y(*)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 J = 1, IREP
         T = 1.0/DBLE(J)
CCDIR$ VECTOR
CCDIR$ SWP
         DO 100 I = 1, N
            Y(I) = Y(I) + T*X(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE WOTVI( N, T, X, INCX, Y, IREP )
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(*), Y(*)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 J = 1, IREP
         T = 1.0/DBLE(J)
CCDIR$ IVDEP
CCDIR$ SWP
         DO 100 I = 1, N*INCX, INCX
            Y(I) = Y(I) + T*X(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
