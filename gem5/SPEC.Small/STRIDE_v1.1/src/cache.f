c
c       This benchmark tests cache/memory speeds in keeping up
c       with vector saxpy memory references.  The parameter that
c       changes is the loop length in a single strided saxpy loop.
c       Times, MFLOPS and speed ratios (larger is worse) with vector
c       length 64 are printed out.
c       The timing for this package uses Cray second() function and a
c       sample 4.3BSD implementation is also provided.  If the system
c       clock is of poor quality this will show up as clumps of the MFLOP
c       and ratio data as plotted against the loop length.  If this occurs
c       simply increase NREP and rerun the benchmark.
c               Send Results to:
c               Dr. Mark K. Seager (seager@llnl.gov)
c               LLNL/Livermore Computing
c               PO BOX 808, L-554
c               Livermore, CA 94551-9911
c               925-423-3141-Voice
c               925-423-8715-Fax
c
      IMPLICIT REAL*8 (A-H, O-Z)
C        The following parameters yield a data size of about 1.0MB
C$$$      PARAMETER(M=256, N=256, M1=M+1,NNREP=102400, NTRIAL=128)
C        The following parameters yield a data size of about 4.0MB
      PARAMETER(M=512, N=512, M1=M+1,NNREP=1024000, NTRIAL=128)
C        The following parameters yield a data size of about 16.0MB
C$$$      PARAMETER(M=1024, N=1024, M1=M+1,NNREP=102400, NTRIAL=128)
      COMMON /ARRAYS/ X(M1*N), Y(M1*N)
      REAL TIC, TOC, TIR, TOR, TTC, TTR, ETIME, TARRAY(2), TMC, TMR
      DATA PI/3.141592653589793/
C
C        Set up the vector and scalar data
C
      DO 10 J = 1, N*M1
         T = 2. * PI * (J-1) / N
         X(J) = SIN(T)
         Y(J) = COS(T)
 10   CONTINUE
      T = 0.5
C
C        Run the reference loop length of 64 for base timing...
C
      TIC = ETIME( TARRAY )
C$$$      TIC = second( tarray )
      CALL WOTV( 64, T, X, Y, NNREP )
      TOC = ETIME( TARRAY )
C$$$      TOC = SECOND( TARRAY )
      TTC = TOC - TIC
      TMC = 2.0*1.0E-6*FLOAT(64)*FLOAT(NNREP)/TTC
      WRITE(6,1000)
      WRITE(6,1010) 64, TTC, TMC, 1.0
C
C       Run the experiment with different vector loop lengths.  We
C       exponentially increase the loop length to minimize output and
C       overall program runtime.  We also monitor the runtime of each
C       trial and decrease the number of repetitions so that the
C       runtime for each trial is about 2.0 seconds.
C
      BASE = EXP(ALOG(FLOAT(N*M))/FLOAT(NTRIAL))
      NREP = NNREP
      NLAST = 0
      DO 30 I = 1, NTRIAL
         NN = MAX( INT(BASE**I+1), NLAST+1 )
         NLAST = NN
         TIR = ETIME( TARRAY )
         CALL WOTV( NN, T, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR = 2.0*1.0E-6*FLOAT(NN)*FLOAT(NREP)/TTR
         WRITE(6,1010) NN, TTR, TMR, TMC/TMR
         IF( TTR.GT.2.0) THEN
            NREP = MAX(5,NREP/2)
CDEBUG            WRITE(6,*) 'NREP = ', NREP
         ENDIF
 30   CONTINUE
      CALL EXIT( 0 )
 1000 FORMAT(1X,'LENGTH',13X,'TIME',11X,'MFLOPS',12X,'RATIO')
 1010 FORMAT(1X,I6,3(1X,1PE16.7))
      END
      SUBROUTINE WOTV( N, T, X, Y, IREP )
C
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
