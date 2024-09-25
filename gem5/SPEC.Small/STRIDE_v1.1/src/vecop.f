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
c               Lawrence Livermore National Laboratory
c               PO BOX 808, L-554
c               Livermore, CA 94550
c               925-423-3141-Voice
c               925-423-6961-Fax
c
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER(N=1024, N1=N+1, NNREP=1024000)
      COMMON /ARRAYS/ X(N1), Y(N1), Z(N1), IDX(N1)
      REAL ETIME, TARRAY(2)
      REAL TIR, TOR, TTR, TMR2, TMR3, TMR4, TMR5, TMR6, TMR7
      DATA PI/3.141592653589793/
C
      DO 10 J = 1, N1
         T = 2. * PI * (J-1) / N
         X(J) = SIN(T)
         Y(J) = COS(T)
         Z(J) = X(J)
         IDX(J) = N1-J
 10   CONTINUE
      T1 = 1.0D-6
      T2 = 1.0D-6
      T3 = 1.0D-6
      NREP = NNREP
      LEN = 4
C
         TIR = ETIME( TARRAY )
         CALL V1S1M3( LEN, T1, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR2 = 2.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S2M3( LEN, T1, T2, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR3 = 3.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S3M3( LEN, T1, T2, T3, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR4 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V2S2M3( LEN, T1, T2, T3, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR5 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V2S2M4( LEN, T1, T2, T3, X, Y, Z, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR6 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S1I3( LEN, T1, X, Y, IDX, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR7 = 2.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
      WRITE(6,1000)
      WRITE(6,1010) LEN, TMR2, TMR3, TMR4, TMR5, TMR6, TMR7
      DO 20 LEN = 4, N, 4
C
         TIR = ETIME( TARRAY )
         CALL V1S1M3( LEN, T1, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR2 = 2.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S2M3( LEN, T1, T2, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR3 = 3.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S3M3( LEN, T1, T2, T3, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR4 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V2S2M3( LEN, T1, T2, T3, X, Y, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR5 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V2S2M4( LEN, T1, T2, T3, X, Y, Z, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR6 = 4.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         TIR = ETIME( TARRAY )
         CALL V1S1I3( LEN, T1, X, Y, IDX, NREP )
         TOR = ETIME( TARRAY )
         TTR = TOR - TIR
         TMR7 = 2.0*1.0E-6*FLOAT(LEN)*FLOAT(NREP)/TTR
C
         WRITE(6,1010) LEN, TMR2, TMR3, TMR4, TMR5, TMR6, TMR7
C
         if( TTR.gt.2.0) THEN
            NREP = max(1,NREP/2)
C$$$            write(6,*) 'NREP = ', NREP
         endif
 20   CONTINUE
      STOP "ALL DONE."
 1000 FORMAT(1X,' LEN',4X,'V1S1M3',4X,'V1S2M3',4X,'V1S3M3',
     $     4X,'V2S2M3',4X,'V2S2M4',4X,'V1S1I3')
 1010 FORMAT(1X,I4,6(1X,1PE9.3))
      END
      SUBROUTINE V1S1M3( N, T1, X, Y, IREP )
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
            Y(I) = Y(I) + T1*X(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE V1S1I3( N, T1, X, Y, IDX, IREP )
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(*), Y(*), IDX(*)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 J = 1, IREP
         T = 1.0/DBLE(J)
CCDIR$ VECTOR
CCDIR$ IVDEP
         DO 100 I = 1, N
            Y(IDX(I)) = Y(IDX(I)) + T1*X(IDX(I))
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE V1S2M3( N, T1, T2, X, Y, IREP )
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
            Y(I) = T1*Y(I) + T2*X(I)
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE V1S3M3( N, T1, T2, T3, X, Y, IREP )
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
            Y(I) = T1*(T2*Y(I) + T3*X(I))
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE V2S2M3( N, T1, T2, T3, X, Y, IREP )
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
            Y(I) = Y(I) + T1*(Y(I) + T3*X(I))
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
      SUBROUTINE V2S2M4( N, T1, T2, T3, X, Y, Z, IREP )
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION X(*), Y(*), Z(*)
C
C       Yo! Make sure the j loop is not optimized away!!!
CCDIR$ SCALAR
      DO 110 J = 1, IREP
         T = 1.0/DBLE(J)
CCDIR$ VECTOR
CCDIR$ SWP
         DO 100 I = 1, N
            Y(I) = Y(I) + T1*(Z(I) + T3*X(I))
 100     CONTINUE
 110  CONTINUE
      RETURN
      END
