/*
       This benchmark tests cache/memory speeds in keeping up
       with vector saxpy memory references.  The parameter that
       changes is the stride in the saxpy loop.  Times, MFLOPS
       and speed ratios (larger is worse) with unit stride are
       printed out.  The timing for this package uses Cray second()
        function and a sample 4.3BSD implementation is also provided.
        If the system clock is of poor quality this will show up
       as clumps of the MFLOP and ratio data as plotted against the
       stride.  If this occurs simply increase NREP and rerun the
       benchmark.
               Send Results to:
               Dr. Mark K. Seager (seager@llnl.gov)
               LLNL
               PO BOX 808, L-60
               Livermore, CA 94550
               925-423-3141
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(x,y) (x<y?y:x)
#define M       4096
#define M1      (M+1)
#define N       64
#define NREP    1024000
#define PI      3.141592654
#define IZERO   0

double x[M1*N], y[M1*N];
double tic, toc, tir, tor, ttc, ttr, tmc, tmr;

main( argc, argv )
     int argc;
     char *argv[];
{
  int j;
  int len, irep;
  register double t;
  double second(), dummy;
  void wot();

  for( j=0; j<N*M1; j++ ) {
    t = 2. * PI * (j-1) / N;
    x[j] = sin(t);
    y[j] = cos(t);
  }
  tic = second( &dummy );
  wot( 64, t, NREP );
  toc = second( &dummy );
  ttc = toc - tic;
  tmc = 2.0*1.0e-6*(double)64*(double)NREP/ttc;
  printf(" #LENGTH             TIME           MFLOPS            RATIO\n");
  printf("  %6d %16e %16e %16e\n", 64, ttc, tmc, 1.0);

  for( len=8; len<= N*M; len+=N ) {
    irep = MAX( 2, 64*NREP/len );
    tir = second( &dummy );
    wot( len, t, irep );
    tor = second( &dummy );
    ttr = tor - tir;
    tmr = 2.0*1.0e-6*(double)len*(double)irep/ttr;
    printf("  %6d %16e %16e %16e\n", len, ttr, tmr,
          (64.0*(double)NREP/((double)(len*irep)))*ttr/ttc );
  }
  exit( IZERO );
}
void wot( n, t, irep )
     int n, irep;
     double t;
{
  int i, j;

  /*
    Yo! Make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t = 1.0/(double)(j+1);
#pragma vector
    for( i=0; i<n; i++ ) {
      y[i] += t*x[i];
    }
  }
}
