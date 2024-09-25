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
               LLNL/Integrated Computing & Communications Department
               PO BOX 808, L-554
               Livermore, CA 94551-9900
               925-423-3141-Voice
               925-423-8911-Fax
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M       1024
#define N       256
#define M1      (M+1)
#define NREP    102400
#define PI      3.141592653589793

double x[M1*N], y[M1*N];

main( argc, argv )
     int argc;
     char *argv[];
{
  double dsecnd( double *), dummy;
  double t, tic, toc, tir, tor, ttc, ttr, tmc, tmr;
  int j, incx;
  void wotv1( int, double, double *, int, double *, int );
  void wotvi( int, double, double *, int, double *, int );

  for( j=0; j<N*M1; j++ ) {
    t = 2. * PI * (double)j / (double)N;
    x[j] = sin(t);
    y[j] = cos(t);
  }
  t = 0.5;
  tic = dsecnd( &dummy );
  wotv1( N, t, x, 1, y, NREP );
  toc = dsecnd( &dummy );
  ttc = toc - tic;
  tmc = 2.0*1.0e-6*(double)N*(double)NREP/ttc;
  printf(" STRIDE             TIME           MFLOPS            RATIO\n");
  printf(" %6d %16.7e %16.7e %16.7e\n", 1, ttc, tmc, 1.0);
  for( incx = 1; incx <= M; incx++ ) {
    tir = dsecnd( &dummy );
    wotvi( N, t, x, incx, y, NREP );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    tmr = 2.0*1.0e-6*(double)N*(double)NREP/ttr;
    printf(" %6d %16.7e %16.7e %16.7e\n", incx, ttr, tmr, ttr/ttc);
  }
  exit( 0 );
}

void wotv1( n, t, x, incx, y, irep )
     int n, incx, irep;
     double t, *x, *y;
{
  register int i, j;
  register double *xx = x, *yy = y;
  /*
    Yo! Make sure the j loop is not optimized away!!!
  */

  for( j=0; j<irep; j++ ) {
    t = 1.0/(double)(j+1);
#pragma ivdep
#pragma swp
#pragma vector
    for( i=0; i<n; i++ ) {
      yy[i] += t*xx[i];
    }
  }
}

void wotvi( n, t, x, incx, y, irep )
     int n, incx, irep;
     double t, *x, *y;
{
  register int i, j;
  double *xx=x, *yy=y;

  /*
      Yo! Make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t = 1.0/(double)(j+1);
#pragma vector
#pragma ivdep
#pragma swp
    for( i=0; i<n*incx; i+=incx ) {
      yy[i] += t*xx[i];
    }
  }
}
