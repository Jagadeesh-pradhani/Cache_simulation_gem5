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
               Lawrence Livermore National Laboratory
               PO BOX 808, L-554
               Livermore, CA 94551-9900
               925-423-3141-Voice
               925-423-8911-Fax

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M       4096
#define N       1024
#define N1      (N+1)
#define NNREP   1024000
#define PI      3.141592653589793
#define IZERO   0
#define MAX(x,y) (x<y?y:x)

double x[M*N1], y[M*N1], z[M*N1];
int idx[N1];
main(argc, argv)
     int argc;
     char *argv[];
{
  int j, nrep, len;
  register double t;
  double t1, t2, t3;
  double tir, tor, ttr, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7, ttmax;
  double sin(double), cos(double);
  double dsecnd( double *), dummy;
  void v1s1m3(), v1s2m3(), v1s3m3(), v2s2m3(), v2s2m4(), v1s1i3();

  printf("Starting Setup...\n");
  for( j = 0; j<N1; j++ ) {
    t = 2. * PI * j / N;
    x[j] = sin(t);
    y[j] = cos(t);
    z[j] = x[j];
    idx[j] = j*4 % N1;
  }
  printf("Setup Complete...\n");
  t1 = 1./t;
  t2 = 1./t;
  t3 = 1./t;
  nrep = NNREP;
  len = 4;

  tir = dsecnd( &dummy );
  v1s1m3( len, t1, x, y, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr2 = 2.0*1.0e-6*(double)len*(double)nrep/ttr;

  tir = dsecnd( &dummy );
  v1s2m3( len, t1, t2, x, y, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr3 = 3.0*1.0e-6*(double)len*(double)nrep/ttr;

  tir = dsecnd( &dummy );
  v1s3m3( len, t1, t2, t3, x, y, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr4 = 4.0*1.0e-6*(double)len*(double)nrep/ttr;

  tir = dsecnd( &dummy );
  v2s2m3( len, t1, t2, t3, x, y, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr5 = 4.0*1.0e-6*(double)len*(double)nrep/ttr;

  tir = dsecnd( &dummy );
  v2s2m4( len, t1, t2, x, y, z, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr6 = 4.0*1.0e-6*(double)len*(double)nrep/ttr;

  tir = dsecnd( &dummy );
  v1s1i3( len, t1, x, y, idx, nrep );
  tor = dsecnd( &dummy );
  ttr = tor - tir;
  tmr7 = 2.0*1.0e-6*(double)len*(double)nrep/ttr;

  printf(" LEN   V1S1M3      V1S2M3      V1S3M3      V2S2M3      V2S2M4      V1S1I3\n");
  printf("%4d %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e\n",
       len, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7 );

  for( len=4; len<N; len+=4 ) {
    ttmax = 1.0;

    tir = dsecnd( &dummy );
    v1s1m3( len, t1, x, y, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr2 = 2.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    tir = dsecnd( &dummy );
    v1s2m3( len, t1, t2, x, y, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr3 = 3.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    tir = dsecnd( &dummy );
    v1s3m3( len, t1, t2, t3, x, y, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr4 = 4.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    tir = dsecnd( &dummy );
    v2s2m3( len, t1, t2, t3, x, y, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr5 = 4.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    tir = dsecnd( &dummy );
    v2s2m4( len, t1, t2, x, y, z, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr6 = 4.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    tir = dsecnd( &dummy );
    v1s1i3( len, t1, x, y, idx, nrep );
    tor = dsecnd( &dummy );
    ttr = tor - tir;
    ttmax = MAX( ttmax, ttr );
    tmr7 = 2.0*1.0e-6*(double)(len)*(double)(nrep)/ttr;

    printf("%4d %9.5e %9.5e %9.5e %9.5e %9.5e %9.5e\n",
          len, tmr2, tmr3, tmr4, tmr5, tmr6, tmr7 );

    if( ttmax>2.0) {
      nrep = MAX(1,nrep/2);
      printf( "TTMAX=%e, NREP = %d\n", ttmax, nrep );
    }
  }
  exit( IZERO );
}
void v1s1m3( n, t1, x, y, irep )
     int n, irep;
     double t1, *x, *y;
{
  register int i, j;
  register double *xx = x, *yy = y;

  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
    for( i=0; i<n; i++ ) {
      yy[i] += t1*xx[i];
    }
  }
}
void v1s1i3( n, t1, x, y, idx, irep )
     int n, irep, *idx;
     double t1, *x, *y;
{
  register int i, j;
  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
#pragma ivdep
    for( i=0; i<n; i++ ) {
      y[idx[i]] += t1*x[idx[i]];
    }
  }
}
void v1s2m3( n, t1, t2, x, y, irep )
     int n, irep;
     double t1, t2, *x, *y;
{
  register int i, j;
  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
    for( i=1; i<n; i++ ) {
      y[i] = t1*y[i] + t2*x[i];
    }
  }
}

void v1s3m3( n, t1, t2, t3, x, y, irep )
     int n, irep;
     double t1, t2, t3, *x, *y;
{
  register int i, j;
  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
    for(i=0; i<n; i++ ) {
      y[i] = t1*(t2*y[i] + t3*x[i]);
    }
  }
}
void v2s2m3( n, t1, t2, t3, x, y, irep )
     int n, irep;
     double t1, t2, t3, *x, *y;
{
  register int i, j;
  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
    for( i=0; i<n; i++ ) {
      y[i] += t1*(y[i] + t3*x[i]);
    }
  }
}
void v2s2m4( n, t1, t2, x, y, z, irep )
     int n, irep;
     double t1, t2, *x, *y, *z;
{
  register int i, j;
  /*
    yo! make sure the j loop is not optimized away!!!
  */
  for( j=0; j<irep; j++ ) {
    t1 = 1.0/(double)(j+1);
#pragma vector
    for( i=0; i<n; i++ ) {
      y[i] += t1*(z[i] + t2*x[i]);
    }
  }
}
