/****************************************************************************
 *              Returns the total cpu time used in seconds.
 *      Emulates the Cray fortran library function of the same name.
 ****************************************************************************/
#ifdef LINUX            /* LINUX */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
double second( arg )
double *arg;
{
    struct rusage buf;
    double t1, t2;
    double t11, t12;
    double temp;

    getrusage( RUSAGE_SELF, &buf );

    /* Get system time and user time in SECONDS.*/
#ifdef DEBUG
    temp = (double)buf.ru_utime.tv_sec + (double)buf.ru_utime.tv_usec*1.0e-6 +
           (double)buf.ru_stime.tv_sec + (double)buf.ru_stime.tv_usec*1.0e-6;
    t11 = (double)buf.ru_utime.tv_sec;
    t12 = (double)(buf.ru_utime.tv_usec);
    printf("SECOND: t11 = %e t12 = %e\n", t11, t12 );
    t12 = (double)(buf.ru_utime.tv_usec)*1.0e-6;
    printf("SECOND: t11 = %e t12 = %e\n", t11, t12 );
#endif
    t1 = (double)buf.ru_utime.tv_sec + (double)(buf.ru_utime.tv_usec)*1.0e-6;
    t2 = (double)buf.ru_stime.tv_sec + (double)(buf.ru_stime.tv_usec)*1.0e-6;
    temp = t1 + t2;

    /* Return the sum of system and user time in SECONDS.*/
    return( temp );
}

double dsecnd( arg )
double *arg;
{
  return( second( arg ) ) ;
}

#elif AIX
#include        <unistd.h>
#include        <sys/time.h>
#include        <sys/times.h>
#include        <sys/resource.h>


/* Returns the current value of the wall clock timer.
 * C entry point.
 */
double
wc_second()
{
  struct timeval s_val;

  gettimeofday(&s_val,0);
  return ((double) s_val.tv_sec + 0.000001 * (double) (s_val.tv_usec));
}

/* Returns the current value of the user+system timer.  C entry point.
 */
double
us_second()
{
  struct        rusage  ru;
  double        tu, ts;

  getrusage(RUSAGE_SELF,&ru);

  tu = (double) (ru.ru_utime.tv_sec) +
         1.0e-6 * (double) (ru.ru_utime.tv_usec);
  ts = (double) ru.ru_stime.tv_sec +
         1.0e-6 * (double) (ru.ru_stime.tv_usec);

  return (tu + ts);
}
/* Returns the current value of the user+system timer.  C entry point.
 */
double
dsecnd( arg )
     double *arg;
{
  return( us_second() );
}
double
second( arg )
     double *arg;
{
  return( us_second() );
}

#elif _BGL
/*----------------------------------------------------------*/
/*    elapsed-time timing functions                         */
/*----------------------------------------------------------*/
extern unsigned long long	rts_get_timebase(void);

static double seconds_per_cycle = 1.4285714285714285714e-9; /* 700 MHz default */

#define WTIME(TB) TB = rts_get_timebase()
#define TCONV(TB1,TB2) seconds_per_cycle*((double) (TB1 - TB2))

double dsecnd(

double	*arg)

{
  static int			first = 1;
  static unsigned long long	tb0;
  unsigned long long		tb;

  if (first) {
    first = 0;
    WTIME(tb0);
  }

  return (TCONV(WTIME(tb), tb0));
}

double second(

double	*arg)

{
  static int			first = 1;
  static unsigned long long	tb0;
  unsigned long long		tb;

  if (first) {
    first = 0;
    WTIME(tb0);
  }

  return (TCONV(WTIME(tb), tb0));
}
#endif
