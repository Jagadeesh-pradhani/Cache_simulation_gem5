#include <sys/types.h>
#include <time.h>

double etime( arg )
double *arg;
{
    clock_t clock();           /* Returns Micro-Seconds. */
    register double foo;

    foo = (double)clock()/(double)CLOCKS_PER_SEC;
    return( foo );
}
double dsecond_( arg )
double *arg;
{
    clock_t clock();           /* Returns Micro-Seconds. */
    register double foo;

    foo = (double)clock()/(double)CLOCKS_PER_SEC;
    return( foo );
}
