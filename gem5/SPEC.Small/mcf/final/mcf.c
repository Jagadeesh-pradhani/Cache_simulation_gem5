/**************************************************************************
DEFINES.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 21:46:48 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: defines.h,v 1.13 2005/02/17 21:43:12 bzfloebe Exp $  */


#ifndef _PROTOTYP_H
#define _PROTOTYP_H


#ifndef _PROTO_
#if    defined(__STDC__) || defined(__cplusplus) \
    || defined(WANT_STDC_PROTO) || defined(SPEC_CPU)
#define _PROTO_( args ) args
#else
#define _PROTO_( args ) 
#endif
#endif

#endif

#ifndef _DEFINES_H
#define _DEFINES_H

#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <assert.h>

#ifdef INTERNAL_TIMING
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#endif


/*#include "prototyp.h"*/


#define UNBOUNDED        1000000000
#define ZERO             0
#define MAX_ART_COST     (long)(100000000L)
#define ARITHMETIC_TYPE "I"



#define FIXED           -1
#define BASIC            0
#define AT_LOWER         1
#define AT_UPPER         2
/* #define AT_ZERO           3  NOT ALLOWED FOR THE SPEC VERSION */
#undef AT_ZERO


#define UP    1
#define DOWN  0



typedef long flow_t;
typedef long cost_t;




#ifndef NULL
#define NULL 0
#endif


#ifndef ABS
#define ABS( x ) ( ((x) >= 0) ? ( x ) : -( x ) )
#endif


#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif


#ifndef SET_ZERO
#define SET_ZERO( vec, n ) if( vec ) memset( (void *)vec, 0, (size_t)n )
#endif


#ifndef FREE
#define FREE( vec ) if( vec ) free( (void *)vec )
#endif


typedef struct node node_t;
typedef struct node *node_p;

typedef struct arc arc_t;
typedef struct arc *arc_p;



struct node
{
  cost_t potential; 
  int orientation;
  node_p child;
  node_p pred;
  node_p sibling;
  node_p sibling_prev;     
  arc_p basic_arc; 
  arc_p firstout, firstin;
  arc_p arc_tmp;
  flow_t flow;
  long depth; 
  int number;
  int time;
};



struct arc
{
  cost_t cost;
  node_p tail, head;
  int ident;
  arc_p nextout, nextin;
  flow_t flow;
  cost_t org_cost;
};



typedef struct network
{
  char inputfile[200];
  char clustfile[200];
  long n, n_trips;
  long max_m, m, m_org, m_impl;
  long max_residual_new_m, max_new_m;
  
  long primal_unbounded;
  long dual_unbounded;
  long perturbed;
  long feasible;
  long eps;
  long opt_tol;
  long feas_tol;
  long pert_val;
  long bigM;
  double optcost;  
  cost_t ignore_impl;
  node_p nodes, stop_nodes;
  arc_p arcs, stop_arcs;
  arc_p dummy_arcs, stop_dummy; 
  long iterations;
  long bound_exchanges;
  long checksum;
} network_t;


#endif
/**************************************************************************
IMPLICIT.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:21:18 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: implicit.h,v 1.11 2005/02/17 19:42:21 bzfloebe Exp $  */


#ifndef _IMPLICIT_H
#define _IMPLICIT_H


/*#include "mcfutil.h"
#include "mcflimit.h"*/


long price_out_impl _PROTO_(( network_t * ));
long suspend_impl _PROTO_(( network_t *, cost_t, long ));


#endif
/**************************************************************************
MCF.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:21:29 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: mcf.h,v 1.9 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _MCF_H
#define _MCF_H

/*
#include "defines.h"
#include "mcfutil.h"
#include "readmin.h"
#include "output.h"  
#include "pstart.h"
#include "psimplex.h"
#include "pbeampp.h"
#include "implicit.h"
#include "limits.h"

*/
#endif
/**************************************************************************
MCFLIMIT.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 22:24:36 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: mcflimit.h,v 1.12 2005/02/17 21:43:12 bzfloebe Exp $  */


#ifndef _MCF_LIMITS_H
#define _MCF_LIMITS_H


#define BIGM 1.0e7
#define STRECHT(x) ((long)(1.25 * (double)(x)))

#define MAX_NB_TRIPS_FOR_SMALL_NET 15000

#define MAX_NEW_ARCS_SMALL_NET 3000000
#define MAX_NEW_ARCS_LARGE_NET 28900000

#define MAX_NB_ITERATIONS_SMALL_NET  5
#define MAX_NB_ITERATIONS_LARGE_NET  5


/*
// Some operating systems and compiler, respectively, do not handle reallocs
// properly. Thus, this program requires a somehow static memory handling
// without reallocation of the main (and quite huge) arc array.
*/
#define SPEC_STATIC


#endif
/**************************************************************************
MCFUTIL.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:21:49 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: mcfutil.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _MCFUTIL_H
#define _MCFUTIL_H


/*#include "defines.h"*/


extern void refresh_neighbour_lists _PROTO_(( network_t * ));
extern long refresh_potential _PROTO_(( network_t * ));
extern double flow_cost _PROTO_(( network_t * ));
extern double flow_org_cost _PROTO_(( network_t * ));
extern long primal_feasible _PROTO_(( network_t * ));
extern long dual_feasible _PROTO_(( network_t * ));
extern long getfree _PROTO_(( network_t * ));


#endif
/**************************************************************************
OUTPUT.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:21:59 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: output.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _OUTPUT_H
#define _OUTPUT_H


/*#include "mcfutil.h"*/


extern long write_circulations _PROTO_(( char *, network_t * ));


#endif
/**************************************************************************
PBEAMPP.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:09 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbeampp.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _PBEAMPP_H
#define _PBEAMPP_H


/*#include "defines.h"*/


extern arc_t *primal_bea_mpp _PROTO_(( long, arc_t*, arc_t*, cost_t* ));


#endif
/**************************************************************************
PBLA.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:19 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbla.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _PBLA_H
#define _PBLA_H


/*#include "defines.h"*/


extern node_t *primal_iminus _PROTO_(( flow_t *, long *, node_t *, 
                                       node_t *, node_t ** ));


#endif
/**************************************************************************
PFLOWUP.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:28 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pflowup.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _PFLOWUP_H
#define _PFLOWUP_H


/*#include "defines.h"*/


extern void primal_update_flow _PROTO_(( node_t *, node_t *, node_t * ));


#endif
/**************************************************************************
PROTOTYPE.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Wed Feb 16 20:24:10 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: prototyp.h,v 1.11 2005/02/17 19:42:21 bzfloebe Exp $  */










/**************************************************************************
PSIMPLEX.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:48 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: psimplex.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _PSIMPLEX_H
#define _PSIMPLEX_H

/*
#include "defines.h"
#include "pbeampp.h"
#include "pbla.h"
#include "pflowup.h"
#include "treeup.h"
#include "mcfutil.h"
*/

extern long primal_net_simplex _PROTO_(( network_t * ));


#endif
/**************************************************************************
PSTART.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:58 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pstart.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _PSTART_H
#define _PSTART_H


/*#include "defines.h"*/


extern long primal_start_artificial _PROTO_(( network_t * ));


#endif
/**************************************************************************
READMIN.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:23:07 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: readmin.h,v 1.11 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _READMIN_H
#define _READMIN_H

/*
#include "defines.h"
#include "mcfutil.h"
#include "mcflimit.h"
*/

extern long read_min _PROTO_(( network_t * ));


#endif
/**************************************************************************
TREEUP.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:23:16 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: treeup.h,v 1.10 2005/02/17 19:42:21 bzfloebe Exp $  */



#ifndef _TREEUP_H
#define _TREEUP_H


/*#include "defines.h"*/


extern void update_tree _PROTO_(( long, long, flow_t, flow_t, node_t *, 
                                  node_t *, node_t *, node_t *, node_t *, 
                                  arc_t *, cost_t, flow_t ));


#endif
/**************************************************************************
IMPLICIT.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 21:45:45 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: implicit.c,v 1.20 2005/02/17 21:43:12 bzfloebe Exp $  */



/*#include "implicit.h"*/



#ifdef _PROTO_
long resize_prob( network_t *net )
#else
long resize_prob( net )
     network_t *net;
#endif
{
    arc_t *arc;
    node_t *node, *stop, *root;
    size_t off;
            
    
    assert( net->max_new_m >= 3 );


    net->max_m += net->max_new_m;
    net->max_residual_new_m += net->max_new_m;
    

#if defined AT_HOME
    printf( "\nresize arcs to %4ld MB (%ld elements a %d B)\n\n",
            net->max_m * sizeof(arc_t) / 0x100000,
            net->max_m,
            sizeof(arc_t) );
    fflush( stdout );
#endif


    arc = (arc_t *) realloc( net->arcs, net->max_m * sizeof(arc_t) );
    if( !arc )
    {
        printf( "network %s: not enough memory\n", net->inputfile );
        fflush( stdout );
        return -1;
    }
    
    off = (size_t)arc - (size_t)net->arcs;
        
    net->arcs = arc;
    net->stop_arcs = arc + net->m;

    root = node = net->nodes;
    for( node++, stop = (void *)net->stop_nodes; node < stop; node++ )
        if( node->pred != root )
            node->basic_arc = (arc_t *)((size_t)node->basic_arc + off);
        
    return 0;
}







#ifdef _PROTO_
void insert_new_arc( arc_t *new, long newpos, node_t *tail, node_t *head,
                     cost_t cost, cost_t red_cost )
#else
void insert_new_arc( new, newpos, tail, head, cost, red_cost )
     arc_t *new;
     long newpos;
     node_t *tail;
     node_t *head;
     cost_t cost;
     cost_t red_cost;
#endif
{
    long pos;

    new[newpos].tail      = tail;
    new[newpos].head      = head;
    new[newpos].org_cost  = cost;
    new[newpos].cost      = cost;
    new[newpos].flow      = (flow_t)red_cost; 
    
    pos = newpos+1;
    while( pos-1 && red_cost > (cost_t)new[pos/2-1].flow )
    {
        new[pos-1].tail     = new[pos/2-1].tail;
        new[pos-1].head     = new[pos/2-1].head;
        new[pos-1].cost     = new[pos/2-1].cost;
        new[pos-1].org_cost = new[pos/2-1].cost;
        new[pos-1].flow     = new[pos/2-1].flow;
        
        pos = pos/2;
        new[pos-1].tail     = tail;
        new[pos-1].head     = head;
        new[pos-1].cost     = cost;
        new[pos-1].org_cost = cost;
        new[pos-1].flow     = (flow_t)red_cost; 
    }
    
    return;
}   






#ifdef _PROTO_
void replace_weaker_arc( network_t *net, arc_t *new, node_t *tail, node_t *head,
                         cost_t cost, cost_t red_cost )
#else
void replace_weaker_arc( net, new, tail, head, cost, red_cost )
     network *net;
     arc_t *new;
     node_t *tail;
     node_t *head;
     cost_t cost;
     cost_t red_cost;
#endif
{
    long pos;
    long cmp;

    new[0].tail     = tail;
    new[0].head     = head;
    new[0].org_cost = cost;
    new[0].cost     = cost;
    new[0].flow     = (flow_t)red_cost; 
                    
    pos = 1;
    cmp = (new[1].flow > new[2].flow) ? 2 : 3;
    while( cmp <= net->max_residual_new_m && red_cost < new[cmp-1].flow )
    {
        new[pos-1].tail = new[cmp-1].tail;
        new[pos-1].head = new[cmp-1].head;
        new[pos-1].cost = new[cmp-1].cost;
        new[pos-1].org_cost = new[cmp-1].cost;
        new[pos-1].flow = new[cmp-1].flow;
        
        new[cmp-1].tail = tail;
        new[cmp-1].head = head;
        new[cmp-1].cost = cost;
        new[cmp-1].org_cost = cost;
        new[cmp-1].flow = (flow_t)red_cost; 
        pos = cmp;
        cmp *= 2;
        if( cmp + 1 <= net->max_residual_new_m )
            if( new[cmp-1].flow < new[cmp].flow )
                cmp++;
    }
    
    return;
}   




#if defined AT_HOME
#include <sys/time.h>
double Get_Time( void  ) 
{
    struct timeval tp;
    struct timezone tzp;
    if( gettimeofday( &tp, &tzp ) == 0 )
        return (double)(tp.tv_sec) + (double)(tp.tv_usec)/1.0e6;
    else
        return 0.0;
}
static double wall_time = 0; 
#endif


#ifdef _PROTO_
long price_out_impl( network_t *net )
#else
long price_out_impl( net )
     network_t *net;
#endif
{
    long i;
    long trips;
    long new_arcs = 0;
    long resized = 0;
    long latest;
    long min_impl_duration = 15;

    register cost_t bigM = net->bigM;
    register cost_t head_potential;
    register cost_t arc_cost = 30;
    register cost_t red_cost;
    register cost_t bigM_minus_min_impl_duration;
        
    register arc_t *arcout, *arcin, *arcnew, *stop;
    register arc_t *first_of_sparse_list;
    register node_t *tail, *head;


#if defined AT_HOME
    wall_time -= Get_Time();
#endif

    
    bigM_minus_min_impl_duration = (cost_t)bigM - min_impl_duration;
    

    
    if( net->n_trips <= MAX_NB_TRIPS_FOR_SMALL_NET )
    {
      if( net->m + net->max_new_m > net->max_m 
          &&
          (net->n_trips*net->n_trips)/2 + net->m > net->max_m
          )
      {
        resized = 1;
        if( resize_prob( net ) )
          return -1;
        
        refresh_neighbour_lists( net );
      }
    }
#if !defined SPEC_STATIC
    else
    {
      if( net->m + net->max_new_m > net->max_m 
          &&
          (net->n_trips*net->n_trips)/2 + net->m > net->max_m
          )
      {
        resized = 1;
        if( resize_prob( net ) )
          return -1;
        
        refresh_neighbour_lists( net );
      }
    }
#endif

        
    arcnew = net->stop_arcs;
    trips = net->n_trips;

    arcout = net->arcs;
    for( i = 0; i < trips && arcout[1].ident == FIXED; i++, arcout += 3 );
    first_of_sparse_list = (arc_t *)NULL;
    for( ; i < trips; i++, arcout += 3 )
    {
        if( arcout[1].ident != FIXED )
        {
            arcout->head->firstout->head->arc_tmp = first_of_sparse_list;
            first_of_sparse_list = arcout + 1;
        }
        
        if( arcout->ident == FIXED )
            continue;
        
        head = arcout->head;
        latest = head->time - arcout->org_cost 
            + (long)bigM_minus_min_impl_duration;
                
        head_potential = head->potential;
        
        arcin = first_of_sparse_list->tail->arc_tmp;
        while( arcin )
        {
            tail = arcin->tail;

            if( tail->time + arcin->org_cost > latest )
            {
                arcin = tail->arc_tmp;
                continue;
            }
            
            red_cost = arc_cost - tail->potential + head->potential;
            
            if( red_cost < 0 )
            {
                if( new_arcs < net->max_residual_new_m )
                {
                    insert_new_arc( arcnew, new_arcs, tail, head, 
                                    arc_cost, red_cost );
                    new_arcs++;                 
                }
                else if( (cost_t)arcnew[0].flow > red_cost )
                    replace_weaker_arc( net, arcnew, tail, head, 
                                        arc_cost, red_cost );
            }

            arcin = tail->arc_tmp;
        }
    }
    
    if( new_arcs )
    {
        arcnew = net->stop_arcs;
        net->stop_arcs += new_arcs;
        stop = (void *)net->stop_arcs;
        if( resized )
        {
            for( ; arcnew != stop; arcnew++ )
            {
                arcnew->flow = (flow_t)0;
                arcnew->ident = AT_LOWER;
            }
        }
        else
        {
            for( ; arcnew != stop; arcnew++ )
            {
                arcnew->flow = (flow_t)0;
                arcnew->ident = AT_LOWER;
                arcnew->nextout = arcnew->tail->firstout;
                arcnew->tail->firstout = arcnew;
                arcnew->nextin = arcnew->head->firstin;
                arcnew->head->firstin = arcnew;
            }
        }
        
        net->m += new_arcs;
        net->m_impl += new_arcs;
        net->max_residual_new_m -= new_arcs;
    }
    

#if defined AT_HOME
    wall_time += Get_Time();
    printf( "total time price_out_impl(): %0.0f\n", wall_time );
#endif


    return new_arcs;
}   






#ifdef _PROTO_
long suspend_impl( network_t *net, cost_t threshold, long all )
#else
long suspend_impl( net, threshold, all )
     network_t *net;
     cost_t threshold;
     long all;
#endif
{
    long susp;
    
    cost_t red_cost;
    arc_t *new_arc, *arc;
    void *stop;

    

    if( all )
        susp = net->m_impl;
    else
    {
        stop = (void *)net->stop_arcs;
        new_arc = &(net->arcs[net->m - net->m_impl]);
        for( susp = 0, arc = new_arc; arc < (arc_t *)stop; arc++ )
        {
            if( arc->ident == AT_LOWER )
                red_cost = arc->cost - arc->tail->potential 
                        + arc->head->potential;
            else
            {
                red_cost = (cost_t)-2;
                
                if( arc->ident == BASIC )
                {
                    if( arc->tail->basic_arc == arc )
                        arc->tail->basic_arc = new_arc;
                    else
                        arc->head->basic_arc = new_arc;
                }
            }
            
            if( red_cost > threshold )
                susp++;
            else
            {
                *new_arc = *arc;
                new_arc++;
            }
        }
    }
    
        
#if defined AT_HOME
    printf( "\nremove %ld arcs\n\n", susp );
    fflush( stdout );
#endif

    if( susp )
    {
        net->m -= susp;
        net->m_impl -= susp;
        net->stop_arcs -= susp;
        net->max_residual_new_m += susp;
        
        refresh_neighbour_lists( net );
    }

    return susp;
}



/**************************************************************************
MCF.H of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 22:10:51 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: mcf.c,v 1.15 2005/02/17 21:43:12 bzfloebe Exp $  */



/*#include "mcf.h"*/

#define REPORT

extern long min_impl_duration;
network_t net;





#ifdef _PROTO_
long global_opt( void )
#else
long global_opt( )
#endif
{
    long new_arcs;
    long residual_nb_it;
    

    new_arcs = -1;
    residual_nb_it = net.n_trips <= MAX_NB_TRIPS_FOR_SMALL_NET ?
        MAX_NB_ITERATIONS_SMALL_NET : MAX_NB_ITERATIONS_LARGE_NET;

    while( new_arcs )
    {
#ifdef REPORT
        printf( "active arcs                : %ld\n", net.m );
#endif

        primal_net_simplex( &net );


#ifdef REPORT
        printf( "simplex iterations         : %ld\n", net.iterations );
        printf( "objective value            : %0.0f\n", flow_cost(&net) );
#endif


#if defined AT_HOME
        printf( "%ld residual iterations\n", residual_nb_it );
#endif

        if( !residual_nb_it )
            break;


        if( net.m_impl )
        {
          new_arcs = suspend_impl( &net, (cost_t)-1, 0 );

#ifdef REPORT
          if( new_arcs )
            printf( "erased arcs                : %ld\n", new_arcs );
#endif
        }


        new_arcs = price_out_impl( &net );

#ifdef REPORT
        if( new_arcs )
            printf( "new implicit arcs          : %ld\n", new_arcs );
#endif
        
        if( new_arcs < 0 )
        {
#ifdef REPORT
            printf( "not enough memory, exit(-1)\n" );
#endif

            exit(-1);
        }

#ifndef REPORT
        printf( "\n" );
#endif


        residual_nb_it--;
    }

    printf( "checksum                   : %ld\n", net.checksum );

    return 0;
}






#ifdef _PROTO_
int main( int argc, char *argv[] )
#else
int main( argc, argv )
    int argc;
    char *argv[];
#endif
{
    if( argc < 2 )
        return -1;


    printf( "\nMCF SPEC CPU2006 version 1.10\n" );
    printf( "Copyright (c) 1998-2000 Zuse Institut Berlin (ZIB)\n" );
    printf( "Copyright (c) 2000-2002 Andreas Loebel & ZIB\n" );
    printf( "Copyright (c) 2003-2005 Andreas Loebel\n" );
    printf( "\n" );



    memset( (void *)(&net), 0, (size_t)sizeof(network_t) );
    net.bigM = (long)BIGM;

    strcpy( net.inputfile, argv[1] );
    
    if( read_min( &net ) )
    {
        printf( "read error, exit\n" );
        getfree( &net );
        return -1;
    }


#ifdef REPORT
    printf( "nodes                      : %ld\n", net.n_trips );
#endif


    primal_start_artificial( &net );
    global_opt( );


#ifdef REPORT
    printf( "done\n" );
#endif

    

    if( write_circulations( "mcf.out", &net ) )
    {
        getfree( &net );
        return -1;    
    }


    getfree( &net );
    return 0;
}
/**************************************************************************
MCFUTIL.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 21:46:06 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: mcfutil.c,v 1.11 2005/02/17 21:43:12 bzfloebe Exp $  */



/*#include "mcfutil.h"*/



#ifdef _PROTO_
void refresh_neighbour_lists( network_t *net )
#else
void refresh_neighbour_lists( net )
    network_t *net;
#endif
{
    node_t *node;
    arc_t *arc;
    void *stop;
        

    node = net->nodes;
    for( stop = (void *)net->stop_nodes; node < (node_t *)stop; node++ )
    {
        node->firstin = (arc_t *)NULL;
        node->firstout = (arc_t *)NULL;
    }
    
    arc = net->arcs;
    for( stop = (void *)net->stop_arcs; arc < (arc_t *)stop; arc++ )
    {
        arc->nextout = arc->tail->firstout;
        arc->tail->firstout = arc;
        arc->nextin = arc->head->firstin;
        arc->head->firstin = arc;
    }
    
    return;
}










#ifdef _PROTO_
long refresh_potential( network_t *net )
#else
long refresh_potential( net )
    network_t *net;
#endif
{
    node_t *node, *tmp;
    node_t *root = net->nodes;
    long checksum = 0;
    

    root->potential = (cost_t) -MAX_ART_COST;
    tmp = node = root->child;
    while( node != root )
    {
        while( node )
        {
            if( node->orientation == UP )
                node->potential = node->basic_arc->cost + node->pred->potential;
            else /* == DOWN */
            {
                node->potential = node->pred->potential - node->basic_arc->cost;
                checksum++;
            }

            tmp = node;
            node = node->child;
        }
        
        node = tmp;

        while( node->pred )
        {
            tmp = node->sibling;
            if( tmp )
            {
                node = tmp;
                break;
            }
            else
                node = node->pred;
        }
    }
    
    return checksum;
}







#ifdef _PROTO_
double flow_cost( network_t *net )
#else
double flow_cost( net )
    network_t *net;
#endif
{
    arc_t *arc;
    node_t *node;
    void *stop;
    
    long fleet = 0;
    cost_t operational_cost = 0;
    

    stop = (void *)net->stop_arcs;
    for( arc = net->arcs; arc != (arc_t *)stop; arc++ )
    {
        if( arc->ident == AT_UPPER )
            arc->flow = (flow_t)1;
        else
            arc->flow = (flow_t)0;
    }

    stop = (void *)net->stop_nodes;
    for( node = net->nodes, node++; node != (node_t *)stop; node++ )
        node->basic_arc->flow = node->flow;
    
    stop = (void *)net->stop_arcs;
    for( arc = net->arcs; arc != (arc_t *)stop; arc++ )
    {
        if( arc->flow )
        {
            if( !(arc->tail->number < 0 && arc->head->number > 0) )
            {
                if( !arc->tail->number )
                {
                    operational_cost += (arc->cost - net->bigM);
                    fleet++;
                }
                else
                    operational_cost += arc->cost;
            }
        }

    }
    
    return (double)fleet * (double)net->bigM + (double)operational_cost;
}










#ifdef _PROTO_
double flow_org_cost( network_t *net )
#else
double flow_org_cost( net )
    network_t *net;
#endif
{
    arc_t *arc;
    node_t *node;
    void *stop;
    
    long fleet = 0;
    cost_t operational_cost = 0;
    

    stop = (void *)net->stop_arcs;
    for( arc = net->arcs; arc != (arc_t *)stop; arc++ )
    {
        if( arc->ident == AT_UPPER )
            arc->flow = (flow_t)1;
        else
            arc->flow = (flow_t)0;
    }

    stop = (void *)net->stop_nodes;
    for( node = net->nodes, node++; node != (node_t *)stop; node++ )
        node->basic_arc->flow = node->flow;
    
    stop = (void *)net->stop_arcs;
    for( arc = net->arcs; arc != (arc_t *)stop; arc++ )
    {
        if( arc->flow )
        {
            if( !(arc->tail->number < 0 && arc->head->number > 0) )
            {
                if( !arc->tail->number )
                {
                    operational_cost += (arc->org_cost - net->bigM);
                    fleet++;
                }
                else
                    operational_cost += arc->org_cost;
            }
        }
    }
    
    return (double)fleet * (double)net->bigM + (double)operational_cost;
}










#ifdef _PROTO_
long primal_feasible( network_t *net )
#else
long primal_feasible( net )
    network_t *net;
#endif
{
    void *stop;
    node_t *node;
    arc_t *dummy = net->dummy_arcs;
    arc_t *stop_dummy = net->stop_dummy;
    arc_t *arc;
    flow_t flow;
    

    node = net->nodes;
    stop = (void *)net->stop_nodes;

    for( node++; node < (node_t *)stop; node++ )
    {
        arc = node->basic_arc;
        flow = node->flow;
        if( arc >= dummy && arc < stop_dummy )
        {
            if( ABS(flow) > (flow_t)net->feas_tol )
            {
                printf( "PRIMAL NETWORK SIMPLEX: " );
                printf( "artificial arc with nonzero flow, node %d (%ld)\n",
                        node->number, flow );
            }
        }
        else
        {
            if( flow < (flow_t)(-net->feas_tol)
               || flow - (flow_t)1 > (flow_t)net->feas_tol )
            {
                printf( "PRIMAL NETWORK SIMPLEX: " );
                printf( "basis primal infeasible (%ld)\n", flow );
                net->feasible = 0;
                return 1;
            }
        }
    }
    
    net->feasible = 1;
    
    return 0;
}










#ifdef _PROTO_
long dual_feasible( network_t *net )
#else
long dual_feasible(  net )
    network_t *net;
#endif
{
    arc_t         *arc;
    arc_t         *stop     = net->stop_arcs;
    cost_t        red_cost;
    
    

    for( arc = net->arcs; arc < stop; arc++ )
    {
        red_cost = arc->cost - arc->tail->potential 
            + arc->head->potential;
        switch( arc->ident )
        {
        case BASIC:
#ifdef AT_ZERO
        case AT_ZERO:
            if( ABS(red_cost) > (cost_t)net->feas_tol )
#ifdef DEBUG
                printf("%d %d %d %ld\n", arc->tail->number, arc->head->number,
                       arc->ident, red_cost );
#else
                goto DUAL_INFEAS;
#endif
            
            break;
#endif
        case AT_LOWER:
            if( red_cost < (cost_t)-net->feas_tol )
#ifdef DEBUG
                printf("%d %d %d %ld\n", arc->tail->number, arc->head->number,
                       arc->ident, red_cost );
#else
                goto DUAL_INFEAS;
#endif

            break;
        case AT_UPPER:
            if( red_cost > (cost_t)net->feas_tol )
#ifdef DEBUG
                printf("%d %d %d %ld\n", arc->tail->number, arc->head->number,
                       arc->ident, red_cost );
#else
                goto DUAL_INFEAS;
#endif

            break;
        case FIXED:
        default:
            break;
        }
    }
    
    return 0;
    
DUAL_INFEAS:
    fprintf( stderr, "DUAL NETWORK SIMPLEX: " );
    fprintf( stderr, "basis dual infeasible\n" );
    return 1;
}







#ifdef _PROTO_
long getfree( 
            network_t *net
            )
#else
long getfree( net )
     network_t *net;
#endif
{  
    FREE( net->nodes );
    FREE( net->arcs );
    FREE( net->dummy_arcs );
    net->nodes = net->stop_nodes = NULL;
    net->arcs = net->stop_arcs = NULL;
    net->dummy_arcs = net->stop_dummy = NULL;

    return 0;
}



/**************************************************************************
OUTPUT.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:21:54 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: output.c,v 1.10 2005/02/17 19:42:33 bzfloebe Exp $  */



/*#include "output.h"*/





#ifdef _PROTO_
long write_circulations(
                   char *outfile,
                   network_t *net
                   )
#else
long write_circulations( outfile, net )
     char *outfile;
     network_t *net;
#endif 
{
    FILE *out = NULL;
    arc_t *block;
    arc_t *arc;
    arc_t *arc2;
    arc_t *first_impl = net->stop_arcs - net->m_impl;

    if(( out = fopen( outfile, "w" )) == NULL )
        return -1;

    refresh_neighbour_lists( net );
    
    for( block = net->nodes[net->n].firstout; block; block = block->nextout )
    {
        if( block->flow )
        {
            fprintf( out, "()\n" );
            
            arc = block;
            while( arc )
            {
                if( arc >= first_impl )
                    fprintf( out, "***\n" );

                fprintf( out, "%d\n", - arc->head->number );
                arc2 = arc->head[net->n_trips].firstout; 
                for( ; arc2; arc2 = arc2->nextout )
                    if( arc2->flow )
                        break;
                if( !arc2 )
                {
                    fclose( out );
                    return -1;
                }
                
                if( arc2->head->number )
                    arc = arc2;
                else
                    arc = NULL;
            }
        }
    }
    


    fclose(out);
    
    return 0;
}
/**************************************************************************
PBEAMPP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:04 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbeampp.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



#define K 300
#define B  50



/*#include "pbeampp.h"*/




#ifdef _PROTO_
int bea_is_dual_infeasible( arc_t *arc, cost_t red_cost )
#else
int bea_is_dual_infeasible( arc, red_cost )
    arc_t *arc;
    cost_t red_cost;
#endif
{
    return(    (red_cost < 0 && arc->ident == AT_LOWER)
            || (red_cost > 0 && arc->ident == AT_UPPER) );
}







typedef struct basket
{
    arc_t *a;
    cost_t cost;
    cost_t abs_cost;
} BASKET;

static long basket_size;
static BASKET basket[B+K+1];
static BASKET *perm[B+K+1];



#ifdef _PROTO_
void sort_basket( long min, long max )
#else
void sort_basket( min, max )
    long min, max;
#endif
{
    long l, r;
    cost_t cut;
    BASKET *xchange;

    l = min; r = max;

    cut = perm[ (long)( (l+r) / 2 ) ]->abs_cost;

    do
    {
        while( perm[l]->abs_cost > cut )
            l++;
        while( cut > perm[r]->abs_cost )
            r--;
            
        if( l < r )
        {
            xchange = perm[l];
            perm[l] = perm[r];
            perm[r] = xchange;
        }
        if( l <= r )
        {
            l++; r--;
        }

    }
    while( l <= r );

    if( min < r )
        sort_basket( min, r );
    if( l < max && l <= B )
        sort_basket( l, max ); 
}






static long nr_group;
static long group_pos;


static long initialize = 1;


#ifdef _PROTO_
arc_t *primal_bea_mpp( long m,  arc_t *arcs, arc_t *stop_arcs, 
                              cost_t *red_cost_of_bea )
#else
arc_t *primal_bea_mpp( m, arcs, stop_arcs, red_cost_of_bea )
    long m;
    arc_t *arcs;
    arc_t *stop_arcs;
    cost_t *red_cost_of_bea;
#endif
{
    long i, next, old_group_pos;
    arc_t *arc;
    cost_t red_cost;

    if( initialize )
    {
        for( i=1; i < K+B+1; i++ )
            perm[i] = &(basket[i]);
        nr_group = ( (m-1) / K ) + 1;
        group_pos = 0;
        basket_size = 0;
        initialize = 0;
    }
    else
    {
        for( i = 2, next = 0; i <= B && i <= basket_size; i++ )
        {
            arc = perm[i]->a;
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if( (red_cost < 0 && arc->ident == AT_LOWER)
                || (red_cost > 0 && arc->ident == AT_UPPER) )
            {
                next++;
                perm[next]->a = arc;
                perm[next]->cost = red_cost;
                perm[next]->abs_cost = ABS(red_cost);
            }
                }   
        basket_size = next;
        }

    old_group_pos = group_pos;

NEXT:
    /* price next group */
    arc = arcs + group_pos;
    for( ; arc < stop_arcs; arc += nr_group )
    {
        if( arc->ident > BASIC )
        {
            /* red_cost = bea_compute_red_cost( arc ); */
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if( bea_is_dual_infeasible( arc, red_cost ) )
            {
                basket_size++;
                perm[basket_size]->a = arc;
                perm[basket_size]->cost = red_cost;
                perm[basket_size]->abs_cost = ABS(red_cost);
            }
        }
        
    }

    if( ++group_pos == nr_group )
        group_pos = 0;

    if( basket_size < B && group_pos != old_group_pos )
        goto NEXT;

    if( basket_size == 0 )
    {
        initialize = 1;
        *red_cost_of_bea = 0; 
        return NULL;
    }
    
    sort_basket( 1, basket_size );
    
    *red_cost_of_bea = perm[1]->cost;
    return( perm[1]->a );
}










/**************************************************************************
PBLA.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:14 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbla.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



/*#include "pbla.h"*/





#define TEST_MIN( nod, ex, comp, op ) \
{ \
      if( *delta op (comp) ) \
      { \
            iminus = nod; \
            *delta = (comp); \
            *xchange = ex; \
      } \
}


#ifdef _PROTO_
node_t *primal_iminus( 
                      flow_t *delta,
                      long *xchange,
                      node_t *iplus, 
                      node_t*jplus,
                      node_t **w
                    )
#else
node_t *primal_iminus( delta, xchange, iplus, jplus, w )
    flow_t *delta;
    long *xchange;
    node_t *iplus, *jplus;
    node_t **w;
#endif
{
    node_t *iminus = NULL;
    

    while( iplus != jplus )
    {
        if( iplus->depth < jplus->depth )
        {
            if( iplus->orientation )
                TEST_MIN( iplus, 0, iplus->flow, > )
            else if( iplus->pred->pred )
                TEST_MIN( iplus, 0, (flow_t)1 - iplus->flow, > )
            iplus = iplus->pred;
        }
        else
        {
            if( !jplus->orientation )
                TEST_MIN( jplus, 1, jplus->flow, >= )
            else if( jplus->pred->pred )
                TEST_MIN( jplus, 1, (flow_t)1 - jplus->flow, >= )
            jplus = jplus->pred;
        }
    } 

    *w = iplus;

    return iminus;
}
/**************************************************************************
PFLOWUP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:23 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pflowup.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



/*#include "pflowup.h"*/




#ifdef _PROTO_
void primal_update_flow( 
                 node_t *iplus,
                 node_t *jplus,
                 node_t *w
                 )
#else
void primal_update_flow( iplus, jplus, w )
    node_t *iplus, *jplus;
    node_t *w; 
#endif
{
    for( ; iplus != w; iplus = iplus->pred )
    {
        if( iplus->orientation )
            iplus->flow = (flow_t)0;
        else
            iplus->flow = (flow_t)1;
    }

    for( ; jplus != w; jplus = jplus->pred )
    {
        if( jplus->orientation )
            jplus->flow = (flow_t)1;
        else
            jplus->flow = (flow_t)0;
    }
}
/**************************************************************************
PSIMPLEX.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:44 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: psimplex.c,v 1.9 2005/02/17 19:44:56 bzfloebe Exp $  */



#undef DEBUG

/*#include "psimplex.h"*/



#ifdef _PROTO_
long primal_net_simplex( network_t *net )
#else
long primal_net_simplex(  net )
    network_t *net;
#endif
{
    flow_t        delta;
    flow_t        new_flow;
    long          opt = 0;
    long          xchange;
    long          new_orientation;
    node_t        *iplus;
    node_t        *jplus; 
    node_t        *iminus;
    node_t        *jminus;
    node_t        *w; 
    arc_t         *bea;
    arc_t         *bla;
    arc_t         *arcs          = net->arcs;
    arc_t         *stop_arcs     = net->stop_arcs;
    node_t        *temp;
    long          m = net->m;
    long          new_set;
    cost_t        red_cost_of_bea;
    long          *iterations = &(net->iterations);
    long          *bound_exchanges = &(net->bound_exchanges);
    long          *checksum = &(net->checksum);


    while( !opt )
    {       
        if( (bea = primal_bea_mpp( m, arcs, stop_arcs, &red_cost_of_bea )) )
        {
            (*iterations)++;

#ifdef DEBUG
            printf( "it %ld: bea = (%ld,%ld), red_cost = %ld\n", 
                    *iterations, bea->tail->number, bea->head->number,
                    red_cost_of_bea );
#endif

            if( red_cost_of_bea > ZERO ) 
            {
                iplus = bea->head;
                jplus = bea->tail;
            }
            else 
            {
                iplus = bea->tail;
                jplus = bea->head;
            }

            delta = (flow_t)1;
            iminus = primal_iminus( &delta, &xchange, iplus, 
                    jplus, &w );

            if( !iminus )
            {
                (*bound_exchanges)++;
                
                if( bea->ident == AT_UPPER)
                    bea->ident = AT_LOWER;
                else
                    bea->ident = AT_UPPER;

                if( delta )
                    primal_update_flow( iplus, jplus, w );
            }
            else 
            {
                if( xchange )
                {
                    temp = jplus;
                    jplus = iplus;
                    iplus = temp;
                }

                jminus = iminus->pred;

                bla = iminus->basic_arc;
                 
                if( xchange != iminus->orientation )
                    new_set = AT_LOWER;
                else
                    new_set = AT_UPPER;

                if( red_cost_of_bea > 0 )
                    new_flow = (flow_t)1 - delta;
                else
                    new_flow = delta;

                if( bea->tail == iplus )
                    new_orientation = UP;
                else
                    new_orientation = DOWN;

                update_tree( !xchange, new_orientation,
                            delta, new_flow, iplus, jplus, iminus, 
                            jminus, w, bea, red_cost_of_bea,
                            (flow_t)net->feas_tol );

                bea->ident = BASIC; 
                bla->ident = new_set;
               
                if( !((*iterations-1) % 200) )
                {
                    *checksum += refresh_potential( net );
#if defined AT_HOME
                    if( *checksum > 2000000000l )
                    {
                        printf( "%ld\n", *checksum );
                        fflush(stdout);
                    }
#endif
                }   
            }
        }
        else
            opt = 1;
    }


    *checksum += refresh_potential( net );
    primal_feasible( net );
    dual_feasible( net );
    
    return 0;
}



/**************************************************************************
PSTART.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:53 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pstart.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



/*#include "pstart.h"*/




#ifdef _PROTO_ 
long primal_start_artificial( network_t *net )
#else
long primal_start_artificial( net )
    network_t *net;
#endif
{      
    node_t *node, *root;
    arc_t *arc;
    void *stop;


    root = node = net->nodes; node++;
    root->basic_arc = NULL;
    root->pred = NULL;
    root->child = node;
    root->sibling = NULL;
    root->sibling_prev = NULL;
    root->depth = (net->n) + 1;
    root->orientation = 0;
    root->potential = (cost_t) -MAX_ART_COST;
    root->flow = ZERO;

    stop = (void *)net->stop_arcs;
    for( arc = net->arcs; arc != (arc_t *)stop; arc++ )
        if( arc->ident != FIXED )
            arc->ident = AT_LOWER;

    arc = net->dummy_arcs;
    for( stop = (void *)net->stop_nodes; node != (node_t *)stop; arc++, node++ )
    {
        node->basic_arc = arc;
        node->pred = root;
        node->child = NULL;
        node->sibling = node + 1; 
        node->sibling_prev = node - 1;
        node->depth = 1;

        arc->cost = (cost_t) MAX_ART_COST;
        arc->ident = BASIC;

        node->orientation = UP; 
        node->potential = ZERO;
        arc->tail = node;
        arc->head = root;                
        node->flow = (flow_t)0;
    }

    node--; root++;
    node->sibling = NULL;
    root->sibling_prev = NULL;

    return 0;
}
/**************************************************************************
READMIN.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Thu Feb 17 20:44:29 2005 by Andreas Loebel (boss.local.de)  */
/*  $Id: readmin.c,v 1.16 2005/02/17 19:44:40 bzfloebe Exp $  */



/*#include "readmin.h"*/




#ifdef _PROTO_
long read_min( network_t *net )
#else
long read_min( net )
     network_t *net;
#endif
{                                       
    FILE *in = NULL;
    char instring[201];
    long t, h, c;
    long i;
    arc_t *arc;
    node_t *node;


    if(( in = fopen( net->inputfile, "r")) == NULL )
        return -1;

    fgets( instring, 200, in );
    if( sscanf( instring, "%ld %ld", &t, &h ) != 2 )
        return -1;
    

    net->n_trips = t;
    net->m_org = h;
    net->n = (t+t+1); 
    net->m = (t+t+t+h);

    if( net->n_trips <= MAX_NB_TRIPS_FOR_SMALL_NET )
    {
      net->max_m = net->m;
      net->max_new_m = MAX_NEW_ARCS_SMALL_NET;
    }
    else
    {
#ifdef SPEC_STATIC
/*
      //net->max_m = 0x1c00000l;
*/
      net->max_m = 0x1a10000l;
#else
      net->max_m = MAX( net->m + MAX_NEW_ARCS, STRECHT(STRECHT(net->m)) );
#endif
      net->max_new_m = MAX_NEW_ARCS_LARGE_NET;
    }

    net->max_residual_new_m = net->max_m - net->m;


    assert( net->max_new_m >= 3 );

    
    net->nodes      = (node_t *) calloc( net->n + 1, sizeof(node_t) );
    net->dummy_arcs = (arc_t *)  calloc( net->n,   sizeof(arc_t) );
    net->arcs       = (arc_t *)  calloc( net->max_m,   sizeof(arc_t) );

    if( !( net->nodes && net->arcs && net->dummy_arcs ) )
    {
      printf( "read_min(): not enough memory\n" );
      getfree( net );
      return -1;
    }


#if defined AT_HOME
    printf( "malloc for nodes       MB %4ld\n", 
            (long)((net->n + 1)*sizeof(node_t) / 0x100000) );
    printf( "malloc for dummy arcs  MB %4ld\n", 
            (long)((net->n)*sizeof(arc_t) / 0x100000) );
    printf( "malloc for arcs        MB %4ld\n", 
            (long)((net->max_m)*sizeof(arc_t) / 0x100000) );
    printf( "------------------------------\n" );
    printf( "heap about             MB %4ld\n\n", 
            (long)((net->n +1)*sizeof(node_t) / 0x100000)
            +(long)((net->n)*sizeof(arc_t) / 0x100000)
            +(long)((net->max_m)*sizeof(arc_t) / 0x100000)
            );
#endif


    net->stop_nodes = net->nodes + net->n + 1; 
    net->stop_arcs  = net->arcs + net->m;
    net->stop_dummy = net->dummy_arcs + net->n;


    node = net->nodes;
    arc = net->arcs;

    for( i = 1; i <= net->n_trips; i++ )
    {
        fgets( instring, 200, in );

        if( sscanf( instring, "%ld %ld", &t, &h ) != 2 || t > h )
            return -1;

        node[i].number = -i;
        node[i].flow = (flow_t)-1;
            
        node[i+net->n_trips].number = i;
        node[i+net->n_trips].flow = (flow_t)1;
        
        node[i].time = t;
        node[i+net->n_trips].time = h;

        arc->tail = &(node[net->n]);
        arc->head = &(node[i]);
        arc->org_cost = arc->cost = (cost_t)(net->bigM+15);
        arc->nextout = arc->tail->firstout;
        arc->tail->firstout = arc;
        arc->nextin = arc->head->firstin;
        arc->head->firstin = arc; 
        arc++;
                                    
        arc->tail = &(node[i+net->n_trips]);
        arc->head = &(node[net->n]);
        arc->org_cost = arc->cost = (cost_t)15;
        arc->nextout = arc->tail->firstout;
        arc->tail->firstout = arc;
        arc->nextin = arc->head->firstin;
        arc->head->firstin = arc; 
        arc++;

        arc->tail = &(node[i]);
        arc->head = &(node[i+net->n_trips]);
        arc->org_cost = arc->cost = (cost_t)(2*MAX(net->bigM,(long)BIGM));
        arc->nextout = arc->tail->firstout;
        arc->tail->firstout = arc;
        arc->nextin = arc->head->firstin;
        arc->head->firstin = arc; 
        arc++;
    }

    
    if( i != net->n_trips + 1 )
        return -1;


    for( i = 0; i < net->m_org; i++, arc++ )
    {
        fgets( instring, 200, in );
        
        if( sscanf( instring, "%ld %ld %ld", &t, &h, &c ) != 3 )
                return -1;

        arc->tail = &(node[t+net->n_trips]);
        arc->head = &(node[h]);
        arc->org_cost = (cost_t)c;
        arc->cost = (cost_t)c;
        arc->nextout = arc->tail->firstout;
        arc->tail->firstout = arc;
        arc->nextin = arc->head->firstin;
        arc->head->firstin = arc; 
    }


    if( net->stop_arcs != arc )
    {
        net->stop_arcs = arc;
        arc = net->arcs;
        for( net->m = 0; arc < net->stop_arcs; arc++ )
            (net->m)++;
        net->m_org = net->m;
    }
    
    fclose( in );


    net->clustfile[0] = (char)0;
        
    for( i = 1; i <= net->n_trips; i++ )
    {
        net->arcs[3*i-1].cost = 
            (cost_t)((-2)*MAX(net->bigM,(long) BIGM));
        net->arcs[3*i-1].org_cost = 
            (cost_t)((-2)*(MAX(net->bigM,(long) BIGM)));
    }
    
    
    return 0;
}
/**************************************************************************
TREEUP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:23:12 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: treeup.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



/*#include "treeup.h"*/




#ifdef _PROTO_
void update_tree( 
                 long cycle_ori,
                 long new_orientation,
                 flow_t delta,
                 flow_t new_flow,
                 node_t *iplus,
                 node_t *jplus,
                 node_t *iminus,
                 node_t *jminus,
                 node_t *w,
                 arc_t *bea,
                 cost_t sigma,
                 flow_t feas_tol
                )
#else
void update_tree( cycle_ori, new_orientation, delta, new_flow, 
                 iplus, jplus, iminus, jminus, w, bea, sigma, feas_tol )
     long cycle_ori;
     long new_orientation;
     flow_t delta; 
     flow_t new_flow;
     node_t *iplus, *jplus;
     node_t *iminus, *jminus;
     node_t *w;
     arc_t *bea;
     cost_t sigma; 
     flow_t feas_tol;
#endif
{
    arc_t    *basic_arc_temp;
    arc_t    *new_basic_arc;  
    node_t   *father;         
    node_t   *temp;           
    node_t   *new_pred;       
    long     orientation_temp;
    long     depth_temp;      
    long     depth_iminus;    
    long     new_depth;       
    flow_t   flow_temp;       


    /**/
    if( (bea->tail == jplus && sigma < 0) ||
        (bea->tail == iplus && sigma > 0) )
        sigma = ABS(sigma);
    else
        sigma = -(ABS(sigma));
    
    father = iminus;
    father->potential += sigma;
 RECURSION:
    temp = father->child;
    if( temp )
    {
    ITERATION:
        temp->potential += sigma;
        father = temp;
        goto RECURSION;
    }
 TEST:
    if( father == iminus )
        goto CONTINUE;
    temp = father->sibling;
    if( temp )
        goto ITERATION;
    father = father->pred;
    goto TEST;
    
 CONTINUE:
    /**/


    temp = iplus;
    father = temp->pred;
    new_depth = depth_iminus = iminus->depth;
    new_pred = jplus;
    new_basic_arc = bea;
    while( temp != jminus )
    {
        if( temp->sibling )
            temp->sibling->sibling_prev = temp->sibling_prev;
        if( temp->sibling_prev )
            temp->sibling_prev->sibling = temp->sibling;
        else father->child = temp->sibling;


        temp->pred = new_pred;
        temp->sibling = new_pred->child;
        if( temp->sibling )
            temp->sibling->sibling_prev = temp;
        new_pred->child = temp;
        temp->sibling_prev = 0;

        orientation_temp = !(temp->orientation); 
        if( orientation_temp == cycle_ori )
            flow_temp = temp->flow + delta;
        else
            flow_temp = temp->flow - delta;
        basic_arc_temp = temp->basic_arc;
        depth_temp = temp->depth;

        temp->orientation = new_orientation;
        temp->flow = new_flow;
        temp->basic_arc = new_basic_arc;
        temp->depth = new_depth;

        new_pred = temp;
        new_orientation = orientation_temp;
        new_flow = flow_temp;
        new_basic_arc = basic_arc_temp;
        new_depth = depth_iminus - depth_temp;      
        temp = father;
        father = temp->pred;
    } 

    if( delta > feas_tol )
    {
        for( temp = jminus; temp != w; temp = temp->pred )
        {
            temp->depth -= depth_iminus;
            if( temp->orientation != cycle_ori )
                temp->flow += delta;
            else
                temp->flow -= delta;
        }
        for( temp = jplus; temp != w; temp = temp->pred )
        {
            temp->depth += depth_iminus;
            if( temp->orientation == cycle_ori )
                temp->flow += delta;
            else
                temp->flow -= delta;
        }
    }
    else
    {
        for( temp = jminus; temp != w; temp = temp->pred )
            temp->depth -= depth_iminus;
        for( temp = jplus; temp != w; temp = temp->pred )
            temp->depth += depth_iminus;
    }

}


