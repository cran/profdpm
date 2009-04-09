#ifndef PDPMLM_H
#define PDPMLM_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "util.h"

// The following should be changed if ported to another interface
// In addition, allocated memory may need to be freed if ported to
// another interface. This is not necessary if allocated using R's
// memory manager.
#define pdpmlm_alloc(count, size)  R_alloc(count, size)
#define pdpmlm_printf Rprintf

#define DEFAULT_ALP    1.000
#define DEFAULT_A0     0.001
#define DEFAULT_B0     0.001
#define DEFAULT_M0     0.000
#define DEFAULT_S0     1.000

#define BAD_CLS INT_MAX

// bit masks for flags
#define FLAG_VERBOSE  1<<0  // should routine be verbose
#define FLAG_OPTCRIT  1<<1  // has optimization criterion been met
#define FLAG_PRICLUS  1<<2  // use Cluster prior instead of Dirichlet prior
#define FLAG_EMPTY_3  1<<3  // not used
#define FLAG_EMPTY_4  1<<4  // not used
#define FLAG_EMPTY_5  1<<5  // not used
#define FLAG_EMPTY_6  1<<6  // not used
#define FLAG_EMPTY_7  1<<7  // not used


typedef struct {

unsigned char   flags;  // some options

double          alp;  // prior alpha parameter
double          s0;   // prior s0 parameter
double        * m0;   // prior m0 parameter
double          a0;   // prior a0 parameter
double          b0;   // prior b0 parameter

// Each of the p entries in y has a corresponding entry in vgr
// that indicates group membership. The values in vgr are begin
// at zero and are increasing. The largest value is ngr-1.
// Hence, y and x are sorted before being passed to the C code.
unsigned int  * vgr;  // group vector {0,0,...,ngr-1,ngr-1}
unsigned int  * pgr;  // number in each group (ngr)
unsigned int    ngr;  // number of groups

// Each of the ngr groups has an entry in vcl indicating
// group membership. The group indicator is used to index
// the values in vcl. Hence, 'vcl[ 0 ]' would yield the 
// cluster indicator for group 0. The values in vcl may 
// range from 0 to ncl-1, but must be less than or equal to
// ngr-1. These values may not be ordered.
unsigned int  * vcl;  // cluster vector {0,...,ncl-1} 
unsigned int  * pcl;  // number in each cluster (ngr)
unsigned int    ncl;  // number of clusters

double        * y;    // y vector (px1)
double        * x;    // x matrix (pxq)
unsigned int    p;    // nrow(x)
unsigned int    q;    // ncol(x)

double       ** xxgr;   // x'x matrix (qxq) for each group
double       ** xygr;   // x'y vector (qx1) for each group
double        * yygr;   // y'y scalar for each group

double       ** xxcl;   // x'x matrix (qxq) for each cluster
double       ** xycl;   // x'y vector (qx1) for each cluster
double        * yycl;   // y'y scalar for each cluster

double        * s;      // storage for an s matrix (qxq)
double        * m;      // storage for an m vector (qx1)
double          a;      // storage for an a scalar
double          b;      // storage for an b scalar

double        * fbuf;   // temporary storage for fortran routines (qx1)
unsigned int  * pbuf;   // temporary storage for pdpmlm routines (ngrx1)
unsigned int    mem;    // memory usage counter

} pdpmlm_t;

void memerror(); 

// Assign the observations/groups according to a simple algorithm
void         pdpmlm_divy( pdpmlm_t * obj );

// Add an observation/group to a cluster cls
void         pdpmlm_add( pdpmlm_t * obj, unsigned grp, unsigned int cls );

// Remove an observation/group from cluster clt 
void         pdpmlm_sub( pdpmlm_t * obj, unsigned grp, unsigned int cls );

// Move an observation/group to cluster clt
void         pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Move an observation/group to cluster cls, return the resulting change in logp
double       pdpmlm_movep( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Get the index of an free (empty) cluster of BAD_CLS if none exist
unsigned int pdpmlm_free( pdpmlm_t * obj );

// Try to split cluster cls
double       pdpmlm_split( pdpmlm_t * obj, unsigned int cls );

// Try to merge cluster cls
double       pdpmlm_merge( pdpmlm_t * obj, unsigned int cls );

// Move an observation/group to the cluster that minimizes the logp
void         pdpmlm_best( pdpmlm_t * obj, unsigned int grp );

// Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmlm_parm( pdpmlm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b );

// Compute part of the log posterior value for a particular cluster
double       pdpmlm_logpcls( pdpmlm_t * obj, unsigned int cls );

// Compute the log posterior value for the model 
double       pdpmlm_logp( pdpmlm_t * obj );

// Perform split merge style optimization
void         pdpmlm_spmer( pdpmlm_t * obj, unsigned int itermax, double crit);

// Perform chunk style optimization
void         pdpmlm_chunk( pdpmlm_t * obj, unsigned int itermax, double crit);
#endif
