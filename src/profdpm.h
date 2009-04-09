#ifndef PROFDPM_H
#define PROFDPM_H

#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

//some R specific macros
#define pdpm_printf Rprintf
#define pdpm_runif  runif

//utility functions
#define DEBUG pdpm_printf("F: %s, C: %u\n", __FUNCTION__, __COUNTER__)
/* log(sqrt(2*pi)) == log(2*pi)/2 */
#define LN_SQRT_2PI 0.918938533204672741780329736406
#define rlcat_runif runif

/* stats */
double       lfactorial(unsigned int x);
unsigned int rlcat(double * logp, unsigned int n);

/* R */
SEXP         getListElementByName(SEXP list, const char * name);

//address upper triangular packed storage matrix
//a00, a01, a11, a02, a12, a22, a03, a13, a23, a33, ...
//0 <= i <= j < n
#define UMAT(i, j) (i + j * ( j + 1 ) / 2)

//address symmetric full storage (by column) matrix
//a00, a10, ..., an0, a01, a11, ..., an1, ...
//0 <= i, j < n
#define FMAT(i, j, n) (i + j * n)

//address lower triangular packed storage matrix
//a00, a10, ... ,an0, a11, a21, ..., an1, a22, a32, ...
//0 <= j <= i < n
#define LMAT(i, j, n) (i + j * ( 2 * n - j + 1 ) / 2)

//absolute value
#define ABS(x) (x < 0 ? -x : x)

//pdpm defaults
#define DEFAULT_LAM    0.000
#define DEFAULT_ALP    1.000
#define BAD_VCL UINT_MAX

//pdpmlm defaults
#define DEFAULT_LM_A0     0.001
#define DEFAULT_LM_B0     0.001
#define DEFAULT_LM_M0     0.000
#define DEFAULT_LM_S0     1.000

//pdpmbm defaults
#define DEFAULT_BM_A0     1.000
#define DEFAULT_BM_B0     1.000

//bit masks for flags
#define FLAG_VERBOSE  1<<0  //should routine be verbose
#define FLAG_OPTCRIT  1<<1  //has optimization criterion been met
#define FLAG_DIRICHL  1<<2  //use cluster Dirichlet prior
#define FLAG_SINGULA  1<<3  //singularities in optimization
#define FLAG_EMPTY_4  1<<4  //not used
#define FLAG_EMPTY_5  1<<5  //not used
#define FLAG_EMPTY_6  1<<6  //not used
#define FLAG_EMPTY_7  1<<7  //not used

typedef struct pdpm_t {

unsigned char   flags;  //some options

double          lam;  //prior lambda parameter
double          alp;  //prior alpha parameter

//ngr - number of groups
//ngr is the number of distinct group values
unsigned int    ngr;  // number of groups

//ncl - number of clusters. 
//Each group is assigned to exactly one cluster. Each 
//cluster is made up of one or more groups. ncl is the current 
//number of clusters. ncl must be less than or equal to ngr.
unsigned int    ncl;

//gcl - number of groups in each cluster
//Hence, gcl[ cls ] is the number of groups in cluster
//cls (length ngr)
unsigned int  * gcl;  

//vcl - cluster membership vector
//Each of the ngr groups has an entry in vcl indicating
//group membership. The group indicator is used to index
//the values in vcl. Hence, 'vcl[ 0 ]' would yield the 
//cluster indicator for group 0. The values in vcl may 
//range from 0 to ngr-1. These values are not ordered.
//However, there will always be ncl distinct values 
//other than BAD_CLS. (length ngr)
unsigned int  * vcl;

//unnormalized log posterior value
double          logpval;

//temporary storage for pdpmlm routines (array of length ngr)
unsigned int  * pbuf;

//memory usage counter
unsigned int    mem;    

//Assign the groups according to a simple algorithm
void         (*divy)( struct pdpm_t * obj);
//Add a group to a cluster cls
void         (*add)( struct pdpm_t * obj, unsigned int grp, unsigned int cls );
//Remove a group from cluster cls
void         (*sub)( struct pdpm_t * obj, unsigned int grp, unsigned int cls );
//Compute part of the log posterior value for a particular cluster
double       (*logpcls)( struct pdpm_t * obj, unsigned int cls );
//Compute the log posterior value for the current partition 
double       (*logp)( struct pdpm_t * obj );
//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       (*logponly)( struct pdpm_t * obj, unsigned int * only, unsigned int size );

//model struct pointer
void        *model;

} pdpm_t;

//optimization methods
#define METHOD_NONE   0
#define METHOD_STOCH  1
#define METHOD_AGGLO  2
#define METHOD_GIBBS  3

//Allocate memory and count usage in obj->mem
void *       pdpm_alloc( pdpm_t * obj, unsigned int count, unsigned int size );
//Allocate and zero memory and count usage in obj->mem
void *       pdpm_zalloc( pdpm_t * obj, unsigned int count, unsigned int size );
//Move group grp to cluster cls
void         method_move( pdpm_t * obj, unsigned int grp, unsigned int cls );
//More group grp to cluster cls, return the resulting change in logp
double       method_movep( pdpm_t * obj, unsigned int grp, unsigned int cls );
//Return the index of an empty cluster, or BAD_VCL if none
unsigned int method_free( pdpm_t * obj );
//Move group grp to the best fitting cluster
void         method_best( pdpm_t * obj, unsigned int grp );
//Merge cluster cls1 into cluster cls2
void         method_merge( pdpm_t * obj, unsigned int cls1, unsigned int cls2 );
//Merge cluster cls1 into cluster cls2, return the resulting change in logp
double       method_mergep( pdpm_t * obj, unsigned int cls1, unsigned int cls2 );
//Merge cluster cls1 into cluster cls2, return the resulting change in logp, unmerge
double       method_testmergep( pdpm_t * obj, unsigned int cls1, unsigned int cls2 );
//Gibbs sampler optimization
void         method_gibbs( pdpm_t * obj, int maxiter, double crit );
//stochastic optimization
void         method_stoch( pdpm_t * obj, int maxiter, double crit );
//agglomerative optimization
void         method_agglo( pdpm_t * obj, int maxiter );

typedef struct {

double          s0;   //prior s0 parameter
double        * m0;   //prior m0 parameter
double          a0;   //prior a0 parameter
double          b0;   //prior b0 parameter

//vgr - group membership vector
//Each of the p entries in y has a corresponding integer in vgr
//that indicates group membership. The values in vgr begin
//at zero and are increasing. The largest value is ngr-1.
//y and x are sorted before being passed to the C code such
//that the values in vgr are in order. (length p)
unsigned int  * vgr;

//pgr - number of observations in each group
//The number of observations assigned to each group is stored 
//in pgr. pgr is a vector of length ngr and contains values between 1 and p.
//The group indicator us used to index pgr. Hence, 'pgr[ 0 ]' would give 
//number of observations in group 0. (length ngr)
unsigned int  * pgr;

//pcl - number of observations in each cluster
//pcl is similar to pgr, but for clusters rather than
//groups. (length ngr)
unsigned int  * pcl;

double        * y;    //y vector (array of length p)
double        * x;    //x matrix (array of length p*q)
unsigned int    p;    //nrow(x)
unsigned int    q;    //ncol(x)

//xxgr, xygr, and yygr store, for each group the matrix
//x'x, vector x'y, and scalar y'y, where x is the matrix of 
//covariates of a particular group and y are the dependent 
//observations of the same group. The matrices pointed to by
//xxgr are upper triangular packed storage
double       ** xxgr;   //x'x matrix (array of ngr arrays of length q*(q+1)/2)
double       ** xygr;   //x'y vector (array of ngr arrays of length q)
double        * yygr;   //y'y scalar (array of length ngr)

// xxcl, xycl, and yycl are similar to the variables above,
// but for each cluster rather than each group.
double       ** xxcl;   //x'x matrix (array of at least ncl arrays of length (q*(q+1))/2)
double       ** xycl;   //x'y matrix (array of at least ncl arrays of length q)
double        * yycl;   //y'y scalar (array of at least length ncl)

//s, m, a, and b are temporary storage variables used to hold
//posterior quantities
double        * s;      //storage for an s matrix (array of length q*(q+1)/2)
double        * m;      //storage for an m vector (array of length q)
double          a;      //storage for an a scalar
double          b;      //storage for an b scalar
double          d;      //storage for an d scalar

double        * fbuf;   //temporary storage for fortran routines (array of length q)

} pdpmlm_t;

//Assign the observations/groups according to a simple algorithm
void         pdpmlm_divy( pdpm_t * obj );
//Add an observation/group to a cluster cls
void         pdpmlm_add( pdpm_t * obj, unsigned int grp, unsigned int cls );
//Remove an observation/group from cluster cls
void         pdpmlm_sub( pdpm_t * obj, unsigned int grp, unsigned int cls );
//Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmlm_parm( pdpm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b, double * d );
//Compute part of the log posterior value for a particular cluster
double       pdpmlm_logpcls( pdpm_t * obj, unsigned int cls );
//Compute the log posterior value for the model 
double       pdpmlm_logp( pdpm_t * obj );
//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       pdpmlm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size );

typedef struct {

double          a0;   //prior a0 parameter
double          b0;   //prior b0 parameter

unsigned int  * y;    //y vector (array of length p*q)
unsigned int    q;    //ncol(y)

//gqcl - number of 1s in the q'th column of the submatrix
//of y corresponding to each cluster, corresponds to nkj
//in model notation (length ngr*q)
unsigned int  * gqcl;

//a and b are temporary storage variables used to hold
//posterior quantities
double        * a;      //storage for a vector (length q)
double        * b;      //storage for a vector (length q)

} pdpmbm_t;

//Assign the observations/groups according to a simple algorithm
void         pdpmbm_divy( pdpm_t * obj );
//Add an observation/group to a cluster cls
void         pdpmbm_add( pdpm_t * obj, unsigned int grp, unsigned int cls );
//Remove an observation/group from cluster cls
void         pdpmbm_sub( pdpm_t * obj, unsigned int grp, unsigned int cls );
//Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmbm_parm( pdpm_t * obj, unsigned int cls, double * a, double * b );
//Compute part of the log posterior value for a particular cluster
double       pdpmbm_logpcls( pdpm_t * obj, unsigned int cls );
//Compute the log posterior value for the model 
double       pdpmbm_logp( pdpm_t * obj );
//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       pdpmbm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size );
#endif
