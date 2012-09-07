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
/* log(2) */
#define LN_2 0.69314718055994530941
#define rlcat_runif runif

/* stats */
//approximate the log factorial
double       lfactorial(unsigned int x);
//draw a categorical random variate, specifying the
//unnormalized log probabilities 
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
#define DEFAULT_LAM    1.000
#define DEFAULT_ALP    0.001  // 1/1000
#define BAD_VCL UINT_MAX

//pdpmlm defaults
#define DEFAULT_LM_A0     0.001
#define DEFAULT_LM_B0     0.001
#define DEFAULT_LM_M0     0.000
#define DEFAULT_LM_S0     0.001 

//pdpmbm defaults
#define DEFAULT_BM_A0     1.000
#define DEFAULT_BM_B0     1.000

//bit masks for flags (options)
#define FLAG_VERBOSE  1<<0  //should routine be verbose
#define FLAG_OPTCRIT  1<<1  //has optimization criterion been met
#define FLAG_DIRICHL  1<<2  //use cluster Dirichlet prior
#define FLAG_SINGULA  1<<3  //singularities in optimization
#define FLAG_PREDICT  1<<4  //use a special prediction posterior
#define FLAG_SAMPLER  1<<5  //run the Gibbs method in sampler mode
#define FLAG_UNUSED3  1<<6
#define FLAG_UNUSED4  1<<7

//pdpm_t is a generic data type representing the elements of 
//a product partition model
typedef struct pdpm_t {

unsigned char   flags;  //a bitmask of options

double          lam;  //prior lambda parameter
double          alp;  //prior alpha parameter

//ngr - number of groups
//ngr is the number of distinct groups, or observations that should always be
//grouped together, may be a single observation or multiple observations
unsigned int    ngr;  // number of groups

//ncl - number of clusters. 
//Each group is assigned to exactly one cluster. Each cluster is made up of one
//or more groups. ncl is the current number of clusters. ncl must be less than
//or equal to ngr.
unsigned int    ncl;

//vcl - cluster membership vector
//Each of the ngr groups has an entry in vcl indicating cluster membership. The
//group indicator indexes the values in vcl. Hence, 'vcl[ 0 ]' indexes the
//cluster number for group 0. The vcl array should be of length ngr. The values
//in vcl may range from 0 to ngr-1, but may also take the value defined by the
//BAD_CLS macro.
unsigned int  * vcl;

//gcl - number of groups in each cluster
//Hence, gcl[ cls ] is the number of groups in cluster cls The gcl array is of
//length ngr, the number of groups.
unsigned int  * gcl;  

//unnormalized log posterior value
//This value is used to store the log posterior value for the current state,
//unless the value was not computed after a change in state. Good practice 
//dictates that this value should be updated at every state change.
double          logpval;

//temporary storage of length ngr, use by optimization methods 
unsigned int  * pbuf;

//memory usage counter
//This value enumerates the number of bytes allocated by the optimization
//routine using pdpm_alloc.
unsigned int    mem;    

//These function pointers are implemented separately for each type of PPM.

//Assign or reassign group grp to cluster cls
void         (*move)( struct pdpm_t * obj, unsigned int grp, unsigned int cls );

//Compute the log posterior value for the current partition 
double       (*logp)( struct pdpm_t * obj );

//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       (*logponly)( struct pdpm_t * obj, unsigned int * only, unsigned int size );

//PPM specific model struct pointer
void        *model;

} pdpm_t;

//optimization methods
#define METHOD_NONE   0
#define METHOD_STOCH  1
#define METHOD_AGGLO  2
#define METHOD_GIBBS  3
#define METHOD_FAST   4

//These functions operate generically on the PPM data structure pdpm_t

//Allocate memory for an instance of pdpm_t
pdpm_t *     pdpm_init(unsigned int ngr);

//Allocate memory and count usage in obj->mem
void *       pdpm_alloc( pdpm_t * obj, unsigned int count, unsigned int size );

//Allocate and zero memory and count usage in obj->mem
void *       pdpm_zalloc( pdpm_t * obj, unsigned int count, unsigned int size );

//More group grp to cluster cls, return the resulting change in logp
double       method_movep( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Return the index of an empty cluster, or BAD_VCL if none
unsigned int method_free( pdpm_t * obj );

//Move group grp to the best fitting cluster
void         method_best( pdpm_t * obj, unsigned int grp );

//Simple split along vcl[grp1], vcl[grp2]
void         method_ssplit( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls);

//Simple split along vcl[grp1], vcl[grp2], return resulting change in logp
double       method_ssplitp( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls);

//Simple split along vcl[grp1], vcl[grp2], return resulting change in logp, unssplit
double       method_testssplitp( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls);

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

//Fast method, adapted from Wang and Dunson, 2010
void         method_fast( pdpm_t * obj );

//This structure is specific to the linear model PPM.
typedef struct {

double          s0;   //prior s0 parameter
double        * m0;   //prior m0 parameter
double          a0;   //prior a0 parameter
double          b0;   //prior b0 parameter

unsigned int    p;    //nrow(x)
unsigned int    q;    //ncol(x)
double        * y;    //y vector (array of length p)
double        * x;    //x matrix (array of length p*q)

//vgr - group membership vector
//Each of the p entries in y has a corresponding integer in vgr that indicates
//group membership. The values in vgr begin at zero and are increasing. The
//largest value is ngr-1. y and x are sorted in R/profLinear.R before being
//passed to the C code, such that the values in vgr are in order. The vgr
//array should be of length p
unsigned int  * vgr;

//pgr - number of observations in each group
//The number of observations assigned to each group is stored in pgr. pgr is a
//vector of length ngr and contains values between 1 and p. The group indicator
//is used to index pgr. Hence, 'pgr[ 0 ]' would give number of observations in
//group 0.
unsigned int  * pgr;

//pcl - number of observations in each cluster
//pcl is similar to pgr, but for clusters rather than groups. The pcl array
//is also of length ngr.
unsigned int  * pcl;


//xxgr, xygr, and yygr store, for each group the matrix  x'x, vector x'y, and
//scalar y'y, where x is the matrix of covariates of a particular group and y
//are the dependent observations of the same group. The matrices pointed to by
//xxgr are upper triangular packed storage.
double       ** xxgr;   //x'x matrix (array of ngr arrays of length q*(q+1)/2)
double       ** xygr;   //x'y vector (array of ngr arrays of length q)
double        * yygr;   //y'y scalar (array of length ngr)

//xxcl, xycl, and yycl are similar to the variables above, but for each cluster
//rather than each group.
double       ** xxcl;   //x'x matrix (array of at least ncl arrays of length (q*(q+1))/2)
double       ** xycl;   //x'y matrix (array of at least ncl arrays of length q)
double        * yycl;   //y'y scalar (array of at least length ncl)

//s, m, a, and b are temporary storage variables used to hold posterior quantities.
//Ideally, these values are updated each time the model state changes.
double        * s;      //storage for an s matrix (array of length q*(q+1)/2)
double        * m;      //storage for an m vector (array of length q)
double          a;      //storage for an a scalar
double          b;      //storage for an b scalar
double          d;      //storage for an d scalar

double        * fbuf;   //temporary storage for fortran routines (array of length q)

} pdpmlm_t;

//Add an observation/group to a cluster cls
void         pdpmlm_add( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Remove an observation/group from cluster cls
void         pdpmlm_sub( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Move an observation/group to cluster cls
void         pdpmlm_move( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmlm_parm( pdpm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b, double * d );

//Compute part of the log posterior value for a particular cluster
double       pdpmlm_logpcls( pdpm_t * obj, unsigned int cls );

//Compute the log posterior value for the model 
double       pdpmlm_logp( pdpm_t * obj );

//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       pdpmlm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size );


//This structure is specific to the linear model PPM.
typedef struct {

double          a0;   //prior a0 parameter
double          b0;   //prior b0 parameter

unsigned int    q;    //ncol(y)
unsigned int  * y;    //y vector (array of length p*q)

//gqcl - number of 1s in the q'th column of the submatrix
//of y corresponding to each cluster, corresponds to nkj
//in model notation (length ngr*q)
unsigned int  * gqcl;

//a and b are temporary storage variables used to hold posterior quantities.
//Ideally, these values are updates each time the model state changes.
double        * a;      //storage for a vector (length q)
double        * b;      //storage for a vector (length q)

} pdpmbm_t;

//Add an observation/group to a cluster cls
void         pdpmbm_add( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Remove an observation/group from cluster cls
void         pdpmbm_sub( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Assign or reassign group grp to cluster cls
void         pdpmbm_move( pdpm_t * obj, unsigned int grp, unsigned int cls );

//Compute part of the log posterior value for a particular cluster
double       pdpmbm_logpcls( pdpm_t * obj, unsigned int cls );

//Compute the log posterior value for the model 
double       pdpmbm_logp( pdpm_t * obj );

//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       pdpmbm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size );

#endif
