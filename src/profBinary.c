#include "profdpm.h"

SEXP profBinary(SEXP y, SEXP clust, SEXP param, SEXP method,\
  SEXP maxiter, SEXP crit, SEXP verbose) {
  SEXP retval, elem, names, class, dim;
  pdpm_t * obj;
  pdpmbm_t * mdl;
  int i, j, k, cls, elt, onei=1; 
  double *xp, *yp, oned=1.0;
  unsigned int *buf;

  //setup the return value 
  PROTECT(retval = allocVector(VECSXP, 6));
  PROTECT(names = allocVector(STRSXP, 6));
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(names, 0, mkChar("y"));
  SET_STRING_ELT(names, 1, mkChar("param"));
  SET_STRING_ELT(names, 2, mkChar("clust"));
  SET_STRING_ELT(names, 3, mkChar("a"));
  SET_STRING_ELT(names, 4, mkChar("b"));
  SET_STRING_ELT(names, 5, mkChar("logp"));
  SET_STRING_ELT(class, 0, mkChar("profBinary"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, y);
  SET_VECTOR_ELT(retval, 1, param);
  dim = getAttrib(y, R_DimSymbol); 
  SET_VECTOR_ELT(retval, 2, allocVector(INTSXP, INTEGER(dim)[ 0 ]));

  //allocate obj, make assignments, check priors
  obj = (pdpm_t *) R_alloc( 1, sizeof(pdpm_t) );
  obj->mem = sizeof(pdpm_t);

  //allocate model, convenience pointer
  obj->model = (pdpmbm_t *) pdpm_alloc( obj, 1, sizeof(pdpmbm_t) );
  mdl = obj->model;

  //set flags
  obj->flags   = 0;
  if( LOGICAL(verbose)[0] )    { obj->flags |= FLAG_VERBOSE; }

  //set pointers to data
  mdl->y = INTEGER(y);
  mdl->q = INTEGER(dim)[ 1 ];
  obj->ngr      = INTEGER(dim)[ 0 ];

  //check values in param list
  elem       = getListElementByName(param, "lambda");
  if( elem == R_NilValue ) { obj->flags != FLAG_DIRICHL; obj->lam = 0; }
  else if( REAL(elem)[0] < 0 || REAL(elem)[0] > 1 ) {
    warning( "list item \"lambda\" must be between zero and one, using default value" );
    obj->lam = DEFAULT_LAM;
  } else { obj->lam = REAL(elem)[0]; }
  elem       = getListElementByName(param, "alpha");
  if( elem == R_NilValue ) {
    obj->alp = DEFAULT_ALP;
  } else if( REAL(elem)[0] <= 0 ) {
    warning( "list item \"alpha\" must be positive, using default value" );
    obj->alp = DEFAULT_ALP;
  } else { obj->alp = REAL(elem)[0]; }
  elem       = getListElementByName(param, "a0");
  if( elem == R_NilValue ) { 
    mdl->a0 = DEFAULT_BM_A0;
  } else if( REAL(elem)[0] <= 0 ) {
    warning( "list item \"a0\" must be positive, using default value" );
    mdl->a0 = DEFAULT_BM_A0;
  } else { mdl->a0 = REAL(elem)[0]; }
  elem       = getListElementByName(param, "b0");
  if( elem == R_NilValue ) {
    mdl->b0 = DEFAULT_BM_B0;
  } else if( REAL(elem)[0] < 0 ) {
    warning( "list item \"b0\" must be nonnegative, using default value" );
    mdl->b0 = DEFAULT_BM_B0;
  } else { mdl->b0 = REAL(elem)[0]; }

  //set pointers to methods
  obj->divy = &pdpmbm_divy;
  obj->add  = &pdpmbm_add;
  obj->sub  = &pdpmbm_sub;
  obj->logpcls = &pdpmbm_logpcls;
  obj->logp = &pdpmbm_logp;
  obj->logponly = &pdpmbm_logponly;

  //allocate and zero memory ncl, vcl, gcl, gqcl, pbuf
  obj->ncl = 0;
  obj->vcl = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof(unsigned int) );
  obj->gcl = (unsigned int *) pdpm_zalloc( obj, obj->ngr, sizeof(unsigned int) );
  obj->pbuf = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof(unsigned int) ); 
  for( i = 0; i < obj->ngr; i++ ) { obj->vcl[ i ] = BAD_VCL; }
  mdl->gqcl = (unsigned int *) pdpm_zalloc( obj, obj->ngr * mdl->q, sizeof(unsigned int) );

  //9. distribute clusters initially and perform optimization
  if( isInteger(clust) ) {
    for( j = 0; j < obj->ngr; j++ ) {
      obj->add( obj, j, INTEGER(clust)[j] );
    }
  } 
  
  //***********************************************************
  if( INTEGER(method)[0] == METHOD_NONE ) {
    if( isLogical(clust) ) { obj->divy( obj ); }
    obj->logpval = obj->logp( obj );
  }
  else if( INTEGER(method)[0] == METHOD_STOCH ) {
    if( isLogical(clust) ) { obj->divy( obj ); }
    GetRNGstate();
    method_stoch( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
    PutRNGstate();
  }
  else if( INTEGER(method)[0] == METHOD_GIBBS ) {
    if( isLogical(clust) ) { obj->divy( obj ); }
    GetRNGstate();
    method_gibbs( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
    PutRNGstate();
  }
  else if( INTEGER(method)[0] == METHOD_AGGLO ) {
    if( isLogical(clust) ) { for( i = 0; i < obj->ngr; i++ ) { obj->add( obj, i, i ); } }
    method_agglo( obj, INTEGER(maxiter)[0] );
  }
  //***********************************************************

  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
  if( obj->flags & FLAG_VERBOSE ) {
    pdpm_printf( "allocated memory: %fMb\n", obj->mem/1000000.0 );
  }

  //10. complete the return value
  SET_VECTOR_ELT(retval, 3, allocVector(VECSXP, obj->ncl)); //a
  SET_VECTOR_ELT(retval, 4, allocVector(VECSXP, obj->ncl)); //b
  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, 1)); //logp
  REAL(VECTOR_ELT(retval, 5))[0] = obj->logpval;

  for( i = 0; i < obj->ngr; i++ ) { obj->pbuf[ i ] = BAD_VCL; }

  //10.1 remap cluster labels 1--ncl
  cls = 1;
  for( i = 0; i < obj->ngr; i++ ) { 
    if( obj->pbuf[ obj->vcl[ i ] ] == BAD_VCL ) { 
      obj->pbuf[ obj->vcl[ i ] ] = cls++;
    }
    INTEGER(VECTOR_ELT(retval, 2))[i] = obj->pbuf[ obj->vcl[ i ] ];
  }

  cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ cls ] == 0 ) { cls++; }
    elt = obj->pbuf[ cls ] - 1; /* remap cls */
    SET_VECTOR_ELT(VECTOR_ELT(retval, 3), elt, allocVector(REALSXP, mdl->q));
    SET_VECTOR_ELT(VECTOR_ELT(retval, 4), elt, allocVector(REALSXP, mdl->q));
    for( j = 0; j < mdl->q; j++ ) {
      REAL(VECTOR_ELT(VECTOR_ELT(retval, 3), elt))[j] = mdl->a0 +\
         mdl->gqcl[ FMAT(cls, j, obj->ngr) ];
      REAL(VECTOR_ELT(VECTOR_ELT(retval, 4), elt))[j] = mdl->b0 +\
         obj->gcl[ cls ] - mdl->gqcl[ FMAT(cls, j, obj->ngr) ];
    }
    cls++;
  }
  UNPROTECT(3);
  return retval;
}
