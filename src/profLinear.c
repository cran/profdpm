#include "profdpm.h"


SEXP profLinear(SEXP y, SEXP x, SEXP group, SEXP clust,\
    SEXP param, SEXP method, SEXP maxiter, SEXP crit, SEXP verbose) {
    SEXP retval, elem, names, class, dim;
    pdpm_t * obj;
    pdpmlm_t * mdl;
    int i, j, k, cls, onei=1, gr; 
    double *xp, *yp, oned=1.0;

    // setup the return value 
    PROTECT(retval = allocVector(VECSXP, 10));
    PROTECT(names = allocVector(STRSXP, 10));
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("y"));
    SET_STRING_ELT(names, 1, mkChar("x"));
    SET_STRING_ELT(names, 2, mkChar("group"));
    SET_STRING_ELT(names, 3, mkChar("param"));
    SET_STRING_ELT(names, 4, mkChar("clust"));
    SET_STRING_ELT(names, 5, mkChar("a"));
    SET_STRING_ELT(names, 6, mkChar("b"));
    SET_STRING_ELT(names, 7, mkChar("m"));
    SET_STRING_ELT(names, 8, mkChar("s"));
    SET_STRING_ELT(names, 9, mkChar("logp"));
    SET_STRING_ELT(class, 0, mkChar("profLinear"));
    setAttrib(retval, R_NamesSymbol, names);
    setAttrib(retval, R_ClassSymbol, class);
    SET_VECTOR_ELT(retval, 0, y);
    SET_VECTOR_ELT(retval, 1, x);
    dim = getAttrib(x, R_DimSymbol); 
    SET_VECTOR_ELT(retval, 2, group);
    SET_VECTOR_ELT(retval, 3, param);
    SET_VECTOR_ELT(retval, 4, allocVector(INTSXP, LENGTH(y)));

    //allocate linear PPM structure
    mdl = (pdpmlm_t *) R_alloc( 1, sizeof(pdpmlm_t) );

    //set pointers to data
    mdl->y   = REAL(y);
    mdl->x   = REAL(x);
    mdl->vgr = (unsigned int *) INTEGER(group);
    mdl->p   = INTEGER(dim)[ 0 ];
    mdl->q   = INTEGER(dim)[ 1 ];

    //allocate memory for pgr
    mdl->pgr = (unsigned int *) R_alloc( mdl->p, sizeof(unsigned int) );

    //compute pgr, ngr
    gr = 0;
    for( i = 0; i < mdl->p; i++ ) { mdl->pgr[ i ] = 0; }
    for( i = 0; i < mdl->p; i++ ) { mdl->pgr[ mdl->vgr[ i ] ]++; }
    for( i = 0; i < mdl->p; i++ ) { if( mdl->pgr[ i ] > 0 ) { gr++; } }

    //allocate and initialize generic PPM object
    obj = pdpm_init(gr);
    obj->mem += mdl->p * sizeof(unsigned int) + sizeof(pdpmlm_t);
 
    //set flags
    if( LOGICAL(verbose)[0] )    { obj->flags |= FLAG_VERBOSE; }

    //set pointers to methods
    obj->move     = &pdpmlm_move;
    obj->logp     = &pdpmlm_logp;
    obj->logponly = &pdpmlm_logponly;
    obj->model    = mdl;

    //allocate and zero memory for xxgr xygr, and yygr
    mdl->pcl  = (unsigned int *) pdpm_zalloc( obj, obj->ngr, sizeof(unsigned int) );
    mdl->xxgr = (double **) pdpm_alloc( obj, obj->ngr, sizeof(double *) );
    mdl->xygr = (double **) pdpm_alloc( obj, obj->ngr, sizeof(double *) );
    mdl->yygr = (double *)  pdpm_alloc( obj, obj->ngr, sizeof(double) );
    for( i = 0; i < obj->ngr; i++ ) {
        //xxgr is symmetric packed
        mdl->xxgr[ i ] = (double *) pdpm_zalloc( obj, ( mdl->q * ( mdl->q + 1 ) ) / 2, sizeof(double) );
        mdl->xygr[ i ] = (double *) pdpm_zalloc( obj, mdl->q, sizeof(double) );
        mdl->yygr[i] = 0.0;
    }
  
    //compute xxgr, xygr, yygr
    for( i = 0; i < mdl->p; i++ ) {
        xp = mdl->x + i;
        yp = mdl->y + i; 
        //xxgr += xx' , symmetric packed
        F77_CALL(dspr)("U", (int *) &mdl->q, &oned, xp, (int *) &mdl->p, mdl->xxgr[ mdl->vgr[ i ] ] );
        //xygr += xy
        F77_CALL(daxpy)((int *) &mdl->q, yp, xp, (int *) &mdl->p, mdl->xygr[ mdl->vgr[ i ] ], &onei); 
        //yygr += yy
        mdl->yygr[ mdl->vgr[ i ] ] += (*yp) * (*yp);
    }
  
    //allocate and zero xxcl, xycl, yycl
    mdl->xxcl = (double **) pdpm_zalloc( obj, obj->ngr, sizeof(double *) );
    mdl->xycl = (double **) pdpm_zalloc( obj, obj->ngr, sizeof(double *) );
    mdl->yycl = (double *)  pdpm_zalloc( obj, obj->ngr, sizeof(double) );

    //allocate s, m, fbuf, and pbuf
    //s is symmetric packed
    mdl->s = (double *) pdpm_alloc( obj, ( mdl->q * ( mdl->q + 1 ) ) / 2, sizeof(double) );
    mdl->m = (double *) pdpm_alloc( obj, mdl->q, sizeof(double) );
    mdl->fbuf = (double *) pdpm_alloc( obj, mdl->q, sizeof(double) );


    //check values in param list
    elem = getListElementByName(param, "lambda");
    if( elem == R_NilValue ) { 
        obj->lam = DEFAULT_LAM; 
    } else { obj->lam = REAL(elem)[0]; }
    elem  = getListElementByName(param, "alpha");
    if( elem == R_NilValue ) {
        obj->alp = DEFAULT_ALP;
    } else if( REAL(elem)[0] <= 0 ) {
        warning( "list item \"alpha\" must be positive, using default value" );
        obj->alp = DEFAULT_ALP;
    } else { obj->alp = REAL(elem)[0]; }
    elem = getListElementByName(param, "s0");
    if( elem == R_NilValue ) {
        mdl->s0 = DEFAULT_LM_S0;
    } else if( REAL(elem)[0] <= 0 ) {
        warning( "list item \"s0\" must be positive, using default value" );
        mdl->s0 = DEFAULT_LM_S0;
    } else { mdl->s0 = REAL(elem)[0]; }
    elem = getListElementByName(param, "m0");
    if( elem == R_NilValue ) {
        mdl->m0 = (double *) pdpm_alloc( obj, mdl->q, sizeof(double) ); 
        for( i = 0; i < mdl->q; i++ ) { mdl->m0[i] = DEFAULT_LM_M0; }
    } else if ( LENGTH(elem) < mdl->q ) {
        warning( "list item \"m0\" should be of length ncol(x), using default values" );
        mdl->m0 = (double *) pdpm_alloc( obj, mdl->q, sizeof(double) );
        for( i = 0; i < mdl->q; i++ ) { mdl->m0[i] = DEFAULT_LM_M0; }
    } else { mdl->m0 = REAL(elem); }
    elem = getListElementByName(param, "a0");
    if( elem == R_NilValue ) { 
        mdl->a0 = DEFAULT_LM_A0;
    } else if( REAL(elem)[0] <= 0 ) {
        warning( "list item \"a0\" must be positive, using default value" );
        mdl->a0 = DEFAULT_LM_A0;
    } else { mdl->a0 = REAL(elem)[0]; }
    elem = getListElementByName(param, "b0");
    if( elem == R_NilValue ) {
        mdl->b0 = DEFAULT_LM_B0;
    } else if( REAL(elem)[0] < 0 ) {
        warning( "list item \"b0\" must be nonnegative, using default value" );
        mdl->b0 = DEFAULT_LM_B0;
    } else { mdl->b0 = REAL(elem)[0]; }
 
    //9. distribute clusters initially and perform optimization
    if( isInteger(clust) ) {
        i = 0;
        for( j = 0; j < obj->ngr; j++ ) {
            obj->move( obj, j, INTEGER(clust)[i] );
            i += mdl->pgr[ j ];
        }
    } 
  
    //branch to specific optimization routines
    if( INTEGER(method)[0] == METHOD_NONE ) {
        if( isLogical(clust) ) { method_fast( obj ); }
        obj->logpval = obj->logp( obj );
    } else if( INTEGER(method)[0] == METHOD_STOCH ) {
        if( isLogical(clust) ) { method_fast( obj ); }
        GetRNGstate();
        method_stoch( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
        PutRNGstate();
    } else if( INTEGER(method)[0] == METHOD_GIBBS ) {
        if( isLogical(clust) ) { method_fast( obj ); }
        GetRNGstate();
        method_gibbs( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
        PutRNGstate();
    } else if( INTEGER(method)[0] == METHOD_AGGLO ) {
        if( isLogical(clust) )
            for( i = 0; i < obj->ngr; i++ )
                obj->move( obj, i, i );
        method_agglo( obj, INTEGER(maxiter)[0] );
    } else if( INTEGER(method)[0] == METHOD_FAST ) {
        if(!isLogical(clust))
            warning("\'clust\' argument ignored for \'fast\' method");
        method_fast( obj );
    }

    //check optimization criterion
    if( !(obj->flags & FLAG_OPTCRIT) )
        warning("optimization criterion not met");
    if( obj->flags & FLAG_SINGULA )
        warning("singularities detected during optimization");
    if( obj->flags & FLAG_VERBOSE )
        pdpm_printf( "allocated memory: %fMb\n", obj->mem/1000000.0 );

    //10. complete the return value
    SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, obj->ncl)); //a
    SET_VECTOR_ELT(retval, 6, allocVector(REALSXP, obj->ncl)); //b
    SET_VECTOR_ELT(retval, 7, allocVector(VECSXP, obj->ncl));  //m
    SET_VECTOR_ELT(retval, 8, allocVector(VECSXP, obj->ncl));  //s
    SET_VECTOR_ELT(retval, 9, allocVector(REALSXP, 1));        //logp
    REAL(VECTOR_ELT(retval, 9))[0] = obj->logpval;

    for( i = 0; i < obj->ngr; i++ ) obj->pbuf[ i ] = BAD_VCL;
    cls = 1;
    for( i = 0; i < mdl->p; i++ ) {
        if( obj->pbuf[ obj->vcl[ mdl->vgr[ i ] ] ] == BAD_VCL )
            obj->pbuf[ obj->vcl[ mdl->vgr[ i ] ] ] = cls++;
        INTEGER(VECTOR_ELT(retval, 4))[i] = obj->pbuf[ obj->vcl[ mdl->vgr[ i ] ] ];
    }
    cls = 0;
    for( i = 0; i < obj->ncl; i++) {
        while( obj->gcl[ cls ] == 0 ) { cls++; }
        SET_VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1, allocVector(REALSXP, mdl->q*mdl->q));
        SET_VECTOR_ELT(VECTOR_ELT(retval, 7), obj->pbuf[ cls ]-1, allocVector(REALSXP, mdl->q));
        PROTECT(dim = allocVector(INTSXP, 2));
        INTEGER(dim)[0] = mdl->q;
        INTEGER(dim)[1] = mdl->q;
        setAttrib(VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1), R_DimSymbol, dim);

        pdpmlm_parm( obj, cls, mdl->s,\
            REAL(VECTOR_ELT(VECTOR_ELT(retval, 7), obj->pbuf[ cls ]-1)),\
            REAL(VECTOR_ELT(retval, 5))+obj->pbuf[ cls ]-1,\
            REAL(VECTOR_ELT(retval, 6))+obj->pbuf[ cls ]-1,\
            &mdl->d\
        );

        //copy obj->s (packed) to R matrix (full)
        for( j = 0; j < mdl->q; j++ ) {
            for( k = 0; k < mdl->q; k++ ) {
                //FIXME this is too complicated
                REAL(VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1))[ FMAT(j, k, mdl->q) ] = mdl->s[ j <= k ? UMAT(j, k) : UMAT(k, j) ];
            }
        }
        cls++;
    }
 
    UNPROTECT(3+obj->ncl);
    return retval;
}
