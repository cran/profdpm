#include "profdpm.h"

void pdpmlm_add( pdpm_t * obj, unsigned int grp, unsigned int cls ) {
    pdpmlm_t * mdl = (pdpmlm_t *) obj->model;
    unsigned int i, j, index;
    if( grp >= obj->ngr || cls >= obj->ngr )
        error( "pdpmlm_add: invalid argument: grp = %u cls = %u", grp, cls ); 
    //set vcl, recompute gcl, pcl, and possibly ncl
    obj->vcl[ grp ] = cls;
    if( obj->gcl[ cls ] == 0 ) obj->ncl++;
    mdl->pcl[ cls ] += mdl->pgr[ grp ];
    obj->gcl[ cls ] += 1;
    //allocate and zero memory for xxcl and xycl if necessary
    if( mdl->xxcl[ cls ] == NULL ) {
        mdl->xxcl[ cls ] = (double *) pdpm_zalloc( obj, ( mdl->q * ( mdl->q + 1 ) ) / 2, sizeof(double) );
        mdl->xycl[ cls ] = (double *) pdpm_zalloc( obj, mdl->q, sizeof(double) );
        mdl->yycl[ cls ] = 0.0;
    }
    //(re)compute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
    //use index here rather then UMAT(i, j) to avoid unnecessary index computation
    index = 0;
    mdl->yycl[ cls ] += mdl->yygr[ grp ];
    for( i = 0; i < mdl->q; i++ ) {
        mdl->xycl[ cls ][ i ] += mdl->xygr[ grp ][ i ];
        for( j = i; j < mdl->q; j++ ) {
            mdl->xxcl[ cls ][ index ] += mdl->xxgr[ grp ][ index ];
            index++;
        }
    }
}

void pdpmlm_sub( pdpm_t * obj, unsigned grp, unsigned int cls ) {
    pdpmlm_t * mdl = (pdpmlm_t *) obj->model;
    unsigned int i, j, index;
    if( grp >= obj->ngr || cls >= obj->ngr )
        error( "pdpmlm_sub: invalid argument: grp = %u", grp );  
    //set vcl, recompute gcl, and possibly ncl
    obj->vcl[ grp ] = BAD_VCL;
    mdl->pcl[ cls ] -= mdl->pgr[ grp ];
    obj->gcl[ cls ] -= 1;
    if( obj->gcl[ cls ] == 0 ) obj->ncl--;
    //recompute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
    //use index here rather then UMAT(i, j) to avoid unnecessary index computation
    index = 0;
    mdl->yycl[ cls ] -= mdl->yygr[ grp ];
    for( i = 0; i < mdl->q; i++ ) {
        mdl->xycl[ cls ][ i ] -= mdl->xygr[ grp ][ i ];
        for( j = i; j < mdl->q; j++ ) {
            mdl->xxcl[ cls ][ index ] -= mdl->xxgr[ grp ][ index ];
            index++;
        }
    }
}

void pdpmlm_move( pdpm_t * obj, unsigned int grp, unsigned int cls ) {
    unsigned int old = obj->vcl[ grp ];
    if( old == cls ) return;
    if( old != BAD_VCL )
        pdpmlm_sub( obj, grp, old );
    pdpmlm_add( obj, grp, cls );
}

double pdpmlm_logpcls( pdpm_t * obj, unsigned int cls ) {
    pdpmlm_t * mdl = (pdpmlm_t *) obj->model;
    double logp = 0.0;
    if( obj->gcl[ cls ] == 0 ) { return logp; }
    //compute posterior statistics
    pdpmlm_parm( obj, cls, mdl->s, mdl->m, &mdl->a, &mdl->b, &mdl->d );
    //compute posterior mass
    logp = lgamma( mdl->a / 2 ) - ( mdl->a / 2 ) * log( mdl->b / 2 ) - mdl->d;
    logp += obj->lam * lgamma( obj->gcl[ cls ] );
    return logp;
}


double pdpmlm_logp( pdpm_t * obj ) {
    unsigned int i, cls = 0;
    double logp, accu = 0;
    logp = obj->ncl * log( obj->alp );
    for( i = 0; i < obj->ncl; i++ ) {
        while( obj->gcl[ cls ] == 0 ) { cls++; }
        logp += pdpmlm_logpcls( obj, cls );
        cls++;
    }
    if(obj->flags & FLAG_PREDICT) {
        cls = 0;
        for( i = 0; i < obj->ncl; i++ ) {
            while( obj->gcl[ cls ] == 0 ) { cls++; }
            accu += pow(obj->gcl[ cls ], obj->lam);
            cls++;
        }
        logp += log(obj->alp + accu);
    }        
    return logp;
}

double pdpmlm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size ) {
    unsigned int i, cls;
    double logp, accu;
    logp = obj->ncl * log( obj->alp );
    for( i = 0; i < size; i++ ) {
        logp += pdpmlm_logpcls( obj, only[ i ] );
    }
    if(obj->flags & FLAG_PREDICT) {
        cls = 0;
        for( i = 0; i < obj->ncl; i++ ) {
            while( obj->gcl[ cls ] == 0 ) { cls++; }
            accu += pow(obj->gcl[ cls ], obj->lam);
            cls++;
        }
        logp += log(obj->alp + accu);
    }
    return logp;
}

void pdpmlm_parm( pdpm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b, double * d ) {
    pdpmlm_t * mdl = (pdpmlm_t *) obj->model;
    int i, j, index, ione=1, info, *ipiv;
    double done=1.0, dzero=0.0;
    if( obj->gcl[ cls ] == 0 ) return;
    if( cls >= obj->ngr ) error( "pdpmlm_parm: invalid argument" );

    //load s with s0I + x'x, load m with s0m0 + x'y
    index = 0;
    for( i = 0; i < mdl->q; i++ ) {
        m[ i ] = mdl->s0 * mdl->m0[ i ] + mdl->xycl[ cls ][ i ];
        for( j = 0; j <= i; j++ ) {
            s[ index ] = mdl->xxcl[ cls ][ index ];
            if( j == i ) s[ index ] += mdl->s0;
            index++;
        }
    }

    //m = (s)^(-1) * m
    //dppsv overwrites s, must reload s aftward.
/*
    ipiv = (int *) mdl->fbuf;
    F77_CALL(dppsv)("U", (int *) &mdl->q, &ione, s, m, (int *) &mdl->q, &info);
    if( info < 0 ) error("dppsv: invalid argument");

    if( info > 0 ) { //not numerically pd 
        //load s with s0I + x'x, load m with s0m0 + x'y
        index = 0;
        for( i = 0; i < mdl->q; i++ ) {
            m[ i ] = mdl->s0 * mdl->m0[ i ] + mdl->xycl[ cls ][ i ];
            for( j = i; j < mdl->q; j++ ) {
                s[ index ] = mdl->xxcl[ cls ][ index ];
                if( j == i ) s[ index ] += mdl->s0;
                index++;
            }
        }
*/
        //m = (s)^(-1) * m
        //dspsv overwrites s, must reload s aftward.
        ipiv = (int *) mdl->fbuf;
        F77_CALL(dspsv)("U", (int *) &mdl->q, &ione, s, ipiv, m, (int *) &mdl->q, &info);
        if( info > 0 ) obj->flags |= FLAG_SINGULA;
        if( info < 0 ) error("dspsv: invalid argument");
        //d = 0.5*log|det(s)| (see dspsv/dsptrf documentation)
        *d = 0.0;
        for( i = 0; i < mdl->q; i++ ) {
            if( ipiv[ i ] > 0 ) {
                *d += log( ABS( s[ UMAT(i, i) ] ) );
            } else if( i > 0 && ipiv[ i ] < 0 && ipiv[ i-1 ] == ipiv[ i ] ) {
                *d += log( ABS(
                    s[ UMAT(i-1,i-1) ] * s[ UMAT(i,i) ] -\
                    s[ UMAT(i-1,i) ] * s[ UMAT(i-1,i) ] )\
                );
            }
        }
        *d *= 0.5;
/*
    } else {
        //d = 0.5*log|det(s)| (see dppsv documentation)
        *d = 0.0;
        for( i = 0; i < mdl->q; i++ )
            *d += log( ABS( s[ UMAT(i, i) ] ) );
        *d *= 0.5;
    }
*/  
    //reload s
    index = 0;
    for( i = 0; i < mdl->q; i++ ) {
        for( j = 0; j <= i; j++ ) {
            s[ index ] = mdl->xxcl[ cls ][ index ];
            if( j == i ) { s[ index ] += mdl->s0; }
            index++;
        }
    }

    //b = b0 + y'y + s0*m0'm0 - m'sm
    //b = b0 + y'y
    *b = mdl->b0 + mdl->yycl[ cls ];   
    //b += s0*m0'm0
    *b += mdl->s0*F77_CALL(ddot)( (int *) &mdl->q, mdl->m0, &ione, mdl->m0, &ione);
    //mdl->fbuf = s*m
    //if info > 0, system is singualr, then dont allow 
    //correction for the mean (b -= m'sm). this could 
    //result in negative b! 
    if( info == 0 ) { 
        F77_CALL(dspmv)("U", (int *) &mdl->q, &done, s, m, &ione, &dzero, mdl->fbuf, &ione);
        //b -= m'mdl->fbuf
        *b -= F77_CALL(ddot)( (int *) &mdl->q, m, &ione, mdl->fbuf, &ione );
    } 

    //a = a0 + #values in cluster;
    *a = mdl->a0 + mdl->pcl[ cls ];
}
