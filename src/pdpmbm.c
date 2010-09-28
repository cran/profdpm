#include "profdpm.h"

void pdpmbm_add( pdpm_t * obj, unsigned int grp, unsigned int cls ) {
    pdpmbm_t * mdl = (pdpmbm_t *) obj->model;
    unsigned int i;
    if( grp >= obj->ngr || cls >= obj->ngr )
        error( "pdpmb_add: invalid argument: grp = %u cls = %u", grp, cls ); 
    //set vcl, recompute gcl, and possibly ncl
    obj->vcl[ grp ] = cls;
    if( obj->gcl[ cls ] == 0 ) obj->ncl++;
    obj->gcl[ cls ] += 1;
    //(re)compute gqcl
    for( i = 0; i < mdl->q; i++ )
        mdl->gqcl[ FMAT(cls, i, obj->ngr) ] += mdl->y[ FMAT(grp, i, obj->ngr) ];
}

void pdpmbm_sub( pdpm_t * obj, unsigned grp, unsigned int cls ) {
    pdpmbm_t * mdl = (pdpmbm_t *) obj->model;
    unsigned int i;
    if( grp >= obj->ngr || cls >= obj->ngr )
        error( "pdpmb_sub: invalid argument: grp = %u", grp );
    //set vcl, recompute gcl, and possibly ncl
    obj->vcl[ grp ] = BAD_VCL;
    obj->gcl[ cls ] -= 1;
    if( obj->gcl[ cls ] == 0 ) { obj->ncl--; }
    //recompute gqcl
    for( i = 0; i < mdl->q; i++ )
        mdl->gqcl[ FMAT(cls, i, obj->ngr) ] -= mdl->y[ FMAT(grp, i, obj->ngr) ];
}

void pdpmbm_move( pdpm_t * obj, unsigned int grp, unsigned int cls ) {
    unsigned int old = obj->vcl[ grp ];
    if( old == cls ) return;
    if( old != BAD_VCL )
        pdpmbm_sub( obj, grp, old );
    pdpmbm_add( obj, grp, cls );
}

double pdpmbm_logpcls( pdpm_t * obj, unsigned int cls ) {
    pdpmbm_t * mdl = (pdpmbm_t *) obj->model;
    unsigned int i;
    double logp = 0.0;
    if( obj->gcl[ cls ] == 0 ) return logp;
    //compute posterior mass
    for( i = 0; i < mdl->q; i++ ) {
        logp += lgamma( mdl->a0 + (double) mdl->gqcl[ FMAT(cls, i, obj->ngr) ] ) +\
                lgamma( mdl->b0 + (double) obj->gcl[ cls ] -\
                (double) mdl->gqcl[ FMAT(cls, i, obj->ngr) ] ) -\
                lgamma( (double) obj->gcl[ cls ] + mdl->a0 + mdl->b0 );
    }
    if( obj->flags & FLAG_DIRICHL ) { logp += lgamma( obj->gcl[ cls ] ); }
    else { logp += obj->lam * lgamma( obj->gcl[ cls ] + 1 ); }
    return logp;
}

double pdpmbm_logp( pdpm_t * obj ) {
    unsigned int i, cls = 0;
    double logp;
    logp = obj->ncl * log( obj->alp );
    for( i = 0; i < obj->ncl; i++ ) {
        while( obj->gcl[ cls ] == 0 ) { cls++; }
        logp += pdpmbm_logpcls( obj, cls );
        cls++;
    }
    return logp;
}

double pdpmbm_logponly( pdpm_t * obj, unsigned int * only, unsigned int size ) {
    unsigned int i;
    double logp;
    logp = obj->ncl * log( obj->alp );
    for( i = 0; i < size; i++ ) {
        logp += pdpmbm_logpcls( obj, only[ i ] );
    }
    return logp;
}
