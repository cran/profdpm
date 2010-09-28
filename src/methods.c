#include "profdpm.h"

pdpm_t *     pdpm_init(unsigned int ngr) {
    unsigned int i;
    pdpm_t * obj = (pdpm_t*) R_alloc(1, sizeof(pdpm_t));
    obj->flags = 0;
    obj->mem = sizeof(pdpm_t);
    obj->vcl = pdpm_alloc(obj, ngr, sizeof(unsigned int));
    obj->gcl = pdpm_zalloc(obj, ngr, sizeof(unsigned int));
    obj->pbuf = pdpm_alloc(obj, ngr, sizeof(unsigned int));
    obj->ngr = ngr;
    obj->ncl = 0;
    obj->move     = NULL;
    obj->logp     = NULL;
    obj->logponly = NULL;
    for( i = 0; i < obj->ngr; i++ )
        obj->vcl[ i ] = BAD_VCL;
    return obj;
}

void * pdpm_alloc( pdpm_t * obj, unsigned int count, unsigned int size ) {
    obj->mem += count * size;
    return R_alloc( count, size );
} 

void * pdpm_zalloc( pdpm_t * obj, unsigned int count, unsigned int size ) {
    void * data = R_alloc( count, size );
    char * start = (char *) data;
    char * end = start + count * size;
    do { *(start++) = 0; } while( start < end );
    obj->mem += count * size;
    return data;
} 

double method_movep( pdpm_t * obj, unsigned int grp, unsigned int cls ) {
    double logp = 0.0;
    unsigned int only[2];
    unsigned int old = obj->vcl[ grp ];
    if( old == cls ) return logp;
    only[0] = old;
    only[1] = cls;
    if(obj->logponly) {
        logp -= obj->logponly( obj, only, 2 );
        obj->move( obj, grp, cls );
        logp += obj->logponly( obj, only, 2 );
    } else {
        logp -= obj->logp( obj );
        obj->move( obj, grp, cls );
        logp += obj->logp( obj );
    }
    return logp;
}

unsigned int method_free( pdpm_t * obj ) {
    unsigned int cls = 0;
    while( cls < obj->ngr && obj->gcl[ cls ] > 0 ) cls++;
    if( cls == obj->ngr ) cls = BAD_VCL;
    return cls;
}

void method_best( pdpm_t * obj, unsigned int grp ) {
    unsigned int i, test_cls, best_cls;
    double test_delp=0, best_delp=0;
    best_cls = obj->vcl[ grp ];

    //try and empty cluster 
    if( obj->gcl[ best_cls ] > 1 ) {
        test_cls = method_free( obj );
        test_delp += method_movep( obj, grp, test_cls );
        if( test_delp > best_delp ) { 
            best_delp = test_delp;
            best_cls  = test_cls;
        }
    }
    test_cls = 0;

    //try existing clusters
    for( i = 0; i < obj->ncl; i++ ) {
        while( obj->gcl[ test_cls ] == 0 ) test_cls++;
        test_delp += method_movep( obj, grp, test_cls );
        if( test_delp > best_delp ) { 
            best_delp = test_delp;
            best_cls  = test_cls;
        }
        test_cls++;
    }

    if( obj->vcl[ grp ] != best_cls )
        obj->move( obj, grp, best_cls );
}

void method_ssplit( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls) {
    unsigned int grp = 0, size, cls_old, i;
    //cannot split the same group
    if( grp1 == grp2 ) return;
    //cannot split already different groups
    if( obj->vcl[grp1] != obj->vcl[grp2] ) return;
    //cannot split to same clusters
    if( obj->vcl[grp1] == cls || obj->vcl[grp2] == cls ) return;
    obj->move( obj, grp1, cls );
    cls_old = obj->vcl[ grp2 ];
    size = obj->gcl[ cls_old ];
    //split the remaining members uniformly at random
    for( i = 0; i < size; i++ ) { 
        while( obj->vcl[ grp ] != cls_old ) grp++;
        if( pdpm_runif( 0.0, 1.0 ) <= 0.5 )
            obj->move( obj, grp, cls );
    }
}
  
double method_ssplitp( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls) {
    unsigned int only[2];
    double del;

    //cannot split the same group
    if( grp1 == grp2 ) return 0.0;
    //cannot split already different groups
    if( obj->vcl[grp1] != obj->vcl[grp2] ) return 0.0;
    //cannot split to same clusters
    if( obj->vcl[grp1] == cls || obj->vcl[grp2] == cls ) return 0.0;

    only[0] = obj->vcl[grp1];
    only[1] = cls;
    if(obj->logponly) {
        del -= obj->logponly( obj, only, 1 );
        method_ssplit( obj, grp1, grp2, cls );
        del += obj->logponly( obj, only, 2 );
    } else {
        del -= obj->logp( obj );
        method_ssplit( obj, grp1, grp2, cls );
        del += obj->logp( obj );
    }
    return del;
} 

double method_testssplitp( pdpm_t * obj, unsigned int grp1, unsigned int grp2, unsigned int cls ) {
    unsigned int grp = 0, testgrp, size, cls_old;
    double del = 0.0;

    //cannot split the same group
    if( grp1 == grp2 ) return 0.0;
    //cannot split already different groups
    if( obj->vcl[grp1] != obj->vcl[grp2] ) return 0.0;
    //cannot split to same clusters
    if( obj->vcl[grp1] == cls || obj->vcl[grp2] == cls ) return 0.0;

    //enumerate groups in cls_old (cannot use pbuf elsewhere!!!)
    cls_old = obj->vcl[ grp2 ];
    size = obj->gcl[ cls_old ];
    for( testgrp = 0; testgrp < size; testgrp++ ) {
        while( obj->vcl[ grp ] != cls_old ) grp++;
        obj->pbuf[ testgrp ] = grp++;
    }

    //ssplit && unssplit
    del = method_ssplitp( obj, grp1, grp2, cls );
    for( testgrp = 0; testgrp < size; testgrp++ )
        obj->move( obj, obj->pbuf[ testgrp ], cls_old );

    return del;
}

void method_merge( pdpm_t * obj, unsigned int cls1, unsigned int cls2 ) {
    unsigned int i, grp = 0, size;
    //cannot merge an empty group
    if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) return;
    size = obj->gcl[ cls1 ];
    //merge the cluster
    for( i = 0; i < size; i++ ) {
        while( obj->vcl[ grp ] != cls1 ) grp++;
        obj->move( obj, grp, cls2 );
    }
}    

double method_mergep( pdpm_t * obj, unsigned int cls1, unsigned int cls2 ) {
    unsigned int only[2];
    double del = 0.0;
    //cannot merge an empty group
    if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) return 0.0;
    //merge the cluster
    only[0] = cls2;
    only[1] = cls1;
    if(obj->logponly) {
        del -= obj->logponly( obj, only, 2 );
        method_merge( obj, cls1, cls2 );
        del += obj->logponly( obj, only, 1 );
    } else {
        del -= obj->logp( obj );
        method_merge( obj, cls1, cls2 );
        del += obj->logp( obj );
    }
    return del;
}    

double method_testmergep( pdpm_t * obj, unsigned int cls1, unsigned int cls2 ) {
    unsigned int grp = 0, testgrp, size;
    double del = 0.0;
    //enumerate groups in cls1 (cannot use pbuf elsewhere!!!)
    size = obj->gcl[ cls1 ];
    for( testgrp = 0; testgrp < size; testgrp++ ) {
        while( obj->vcl[ grp ] != cls1 ) { grp++; }
        obj->pbuf[ testgrp ] = grp++;
    }
    //merge clusters 
    del = method_mergep( obj, cls1, cls2 );
    //unmerge clusters
    for( testgrp = 0; testgrp < size; testgrp++ )
        obj->move( obj, obj->pbuf[ testgrp ], cls1 );
    return del;
}

void method_gibbs( pdpm_t * obj, int maxiter, double crit) {
    unsigned int i, *vcl_best, cls, grp, iter = 0;
    unsigned int cls_old, cls_new, *proposal_cls, proposal_ncl;
    double stopcrit = 1.0, logp_best, *proposal_logp;

    //compute initial logp, save initial partition
    obj->logpval = obj->logp( obj );
    logp_best = obj->logpval;
    vcl_best = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
    for( i = 0; i < obj->ngr; i++ ) vcl_best[ i ] = obj->vcl[ i ];
    proposal_logp = (double*) pdpm_alloc( obj, obj->ngr, sizeof( double ) );
    proposal_cls  = (unsigned int*) pdpm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

    while( iter++ < maxiter && stopcrit > crit ) {
        R_CheckUserInterrupt();
        for( grp = 0; grp < obj->ngr; grp++ ) {
            cls_old = obj->vcl[ grp ];
            //compute change in logp for possible moves
            cls = 0;
            proposal_ncl = 0;
            for( i = 0; i < obj->ncl; i++ ) {
                while( obj->gcl[ cls ] == 0 ) { cls++; }
                proposal_logp[ i ] = method_movep( obj, grp, cls );
                proposal_cls[ i ] = cls;
                proposal_ncl++;
                obj->move( obj, grp, cls_old );
                cls++;
            }
            if( obj->gcl[ obj->vcl[ grp ] ] > 1 ) { 
                proposal_cls[ proposal_ncl ] = method_free( obj );
                proposal_logp[ proposal_ncl ] = method_movep( obj, grp, proposal_cls[ proposal_ncl ] );
                proposal_ncl++;
                obj->move( obj, grp, cls_old );
            }
            //draw from the conditional
            cls_new = proposal_cls[ rlcat( proposal_logp, proposal_ncl ) ];
            obj->logpval += method_movep( obj, grp, cls_new );
        }
        //save
        if( obj->logpval > logp_best ) {
            stopcrit += obj->logpval - logp_best;
            logp_best = obj->logpval;
            for( i = 0; i < obj->ngr; i++ ) vcl_best[ i ] = obj->vcl[ i ];
        }
        //print summary if requested every iteration
        if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) 
            pdpm_printf("iter: %u, ncl: %u, logp: %f, best: %f, crit: %f\n",\
                     iter, obj->ncl, obj->logpval, logp_best, stopcrit );
    
        //update stopping criterion
        if(stopcrit != DBL_MAX) stopcrit *= 0.95;
    }
    if(stopcrit <= crit) obj->flags |= FLAG_OPTCRIT;
    //restore
    obj->logpval = logp_best;
    for( grp = 0; grp < obj->ngr; grp++ )
        obj->move( obj, grp, vcl_best[ grp ] );  
}

void method_stoch( pdpm_t * obj, int maxiter, double crit) {
    unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0;
    double logp_old, pdel = 1, pcum = 0;

    //allocate memory for vcl_old, grps
    vcl_old = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
    grps    = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

    //compute initial logp
    obj->logpval = obj->logp( obj );
    while( iter++ < maxiter ) {
        R_CheckUserInterrupt();

        //select the number of groups to shuffle
        ngrps = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0, 1.0 ) );
        ngrps = ngrps == 0 ? 1 : ngrps;

        //randomly select ngrps groups to shuffle, save indicators
        for( i = 0; i < ngrps; i++ ) {
            grps[ i ] = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0, 1.0 ) );
            vcl_old[ i ] = obj->vcl[ grps[ i ] ];
        }

        //compute old logp, move groups to random cluster 
        logp_old = obj->logpval;
        cls = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0, 1.0 ) );
        for( i = 0; i < ngrps; i++ ) obj->move( obj, grps[ i ], cls );
        obj->logpval = obj->logp( obj );

        if( obj->logpval <= logp_old ) {  
            //move each group to best cluster
            for( i = 0; i < ngrps; i++ ) method_best( obj, grps[ i ] );
            //compute logp, keep new clustering if better, else revert to old
            obj->logpval = obj->logp( obj );
            if( obj->logpval <= logp_old ) {    
                for( i = 0; i < ngrps; i++ ) obj->move( obj, grps[ i ], vcl_old[ i ] );
                pdel *= 0.9;
                obj->logpval = logp_old;
            }
        }

        //update the stopping criterion
        else { 
            pdel = 0.5 * (obj->logpval - logp_old) + 0.5 * pdel;
            logp_old = obj->logpval;
        }
        pcum += pdel;

        //print summary if requested every 20 iterations
        if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 )
            pdpm_printf("iter: %u, ncl: %u, logp: %f, exp: %u, crit: %f\n",\
                iter, obj->ncl, logp_old, ngrps, pdel / pcum );

        //check stopping criterion, break the optimization loop if met
        if( pcum > 0 && ( pdel / pcum ) < crit ) {
            obj->flags |= FLAG_OPTCRIT;
            break;
        }
    }
}

void method_agglo( pdpm_t * obj, int maxiter ) {
    double *del, del_best, logp_best = -DBL_MAX;
    unsigned int *vcl_best, i, j, index;
    unsigned int icls, jcls, icls_best = BAD_VCL, jcls_best = BAD_VCL;
    unsigned int icls_last = BAD_VCL, jcls_last = BAD_VCL;
    int calcs = 0, cent;

    //allocate some additional memory
    //del[ i, j ] - upper triangular packed storage
    del = (double *) pdpm_alloc( obj, ( obj->ngr * ( obj->ngr + 1 ) / 2 ), sizeof( double ) );
    vcl_best = (unsigned int *) pdpm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

    //compute initial logp 
    obj->logpval = obj->logp( obj );
    //compute required calculations (and 1%)
    cent = obj->ngr * ( obj->ngr - 1 ) + 1;
    cent = cent > 100 ? cent / 100 : 1;

    //repeat until all clusters are merged into one
    //while( obj->ncl > 1 && maxiter-- != 0 ) {
    while( obj->ncl > 1 ) {
        R_CheckUserInterrupt();
        //compute best merge
        del_best = -DBL_MAX;
        icls = 0;
        for( i = 0; i < obj->ncl - 1; i++ ) {
            while( obj->gcl[ icls ] == 0 ) icls++;
            jcls = icls + 1;
            for( j = 0; j < obj->ncl - i - 1; j++ ) {
                while( obj->gcl[ jcls ] == 0 ) jcls++;
                index = UMAT(icls, jcls);
                if( obj->ncl == obj->ngr || icls == jcls_last || jcls == jcls_last ) { 
                    del[ index ] = method_testmergep( obj, icls, jcls ); 
                    calcs++;
                    if( (obj->flags & FLAG_VERBOSE) && (calcs % cent == 0) )
                        pdpm_printf("\rpercent complete: %d%", calcs / cent);
                }
                if( del[ index ] >= del_best ) {
                    del_best = del[ index ];
                    icls_best = icls;
                    jcls_best = jcls;
                }
                jcls++;
            }
            icls++;
        }

        //merge
        method_merge( obj, icls_best, jcls_best );
        icls_last = icls_best;
        jcls_last = jcls_best;
        obj->logpval += del_best;

        //save
        if( obj->logpval > logp_best ) {
            logp_best = obj->logpval;
            for( i = 0; i < obj->ngr; i++ ) vcl_best[ i ] = obj->vcl[ i ];
        }  
    }

    //restore
    obj->logpval = logp_best;
    obj->flags |= FLAG_OPTCRIT;
    for( i = 0; i < obj->ngr; i++ )
        obj->move( obj, i, vcl_best[ i ] );

    if( obj->flags & FLAG_VERBOSE )
        pdpm_printf("\rpercent complete: 100%\nncl: %u logp: %f\n",\
            obj->ncl, obj->logpval);
}

/* Fast method, adapted from SUGS, Wang & Dunson (2010) */
void method_fast(pdpm_t * obj) {
    //strategy: 
    // 1. Select an ordering at random
    // 2. Assign first observation to first cluster
    // 3. Assign subsequent observations to new or existing clusters
    //    such that conditional posterior is largest.
    unsigned int tmp, swp, grp, *grp_list;
    unsigned int cls, cls_best, ncl;
    double del, accu, cent, hund;
    accu = 0.0;
    cent = 100.0/obj->ngr;
    hund = 1.0;

    grp_list = pdpm_alloc( obj, obj->ngr, sizeof(unsigned int) );
    for( grp = 0; grp < obj->ngr; grp++ )
        grp_list[ grp ] = grp;
    for( grp = 0; grp < obj->ngr; grp++ ) {
        swp = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0, 1.0 ) );
        tmp = grp_list[ grp ];
        grp_list[ grp ] = grp_list[ swp ];
        grp_list[ swp ] = tmp;
    }

    obj->move( obj, grp_list[ 0 ], 0 );
    obj->logpval = obj->logp( obj );
    for( grp = 1; grp < obj->ngr; grp++ ) {
        R_CheckUserInterrupt();
        obj->move( obj, grp_list[ grp ], 0 );
        cls_best = 0;
        ncl = obj->ncl;
        for( cls = 0; cls < ncl + 1; cls++ ) {
            del = method_movep( obj, grp_list[ grp ], cls );
            if(del > 0) cls_best = cls;
        }
        obj->move( obj, grp_list[ grp ], cls_best );
        accu += cent;
        if( (obj->flags & FLAG_VERBOSE) && ( accu > hund ) ) {
            pdpm_printf("\rpercent complete: %d%", (int)hund); 
            hund += 1.0;
        }
    }
    if( obj->flags & FLAG_VERBOSE )
        pdpm_printf("\rpercent complete: 100%\n"); 
    obj->logpval = obj->logp( obj );
    obj->flags |= FLAG_OPTCRIT;
}

/* method_spmer (simple and restricted split-merge)
   Jain S and Neil R (2004) A Split-Merge Markov Chain Monte 
   Carlo Procedure for the Dirichlet Process Mixture Model.
   JCGS 13(1) 158-182. 
   
void method_spmer(pdpm_t * obj, int maxiter, double crit, int simple) {
    //strategy: iterate over the following steps:
    // 1. select grp1, grp2, uniformly at random from {0, ngr-1}
    // 2. split if(grp1 == grp2) or merge if(grp1 != grp2) (simple or restricted)
    // 3. evaluate MH acceptance probability, accept or regect move

    unsigned int grp1, grp2, i, split = 0, n1, n2, cls_old, cls_new;
    //unsigned int testgrp, grp, size;
    double logp_old, del, mhlogp;
    logp_old = obj->logp( obj );
    for( i = 0; i < maxiter; i++ ) {
        grp1 = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0 , 1.0 ) );
        do grp2 = (unsigned int) floor( obj->ngr * pdpm_runif( 0.0 , 1.0 ) );
        while( grp2 == grp1 && (obj->ngr > 1) );
        split = ( obj->vcl[ grp1 ] == obj->vcl[ grp2 ] );

        //store group data
        //cls_old = obj->vcl[ grp1 ];
        //size = obj->gcl[ cls_old ];
        //for( testgrp = 0; testgrp < size; testgrp++ ) {
        //    while( obj->vcl[ grp ] != cls_old ) grp++;
        //    obj->pbuf[ testgrp ] = grp++;
        //}

        if( split ) {
            cls_old = obj->vcl[ grp1 ];
            cls_new = method_free( obj );
            del = simple ? method_testssplitp( obj, grp1, grp2, cls_new )\
                         : method_testrsplitp( obj, grp1, grp2, cls_new );
            n1 = obj->gcl[ cls_old ];
            n2 = obj->gcl[ cls_new ];
            mhlogp = simple ? del + (n1 + n2 - 2) * LN_2\
                            : del + (n1 + n2 - 2) * LN_2; //FIXME need restricted mhlogp here
        } else {
            n1 = obj->gcl[ obj->vcl[ grp1 ];
            n2 = obj->gcl[ obj->vcl[ grp2 ];
            del = method_testmergep( obj, obj->vcl[ grp1 ], obj->vcl[ grp2 ] );
            mhlogp = del - (n1 + n2 - 2) * LN_2;
        }

        //accept MH move
        if( mklogp >= 0 || mklogp >= log( pdpm_runif( 0.0, 1.0 ) ) ) {
            if( split ) {
                if( simple ) method_ssplit( obj, grp1, grp2, cls_new );
                else method_rsplit( obj, grp1, grp2, cls_new );
            } else {
                method_merge( obj, obj->vcl[ grp1 ], obj->vcl[ grp2 ] );
            }
            obj->logp += del;
        }
    }
}

*/
