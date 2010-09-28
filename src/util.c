#include "profdpm.h"

/* stats */

/* 
    The lfactorial function computes the log (base e) factorial for an integer
    through a series expansion approximation given by Jolly (1961).
    Jolley, L.B.W. (1961) Summation of Series. pp. 28 ISBN 0-486-60023-8 
*/
double lfactorial(unsigned int x) {
    if( x == 0 ) { return 0; }
    return ( LN_SQRT_2PI + (0.5 + x)*log(x) - x );
}

/* 
    Draw an integer at random from 0..(n-1), where the log probability of
    selecting integer k is stored in logp[k]
*/
unsigned int rlcat( double * logp, unsigned int n ) {
    unsigned int i, j, ret = n-1;
    double low=0, high, u;
    u = rlcat_runif(0.0, 1.0);
    for( i = 0; i < n; high=0, i++ ) {
        for( j = 0; j < n; j++ )
            high += exp( logp[ j ] - logp[ i ] );
        high = low + 1/high;
        if( u >= low && u < high ) { 
            ret = i;
            break; 
        }
        low = high;
    }
    return ret;
}

/* From 'Writing R Extenstions' */

SEXP getListElementByName(SEXP list, const char * name) {
    SEXP elem = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    unsigned int i;
    for (i = 0; i < LENGTH(list); i++) {
        if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            elem = VECTOR_ELT(list, i);
            break;
        }
    }
    return(elem);
}
