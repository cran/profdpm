#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
    pci - Partition Comparison Indices
    The pci function computes some statstics used to compare two cluster
    partitions. The statistics are the Rand index, Wallace coefficients, and 
    the Fowlkes and Mallows index.
 
    1. Rand, W. (1971). Objective Criteria for the Evaluation of Clustering
       Methods. JASA 66:846-850.
    2. Fowlkes, E. B. and Mallows, C. L. (1983). A Method for Comparing Two
       Hierarchical Clusterings. JASA 78:553-569.
    3. Wallace, D. L. (1983). A Method for Comparing Two Hierarchical
       Clusterings: Comment. JASA 78:569-576. 
*/

SEXP pci( SEXP c1, SEXP c2 ) {
    //n11 - number of pairs in same cluster in c1 and c2
    //n00 - number of pairt in different clusters in c1 and c2
    //n10 - number of pairs in same cluster in c1 but different in c2
    //n01 - number of pairs in same cluster in c2 but different in c1
    double n11=0.0, n00=0.0, n10=0.0, n01=0.0;
    unsigned int i, j, n;
    unsigned int *v1 = (unsigned int *) INTEGER( c1 );
    unsigned int *v2 = (unsigned int *) INTEGER( c2 );
    SEXP ret, names;
    n = (unsigned int) LENGTH( c1 );
    for( i = 0; i < n - 1; i++ ) {
        for( j = i + 1; j < n; j++ ) {
            if( v1[ i ] == v1[ j ] ) {
                if( v2[ i ] == v2[ j ] ) { n11++; }
                else { n10++; }
            } 
            else {
                if( v2[ i ] != v2[ j ] ) { n00++; }
                else { n01++; }
            }
        }
    }
    PROTECT( ret = allocVector( REALSXP, 5 ) );
    PROTECT( names = allocVector( STRSXP, 5 ) );
    SET_STRING_ELT(names, 0, mkChar("R"));   //Rand
    SET_STRING_ELT(names, 1, mkChar("FM"));  //Fowlkes and Mallows
    SET_STRING_ELT(names, 2, mkChar("W10")); //Wallace 10
    SET_STRING_ELT(names, 3, mkChar("W01")); //Wallace 01
    SET_STRING_ELT(names, 4, mkChar("J"));   //Jaccard
    setAttrib(ret, R_NamesSymbol, names);
    REAL( ret )[ 0 ] =  (n11+n00)/(n11+n00+n10+n01);
    REAL( ret )[ 1 ] = n11 == 0.0 ? n11 : n11/sqrt((n11+n01)*(n11+n10));
    REAL( ret )[ 2 ] = n11 == 0.0 ? n11 : n11/(n11+n10);
    REAL( ret )[ 3 ] = n11 == 0.0 ? n11 : n11/(n11+n01);
    REAL( ret )[ 4 ] = n11 == 0.0 ? n11 : n11/(n11+n01+n10);
    UNPROTECT( 2 );
    return ret;
}
        

