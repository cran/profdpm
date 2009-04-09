#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
   The profSimilarity function computes the Rand statistic
   for two vectors of integers. The computing time of this 
   function increases with the square of the length of the 
   two vectors. The use of bitwise comparrisons
   significantly improves performance.

   Rand, W. (1971). Objective Criteria for the Evaluation
   of Clustering Methods. JASA 66:846-850.
 */

SEXP profSimilarity( SEXP c1, SEXP c2 ) {
  unsigned int i, j, tot = 0, con = 0;
  unsigned int *v1 = (unsigned int *) INTEGER( c1 ), *v2 = (unsigned int *) INTEGER( c2 );
  SEXP ret;
  for( i = 0; i < LENGTH( c1 ) - 1; i++ ) {
    for( j = i + 1; j < LENGTH( c1 ); j++ ) {
      tot++;
      if( !( (v1[ i ] ^ v1[ j ] ? 1 : 0) ^ (v2[ i ] ^ v2[ j ] ? 1 : 0) ) ) { con++; }
    }
  }
  PROTECT( ret = allocVector( REALSXP, 1 ) );
  REAL( ret )[ 0 ] = ( (double) con ) / tot;
  UNPROTECT( 1 );
  return ret;
}
        

