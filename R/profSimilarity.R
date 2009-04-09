profSimilarity <- function( p1, p2 ) {
  ###################################################
  #do some argument checking
  if( !( "profLinear" %in% is( p1 ) ) ) { stop("p1 not of class \'profLinear\'") }
  if( !( "profLinear" %in% is( p1 ) ) ) { stop("p2 not of class \'profLinear\'") }  
  if( length( p1$clust) != length( p2$clust ) ) { stop("p1 and p2 lengths differ") }
  if( length( p1$group) != length( p2$group ) ) { stop("p1 and p2 lengths differ") }
  if( !all.equal( p1$group, p2$group ) ) { stop("p1 and p2 groups differ") }

  uni1 <- unique( cbind( p1$group, p1$clust ) )[,2]
  uni2 <- unique( cbind( p2$group, p2$clust ) )[,2]
  ret <- .Call("profSimilarity", as.integer(uni1), as.integer(uni2), PACKAGE="profdpm")
  return( ret )
}
