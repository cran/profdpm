profLinear <- function(y, x, group, param, maxiter=1000, crit=1e-5, prior="Dirichlet", verbose=FALSE) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { stop("y must be numeric") }
  if(!is.matrix(x)|!is.numeric(x)) { stop("x must be a numeric matrix") }
  if(missing(group)) { group <- seq(1, length(y)) }
  if(length(y) != length(group)) { stop("length(y) must equal length(group)") }
  if(length(y) != nrow(x)) { stop("length(y) must equal nrow(x)") }
  if(missing(param)) { param <- list(alpha=1,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000) }

  ###################################################
  #remove missing observations, issue warning
  miss <- apply( is.na( cbind( y, x ) ), 1, any ) 
  ry <- y[!miss]
  rx <- x[!miss,]
  rg <- group[!miss]
  if( any( miss ) ) {
    warning( "removed observations with missing values: ", paste(" ", (1:length(y))[miss], sep="") )
  }

  ###################################################
  #order the data according to group
  #convert ordered y to double
  rg <- factor(rg)
  ord <- order(rg)
  ry <- as.double(ry[ord])

  ###################################################
  #transpose ordered x to simplify BLAS/LAPACK calls
  #matrices are double storage by column
  rx <- t(as.matrix(x[ord,]))

  ###################################################
  #convert ordered group to integers from 0,1,...
  rg <- as.integer(unclass(rg[ord])-1)

  ###################################################
  #convert prior to integer
  if( prior == "Dirichlet" ) { prior <- 0 }
  else if( prior == "cluster" ) { prior <- 1 }
  else{ 
    prior <- 0
    warning("prior must be \'Dirichlet\' or \'cluster\', using \'Dirichlet\'")  
  }

  ###################################################
  #call the C function
  ret <- .Call("profLinear", ry, rx, rg, as.list(param), as.integer(maxiter), as.numeric(crit), as.integer(prior), as.logical(verbose), PACKAGE="profdpm")

  ###################################################
  #undo ordering
  ret$y[ord] <- ret$y
  ret$x <- t(ret$x)
  ret$x[ord,] <- ret$x
  ret$group[ord] <- ret$group
  ret$clust[ord] <- ret$clust
  ret$clust <- unclass(as.factor(ret$clust))
  attributes(ret$clust) <- NULL
  return(ret)  
}


