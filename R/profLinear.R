profLinear <- function(y, x, group, clust, param, method="stochastic", 
                       maxiter=1000, crit=1e-5, verbose=FALSE) {
  ###################################################
  #do some argument checking
  if(!is.numeric(y)) { 
    stop("y must be numeric") 
  }
  if(!is.matrix(x) | !is.numeric(x)) { 
    stop("x must be a numeric matrix") 
  }
  if(missing(group)) { 
    group <- seq(1, length(y)) 
  }
  else {
    group <- as.factor(group)
  }
  if(length(y) != length(group)) { 
    stop("length(y) must equal length(group)") 
  }
  if(length(y) != nrow(x)) {
    stop("length(y) must equal nrow(x)") 
  }
  if(missing(clust)) { 
    clust <- FALSE 
  }
  else { 
    clust <- as.factor(clust)
    if(length(y) != length(clust)) { 
      stop("length(y) must equal length(clust)") 
    } 
    for(grp in unique(group)) {
      if( length(unique(clust[group==grp])) > 1 ) { 
        stop("clust and group are conflicting") 
      }
    }
  }
  if(missing(param)) { 
    param <- list(alpha=1,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000) 
  }
  else if(!is.list(param)) {
    warning("param must be a list, using defaults")
    param <- list(alpha=1,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000)
  } 
  else{
    if(length(names(param)) == 0) {
      warning("param argument does not include any named items, using defaults")
      param <- list(alpha=1,a0=0.001,b0=0.001,m0=rep(0,ncol(x)),s0=1.000)
    }
    else if(length(param) > length(names(param))) {
      warning("param contains unnamed items")
    }
  }
  if(!is.character(method)) {
    warning("method must be a character string, using default")
    method <- "stochastic"
  }
  if(!is.numeric(maxiter) | maxiter < 0) {
    warning("maxiter must be numeric and non-negative, using default")
    maxiter <- 1000
  }
  if(!is.numeric(crit) | crit < 0) { 
    warning("crit must be numeric and non-negative, using default") 
    crit <- 1e-5
  }
  if(!is.logical(verbose)) {
    warning("verbose must be a logical, using default")
    verbose <- FALSE
  }

  ###################################################
  #remove missing observations, issue warning
  miss <- apply( is.na( cbind( y, x ) ), 1, any ) 
  ry <- y[!miss]
  rx <- x[!miss,]
  rg <- group[!miss]
  if( is.logical(clust) ) { rc <- FALSE }
  else{ rc <- clust[!miss] }
  if( any( miss ) ) {
    warning( "removed observations with missing values: ", 
      paste(" ", which(miss), sep="") )
  }

  ###################################################
  #order the data according to group
  #convert ordered y to double
  #convert ordered group to integers from 0,1,...
  #convert ordered clust to integers from 0,1,...
  rg <- factor(rg)
  ord <- order(rg)
  ry <- as.double(ry[ord])
  rg <- as.integer(unclass(rg[ord])-1)
  if( !is.logical(rc) ) { rc <- as.integer(unclass(rc[ord])-1) }

  ###################################################
  rx <- as.matrix(x[ord,])

  ###################################################
  #convert method to integer
  if(      method == "none" )          { method <- 0 }
  else if( method == "stochastic" )    { method <- 1 }
  else if( method == "agglomerative" ) { method <- 2 }
  else if( method == "gibbs" )         { method <- 3 }
  else {
    method <- 1 #default is "stochastic"
    warning("method must be \'stochastic\', \'agglomerative\', or \'none\'", )
  }

  ###################################################
  #call the C function
  ret <- .Call("profLinear", ry, rx, rg, rc, as.list(param), as.integer(method),
                as.integer(maxiter), as.double(crit), as.logical(verbose), PACKAGE="profdpm")

  ###################################################
  #undo ordering
  ret$y[ord] <- ret$y
  ret$x[ord,] <- ret$x
  ret$group[ord] <- ret$group
  ret$clust[ord] <- ret$clust
  #ret$clust <- unclass(as.factor(ret$clust))
  #attributes(ret$clust) <- NULL
  return(ret)  
}
