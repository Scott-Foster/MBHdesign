# This is package MBHdesign 

"alterInclProbs" <- function (legacy.sites, potential.sites = NULL, n = NULL, inclusion.probs = NULL, 
    mc.cores = 1, sigma = NULL) {
    if (is.null(n) & is.null(inclusion.probs)) 
        stop("One of the arguments n or inclusion.probs should be non-NULL")
    if (!is.null(n) & !is.null(inclusion.probs)) 
        warning("\n\tBoth n and inclusion.probs specified!\n\tUsing inclusion.probs but scaling them to sum to n")
    dimen <- ncol(legacy.sites)
    if (is.null(potential.sites)) {
        N <- 100
        potential.sites <- getGridInHull(legacy.sites, N = N)
        colnames(potential.sites) <- colnames(legacy.sites)
    }
    N <- nrow(potential.sites)
    if (ncol(legacy.sites) != ncol(potential.sites)) 
        stop("Dimension of legacy points does not match dimension of potential sites")
    if (is.null(inclusion.probs)) {
        message("No inclusion.probs supplied, assuming uniform")
        inclusion.probs <- n * rep(1/N, N)
    }
    if (length(legacy.sites) == 0) {
        warning("No legacy sites provided -- nothing to adjust for.  Returning inclusion probabilities.")
        return(inclusion.probs)
    }
    distFun4 <- function(ii) {
        lpt <- legacy.sites[ii, ]
        disty <- sweep(potential.sites, 2, FUN = "-", STATS = lpt)
        disty <- disty^2
        disty <- rowSums(disty, na.rm=TRUE)
        disty <- sqrt(disty)
        return(disty)
    }
    if( mc.cores >1){
      cl <- parallel::makeCluster(mc.cores)
      parallel::clusterExport(cl, c("potential.sites", "legacy.sites", 
	  "distFun4"), envir = environment())
      tmp <- parallel::parLapply(cl, 1:nrow(legacy.sites), distFun4)
      parallel::stopCluster(cl)
    }
    else
      tmp <- lapply( 1:nrow( legacy.sites), distFun4)
    tmp <- do.call("cbind", tmp)
    if( !is.null( n))
      inclusion.probs <- n * inclusion.probs / sum( inclusion.probs, na.rm=TRUE)
    else
      n <- sum(inclusion.probs, na.rm=TRUE)
    if (is.null(sigma)) 
        sigma <- find.sigma(n = n, nL = nrow(legacy.sites), potSites = potential.sites)
    legacyLocs <- apply(tmp, 2, which.min)
    d <- exp(-(tmp^2)/(2 * sigma * sigma))
    adj.probs.ret <- inclusion.probs
    for (jj in 1:ncol(tmp)) {
        adj <- (1 - inclusion.probs[legacyLocs[jj]])
        adj.probs.ret <- adj.probs.ret - adj * adj.probs.ret * 
            d[, jj]
    }
    tmp1 <- n * adj.probs.ret/sum(adj.probs.ret, na.rm=TRUE)
    return(tmp1)
}

{#"calcExactInclusionProbs2" <- function( potential.sites, sampled.sites, legacy.sites=NULL, 
#                                        inclusion.probs=NULL, times=1e7, mc.cores=4) {
#  #inclusion probs should be unaltered probs?
#  n <- nrow( sampled.sites)
#  N <- nrow( potential.sites)
#  dimension <- ncol( potential.sites)
#  if( ncol( sampled.sites) != dimension)
#    stop( "The dimensions of the potential sites and the sampled sites are different! Please check. Only pass those columns of each matrix that are coordinates")
#  if( is.null( inclusion.probs)){
#    message( "No inclusion.probs supplied, assuming uniform")
#    inclusion.probs <- n*rep( 1/N, N) #even probs
#  }
#  if( n != sum( inclusion.probs)){
#    message( "Sum of inclusion probabilities not equal to the number of sampled.sites. Rescaling inclusion probabilities")
#    inclusion.probs <- inclusion.probs * n / sum( inclusion.probs)
#  }  
#  if( is.null( legacy.sites)){
#    message( "No legacy sites supplied")
#    all.sampled.sites <- sampled.sites
#  }
#  else{
#    all.sampled.sites <- cbind( sampled.sites, legacy.sites)
#    L <- nrow( legacy.sites)
#  }
#  #Whittle down the potential.sites to include only the all.sampled.sites and their immediate M neighbours
#  M <- 100
#  
#  funny <- function(x){
#    #this used to be qasuiSamp1()  !?!?!
#    tmp <- quasiSamp1( n=n, dimension=dimension, potential.sites=potential.sites, inclusion.probs=inclusion.probs, return.potential.sites=FALSE)
#    return( tmp$ID)
#  }
#  cl <- parallel::makeCluster( mc.cores)
#  tmppy <- .libPaths()
#  parallel::clusterExport( cl, c( "dimension", "potential.sites", "inclusion.probs", "funny", "quasiSamp", "tmppy"), envir=environment())
#  parallel::clusterEvalQ( cl, {.libPaths( tmppy); library( randtoolbox); library( class)})
#  tmp <- parallel::parLapply( cl, 1:times, funny)
#  parallel::stopCluster( cl)
#
#  tmp <- factor( unlist( tmp), levels=1:nrow( potential.sites))
#  tmp <- table( tmp) / times
#  
#  return( tmp)
#}
}

"find.sigma" <- function( n, nL, potSites) {
  dim <- ncol( potSites)
  hully <- geometry::convhulln( potSites, options="FA")
  A <- hully$vol
  A.dashed <- A / (n+nL)
  r <- ((A.dashed*gamma( dim/2 + 1))/(pi^(dim/2)))^(1/dim)
#  r <- sqrt( A.dashed/pi)
  alpha <- 0.95
  
  fun <- function(sig)
    stats::pnorm(r, mean=0, sd=sig) - stats::pnorm( -r, mean=0, sd=sig) - alpha
  sigma <- stats::uniroot( fun, interval=c(1e-8,r))$root

  return( sigma)
}

"getControl" <- function( contr) {
  if( !"k" %in% names( contr))
    contr$k <- 3
  if( !"N" %in% names( contr))
    contr$N <- 100
  if( !"B" %in% names( contr))
    contr$B <- 1e3
  if( !"mc.cores" %in% names( contr))
    contr$mc.cores <- 1
    
  return( contr)
}

"getGAMpred" <- function( fm, predPts, B=1e3, seMethod="parameter", control=NULL) {
	Xp <- mgcv::predict.gam( fm, newdata=as.data.frame( predPts), type='lpmatrix')
	if( seMethod=="parameter")
  	betas <- mvtnorm::rmvnorm( n=B, mean=coef( fm), sigma=fm$Vp)
  if( seMethod=="BayesianBoot"){
    control <- getControl( control)
    funny <- function(x){
      wts <- stats::rexp( length(fm$y), 1)
      wts <- wts / sum( wts)
      fm1 <- stats::update( fm, weights=wts)
      return( stats::coef( fm1))
    }
    betas <- do.call( "rbind", lapply( 1:B, funny))
  }    
	tmp <- try( etas <- Xp %*% t( betas), silent=FALSE)
	if( inherits( tmp, "try-error"))
	  stop("Crash. Probably due to matrix of prediction sites by bootstrap samples being too large for your machine.  Try reducing either or both")
	preds <- fm$family$linkinv( etas)
  preds <- colMeans( preds)
  ret <- list( mean=mean( preds), se=stats::sd( preds), CI=stats::quantile( preds, c(0.025, 0.975)))
  return( ret)
}

"getGridInHull" <- function( locations, N=100) {
  hully <- geometry::convhulln( locations)
	predPts <- data.frame( tmp=rep( NA, N))
	for( ii in 1:ncol( locations))
	  predPts <- cbind( predPts, seq( from=min( locations[,ii], na.rm=TRUE), to=max( locations[,ii], na.rm=TRUE), length=N))
  predPts <- predPts[,-1]
  if( is.null( colnames( locations)))
    colnames( locations) <- paste0( "Var", 1:ncol( locations))
  colnames( predPts) <- colnames( locations)
  predPts <- expand.grid( predPts)
  phull2 <- geometry::convhulln(rbind(predPts,locations))
  nrp <- nrow( locations)
  nrx <- nrow( predPts)
  outside <- unique( phull2[phull2>nrp])-nrp
  done <- FALSE
  while(!done){
    phull3 <- geometry::convhulln( rbind( locations, predPts[-(outside),]))
    also.outside <- (1:nrx)[-outside][unique(phull3[phull3>nrp])-nrp]
    outside <- c(outside,also.outside)
    done <- length(also.outside)==0
  }
  predPts <- predPts[-(outside),]
  return( predPts)
}

"modEsti" <- function( y, locations, includeLegacyLocation=TRUE, legacyIDs=NULL, predPts=NULL, 
                       family=stats::gaussian(), offset=rep(0,length(y)), control=list()) {
  control <- getControl( control)
  locations <- as.matrix( locations)
	if( is.null( predPts)){ #if not provided then make up a grid
	  message( "No prediction points provided. Using a dense grid on the convex hull of the sampling locations.")
    predPts <- getGridInHull( locations, N=control$N)
  }
  formPart1 <- "y~"
  formPart2 <- paste0("te(",paste( colnames( locations), collapse=','),",k=control$k)")
  my.form <- paste0(formPart1, formPart2)
  if( includeLegacyLocation){
    if( is.null( legacyIDs))
      stop("need to provide the rownumbers (of locations argument) that are legacy sites in the legacyIDs argument")
    formPart3 <- paste0("+combinedDistToLegacy")#paste0("+s(distToNearLegacy,k=control$k,bs='cr')")
    my.form <- paste0(my.form, formPart3)  
    disty <- as.matrix( stats::dist( locations, diag=TRUE, upper=TRUE)) #all locations
    disty <- disty[,legacyIDs]
    sigma <- find.sigma( n=nrow( locations)-length( legacyIDs), nL=length( legacyIDs), potSites=predPts)
    combinedDistToLegacy <- rowSums( exp( -(disty^2)/(sigma*sigma) ) )
  }
  else
    combinedDistToLegacy <- rep( NA, length( y)) #distToNearLegacy <- rep( NA, length( y))
  my.form <- stats::as.formula( my.form)
  
  mod.dat <- as.data.frame( cbind( y, locations, combinedDistToLegacy))
	fm <- mgcv::gam( my.form, data=mod.dat, family=family, offset=offset)
  #####getting predictions
  #first find distance to legacy points, if needed.
  if( includeLegacyLocation){
    distFun4 <- function( ii){#, potential.sites, legacy.sites){
      lpt <- locations[legacyIDs[ii],]
      disty <- sweep( predPts, 2, FUN='-', STATS=lpt) ####Could have nasty implications when dealing with 3+ Dims
      disty <- disty^2
      disty <- rowSums( disty)
      disty <- sqrt( disty)
      return( disty)
    }
    if( control$mc.cores > 1){
      cl <- parallel::makeCluster( control$mc.cores)
      parallel::clusterExport( cl, c("predPts", "locations", "legacyIDs", "distFun4"), envir=environment())
      tmp <- parallel::parLapply( cl, 1:length( legacyIDs), distFun4)
      parallel::stopCluster( cl)
    }
    else
      tmp <- lapply( 1:length( legacyIDs), distFun4)
    tmp1 <- do.call( "cbind", tmp)
    tmp1 <- exp( -(tmp1^2)/(sigma^2))
    combinedDistToLegacy <- rowSums( tmp1)
    updatedPredPts <- cbind( predPts, combinedDistToLegacy)
  }
  else
    updatedPredPts <- predPts
 	Xp <- mgcv::predict.gam( fm, newdata=as.data.frame( updatedPredPts), type='lpmatrix')
 	betas <- mvtnorm::rmvnorm( n=control$B, mean=coef( fm), sigma=fm$Vp)
	tmp <- try( etas <- Xp %*% t( betas), silent=FALSE)
	if( inherits( tmp, "try-error"))
	  stop("Crash. Probably due to matrix of prediction sites by B samples being too large for your machine.  Try reducing either or both (or getting more memory on your machine)")
	preds <- fm$family$linkinv( etas)
  preds <- colMeans( preds)
  ret <- list( mean=mean( preds), se=stats::sd( preds), CI=quantile( preds, c(0.025, 0.975)))

  return( ret)
}

"setDesignParams" <- function( dimension, study.area, potential.sites, inclusion.probs, n) {
  if( is.null( study.area)){
    if( is.null( potential.sites)){
      message( "No study.area defined and no potential.sites given. Using unit interval/square/cube/hyper-cube (as dimension dictates)")
      N <- 100
      potential.sites <- as.matrix( expand.grid( as.data.frame( matrix( rep( 1:N, times=dimension), ncol=dimension)))/N - 1/(2*N)) #the minus is so that the pts are centres of cells
      study.area <- as.matrix( expand.grid( as.data.frame( matrix( c( rep( 0, dimension), rep( 1, dimension)), nrow=2, byrow=TRUE))))
      colnames( potential.sites) <- colnames( study.area) <- paste0("dimension",1:dimension)
    }
    else{
      message( "No study.area defined. Using an interval/rectangle/hyper-rectangle based on the extreme locations in potential.sites")
      #if( !class( potential.sites) %in% c("matrix","data.frame"))
      if( !( "matrix" %in% class( potential.sites) | "data.frame" %in% class( potential.sites)))
        potential.sites <- as.matrix( potential.sites)
      if( is.null( colnames( potential.sites)))
        colnames( potential.sites) <- paste0( "dimension",1:dimension)
      study.area <- as.matrix( expand.grid( as.data.frame( apply( potential.sites, -1, range))))
      colnames( study.area) <- colnames( potential.sites)
    }
    if( dimension==2)
      study.area <- study.area[c(1,3,4,2),] #just so, in this case, it forms a polygon (for checking in boundness)
  }
  #if( !class( study.area) %in% c("matrix","data.frame"))
  if( !( "matrix" %in% class( potential.sites) | "data.frame" %in% class( potential.sites)))
    study.area <- as.matrix( study.area)
  if( is.null( colnames( study.area)))
    colnames( study.area) <- paste0( "dimension",1:dimension)
  if( is.null( potential.sites)){
    #use study.area to sample in
    message( "Sampling from a 100x100 grid (max -- will be cropped) within the study area (well hyper-rectangle if dimension>2)")
    N <- 100
    potential.sites <- as.matrix( expand.grid( as.data.frame( matrix( rep( 1:N, times=dimension), ncol=dimension)))/N - 1/(2*N)) #the minus is so that the pts are centres of cells
    colnames( potential.sites) <- paste0("dimension",1:dimension)
    myRange <- apply( study.area, -1, range)
    for( ii in 1:dimension){
      potential.sites[,ii] <- myRange[1,ii] + (myRange[2,ii]-myRange[1,ii]) * potential.sites[,ii]
    }
    if( dimension==2){
      tmp <- mgcv::in.out( study.area, potential.sites)
      potential.sites <- potential.sites[tmp,]
    }
    if( dimension>2)
      message( "Sampling grid not checked for interior of study area (please do so yourself)")
  }
  N <- nrow( potential.sites)
  if( dimension != ncol( potential.sites))
    stop( "The dimension supplied does not match the dimension that the potential.sites implies")
  if( is.null( inclusion.probs)){
    message( "No inclusion.probs supplied, assuming uniform")
    inclusion.probs <- rep( 1/N, N) #even probs
  }
  #standardise properly
  inclusion.probs <- n * inclusion.probs / sum( inclusion.probs, na.rm=TRUE)
  #standardise inclusion probabilities (for efficient sampling)
  inclusion.probs1 <- inclusion.probs / max( inclusion.probs, na.rm=TRUE)
  
  ret <- list( dimension=dimension, study.area=study.area, potential.sites=potential.sites, inclusion.probs=inclusion.probs, inclusion.probs1=inclusion.probs1, N=N)
  return( ret)
}

"quasiSamp_fromhyperRect" <- function( nSampsToConsider, randStartType=3, designParams) {
  #Generate quasi random numbers within survey area
  
  #initialise the sequence and subsample from it
  samp <- randtoolbox::halton( nSampsToConsider*2, dim=designParams$dimension+1, init=TRUE) #The big sequence of quasi random numbers
  if( randStartType==1)
    skips <- rep( sample( 1:nSampsToConsider, size=1, replace=TRUE), designParams$dimension+1)
  if( randStartType %in% 2:3)
    skips <- sample( 1:nSampsToConsider, size=designParams$dimension+1, replace=TRUE) #the start points
  samp <- do.call( "cbind", lapply( 1:(designParams$dimension+1), function(x) samp[skips[x]+0:(nSampsToConsider-1),x]))  #a tedious way to paste it all together?  
  
  #convert to scale of study region
  myRange <- apply( designParams$study.area, -1, range)  
  for( ii in 1:designParams$dimension)
    samp[,ii] <- myRange[1,ii] + (myRange[2,ii]-myRange[1,ii]) * samp[,ii]
  #select those points within study area
  if( designParams$dimension==2){
    tmp <- mgcv::in.out( designParams$study.area, samp[,1:designParams$dimension])
    samp <- samp[tmp,] #get rid of samples that are outside the study region
  }
  return( samp)
}

quasiSamp <- function (n, dimension = 2, study.area = NULL, potential.sites = NULL, inclusion.probs = NULL, randStartType = 3, nSampsToConsider = 25*n, nStartsToConsider=100*n) 
{
  if (inherits(inclusion.probs, c( "RasterLayer", "SpatRaster"))) {
    samp <- quasiSamp.raster(n = n, inclusion.probs = inclusion.probs, randStartType = randStartType, nSampsToConsider = nSampsToConsider)
    return(samp)
  }
  designParams <- setDesignParams(dimension, study.area, potential.sites, inclusion.probs, n=n)
  c2cAll <- lapply(as.data.frame(designParams$potential.sites), function(xx) sort(unique(xx)))
  c2c <- sapply(c2cAll, function(xx) min(diff(xx)))

  kount <- 1
  repeat{
    samp <- quasiSamp_fromhyperRect(nSampsToConsider, randStartType = randStartType, designParams) #randStartType was hardwired to be 2 (19 May 2023)
    sampIDs <- class::knn1(designParams$potential.sites, samp[, 1:designParams$dimension, drop = FALSE], 1:nrow(designParams$potential.sites))
    closePotential <- designParams$potential.sites[sampIDs, ]
    dist.1d <- abs(closePotential - samp[, 1:designParams$dimension, drop = FALSE])
    within.cells <- matrix(NA, nrow = nrow(dist.1d), ncol = ncol(dist.1d))
    for (ii in 1:designParams$dimension) 
      within.cells[, ii] <- dist.1d[,ii] <= c2c[ii]
    inArea <- rowSums(within.cells) == designParams$dimension
    samp <- samp[inArea, ]
    sampIDs <- sampIDs[inArea]
    sampIDs.2 <- which(samp[, designParams$dimension + 1] < designParams$inclusion.probs1[sampIDs])
    if( (length( sampIDs.2) > 0) & (randStartType!=3 | sampIDs.2[1]==1))
      break
    kount <- kount + 1
    if( kount > nStartsToConsider)
      stop( "Failed to find a design. Too many starts when using the randomisation in Robertson et al (2017). Please consider 1) increasing nStartsToConsider (computationally more demanding), 2) try again (you may get lucky), 3) make inclusion probabilities more even (more likely but possibly undesireable).")
  }
  if (length(sampIDs.2) >= n) 
    sampIDs <- sampIDs[sampIDs.2][1:n]
  else 
    stop("Failed to find a design. Too few samples considered for BAS. It is possible that the inclusion probabilities are very low and uneven OR that the sampling area is very irregular (e.g. long and skinny) OR something else. Please try again (less likely to work) OR make inclusion probabilities more even (more likely but possibly undesireable) OR increase the number of sites considered (likely but computationally expensive).")  
#    stop("Failed to find a design. It is possible that the inclusion probabilities are very low and uneven OR that the sampling area is very irregular (e.g. long and skinny) OR something else. Please try again (less likely to work) OR make inclusion probabilities more even (more likely but possibly undesireable) OR increase the number of sites considered (likely but computationally expensive).")
  samp <- as.data.frame(cbind(designParams$potential.sites[sampIDs,, drop = FALSE], designParams$inclusion.probs[sampIDs],sampIDs))
  colnames(samp) <- c(colnames(designParams$potential.sites), "inclusion.probabilities", "ID")
  return(samp)
}

quasiSamp.raster <- function (n, inclusion.probs, randStartType = 3, nSampsToConsider = 25*n, nStartsToConsider=100*n) 
{
#  if (!inherits(inclusion.probs, "RasterLayer")) 
#    inclusion.probs <- try( terra::rast( inclusion.probs), silent=TRUE)
#  if (!inherits(inclusion.probs, "RasterLayer")) 
#    stop("RasterLayer must be passed as input (argument inclusion.probs). It must be a RasterLayer, or something that can be coerced to a RasterLayer using raster::raster. Please revise, or use quasiSamp (not quasiSamp.raster)")
  terra::values(inclusion.probs) <- terra::values(inclusion.probs)/max(terra::values(inclusion.probs), na.rm = TRUE)
  #saving and scaling the defined IPs
  IP.orig <- inclusion.probs
  terra::values(IP.orig) <- terra::values(IP.orig)/sum(terra::values(IP.orig), na.rm = TRUE)
  tmp1 <- terra::ext(inclusion.probs)
  tmp <- matrix( tmp1[], ncol=2, byrow=TRUE)  #matrix(c(tmp1@xmin, tmp1@xmax, tmp1@ymin, tmp1@ymax), ncol = 2, byrow = TRUE)
  designParams.short <- list(dimension = 2, study.area = matrix(c(tmp[1,1], tmp[2, 1], tmp[1, 2], tmp[2, 1], tmp[1, 2], tmp[2,2], tmp[1, 1], tmp[2, 2]), ncol = 2, byrow = TRUE))
  
  kount <- 1
  repeat{
    samp <- quasiSamp_fromhyperRect(nSampsToConsider, randStartType = randStartType, designParams = designParams.short)
    sampIDs <- terra::extract( x=inclusion.probs, y=samp[, 1:2], cells = TRUE)
    lotsOvals <- sampIDs[, 2]
    sampIDs <- sampIDs[, 1]
    sampIDs.2 <- which(samp[, 3] < lotsOvals)
    if( (length( sampIDs.2) > 0) & (randStartType!=3 | sampIDs.2[1]==1))
      break
    kount <- kount + 1
    if( kount > nStartsToConsider)
      stop( "Failed to find a design. Too many starts when using the randomisation in Robertson et al (2017). Please consider 1) increasing nStartsToConsider (computationally more demanding), 2) try again (you may get lucky), 3) make inclusion probabilities more even (more likely but possibly undesireable).")
  }
  if (length(sampIDs.2) >= n) {
    sampIDs <- sampIDs[sampIDs.2][1:n]
    lotsOvals <- lotsOvals[sampIDs.2][1:n]
  }
  else 
    stop("Failed to find a design. Too few samples considered for BAS. It is possible that the inclusion probabilities are very low and uneven OR that the sampling area is very irregular (e.g. long and skinny) OR something else. Please try again (less likely to work) OR make inclusion probabilities more even (more likely but possibly undesireable) OR increase the number of sites considered (likely but computationally expensive).")
  samp <- as.data.frame(cbind(terra::crds(inclusion.probs, na.rm=FALSE)[sampIDs,], terra::extract( x=IP.orig, y=terra::crds( inclusion.probs, na.rm=FALSE)[sampIDs,]), sampIDs))
  colnames(samp) <- c(colnames(samp)[1:2], "inclusion.probabilities", "ID")
  return(samp)
}
