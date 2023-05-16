
"transectSamp" <- function( n, study.area=NULL, potential.sites=NULL,
                            inclusion.probs=NULL, control=NULL, constrainedSet=NULL){
  #This is the public face of transectSamp.  The other version *.internal just has slightly
  #different argument list, for doing simulation and the like (to avoid recalculating
  #things that do not need to be)
  tmp <- transectSamp.internal( n=n, study.area=study.area, potential.sites=potential.sites,
                         inclusion.probs=inclusion.probs, control=control,
                         constrainedSet=constrainedSet)
  return( tmp)
}

"transectSamp.internal" <- function( n, study.area=NULL, potential.sites=NULL, 
                            inclusion.probs=NULL, control=NULL, constrainedSet, ...) {
  #set the control parameters
  control <- set.transect.control( control)
  #computationally expensive bits from previous runs
  id.list <- getIDs( ...)
  #set up a cluster for parallel computing (only if needed)
  if( is.null( id.list$IDs) | is.null( id.list$adjustedSpecified)){
    cl <- parallel::makeCluster(control$mc.cores)
    parallel::clusterExport(cl, c( "distFun4", "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir = environment())
    parallel::clusterExport( cl, "in.out", envir=.getNamespace("mgcv"))
    cl.flag <- TRUE
  }
  else{
    cl <- NA
    cl.flag <- FALSE
  }
  #set the parameters of the design (survey area, potential sites, inclusion probs, all the different n's, ...)
  designParams <- setDesignParams_transect( study.area=study.area, potential.sites=potential.sites, inclusion.probs=inclusion.probs, control=control)
  #make sure that the potential.sites are appropriately orderd and comprehensive
  check.potential.sites( designParams$potential.sites)
  #set parameters of the individual transects (pattern, direction, point representation, ...)
  transectParams <- setTransectParams2( designParams=designParams, control=control)
  #cell/transect inclusion probs and transect composition (which cells in which transect)
  locProbs.raw <- getProbs2( designParams, transectParams, constraints=constrainedSet, mc.cores=control$mc.cores, prev.ids=id.list$IDs, cl=cl)
  #Edge adjusted inclusion probabilities for transects (nSites X nRotates) and cells (nSites X 1)
  locProbs.edge <- adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=constrainedSet, control, cl=cl, adjustedSpecified=id.list$adjustedSpecified)
  locProbs.edge$transects[is.na( locProbs.edge$transects)] <- 0 #so that non-existant transects are not sampled
  
  transIDs <- list()
  locProbs.scaled.adjusted <- locProbs.edge$transects   #scaled
  
  #get sample sites if control$spat.random.type is quasi or pseudo
  transIDs <- NULL
  if( control$spat.random.type %in% c("quasi","pseudo"))
    transIDs <- getTransects_2step( n, locProbs.scaled.adjusted, designParams, transectParams, control)
  else
    stop( "Unknown randomisation type. It should be either 'quasi' or 'pseudo'")
#  if( control$spat.random.type == "sequential")
#    transIDs <- getTransects_sequential( n=n, locProbs.edge=locProbs.edge, constrainedSet=constrainedSet, designParams=designParams,
#                                         transectParams=transectParams, control=control, cl=cl, id.list=list( IDs=locProbs.raw$IDs, adjustedSpecified=NULL))
  if( is.null( transIDs))
    stop("Unknown randomisation type -- check control list")
  
  ret <- list()
  ret$transect <- cbind( 1:n, do.call( "rbind", lapply( transIDs, function(x) c(x$start_locat, x$direction, x$inclusion.prob.transect*n))))
  colnames( ret$transect) <- c("transect", names( transIDs[[1]]$start_locat), "direction", "transect_inclusionProb")
  ret$points <- lapply( transIDs, function(x) x$locats)
  ret$points <- lapply( 1:length( ret$points), function(x) cbind( x, ret$transect[x, names( transIDs[[1]]$start_locat)[1]], ret$transect[x, names( transIDs[[1]]$start_locat)[2]], ret$transect[x, "direction"], ret$points[[x]]))
  ret$points <- do.call( "rbind", ret$points)
  ret$points <- cbind( ret$points, do.call( 'c', lapply( transIDs, function(x) locProbs.raw$cell.probs[x$site, x$rotate, ])))
  ret$points <- cbind( ret$points, locProbs.edge$cells[do.call( 'c', lapply( transIDs, function(x) locProbs.raw$IDs[x$site, x$rotate, ]))])
  colnames( ret$points) <- c( "transect", paste0( "mid_", names( transIDs[[1]]$start_locat)), "direction", names( transIDs[[1]]$start_locat), "specifiedInclProb", "AdjustedInclProb")
  ret$points <- as.data.frame( ret$points)
  ret$transect <- as.data.frame( ret$transect)
  
  if( control$return.index)
    ret$index <- list( IDs=locProbs.raw$IDs, adjustedSpecified=locProbs.edge)
  
  if( cl.flag)
    parallel::stopCluster( cl)
  
  return( ret)
}

"set.transect.control" <- function( control){
  if( is.null( control))
    ret <- list()
  else
    ret <- control
  if( !"transect.pattern" %in% names( ret))#| ret$transect.pattern == "line")
    ret$transect.pattern <- "line"
  if( !all(ret$transect.pattern == "line"))
    ret$transect.nPts <- nrow( control$transect.pattern)
  if( !"transect.nPts" %in% names( ret))
    ret$transect.nPts <- 20
  if( !"line.length" %in% names( ret))
    ret$line.length <- NULL   #this will be filled out in setTransectParams
  if( !"nRotate" %in% names( ret))
    ret$nRotate <- 21
  if( !"mc.cores" %in% names( ret))
    ret$mc.cores <- 1#max( 1, parallel::detectCores() - 1)
  if( !"transAdjust" %in% names( ret))
    ret$transAdjust <- TRUE
  if( !"edge.max.iter" %in% names( ret))
    ret$edge.max.iter <- 25
  if( !"conv.tol.diff" %in% names( ret))
    ret$conv.tol.diff <- 0.0001
  if( !"gamma" %in% names( ret))
    ret$gamma <- 0.8
  if( !"calcObsProbs" %in% names( ret))
    ret$calcObsProbs <- FALSE
  if( !ret$calcObsProbs & ret$transAdjust)
    ret$calcObsProbs <- TRUE
  if( !"return.index" %in% names( ret))
    ret$return.index <- FALSE
  if( !"spat.random.type" %in% names( ret))
    ret$spat.random.type <- "quasi"#other option is "pseudo" or "sequential"
  if( !"nSampsToConsider" %in% names( ret))
    ret$nSampsToConsider <- 5000
  if( !"nCells" %in% names( ret))
    ret$nCells <- 100
  
  

    
  
  if( !"sigma.alpha" %in% names( ret))
    ret$sigma.alpha <- 0.95

#  if( !"stop.crit" %in% names( ret))
#    ret$stop.crit <- "percent.diff"
##  if( !"conv.perc" %in% names( ret))
##    ret$conv.perc <- 0.01
##  if( !"cell.conv.tol" %in% names( ret))
##    ret$cell.conv.tol <- 0.05
##  if( !"prop.conv.tol" %in% names( ret))
##    ret$prop.conv.tol <- 0.9

  return( ret)
}

"setDesignParams_transect" <- function( study.area, potential.sites, inclusion.probs, control){#Nrotate=10, npts){
  dimension <- 2
  #set the number of sites (in each direction)
  if( is.null( potential.sites)){
    N <- list( dimension1=control$nCells, dimension2=control$nCells, rotate=control$nRotate)
    N$Tot.xy <- N$dimension1*N$dimension2
  }
  else{
    N <- list()
    N$dimension1 <- length( unique( potential.sites[,1]))
    N$dimension2 <- length( unique( potential.sites[,2]))
    N$rotate <- control$nRotate
    N$Tot.xy <- NA
  }
  N$transect.nPts <- control$transect.nPts
  if( !is.null( study.area) & !is.null( potential.sites)){
    message( "User supplied study.area and potential.sites")
  }
  potential.sites.list <- list()
  if( is.null( study.area)){
    if( is.null( potential.sites)){
      message( "No study.area defined and no potential.sites given. Using unit square")
      #the minus is so that the pts are centres of cells
      potential.sites.list <- list( dimension1=1:N$dimension1 / N$dimension1 - 1/(2*N$dimension1), dimension2=1:N$dimension2 / N$dimension2 - 1/(2*N$dimension2), direction=seq( from=0, to=360-360/N$rotate, length=N$rotate))
      potential.sites <- as.matrix( expand.grid( list( dimension1=potential.sites.list$dimension1, dimension2=potential.sites.list$dimension2)))
      study.area <- as.matrix( expand.grid( as.data.frame( matrix( c( rep( 0, dimension), rep( 1, dimension)), nrow=2, byrow=TRUE))))
      study.area <- study.area[c(1,2,4,3),]
      colnames( potential.sites) <- colnames( study.area) <- paste0("dimension",1:(dimension))
    }
    else{
      message( "No study.area defined. Using a rectangle based on the extreme locations in potential.sites")
      #if( !class( potential.sites) %in% c("matrix", "array"))
      if( !("matrix" %in% class( potential.sites)))
        potential.sites <- as.matrix( potential.sites)
      if( ncol( potential.sites) != 2)
        stop( "The argument potential.sites contains more than 2 dimensions. For transect designs, this is not supported. Please revise.")
      if( is.null( colnames( potential.sites)))
        colnames( potential.sites) <- paste0( "dimension",1:dimension)
      potential.sites.list <- list( dimension1=sort( unique( potential.sites[,1])), dimension2=sort( unique( potential.sites[,2])), direction=seq( from=0, to=360-360/N$rotate, length=N$rotate))
      names( potential.sites.list) <- c( colnames( potential.sites), "direction")
      study.area <- as.matrix( expand.grid( as.data.frame( apply( potential.sites, -1, range))))
      colnames( study.area) <- colnames( potential.sites)
      study.area <- study.area[c(1,2,4,3),]
      N$Tot.xy <- nrow( potential.sites)
    }
  }
  else{ 
    if( is.null( potential.sites)){
      message( "No potential.sites given. Using grid within supplied study.area")
      potential.sites.list <- list()
      for( ii in 1:2){
        rang <- range( study.area[,ii])
        stepsize <- diff( rang) / N[[ii]]
        potential.sites.list[[ii]] <- seq( from=rang[1]+stepsize, to=rang[2], length=N[[ii]])
        potential.sites.list[[ii]] <- potential.sites.list[[ii]] - stepsize/2
      }
      potential.sites.list[[3]] <- seq( from=0, to=360-360/N$rotate, length=N$rotate)
      potential.sites <- as.matrix( expand.grid( list( dimension1=potential.sites.list[[1]], dimension2=potential.sites.list[[2]])))
      if( !is.null( colnames( study.area)))
        colnames( potential.sites) <- colnames( study.area)
      else 
        colnames( potential.sites) <- paste0( "dimension",1:dimension)
      names( potential.sites.list) <- c( colnames( study.area), "direction")
      tmp <- mgcv::in.out( study.area, potential.sites)
      potential.sites <- potential.sites[tmp,]
      N$Tot.xy <- nrow( potential.sites)
    }
  }
  if( is.na( N$Tot.xy))
    N$Tot.xy <- nrow( potential.sites)
  if( length( potential.sites.list) == 0)
    potential.sites.list <- list( sort( unique( potential.sites[,1])), sort( unique( potential.sites[,2])), seq( from=0, to=360-360/N$rotate, length=N$rotate))
  names( potential.sites.list) <- c( colnames( potential.sites), "direction")
  averageCellDist <- c( mean( diff( potential.sites.list[[1]])), mean( diff( potential.sites.list[[2]])))
  #if( !class( study.area) %in% c("matrix","data.frame"))
  if( !( "matrix" %in% class( study.area)))
    study.area <- as.matrix( study.area)
  if( is.null( colnames( study.area)))
    colnames( study.area) <- paste0( "dimension",1:dimension)
  if( is.null( colnames( potential.sites)))
    colnames( potential.sites) <- paste0( "dimension",1:dimension)
  if( is.null( inclusion.probs)){
    message( "No inclusion.probs supplied, assuming uniform throughout study.area")
    inclusion.probs <- rep( 1/N$Tot.xy, N$Tot.xy) #even probs
  }
  else
    inclusion.probs <- inclusion.probs / sum( inclusion.probs, na.rm=TRUE)
  
  #adjust limits of study.area to be ever-so-slightly larger than grid (for correct mgcv::in.out behaviour)
  study.area[study.area[,1]==max(study.area[,1]),1] <- max(study.area[,1]) + 5*.Machine$double.eps
  study.area[study.area[,2]==max(study.area[,2]),2] <- max(study.area[,2]) + 5*.Machine$double.eps
  study.area[study.area[,1]==min(study.area[,1]),1] <- min(study.area[,1]) - 5*.Machine$double.eps
  study.area[study.area[,2]==min(study.area[,2]),2] <- min(study.area[,2]) - 5*.Machine$double.eps
  
  ret <- list( study.area=study.area, potential.sites=potential.sites, potential.sites.list=potential.sites.list, inclusion.probs=inclusion.probs, inclusion.probs.original=inclusion.probs, N=N, averageCellDist=averageCellDist)
  return( ret)
}

"check.potential.sites" <- function( pots){
  #check to see if regular grid
  n.dim1 <- length( unique( pots[,1]))
  n.dim2 <- length( unique( pots[,2]))
  if( nrow( pots) != n.dim1*n.dim2)
    stop( "potential.sites needs to specify a grid. Hint: mark locations outside the survey area with an NA.")
  #check to see if ordering is OK
  unique.dim1 <- sort( unique( pots[,1]))
  unique.dim2 <- sort( unique( pots[,2]))
  test.pot.sites <- expand.grid( unique.dim1, unique.dim2)
  if( !all( test.pot.sites==pots))
    stop( "potential.sites needs to be ordered so that the first dimension changes more rapidly than the second. As a concrete example try running expand.grid(1:4,1:5)")
}

"setTransectParams2" <- function( designParams, control){
  transParams <- list()
  if( all( control$transect.pattern=='line')){
    if( is.null( control$line.length)){ #!"line.length" %in% names( control)){
      message( "Linear transect. No length given though -- will assume 10% of maximal boundary of study area.")
      maxRange <- apply( designParams$study.area, MARGIN=2, FUN=function(x) diff( range(x)))
      transParams$line.length <- 0.1 * max( maxRange)
    }
    else
      transParams$line.length <- control$line.length
    if( !"transect.nPts" %in% names( control)){
      cellDensities <- designParams$N[c("x","y")] / apply( designParams$study.area, MARGIN=2, FUN=function(x) diff( range(x)))
      transParams$npts <- ceiling( max( cellDensities) * transParams$line.length)
    }
    else
      transParams$npts <- control$transect.nPts
    transParams$shape <- cbind( 0, seq( from=0, to=transParams$line.length, length=transParams$npts)) #line facing directly north
  }
  else{
    message( "User supplied transect shape.")
    transParams$npts <- nrow( control$transect.pattern)
    transParams$shape <- control$transect.pattern
    #total distance of the transect (approx)
    transParams$line.length <- sum( sqrt( rowSums( (transParams$shape[ -( nrow( transParams$shape)),] - transParams$shape[ -1,])^2)))
  }
  transParams$shape <- scale( transParams$shape, center=TRUE, scale=FALSE)
  #determine points along each of the potential transects for each direction (from origin)
  angle.rad <- degrees2radians(designParams$potential.sites.list$direction)
  transParams$potential.transects <- lapply( angle.rad, FUN=function(x) transParams$shape %*% matrix( c(cos(x), -sin(x), sin(x), cos(x)), byrow=TRUE, ncol=2))
  
  return( transParams)
}

"degrees2radians" <- function(deg){
  rad <- 2*pi*deg / 360
  return( rad)
}

"getProbs2" <- function( designParams, transectParams, constraints, mc.cores, prev.ids=NULL, cl=NULL){
  if( is.null( cl)){
    cl.flag <- TRUE
    if( is.null( mc.cores))
      mc.cores <- max( 1, parallel::detectCores()-1) #want to be greedy but not too greedy
    cl <- parallel::makeCluster(mc.cores)
    parallel::clusterExport(cl, c( "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir = environment())
    parallel::clusterExport( cl, "in.out", envir=.getNamespace("mgcv"))
  }
  else
    cl.flag <- FALSE
  if( is.null( prev.ids)){
    tmp <- parallel::parLapply(cl, 1:nrow( designParams$potential.sites), getSingleSiteLocProbs, designParams=designParams, transectParams=transectParams, constraints=constraints)
    tmp.combined <- parallel::parLapply( cl, tmp, function(x) sapply( x, function( y) if( all( is.na( y$IDs))) matrix( NA, ncol=2, nrow=transectParams$npts) else cbind( y$IDs, y$inclProbs)))
    tmp.combined <- simplify2array( tmp.combined)
    tmp.combined <- aperm( tmp.combined, c(3,2,1))
    tmp.IDs <- tmp.combined[,,1:designParams$N$transect.nPts]
    tmp.probs <- tmp.combined[,,1:designParams$N$transect.nPts+designParams$N$transect.nPts]
  }
  else{
    tmp.probs <- array( designParams$inclusion.probs[prev.ids], dim=dim( prev.ids))
    tmp.IDs <- prev.ids
  }
  TransProbs <- apply( tmp.probs, c(1,2), sum, na.rm=FALSE) #The inclusion probability of the transect
  TransProbs <- TransProbs / sum( TransProbs, na.rm=TRUE)

  return( list( cell.probs=tmp.probs, transect.probs=TransProbs, IDs=tmp.IDs))
}

"getSingleSiteLocProbs" <- function( startLoc, transectParams, designParams, constraints){
  transectInclProbs <- lapply( 1:designParams$N$rotate, FUN=getSingleTransLocProbs, 
                               startLoc=designParams$potential.sites[startLoc,], transectParams=transectParams, 
                               designParams=designParams, constraint.transect=constraints[startLoc,])
  
  return( transectInclProbs)
}

"getSingleTransLocProbs" <- function( currentRotation, startLoc, transectParams, designParams, constraint.transect){
  if( !is.null( constraint.transect))
    if( !constraint.transect[currentRotation])
      return( list( inclProbs=NA, IDs=NA))
  transect <- matrix( rep( startLoc, each=transectParams$npts), ncol=2) + transectParams$potential.transects[[currentRotation]]
  inOut <- mgcv::in.out( designParams$study.area, transect)
  if( !all( inOut)) #transect is not wholly contained within study area
    return( list( inclProbs=NA, IDs=NA))
  xs <- apply( transect[,1,drop=FALSE], 1, function(x) which.min( abs( x-designParams$potential.sites.list[[1]])))
  ys <- apply( transect[,2,drop=FALSE], 1, function(x) which.min( abs( x-designParams$potential.sites.list[[2]])))
  nearID <- xs + (ys-1)*designParams$N$dimension1
  ret <- list()
  ret$inclProbs <- designParams$inclusion.probs[nearID]
  ret$IDs <- nearID

  return( ret)
}

"transectsOverCells" <- function( cell, IDs){
  #y is a cell that we want to discover which transects pass over
#  tmpy <- apply( IDs, 1:2, function( x) cell %in% x)
#  return( which( tmpy, arr.ind=TRUE))
  locats <- which( IDs==cell, arr.ind=TRUE)
  locats <- locats[,-3]
  return( locats)
}

"obsCellProb2" <- function( x, locProbs.raw, toc.list){
  if( length(toc.list[[x]])>0){
    tmp <- apply( matrix( toc.list[[x]], ncol=2), 1, function(x) locProbs.raw$transect.probs[x[1],x[2]])
    return( sum( tmp, na.rm=TRUE))
#    return( sum( locProbs.raw$inclusion.probs[toc.list[[x]]]))#which( locProbs.raw1[[2]]==x)])
  }
  else
    return( 0)
}

"adjustEdge4" <- function( locProbs.raw, designParams, transectParams, constrainedSet, control, cl=NULL, adjustedSpecified=NULL){
  if( is.null( adjustedSpecified)){
    if( is.null( cl)){
      cl <- parallel::makeCluster( control$mc.cores)
      #next three lines for getProbs2 (and obsCellProb for adjustEdge3)
      parallel::clusterExport(cl, c( "obsCellProb2", "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir = environment())
      parallel::clusterExport( cl, "in.out", envir=.getNamespace("mgcv"))
      cl.flag <- TRUE
    }
    if( control$calcObsProbs)
      tranOverCell.list <- parallel::parLapply( cl, 1:designParams$N$Tot.xy, transectsOverCells, IDs=locProbs.raw$IDs)
    if (control$transAdjust == FALSE) {
      if( control$calcObsProbs){
        obsProbs <- parallel::parLapply(cl, 1:designParams$N$Tot.xy, obsCellProb2, locProbs.raw=locProbs.raw, toc.list=tranOverCell.list)
        tmp1 <- simplify2array( obsProbs)
        tmp1 <- tmp1 / sum( tmp1, na.rm=TRUE)
      }
      else
        tmp1 <- NA
      kount <- -99
    }
    else {
      stop.crit <- Inf
#      stop.crit.seamount <- 0#-Inf
      kount <- 1
      message( "Adjusting Inclusion Probabilities for Transect Pattern.", appendLF=TRUE)
#      if( control$stop.crit=="Seamount")
#        message( "Iteration : Seamount-Criterion : Agreement-Criterion")
#      else
        message( "Iteration : Agreement-Criterion", appendLF=TRUE)
#      message( "EXPERIMENTAL VERSION: saving observed and updated probabilities", appendLF=TRUE)
      tol.limit <- control$conv.tol.diff
#      if( control$stop.crit=="Seamount")
#        tol.limit <- control$prop.conv.tol
      flag <- TRUE
      while( flag & kount < control$edge.max.iter){
#      while( !(stop.crit.diff < control$conv.tol.diff) & kount < control$edge.max.iter){
#        if( control$stop.crit=="Seamount")
#          message(kount, " : ", round( stop.crit.seamount, 4), " : ", round( stop.crit, 4), appendLF=TRUE)
#        else
          message( kount, " : ", round( stop.crit, 4), appendLF=TRUE)
        #find the observed inclusion probabilities
        obsProbs <- parallel::parLapply(cl, 1:designParams$N$Tot.xy, obsCellProb2, locProbs.raw=locProbs.raw, toc.list=tranOverCell.list)
        tmp1 <- simplify2array( obsProbs)
        tmp1 <- tmp1 / sum( tmp1, na.rm=TRUE)
        #the modified inclusion probs (target + (target-obs))
        delta <- designParams$inclusion.probs.original - tmp1
        designParams$inclusion.probs <- designParams$inclusion.probs + control$gamma * delta
        designParams$inclusion.probs[designParams$inclusion.probs < 0] <- 0
        designParams$inclusion.probs[designParams$inclusion.probs > 1] <- 1
        #temporary saving for diagnostics for seamounts
        updateProbs <- designParams$inclusion.probs
#        save( obsProbs, updateProbs, file=paste0("iter_",kount,".RData"))
        
        #find the new observed inclusion probabilities
        locProbs.raw <- getProbs2(designParams, transectParams, constraints=constrainedSet, mc.cores = control$mc.cores, prev.ids=locProbs.raw$IDs, cl=cl)
        allNonNA <- which( designParams$inclusion.probs.original != 0 & !is.na( designParams$inclusion.probs.original))
        CellGood <- abs( tmp1[allNonNA] - designParams$inclusion.probs.original[allNonNA]) / designParams$inclusion.probs.original[allNonNA]
        stop.crit.old <- stop.crit 
        stop.crit <- mean( CellGood)
        stop.crit.diff <- abs( stop.crit.old-stop.crit)
        flag <- !( stop.crit.diff < tol.limit)
#        if(control$stop.crit=="Seamount"){
#          stop.crit.seamount <- mean( CellGood < control$cell.conv.tol)
#          flag <- !( stop.crit.seamount > tol.limit)
#        }
        kount <- kount + 1
      }
#      if( control$stop.crit=="Seamount")
#        message(kount, " : ", round( stop.crit.seamount, 4), " : ", round( stop.crit, 4), appendLF=TRUE)
#      else
        message(kount, " : ", round( stop.crit, 4), appendLF=TRUE)
    }
    TranProbs <- apply( locProbs.raw$transect.probs, c(1,2), sum, na.rm=FALSE) #The inclusion probability of the transect
    TranProbs <- TranProbs / sum( TranProbs, na.rm=TRUE)
  
    ret <- list()
    ret$transects <- TranProbs
    ret$cells <- designParams$inclusion.probs
    attr( ret, "conv") <- (kount < control$edge.max.iter)
    attr( ret, "observed") <- tmp1
    attr( ret, "delta") <- designParams$inclusion.probs.original - designParams$inclusion.probs
  }
  else{
    ret <- adjustedSpecified
  }
  
  return( ret)
}  

"rescaleProbs2" <- function( locProbs){
  maxxy <- max( locProbs, na.rm=TRUE)
  return( locProbs / maxxy)
}

{#"find.sigma_transect" <- function( n.tot, survey.area, l, transectParams, control){
#  dim <- ncol( survey.area)
#  if( dim != 2)
#    stop( "Too many dimensions for finding area-of-influence. This is a programming error and unless you are the package maintainer, it almost is surely not your problem. Pleas contact the maintainer.")
#  dim <- 2
#  A <- geometry::polyarea( survey.area[,1], survey.area[,2])
#  A.dashed <- A / (n.tot)
#  if( all( control$transect.pattern=="line"))
#    r <- (-l + sqrt( l^2 + pi*A.dashed)) / pi
#  else{
#    hully <- geometry::convhulln( transectParams$shape, options="FA")
#    perim <- hully$area #inappropriate name for 2D polygons
#    area <- hully$vol #inappropriate name (from conhulln) for 2d polygons
#    r <- ( -perim + sqrt( perim^2 - 4*pi*(area-A.dashed))) / ( 2*pi)
#    
#  }
#  alpha <- control$sigma.alpha#0.9
#  fun <- function(sig)
#    stats::pnorm(r, mean=0, sd=sig) - alpha#note that this has a normalising constant which can 
#          #be ignored / cancelled from the distance function.
#  sigma <- stats::uniroot( fun, interval=c(1e-8, 10*r))$root
#  
#  return( sigma)
#}
}

"getTransects_2step" <- function( n, locProbs, designParams, transectParams, control){
  site.probs <- rowSums( locProbs)
  if( control$spat.random.type=="quasi"){
    tmp <- quasiSamp( n=n, dimension=2, study.area=designParams$study.area, potential.sites=designParams$potential.sites, 
                      inclusion.probs=site.probs, randStartType=2, nSampsToConsider=control$nSampsToConsider)
    tmp <- tmp[,ncol( tmp)]
  }
  if( control$spat.random.type=="pseudo"){
    tmp <- sample( x=1:nrow( locProbs), size=n, prob=site.probs)
  }
  rotat.probs <- apply( locProbs[tmp,], 1, function(x) x/sum( x, na.rm=TRUE))
  rotations <- apply( rotat.probs, 2, function(x) sample( x=1:length(x), size=1, prob=x))
  ret <- list()
  for( ii in 1:n){
    ret[[ii]] <- list()
    ret[[ii]]$longID <- NA
    ret[[ii]]$site <- tmp[ii]
    ret[[ii]]$rotate <- rotations[ii]
    ret[[ii]]$inclusion.prob.transect <- locProbs[tmp[ii], rotations[ii]]
    ret[[ii]]$start_locat <- designParams$potential.sites[tmp[ii],]
    ret[[ii]]$direction <- designParams$potential.sites.list$direction[rotations[ii]]
    ret[[ii]]$locats <- rep( ret[[ii]]$start_locat, each=designParams$N$transect.nPts) + transectParams$potential.transects[[ret[[ii]]$rotate]]
  }
  return( ret)  
}

"getTransect" <- function( locProbs, designParams){
  ret <- list()
  tmp <- sample( x=1:prod( dim( locProbs)), size=1, replace=FALSE, prob=locProbs)
  ret$longID <- tmp
  ret$site <- (tmp-1) %% dim( locProbs)[1] + 1  #R indexes from 1 (correct for first row)
  ret$rotate <- (tmp-1) %/% dim( locProbs)[1] + 1  #+1 needed to correct for first nrow().
  ret$inclusion.prob.transect <- locProbs[ret$longID]
  ret$start_locat <- designParams$potential.sites[ret$site,]
  ret$direction <- designParams$potential.sites.list$direction[ret$rotate]
  
  return( ret)
}

{#"getTransects_sequential" <- function( n, locProbs.edge, constrainedSet, designParams, transectParams, control, cl, id.list){
#  #Area of Incluence parameter -- may need revisiting for non-linear transects
#  sigma <- find.sigma_transect(n.tot=n, survey.area=designParams$study.area, l=transectParams$line.length, 
#                               control=control, transectParams=transectParams)
#  transIDs <- list()
#  for( ii in 1:n){
#    #Choose the transect
#    transIDs[[ii]] <- getTransect( locProbs.edge$transects, designParams)
#    transIDs[[ii]]$locats <- rep( transIDs[[ii]]$start_locat, each=transectParams$npts) + 
#                                transectParams$potential.transects[[transIDs[[ii]]$rotate]]
#    #update for the latest chosen transect
#    tmp <- alterInclProbs_transect( legacy.sites=transIDs[[ii]]$locats, 
#                                  potential.sites=designParams$potential.sites, 
#                                  n=n, inclusion.probs=designParams$inclusion.probs, 
#                                  mc.cores=control$mc.cores, sigma=sigma, cl=cl)
#    designParams$inclusion.probs <- designParams$inclusion.probs.original <- tmp / sum( tmp, na.rm=TRUE)
#    ##refind the transect probs
#    #cell/transect inclusion probs and transect composition (which cells in which transect)
#    locProbs.raw <- getProbs2( designParams, transectParams, constraints=constrainedSet, mc.cores=control$mc.cores, prev.ids=id.list$IDs, cl=cl)
#    
#    #Edge adjusted inclusion probabilities for transects (nSites X nRotates) and cells (nSites X 1)
#    locProbs.edge <- adjustEdge4( locProbs.raw, designParams, transectParams, constrainedSet=constrainedSet, control, cl=cl, 
#                                  adjustedSpecified=id.list$adjustedSpecified)
#    locProbs.edge$transects[is.na( locProbs.edge$transects)] <- 0 #so that non-existant transects are not sampled
#  }
#  return( transIDs)
#  
#}
}

"getIDs" <- function( ...){
  tmp <- list( ...)
  if( !"IDs" %in% names( tmp))
    tmp$IDs <- NULL
  if( !"transectsOverCells" %in% names( tmp))
    tmp$transectOverCells <- NULL
  if( !"adjustedSpecified" %in% names( tmp))
    tmp$adjustedSpecified <- NULL
  return( list( IDs=tmp$IDs, adjustedSpecified=tmp$adjustedSpecified))#transectsOverCells=tmp$transectsOverCells))
}

"distFun4" <- function(ii, legacy.sites, potential.sites){
  lpt <- legacy.sites[ii, ]
  disty <- sweep(potential.sites, 2, FUN = "-", STATS = lpt)
  disty <- disty^2
  disty <- rowSums(disty, na.rm=TRUE)
  disty <- sqrt(disty)
  return(disty)
}

{#"alterInclProbs_transect" <- function (legacy.sites, potential.sites, n, inclusion.probs, 
#                              mc.cores, sigma, cl) {
#  #methodologically, this is the same as alterInclProbs().  However, there are different
#  #arguments as this will have to be called repeatedly.
#  N <- nrow(potential.sites)
##  cl <- parallel::makeCluster(mc.cores)
##  parallel::clusterExport(cl, c("potential.sites", "legacy.sites", 
##                                "distFun4"), envir = environment())
#  tmp <- parallel::parLapply(cl, 1:nrow(legacy.sites), distFun4, 
#                             legacy.sites=legacy.sites, potential.sites=potential.sites)
##  parallel::stopCluster(cl)
#  tmp <- do.call("cbind", tmp)
#  if( !is.null( n))
#    inclusion.probs <- n * inclusion.probs / sum( inclusion.probs, na.rm=TRUE)
#  else
#    n <- sum(inclusion.probs, na.rm=TRUE)
#  if (is.null(sigma)) 
#    sigma <- find.sigma(n = n, nL = nrow(legacy.sites), potSites = potential.sites)
#  legacyLocs <- apply(tmp, 2, which.min)
#  d <- exp(-(tmp^2)/(2 * sigma * sigma))
#  adj.probs.ret <- inclusion.probs
#  for (jj in 1:ncol(tmp)) {
#    adj <- (1 - inclusion.probs[legacyLocs[jj]])
#    adj.probs.ret <- adj.probs.ret - adj * adj.probs.ret * 
#      d[, jj]
#  }
#  tmp1 <- n * adj.probs.ret/sum(adj.probs.ret, na.rm=TRUE)
#  return(tmp1)
#}
}

"cellFindDesend" <- function( X, IDs, bathy, in.area, designParams, descend.cutoff=0){
  #Another support for finding descending transects
  
  inner.in.area <- in.area[IDs[X,,]]
  inner.in.area <- matrix( inner.in.area, nrow=designParams$N$rotate,	ncol=designParams$N$transect.nPts)
  tmptmp <- apply( inner.in.area, 1, function(y) all( y))
  
  curBathy <- bathy[IDs[X,,]]
  curBathy <- matrix( curBathy, nrow=designParams$N$rotate, ncol=designParams$N$transect.nPts)
#  tmp <- sign( t( apply( curBathy, 1, diff)))
  tmp <- t( apply( curBathy, 1, diff))
  tmp <- ifelse( tmp < descend.cutoff, -1, 
                 ifelse( tmp > descend.cutoff, 1,
                         ifelse( tmp==descend.cutoff, 0, NA)))
  #supporting function for finding descending transects
#  innerFunny <- function(y){
#    if( all( y %in% c(0,1)) | all( y %in% c(0,-1)))
#      return( "descend")
#    if( all( y %in% c(NA, NaN)))
#      return( "allNA_NaN")
#    y1 <- na.omit( y)
#    if( all( y1 %in% c(0,1)) | all( y1 %in% c(0,-1)))
#      return( "descendAndNA_NaN")
#    return( "upAndDown")  #default case capturing -1, 0, 1 (with/without NAs)
#  }
  innerFunny <- function( y){
    if( all( y %in% c(0,1)) | all( y %in% c(0,-1)))
      return( "descend")
    if( all( y %in% c(0,1,NA,NaN)) | all( y %in% c(0,-1,NA,NaN)))
      return( "descendAndNA_NaN")
    if( all( y %in% c(NA,NaN)))
      return( "allNA_NaN")
    if( all( y %in% c(0,1,-1)))
      return( "upAndDown")
    if( all( y %in% c(0,1,-1,NA,NaN)))
      return( "upAndDownAndNA_NaN")
    return( "Dunno_PROBLEM")  #catch
  }
  tmp <- apply( tmp, 1, innerFunny)
  tmp[!tmptmp] <- "outsideArea"
  
  return( tmp)
}

"findDescendingTrans" <- function( potential.sites, bathy, in.area, descend.cutoff=0, control=NULL){
  
  #flag if there are missing values in bathy -- user will need to decide what to do
  if( any( is.na( bathy))){
    message("There are NA values in the bathymetry variable!")
    message("These will be returned as FALSE to being descending transects")
    message("Please decide if this is appropriate and correct those that are not\n")
  }
  
  #set control, or augment as needed
  control <- set.transect.control( control)
  
  #set up the cluster for parallel
  cl <- parallel::makeCluster(control$mc.cores)
  #will need to change envir when going back into MBHdesign
  parallel::clusterExport(cl, c("distFun4", "getSingleTransLocProbs", "getSingleSiteLocProbs"), envir = environment())#.getNamespace("MBHdesign"))
  parallel::clusterExport(cl, "in.out", envir = .getNamespace("mgcv"))
  
  #set up the design stuff, as you would within transectSamp
  designParams <- setDesignParams_transect(study.area = NULL, potential.sites = potential.sites, 
                                                       inclusion.probs = rep( 1, nrow( potential.sites)), control = control)
  check.potential.sites(designParams$potential.sites)
  transectParams <- setTransectParams2(designParams = designParams, control = control)
  constrainedSet <- matrix( TRUE, nrow=nrow( potential.sites), ncol=designParams$N$rotate)
  locProbs.raw <- getProbs2(designParams, transectParams, constraints=constrainedSet, mc.cores = control$mc.cores, prev.ids = NULL, cl = cl)
  IDs <- locProbs.raw$IDs
  
  #actually flag the descending transects
  #will need to be updated for environment?
  IDs.down <- parallel::parLapply( cl, X=1:nrow(potential.sites), cellFindDesend, IDs=IDs, bathy=bathy, in.area=in.area, descend.cutoff=descend.cutoff, designParams=designParams)
  tmp <- do.call( "rbind", IDs.down)
  
  parallel::stopCluster( cl)
  
  return( tmp)
}

"findTransFromPoint" <- function( potential.sites, originPoints, in.area, control = NULL){
  control <- set.transect.control(control)
  cl <- parallel::makeCluster(control$mc.cores)
  parallel::clusterExport(cl, c("distFun4", "getSingleTransLocProbs", 
                                "getSingleSiteLocProbs"), envir = environment())
  parallel::clusterExport(cl, "in.out", envir = .getNamespace("mgcv"))
  designParams <- setDesignParams_transect(study.area = NULL, 
                                           potential.sites = potential.sites, inclusion.probs = rep(1, 
                                           nrow(potential.sites)), control = control)
  check.potential.sites(designParams$potential.sites)
  transectParams <- setTransectParams2(designParams = designParams, 
                                       control = control)
  constrainedSet <- matrix(TRUE, nrow = nrow(potential.sites), 
                           ncol = designParams$N$rotate)
  locProbs.raw <- getProbs2(designParams, transectParams, constraints = constrainedSet, 
                            mc.cores = control$mc.cores, prev.ids = NULL, cl = cl)
  IDs <- locProbs.raw$IDs
  IDs.startEnd <- IDs[,,c(1,control$transect.nPts)]
  max.IDs <- apply( originPoints, 1, function(xx) which.min( (potential.sites[,1] - xx[1])^2 + (potential.sites[,2]-xx[2])^2))
  
  tmp <- apply( IDs.startEnd, c(1,2), function(y) any( y %in% max.IDs))
  startFromPoint <- array( NA, dim=dim( tmp))
  startFromPoint[tmp==TRUE] <- "startsFromPoint"
  startFromPoint[tmp==FALSE] <- "badStart"
  
  #making sure that transects are within area
  inner.in.area <- matrix( NA, nrow=length( in.area), ncol=designParams$N$rotate)
  for( X in 1:length( in.area)){
    tmp.inner <- in.area[IDs[X,,]]
    tmp.inner <- matrix( tmp.inner, nrow=designParams$N$rotate,	ncol=designParams$N$transect.nPts)
    inner.in.area[X,] <- apply( tmp.inner, 1, all)
  }
  startFromPoint[!inner.in.area] <- "outsideArea" 
  
  parallel::stopCluster(cl)
  return(startFromPoint)
}



