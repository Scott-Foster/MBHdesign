
"countFromRow" <- function( x, wtMat, inclusion.probs){
  inclusion.probs <- terra::unwrap( inclusion.probs)
  tmpRow <- terra::focalValues( x=inclusion.probs, w=wtMat, row=x, nrows=1, fill=NA)
  tmpRow <- tmpRow[,as.numeric( t( wtMat)) == 1]
  
  tmpCount <- apply( tmpRow, 1, FUN=function(xx) sum( !is.na( xx) & xx!=0))
  
  return( tmpCount)
}

"IP.barFromRow" <- function( x, ras, wtMat){
  ras <- terra::unwrap( ras)
  tmpRow <- terra::focalValues( x=ras, w=wtMat, row=x, nrows=1)
  tmpRow <- tmpRow[,as.numeric( t( wtMat)) == 1]
  
  tmpCount <- apply( tmpRow, 1, FUN=function(xx) sum( xx, na.rm=TRUE))
  
  return( tmpCount)
}

"getQcell" <- function( x, cellRowIDs, cellIDsForRows, IP.bar, IP.w, clusterSize){
  local.RowID <- cellIDsForRows[x,]
  local.RowID <- local.RowID[!is.na( local.RowID)]
  local.bars <- IP.bar[local.RowID][,1]
  
  local.condP <- as.numeric( IP.w[cellRowIDs[x]]) / (local.bars)
  tmpStruct1 <- ( 1 - local.condP)^clusterSize
  res <- rep( NA, 2)
  #observed prob
  res[1] <- sum( local.bars * (1-tmpStruct1), na.rm=TRUE)
  res[2] <- sum( 1 + (clusterSize-1) * tmpStruct1, na.rm=TRUE)
  
  return( res)
}

"getQrow" <- function(x, IP.bar, IP.w, clusterSize, wtMat, myCellIDs){
  IP.bar <- terra::unwrap( IP.bar)
  IP.w <- terra::unwrap( IP.w)
  myCellIDs <- terra::unwrap( myCellIDs)
  cellIDsForRows <- terra::focalValues( x=myCellIDs, w=rev( dim( wtMat)), row=x, nrows=1)
  cellIDsForRows <- cellIDsForRows[, as.numeric( t( wtMat))==1] #the cells that make up the focal areas for that row ncol x focus_size
  cellRowIDs <- terra::cellFromRowCol( object=myCellIDs, row=x, col=1:terra::ncol( myCellIDs)) #the cells that make up the row (ncol x 1)
  
  tmp <- sapply( 1:length( cellRowIDs), getQcell, cellRowIDs=cellRowIDs, cellIDsForRows=cellIDsForRows, IP.bar=IP.bar, IP.w=IP.w, clusterSize=clusterSize)
  
  return( tmp)
}

"scaleQ" <- function( Q, tot, maskky){
  scaleFac <- tot / terra::global( Q[[1]], sum, na.rm=TRUE)[1,1]
#  scaleFac <- tot / sum( Q[,1], na.rm=TRUE)
  Q <- Q / scaleFac
#  Q[,1] <- Q[,1] / scaleFac
#  Q[,2] <- Q[,2] / scaleFac
  Q[is.na( Q)] <- 0
#  Q[is.na( Q[,1]),] <- 0
  Q <- terra::mask( Q, maskky)
#  Q[maskky,] <- NA
  return( Q)
}

"getCrits" <- function( Q, IP, IP.w, doPlot, maskky, kount, clusterRadius){
  crit.ras <- ( Q - IP)
#  crit.ras <- ( Q - terra::values( IP, na.rm=FALSE))
  crit.ras <- mask( crit.ras, maskky)
#  crit.ras[maskky] <- NA
  crit.ras[crit.ras==0] <- NA #remove those cells that are zero.  They don't tell us much about convergence.
  crit <- c( terra::global( abs( crit.ras), quantile, probs=0.999, na.rm=TRUE)[1,1], terra::global( abs( crit.ras), median, na.rm=TRUE)[1,1], terra::global( abs( crit.ras), max, na.rm=TRUE)[1,1])
#  crit <- c( stats::quantile( abs( crit.ras), probs=0.999, na.rm=TRUE), stats::median( abs( crit.ras), na.rm=TRUE), max( abs( crit.ras), na.rm=TRUE))
  
  if( doPlot){
    tmpRas <- crit.ras#terra::rast( cbind( terra::crds( IP, na.rm=FALSE), crit.ras), type='xyz')
    graphics::par( mfrow=c(1,3))
    terra::plot( tmpRas, zlim=crit[1]*c(-1,1), main="Observed - Specified")
    worst <- terra::crds( tmpRas)[which( abs( terra::values( tmpRas, na.rm=FALSE))>=crit[3]),,drop=FALSE]  #extra subsetting just in case there is more than one point
#    worst <- terra::crds( tmpRas)[which( abs( tmpRas))>=crit[3],,drop=FALSE]  #extra subsetting just in case there is more than one point
    graphics::points( x=worst[,1], y=worst[,2], col=grDevices::grey(0.5), pch=20, cex=0.1)
    graphics::points( x=worst[1,1], y=worst[1,2], col='black', pch=20, cex=0.25)
    terra::plot( IP.w, xlim=worst[1,1]+clusterRadius*c(-1,1), ylim=worst[1,2]+clusterRadius*c(-1,1), main="Working")
    graphics::points( worst[,1], worst[,2], pch=20, col='black')
    terra::plot( tmpRas, xlim=worst[1,1]+clusterRadius*c(-1,1), ylim=worst[1,2]+clusterRadius*c(-1,1), main="Observed - Specified")
    graphics::points( worst[,1], worst[,2], pch=20, col='black')
    graphics::mtext( paste0("Iter: ",kount), side=1, line=-2, outer=TRUE)
  }
  return( list( ras=crit.ras, max=crit[1], median=crit[2]))
}

"alterInclProbs.cluster" <- function( nCluster, clusterSize, clusterRadius, inclusion.probs, 
                                    maxIter=50, tolerance=NULL, mc.cores=parallel::detectCores()-1, doPlot=FALSE){
  #no precalculation (too memory hungry) but parallelisation (per row)
  
  #this function is riddled with the use of terra::values. Rob Hijmans suggests that this is a potential error-producing problem for large rasters
  
  #rescale to clusterSize
  inclusion.probs <- nCluster * clusterSize * inclusion.probs / (terra::global( inclusion.probs, sum, na.rm=TRUE)[1,1])
    
  if( is.null( tolerance)){
    tmpRas <- terra::ifel( inclusion.probs>0, inclusion.probs, NA)
    tmp <- terra::global( tmpRas, median, na.rm=TRUE)[1,1]
    tolerance <- 0.0001 * tmp
  }
  
  #message about task.
  message( "Optimisation summary.")
  message( "Maximum Iterations: ",maxIter,". Absolute tolerance: ",signif( tolerance,4),". Number of cores: ",mc.cores,".\n")
  
  #vectorised version of the spatial mask
#  myMask <- !is.na( inclusion.probs)
#  myMask <- is.na( terra::values( inclusion.probs))
#  myMaskID <- which( myMask)
  #raster of cell numbers for extraction
  myCellIDs <- terra::rast( inclusion.probs, nlyrs=1, names="myCellIDs", vals=1:terra::ncell( inclusion.probs))
  myCellIDs <- mask( x=myCellIDs, mask=inclusion.probs)
  myCellIDs.wrap <- terra::wrap( myCellIDs)
  
  #weight matrix for focal calculations
  wtMat <- terra::focalMat( x=inclusion.probs, d=clusterRadius, type="circle")
  if( !max( dim( wtMat), na.rm=TRUE) > 1)
    stop( "clusterRadius is smaller than a cell. Please re-consider.")
  wtMat[wtMat!=0] <- 1
  wtMatMask <- which( t( wtMat)==1)
  wtDim <- dim( wtMat)

  #setup cluster
  cl <- parallel::makeCluster(mc.cores)
  parallel::clusterExport( cl, c( "countFromRow", "IP.barFromRow", "getQrow", "getQcell"), envir = .getNamespace( "MBHdesign"))
  parallel::clusterExport( cl, c( "focalValues", "cellFromRowCol", "unwrap", "ncol"), envir=.getNamespace("terra"))
  inclusion.probs.wrap <- terra::wrap( inclusion.probs)
  
  #count cells in neighbourhoods
#  countty <- terra::focalValues( x=inclusion.probs, w=wtMat, row=1, nrows=terra::nrow( inclusion.probs), fill=NA)
#  countty <- apply( countty, 1, function(xx) sum( !is.na( xx)))
  countty <- parallel::parSapply( cl, 1:terra::nrow( inclusion.probs), countFromRow, wtMat=wtMat, inclusion.probs=inclusion.probs.wrap)
  countty <- terra::rast( t(countty), crs=terra::crs(inclusion.probs), ext=terra::ext( inclusion.probs))
  countty <- terra::mask( countty, inclusion.probs)
#  countty[myMaskID] <- NA
  
  #Initialise the working values from linear approx
  IP.w <- inclusion.probs / (clusterSize*countty)
#  IP.w <- inclusion.probs / terra::rast( cbind( terra::crds(inclusion.probs, na.rm=FALSE), clusterSize*terra::values( countty)), type='xyz', crs=terra::crs( inclusion.probs))
  IP.w.wrap <- terra::wrap( IP.w)
  #IP.w <- inclusion.probs / ( clusterSize * countty)#sum( wtMat==1))
  IP.w <- ifel( is.na( IP.w), 0, IP.w)
#  terra::values( IP.w)[is.na( terra::values( IP.w))] <- 0
  IP.w <- ifel( inclusion.probs==0, 0, IP.w)
#  terra::values( IP.w)[terra::values( inclusion.probs)==0] <- 0
  IP.w <- mask( x=IP.w, mask=inclusion.probs)
#  terra::values( IP.w)[myMaskID] <- NA
  IP.w.old <- IP.w
  #the first (pre-)iteration
  kount <- 0
  #take local averages
#  IP.bar <- terra::focalValues( x=IP.w, w=wtMat, row=1, nrows=terra::nrow( IP.bar), fill=NA)
#  wtMat.numeric <- as.numeric( wtMat)
#  IP.bar <- apply( IP.bar, 1, function(xx) sum( xx * wtMat.numeric, na.rm=TRUE))
#  IP.bar <- terra::rast( cbind( terra::crds( inclusion.probs), IP.bar), type='xyz')
  IP.bar <- terra::rast( t( parallel::parSapply( cl, 1:terra::nrow( IP.w), IP.barFromRow, ras=IP.w.wrap, wtMat=wtMat)), crs=terra::crs( IP.w), ext=terra::ext( IP.w))
  IP.bar <- mask( x=IP.bar, mask=inclusion.probs)
#  terra::values( IP.bar)[myMaskID] <- NA
  IP.bar.wrap <- terra::wrap( IP.bar)
  
  QandDeriv <- parallel::parLapply( cl, 1:nrow( inclusion.probs), getQrow, IP.bar=IP.bar.wrap, IP.w=IP.w.wrap, clusterSize=clusterSize, wtMat=wtMat, myCellIDs=myCellIDs.wrap)
  QandDeriv <- lapply( QandDeriv, function(xx) t(xx))
  QandDeriv <- do.call( "rbind", QandDeriv)
  QandDeriv <- terra::rast( aperm( array( QandDeriv, dim=c( dim( inclusion.probs)[2:1],2)), c(2,1,3)), crs=terra::crs( inclusion.probs), ext=terra::ext( inclusion.probs))

  #avoid numerical errors, which when summed up can be a gotcha.
  #More importantly gives context to the working probs.
  QandDeriv <- scaleQ( Q=QandDeriv, tot=terra::global( inclusion.probs, sum, na.rm=TRUE)[1,1], maskky=inclusion.probs)
#  QandDeriv <- scaleQ( Q=QandDeriv, tot=terra::global( inclusion.probs, sum, na.rm=TRUE)[1,1], maskky=myMaskID)

  #create and possibly plot crit values
  crit <- getCrits( QandDeriv[[1]], inclusion.probs, IP.w, doPlot, inclusion.probs, kount, clusterRadius)
  crit.old <- crit
  crit.old$ras[] <- Inf
  crit.old$ras <- terra::mask( crit.old$ras, inclusion.probs)
#  crit.old$ras <- rep( Inf, length( crit.old$ras))
  message( "Iteration: ", kount, ". Step Length: ", NA, ". Tolerance: ", signif( crit$max, 4), " (99.9 perc.) ", signif( crit$median, 4), " (median).")
  kount <- kount+1
  
  numNonZero <- terra::global( !is.na( crit.old$ras) & (crit.old$ras !=0), sum, na.rm=TRUE)[1,1]
#  numNonZero <- sum( !is.na( crit.old$ras) & (crit.old$ras !=0), na.rm=TRUE)
  
  #the actual iteration now  
  while( kount <= maxIter & crit$max>tolerance){
    
    #adjust the previous update if necessary
    setWorse <- !is.na( crit.old$ras) & (crit.old$ras !=0) & (abs( crit.old$ras) < abs( crit$ras))
#    setWorse <- which( !is.na( crit.old$ras) & (crit.old$ras !=0) & abs( crit.old$ras) < abs( crit$ras))
    numworse <- terra::global( setWorse, sum, na.rm=TRUE)[1,1]
    if( terra::global( setWorse, any, na.rm=TRUE)[1,1] > 0){#length( setWorse) > 0){
      message( "Previous step increases tolerance of ", numworse, " (",signif( 100*numworse/numNonZero, 3), "% of all non-zero) cells. Reversing the update of those cells.")
      #again, not updating, which may cause issues in future
      IP.w[setWorse] <- IP.w.old[setWorse]
#      terra::values( IP.w)[setWorse] <- terra::values( IP.w.old)[setWorse]
    }
    
    #step size (increasing as iterations progress)
    alpha <- 0.9 * ( 1 - exp( - kount / 5))
    
    #Pieces of the pie for updating #convert rasters to numeric
    Q <- ( QandDeriv[[1]] - terra::values( inclusion.probs))^2
#    Q <- ( QandDeriv[,1] - terra::values( inclusion.probs))^2
    Q.dash <- 2 * ( QandDeriv[[1]] - terra::values( inclusion.probs)) * QandDeriv[[2]]
#    Q.dash <- 2 * ( QandDeriv[,1] - terra::values( inclusion.probs)) * QandDeriv[,2]
    adj <- (Q/Q.dash)*alpha
  
    #catching those updates that are obviously stupid.
    tmpzero <- ifel( IP.w.old < adj, 1, NA)
    numzero <- terra::global( tmpzero, sum, na.rm=TRUE)[1,1]
    if( !is.na( numzero) & numzero > 0){
      message( terra::global( tmpzero, sum, na.rm=TRUE), " Updates less than 0. Reducing those probabilities towards zero (by 0.5*step_to_zero).")
      adj[tmpzero==1] <- 0.5* IP.w.old[tmpzero==1]
    }
##    tmpzero <- which( terra::values( IP.w.old) < adj)
#    if( length( tmpzero) > 0){
#      message( length( tmpzero), " Updates less than 0. Reducing those probabilities towards zero (by 0.5*step_to_zero).")
#      adj[tmpzero] <- 0.5*terra::values( IP.w.old)[tmpzero]
#    }
    tmpone <- ifel( adj > 1-IP.w.old, 1, NA)
    numone <- terra::global( tmpone, sum, na.rm=TRUE)[1,1]
    if( !is.na( numone) & numone > 0){
      message( terra::global( tmpone, sum, na.rm=TRUE), " Updates larger than 1. Increasing those probabilities towards 1 by 0.5*step_to_one.")
      adj[tmpone==1] <- 0.5* (1-IP.w.old[tmpone==1])
    }
#    tmpone <- which( adj > 1-terra::values( IP.w.old))
#    if( length( tmpone) > 0){
#      message( length( tmpone), " Updates larger than 1. Increasing those probabilities towards 1 by 0.5*step_to_one.")
#      adj[tmpone] <- 0.5 * (1-terra::values( IP.w.old)[tmpone])
#    }

    #actually do the update
    IP.w.old <- IP.w
    IP.w <- IP.w.old - adj
#    terra::values( IP.w) <- terra::values( IP.w.old) - adj
    IP.w.wrap <- terra::wrap( IP.w)
    
    #take local averages
    IP.bar <- terra::rast( t( parallel::parSapply( cl, 1:terra::nrow( IP.w), IP.barFromRow, ras=IP.w.wrap, wtMat=wtMat)), crs=terra::crs( IP.w), ext=terra::ext( IP.w))
    IP.bar <- terra::mask( IP.bar, inclusion.probs)
#    terra::values( IP.bar)[myMaskID] <- NA #mixed casting...
    IP.bar.wrap <- terra::wrap( IP.bar)
    
    QandDeriv <- parallel::parLapply( cl, 1:nrow( inclusion.probs), getQrow, IP.bar=IP.bar.wrap, IP.w=IP.w.wrap, clusterSize=clusterSize, wtMat=wtMat, myCellIDs=myCellIDs.wrap)
    QandDeriv <- lapply( QandDeriv, function(xx) t(xx))
    QandDeriv <- do.call( "rbind", QandDeriv)
    QandDeriv <- terra::rast( aperm( array( QandDeriv, dim=c( dim( inclusion.probs)[2:1],2)), c(2,1,3)), crs=terra::crs( inclusion.probs), ext=terra::ext( inclusion.probs))
    
    #avoid numerical errors, which when summed up can be a gotcha.
    #More importantly gives context to the working probs.
#    QandDeriv <- scaleQ( Q=QandDeriv, tot=sum( terra::values( inclusion.probs), na.rm=TRUE), maskky=myMaskID)
    QandDeriv <- scaleQ( Q=QandDeriv, tot=terra::global( inclusion.probs, sum, na.rm=TRUE)[1,1], maskky=inclusion.probs)
    
    #create and possibly plot crit values
    crit.old <- crit  #update irrespective of improvement or not
    crit <- getCrits( QandDeriv[[1]], inclusion.probs, IP.w, doPlot, inclusion.probs, kount, clusterRadius)
    
    #output to console (or whereever)
    message( "Iteration: ", kount, ". Step Length: ", signif( alpha, 2), ". Tolerance: ", signif( crit$max, 4), " (99.9 perc.) ", signif( crit$median, 4), " (median).")
    kount <- kount+1    
  }
  if( kount == maxIter)
    message( "Failed to Converge.")
  else
    message( "Converged!")
  
  #form the return object 'n' stuff
  IP.w[is.na( IP.w)] <- 0
  IP.w <- terra::mask( IP.w, inclusion.probs)
#  IP.w[myMask] <- NA
  IP.w.wrap <- terra::wrap( IP.w)
  QandDeriv[[1]][is.na( QandDeriv[[1]])] <- 0
#  QandDeriv[is.na( QandDeriv[,1]),1] <- 0
  QandDeriv <- terra::mask( QandDeriv, inclusion.probs)
#  QandDeriv[myMask,1] <- NA
  IP.o <- QandDeriv[[1]]
#  IP.o <- terra::rast( cbind( terra::crds( inclusion.probs, na.rm=FALSE), QandDeriv[,1]), type='xyz', crs=terra::crs(inclusion.probs))
  IP.bar <- terra::rast( t( parallel::parSapply( cl, 1:terra::nrow( IP.w), IP.barFromRow, ras=IP.w.wrap, wtMat=wtMat)), crs=terra::crs( IP.w), ext=terra::ext( IP.w))
  IP.bar <- terra::mask( IP.bar, inclusion.probs)
#  terra::values( IP.bar)[myMaskID] <- NA
  
  parallel::stopCluster(cl)
  
  res <- c( inclusion.probs, IP.o, IP.w, IP.bar)
  names( res) <- c("IP.s","IP.o","IP.w","IP.bar")
  
  return( res)
}
