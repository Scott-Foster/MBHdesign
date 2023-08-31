"countFromRow" <- function( x, wtMat, inclusion.probs){
  tmpRow <- raster::getValuesFocal( inclusion.probs, row=x, nrows=1, ngb=rev( dim( wtMat)))
  tmpRow <- tmpRow[,as.numeric( t( wtMat)) == 1]
  
  tmpCount <- apply( tmpRow, 1, FUN=function(xx) sum( !is.na( xx) & xx!=0))
  
  return( tmpCount)
}

"IP.barFromRow" <- function( x, ras, wtMat){
  tmpRow <- raster::getValuesFocal( ras, row=x, nrows=1, ngb=rev( dim( wtMat)))
  tmpRow <- tmpRow[,as.numeric( t( wtMat)) == 1]
  
  tmpCount <- apply( tmpRow, 1, FUN=function(xx) sum( xx, na.rm=TRUE))
  
  return( tmpCount)
}

"getQcell" <- function( x, cellRowIDs, cellIDsForRows, IP.bar, IP.w, clusterSize){
  local.RowID <- cellIDsForRows[x,]
  local.RowID <- local.RowID[!is.na( local.RowID)]
  local.bars <- IP.bar[local.RowID]
  
  local.condP <- IP.w[cellRowIDs[x]] / (local.bars)
  tmpStruct1 <- ( 1 - local.condP)^clusterSize
  res <- rep( NA, 2)
  #observed prob
  #  res[1] <- sum( local.counts * local.bars * local.probSampled, na.rm=TRUE)
  res[1] <- sum( local.bars * (1-tmpStruct1), na.rm=TRUE)
  #deriv of observed prob (local.bars not considered constant.)
  res[2] <- sum( 1 + (clusterSize-1) * tmpStruct1, na.rm=TRUE)
  
  return( res)
}

"getQrow" <- function(x, IP.bar, IP.w, clusterSize, wtMat, myCellIDs){
  cellIDsForRows <- raster::getValuesFocal( myCellIDs, row=x, nrows=1, ngb=rev( dim( wtMat)))
  cellIDsForRows <- cellIDsForRows[, as.numeric( t( wtMat))==1] #the cells that make up the focal areas for that row ncol x focus_size
  cellRowIDs <- raster::cellFromRow( object=myCellIDs, rownr=x) #the cells that make up the row (ncol x 1)
  
  tmp <- sapply( 1:length( cellRowIDs), getQcell, cellRowIDs=cellRowIDs, cellIDsForRows=cellIDsForRows, IP.bar=IP.bar, IP.w=IP.w, clusterSize=clusterSize)
  
  return( tmp)
}

"scaleQ" <- function( Q, tot, maskky){
  scaleFac <- tot / sum( Q[,1], na.rm=TRUE)
  Q[,1] <- Q[,1] / scaleFac
  Q[,2] <- Q[,2] / scaleFac
  Q[is.na( Q[,1]),] <- 0
  Q[maskky,] <- NA
  return( Q)
}

"getCrits" <- function( Q, IP, IP.w, doPlot, maskky, kount, clusterRadius){
  crit.ras <- ( Q - raster::values( IP))
  crit.ras[maskky] <- NA
  crit.ras[crit.ras==0] <- NA #remove those cells that are zero.  They don't tell us much about convergence.
  crit <- c( stats::quantile( abs( crit.ras), probs=0.999, na.rm=TRUE), stats::median( abs( crit.ras), na.rm=TRUE), max( abs( crit.ras), na.rm=TRUE))
  
  if( doPlot){
    tmpRas <- raster::rasterFromXYZ( cbind( coordinates( IP), crit.ras))
    graphics::par( mfrow=c(1,3))
    raster::plot( tmpRas, zlim=crit[1]*c(-1,1), main="Observed - Specified")
    worst <- raster::coordinates( tmpRas)[which( abs( raster::values( tmpRas))>=crit[3]),,drop=FALSE]  #extra subsetting just in case there is more than one point
    graphics::points( x=worst[,1], y=worst[,2], col=grDevices::grey(0.5), pch=20, cex=0.1)
    graphics::points( x=worst[1,1], y=worst[1,2], col='black', pch=20, cex=0.25)
    raster::plot( IP.w, xlim=worst[1,1]+clusterRadius*c(-1,1), ylim=worst[1,2]+clusterRadius*c(-1,1), main="Working")
    graphics::points( worst[,1], worst[,2], pch=20, col='black')
    raster::plot( tmpRas, xlim=worst[1,1]+clusterRadius*c(-1,1), ylim=worst[1,2]+clusterRadius*c(-1,1), main="Observed - Specified")
    graphics::points( worst[,1], worst[,2], pch=20, col='black')
    graphics::mtext( paste0("Iter: ",kount), side=1, line=-2, outer=TRUE)
  }
  return( list( ras=crit.ras, max=crit[1], median=crit[2]))
}

"alterInclProbs.cluster" <- function( nCluster, clusterSize, clusterRadius, inclusion.probs, 
                                    maxIter=50, tolerance=NULL, mc.cores=parallel::detectCores()-1, doPlot=FALSE){
  #no precalculation (too memory hungry) but parallelisation (per row)
  
  #rescale to clusterSize
  raster::values( inclusion.probs) <- nCluster * clusterSize * raster::values( inclusion.probs) / sum( raster::values( inclusion.probs), na.rm=TRUE)
  
  if( is.null( tolerance)){
    tmp <- stats::median( raster::values( inclusion.probs)[raster::values( inclusion.probs)>0], na.rm=TRUE)
    tolerance <- 0.0001 * tmp
  }
  
  #message about task.
  message( "Optimisation summary.")
  message( "Maximum Iterations: ",maxIter,". Absolute tolerance: ",signif( tolerance,4),". Number of cores: ",mc.cores,".\n")
  
  #vectorised version of the spatial mask
  myMask <- is.na( raster::values( inclusion.probs))
  myMaskID <- which( myMask)
  #raster of cell numbers for extraction
  myCellIDs <- inclusion.probs
  raster::values( myCellIDs) <- 1:length( myCellIDs)
  raster::values( myCellIDs)[myMaskID] <- NA
  
  #weight matrix for focal calculations
  wtMat <- raster::focalWeight( x=inclusion.probs, d=clusterRadius, type="circle")
  wtMat[wtMat!=0] <- 1
  wtMatMask <- which( t( wtMat)==1)
  wtDim <- dim( wtMat)
  
  #setup cluster and count cells in neighbourhoods
  cl <- parallel::makeCluster(mc.cores)
  parallel::clusterExport( cl, c( "countFromRow", "getQrow", "getQcell"), envir = .getNamespace( "MBHdesign"))
  parallel::clusterExport( cl, c( "getValuesFocal"), envir=.getNamespace("raster"))
  countty <- as.numeric( parallel::parSapply( cl, 1:raster::nrow( inclusion.probs), countFromRow, wtMat=wtMat, inclusion.probs=inclusion.probs))
  countty[myMaskID] <- NA
  
  #Initialise the working values from linear approx
  IP.w <- inclusion.probs / ( clusterSize * countty)#sum( wtMat==1))
  raster::values( IP.w)[is.na( raster::values( IP.w))] <- 0
  raster::values( IP.w)[raster::values( inclusion.probs)==0] <- 0
  raster::values( IP.w)[myMaskID] <- NA
  IP.w.old <- IP.w
  #the first (pre-)iteration
  kount <- 0
  #take local averages
  IP.bar <- raster::raster( t( parallel::parSapply( cl, 1:raster::nrow( IP.w), IP.barFromRow, ras=IP.w, wtMat=wtMat)))
  raster::values( IP.bar)[myMaskID] <- NA
  
  IP.o <- rep( 0, length( IP.w))
  deriv.IP.o <- rep( NA, length( IP.w))
  
  QandDeriv <- parallel::parLapply( cl, 1:nrow( inclusion.probs), getQrow, IP.bar=IP.bar, IP.w=IP.w, clusterSize=clusterSize, wtMat=wtMat, myCellIDs=myCellIDs)
  QandDeriv <- lapply( QandDeriv, function(xx) t(xx))
  QandDeriv <- do.call( "rbind", QandDeriv)
  
  #avoid numerical errors, which when summed up can be a gotcha.
  #More importantly gives context to the working probs.
  QandDeriv <- scaleQ( Q=QandDeriv, tot=sum( raster::values( inclusion.probs), na.rm=TRUE), maskky=myMaskID)
  
  #create and possibly plot crit values
  crit <- getCrits( QandDeriv[,1], inclusion.probs, IP.w, doPlot, myMaskID, kount, clusterRadius)
  crit.old <- crit
  crit.old$ras <- rep( Inf, length( crit.old$ras))
  message( "Iteration: ", kount, ". Step Length: ", NA, ". Tolerance: ", signif( crit$max, 4), " (99.9 perc.) ", signif( crit$median, 4), " (median).")
  kount <- kount+1
  
  numNonZero <- sum( !is.na( crit.old$ras) & (crit.old$ras !=0), na.rm=TRUE)
  
  #the actual iteration now  
  while( kount <= maxIter & crit$max>tolerance){
    
    #adjust the previous update if necessary
    setWorse <- which( !is.na( crit.old$ras) & (crit.old$ras !=0) & abs( crit.old$ras) < abs( crit$ras))
    if( length( setWorse) > 0){
      message( "Previous step increases tolerance of ", length( setWorse), " (",signif( 100*length(setWorse)/numNonZero, 3), "% of all non-zero) cells. Reversing the update of those cells.")
      raster::values( IP.w)[setWorse] <- raster::values( IP.w.old)[setWorse]
    }
   
    
    #update from previous
    
    #step size (increasing as iterations progress)
    alpha <- 0.9 * ( 1 - exp( - kount / 5))
    
    #Pieces of the pie for updating
    Q <- ( QandDeriv[,1] - raster::values( inclusion.probs))^2
    Q.dash <- 2 * ( QandDeriv[,1] - raster::values( inclusion.probs)) * QandDeriv[,2]
    adj <- (Q/Q.dash)*alpha
  
    #catching those updates that are obviously stupid.
    tmpzero <- which( raster::values( IP.w.old) < adj)
    if( length( tmpzero) > 0){
      message( length( tmpzero), " Updates less than 0. Reducing those probabilities towards zero (by 0.5*step_to_zero).")
      adj[tmpzero] <- 0.5*raster::values( IP.w.old)[tmpzero]
    }
    tmpone <- which( adj > 1-raster::values( IP.w.old))
    if( length( tmpone) > 0){
      message( length( tmpone), " Updates larger than 1. Increasing those probabilities towards 1 by 0.5*step_to_one.")
      adj[tmpone] <- 0.5 * (1-raster::values( IP.w.old)[tmpone])
    }

    #actually do the update
    IP.w.old <- IP.w
    raster::values( IP.w) <- raster::values( IP.w.old) - adj
    
    #take local averages
    IP.bar <- raster::raster( t( parallel::parSapply( cl, 1:raster::nrow( IP.w), IP.barFromRow, ras=IP.w, wtMat=wtMat)))
    raster::values( IP.bar)[myMaskID] <- NA
    
    IP.o <- rep( 0, length( IP.w))
    deriv.IP.o <- rep( NA, length( IP.w))
    
    QandDeriv <- parallel::parLapply( cl, 1:nrow( inclusion.probs), getQrow, IP.bar=IP.bar, IP.w=IP.w, clusterSize=clusterSize, wtMat=wtMat, myCellIDs=myCellIDs)
    QandDeriv <- lapply( QandDeriv, function(xx) t(xx))
    QandDeriv <- do.call( "rbind", QandDeriv)
    
    #avoid numerical errors, which when summed up can be a gotcha.
    #More importantly gives context to the working probs.
    QandDeriv <- scaleQ( Q=QandDeriv, tot=sum( raster::values( inclusion.probs), na.rm=TRUE), maskky=myMaskID)
    
    #create and possibly plot crit values
    crit.old <- crit  #update irrespective of improvement or not
    crit <- getCrits( QandDeriv[,1], inclusion.probs, IP.w, doPlot, myMaskID, kount, clusterRadius)
    
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
  IP.w[myMask] <- NA
  QandDeriv[is.na( QandDeriv[,1]),1] <- 0
  QandDeriv[myMask,1] <- NA
  IP.o <- raster::rasterFromXYZ( cbind( coordinates( inclusion.probs), QandDeriv[,1]))
  IP.bar <- raster::raster( t( parallel::parSapply( cl, 1:raster::nrow( IP.w), IP.barFromRow, ras=IP.w, wtMat=wtMat)))
  raster::values( IP.bar)[myMaskID] <- NA
  tmp.bar <- IP.w
  raster::values( tmp.bar) <- raster::values( IP.bar)
  
  parallel::stopCluster(cl)
  
  res <- raster::stack(inclusion.probs, IP.o, IP.w, tmp.bar)
  names( res) <- c("IP.s","IP.o","IP.w","IP.bar")
  
  return( res)
}
