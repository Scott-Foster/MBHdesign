"quasiSamp.cluster"  <- function( nCluster, clusterSize, clusterRadius, 
                                inclusion.probs=NULL, working.inclusion.probs=NULL,
                                nSampsToConsider=c(25*nCluster,25*clusterSize), 
				nStartsToConsider=100*c(nCluster, clusterSize)){
  #inclusion.probs are not working inclusion probs, rather the raw target inclusion probs
  #working inclusion.probs are an object (raster stack) from alterInclProbs.cluster
  
  if( is.null( inclusion.probs) & is.null( working.inclusion.probs))
    stop( "One of the following arguments MUST be supplied: inclusion.probs working.inclusion.probs.  Please specify one (inclusion.probs will be ignored if both are given).")
  if( !is.null( inclusion.probs) & !is.null( working.inclusion.probs))
    message( "Both inclusion.probs AND working.inclusion.probs specified. Ignoring inclusion.probs.")
  if( !is.null( inclusion.probs) & is.na( raster::projection( inclusion.probs)))
    message( "No projection specified for inclusion.probs, assuming lat/long (could cause an error later")
  if( !is.null( working.inclusion.probs) & is.na( raster::projection( working.inclusion.probs)))
    message( "No projection specified for working.inclusion.probs, assuming lat/long (could cause an error later")
    
  if( is.null( working.inclusion.probs)){
    message( "No working.inclusion.probs specified. Calculating now (with default computational options.")
    working.inclusion.probs <- alterInclProbs.cluster( nCluster=nCluster, clusterSize=clusterSize, clusterRadius=clusterRadius, inclusion.probs=inclusion.probs, maxIter=500, tolerance=NULL, mc.cores=2, doPlot=FALSE)
  }

  #sample the cluster centres
  clusterDes <- quasiSamp.raster( nCluster, inclusion.probs=working.inclusion.probs$IP.bar, randStartType=3, nSampsToConsider=nSampsToConsider[1], nStartsToConsider=nStartsToConsider[1])
  
  #the data for each cluster "swatch"
  clusterIPs <- raster::extract( working.inclusion.probs, clusterDes[,c("x","y")], buffer=clusterRadius, method='simple', cellnumbers=TRUE)
  #cell IDs and conditional probs
  tmp <- lapply( 1:nCluster, function(xx) cbind( cell=clusterIPs[[xx]][,"cell"], IP.s=clusterIPs[[xx]][,"IP.s"], IP.bar=clusterIPs[[xx]][,"IP.bar"], condIP=clusterIPs[[xx]][,"IP.w"] / sum( clusterIPs[[xx]][,"IP.w"], na.rm=TRUE), IP.w=clusterIPs[[xx]][,"IP.w"]))
  #coordinates and conditional probs
  tmp1 <- lapply( 1:nCluster, function(xx) cbind( raster::xyFromCell( working.inclusion.probs$IP.w, tmp[[xx]][,1]), tmp[[xx]][,c("cell","IP.s","IP.bar","condIP","IP.w")]))
  #turned into a list of rasters
  tmp1Ras <- lapply( tmp1, function(xx) raster::rasterFromXYZ( cbind( xx[,c('x','y','cell','IP.s','IP.bar','condIP','IP.w')])))
  #spatial sample in each cluster
  tmp2 <- lapply( 1:nCluster, function(xx) quasiSamp.raster(n=clusterSize, inclusion.probs = tmp1Ras[[xx]]$condIP, randStartType=3, nSampsToConsider=nSampsToConsider[2], nStartsToConsider=nStartsToConsider[2]))
  #re-append cell information
  tmp2a <- lapply( 1:nCluster, function(xx) cbind( raster::coordinates( tmp1Ras[[xx]]), raster::values( tmp1Ras[[xx]]))[tmp2[[xx]][,"ID"],])
  #append indexing
  tmp2a <- lapply( 1:nCluster, function(xx) cbind( tmp2a[[xx]], 1:clusterSize, xx))
  #flattening
  tmp3 <- do.call( "rbind", tmp2a)
  colnames( tmp3) <- c("x","y", "cellID", "IP.s", "IP.bar", "IP.cond", "IP.w", "point", "cluster")
  tmp3 <- as.data.frame( tmp3[,c("x","y","cluster","point","cellID","IP.s","IP.bar","IP.cond","IP.w")])
  tmp4 <- sp::SpatialPointsDataFrame( coords = as.matrix( tmp3[,c("x","y")]), data=tmp3[,-(1:2)], proj4string=sp::CRS( proj4string(working.inclusion.probs)))
  
  tmpClust <- sp::SpatialPointsDataFrame( coords=as.matrix( clusterDes[,1:2]), data=clusterDes[,-(1:2)], proj4string=sp::CRS( proj4string( working.inclusion.probs)))
  attr( tmp4, "clusterDesign") <- tmpClust

  return( tmp4)
}

