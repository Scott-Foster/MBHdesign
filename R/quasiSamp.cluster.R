"extractWithBuffer" <- function( ras, cellIDs, width){
  pnts <- terra::as.points( ras, na.rm=FALSE)[cellIDs]
  buffs <- terra::buffer( pnts, width=width)
  extry <- terra::extract( ras, buffs, cells=TRUE)

  return( extry)
}

"quasiSamp.cluster"  <- function( nCluster, clusterSize, clusterRadius, 
                                inclusion.probs=NULL, working.inclusion.probs=NULL,
                                nSampsToConsider=c(25*nCluster,25*clusterSize), 
				nStartsToConsider=100*c(nCluster, clusterSize),
				randStartType=c(3,3),
				mc.cores=parallel::detectCores()-1){
  #inclusion.probs are not working inclusion probs, rather the raw target inclusion probs
  #working inclusion.probs are an object (raster stack) from alterInclProbs.cluster
  
  if( is.null( inclusion.probs) & is.null( working.inclusion.probs))
    stop( "One of the following arguments MUST be supplied: inclusion.probs working.inclusion.probs.  Please specify one (inclusion.probs will be ignored if both are given).")
  if( !is.null( inclusion.probs) & !is.null( working.inclusion.probs))
    message( "Both inclusion.probs AND working.inclusion.probs specified. Ignoring inclusion.probs.")
  if( !is.null( inclusion.probs))
    if( is.na( terra::crs( inclusion.probs)))
      message( "No projection specified for inclusion.probs, assuming lat/long (could cause an error later")
  if( !is.null( working.inclusion.probs))
    if( is.na( terra::crs( working.inclusion.probs)))
      message( "No projection specified for working.inclusion.probs, assuming lat/long (could cause an error later")
  if( length( nStartsToConsider) != 2){
    if( length( nStartsToConsider) == 1){
      message( "nStartsToConsider is a scalar when it should be a vector of length 2. Please consider, but continuing with the specified scalar replicated (for now)")
      nStartsToConsider <- rep( nStartsToConsider, 2)
    }
    else
      stop( "nStartsToConsider should be a integer vector of length 2.  See ?quasiSamp.cluster.")
  }
  if( length( nSampsToConsider) != 2){
    if( length( nSampsToConsider) == 1){
      message( "nSampsToConsider is a scalar when it should be a vector of length 2. Please consider, but continuing with the specified scalar replicated (for now)")
      nSampsToConsider <- rep( nSampsToConsider, 2)
    }
    else
      stop( "nSampsToConsider should be a integer vector of length 2.  See ?quasiSamp.cluster.")
  }
    
  if( is.null( working.inclusion.probs)){
    message( "No working.inclusion.probs specified. Calculating now (with default computational options.")
    working.inclusion.probs <- alterInclProbs.cluster( nCluster=nCluster, clusterSize=clusterSize, clusterRadius=clusterRadius, inclusion.probs=inclusion.probs, maxIter=500, tolerance=NULL, mc.cores=mc.cores, doPlot=FALSE)
  }

  #sample the cluster centres
  clusterDes <- quasiSamp.raster( nCluster, inclusion.probs=working.inclusion.probs$IP.bar, randStartType=randStartType[1], nSampsToConsider=nSampsToConsider[1], nStartsToConsider=nStartsToConsider[1])
  
  #the data for each cluster "swatch"
  clusterIPs <- terra::extract( working.inclusion.probs, clusterDes[,c("x","y")], buffer=clusterRadius, method='simple', cells=TRUE)
  tmp <- extractWithBuffer( ras=working.inclusion.probs, cellIDs=clusterIPs$cell, width=clusterRadius)
  #conditional probs
  tmp1 <- tapply( X=tmp$IP.w, INDEX=tmp$ID, FUN=function(xx) xx / sum( xx, na.rm=TRUE))
  tmp <- cbind( tmp, IP.cond=unlist( tmp1))
  #coordinates
  tmp <- cbind( terra::xyFromCell( working.inclusion.probs$IP.w, tmp$cell), tmp[,c("cell","ID", "IP.s","IP.bar","IP.cond","IP.w")])
  #turn into a list of SpatRast (first try is down-grading to previous versions of R, second (tapply) seems to only work in later versions. Dumbing down for now.
  tmpRast <- list()
  for( ii in 1:length( unique( tmp$ID)))
    tmpRast[[ii]] <- terra::rast( tmp[tmp$ID==unique( tmp$ID)[ii],c("x","y","cell","IP.s","IP.bar","IP.cond","IP.w")], type='xyz')
#  tmpRast <- tapply( X=tmp, INDEX=tmp$ID, FUN=function(xx) terra::rast( xx[,c("x","y","cell","IP.s","IP.bar","IP.cond","IP.w")], type='xyz'))
  #spatial sample in each cluster
  tmp2 <- lapply( 1:nCluster, function(xx) quasiSamp.raster(n=clusterSize, inclusion.probs = tmpRast[[xx]]$IP.cond, randStartType=randStartType[2], nSampsToConsider=nSampsToConsider[2], nStartsToConsider=nStartsToConsider[2]))
   
  #cell IDs and conditional probs
#  tmp <- lapply( 1:nCluster, function(xx) cbind( cell=clusterIPs[[xx]][,"cell"], IP.s=clusterIPs[[xx]][,"IP.s"], IP.bar=clusterIPs[[xx]][,"IP.bar"], condIP=clusterIPs[[xx]][,"IP.w"] / sum( clusterIPs[[xx]][,"IP.w"], na.rm=TRUE), IP.w=clusterIPs[[xx]][,"IP.w"]))
  
  #coordinates and conditional probs
#  tmp1 <- lapply( 1:nCluster, function(xx) cbind( terra::xyFromCell( working.inclusion.probs$IP.w, tmp[[xx]][,1]), tmp[[xx]][,c("cell","IP.s","IP.bar","condIP","IP.w")]))
  #turned into a list of rasters
#  tmp1Ras <- lapply( tmp1, function(xx) terra::rast( cbind( xx[,c('x','y','cell','IP.s','IP.bar','condIP','IP.w')]), type='xyz'))
  #spatial sample in each cluster
#  tmp2 <- lapply( 1:nCluster, function(xx) quasiSamp.raster(n=clusterSize, inclusion.probs = tmp1Ras[[xx]]$condIP, randStartType=3, nSampsToConsider=nSampsToConsider[2], nStartsToConsider=nStartsToConsider[2]))
  #re-append cell information
  tmp2a <- lapply( 1:nCluster, function(xx) cbind( terra::crds( tmpRast[[xx]], na.rm=FALSE), terra::values( tmpRast[[xx]], na.rm=FALSE))[tmp2[[xx]][,"ID"],])
  #append indexing
  tmp2a <- lapply( 1:nCluster, function(xx) cbind( tmp2a[[xx]], 1:clusterSize, xx))
  #flattening
  tmp3 <- do.call( "rbind", tmp2a)
  colnames( tmp3) <- c("x","y", "cellID", "IP.s", "IP.bar", "IP.cond", "IP.w", "point", "cluster")
  tmp3 <- as.data.frame( tmp3[,c("x","y","cluster","point","cellID","IP.s","IP.bar","IP.cond","IP.w")])
  
#  #converting to spatialPoints
#  tmp4 <- sp::SpatialPointsDataFrame( coords = as.matrix( tmp3[,c("x","y")]), data=tmp3[,-(1:2)], proj4string=sp::CRS( raster::proj4string(working.inclusion.probs)))
#  
#  tmpClust <- sp::SpatialPointsDataFrame( coords=as.matrix( clusterDes[,1:2]), data=clusterDes[,-(1:2)], proj4string=sp::CRS( raster::proj4string( working.inclusion.probs)))
#  attr( tmp4, "clusterDesign") <- tmpClust

  attr( tmp3, "clusterDesign") <- clusterDes

  return( tmp3)
}

