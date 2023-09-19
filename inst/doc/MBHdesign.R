## ----prelim, echo = FALSE, results="hide"-------------------------------------
library( knitr)
opts_chunk$set(cache=TRUE, message = FALSE, comment = "", dev="pdf",
                      dpi=300, fig.show = "hold", fig.align = "center")

## ----setup1, eval=FALSE-------------------------------------------------------
#  install.packages( "MBHdesign")
#  ## or ##
#  devtools::install_github( repo="Scott-Foster/MBHdesign", build_vignettes=FALSE)

## ----setup2, eval=TRUE--------------------------------------------------------
library( MBHdesign)

## ----setSeed------------------------------------------------------------------
set.seed( 747)  #a 747 is a big plane

## ----precompiled--------------------------------------------------------------
#TRUE if you want some of the examples shortened by using pre-saved output
usePrecompiledData <- FALSE

## ----legacySites--------------------------------------------------------------
#number of samples
n <- 10
#number of points to sample from
N <- 100^2
#the sampling grid (offset so that the edge locations have same area)
offsetX <- 1/(2*sqrt( N))
my.seq <- seq( from=offsetX, to=1-offsetX, length=sqrt(N))
X <- expand.grid( my.seq, my.seq)
#the legacy sites (three of them)
legacySites <- matrix( runif( 6), ncol=2, byrow=TRUE)
#names can be useful
colnames( X) <- colnames( legacySites) <- c("X1","X2")

## ----inclProbs, dpi=300, out.width='60%'--------------------------------------
#non-uniform inclusion probabilities
inclProbs <- 1-exp(-X[,1])
#scaling to enforce summation to n
inclProbs <- n * inclProbs / sum( inclProbs)
#uniform inclusion probabilities would be inclProbs <- rep( n/N, times=N)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]), 
    z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))), 
    main="(Undadjusted) Inclusion Probabilities", 
    ylab=colnames( X)[2], xlab=colnames( X)[1])
#The legacy locations
points( legacySites, pch=21, bg=grey(0.75), cex=1.5)

## ----alterInclProbs, dpi=300, out.width='60%'---------------------------------
#alter inclusion probabilities 
#   so that new samples should be well-spaced from legacy
altInclProbs <- alterInclProbs( legacy.sites=legacySites, 
		potential.sites=X, inclusion.probs = inclProbs)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]), 
    z=matrix( altInclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))), 
    main="Adjusted Inclusion Probabilities",
    ylab=colnames( X)[2], xlab=colnames( X)[1])
#The legacy locations
points( legacySites, pch=21, bg=grey(0.75), cex=1.5)

## ----GetDesign, dpi=300, out.width='60%'--------------------------------------
#generate the design according to the altered inclusion probabilities.
samp <- quasiSamp( n=n, dimension=2, 
	study.area=matrix( c(0,0, 0,1, 1,1, 1,0),ncol=2,  byrow=TRUE), 
	potential.sites=X, inclusion.probs=altInclProbs)
#for faster sampling (large problems), 
#	consider using quasiSamp.raster, which utilises SpatRaster formats
#visualise
image( x=unique( X[,1]), y=unique( X[,2]), 
    z=matrix( altInclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))), 
    main="Adjusted Inclusion Probabilities",
    ylab=colnames( X)[2], xlab=colnames( X)[1])
#The legacy locations
points( legacySites, pch=21, bg=grey(0.75), cex=1.5)
points( samp[,1:2], pch=21)

## ----ShowDesign---------------------------------------------------------------
print( samp, row.names=FALSE)

## ----GetData------------------------------------------------------------------
#generate some `observations' for the new sites
Z <- 3*( X[samp$ID,1]+X[samp$ID,2]) + 
			sin( 6*( X[samp$ID,1]+X[samp$ID,2]))
#and some for the legacy sites
Zlegacy <- 3*( legacySites[,1]+legacySites[,2]) + 
			sin( 6*( legacySites[,1]+legacySites[,2]))

## ----HTestimate---------------------------------------------------------------
#the proportion of legacy sites in the whole sample
fracLegacy <- nrow( legacySites) / (n+nrow( legacySites))
#inclusion probabilities for legacy sites
#   (these are just made up, from uniform)
LegInclProbs <- rep( nrow( legacySites) / N, nrow( legacySites))
#estimator based on legacy sites only
legacyHT <- (1/N) * sum( Zlegacy / LegInclProbs)
#estimator based on new sites only
newHT <- (1/N) * sum( Z / inclProbs[samp$ID])
mean.estimator <- fracLegacy * legacyHT + (1-fracLegacy) * newHT
#print the mean
print( mean.estimator)

## ----NNestimate---------------------------------------------------------------
#load the spsurvey package
library( spsurvey)
#rescale the inclusion probs
#   (the sample frames are the same in legacy and new sites)
tmpInclProbs <- ( c( inclProbs[samp$ID], LegInclProbs) / n) *
						(n+nrow(legacySites))
#create a temporary data frame
tmpDat <- data.frame( siteID=
		      c( samp$ID, paste0( "legacy", 1:nrow(legacySites))),
                      wgt=1/tmpInclProbs,
                      xcoord=c(samp$X1,legacySites[,1]),
                      ycoord=c(samp$X2,legacySites[,2]), Z=c(Z,Zlegacy))
#calculate the standard error
se.estimator <- cont_analysis( tmpDat, vars="Z",
                      weight="wgt",
                      xcoord="xcoord",
                      ycoord="ycoord")$Mean$StdError
#print it
print( se.estimator)

## ----ModEstimate--------------------------------------------------------------
tmp <- modEsti( y=c( Z, Zlegacy), locations=rbind( X[samp$ID,], legacySites),
	includeLegacyLocation=TRUE, legacyIDs=n + 1:nrow( legacySites),
	predPts=X, control=list(B=1000))
print( tmp)

## ----Tidy---------------------------------------------------------------------
##write csv if wanted
#write.csv( samp, file="pointSample1.csv", row.names=FALSE)
#tidy
rm( list=c( "samp", "tmp", "se.estimator","tmpDat","tmpInclProbs",
	    "mean.estimator","newHT","legacyHT","LegInclProbs",
	    "fracLegacy","Zlegacy","Z","X","altInclProbs", "n", "N", 
	    "my.seq","inclProbs", "offsetX", "legacySites"))

## ----transSetup---------------------------------------------------------------

set.seed( 747)  #I'm currently on a 787, so it *almost* seems appropriate
#number of transects
n <- 10
#number of points to sample from
N <- 100^2
#the sampling grid (offset so that the edge locations have same area)
offsetX <- 1/(2*sqrt( N))
my.seq <- seq( from=offsetX, to=1-offsetX, length=sqrt(N))
X <- expand.grid( my.seq, my.seq)
colnames( X) <- c("X1","X2")

## ----transIinclProbs, dpi=300, out.width='60%'--------------------------------
#non-uniform inclusion probabilities
inclProbs <- 1-exp(-X[,1])
#scaling to enforce summation to n
inclProbs <- n * inclProbs / sum( inclProbs)
#uniform inclusion probabilities would be inclProbs <- rep( n/N, times=N)
#visualise
image( x=unique( X[,1]), y=unique( X[,2]),
    z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
    main="(Undadjusted) Inclusion Probabilities",
    ylab=colnames( X)[2], xlab=colnames( X)[1])

## ----transSetControl----------------------------------------------------------
#my.control is a list that contains
my.control <- list(
  #the type of transect
  transect.pattern="line",
  #the length of transect
  line.length=0.15,
  #the number of points that define the transect
  transect.nPts=15,
  #the number of putative directions that a transect can take
  nRotate=9
)

## ----callTransectSamp---------------------------------------------------------
#take the transect sample
if( !usePrecompiledData){
  samp <- transectSamp( n=n, potential.sites=X, 
                inclusion.probs=inclProbs, control=my.control)
} else{
  samp <- readRDS( system.file(
		    "extdata", "transectSamp1.RDS", 
		    package="MBHdesign"))}
image( x=unique( X[,1]), y=unique( X[,2]),
    z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
    main="(Undadjusted) Inclusion Probabilities",
    sub="10 Transects",
    ylab=colnames( X)[2], xlab=colnames( X)[1])
points( samp$points[,5:6], pch=20, cex=0.6)

## ----transTidy----------------------------------------------------------------
##write csv
#write.csv( samp$transect, file="transectSample1.csv", row.names=FALSE)
#tidy
rm( list=c( "X", "inclProbs", "samp", "my.control", 
	    "my.seq", "offsetX", "n", "N"))

## ----volSetup, fig.width=9.43-------------------------------------------------
library( MASS)  #for the data
library( fields)  #for image.plot
#library( MBHdesign) #for the spatial design and constraints
set.seed( 717)  #Last plan I was on
#number of transects
n <- 20
#load the altitude data
data( volcano)  #this is a matrix
n.x <- nrow( volcano)
n.y <- ncol( volcano)
image.plot( x=1:n.x, y=1:n.y, z=volcano, main="Mountain Height (m)", asp=1)
#format for MBHdesign functions
pot.sites <- expand.grid( x=1:n.x, y=1:n.y)
pot.sites$height <- as.vector( volcano)
#details of the transects (see Details section in ?transectSamp)
vol.control <- list( transect.pattern="line", transect.nPts=10,
                     line.length=7, nRotate=11, mc.cores=1)
#In a real application, transect.nPts and nRotate may need to be increased
#1 cores have been used to ensure generality for all computers. 
#		Use more to speed things up

## ----volConstraint, fig.width=9.43--------------------------------------------
if( !usePrecompiledData){
  vol.constraints <- findDescendingTrans(
        potential.sites = pot.sites[,c("x","y")], 
        bathy=pot.sites$height, in.area=rep( TRUE, nrow( pot.sites)), 
        control=vol.control)
} else{
  vol.constraints <- readRDS( system.file(
		    "extdata", "transectConstraints1.RDS", 
		    package="MBHdesign"))}
#this is a matrix with nrow given by the number of sites and ncol by
#   the number of rotations around each site
print( dim( vol.constraints))
#The contents describe how the transect lays over the landscape
#So, there are 15592 putative transects that ascend and descend
#   (and can't be used in the sample)
table( as.vector( vol.constraints))
#convert to TRUE/FALSE
#Note that the final possible transect type ('descendAndNA') is
#   not present in these data
#If present, we would have to decide to sample these or not
vol.constraints.bool <- matrix( FALSE, nrow=nrow( vol.constraints),
                                ncol=ncol( vol.constraints))
vol.constraints.bool[vol.constraints %in% c("descend")] <- TRUE
#Let's get a visual to see what has just been done.
tmpMat <- matrix( apply( vol.constraints.bool, 1, sum), nrow=n.x, ncol=n.y)
image.plot( x=1:n.x, y=1:n.y, z=tmpMat,
            main="Number of Transects",
            sub="Transects centered at cell (max 11)", asp=1)
#There aren't any transects that are centred on ridges or depressions.

## ----volSample, fig.width=9.43------------------------------------------------
#take the sample
if( !usePrecompiledData){
  volSamp <- transectSamp( n=n, potential.sites=pot.sites[,c("x","y")],
                         control=vol.control,
                         constrainedSet=vol.constraints.bool)
} else{
  volSamp <- readRDS( system.file(
		    "extdata", "transectSamp2.RDS", 
		    package="MBHdesign"))}
#visualise the sample
image.plot( x=1:n.x, y=1:n.y, z=volcano,
            main="Uniform Probability Transect Sample", asp=1)
points( volSamp$points[,c("x","y")], pch=20)

## ----volTidy------------------------------------------------------------------
##write csv
#write.csv( volSamp$transect, file="volcanoSample1.csv", row.names=FALSE)
#tidy
rm( list=c( "volSamp", "tmpMat", "vol.constraints.bool", 
   "vol.constraints", "vol.control", "pot.sites", 
   "n", "n.x", "n.y", "volcano"))

## ----clusSetup----------------------------------------------------------------
#need raster functions
library( terra) 
#import example data
egDat <- rast(system.file(
		    "extdata", "ACT_DemoData.grd", 
		    package="MBHdesign"))$soilMoisture
values( egDat) <- ( values( egDat) - min( values( egDat), na.rm=TRUE)) * 5

## ----clusterSamp--------------------------------------------------------------
set.seed( 727)  
#take the cluster sample
#increase mc.cores for faster processing
if( !usePrecompiledData){
  samp <- quasiSamp.cluster( nCluster=10, clusterSize=5, clusterRadius=5, 
				inclusion.probs = egDat, mc.cores=1)
} else{
  samp <- readRDS( system.file(
		    "extdata", "clusterSamp1.RDS", 
		    package="MBHdesign"))}
#plot it over the egData data
plot( egDat)
#the sample points
points( samp$x, samp$y, pch=20, cex=0.5)
#the centres of the clusters 
#		(not sample points but potentially useful nevertheless)
points( attr( samp, "clusterDes")[,c("x","y")], pch=1, col='red', cex=0.5)

## ----clusterOverSamp----------------------------------------------------------
#Create the working probabilties for the correct sized cluster.
if( !usePrecompiledData){
  workProbs <- alterInclProbs.cluster( nCluster=15, clusterSize=5, 
                mc.cores=1, clusterRadius=5, inclusion.probs=egDat)
} else{
  workProbs <- readRDS( system.file(
		    "extdata", "clusterWorkProbs1.RDS", 
		    package="MBHdesign"))}
#take the (over-sample)
set.seed( 747)
overSamp <- quasiSamp.cluster( nCluster=15, clusterSize=10, 
		clusterRadius=5, working.inclusion.probs = workProbs)
#plot the results
par( mfrow=c(1,2))
plot( egDat, main="Planned and Spare points")
#the planned sample
points( overSamp[overSamp$cluster<=10 & overSamp$point<=5,c("x","y")], cex=0.5)
#the over-sample (within clusters 1:10)
points( overSamp[overSamp$cluster<=10 & overSamp$point>5,c("x","y")], 
		cex=0.5, col='red')
plot( egDat, main="Over-sampled clusters")
#the overs-sampled clusters (themselves oversampled)
points( overSamp[overSamp$cluster>10 & overSamp$point<=5,c("x","y")], cex=0.5)
points( overSamp[overSamp$cluster>10 & overSamp$point>5,c("x","y")], 
		cex=0.5, col='red')

## ----clusterTidy--------------------------------------------------------------
##write csv
#write.csv( as.data.frame( overSamp), 
#      file="clusterSamp1.csv", row.names=FALSE)
#tidy
rm( list=c("egDat","overSamp","workProbs","samp","usePrecompiledData"))

## ----sessionInfo, results = "asis", echo = FALSE------------------------------
toLatex(sessionInfo())

