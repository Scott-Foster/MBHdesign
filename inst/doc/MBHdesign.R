## ----prelim, echo = FALSE, results="hide"-------------------------------------
library( knitr)
opts_chunk$set(cache=TRUE, message = FALSE, comment = "", dev="pdf",
                      dpi=300, fig.show = "hold", fig.align = "center")

## ----setup1, eval=FALSE-------------------------------------------------------
#  install.packages( "MBHdesign")
#  ## or ##
#  devtools::install_github( repo="Scott-Foster/MBHdesign", build_vignettes=FALSE)

## ----setup2-------------------------------------------------------------------
library( MBHdesign)

## ----setSeed------------------------------------------------------------------
set.seed( 747)  #a 747 is a big plane

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
#write csv
write.csv( samp, file="pointSample1.csv", row.names=FALSE)
#tidy
rm( list=ls())

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
samp <- transectSamp( n=n, potential.sites=X, inclusion.probs=inclProbs,
		    control=my.control)
image( x=unique( X[,1]), y=unique( X[,2]),
    z=matrix( inclProbs, nrow=sqrt(nrow(X)), ncol=sqrt(nrow( X))),
    main="(Undadjusted) Inclusion Probabilities",
    sub="10 Transects",
    ylab=colnames( X)[2], xlab=colnames( X)[1])
points( samp$points[,5:6], pch=20, cex=0.6)

## ----transTidy----------------------------------------------------------------
#write csv
write.csv( samp$transect, file="transectSample1.csv", row.names=FALSE)
#tidy
rm( list=ls())

## ----volSetup, fig.width=9.43-------------------------------------------------
library( MASS)  #for the data
library( fields)  #for image.plot
library( MBHdesign) #for the spatial design and constraints
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
vol.constraints <- findDescendingTrans(
  potential.sites = pot.sites[,c("x","y")], bathy=pot.sites$height,
  in.area=rep( TRUE, nrow( pot.sites)), control=vol.control)
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
volSamp <- transectSamp( n=n, potential.sites=pot.sites[,c("x","y")],
                         control=vol.control,
                         constrainedSet=vol.constraints.bool)
#visualise the sample
image.plot( x=1:n.x, y=1:n.y, z=volcano,
            main="Uniform Probability Transect Sample", asp=1)
points( volSamp$points[,c("x","y")], pch=20)

## ----volTidy------------------------------------------------------------------
#write csv
write.csv( volSamp$transect, file="volcanoSample1.csv", row.names=FALSE)
#tidy
rm( list=ls())

## ----clusSetup----------------------------------------------------------------
#need raster functions
library( raster) 
#I'm currently on a Embraer 190, so it *almost* seems appropriate
set.seed( 747)  
#import the data and flip it so that latitudes index south to north
volcano <- datasets::volcano[nrow(datasets::volcano):1,]
#re-align the minimum altitude
volcano <- volcano - min( volcano)
#exaggerate altitude
volcano <- 50*volcano
#rescale for inclusion probabilities
#10 clusters of 5 sites
volcano <- 10*5*volcano / sum( volcano)
#cast to a raster
volcano <- raster( volcano, crs="+proj=utm +datum=WGS84 +units=km")

## ----clusterSamp--------------------------------------------------------------
#take the cluster sample
samp <- quasiSamp.cluster( nCluster=10, clusterSize=5, clusterRadius=0.05, 
				inclusion.probs = volcano)
#plot it over the volcano data
plot( volcano)
#the sample points
plot( samp, add=TRUE)
#the centres of the clusters 
#		(not sample points but potentially useful nevertheless)
plot( attr( samp, "clusterDes"), add=TRUE, pch=1, col='red')

## ----clusterOverSamp----------------------------------------------------------
#Create the working probabilties for the correct sized cluster.
workProbs <- alterInclProbs.cluster( nCluster=15, clusterSize=5, mc.cores=1,
			    clusterRadius=0.05, inclusion.probs=volcano)
#take the (over-sample)
overSamp <- quasiSamp.cluster( nCluster=15, clusterSize=10, 
		clusterRadius=0.05, working.inclusion.probs = workProbs)
#plot the results
par( mfrow=c(1,2))
plot( volcano, main="Planned sample and spare points within clusters")
#the planned sample
plot( overSamp[overSamp$cluster<=10 & overSamp$point<=5,], add=TRUE, cex=0.5)
#the over-sample (within clusters 1:10)
plot( overSamp[overSamp$cluster<=10 & overSamp$point>5,], add=TRUE, cex=0.5, col='red')
plot( volcano, main="Over-sampled clusters")
#the overs-sampled clusters (themselves oversampled)
plot( overSamp[overSamp$cluster>10 & overSamp$point<=5,], add=TRUE, cex=0.5)
plot( overSamp[overSamp$cluster>10 & overSamp$point>5,], add=TRUE, cex=0.5, col='red')

## ----clusterTidy--------------------------------------------------------------
#write csv
write.csv( as.data.frame( overSamp), file="clusterSamp1.csv", row.names=FALSE)
#tidy
rm( list=ls())

## ----sessionInfo, results = "asis", echo = FALSE------------------------------
toLatex(sessionInfo())

