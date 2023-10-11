\name{quasiSamp}
\alias{quasiSamp}
\alias{quasiSamp.raster}
\alias{quasiSamp.cluster}
\title{Generates a spatial design using Quasi-random numbers}
\description{ Generates a spatially balanced design for given inclusion probabilities over a grid of potential sampling locations
}
\usage{
 quasiSamp( n, dimension=2, study.area=NULL, potential.sites=NULL, inclusion.probs=NULL,
                                randStartType=3, nSampsToConsider=25*n, 
				nStartsToConsider=100*n)
 quasiSamp.raster( n, inclusion.probs, randStartType=3, 
				nSampsToConsider=25*n, nStartsToConsider=100*n)
 quasiSamp.cluster( nCluster, clusterSize, clusterRadius, inclusion.probs=NULL, 
                                working.inclusion.probs=NULL, 
				nSampsToConsider=c(25*nCluster,25*clusterSize),
				nStartsToConsider=100*c(nCluster, clusterSize),
				randStartType=c(3,3),
				mc.cores=parallel::detectCores()-1)
}
\arguments{
\item{ n}{the number of samples to take.}
\item{ nCluster}{the number of clusters to sample. Only used when sampling sites in cluster groups.}
\item{ clusterSize}{the number of sites within each cluster to sample.}
\item{clusterRadius}{the radius of each individual clusters. Sites within a cluster will be chosen so that they are at most a distance of clusterRadius from the cluster centre. Must be specified in metres (see \code{\link[raster]{extract}} and \code{\link[raster]{buffer}}).}
\item{ dimension}{the number of dimensions that the samples are located in. Equal to 2 for areal sampling. Care should be taken with large dimensions as:1) the number of potential sampling sites needed for effective coverage starts to explode (curse of dimensionality); and 2) the well-spaced behaviour of the Halton sequence starts to deteriorate (but this requires very very many dimensions to be problematic -- included as a warning here for largely academic reasons).}
\item{ study.area}{a numeric matrix with dimension columns. This defines the sampling area from where the sites are selected -- each row defines a vertex of the sampling area and the order of rows is such that the vertices are joined in order. The last vertex is joined to the first. If NULL (default), the study.area is defined to be the smallest (hyper-)rectangle that bounds the potential.sites. If potential.sites is also NULL (default), then the study area is taken to be the unit (hyper-)square. This argument is closely related to potential.sites.}
\item{ potential.sites}{a matrix (of size Nxdimension) of the spatial coordinates of the N sampling locations, of which n<<N are taken as the sample. If NULL (default) N=10000 samples are placed on a regular grid. If study.area is defined, then this grid is over the smallest bounding (hyper-)rectangle for the study.area. If study.area is NULL, the grid is over the unit (hyper-)square.}
\item{ inclusion.probs}{either a vector or a SpatRaster specifying the inclusion probability for each of the N potential sampling sites. For quasiSamp.raster and for quasiSamp.cluster, inclusion.probs must be a raster. The values contained are the probability that each site will be included in the final sample. Note that inclusion.probs will be scaled internally so that they sum to the number of sites to be sampled. If a vector, then the locations must be ordered the same as the potential.sites argument. If NULL and a vector (default) equal inclusion probabilities are specified. Must be specified for quasiSamp.raster. If quasiSamp is called with a SpatRaster value, then quasiSamp will internally call quasiSamp.raster and ignore the dimension, study.area and potential.sites arguments. Argument is ignored in quasiSamp.cluster if (and only if) working.inclusion.probs is also supplied. If not ignored in quasiSamp.cluster, then \code{\link{alterInclProbs}} is called internally.}
\item{ working.inclusion.probs}{a SpatRaster of relevant inclusion probabilities. In particular, the working inclusion probabilities and their local (geogrphic) sums, which are uased in place of specified inclusion probabilities to ensure that cluster sampling respects the specified inclusion probabilities. See Foster et al (in review) for technical details. It is expected that, by far, the easiest way to generate this object is via a call to \code{\link{alterInclProbs.cluster}}. If working.inclusion.probs is not supplied then inclusion.probs must be. The function \code{\link{alterInclProbs.cluster}} internally performs the desired scaling of the inclusion probabilties, so this function (\code{\link{quasiSamp.cluster}}) does not perform any extra checks.}
\item{ randStartType}{the type of random start Halton sequence to use. For quasiSamp and quasiSamp.raster, the choices are 3 (default) as recommended in Robertson et al (2017), which improves the match between observed and specified inclusion probabilities (i.e. you get closer to what you want). Other options are 2 which gives the process in Robertson et al (2013), and 1 which is a mis-interpretation of method 2 (constrained so that the size of the skip in each dimension is equal). Note that randStartType=1 is used in Foster et al (2017). For quasiSamp.cluster, this argument is a integer vector of length two providing the random start type for both stages of cluster sampling (choosing clusters in first position and choosing sites within clusters in second position).}
\item{ nSampsToConsider}{the total number of samples to consider in the BAS step (rejection sampling). The default is 25*n, which means that 25*n halton numbers are drawn and then thinned according to the inclusion probabilities. Users will want to increase this number if inclusion probabilities are extremely unbalanced or if the number of samples required is close to, or exceeds, 25*n. Reduce if you want the code to run quicker and are confident that a sample will be found using less. For quasiSamp.cluster, nsampsToConsider is a vector of length 2. The first element specifies the number of samples to consider for sampling clusters. The second element specifies the number of samples for sampling within each cluster.}
\item{ nStartsToConsider}{(only used when randStartType=3). For quasiSamp and quasiSamp.raster: The maximum number of times the randomisation process should be performed before giving up. Default is 100*n. If this is not enough, then consider increasing it (and probably waiting longer for your computer to finish). For quasiSamp.cluster: a two element numeric vector giving the number of starts to attempt at each level of the randomisation process.}
\item{ mc.cores}{ When quasiSamp.cluster is called without a working.inclusion.probs argument (NULL), then \code{ \link{alterInclProbs.cluster}} is called with this many cores used.}
}
\details{ These function are an implementation of the balanced adaptive sampling (BAS) designs presented in Robertson et al. (2013) and Robertson et al. (2017). The former forms the basis for the methods in Foster et al (2017) and the latter is a modification of the former. The BAS approach uses Halton sequences of quasi-random numbers, which are evenly spread over space, as the basis for generating spatially balanced designs. In this implementation, we requrie that the inclusion probabilities be given as points in space and the BAS design is the set of these points that lie closest to a continuous-space Halton sequence. Computational speed has been rudimentily optimised, but (of course) it could be done better.

In an updated version of the package (Version 2.2.1 onwards) a raster can be passed to the function. Post 2.3.14, this must be a SpatRaster object from pacakge terra (prior a RasterLayer from pacakge raster). This may be both more convenient and it will be faster for very large design problems. Note though that the underlying algorithm, and much of the code, remains unchanged between the two different versions.

From version 2.3.0 onwards, the spatial cluster sampling approach of Foster et al (in review) is implemented in quasiSamp.cluster. This method proceeds in a two-staged fashion: cluster centres are chosen and then sites are chosen within clusterRadius of these centres. Both stages are chosen using quasi random numbers in BAS (Robertson et al; 2013).

In the edge case, where the number of samples is larger than the number of potential sampling points (or raster cells), the quasiSamp functions will simply sample sites multiple times. This behaviour may also be exhibited for cells with very high inclusion probabilties too, even when the sample size is larger than the number of potential sample sites. For cluster sampling, using \code{quasiSamp.cluster}, this also applies to within cluster sampling, as well as between cluster sampling.
}

\value{The quasiSamp and quasiSamp.raster functions returns a matrix of (dimension+2) columns. The first columns (of number dimension) are the sampled sites locations. The second to last column contains the inclusion probabilities for the sampled locations. The last column is the row number (of potential.sites) that corresponds to that sampled site.

The quasiSamp.cluster function returns a SpatialPointsDataFrame. It contains an identifier for cluster and site within cluster, the cellID from the original (inclusion.probs or working.inclusion.probs) raster, the specified inclusion probability for the cell, the cluster probability for a cluster centred at that cell, the conditional probability of sampling each cell within that cluster, and the working inclusion probabilities. The return object also contains a SpatialPointsDataFrame containing the design for the cluster.
}
\seealso{\code{\link{alterInclProbs}}, \code{\link{modEsti}}, \code{\link{alterInclProbs.cluster}}}

\author{Scott D. Foster}

\references{
 Robertson, B. L., Brown, J. A., McDonald, T. and Jaksons, P. (2013) BAS: Balanced Acceptance Sampling of Natural Resources. Biometrics \emph{69}: 776--784.
 
 Robertson, B.; McDonald, T.; Price, C. & Brown, J. (2017) A modification of balanced acceptance sampling Statistics and Probability Letters, \emph{129}, 107--112
 
 Foster, S.D., Hosack, G.R., Lawrence, E., Przeslawski, R., Hedge,P., Caley, M.J., Barrett, N.S., Williams, A., Li, J., Lynch, T., Dambacher, J.M., Sweatman, H.P.A, and Hayes, K.R. (2017) Spatially-Balanced Designs that Incorporate Legacy Sites. Methods in Ecology and Evolution \emph{8}:1433--1442.

 Foster, S.D., Lawrence, E., and Hoskins, A. (2023). Spatially Clustered Survey Designs. Journal of Agriculural, Biological and Environmental Statistics. 
}
\examples{
#generate samples on a 100 x 100 grid
#Note that, although the random number is set, there may be differences between versions of R. 
#In particular, post R/3.6 might be different to R/3.5 and before
#jet plane
set.seed(707)
#the number of potential sampling locations
N <- 100^2
#number of samples
n <- 10
#the grid on unit square
X <- as.matrix( expand.grid( 1:sqrt( N), 1:sqrt(N)) / sqrt(N) - 1/(2*sqrt(N)))
#the inclusion probabiltiies with gradient according to non-linear function of X[,1]
p <- 1-exp(-X[,1])
#standardise to get n samples
p <- n * p / sum( p)
#get the sample
samp <- quasiSamp( n=n, dimension=2, potential.sites=X, inclusion.probs=p)
par( mfrow=c(1,3), ask=FALSE)
plot( samp[,1:2], main="n=10")
#now let's get sillier
n <- 250
#get the sample
samp <- quasiSamp( n=n, dimension=2, potential.sites=X, inclusion.probs=p)
plot( samp[,1:2], main="n=250")
#silly or sublime?
n <- 1000
#get the sample
samp <- quasiSamp( n=n, dimension=2, potential.sites=X, inclusion.probs=p, nSampsToConsider=5000)
plot( samp[,1:2], main="n=1000")
#I'm sure that you get the idea now.

##The same for raster inclusion probabilities (just for illustration)
#Xp <- terra::rast( cbind( X,p), type='xyz')
#samp <- quasiSamp( n=10, inclusion.probs=Xp)
#plot( samp[,1:2], main="n=10 (raster)")
#samp <- quasiSamp( n=250, inclusion.probs=Xp)
#plot( samp[,1:2], main="n=250 (raster)")
#samp <- quasiSamp( n=1000, inclusion.probs=Xp, nSampsToConsider=5000)
#plot( samp[,1:2], main="n=1000 (raster)")

#tidy
rm( N, n, X, p, samp, Xp)
}
\keyword{misc}
