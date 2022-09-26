\name{quasiSamp}
\alias{quasiSamp}
\alias{quasiSamp.raster}
\title{Generates a spatial design using Quasi-random numbers}
\description{ Generates a spatially balanced design for given inclusion probabilities over a grid of potential sampling locations
}
\usage{
 quasiSamp( n, dimension=2, study.area=NULL, potential.sites=NULL, inclusion.probs=NULL,
                                              randStartType=2, nSampsToConsider=5000)
 quasiSamp.raster( n, inclusion.probs, randStartType=2, nSampsToConsider=5000)
}
\arguments{
\item{ n}{the number of samples to take}
\item{ dimension}{the number of dimensions that the samples are located in. Equal to 2 for areal sampling. Care should be taken with large dimensions as:1) the number of potential sampling sites needed for effective coverage starts to explode (curse of dimensionality); and 2) the well-spaced behaviour of the Halton sequence starts to deteriorate (but this requires very very many dimensions to be problematic -- included as a warning here for largely academic reasons).}
\item{ study.area}{a numeric matrix with dimension columns. This defines the sampling area from where the sites are selected -- each row defines a vertex of the sampling area and the order of rows is such that the vertices are joined in order. The last vertex is joined to the first. If NULL (default), the study.area is defined to be the smallest (hyper-)rectangle that bounds the potential.sites. If potential.sites is also NULL (default), then the study area is taken to be the unit (hyper-)square. This argument is closely related to potential.sites.}
\item{ potential.sites}{a matrix (of size Nxdimension) of the spatial coordinates of the N sampling locations, of which n<<N are taken as the sample. If NULL (default) N=10000 samples are placed on a regular grid. If study.area is defined, then this grid is over the smallest bounding (hyper-)rectangle for the study.area. If study.area is NULL, the grid is over the unit (hyper-)square.}
\item{ inclusion.probs}{either a vector or a RasterLayer specifying the inclusion probability for each of the N potential sampling sites. For quasiSamp.raster inclusion.probs must be a raster. The values contained are the probability that each site will be included in the final sample. If a vector, then the locations must be ordered the same as the potential.sites argument. If NULL and a vector (default) equal inclusion probabilities are specified. Must be specified for quasiSamp.raster. If quasiSamp is called with a RasterLayer value, then quasiSamp will internally call quasiSamp.raster and ignore the dimension, study.area and potential.sites arguments.}
\item{ randStartType}{the type of random start Halton sequence to use. The choices are 2 (default) as described in Robertson et al (2013), and 1 which is a mis-interpretation of that method (constrained so that the size of the skip in each dimension is equal). Note that randStartType=1 is used in Foster et al (2017).}
\item{ nSampsToConsider}{the total number of samples to consider in the BAS step (rejection sampling). The default is 5000, which means that 5000 halton numbers are drawn and then thinned according to the inclusion probabilities. You may want to increase this number if your inclusion probabilities are extremely unbalanced or if the number of samples required is close to 5000. Reduce if you want the code to run quicker and are confident that a sample will be found using less.}
}
\details{ These function are an implementation of the balanced adaptive sampling (BAS) designs presented in Robertson et al. (2013), which forms the basis for the methods in Foster et al (2017). The BAS approach uses Halton sequences of quasi-random numbers, which are evenly spread over space, as the basis for generating spatially balanced designs. In this implementation, we requrie that the inclusion probabilities be given as points in space and the BAS design is the set of these points that lie closest to a continuous-space Halton sequence. Computational speed has been rudimentily optimised, but (of course) it could be done better -- like coding outside of R{}, for example.

In an updated version of the package (Version 2.2.1 onwards) a RasterLayer can be passed to the function. This may be both more convenient and it will be faster for very large design problems. Note though that the underlying algorithm, and much of the code, remains unchanged between the two different versions.

}
\value{\item{ The quasiSamp and quasiSamp.raster functions returns a matrix of (dimension+2) columns. The first columns (of number dimension) are the sampled sites locations. The second to last column contains the inclusion probabilities for the sampled locations. The last column is the row number (of potential.sites) that corresponds to that sampled site.}{}
}
\seealso{\code{\link{alterInclProbs}}, \code{\link{modEsti}}
}
\author{Scott D. Foster}
\references{
 Robertson, B. L., Brown, J. A., McDonald, T. and Jaksons, P. (2013) BAS: Balanced Acceptance Sampling of Natural Resources. Biometrics \emph{69}: 776--784.
 
 Foster, S.D., Hosack, G.R., Lawrence, E., Przeslawski, R., Hedge,P., Caley, M.J., Barrett, N.S., Williams, A., Li, J., Lynch, T., Dambacher, J.M., Sweatman, H.P.A, and Hayes, K.R. (2017) Spatially-Balanced Designs that Incorporate Legacy Sites. Methods in Ecology and Evolution \emph{8}:1433--1442.
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
#Xp <- raster::rasterFromXYZ( cbind( X,p))
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
