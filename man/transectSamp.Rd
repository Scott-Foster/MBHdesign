\name{transectSamp}
\alias{transectSamp}
\alias{transectSamp.internal}
\title{Generates a spatial design for transects}
\description{For arbitrary transect patterns, generates a spatially balanced design for given inclusion probabilities over a grid of potential sampling locations}
\usage{
 transectSamp( n, study.area=NULL, potential.sites=NULL, 
		    inclusion.probs=NULL, control=NULL, constrainedSet=NULL)
 transectSamp.internal( n, study.area=NULL, potential.sites=NULL, 
		    inclusion.probs=NULL, control=NULL, constrainedSet, \dots)
}
\arguments{
\item{n}{the number of transect to provide a sample for}
\item{study.area}{a numeric matrix with dimension columns. This defines the sampling area from where the sites are selected -- each row defines a vertex of the sampling area and the order of rows is such that the vertices are joined in order. The last vertex is joined to the first. If NULL (default), the study.area is defined to be the smallest (hyper-)rectangle that bounds the potential.sites. If potential.sites is also NULL (default), then the study area is taken to be the unit (hyper-)square. This argument is closely related to potential.sites.}
\item{potential.sites}{a matrix (of size Nxdimension) of the spatial coordinates of the N sampling locations, of which n<<N are taken as the sample. If NULL (default) N=10000 samples are placed on a regular grid. If study.area is defined, then this grid is over the smallest bounding (hyper-)rectangle for the study.area. If study.area is NULL, the grid is over the unit (hyper-)square.}
\item{inclusion.probs}{a vector specifying the inclusion probability for each of the N potential sampling sites. This is the probability that each site will be included in the final sample. Locations are ordered the same as the potential.sites argument. If NULL (default) equal inclusion probabilities are specified. Note that there is no sum constraint to this vector -- inclusion probabilities are relative.}
\item{control}{A list containing details of the transects, their computational representation, and decisions about the design. See the Details for more information.}
\item{constrainedSet}{An N x nRotate Boolean matrix specifying whether the putative transect, starting at each location and orientated along each rotation, should be considered for sampling. See \code{\link{findDescendingTrans}} for an example.}
\item{\ldots}{Further arguments mostly for reducing computation for when certain quantities have been pre-calculated. Often not useful for general purpose, unless extra constraints are required for certain types of transects. The arguments include the computationally expensive objects (where repeated survey generation may be sped up): IDs -- the information about which grid cell each and every sampling location (along all transects) falls; transectsOverCells -- information about which sampling locations (from any transect) falls in all cells, and; adjustedSpecified -- The specified inclusion probabilities adjusted for edge effects.}
}
\details{This function aims to provide a randomised design for transects that: 1) adheres to inclusion probabilities for each grid cell; 2) is spatially coherent/balanced so that transects are well spread throughout the study area; and 3) is not too computationally expensive (or too devilishly difficult) to use effectively. However, fine grids over large areas will be demanding nevertheless.

The methods implemented by the function are outlined in Foster et al (2020), but are based on the ideas presented in Robertson et al (2013) and as implemented for Foster et al (2017). Briefly, the locations of the centre of the transects are chosen according to BAS sampling (Robertson et al; 2013) and the direction of the sample is then taken randomly. See Foster et al (2020) for further details.

The control list contains the following elements.
\itemize{
 \item transect.pattern: a set of points describing the transect pattern (a line, or an ``S'' pattern, or a dense grid, or a quadrat, or ...). This is a generic pattern, which will be rotated as required by the function for the randomisation. The pattern will be centred in both the x and y directions, so this does not need to be performed by the user (unless they choose to). Please specify as a 2 column matrix full of X- and Y-coordinates. A special value, and the default, is not a matrix at all. Rather the specific character string ``line'', in whcih case the function will produce its own generic pattern (for a line of length control$line.length with control$transect.nPts points along it). All values that are not ``line'' will be interpreted as a matrix. For symmetric patterns (around the centre point), and to avoid excess computations only, please make sure that an odd number of nRotate is used (see below).
 \item transect.nPts: The number of points to represent the (linear) transect. If control$transect.pattern is specified to something apart from ``line'' then this argument is ignored. A large value of transect.nPts will lead to a more accurate representation of the transects, but it will also cost more computation. The default is the arbitrary value of 20.
 \item line.length: If control$transect.pattern==``line'', then this control parameter gives the length of that line. By default, it is arbitrarily taken to be 0.1 times the maximal edge of the study.area.
 \item nRotate: The number of different bearings (through North, East, South, West) to consider at each sampling location for directions of the transect. The default is 11 (an odd number), meaning that transects laid out at 0 degrees, 32.7 degrees, 65.5 degrees, 98.2 degrees, ... (from north) will be considered. For symmetric a transect.pattern, such as a linear transect, nRotate should be an odd number (like the default) to avoid excess computation. This is because directions around a centre point can be inadvertently repeated (e.g. a linear transect in the 10 and 190 degree directions lead to the same trasnect). Increasing the number of directions to consider also increases computational time -- probably at very modest increase in efficiency or scientific validity.
 \item mc.cores: The number of computer cores to spread the (parallel) computation to. More decreases the real-world time that the function will use to execute, but this occurs at the expense of decreasing the performance of the computer to run other jobs. It is not recommended that the number of cores on the machine is exceeded. If NULL (the default), the function uses just 1 core. This is not generally recommend as it is the slowest choice.
 \item transAdjust: Should the user-specified inclusion probabilities be altered so that, in a transect design, each cell has that inclusion probability of being sampled. The default is TRUE -- the inclusion probabilities should be altered. See Foster et al (2020) for details.
 \item edge.max.iter: The maximum number of iterations to perform whilst adjusting inclusion probabilities for transect sampling. Interger >= 1. The default is 25.
 \item conv.tol.diff: the tolerance for convergence for adjusting inclusion probabilities for transect sampling. Default is 0.01, meaning that successive iterations have to have proportional performance criteria less than 0.0001 apart. 
 \item gamma: the dampening multiplier for successive iterations of the transect adjustment for inclusion probabilities. Default is 0.8. Values should be in [0,1] with larger values likely to be quicker and more unstable.
 \item calcObsProbs: Should the observed inclusion probabilities be calculated for the design. Default is TRUE if inclusion probabilities are being adjusted (transAdjust is TRUE), and FALSE otherwise (transAdjust is FALSE).
 \item return.index: Should the (computationally demanding) indexes be returned? These may be useful if multiple surveys designs are being generated; as done in a simulation study. Default is FALSE, and these large objects are not returned.
 \item spat.random.type: The type of randomisation used to select locations of transects. Default is ``quasi'', after Robertson et al (2013) which gives spatially-balanced designs. The other option is ``pseudo'', a random sample, which does no spatial balancing.
 \item nSampsToConsider: A parameter passed to \code{\link{quasiSamp}} controlling the acceptance sampling. Larger numbers are more computationally expensive, but are also more likely to find a design (important for very unequal inclusion probabilities and/or challengingly shaped survey areas). Default is 5000.
 \item nCells: If no potential sites are given, then some will be created on the unit square. There will be a nCells X nCells grid.
}
}
\value{The transectSamp function returns a list of two data.frames. The first data.frame (named ``transect'') contains the chosen transects' middle locations and their directions. The second data.frame (named ``points'') contains the set of points that the survey's transects are planned to pass through. These are points on each of the transects and are the points used to represent the transect during the design.}

\seealso{\code{\link{alterInclProbs}},\code{\link{quasiSamp},}, \code{\link{findDescendingTrans}}, \code{\link{findTransFromPoint}}}

\author{Scott D. Foster}
\references{
 Foster, S.D., Hosack, G.R., J. Monk, Lawrence, E., Barrett, N.S., Williams, A. and Przeslawski, P. (2020) Spatially-Balanced Designs for Transect-Based Surveys. Methods in Ecology and Evolution \emph{8}: 95--105.
 
 Foster, S.D., Hosack, G.R., Lawrence, E., Przeslawski, R., Hedge,P., Caley, M.J., Barrett, N.S., Williams, A., Li, J., Lynch, T., Dambacher, J.M., Sweatman, H.P.A, and Hayes, K.R. (2017) Spatially-Balanced Designs that Incorporate Legacy Sites. Methods in Ecology and Evolution \emph{8}:1433--1442.
 
 Robertson, B. L., Brown, J. A., McDonald, T. and Jaksons, P. (2013) BAS: Balanced Acceptance Sampling of Natural Resources. Biometrics \emph{69}: 776--784.
}
\examples{
\donttest{
#generate samples on a 50 x 50 grid
#Note that, although the random number is set, there may be differences between versions of R. 
#In particular, post R/3.6 might be different to R/3.5 and before
#jet plane
set.seed( 767)
#the number of potential sampling locations
N <- 50^2
#number of samples
n <- 10
#the grid on unit square
X <- as.matrix( expand.grid( 1:sqrt( N), 1:sqrt(N)) / sqrt(N) - 1/(2*sqrt(N)))
#the inclusion probabilities with gradient according to non-linear function of X[,1]
p <- 1-exp(-X[,1])
#standardise to get n samples
p <- n * p / sum( p)
#get the sample
#note that 5 points on the transect line and 5 directions considered is a bit thin.
#This low-definition is done to avoid trouble with CRAN's checks 
#(no example should take a long time to run).
samp <- transectSamp( n, potential.sites=X, inclusion.probs=p, 
		    control=list( transect.patter="line", nRotate=5, transect.nPts=5, mc.cores=7))
plot( samp$points[,5:6], main="n=10 TRANSECTS")
#tidy
rm( N, n, X, p, samp)
}
}
\keyword{misc}
