\name{modEsti}
\alias{modEsti}
\title{Get a model-based estimate of mean of a sampled area}
\description{ For a given survey design in any number of dimensions, calculate the mean prediction (plus SE plus 95\% CI) for the area.
}
\usage{
 modEsti( y, locations, includeLegacyLocation=TRUE, legacyIDs=NULL, predPts=NULL,
                               family=gaussian(), offset=rep(0,length(y)), control=list())
}
\arguments{
\item{ y}{a numeric vector of all observations (at new sites and legacy sites) collected in the survey}
\item{ locations}{a matrix (or something that can be coerced) containing the set of locations where observations were collected. Note that nrow( locations) == length( y).}
\item{ includeLegacyLocation}{a boolean indicating whether an extra term should be included into the model that corresponds to distance from legacy sites. See Foster et al (2017) for details. If TRUE (default) the extra term is included. If FALSE (so that model is purely spatial), then the extra term is discarded.}
\item{ legacyIDs}{the indexes, for the rows of y and locations, that correspond to legacy sites. For example, if the first, third and sixth rows were legacy sites in y and locations, then legacyIDs would be c(1,3,6).}
\item{ predPts}{A data.frame (or something that can be coerced to a data.frame) containing set of points to do the predictions at. Typically predPts is a dense grid (or similar) of points over the spatial area of interest. Note that the number of columns defines the number of dimensions. Do not include surplus variables in extra columns. If NULL (default), then a dense grid within the convex hull bounding "locations" will be used.}
\item{ family}{A family object giving the distribution of the data. Default is \code{gaussian()} but other sensible choices in ecology include \code{nb(link="log")} and \code{binomial()}. For more details see \code{?family.mgcv} from the mgcv package.}
\item{ offset}{A numeric vector of length equal to "length( y)". This gives any offset to the linear predictor of the model.}
\item{ control}{A list of control parameters (see details)}
}
\value{A list of three elements: 1) a point prediction of mean, 2) standard error of mean (obtained by parametric bootstrap), and 3) 95\% confidence interval of the mean.
}
\details{
This function works by fitting a generalised additive model (see \code{gam()}), predicting at the points predPts, and then averaging. Well, that is the general idea. The actual implementation uses a Monte Carlo routine to account for parameter uncertainty. This is done mirroring the helpfile of predict.gam. Basically, lots of sets of parameters are drawn from the parameters (asymptotic) distribution and then predictions are made for each draw. The overall estimate is then the mean (over parameter draws) of the mean (over prediction locations) of the prediction. Standard errors and confidence intervals are likewise calculated.

The control list contains elements with names:

\describe{
\item{k}{the number of knot points used in each dimension of locations}
\item{N}{the number of prediction points (in each dimension) for the grid, argument not used if "predPts!=NULL"}
\item{B}{the number of bootstrap samples to take of the parameter estimates}
\item{mc.cores}{the number of computer cores to spread the calculation of distances over (only used if includeLegacyLocation==TRUE)}
}
}
\seealso{\code{\link{quasiSamp}}, \code{\link{alterInclProbs}}, and \code{total.est} from the \pkg{spsurvey} package. The total.est function provides design-based estimation of both mean and variance.
}
\author{Scott D. Foster}
\references{
 Foster, S.D., Hosack, G.R., Lawrence, E., Przeslawski, R., Hedge,P., Caley, M.J., Barrett, N.S., Williams, A., Li, J., Lynch, T., Dambacher, J.M., Sweatman, H.P.A, and Hayes, K.R. (2017) Spatially-Balanced Designs that Incorporate Legacy Sites. Methods in Ecology and Evolution \emph{8}:1433--1442.
}
\examples{
#set up design parameters
#taken from the example in alterInclProbs()
#big plane today
set.seed(747)
#the number of potential sampling locations
N <- 50^2
#number of samples
n <- 27
#number of legacy sites
nLegacy <- 3
#the grid
X <- as.matrix( expand.grid( 1:sqrt( N), 1:sqrt(N)) / sqrt(N) - 1/(2*sqrt(N)))
#the inclusion probabiltiies with gradient according to non-linear function of X[,1]
p <- 1-exp(-X[,1])
#standardise to get n samples
p <- n * p / sum( p)
#randomly choose legacy sites
legacySites <- sample( 1:N, nLegacy, prob=p)
#alter inclusion probabilities for legacy sites
p2 <- alterInclProbs( legacy.sites=X[legacySites,], potential.sites=X, inclusion.probs=p)
#get the sample
samp <- quasiSamp( n=n, dimension=2, potential.sites=X, inclusion.probs=p2)
samp <- rbind( cbind( X[legacySites,], inclusion.probabilities=NA, ID=NA), samp)
#generate some fake data
samp$outcomes <- rnorm( nrow( samp))
#get the estimate
esti <- modEsti( y=samp$outcomes, locations=samp[,1:2], includeLegacyLocation=TRUE,
  legacyIDs=1:3, predPts=NULL, family=gaussian(), control=list(mc.cores=1, B=100))
#in real applications the number of bootstrap samples (B), and mc.cores, should be larger
print( esti)
#tidy
rm( esti, legacySites, n, N, nLegacy, p, p2, samp, X)
}
\keyword{misc}
