\name{alterInclProbs.cluster}
\alias{alterInclProbs.cluster}
\title{Alters inclusion probabilities to accommodate a 2-stage cluster sampling process}
\description{ Alters inclusion probabilities so that the cluster sampling approach respects specified inclusion probabilities.}
\usage{
 alterInclProbs.cluster( clusterSize, clusterRadius, inclusion.probs, 
                                    maxIter=50, tolerance=NULL, 
                                    mc.cores=parallel::detectCores()-1, doPlot=FALSE)
}
\arguments{
\item{ clusterSize}{an integer specifying the number of sites within a cluster.}
\item{ clusterRadius}{a numeric scalar giving the geographical extent of each cluster (from cluster centre to most extreme). Must be specified in metres (see \code{\link[raster]{extract}} and \code{\link[raster]{buffer}}).}
\item{ inclusion.probs}{a RasterLayer specifying the inclusion probability for each of the potential sampling sites. The values contained are the probability that each site will be included in the final sample. Must be specified.}
\item{maxIter}{an integer giving the maximum number of iterations that the algorithm is allowed to perform before termination.}
\item{tolerance}{a scalar numeric specifying the desired tolerance between observed inclusion probabilities (from 2-stage cluster sampling process) and the specified inclusion probabilities (from user). This tolerance specifies the accuracy with which the working inclusion probabilities are calculated.}
\item{ mc.cores}{an integer defining the number of processes to for some of the calculations on (in particular the calculation of distances to legacy sites). See parLapply(qv) in package parallel. Defaut is quite greedy, all but one of the cores available.}
\item{ doPlot}{a boolean indicating if a set of diagnostic plots should be produced within each iteration. Can be useful to help identify problems with algorithm convergence.}
}
\details{ Using the spatially-balanced 2-stage clustering process of Foster et al (in prep) will not, without careful consideration and adjustments, respect the inclusion probabilities specified by the user. This function provides the means to create _working_ inclusion probabilities, which when used in the 2-stage process _will_ produce the specified inclusion probabilies.
}
\value{The output of this function is designed for direct input into \code{\link{quasiSamp.cluster}}. It contains a rasterStack with four layers: the user-specified inclusion probabilities (IP.s), the observed inclusion probabilities (IP.o, those inclusion probabilities obtained when using 2-stage sampling with the working inclusion probabilities), the working inclusion probabilities and the area-summed working inclusion probabilities (IP.bar). The IP.bar are based on the working inclusion probabilities and are used to sample the centre of clusters.
}
\seealso{\code{\link{quasiSamp.cluster}}
}
\author{ Scott D. Foster}
\references{
 Foster, S.D., Lawrence, E., Others, TBA, and Hoskins, A. (in prep). Spatially-Balanced Cluster Sampling.
}
\keyword{misc}
