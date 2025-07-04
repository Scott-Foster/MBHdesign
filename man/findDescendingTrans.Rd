\name{ConstrainTransects}
\alias{findDescendingTrans}
\alias{findTransFromPoint}
\title{Finds transects that are descending over complex terrain, or are stemming from a set of locations.}
\description{Creates a 'mask' for the extended sampling frame (sites by rotations) for two reasons: Some sampling gears that can only be deployed down-slope, and those designs where a 'wagon-wheel' is required from a small set of points (e.g. the top of a seamount).}
\usage{
  findDescendingTrans( potential.sites, bathy, in.area, descend.cutoff=0, control=NULL)
  findTransFromPoint( potential.sites, originPoints, in.area, control = NULL)
}
\arguments{
  \item{potential.sites}{a matrix (of size Nxdimension) of the spatial coordinates of the N sampling locations, of which n<<N are taken as the sample. If NULL (default) N=10000 samples are placed on a regular grid. If study.area is defined, then this grid is over the smallest bounding (hyper-)rectangle for the study.area. If study.area is NULL, the grid is over the unit (hyper-)square.}
  \item{bathy}{a vector of length N, giving the depth (or height) at each of the N potential.sites. Ordering must be the same as potential.sites.}
  \item{in.area}{a vector describing which cells are in the study region (TRUE) and which are not (FALSE). Ordering must be the same as potential.sites.}
  \item{descend.cutoff}{the value describing how down a descent actually is (in some cases slight rises are still manageable). Default is 0, meaning that strict descent is required.}
  \item{control}{A list containing details of the transects, their computational representation, and decisions about the design. See the Details Section of \code{\link{transectSamp}} for more information.}
  \item{originPoints}{a matrix with two columns giving the locations of the centre cells. There can be as many central cells as desired (more rows).}
}
\details{The function findDescendingTrans aims to provide a set of transects, which all can be performed down slope. The function findTransFromPoint aims to provide a set of transects, which all start from a particular point (or set of points). The output of these functions is formatted so that it can be directly used as an argument to \code{\link{transectSamp}}. In this case, and when supplied to \code{\link{transectSamp}}, a set of descending or radially stemming, spatially-balanced transects are generated.

This function implements a type of constraint on the randomisation process.  Discussed in Foster et al (2020).}

\value{The findDescendingTrans function returns a matrix giving a character matrix of size N x nRotate that specifies which transects are 'descend'ing (with N centre locations and nRotate different orientations), which transects are 'upAndDown' slope, which contains only NAs (where bathymetry data is not available), which are descending over the available bathymetry data ('descendAndNA'), which are up and down over the available data ('upAndDownAndNA_NaN').

The findTransFromPoint function similarly returns a character matrix containing the strings "startsFromPoint" meaning that the transect starts form one of the specified points, and "badStart" meaning that the transect does not start from one of the specified points. The developer realises that "badStart" is not an extremely informative string -- it must have seemed like a good idea late one night.}
  
\seealso{\code{\link{transectSamp}}}
\author{Scott D. Foster}
\references{
   Foster, S.D., Hosack, G.R., J. Monk, Lawrence, E., Barrett, N.S., Williams, A. and Przeslawski, P. (2020) Spatially-Balanced Designs for Transect-Based Surveys. Methods in Ecology and Evolution. 11: 95--105
}
