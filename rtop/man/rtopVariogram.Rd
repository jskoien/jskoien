\name{rtopVariogram}
\alias{rtopVariogram}
\alias{rtopVariogram.rtop}
\alias{rtopVariogram.SpatialPolygonsDataFrame}
\alias{rtopVariogram.SpatialPointsDataFrame}
\alias{rtopVariogram.STSDF}
\title{ create variogram for data with spatial support }
\description{
rtopVariogram will create binned variogram or cloud variogram of data
with an areal support.
}
\usage{
\method{rtopVariogram}{rtop}(object, params = list(), ...) 
\method{rtopVariogram}{SpatialPolygonsDataFrame}(object, ...) 
\method{rtopVariogram}{SpatialPointsDataFrame}(object, formulaString, params=list(), cloud, 
                       abins, dbins, ...) 
\method{rtopVariogram}{STSDF}(object, formulaString, params=list(), cloud, 
                       abins, dbins, data.table = FALSE, ...) 
}
\arguments{
\item{object}{ object of class \code{rtop} (see \link{rtop-package}) or a \cr
\code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}
or \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}}
with information about observations. If \cr
\code{object} is a \cr
\code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}},
it must have a column with name \code{area}.
}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables; suppose the dependent variable has name \code{z}, 
      for ordinary and simple kriging use the formula \code{z~1}; 
      for universal kriging, suppose \code{z} is linearly dependent on 
      \code{x} and \code{y}, use the formula \code{z~x+y}. The formulaString defaults
      to \code{"value~1"} if \code{value} is a part of the data set. 
      If not, the first column of the data set is used.}
\item{params}{a set of parameters, used to modify the default parameters for 
      the \code{rtop} package, set in \code{\link{getRtopParams}}. }
\item{cloud}{logical; if TRUE, calculate the semivariogram cloud, can be used
      to overrule the cloud parameter in params. }
\item{abins}{possibility to set areal bins (not yet implemented)}
\item{dbins}{possibility to set distance bins (not yet implemented)}
\item{data.table}{an option to use data.table internally for the variogram computation for 
        \code{\link[spacetime]{STSDF}}-objects}
\item{...}{parameters to other functions called, e.g. gstat's 
      \code{\link[gstat]{variogram}}-function and to \code{\link{rtopVariogram.SpatialPointsDataFrame}}
      when the method is called with an object of a different class}

} 

\value{ 
The function creates a variogram, either of type \code{rtopVariogram} or 
\code{rtopVariogramCloud}. This variogram is based on the \code{\link[gstat]{variogram}}
function from gstat, but has additional information about the spatial size or
length of the observations. An rtop-object with the variogram added is 
returned if the function is called with an rtop-object as argument.

For spatio-temporal objects (\code{\link[spacetime]{STSDF}}), the variogram is the spatially variogram, averaged for all time steps. There is a possibility to use data.table internally in this function, which can improve computation time for some cases.
}

\note{
The variogram cloud is similar to the variogram cloud from \code{\link[gstat]{gstat}},
with the area/length added to the resulting data.frame. The binned variogram is 
also based on the area or length, in addition to the distance between observations.
The bins equally distanced in the log10-space of the distances and areas (lengths).
The size of the bins is decided from the parameters \code{amul} and \code{dmul}, 
defining the number of bins per order of magnitude (1:10, 10:100, and so on).

The distances between areas are in this function based on the centre of gravity.

}

\references{ 
Skoien J. O., R. Merz, and G. Bloschl. Top-kriging - geostatistics on stream networks.
Hydrology and Earth System Sciences, 10:277-287, 2006.

Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione, A., 2014. Rtop: An R package for interpolation of data with a variable spatial support, with an example from river networks. Computers & Geosciences, 67.
}


\author{ Jon Olav Skoien }

\examples{
\dontrun{
library(sf)
rpath = system.file("extdata",package="rtop")
observations = st_read(rpath,"observations")
# Create a column with the specific runoff:
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

vario = rtopVariogram(observations, cloud = TRUE)
}
}
\keyword{spatial}
