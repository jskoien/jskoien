\name{optimizeNetwork}
\alias{optimizeNetwork}
\title{Optimization of networks}
\description{
Optimizes the sampling design of observation point locations using a varity of methods including spatial coverage 
by \code{k} means (as described in \code{\link[spcosa:spcosa-package]{spcosa}}) or by maximizing nearest neighbour distances 
and spatial simulated annealing (SSA, as described in \code{\link{ssaOptim}}) using MUKV as a criterion (\code{\link{calculateMukv}}) .
}
\usage{
optimizeNetwork(observations, predGrid, candidates, method, action,
                nDiff, model, criterion = "MUKV", plotOptim = TRUE, nGridCells, 
                nTry, nr_iterations = 10000, formulaString, fun, ...)}

\arguments{
\item{observations}{ object of class \code{\link[sp:Spatial-class]{Spatial}}*
               with coordinates and possible covariates }
\item{predGrid}{ object of class \code{\link[sp:Spatial-class]{Spatial}}*
               used when \code{method = "ssa"}. \code{predGrid} should contain the coordinates 
               of prediction locations for optimization. Usually predGrid is a 
               \code{\link[sp:SpatialGrid]{SpatialGrid}} / \code{\link[sp:SpatialGrid]{SpatialPixels}} 
               or a \code{\link[sp:SpatialGrid]{SpatialGridDataFrame}} /
               \code{\link[sp:SpatialGrid]{SpatialPixelsDataFrame}} when 
               independent covariate predictor variables are used }
\item{candidates}{ when \code{method = "manual"} or \code{method = "ssa"}, candidates is the study area of class 
              \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}; for other methods, 
              when \code{action = add}, candidates are points or polygons of class 
              \code{\link[sp:Spatial-class]{Spatial}}* }
\item{method}{ \code{"spcov"} for spatial coverage, 
              \code{"ssa"} for spatial simulated annealing or 
              \code{"manual"} for manual processing of the network }
\item{action}{ character string indicating which action to perform:
              \code{"add"} to add new measurement stations to the existing 
              network or \code{"del"} to delete existing stations }
\item{nDiff}{ number of stations to add or delete }
\item{model}{ variogram model to consider when \code{method = "ssa"} and \code{criterion = "mukv"}; 
              object of class \code{variogramModel},
              as generated by \code{\link[gstat:vgm]{vgm}} }
\item{criterion}{ Only in use for method \code{"ssa"}:  character string, \code{"mukv"} }
\item{plotOptim}{ logical; if TRUE, creates a plot of the result as optimization progresses; TRUE by default }
\item{nGridCells}{ when method is \code{"spcov"} and action is \code{"add"}: the approximate number gridcells 
             to explore within the candidate map as locations for new observations }
\item{nTry}{ when method is \code{"spcov"} and action is \code{"add"}: 
             \code{nTry} is the number of initial configurations to try. The method will keep the best solution in order 
             to reduce the risk of ending up with an unfavorable solution  }
\item{nr_iterations}{ number of iterations to process before stoping. The default coolingFactor in \code{\link{ssaOptim}}
             is also a function of number of iterations. Refer to \code{\link{ssaOptim}} for more details }
\item{formulaString}{ When \code{method = "ssa"}, this formula defines the dependent 
            variable as a linear model 
            of independent variables; suppose the dependent variable has name \code{z}, 
            for ordinary and simple kriging use the formula \code{z~1}; 
            for universal kriging, suppose \code{z} is linearly dependent on 
            \code{x} and \code{y}, use the formula \code{z~x+y}. The formulaString defaults
            to \code{"value~1"} if \code{value} is a part of the data set. 
            If not, the first column of the data set in \code{observations} is used. }
\item{fun}{Alternative objective function for optimization, the input and output should match
            the ones of  (\code{\link{calculateMukv}} (except for \code{fun})}
\item{...}{ other arguments to be passed on to lower level functions }
}

\details{
This function contains different methods to optimally add or remove point locations 
to or from a measurement network (Baume et al. 2011). 
Points can be added or deleted in the following ways: 
\enumerate{
\item{manually} 
\item{using a spatial coverage approach by \code{k} means to add stations
(as described in \code{\link[spcosa:spcosa-package]{spcosa}}, Brus et al. 2006) 
using a spatial coverage approach by maximizing mean nearest 
neighbour distances to remove stations (as described in \code{\link{spCovDel}})}
\item{or using spatial simulated annealing 
with mean universal kriging variance as a criterion (\code{\link{calculateMukv}}, 
Brus & Heuvelink 2007, Melles et al. 2011)}
}

The results of different methods can be checked using the function \code{\link{calculateMukv}}, 
which returns mean universal kriging variance for an optimized network. 
                                                         
The user should be aware of the following limitations:
\enumerate{
\item{ \code{method = "ssa"} is only implemented for \code{criterion = "mukv"} }
\item{ Input \code{candidates} should preferably be a continuous domain such 
           as \code{\link[sp:SpatialPolygons]{SpatialPolygons}} }
\item{ \code{method = "ssa"} with \code{criterion = "mukv"} makes it possible to assume a linear relationship between 
          independent variables in predGrid and dependent variables at observation locations using
          universal kriging (\code{\link[gstat:krige]{krige}}). However, a correct estimate of 
          mean universal kriging variance requires that the \code{independent} 
          covariate variables be known 
          at candidate locations. Thus it is necessary to have complete spatial 
          coverage for all covariate predictors
          in the linear model. Covariate information must be available at both 
          new candidate measurement locations and 
          prediction locations. This information is acquired (or sampled) from predGrid at 
          candidate locations during SSA using a call 
          to \code{\link[sp]{over}} by default. But see \code{\link{ssaOptim}} 
          for more details and an option to interpolate 
          these values for candidate locations from predGrid. }
\item{ Note that it is not recommended to use independent variables which differ strongly 
          in magnitude (as for traditional universal kriging) }
\item{ If no \code{formulaString} is supplied, an ordinary kriging formula is assumed, and 
          optimization will proceed using mean ordinary kriging variance }
}
}

\value{
Object of class \code{\link[sp:SpatialPoints]{SpatialPoints}}* with spatial coordinates 
of optimized locations (including observation locations when \code{action = "add"}) 
}

\references{

O. P. Baume, A. Gebhardt, C. Gebhardt, G. B. M. Heuvelink, J. Pilz (2011). Network 
optimization algorithms and scenarios in the context of automatic mapping, Computers and Geosciences, 
37: 289-294 (2011). 

S. J. Melles, G. B. M. Heuvelink, C. J. W. Twenhofel, U. Stohlker (2011).	Optimizing 
the spatial pattern of networks for monitoring radioactive releases, Computers and Geosciences, 
37: 280-288 (2011). 

D. J. Brus, G. B. M. Heuvelink (2007). Optimization of sample patterns for universal 
kriging of environmental variables, Geoderma, 138: 86-95 (2007).

D. J. Brus, J. de Gruijter, J. van Groenigen (2006). Designing spatial coverage samples using
the k-means clustering algorithm. In A. McBratney M. Voltz and P. Lagacherie,
editor, Digital Soil Mapping: An Introductory Perspective, Developments in Soil
Science, vol. 3., Elsevier, Amsterdam.

}

\author{ O. Baume, S.J. Melles, J. Skoien }

\seealso{ \code{\link{ssaOptim}}, \code{\link{spCovDel}}, \code{\link{spCovAdd}}, \code{\link{calculateMukv}}, 
\code{\link[spcosa:stratify-methods]{stratify}} }

\examples{
\donttest{
# load data:
library(gstat)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid) = TRUE
predGrid = meuse.grid

# estimate variograms (OK/UK):
vfitOK = fit.variogram(variogram(zinc~1, meuse), vgm(1, "Exp", 300, 1))
vfitUK = fit.variogram(variogram(zinc~x+y, meuse), vgm(1, "Exp", 300, 1))
vfitRK = fit.variogram(variogram(zinc~dist+ffreq+soil, meuse), vgm(1, "Exp", 300, 1))

# study area of interest:
bb = bbox(predGrid)
boun = SpatialPoints(data.frame(x=c(bb[1,1],bb[1,2],bb[1,2],bb[1,1],bb[1,1]),
                                y=c(bb[2,1],bb[2,1],bb[2,2],bb[2,2],bb[2,1])))
Srl = Polygons(list(Polygon(boun)),ID = as.character(1))
candidates = SpatialPolygonsDataFrame(SpatialPolygons(list(Srl)),
                                      data = data.frame(ID=1))


# add 20 more points assuming OK model (SSA method):
 optimOK <- optimizeNetwork(meuse, meuse.grid, candidates = candidates,
  method= "ssa", action= "add", nDiff = 20, model = vfitOK, criterion="MUKV",
  nr_iterations=10000, nmax=40)


# add 20 more points assuming UK model (SSA method):
optimUK <- optimizeNetwork(meuse, meuse.grid, candidates = candidates,
   method = "ssa", action = "add", nDiff = 20, model=vfitUK, criterion="MUKV",
   nr_iterations = 10000, nmax = 40, formulaString = zinc~x+y)

# add 20 more points with auxiliary variables (SSA method):
optimRK <- optimizeNetwork(meuse, meuse.grid, candidates=candidates,
   method="ssa", action="add", nDiff=4, model=vfitRK, criterion="MUKV",
   nr_iterations=10000, formula=zinc~dist+ffreq+soil, nmax=200)

# add optimally 20 stations from current network with method "spcov"
# (spatial coverage method)
optimSC = optimizeNetwork(meuse, meuse.grid, candidates, method = "spcov",
            action = "add", nDiff = 10, model = model, criterion = "MUKV", plotOptim = TRUE,
            nGridCells = 10000,nTry = 100 )

# delete optimally 10 stations from current network with method "manual"
if (interactive()) optimMAN = optimizeNetwork(meuse, meuse.grid, candidates, method = "manual",
            action = "del", nDiff = 10, model = model, criterion = "MUKV", plotOptim = TRUE )


# comparison of results with ordinary kriging variogram, otherwise add formulaString
# ssa method, assuming ordinary kriging
calculateMukv(optimOK, predGrid, vfitOK) 

# ssa method, using spatial location as covariates
calculateMukv(optimUK, predGrid, vfitUK, zinc~x+y)
 
# ssa method, using other variables as covariates
calculateMukv(optimRK, predGrid, vfitRK, zinc~dist+ffreq+soil) 

# spcov method
calculateMukv(optimSC, predGrid, vfitOK) 

# 10 stations manually deleted
if (interactive()) calculateMukv(optimMAN, predGrid, vfitOK, zinc~1) 

}
}

\keyword{spatial}


