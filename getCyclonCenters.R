#     getCyclonCenters.R Algorithm to obtain cyclon centers based on: 
#                        http://www.europeanwindstorms.org/method/track/
#
#     Copyright (C) 2018 Santander Meteorology Group (http://www.meteo.unican.es)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' @title Cyclon Centers Detection
#' @description Implementation of an algorithm to detect the cyclon centers
#'
#' @param slp Sea Level Pressure (hPa).
#' @param vo  Relative Vorticity (s-1).
#' @param seek.radius Maximum distance in degrees to look for a new cyclon center to be included in the same trajectory. The default value is 10.
#' @param slp.diff.threshold Minimum difference between the slp of a grid box and its neighbours in order to 
#' be considered a local minimum and, then, a possible cyclon center. The default value is 0.
#' @param vo.diff.threshold Minimum difference between the vorticity of a grid box and its neighbours in order
#' to be considered a local minimum and, then, a possible cyclon center. The default value is 0.
#' @param lap.diff.threshold Minimum difference between the laplacian of the sea level pressure of a grid box 
#' and its neighbours in order to be considered a local minimum and, then, a possible cyclon center. The default value is 20.
#' @param ndr.threshold Minimum value of the Normalized Deepeningun Ratio to be considered a cyclon center as explosive. The default value is 2.5.
#' @param vo.threshold Minimum value of vorticity to be considered a candidate to cyclon center. The default value is 1e-4.
#' @param wss Wind speed at the same period and geographical domain to be included in the output. The default value is \code{NULL}.
#' @param criteria Criteria considered to select the cyclon center between the possible candidates. Options implemented are:
#' \code{"global"} (default): maximum of the Euclidean norm of a vector containing the vorticity, the NDR and the slp normalized;
#' \code{"max.vo"}: point with the maximum vorticity; \code{"max.ndr"}: point with the maximum NDR; \code{"min.slp"}: point with the minimum slp.
#'     
#' @seealso \code{\link{getCyclonTrack}} for the algorithm to define the cyclon tracks from the cyclon centers detected. 
#' @return A list with the cyclon centers detected for each time step. Each row corresponds to the NDR, slp, vorticity, 
#' laplacian of the slp, longitude, latitude, windspeed and explosive character of the cyclon center.
#' @family cyclonTrack
#'
#' @references
#'
#' \itemize{
#' \item I. Iparragirre (2018) Climate change projections of explosive cyclogenesis events frequency in the Iberian Peninsula. Master's Thesis: Physics, Instrumentation and the Environment - Universidad de Cantabria, http://meteo.unican.es/en/theses/2018_Itsasne
#' 
#' }
#' @author I. Iparragirre, M. D. Frías, S. Herrera
#' @export
#' @examples {
#' ## We define the geographical and temporal domain to be loaded.
#' lonLim <- c(-50,40)
#' latLim <- c(15,75)
#' r <- 6378.1 # Parameter: Earth radius
#' season <- 1:12
#' years <- c(2009:2010)
#' 
#' ## We consider the ERA-Interim dataset from the Santander User Data Gateway
#' dataset <- "http://meteo.unican.es/tds5/dodsC/interim/daily/interim20_daily.ncml"
#' di <- dataInventory(dataset = dataset)
#' names(di)
#' 
#' ## Sea Level Pressure
#' slp <- loadGridData(dataset = dataset, var = "SLP", season = season, years = years,
#'                     lonLim = lonLim, latLim = latLim, time = "none", aggr.d = "none")
#' ## The input of the algorithm should be expressed in hPa instead Pa:
#' ao <- c(0,0.01)
#' slp$Data <- ao[2]*(ao[1]+slp$Data)
#' 
#' ## Vorticity at 850 hPa
#' zg <- loadGridData(dataset = dataset, var = "Z850", season = season, years = years,
#'                    lonLim = lonLim, latLim = latLim, time = "none", aggr.d = "none")
#' vo <- laplacian(zg)
#' rm("zg")
#' 
#' ## Detection and Tracking parameters
#' seek.radius <- 6
#' slp.diff.threshold <- 10
#' vo.diff.threshold <- 1e-6
#' lap.diff.threshold <- 20
#' ndr.threshold <- 1
#' vo.threshold <- 1e-5
#' criteria <- "global"
#' 
#' ## Detection and Tracking
#' varOutputCenters <- getCyclonCenters(slp,vo, seek.radius = seek.radius, slp.diff.threshold = slp.diff.threshold,
#'                                      vo.diff.threshold = vo.diff.threshold, lap.diff.threshold = lap.diff.threshold,
#'                                      ndr.threshold = ndr.threshold, vo.threshold = vo.threshold, criteria = criteria)
#' 
#' ## Plotting the results
#' data("wrld")
#' po <- SpatialPoints(cbind(cyclonTrack[[1]][[2]][,5], cyclonTrack[[1]][[2]][,6]))
#' dat <- as.data.frame(cyclonTrack[[1]][[2]][,2])
#' colnames(dat) <- "y"
#' kl <- SpatialPointsDataFrame(po, data = dat)
#' spplot(kl, zcol = "y", sp.layout = list(wrld, first = F), colorkey = TRUE, 
#'               xlim = lonLim, ylim = latLim, xlab = slp$Dates$start[1])
#'}


getCyclonCenters <- function(slp,
                              vo,
                              seek.radius = 10,
                              slp.diff.threshold = 0,
                              vo.diff.threshold = 0,
                              lap.diff.threshold = 20,
                              ndr.threshold = 2.5,
                              vo.threshold = 1e-4,
                              criteria = "global",
                              wss = NULL){

  ## Definition of the parameters:
  slp.threshold <- 2000 # bound for the SLP to apply a maximization problem.
  r <- 6378.1 # Earth radius
  g <- 9.8 # Earth gravity
  omega <- 7.2921e-5 # Earth angular velocity
  pdiff <- 1e-6 ## pressure threshold (mb) for identifying cyclones.
  maxdist <- 800 ## maximum search distance (km) for tracking storms between time steps. 
  pchange <- 2*1e-5 ## limit on allowable absolute pressure tendency (mb) between time steps. 
  #######################################################################
  ## lat/lon positions of the EASE grid to be interpolated to.
  ieqlon <- length(slp$xyCoords$x)
  jeqlat <- length(slp$xyCoords$y)
  ## Offset and scale of the data
  lat <- matrix(data = slp$xyCoords$y, nrow = length(slp$xyCoords$y), ncol = length(slp$xyCoords$x))
  long <- t(matrix(data = slp$xyCoords$x, nrow = length(slp$xyCoords$x), ncol = length(slp$xyCoords$y)))
  lap <- laplacian(slp)
  ao <- c(0,1e5) ## Change the units to hPa
  lap <- ao[2]*(ao[1]+lap$Data)
  firstDate <- strsplit(vo$Dates$start[1],' ')
  yymmdd <- strsplit(firstDate[[1]][1],'-')
  hhmmss <- strsplit(firstDate[[1]][2],':')
  firstDate <- c(as.double(yymmdd[[1]][1]),as.double(yymmdd[[1]][2]),as.double(yymmdd[[1]][3]),as.double(hhmmss[[1]][1]),as.double(hhmmss[[1]][2]),as.double(hhmmss[[1]][3]))
  secondDate <- strsplit(vo$Dates$start[2],' ')
  yymmdd <- strsplit(secondDate[[1]][1],'-')
  hhmmss <- strsplit(secondDate[[1]][2],':')
  secondDate <- c(as.double(yymmdd[[1]][1]),as.double(yymmdd[[1]][2]),as.double(yymmdd[[1]][3]),as.double(hhmmss[[1]][1]),as.double(hhmmss[[1]][2]),as.double(hhmmss[[1]][3]))
  timeStep <- sum(c(24*365*(secondDate[1] - firstDate[1]), 24*30*(secondDate[2] - firstDate[2]), 24*(secondDate[3] - firstDate[3]), (secondDate[4] - firstDate[4])), na.rm = TRUE)
  # timeStep <- 24*365*(secondDate[1] - firstDate[1]) + 24*30*(secondDate[2] - firstDate[2]) + 24*(secondDate[3] - firstDate[3]) + (secondDate[4] - firstDate[4])
  ndr <- array(data = NA, dim = c(dim(slp$Data)[1], jeqlat, ieqlon))
  for (i in 1:dim(slp$Data)[2]){
    ndr[1:(dim(slp$Data)[1]-1),i,] <- (sin(60*2*pi/360)/sin(slp$xyCoords$y[i]*2*pi/360))*(slp$Data[2:dim(slp$Data)[1],i,]-slp$Data[1:(dim(slp$Data)[1]-1),i,])/timeStep
  }
  cyclonCenter <- vector("list", dim(slp$Data)[1])
  icount <- array(data = 0, dim = c(dim(slp$Data)[1], 1))
  ################################################################################################################
  for (t in 1:dim(vo$Data)[1]){
    auxDate <- strsplit(vo$Dates$start[t],' ')
    yymmdd <- strsplit(auxDate[[1]][1],'-')
    hhmmss <- strsplit(auxDate[[1]][2],':')
    ymdhms <- c(as.double(yymmdd[[1]][1]),as.double(yymmdd[[1]][2]),as.double(yymmdd[[1]][3]),as.double(hhmmss[[1]][1]),as.double(hhmmss[[1]][2]),as.double(hhmmss[[1]][3]))
    nstorms <- 0
    ## Initialize variables
    # tdiff <- matrix(data = 0, nrow = 60, ncol = 1)
    # iswitch <- matrix(data = 0, nrow = jeqlat, ncol = ieqlon)
    # ibad <- matrix(data = 0, nrow = jeqlat, ncol = ieqlon)
    # ncenter <- matrix(data = 0, nrow = jeqlat, ncol = ieqlon)
    k <- 0
    auxInd <- NULL    
    for (j in 4:(ieqlon-3)){
      for (i in 4:(jeqlat-3)){
        ibad1 <- matrix(as.numeric(vo$Data[t,i,j] - vo$Data[t,(i-3):(i+3),(j-3):(j+3)] > vo.diff.threshold & vo$Data[t,i,j] > vo.threshold), nrow = 7, ncol =7) ## & rv$Data[t,i,j]>1e-5
        ibad2 <- matrix(as.numeric(lap[t,i,j] - lap[t,(i-3):(i+3),(j-3):(j+3)] > lap.diff.threshold), nrow = 7, ncol =7) ## & rv$Data[t,i,j]>1e-5
        ibad3 <- matrix(as.numeric(slp$Data[t,i,j] - slp$Data[t,(i-3):(i+3),(j-3):(j+3)] < slp.diff.threshold), nrow = 7, ncol =7) ## & rv$Data[t,i,j]>1e-5
        if (sum(ibad1[3:5,3:5], na.rm=TRUE) == 8 | sum(ibad2[3:5,3:5], na.rm=TRUE) == 8 | sum(ibad3[3:5,3:5], na.rm=TRUE) == 8){
          auxInd <- c(auxInd,(j-1)*jeqlat+i)
        }
      }
    }
    if (any(ndr[t,,] > ndr.threshold | vo$Data[t,,] > vo.threshold, na.rm = TRUE) | !is.null(auxInd)){
      nstorms <- 1
      auxInd <- c(auxInd,which(lap[t,,] > lap.diff.threshold & vo$Data[t,,] > vo.threshold))
      auxInd <- unique(auxInd)
      normDimensions <- c(3,2)
      cycloGenesis <- 0
      a1 <- ndr[t,,]
      if (any(a1[auxInd] > ndr.threshold, na.rm = TRUE)){
        print(paste0("There is an Explosive Cyclogenesis at ",slp$Dates$start[t]))
        auxInd <- c(auxInd,which(lap[t,,] > lap.diff.threshold & vo$Data[t,,] > vo.threshold & ndr[t,,] > ndr.threshold))
        auxInd <- unique(auxInd)
        cycloGenesis <- 1
        normDimensions <- c(1,3,2)
      }
      auxCyclos <- a1[auxInd]
      a1 <- slp$Data[t,,]
      auxCyclos <- cbind(auxCyclos, a1[auxInd])
      a1 <- vo$Data[t,,]
      auxCyclos <- cbind(auxCyclos, a1[auxInd])
      a1 <- lap[t,,]
      auxCyclos <- cbind(auxCyclos, a1[auxInd])
      auxCyclos <- cbind(auxCyclos, long[auxInd], lat[auxInd])
      if (is.null(wss)){
        a1 <- matrix(data = NA, nrow = length(slp$xyCoords$y), ncol = length(slp$xyCoords$x))
      }else{
        a1 <- wss$Data[t,,]
      }
      auxCyclos <- cbind(auxCyclos, a1[auxInd])
      if (length(auxInd) > 1){
        auxCyclos1 <- auxCyclos
        auxCyclos1[which(auxCyclos[,1]<0.5),1] <- 0
        auxCyclos1[,2] <- slp.threshold-auxCyclos[,2]
        auxCyclos1 <- scale(auxCyclos1, center = TRUE, scale = TRUE)
        auxCyclos1 <- sqrt(apply(auxCyclos1[,normDimensions]^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
        if (criteria == "global") {## Global criterium:
          maxGlobal <- which.max(auxCyclos1)
        } else if (criteria == "max.vo") {## Maximum Vorticity:
          maxGlobal <- which.max(auxCyclos[,3])
        } else if (criteria == "max.ndr") {## Maximum NDR:
          maxGlobal <- which.max(auxCyclos[,1])
        } else if (criteria == "min.slp") {## Minimum SLP:
          maxGlobal <- which.min(auxCyclos[,2])
        }
        varOutput1 <- c(auxCyclos[maxGlobal,], cycloGenesis)
        out <- 0
        while (out == 0 && any(sqrt((auxCyclos[,5]-auxCyclos[maxGlobal,5])^2+(auxCyclos[,6]-auxCyclos[maxGlobal,6])^2)>seek.radius, na.rm = TRUE)){
          nstorms <- 1+nstorms
          if (length(which(sqrt((auxCyclos[,5]-auxCyclos[maxGlobal,5])^2+(auxCyclos[,6]-auxCyclos[maxGlobal,6])^2)>seek.radius))>1){
            auxCyclos <- auxCyclos[which(sqrt((auxCyclos[,5]-auxCyclos[maxGlobal,5])^2+(auxCyclos[,6]-auxCyclos[maxGlobal,6])^2)>seek.radius),]
            auxCyclos1 <- auxCyclos
            auxCyclos1[,2] <- slp.threshold-auxCyclos[,2]
            auxCyclos1 <- scale(auxCyclos1, center = TRUE, scale = TRUE)
            auxCyclos1 <- sqrt(apply(auxCyclos1[,normDimensions]^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
            if (criteria == "global") {## Global criterium:
              maxGlobal <- which.max(auxCyclos1)
            } else if (criteria == "max.vo") {## Maximum Vorticity:
              maxGlobal <- which.max(auxCyclos[,3])
            } else if (criteria == "max.ndr") {## Maximum NDR:
              maxGlobal <- which.max(auxCyclos[,1])
            } else if (criteria == "min.slp") {## Minimum SLP:
              maxGlobal <- which.min(auxCyclos[,2])
            }
            if (!is.na(auxCyclos[maxGlobal,1]) && auxCyclos[maxGlobal,1] > ndr.threshold){
              cycloGenesis <- 1
            } else {
              cycloGenesis <- 0
            }
            varOutput1 <- rbind(varOutput1,c(auxCyclos[maxGlobal,], cycloGenesis))
          }else{
            auxCyclos <- auxCyclos[which(sqrt((auxCyclos[,5]-auxCyclos[maxGlobal,5])^2+(auxCyclos[,6]-auxCyclos[maxGlobal,6])^2)>seek.radius),]
            if (!is.na(auxCyclos[1]) && auxCyclos[1] > ndr.threshold){
              cycloGenesis <- 1
            } else {
              cycloGenesis <- 0
            }
            varOutput1 <- rbind(varOutput1,c(auxCyclos, cycloGenesis))
            maxGlobal <- 1
            out <- 1
          }
        }
      }else{
        varOutput1 <- c(auxCyclos, cycloGenesis)
        nstorms <- 1
      }
      cyclonCenter[[t]] <- varOutput1
      icount[t] <- nstorms
    }
  }
  return(cyclonCenter)
}

#######################################################################
#########                DIFERENCIA DE LONGITUD               #########
#######################################################################
difLong <- function (xlon1,xlon2){
  if (is.na(xlon1) | is.na(xlon2)){
    diflong <- NA
  }else if (abs(xlon1-xlon2) <= 180.0){
    diflong <- abs(xlon2 - xlon1)
  }else if(xlon2 < xlon1){
    diflong <- abs(xlon2 - xlon1 + 360.0)
  }else{
    diflong <- abs(xlon2 - xlon1 - 360.0)
  }
  return(diflong)
}
######################################################################
#########                     LAPLACIANO                     #########
######################################################################
laplacian <- function (obj){
  r <- 6378.1 # Earth radius
  ieqlon <- length(obj$xyCoords$x)
  jeqlat <- length(obj$xyCoords$y)
  ## Offset and scale of the data
  lat <- matrix(data = obj$xyCoords$y, nrow = jeqlat, ncol = ieqlon)
  long <- t(matrix(data = obj$xyCoords$x, nrow = ieqlon, ncol = jeqlat))
  ## Auxiliary matrix to obtain the laplacian:
  meters <- array(data = NA, dim=c(jeqlat,ieqlon,4))
  for (j in 2:(dim(lat)[2]-1)) {
    for (i in 2:(dim(lat)[1]-1)) {
      ## d1 is the lat/lon distance between (i+1,j) and (i,j)
      d1lon <- difLong(long[i+1,j],long[i,j])*pi*r*cos(lat[i,j]*(2*pi)/360)/180
      d1lat <- abs((lat[i+1,j]-lat[i,j])*r*(2*pi)/360)
      ## d2 is the lat/lon distance between (i-1,j) and (i,j)
      d2lon <- difLong(long[i,j],long[i-1,j])*pi*r*cos(lat[i,j]*(2*pi)/360)/180
      d2lat <- abs((lat[i,j]-lat[i-1,j])*r*(2*pi)/360)
      ## d3 is the lat/lon distance between (i,j-1) and (i,j)
      d3lon <- difLong(long[i,j-1],long[i,j])*pi*r*cos(lat[i,j]*(2*pi)/360)/180
      d3lat <- abs((lat[i,j-1]-lat[i,j])*r*(2*pi)/360)
      ## d4 is the lat/lon distance between (i,j+1) and (i,j)
      d4lon <- difLong(long[i,j],long[i,j+1])*pi*r*cos(lat[i,j]*(2*pi)/360)/180
      d4lat <- abs((lat[i,j]-lat[i,j+1])*r*(2*pi)/360)
      ## meters(i,j, ) is distance between the points
      meters[i,j,1] <- sqrt(d1lat^2 + d1lon^2)
      meters[i,j,2] <- sqrt(d2lat^2 + d2lon^2)
      meters[i,j,3] <- sqrt(d3lat^2 + d3lon^2)
      meters[i,j,4] <- sqrt(d4lat^2 + d4lon^2)
    }
  }
  d1 <- meters[,,1]*(meters[,,1]+meters[,,2])
  d2 <- meters[,,2]*(meters[,,1]+meters[,,2])
  d3 <- meters[,,3]*(meters[,,3]+meters[,,4])
  d4 <- meters[,,4]*(meters[,,3]+meters[,,4])
  lap <- obj
  lap$Data[,c(1:2,jeqlat-1,jeqlat),] <- NA
  lap$Data[,,c(1:2,ieqlon-1,ieqlon)] <- NA
  zdiff1 <- obj$Data[,5:(jeqlat-2),4:(ieqlon-3)] - obj$Data[,4:(jeqlat-3),4:(ieqlon-3)]
  zdiff2 <- obj$Data[,3:(jeqlat-4),4:(ieqlon-3)] - obj$Data[,4:(jeqlat-3),4:(ieqlon-3)]
  zdiff3 <- obj$Data[,4:(jeqlat-3),3:(ieqlon-4)] - obj$Data[,4:(jeqlat-3),4:(ieqlon-3)]
  zdiff4 <- obj$Data[,4:(jeqlat-3),5:(ieqlon-2)] - obj$Data[,4:(jeqlat-3),4:(ieqlon-3)]
  for (j in 4:(jeqlat-3)){
    lap$Data[,j,4:(ieqlon-3)] <- 2*((zdiff1[,j-3,]/d1[j,4:(ieqlon-3)])+(zdiff2[,j-3,]/d2[j,4:(ieqlon-3)])+(zdiff3[,j-3,]/d3[j,4:(ieqlon-3)])+(zdiff4[,j-3,]/d4[j,4:(ieqlon-3)]))
  }
  lap$Data[,c(1:3,jeqlat-2,jeqlat-1,jeqlat),] <- NA
  lap$Data[,,c(1:3,ieqlon-2,ieqlon-1,ieqlon)] <- NA
  return(lap)
}
