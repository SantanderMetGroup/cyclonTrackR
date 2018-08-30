#     getCyclonTrack.R Algorithm to obtain the cyclon tracks from the cyclon centers
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

#' @title Cyclon Track Detection
#' @description Implementation of an algorithm to obtain the cyclon track from the cyclon centers
#'
#' @param cyclonCenter List of cyclon centers detected for each day.
#' @param seek.radius Maximum distance in degrees to look for a new cyclon center to be included in the same trajectory. The default value is 10.
#' @param ndr.threshold Minimum value of the Normalized Deepeningun Ratio to be considered a cyclon center as explosive. The default value is 1.
#' @param vo.threshold Minimum value of vorticity to be considered a candidate to cyclon center. The default value is 1e-4.
#' @param cut Maximum length, in time steps, to be considered for the track. The default value is 4.
#' @param cyclon.date Argument to obtain all the tracks occurred in an specific date. The default value is \code{NULL}.
#' @param list.date Argument containing the dates corresponding to the cyclon centers. The default value is \code{NULL}.
#' @param criteria Criteria considered to select the cyclon center between the possible candidates. Options implemented are:
#' \code{"global"} (default): maximum of the Euclidean norm of a vector containing the vorticity, the NDR and the slp normalized;
#' \code{"max.vo"}: point with the maximum vorticity; \code{"max.ndr"}: point with the maximum NDR; \code{"min.slp"}: point with the minimum slp.
#'     
#' @seealso \code{\link{getCyclonCenters}} for the algorithm obtaining the cyclon centers. 
#' @return A list with the cyclon tracks detected. Each element contains the dates of the track and the matrix with the corresponding
#' NDR, slp, vorticity, laplacian of the slp, longitude, latitude, windspeed and explosive character of the cyclon center.
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
#' list.date<-slp$Dates$start
#' cyclonTrack <- getCyclonTrack(varOutputCenters,
#'                               seek.radius = seek.radius,
#'                               ndr.threshold = ndr.threshold,
#'                               vo.threshold = vo.threshold,
#'                               max.length = 20,
#'                               cyclon.date = NULL,
#'                               list.date = list.date,
#'                               criteria = criteria, cut = 20)
#' 
#' ## Plotting the results
#' data("wrld")
#' po <- SpatialPoints(cbind(cyclonTrack[[1]][[2]][,5], cyclonTrack[[1]][[2]][,6]))
#' dat <- as.data.frame(cyclonTrack[[1]][[2]][,2])
#' colnames(dat) <- "y"
#' kl <- SpatialPointsDataFrame(po, data = dat)
#' spplot(kl, zcol = "y", sp.layout = list(wrld, first = F), colorkey = TRUE, 
#'               xlim = lonLim, ylim = latLim, xlab = slp$Dates$start[1])
#'
#' ## Considering an specific date
#' cyclonTrack <- getCyclonTrack(varOutputCenters,
#'                               seek.radius = seek.radius,
#'                               ndr.threshold = ndr.threshold,
#'                               vo.threshold = vo.threshold,
#'                               max.length = 20,
#'                               cyclon.date = "2010-02-26 12:00:00 GMT",
#'                               list.date = list.date,
#'                               criteria = criteria, cut = 20)
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


getCyclonTrack <- function(cyclonCenter,
                           seek.radius = 10,
                           ndr.threshold = 1,
                           vo.threshold = 1e-4,
                           max.length = 4,
                           cyclon.date = NULL,
                           list.date = NULL,
                           criteria = "global",
                           cut=4){
  ## Definition of the parameters:
  r <- 6378.1 # Earth radius
  g <- 9.8 # Earth gravity
  omega <- 7.2921e-5 # Earth angular velocity
  pdiff <- 1e-6 ## pressure threshold (mb) for identifying cyclones.
  maxdist <- 800 ## maximum search distance (km) for tracking storms between time steps. 
  pchange <- 2*1e-5 ## limit on allowable absolute pressure tendency (mb) between time steps. 
  #######################################################################
  if (is.null(cyclon.date)){
    cyclonTrack <- NULL
    Ncyclon <- 0
    Nassigned <- c(0,0)
    for (i in 1:length(cyclonCenter)){
      Nassigned[1] <- Nassigned[1] + dim(cyclonCenter[[i]])[1]
    }
    ######### Time Step #########
    t <- 1
    while (Nassigned[2] < Nassigned[1] & t < length(cyclonCenter)){
      ######### Storm in time step#########
      n1 <- 0
      while (any(!is.na(cyclonCenter[[t]][,1])) & n1 < dim(cyclonCenter[[t]])[1]){
        auxTT <- cyclonCenter[[t]][,c(3,1,2)]
        auxTT[which(auxTT[,2]<0.5),2] <- 0
        auxTT[,3] <- 2000 - auxTT[,3]
        no.na<-which(!is.na(cyclonCenter[[t]][,1]))
        if(length(no.na)>1){
          auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
          if (criteria == "global") {## Global criterium:
            auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
            maxGlobal <- which.max(auxTT)
          } else if (criteria == "max.ndr") {## Maximum NDR:
            maxGlobal <- which.max(auxTT[,2])
          } else if (criteria == "max.vo") {## Maximum Vorticity:
            maxGlobal <- which.max(auxTT[,1])
          } else if (criteria == "min.slp") {## Minimum SLP:
            maxGlobal <- which.min(auxTT[,3]-2000)
          }
        }else if (length(no.na)>0) {
          maxGlobal<-no.na
        }
        ### Add date ###
        dates<-NULL
        dates<-rbind(dates, list.date[t])
        ### Add center ###
        xyv <- NULL
        xyv <- rbind(xyv,cyclonCenter[[t]][maxGlobal,])
        ### Remove center ###
        cyclonCenter[[t]][maxGlobal,] <- NA
        Nassigned[2] <- Nassigned[2] + 1
        cut_forward<-0
        ######### Track #########
        for (tt in (t+1):min(c(t+max.length-1,length(cyclonCenter)))){
          if (cut_forward<cut){
            indClose <- which(sqrt((cyclonCenter[[tt]][,6] - xyv[dim(xyv)[1],6])^2+(cyclonCenter[[tt]][,5] - xyv[dim(xyv)[1],5])^2) < seek.radius & 
                              cyclonCenter[[tt]][,1] > ndr.threshold & cyclonCenter[[tt]][,3] > vo.threshold)
          
            if (length(indClose) > 1){
              auxTT <- cyclonCenter[[tt]][indClose,c(3,1,2)]
              auxTT[which(auxTT[,2]<0.5),2] <- 0
              auxTT[,3] <- 2000 - auxTT[,3]
              auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
              if (criteria == "global") {## Global criterium:
                auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
                maxGlobal <- which.max(auxTT)
              } else if (criteria == "max.ndr") {## Maximum NDR:
                maxGlobal <- which.max(auxTT[,2])
              } else if (criteria == "max.vo") {## Maximum Vorticity:
                maxGlobal <- which.max(auxTT[,1])
              } else if (criteria == "min.slp") {## Minimum SLP:
                maxGlobal <- which.min(auxTT[,3]-2000)
              }
              ### Add date ###
              dates<-rbind(dates, list.date[tt])
              ### Add center ###
              xyv <- rbind(xyv,cyclonCenter[[tt]][indClose[maxGlobal],])
              ### Remove center ###
              cyclonCenter[[tt]][indClose[maxGlobal],] <- NA
              Nassigned[2] <- Nassigned[2] + 1
              
            } else if (length(indClose) > 0){
              ### Add date ###
              dates<-rbind(dates, list.date[tt])
              ### Add center ###
              xyv <- rbind(xyv,cyclonCenter[[tt]][indClose,])
              ### Remove center ###
              cyclonCenter[[tt]][indClose,] <- NA
              Nassigned[2] <- Nassigned[2] + 1
            } else {
              cut_forward<-cut_forward+1
            }
          }
          
        }
        n1 <- n1 + 1
        Ncyclon <- Ncyclon + 1
        cyclonTrack[[Ncyclon]] <- list(dates,xyv)
      }
      t <- t+1
    }
  }else{
    t <- which(as.Date(list.date) == cyclon.date)
    if (length(t)>1){
      txyv <- NULL
      for (t1 in 1:length(t)){
        auxTT <- cyclonCenter[[t[t1]]][,c(3,1,2)]
        auxTT[which(auxTT[,2]<0.5),2] <- 0
        auxTT[,3] <- 2000 - auxTT[,3]
        auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
        if (criteria == "global") {## Global criterium:
          auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
          maxGlobal <- which.max(round(auxTT,6))
        } else if (criteria == "max.ndr") {## Maximum NDR:
          maxGlobal <- which.max(auxTT[,2])
        } else if (criteria == "max.vo") {## Maximum Vorticity:
          maxGlobal <- which.max(auxTT[,1])
        } else if (criteria == "min.slp") {## Minimum SLP:
          maxGlobal <- which.min(cyclonCenter[[t[t1]]][,2])
        }
        txyv <- rbind(txyv,cyclonCenter[[t[t1]]][maxGlobal,])
      }
      auxTT <- txyv[,c(3,1,2)]
      auxTT[which(auxTT[,2]<0.5),2] <- 0
      auxTT[,3] <- 2000 - auxTT[,3]
      auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
      if (criteria == "global") {## Global criterium:
        auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
        maxGlobal <- which.max(round(auxTT,6))
      } else if (criteria == "max.ndr") {## Maximum NDR:
        maxGlobal <- which.max(auxTT[,2])
      } else if (criteria == "max.vo") {## Maximum Vorticity:
        maxGlobal <- which.max(auxTT[,1])
      } else if (criteria == "min.slp") {## Minimum SLP:
        maxGlobal <- which.min(auxTT[,3]-2000)
      }
      t <- t[maxGlobal]
    }
    if (is.matrix(cyclonCenter[[t]])){
      Ncyclon <- dim(cyclonCenter[[t]])[1]
      cyclonTrack <- vector("list", Ncyclon)
    } else if (is.value(cyclonCenter[[t]])){
      Ncyclon <- 1
      cyclonTrack <- vector("list", Ncyclon)
    } else {
      Ncyclon <- 0
      cyclonTrack <- NULL
    }
    i <- 0
    while (i < Ncyclon){
      xyv <- NULL
      auxTT <- cyclonCenter[[t]][,c(3,1,2)]
      auxTT[which(auxTT[,2]<0.5),2] <- 0
      auxTT[,3] <- 2000 - auxTT[,3]
      if (i< (Ncyclon-1)){
        auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
        if (criteria == "global") {## Global criterium:
          auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
          maxGlobal <- which.max(round(auxTT,6))
        } else if (criteria == "max.ndr") {## Maximum NDR:
          maxGlobal <- which.max(auxTT[,2])
        } else if (criteria == "max.vo") {## Maximum Vorticity:
          maxGlobal <- which.max(auxTT[,1])
        } else if (criteria == "min.slp") {## Minimum SLP:
          maxGlobal <- which.min(auxTT[,3]-2000)
        }
      }else{
        maxGlobal<-which(!is.na(cyclonCenter[[t]][,1]))
      }
      ### Add date ###
      dates<-NULL
      dates<-rbind(dates, list.date[t])
      ### Add center ###
      xyv <- rbind(xyv,cyclonCenter[[t]][maxGlobal,])
      ### Remove center ###
      cyclonCenter[[t]][maxGlobal,] <- NA
      nn <- 1
      
      cut_backward<-0
      cut_forward<-0
      for (tt in (t+1):min(c(t+max.length,length(cyclonCenter)))){
        if(cut_forward<cut){
          indClose <- which(sqrt((cyclonCenter[[tt]][,6] - xyv[nn,6])^2+(cyclonCenter[[tt]][,5] - xyv[nn,5])^2) < seek.radius & 
                            cyclonCenter[[tt]][,1] > ndr.threshold & 
                            cyclonCenter[[tt]][,3] > vo.threshold)
          if (length(indClose) > 1){
            auxTT <- cyclonCenter[[tt]][indClose,c(3,1,2)]
            auxTT[which(auxTT[,2]<0.5),2] <- 0
            auxTT[,3] <- 2000 - auxTT[,3]
            auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
            if (criteria == "global") {## Global criterium:
              auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
              maxGlobal <- which.max(round(auxTT,6))
            } else if (criteria == "max.ndr") {## Maximum NDR:
              maxGlobal <- which.max(auxTT[,2])
            } else if (criteria == "max.vo") {## Maximum Vorticity:
              maxGlobal <- which.max(auxTT[,1])
            } else if (criteria == "min.slp") {## Minimum SLP:
              maxGlobal <- which.min(auxTT[,3]-2000)
            }
            ### Add date ###
            dates<-rbind(dates, list.date[tt])
            ### Add center ###
            xyv <- rbind(xyv,cyclonCenter[[tt]][indClose[maxGlobal],])
            ### Remove center ###
            cyclonCenter[[tt]][indClose[maxGlobal],] <- NA
            
          } else if (length(indClose) > 0){
            ### Add date ###
            dates<-rbind(dates, list.date[tt])
            ### Add center ###
            xyv <- rbind(xyv,cyclonCenter[[tt]][indClose,])
            ### Remove center ###
            cyclonCenter[[tt]][indClose,] <- NA
            
          }else{
            cut_forward<-cut_forward+1
          }
        }
        
        if (t-(tt-t) > 0 & cut_backward<cut){
          indClose <- which(sqrt((cyclonCenter[[t-(tt-t)]][,6] - xyv[1,6])^2+(cyclonCenter[[t-(tt-t)]][,5] - xyv[1,5])^2) < seek.radius &
                              cyclonCenter[[t-(tt-t)]][,1] > ndr.threshold & cyclonCenter[[t-(tt-t)]][,3] > vo.threshold)
          if (length(indClose) > 1){
            auxTT <- cyclonCenter[[t-(tt-t)]][indClose,c(3,1,2)]
            auxTT[which(auxTT[,2]<0.5),2] <- 0
            auxTT[,3] <- 2000 - auxTT[,3]
            auxTT <- scale(auxTT, center = TRUE, scale = TRUE)
            if (criteria == "global") {## Global criterium:
              auxTT <- sqrt(apply(auxTT^2, FUN = "sum", na.rm = TRUE, MARGIN = 1))
              maxGlobal <- which.max(round(auxTT,6))
            } else if (criteria == "max.ndr") {## Maximum NDR:
              maxGlobal <- which.max(auxTT[,2])
            } else if (criteria == "max.vo") {## Maximum Vorticity:
              maxGlobal <- which.max(auxTT[,1])
            } else if (criteria == "min.slp") {## Minimum SLP:
              maxGlobal <- which.min(auxTT[,3]-2000)
            }
            ### Add date ###
            dates<-rbind(list.date[t-(tt-t)], dates)
            ### Add center ###
            xyv <- rbind(cyclonCenter[[t-(tt-t)]][indClose[maxGlobal],],xyv)
            ### Remove center ###
            cyclonCenter[[t-(tt-t)]][indClose[maxGlobal],] <- NA
            
          } else if (length(indClose) > 0){
            ### Add date ###
            dates<-rbind(list.date[t-(tt-t)], dates)
            ### Add center ###
            xyv <- rbind(cyclonCenter[[t-(tt-t)]][indClose,],xyv)
            ### Remove center ###
            cyclonCenter[[t-(tt-t)]][indClose,] <- NA
          }else{
            cut_backward<-cut_backward+1
          }
        }
        if (is.matrix(xyv)){
          nn <- dim(xyv)[1]
        }
      }
      i <- i + 1
      cyclonTrack[[i]] <- list(dates,xyv)
    }
  }
  return(cyclonTrack)
}
