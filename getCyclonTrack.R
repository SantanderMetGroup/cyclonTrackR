#######################################################################
##                    CYCLONE TRACKING DETECTION                     ##
#######################################################################
## based on: http://www.europeanwindstorms.org/method/track/
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
