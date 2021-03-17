library(Gmisc)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)
library(rgdal)
library(adehabitatHR)
library(ggsn)
library(ggspatial)
library(sp)
library(plyr)
library(dplyr)
library(mapdata)
library(rerddap)
library(data.table)
library(ggplot2)
library(viridis)
library(MASS)
library(diagram)
library(ggthemes)
library(extrafont)

#################################################################################
######################## BRING IN THE FORAGING ESTIMATES ########################
#################################################################################

if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/Temp2018.RData")
    # outloc <- "/Volumes/GoogleDrive/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("F:/UTokyoDrive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("F:/UTokyoDrive/PhD/Data/Temp2018.RData")
    # outloc <- "F:/UTokyoDrive/PhD/Manuscripts/BehaviourIdentification/Figures/"
}
# D18 <- bind_rows(Dat)
###############################################################################
######################## BRING IN THE WIND ESTIMATIONS ########################
###############################################################################

if(Sys.info()['sysname'] == "Darwin"){
    windLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
} else {
    windLoc <- 'F:/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/MinDat/'
}
windFiles <- dir(windLoc)

for(b in 1:length(windFiles)){
    if(b == 1){
        WindDat <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = F)
        WindDat$ID <- sub("*.csv", "", windFiles[b])
    } else {
        toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = F)
        toAdd$ID <- sub("*.csv", "", windFiles[b])
        WindDat <- rbind(WindDat, toAdd)
    }
}
colnames(WindDat) <- c("DT","lat","lon","head","X","Y","ID")
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y-%m-%d %H:%M:%OS", tz = "")

WindDat$WHead <- atan2(WindDat$Y, WindDat$X)
WindDat$WSpeed <- sqrt(WindDat$X^2 + WindDat$Y^2)

ggplot(WindDat[WindDat$ID == "1_S2",], aes(x = lon, y = lat)) +
    geom_point() +
    geom_spoke(data = WindDat[WindDat$ID == "1_S2",], arrow = arrow(length = unit(WindDat$WSpeed[WindDat$ID == "1_S2"]/max(
        WindDat$WSpeed[WindDat$ID == "1_S2"])*0.15, 'inches')),
        aes(x = lon, y = lat, angle = WindDat$WHead[WindDat$ID == "1_S2"], col = WindDat$WSpeed[WindDat$ID == "1_S2"], radius = scales::rescale(WindDat$WSpeed[WindDat$ID == "1_S2"], c(.1, .5)))) +
    scale_colour_gradient("Wind speed", low = "yellow", high = "red")

#################################################################################################################
############################  FINDING FORAGING WITH WIND CALCULATED BEFORE (30 MINS) ############################
#################################################################################################################

# deal with on tag-by-tag basis
# wind data (WindDat)
# foraging data (Dat[[tag#]])
tg <- 1
DatSel <- Dat[[tg]]
DatSel$forageLong[DatSel$tkb == 1 | DatSel$dv == 1] = 1
DatSel$forageLong[is.na(DatSel$forageLong)] = 0
WindSel <- WindDat[WindDat$ID == DatSel$tagID[tg],]
#find shared latlons of tag and wind data
dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
splits <- which(dt > 70)
splits <- c(0, splits, nrow(WindSel))
DatSel$windCal <- 0
for(b in 1:(length(splits) - 1)){
    DatSel$windCal[DatSel$DT > (WindSel$DT[splits[b] + 1] - 30) & DatSel$DT < (WindSel$DT[splits[b + 1] -1] + 30)] = 1
}
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$forageLong[b:nrow(DatSel)] != 1)){
        DatSel$tFromFor[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromFor[b] <- difftime(DatSel$DT[b + min(which(DatSel$forageLong[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# find points where there a wind calculations and <1 hour to foraging
DatSel$windNear <- DatSel$windCal == 1 & (DatSel$tFromFor < (60*60))

# # find turning points (0 xings of UTME and UTMN)
# pos <- which(diff(diff(DatSel$UTME) > 0) != 0) + 2
# DatSel$turn[pos] <- 1
# # timelag until turn
# b=1
# while(b<nrow(DatSel)){
#     if(all(DatSel$turn[b:nrow(DatSel)] != 1)){
#         DatSel$tFromTurn[b:nrow(DatSel)] <- NA
#         b <- nrow(DatSel)
#     } else {
#         DatSel$tFromTurn[b] <- difftime(DatSel$DT[b + min(which(DatSel$turn[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
#         b = b+1
#     }
# }
# where nearby wind calculations (5 mins) and nearby turns (<5 mins) overlap
# DatSel$turnNear <- DatSel$tFromTurn < (5*60) & DatSel$windNear == 1
# alternatively, find turning points of > 30 degrees
DatSel$dx = c(NA,diff(DatSel$UTME))
DatSel$dy = c(NA,diff(DatSel$UTMN))
DatSel$angle = atan2(DatSel$dy, DatSel$dx)
DatSel$lgTurn = DatSel$angle*(180/pi) > 30
# time to large turn
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$lgTurn[b:nrow(DatSel)] != 1)){
        DatSel$tFromLgTurn[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromLgTurn[b] <- difftime(DatSel$DT[b + min(which(DatSel$lgTurn[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# where nearby wind calculations (5 mins) and nearby large turns (>30 deg, <5 mins) overlap
DatSel$lgTurnNear <- DatSel$tFromLgTurn < (5*60) & DatSel$windNear == 1

# create a mean vector for each selected period
st = which(diff(DatSel$turnNear) == 1) + 1
ed = which(diff(DatSel$turnNear) == -1) + 1
if(length(st) < length(ed)){
    st <- c(1, st)
}
if(length(ed) < length(st)){
    ed <- c(ed, nrow(DatSel))
}
for(b in 1:length(st)){
    DatSel$mnWindDx <- mean
}

######################################################################################################################
########################### TAKE NEAR FORAGING AND FIND AVE HEADING 30 MINS PRIOR ####################################
######################################################################################################################

for(tg in 1:length(Dat)){
    DatSel = Dat[[tg]]
    DatSel$forageLong[DatSel$tkb == 1 | DatSel$dv == 1] = 1
    DatSel$forageLong[is.na(DatSel$forageLong)] = 0
    WindSel = WindDat[WindDat$ID == unique(DatSel$tagID),]
    dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
    splits <- which(dt > 70)
    splits <- c(0, splits, nrow(WindSel))
    WindSel$nrFor = NA
    for(b in 1:nrow(WindSel)){
        WindSel$nrFor[b] = any(DatSel$forageLong[DatSel$DT >= WindSel$DT[b] & DatSel$DT < (WindSel$DT[b] + 3600)] == 1)
    }
    WindSel$BSpd = NA
    for(b in 1:(length(splits) - 1)){
        cord.dec <- SpatialPoints(cbind(WindSel$lon[(splits[b]+1):splits[b+1]], WindSel$lat[(splits[b]+1):splits[b+1]]),
            proj4string=CRS("+proj=longlat"))
        cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

        X <- cord.UTM$coords.x1
        Y <- cord.UTM$coords.x1
        vg_x_obs <- X[2:length(X)] - X[1:(length(X) - 1)]
        vg_y_obs <- Y[2:length(Y)] - Y[1:(length(Y) - 1)]
        times = WindSel$DT[(splits[b]+1):splits[b+1]]
        dt = as.numeric(difftime(times[2:length(times)],times[1:(length(times)-1)], units = 'secs'))
        g_speed <- sqrt(vg_y_obs^2 + vg_x_obs^2)/dt
        WindSel$BSpd[(splits[b]+1):splits[b+1]] = c(NA, g_speed)
    }
    if(tg == 1){
        WFor <- WindSel[WindSel$nrFor == 1,]
    } else {
        toAdd <- WindSel[WindSel$nrFor == 1,]
        WFor <- rbind(WFor, toAdd)
    }
}
# plot the difference between the two
ggplot(WFor) +
    geom_point(aes(x = (WHead - head)*(180/pi), y = WSpeed)) +
    coord_polar() + scale_x_continuous(limits = c(0,180))

Wplot <- ggplot(WFor) +
    geom_point(aes(x = (WHead), y = WSpeed), pch = 2) +
    coord_polar() + scale_x_continuous(limits = c(-3,3))

Bplot <- ggplot(WFor) +
    geom_point(aes(x = head, y = BSpd), pch = 4)+
    coord_polar() + scale_x_continuous(limits = c(-3,3))

ggarrange(Wplot, Bplot)


dt <- difftime(WindDat$DT[2:nrow(WindDat)],WindDat$DT[1:(nrow(WindDat) - 1)], "units", "secs")
splits <- which(dt > 0 & dt > 70)
splits <- c(0, splits, nrow(WindDat))
WindDat$BSpd = NA
for(b in 1:(length(splits) - 1)){
    cord.dec <- SpatialPoints(cbind(WindDat$lon[(splits[b]+1):splits[b+1]], WindDat$lat[(splits[b]+1):splits[b+1]]),
        proj4string=CRS("+proj=longlat"))
    cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

    X <- cord.UTM$coords.x1
    Y <- cord.UTM$coords.x1
    vg_x_obs <- X[2:length(X)] - X[1:(length(X) - 1)]
    vg_y_obs <- Y[2:length(Y)] - Y[1:(length(Y) - 1)]
    times = WindDat$DT[(splits[b]+1):splits[b+1]]
    dt = as.numeric(difftime(times[2:length(times)],times[1:(length(times)-1)], units = 'secs'))
    g_speed <- sqrt(vg_y_obs^2 + vg_x_obs^2)/dt
    WindDat$BSpd[(splits[b]+1):splits[b+1]] = c(NA, g_speed)
}

# plot bird speed against wind heading
ggplot(WindDat, aes(x = WHead, y = BSpd)) +
    geom_point() + coord_polar()





# BWcomp <- data.frame(aWh=0,aWs=0,aBh=0,aBs=0,nrFor=0)
# for(b = 1:(length(splits)-1)){
#     BWcomp$aWh[b] = atan2(mean(WindSel$X[(splits[b] + 1):splits[b+1]]), mean(WindSel$Y[(splits[b] + 1):splits[b+1]]))
#     BWcomp$aWs[b] = mean(sqrt(WindSel$X[(splits[b] + 1):splits[b+1]]^2 + WindSel$Y[(splits[b] + 1):splits[b+1]]^2))
#     BWcomp$aBh[b] = 
# }



# dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
# splits <- which(dt > 70)
# splits <- c(0, splits, nrow(WindSel))
# DatSel$windCal <- 0
# for(b in 1:(length(splits) - 1)){
#     DatSel$windCal[DatSel$DT > (WindSel$DT[splits[b] + 1] - 30) & DatSel$DT < (WindSel$DT[splits[b + 1] -1] + 30)] = 1
# }
# b=1
# while(b<nrow(DatSel)){
#     if(all(DatSel$forageLong[b:nrow(DatSel)] != 1)){
#         DatSel$tFromFor[b:nrow(DatSel)] <- NA
#         b <- nrow(DatSel)
#     } else {
#         DatSel$tFromFor[b] <- difftime(DatSel$DT[b + min(which(DatSel$forageLong[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
#         b = b+1
#     }
# }
# windDiff <- as.numeric(difftime(tail(WindSel$DT, -1), head(WindSel$DT, -1), units = 'secs'))
# splits <- which(windDiff > (60*2))
# WindSel$aveX <- NA
# WindSel$aveY <- NA
# for(b in 1:length(splits)){
#     if(b == 1){
#         WindSel$aveX[1:splits[b]] <- sum((WindSel$X[1:splits[b]]), na.omit = T)
#         WindSel$aveY[1:splits[b]] <- sum((WindSel$Y[1:splits[b]]), na.omit = T)
#     } else if(b == length(splits)){
#         WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveX[splits[b]:nrow(WindSel)] <- sum((WindSel$X[splits[b]:nrow(WindSel)]), na.omit = T)
#         WindSel$aveY[splits[b]:nrow(WindSel)] <- sum((WindSel$Y[splits[b]:nrow(WindSel)]), na.omit = T)
#     } else {
#         WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
#     }
# }

# WindSel$headAve <- atan2(WindSel$aveY, WindSel$aveX)

# plot(WindSel$headAve)

# ggplot(data = DatSel, aes(x = Lon, y = Lat)) +
# geom_path() +
# geom_point(data = DatSel[DatSel$windCal == 1, ], aes(x = Lon, y = Lat), pch = 1) + 
#     geom_spoke(data = WindSel, arrow = arrow(length = unit(WindSel$FlSpeed/max(WindSel$FlSpeed)*0.15, 'inches')),
#         aes(x = lon, y = lat, angle = headAve, col = FlSpeed, radius = scales::rescale(FlSpeed, c(.1, .5)))) +
#     scale_colour_gradient("Wind speed", low = "yellow", high = "red")
#     #geom_point(data = DatSel[DatSel$windCal == T, ], aes(x = Lon, y = Lat), pch = 1)



# ggplot(data = DatSel) +
# geom_path(aes(x = Lon, y = Lat)) + 
# geom_point(data = DatSel[DatSel$turnNear == 1, ], aes(x = Lon, y = Lat))
# geom_point(aes(x = Lon, y = Lat, col = tFromFor)) +
# geom_point(data = DatSel[DatSel$Forage == 1,], aes(x = Lon, y = Lat), col = 'red')
