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
WindDat$FlSpeed <- sqrt(WindDat$X^2 + WindDat$Y^2)

ggplot(WindDat[WindDat$ID == "1_S2",], aes(x = lon, y = lat)) +
    geom_point() +
    geom_spoke(data = WindDat[WindDat$ID == "1_S2",], arrow = arrow(length = unit(WindDat$FlSpeed[WindDat$ID == "1_S2"]/max(
        WindDat$FlSpeed[WindDat$ID == "1_S2"])*0.15, 'inches')),
        aes(x = lon, y = lat, angle = WindDat$WHead[WindDat$ID == "1_S2"], col = WindDat$FlSpeed[WindDat$ID == "1_S2"], radius = scales::rescale(WindDat$FlSpeed[WindDat$ID == "1_S2"], c(.1, .5)))) +
    scale_colour_gradient("Wind speed", low = "yellow", high = "red")

#################################################################################################################
############################  FINDING FORAGING WITH WIND CALCULATED BEFORE (30 MINS) ############################
#################################################################################################################

# deal with on tag-by-tag basis
# wind data (WindDat)
# foraging data (Dat[[tag#]])
tg <- 1
DatSel <- Dat[[tg]]
WindSel <- WindDat[WindDat$ID == DatSel$tagID[tg],]
#find shared latlons of tag and wind data
DatSel$windCal <- DatSel$Lat %in% WindSel$Lat & DatSel$Lon %in% WindSel$Lon
# find time lag since foraging
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$Forage[b:nrow(DatSel)] != 1)){
        DatSel$tFromFor[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromFor[b] <- difftime(DatSel$DT[b + min(which(DatSel$Forage[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# find points where there a wind calculations and <5 mins to foraging
DatSel$windNear <- DatSel$windCal == 1 & DatSel$tFromFor < (5*60*60)

# find turning points (0 xings of UTME and UTMN)
pos <- which(diff(diff(DatSel$UTME) > 0) != 0) + 2
DatSel$turn[pos] <- 1
# timelag until turn
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$turn[b:nrow(DatSel)] != 1)){
        DatSel$tFromTurn[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromTurn[b] <- difftime(DatSel$DT[b + min(which(DatSel$turn[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# where nearby wind calculations (5 mins) and nearby turns (<5 mins) overlap
DatSel$turnNear <- DatSel$tFromTurn < (5*60*60) & DatSel$windNear == 1
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
DatSel$lgTurnNear <- DatSel$tFromLgTurn < (5*60*60) & DatSel$windNear == 1

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

windDiff <- as.numeric(difftime(tail(WindSel$DT, -1), head(WindSel$DT, -1), units = 'secs'))
splits <- which(windDiff > (60*2))
WindSel$aveX <- NA
WindSel$aveY <- NA
for(b in 1:length(splits)){
    if(b == 1){
        WindSel$aveX[1:splits[b]] <- sum((WindSel$X[1:splits[b]]), na.omit = T)
        WindSel$aveY[1:splits[b]] <- sum((WindSel$Y[1:splits[b]]), na.omit = T)
    } else if(b == length(splits)){
        WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
        WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
        WindSel$aveX[splits[b]:nrow(WindSel)] <- sum((WindSel$X[splits[b]:nrow(WindSel)]), na.omit = T)
        WindSel$aveY[splits[b]:nrow(WindSel)] <- sum((WindSel$Y[splits[b]:nrow(WindSel)]), na.omit = T)
    } else {
        WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
        WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
    }
}

WindSel$headAve <- atan2(WindSel$aveY, WindSel$aveX)

plot(WindSel$headAve)

ggplot(data = DatSel, aes(x = Lon, y = Lat)) +
geom_path() +
geom_point(data = DatSel[DatSel$windCal == 1, ], aes(x = Lon, y = Lat), pch = 1) + 
    geom_spoke(data = WindSel, arrow = arrow(length = unit(WindSel$FlSpeed/max(WindSel$FlSpeed)*0.15, 'inches')),
        aes(x = lon, y = lat, angle = headAve, col = FlSpeed, radius = scales::rescale(FlSpeed, c(.1, .5)))) +
    scale_colour_gradient("Wind speed", low = "yellow", high = "red")
    #geom_point(data = DatSel[DatSel$windCal == T, ], aes(x = Lon, y = Lat), pch = 1)



ggplot(data = DatSel) +
geom_path(aes(x = Lon, y = Lat)) + 
geom_point(data = DatSel[DatSel$turnNear == 1, ], aes(x = Lon, y = Lat))
geom_point(aes(x = Lon, y = Lat, col = tFromFor)) +
geom_point(data = DatSel[DatSel$Forage == 1,], aes(x = Lon, y = Lat), col = 'red')
