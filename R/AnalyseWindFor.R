install.packages("Gmisc")
install.packages("sf")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("lubridate")
install.packages("rgdal")
install.packages("adehabitatHR")
install.packages("ggsn")
install.packages("ggspatial")
install.packages("sp")
install.packages("plyr")
install.packages("dplyr")
install.packages("mapdata")
install.packages("rerddap")
install.packages("data.table")
install.packages("ggplot2")
install.packages("viridis")
install.packages("MASS")
install.packages("diagram")
install.packages("ggthemes")
install.packages("extrafont")
install.packages("bpnreg")

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
library(rgeos)
library(bpnreg)

#################################################################################
######################## BRING IN THE FORAGING ESTIMATES ########################
#################################################################################

if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/DatEth2018.RData")
    load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/DatEth2019.RData")
    outloc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("E:/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("E:/My Drive/PhD/Data/DatEth2018.RData")
    load("E:/My Drive/PhD/Data/DatEth2019.RData")
    outloc <- "E:/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
}
D18 <- bind_rows(Dat)
D19 <- bind_rows(Dat19)
allD <- data.frame(DT=c(D18$DT, D19$DT),
    lat = c(D18$Lat, D19$Lat),
    lon = c(D18$Lon, D19$Lon),
    forage = c(D18$Forage,D19$Forage),
    tagID = c(D18$tagID, D19$tagID),
    Day = c(D18$Day, D19$Day),
    Sex = c(D18$Sex, D19$Sex),
    distTrav = c(D18$recalDist, D19$recalDist),
    spTrav = c(D18$spTrav, D19$spTrav),
    recalSp = c(D18$recalSp, D19$recalSp),
    distFk = c(D18$distFromFk, D19$distFromFk),
    tripN = c(D18$tripN, D19$tripN),
    tripL = c(D18$tripL, D19$tripL),
    tkb = c(D18$tkb, D19$tkb),
    dv = c(D18$dv, D19$dv),
    UTME = c(D18$UTME, D19$UTME),
    UTMN = c(D18$UTMN, D19$UTMN))
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
allD$yrID <- paste(format(allD$DT,"%Y"),sub('\\_S.*','',allD$tagID),sep="_")
###############################################################################
######################## BRING IN THE WIND ESTIMATIONS ########################
###############################################################################

if(Sys.info()['sysname'] == "Darwin"){
    windLoc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
} else {
    windLoc <- 'E:/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/'
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

ggplot(WindDat[WindDat$ID == "1_S2",], aes(x = DT, y = head-WHead)) +
    geom_point()

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
WFor$dir <- NA
WFor$dir[(WFor$head - WFor$WHead) > (-pi/4) & (WFor$head - WFor$WHead) < (pi/4)] = "Same direction"
WFor$dir[(WFor$head - WFor$WHead) < (-pi/4) & (WFor$head - WFor$WHead) > (-3*pi/4)] = "L-Side"
WFor$dir[(WFor$head - WFor$WHead) < (-3*(pi/4)) | (WFor$head - WFor$WHead) > (3*(pi/4))] = "Toward"
WFor$dir[(WFor$head - WFor$WHead) < (3*pi/4) & (WFor$head - WFor$WHead) > (pi/4)] = "R-Side"
# plot the difference between the two
ggplot(WFor) +
    geom_point(aes(x = (head - WHead)*(180/pi), y = WSpeed)) +
    coord_polar() + scale_x_continuous(limits = c(0,180))

ggplot(WFor) +
    geom_point(aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed)) +
    coord_polar() + scale_x_continuous(limits = c(0,360))

ggplot() +
    geom_point(data = WFor[(WFor$head - WFor$WHead) < pi/2 & (WFor$head - WFor$WHead) >(-pi/2),],
        mapping=aes(x = (head - WHead), y = WSpeed), colour = "green") +
        geom_point(data = WFor[(WFor$head - WFor$WHead) < (-pi/2) | (WFor$head - WFor$WHead) >pi/2,],
        mapping=aes(x = (head - WHead), y = WSpeed), colour = "red") +
    coord_polar() + scale_x_continuous(limits = c(-pi,pi))



ggplot() +
    geom_point(data = WFor[((WFor$head+3) - (WFor$WHead + 3))*(180/pi) < 90,],
        mapping=aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed), colour = "green") +
        geom_point(data = WFor[((WFor$head+3) - (WFor$WHead + 3))*(180/pi) > 90,],
        mapping=aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed), colour = "red") +
    coord_polar() + scale_x_continuous(limits = c(0, 360))

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed, fill = dir), pch = 21) +
    coord_polar() + scale_x_continuous(limits = c(-pi, pi))

Wplot <- ggplot(WFor) +
    geom_point(aes(x = (WHead), y = WSpeed), pch = 2) +
    coord_polar() + scale_x_continuous(limits = c(-pi,pi))

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

sel <- read.delim("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/gribSelected.csv", sep = ",", header = T)
sel$EHead <- atan2(sel$X,sel$Y)
sel$WHead <- atan2(sel$U,sel$V)
ggplot() +
    geom_point(aes(y = EHead, x = WHead), data = sel, pch = 21, fill = "deepskyblue1") +
    geom_line(aes(x = (-pi:pi), y = (-pi:pi))) +
    theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated headings") +
    scale_x_continuous("JMA headings")
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Conferences/2021SeabirdSymposium/WindCorr.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()
library(circular)
res <- cor.circular(sel$WHead, sel$EHead, test = T)
res
ggplot() +
    geom_point(aes(x = WSpd, y = ESpd), data= sel, pch = 21, fill = "orangered2") +
    geom_line(aes(x = (0:12), y =(0:12))) +
    theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)") +
    scale_x_continuous("JMA wind speed (m/s)")
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Conferences/2021SeabirdSymposium/SpeedCorr.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed, fill = dir), pch = 21, size = 2) +
    coord_polar() + scale_x_continuous("Wind heading relative to bird (rad)",limits = c(-pi, pi)) +
    scale_fill_discrete("Wind heading") +
    theme_bw() + 
    theme(panel.border = element_blank(), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)")
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindCol.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed), pch = 21, fill = "black", size = 1) +
    coord_polar() + scale_x_continuous("Wind heading relative to bird (rad)",limits = c(-pi, pi)) +
    theme_bw() + 
    theme(panel.border = element_blank(), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)")
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindBW.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

cor.test(sel$ESpd, sel$WSpd)


    scale_fill_viridis_d(option="magma")
    scale_fill_manual(values=c("#ffffbf","#ffffbf","#2c7bb6","#d7191c"))


###################################################################################################
####################################### HYPLIT ####################################################
###################################################################################################

library(devtools)
devtools::install_github("rich-iannone/splitr")
library(splitr)
install.packages("here",,"https://mac.R-project.org")
library(here)
library(dplyr)
library(ggplot2)
library(plyr)

setwd("~/Documents/SplitR_wd")

trajectory <- 
  hysplit_trajectory(
    lat = 42.83752,
    lon = -80.30364,
    height = 10,
    duration = 24,
    daily_hours = c(0, 6, 12, 18),
    direction = "forward",
    met_type = "gdas1",
    extended_met = TRUE)

if(Sys.info()['sysname'] == "Darwin"){
    exec_loc <- "/Users/aran/hysplit/"
} else {
    exec_loc <- "C:/hysplit/"
}

trajectory <- 
  hysplit_trajectory(
    lat = 39,
    lon = 143,
    height = 10,
    duration = 4,
    days = "2018-09-01",
    daily_hours = c(14),
    direction = "backward",
    met_type = "reanalysis",
    extended_met = FALSE,
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  )

trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = 39,
    lon = 143,
    height = 10,
    duration = 6,
    days = "2018-09-01",
    daily_hours = c(14),
    direction = "backward",
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
  run_model()

trajectory_tbl <- trajectory_model %>% get_output_tbl()

trajectory_tbl
trajectory_tbl %>% trajectory_plot()
trajectory_model %>% trajectory_plot()


dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = 41.0, lon = 143.0, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm("2018-09-01 14:31"),
    release_end = lubridate::ymd_hm("2018-09-01 14:31") + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm("2018-09-01 14:31"),
    end_time = lubridate::ymd_hm("2018-09-01 14:31") + lubridate::hours(6),
    direction = "backward", 
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()

dispersion_tbl <- dispersion_model %>% get_output_tbl()
dispersion_model %>% dispersion_plot()


find_hull <- function(df) df[chull(df$lon, df$lat), ]
hulls <- ddply(dispersion_tbl, "hour", find_hull)

plot <- ggplot(data = dispersion_tbl, aes(x = lon, y = lat, colour=as.factor(hour), fill = as.factor(hour), group = as.factor(hour))) +
# geom_point() + 
geom_polygon(data = hulls, alpha = 0.5) +
labs(x = "Lon", y = "Lat")
plot

#################################################################
#################### ANALYSE 2019 YONE DATA #####################
#################################################################

# bring in wind and foraging data
if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
} else {
    fileloc <- "E:/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
    forLoc <- "E:/My Drive/PhD/Data/2019Shearwater/TxtDat/AxyTrek/AlgorithmOutput/PredictedForage/"
}

files <- dir(fileloc)
forFiles <- dir(forLoc, pattern = "*ForageGPS.txt")

tags = unique(sub("_S.*", "",forFiles))

# read in the foraging data
# forD <- vector(mode="list", length=length(tags))
# for (t in 1:length(tags)) {
#     dayFiles = forFiles[grepl(paste("^",tags[t],"_S.*",sep=""),forFiles)]
#     ds <- data.frame(DT=character(),lat=numeric(),lon=numeric(),Forage=integer())
#     outpt <- ds
#     for (d in 1:length(dayFiles)){
#         outpt <- rbind(outpt,read.delim(paste(forLoc,dayFiles[d],sep=""),sep=",",header=T))
#     }
#     outpt$DT <- as.POSIXct(outpt$DT,"%d-%b-%Y %H:%M:%S",tz="")
#     outpt$ID <- tags[t]
#     w.dec <- SpatialPoints(cbind(outpt$Lon,outpt$Lat),proj4string = CRS("+proj=longlat"))
#     UTMdat <- spTransform(w.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
#     outpt$UTME <- UTMdat$coords.x1
#     outpt$UTMN <- UTMdat$coords.x2
#     forD[[t]] <- outpt
# }
# read in wind data
wind <- vector(mode="list",length=length(tags))
for(b in 1:length(files)){
    wind[[b]] <- read.delim(paste(fileloc,files[b],sep=""),sep=",",header=T)
    wind[[b]]$time <- as.POSIXct(wind[[b]]$time,"%Y/%m/%d,%H:%M:%S",tz="")
    w.dec <- SpatialPoints(cbind(wind[[b]]$lon,wind[[b]]$lat),proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(w.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    wind[[b]]$UTME <- UTMdat$coords.x1
    wind[[b]]$UTMN <- UTMdat$coords.x2
    wind[[b]] <- wind[[b]][!is.nan(wind[[b]]$Resnorm),]
    wind[[b]]$timeTo <- NA
    wind[[b]]$distTo <- NA
    for (g in 1:nrow(wind[[b]])){
        if (any(forD[[b]]$Forage[forD[[b]]$DT > wind[[b]]$time[g]] == 1)) {
            point <- which(abs(forD[[b]]$DT - wind[[b]]$time[g]) == min(abs(forD[[b]]$DT - wind[[b]]$time[g])))
            forInd <- min(which(forD[[b]]$Forage[point:nrow(forD[[b]])] == 1)) + point - 1 # select the minimum (i.e. the next one)
            wind[[b]]$timeTo[g] <- as.numeric(difftime(forD[[b]]$DT[forInd],wind[[b]]$time[g],units="secs"))
            wind[[b]]$distTo[g] <- sqrt( (forD[[b]]$UTMN[forInd] - forD[[b]]$UTMN[point])^2 + 
                (forD[[b]]$UTME[forInd] - forD[[b]]$UTME[point])^2 ) * 10^(-3)
        }
    }
    wind[[b]]$rwh <- wind[[b]]$aveDir - wind[[b]]$wDir
    wind[[b]]$ID <- tags[b]
}
dat$time <- as.POSIXct(dat$time, format = "%Y/%m/%d,%H:%M:%S", tz = "")
dat$rwh <- dat$aveDir - dat$wDir
dat$rwh[dat$aveDir == 0] <- NA
dat$U <- dat$wSp*cos(dat$wDir)
dat$V <- dat$wSp*sin(dat$wDir)

selectDat <- dat[!is.na(dat$rwh),]

colnames(dat)
unique()
ggplot() + 
    geom_point(data = dat[dat$ID == unique(dat$ID)[1],],
        aes(x=time,y=rwh))

ggplot(dat[dat$distTo < 10 & dat$forage != 1,]) + geom_point(aes(y = rwh, x = distTo)) + scale_x_continuous(limits=c(-180,180))
hist(dat$distTo)

ggplot(windAll) + 
    geom_point(aes(x = rwh*(pi/180), y = distTo)) + coord_polar() + scale_x_continuous(limits = c(-pi,pi))

FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
FkOshi.dec <- SpatialPoints(cbind(FkOshi$Long,FkOshi$Lat),proj4string=CRS('+proj=longlat'))
FkOshi.UTM <- spTransform(FkOshi.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
FkUTME <- FkOshi.UTM$coords.x1
FkUTMN <- FkOshi.UTM$coords.x2
for (b in 1:length(wind)){
    wind[[b]]$distFromFk <- sqrt((wind[[b]]$UTMN - FkUTMN)^2 + (wind[[b]]$UTME - FkUTME)^2) * 10^(-3)
    wind[[b]]$tSinceStart <- difftime(wind[[b]]$time, wind[[b]]$time[1], units="secs")
}

# go through each individuals tracks
ggplot(wind[[6]]) +
    geom_point(aes(x=rwh,y=timeTo/60)) + coord_polar() #+
    #geom_line(aes(x=rwh,y=distFromFk))
    # + facet_grid(~ ID) 
colnames(windAll)

ggplot(windAll[windAll$distTo < 30,]) +
    geom_point(aes(x = rwh*(pi/180), y = distTo, fill = ID),pch=21, alpha=.3) + coord_polar(start=pi)


dloadLoc  = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
if(any(gribFls == "Z__C_RJTD_20190824150000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"))

###################################################################
#################### 2019 YONE WIND FORAGING ######################
###################################################################


if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("E:/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("E:/My Drive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "E:/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
}
D19 <- bind_rows(Dat19)
allD <- data.frame(DT=D19$DT,
    lat = D19$Lat,
    lon = D19$Lon,
    tagID = D19$tagID,
    Day = D19$Day,
    Sex = D19$Sex,
    distTrav = D19$recalDist,
    spTrav = D19$spTrav,
    recalSp = D19$recalSp,
    distFk = D19$distFromFk,
    tripN = D19$tripN,
    tripL = D19$tripL,
    tkb = D19$tkb,
    dv = D19$dv,
    UTME = D19$UTME,
    UTMN =  D19$UTMN)
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
# allD$forBeh <- NA
# allD$forBeh[allD$dv == 1] <- "Dive"
# allD$forBeh[allD$tkb == 1] <- "Surf"
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")

# CALCULATE RELATIVE WIND CONDITIONS FROM ESTIMATES WITH TIME/DISTANCE TO FORAGING
if(Sys.info()['sysname'] == "Darwin"){
    windLoc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
} else {
    windLoc = "G:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
}
windFiles <- dir(windLoc)
for(b in 1:length(windFiles)){
  if(b == 1){
    WindDat <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
    colnames(WindDat) <- c("DT","Lat","Lon","Head","wDir","wSp","Resnorm")
    WindDat$ID <- sub("*WindYone.txt", "", windFiles[b])
    Wind.dec <- SpatialPoints(cbind(WindDat$Lon,WindDat$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    WindDat$UTME <- UTMdat$coords.x1
    WindDat$UTMN <- UTMdat$coords.x2
  } else {
    toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
    colnames(toAdd) <- c("DT","Lat","Lon","Head","wDir","wSp","Resnorm")
    toAdd$ID <- sub("*WindYone.txt", "", windFiles[b])
    Add.dec <- SpatialPoints(cbind(toAdd$Lon,toAdd$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Add.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    toAdd$UTME <- UTMdat$coords.x1
    toAdd$UTMN <- UTMdat$coords.x2 
    WindDat <- rbind(WindDat, toAdd)
  }
}
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y/%m/%d,%H:%M:%OS",tz="UTC")
# remove missing wind values
WindDat <- WindDat[!is.nan(WindDat$Resnorm),]
# remove points within 5km of FkOshi
FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
FkOshi.dec <- SpatialPoints(cbind(FkOshi$Long,FkOshi$Lat),proj4string=CRS('+proj=longlat'))
FkOshi.UTM <- spTransform(FkOshi.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
FkUTME <- FkOshi.UTM$coords.x1
FkUTMN <- FkOshi.UTM$coords.x2

WindDat$timeTo <- NA
WindDat$distTo <- NA
allD$forage[is.na(allD$forage)] <- 0
tags <- unique(WindDat$ID)
for(b in 1:nrow(WindDat)){ # for each row in WindDat
  if(any(allD$forage[allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y") & allD$DT > WindDat$DT[b]] == 1)){ # if there are any foraging points after timepoint
    point <- which(allD$lat == WindDat$Lat[b] & allD$lon == WindDat$Lon[b] & allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y")) # find the timepoint of the next foraging points
    forInd <- min(which(allD$forage[point:max(which(allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y")))] == 1)) + point - 1 # select the minimum (i.e. the next one)
    WindDat$timeTo[b] <- as.numeric(difftime(allD$DT[point+forInd],WindDat$DT[b], units="secs"))
    WindDat$distTo[b] <- sqrt((allD$UTMN[forInd] - allD$UTMN[point])^2 + (allD$UTME[forInd] - allD$UTME[point])^2)*10^-3
  }
}

#####################################################################################################
############################### BAYESIAN CIRCULAR MIXED EFFECTS MODEL ###############################
#####################################################################################################

fit.Motor <- bpnr(pred.I = Phaserad ~ 1 + Cond, data = Motor, its = 10000, burn = 100, n.lag = 3, seed = 101)


traceplot(fit.Motor,parameter="beta1")

#####################################################################################################
####################################### STRAIGHTNESS OF TRACK #######################################
#####################################################################################################

library(geosphere)
HdistTrav <- lapply(unique(allD$yrID), function(x) c(NA,distHaversine(cbind(allD$lon[allD$yrID==x][1:(length(allD$lon[allD$yrID==x])-1)],allD$lat[allD$yrID==x][1:(length(allD$lat[allD$yrID==x])-1)]),
    cbind(allD$lon[allD$yrID==x][2:(length(allD$lon[allD$yrID==x]))],allD$lat[allD$yrID==x][2:(length(allD$lat[allD$yrID==x]))]))))

length(HdistTrav[[1]])

nrow(allD)

distHaversine(cbind(allD$lon[allD$yrID==x][1:(length(allD$lon[allD$yrID==x])-1)],allD$lat[allD$yrID==x][1:(length(allD$lat[allD$yrID==x])-1)]),
    cbind(allD$lon[allD$yrID==x][2:(length(allD$lon[allD$yrID==x]))],allD$lat[allD$yrID==x][2:(length(allD$lat[allD$yrID==x]))]))

 SpatialPoints(cbind(WindSel$lon[(splits[b]+1):splits[b+1]], WindSel$lat[(splits[b]+1):splits[b+1]]),
            proj4string=CRS("+proj=longlat"))
