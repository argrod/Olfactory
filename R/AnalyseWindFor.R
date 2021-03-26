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

sel <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/gribSelected.csv", sep = ",", header = T)
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
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/WindCorr.png", , dpi = 300, height = 6,
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
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/SpeedCorr.png", , dpi = 300, height = 6,
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
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindCol.png", , dpi = 300, height = 6,
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
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindBW.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

cor.test(sel$ESpd, sel$WSpd)


    scale_fill_viridis_d(option="magma")
    scale_fill_manual(values=c("#ffffbf","#ffffbf","#2c7bb6","#d7191c"))


###################################################################################################
####################################### HYPLIT ####################################################
###################################################################################################

library(devtools)
devtools::install_github("rich-iannone/SplitR")
library(SplitR)
install.packages("here",,"https://mac.R-project.org")
library(here)
library(dplyr)

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

trajectory <- 
  hysplit_trajectory(
    lat = 39,
    lon = 143,
    height = 10,
    duration = 4,
    days = "2018-09-01",
    daily_hours = c(14),
    direction = "backward",
    met_type = "gdas1",
    extended_met = FALSE
    # met_dir = here::here("met"),
    # exec_dir = here::here("out")
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
    met_type = "reanalysis"
    # met_dir = here::here("met"),
    # exec_dir = here::here("out")
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
    met_type = "reanalysis"
    # met_dir = here::here("met"),
    # exec_dir = here::here("out")
  ) %>%
    run_model()

dispersion_tbl <- dispersion_model %>% get_output_tbl()
dispersion_model %>% dispersion_plot()
