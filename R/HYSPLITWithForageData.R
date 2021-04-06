library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
options(timeout = 200)
if(Sys.info()['sysname'] == "Darwin"){
    exec_loc <- "/Users/aran/hysplit/"
} else {
    exec_loc <- "C:/hysplit/"
}
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/Temp2018.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("F:/UTokyoDrive/PhD/Data/Temp2018.RData")
}
D18 <- bind_rows(Dat)
D19 <- bind_rows(Dat19)
allD <- data.frame(DT=c(D18$DT, D19$DT),
    lat = c(D18$Lat, D19$Lat),
    lon = c(D18$Lon, D19$Lon),
    tagID = c(D18$tagID, D19$tagID),
    Day = c(D18$Day, D19$Day),
    Sex = c(D18$Sex, D19$Sex),
    distTrav = c(D18$recalDist, D19$recalDist),
    recalSp = c(D18$recalSp, D19$recalSp),
    distFk = c(D18$distFromFk, D19$distFromFk),
    tripL = c(D18$tripL, D19$tripL),
    tkb = c(D18$tkb, D19$tkb),
    dv = c(D18$dv, D19$dv))
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")

find_hull <- function(df) df[chull(df$lon, df$lat), ]
TrackDisp <- function(DT, lat, lon, hrs){
  time <- format(DT, "%Y-%m-%d %H:%M")
  dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = lat, lon = lon, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm(time) - lubridate::hours(hrs),
    release_end = lubridate::ymd_hm(time) - lubridate::hours(hrs) + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm(time) - lubridate::hours(hrs),
    end_time = lubridate::ymd_hm(time),
    direction = "forward", 
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()
  dispersion_tbl <- dispersion_model %>% get_output_tbl()
  hulls <- ddply(dispersion_tbl, "hour", find_hull)
  outpt <- list("DispModel" = dispersion_model,"partDisp" = dispersion_tbl, "partPoly" = hulls)
}
TrackTraj <- function(DT, lat, lon, hrs){
  time <- format(DT, "%Y-%m-%d")
  hr <- format(DT, "%H")
  trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = lat,
    lon = lon,
    height = 50,
    duration = hrs,
    days = time,
    daily_hours = hr,
    direction = "backward",
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()
  return(trajectory_model %>% get_output_tbl())
}

# testing ------------------------------------------------------------------------------------------------------------
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
    lat = 42.42768,
    lon = 143.8522,
    height = 10,
    duration = 6,
    days = "2018-08-31",
    daily_hours = c(19),
    direction = "backward",
    met_type = "gdas1",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
  run_model()

trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = 42,
    lon = 143,
    height = 50,
    duration = 6,
    days = "2018-09-01",
    daily_hours = c(12,14),
    direction = "backward",
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()


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
    direction = "forward", 
    met_type = "gdas1",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()

dispersion_tbl <- dispersion_model %>% get_output_tbl()
ab <- dispersion_model %>% dispersion_plot()


hulls <- ddply(dispersion_tbl, "hour", find_hull)

plot <- ggplot(data = dispersion_tbl, aes(x = lon, y = lat, colour=as.factor(hour), fill = as.factor(hour), group = as.factor(hour))) +
# geom_point() + 
  geom_polygon(data = hulls, alpha = 0.5) +
  labs(x = "Lon", y = "Lat") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA))
plot

ggplot() +
  geom_point(dispersion_tbl, mapping = aes(x = lon, y = lat, colour = as.factor(hour), stroke = 0), alpha = .5) +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  coord_sf(xlim = c(141, 143.5), ylim = c(40, 41.5)) +
  labs(x = "Lon", y = "Lat") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA))

ggplot() +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  coord_sf(xlim = c(141, 144), ylim = c(40, 43)) +
  geom_polygon(data = hulls, aes(x = lon, y = lat, colour=as.factor(hour), fill = as.factor(hour), group = as.factor(hour)), alpha = 0.5) +
  labs(x = "Lon", y = "Lat")







for1 <- sample(which(allD$forage == 1),1)
test <- TrackDisp(allD$DT[for1], allD$lat[for1], allD$lon[for1],6)

dispplot <- ggplot() +
  theme(panel.border = element_rect(colour = 'black', fill = NA)) +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  coord_sf(xlim = c(139, 145), ylim = c(39, 45)) +
  # geom_point(data = test$partDisp, aes(x = lon, y = lat, colour=as.factor(hour), group = as.factor(hour))) + 
  geom_polygon(data = test$partPoly, alpha = 0.5, aes(x = lon, y = lat, fill = as.factor(hour))) +
  scale_fill_discrete(name = "Hours before") +
  scale_colour_manual(name = "Hours before", labels = c("1","2","3","4","5","6"), values = (1:6)) +
  labs(x = "Lon", y = "Lat") + theme_bw() + theme(panel.grid = element_blank()) +
  geom_path(data = allD[allD$DT > (allD$DT[for1] - lubridate::hours(6)) & allD$DT <= allD$DT[for1] & allD$tagID == allD$tagID[for1],], aes(x = lon, y = lat))
dispplot

test$DispModel %>% dispersion_plot()


png("/Volumes/GoogleDrive/My Drive/PhD/Notes/HYSPLIT/DispPlotZoom.png",
  width=600, height = 600)
dispplot
dev.off()

test <- TrackTraj(allD$DT[for1], allD$lat[for1], allD$lon[for1],6)
test

trajplot <- ggplot() +
  theme(panel.border = element_rect(colour = 'black', fill = NA)) +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  coord_sf(xlim = c(144, 146), ylim = c(41, 43.5)) +
  geom_point(data = test, aes(x = lon, y = lat, colour = as.factor(hour_along))) +
  scale_fill_discrete(name = "Hours before") +
  labs(x = "Lon", y = "Lat") + theme_bw() + theme(panel.grid = element_blank()) +
  geom_path(data = allD[allD$DT > (allD$DT[for1] - lubridate::hours(6)) & allD$DT <= allD$DT[for1] & allD$tagID == allD$tagID[for1],], aes(x = lon, y = lat)) +
  scale_colour_manual(name = "Hours before", labels = c("-6","-5","-4","-3","-2","-1","0"), values = as.factor(-6:0))
trajplot
png("/Volumes/GoogleDrive/My Drive/PhD/Notes/HYSPLIT/TrajPlot.png",
  width=600, height = 600)
trajplot
dev.off()
######################################################################################################

# trajectory analysis ------------------------------------------------------------------------------------------

######################################################################################################
###################################### TRAJECTORY EVERY 6 HOURS ######################################
######################################################################################################

# go through each tag (each list value within Dat or Dat19)
AllD <- c(Dat,Dat19)
hrlyDat <- vector(mode = "list", length = length(AllD))
for(h in 1:length(AllD)){
  sel <- AllD[[h]]
  # remove non-flight values
  sel <- sel[which(sel$spTrav > 15),]
  # colnames(sel)
  sel$tdiff <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel) - 1)], units = 'secs'))
  # decide outgoing/incoming
  sel$appSpd <- c(NA, diff(sel$distFromFk)) # approach speed to Fk Island
  sel$rtChg <- NA
  for(c in 1:nrow(sel)){
    sel$rtChg[c] <- mean(sel$appSpd[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]/as.numeric(sel$tdiff[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]))
  }
  tsel <- seq(sel$DT[1], sel$DT[length(sel$DT)], by = 3600)
  dtFull <- as.numeric(diff(sel$DT))
  # find data that line up to the new timepoints
  select <- NA
  for(b in 1:length(tsel)){
      choose <- which(sel$DT >= (tsel[b] - (median(dtFull))) & sel$DT <= (tsel[b] + (median(dtFull))))
      if(length(choose) != 0){
          if(length(choose) > 1){
              select[b] <- choose[which.min(abs(tsel[b] - sel$DT [choose]))]
          } else {
              select[b] <- choose
          }
      }
  }
  select <- select[!is.na(select)]
  # select hourly subsetted data
  hrSel <- sel[select,]
  trun <- vector(mode = "list", length = nrow(hrSel))
  for(d in 1:length(select)){
    tryCatch({
      trun[[d]] <- TrackTraj(hrSel$DT[d], hrSel$Lat[d], hrSel$Lon[d], 12)
    }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
  }
  aveHead <- NA
  aveTrHead <- NA
  aveTrWSpd <- NA
  for(c in 1:(length(trun) - 1)){
    cord.dec <- SpatialPoints(cbind(sel$Lon[select[c]:select[c+1]], sel$Lat[select[c]:select[c+1]]),proj4string=CRS("+proj=longlat"))
    cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    utmDiff <- cbind(diff(cord.UTM$coords.x1), diff(cord.UTM$coords.x2))
    aveHead[c] <- atan2(mean(utmDiff[,1]),mean(utmDiff[,2]))

    tryCatch({
      trCord.dec <- SpatialPoints(cbind(trun[[c]]$lon, trun[[c]]$lat),
                  proj4string=CRS("+proj=longlat"))
      trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
      trutmDiff <- cbind(diff(trCord.utm$coords.x1), diff(trCord.utm$coords.x2))
      aveTrHead[c] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
      aveTrWSpd[c] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)
    }, error = function(e){c(NA)})
  }
  hrSel$aveHead <- c(aveHead,NA)
  hrSel$aveTrHead <- c(aveTrHead,NA)
  hrSel$aveTrWSpd <- c(aveTrWSpd,NA)
  hrlyDat[[h]] <- hrSel
}
allHrly <- bind_rows(hrlyDat)
colnames(allHrly)

testDf <- data.frame(avHd = allHrly$aveHead[!is.na(allHrly$aveTrHead)]*(180/pi) + 180,
  avTrHd = allHrly$aveTrHead[!is.na(allHrly$aveTrHead)]*(180/pi) + 180, outRet = allHrly$rtChg[!is.na(allHrly$aveTrHead)] < 0)

ggplot(testDf[!is.na(testDf$avTrHd),], aes(x = avHd, y = avTrHd, colour = outRet)) +
  geom_point() + scale_x_continuous("Average bird heading (1 hr)") + 
  scale_y_continuous("Average air trajectory heading (12 hr)") +
  scale_colour_manual("Direction from nest", values = c("#41b6c4","#ec7014"), breaks = c(FALSE, TRUE),
    labels = c("Return", "Leave"))

ggplot(testDf[!is.na(testDf$avTrHd)], aes(x = outRet, y  = avHd - avTrHd)) + 
  geom_point()
plot(testDf$avHd[testDf$outRet == F])
points(testDf$avTrHd[testDf$outRet == F], pch=2)

sel <- allHrly[!is.na(allHrly$aveTrHead),]
sel$relHead <- sel$aveHead - sel$aveTrHead
sel$relHead[sel$relHead < -pi] = sel$relHead[sel$relHead < -pi] + 2*pi
sel$relHead[sel$relHead > pi] = sel$relHead[sel$relHead > pi] - 2*pi
sel$RelWind <- NA
sel$RelWind[sel$relHead*(180/pi) < 45 & sel$relHead*(180/pi) > -45] <- "Tail"
sel$RelWind[sel$relHead*(180/pi) > 45 & sel$relHead*(180/pi) < 135] <- "RSide"
sel$Rel Wind[sel$relHead*(180/pi) < -45 & sel$relHead*(180/pi) > -135] <- "LSide"
sel$RelWind[sel$relHead*(180/pi) < -135 | sel$relHead*(180/pi) > 135] <- "Head"
sel$Yr <- as.numeric(format(sel$DT, "%Y"))

tags = unique(sel$tagID)
ggplot(sel[sel$tagID == tags[7] & sel$Yr == 2018 & sel$rtChg < 0,]) + 
  # geom_path(aes(x=Lon, y = Lat)) +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  coord_sf(xlim = c(139, 148), ylim = c(39, 45)) +
  geom_point(aes(x = Lon, y = Lat, colour = RelWind)) +
  ggtitle(paste(max()))
  
  scale_colour_manual(palette = "Temperature")
  scale_colour_gradientn(colours = terrain.colors(10))
  (rtChg > 0)))
  (aveHead*(180/pi) - aveTrHead*(180/pi))))

plot(sel$DT[sel$tagID == sel$tagID[1]], sel$rtChg[sel$tagID == sel$tagID[1]], type = 'l')

testDf$relHead <- testDf$avHd - testDf$avTrHd
testDf$relHead[testDf$relHead < -pi] = testDf$relHead[testDf$relHead < -pi] + 2*pi
testDf$relHead[testDf$relHead > pi] = testDf$relHead[testDf$relHead > pi] - 2*pi
testDf$spTrav <- allHrly$spTrav[!is.na(allHrly$aveTrHead)]

ggplot() +
  geom_path(aes(x = cord.dec$coords.x1, y = cord.dec$coords.x2)) +
  geom_point(aes(x = trun[[c]]$lon, y = trun[[c]]$lat, colour = trun[[c]]$hour_along))


ggplot(testDf[testDf$outRet == 1 & testDf$spTrav < 80,]) +
  geom_point(aes(x = relHead, y = spTrav)) +
  coord_polar() + scale_x_continuous(limits = c(-pi, pi))

ggplot() +
  geom_point(aes(x = (allHrly$aveHead[!is.na(allHrly$aveTrHead)] - allHrly$aveTrHead[!is.na(allHrly$aveTrHead)]),
    y = allHrly$spTrav[!is.na(allHrly$aveTrHead)], colour = allHrly$rtChg[!is.na(allHrly$aveTrHead)] > 0)) +
  coord_polar() + scale_x_continuous(limits = c(-2*pi,2*pi))
ggplot(allHrly[allHrly$rtChg < 0,]) +
  geom_point(aes(y = recalSp, x = (aveHead*(180/pi) - aveTrHead*(180/pi)))) +
  coord_polar() + scale_x_continuous(limits = c(0,360))


# testing ------------------------------------------------------------------------------------
sel <- allD[allD$tagID == allD$tagID[1] & allD$Year == allD$Year[1],]
tsel <- seq(sel$DT[1], sel$DT[length(sel$DT)], by = 6*3600)
dtFull <- as.numeric(diff(sel$DT))
# find data that line up to the new timepoints
select <- NA
for(b in 1:length(tsel)){
    choose <- which(sel$DT >= (tsel[b] - (median(dtFull))) & sel$DT <= (tsel[b] + (median(dtFull))))
    if(length(choose) != 0){
        if(length(choose) > 1){
            select[b] <- choose[which.min(abs(tsel[b] - sel$DT [choose]))]
        } else {
            select[b] <- choose
        }
    }
}
select <- select[!is.na(select)]

iso <- sel[select,]
iso <- iso[as.numeric(format(iso$DT, "%H")) < 18,]

trun <- vector(mode = "list", length = nrow(iso))
for(b in 1:length(select)){
  tryCatch({
    trun[[b]] <- TrackTraj(sel$DT[select[b]], sel$lat[select[b]], sel$lon[select[b]], 6)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# calculate the mean headings of the birds between each point
# compare between the wind trajectories and the bird - subtraction?
aveHead <- NA
aveTrHead <- NA
for(c in 1:(length(trun) - 1)){
  cord.dec <- SpatialPoints(cbind(sel$lon[select[c]:select[c+1]], sel$lat[select[c]:select[c+1]]),
              proj4string=CRS("+proj=longlat"))
  cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
  utmDiff <- cbind(diff(cord.UTM$coords.x1), diff(cord.UTM$coords.x2))
  aveHead[c] <- atan2(mean(utmDiff[,1]),mean(utmDiff[,2]))

  trCord.dec <- SpatialPoints(cbind(trun[[c+1]]$lon, trun[[c+1]]$lat),
              proj4string=CRS("+proj=longlat"))
  trCord.utm <- spTransform(rev(trCord.dec), CRS("+proj=utm +zone=54 +datum=WGS84"))
  trutmDiff <- cbind(diff(trCord.utm$coords.x1), diff(trCord.utm$coords.x2))
  aveTrHead[c] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
}
plot((aveHead - aveTrHead)*(180/pi))



tst <- TrackTraj(iso$DT[b], iso$lat[b], iso$lon[b], 6)

DT <- iso$DT[b]
lat <- iso$lat[b]
lon <- iso$lon[b]
hrs <- 6
time <- format(DT, "%Y-%m-%d")
  hr <- format(DT, "%H")
  trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = lat,
    lon = lon,
    height = 50,
    duration = hrs,
    days = time,
    daily_hours = hr,
    direction = "forward",
    met_type = "gdas1",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()
  return(trajectory_model %>% get_output_tbl())


  trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = 42.42768,
    lon = 143.8522,
    height = 10,
    duration = 6,
    days = "2018-08-29",
    daily_hours = 22,
    direction = "forward",
    met_type = "gdas1",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()

trun <- lapply(sel, function(x) TrackTraj(x$DT,x$lat,x$lon,6))

testRun <- lapply(sel[select,], function(x) TrackTraj(x[,1], x[,2], x[,3], 6))
TrackTraj(allD$DT[for1], allD$lat[for1], allD$lon[for1],6)



testRun <- TrackTraj(sel[select,]$DT,sel[select,]$Lat,sel[select,]$Lon,6)
