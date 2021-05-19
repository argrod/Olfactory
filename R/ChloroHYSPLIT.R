# remotes::install_github("ropensci/rerddap")
library(rerddap)
library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
# devtools::install_github('ropensci/plotdap')
library(plotdap)
library(scales)
# install.packages("StreamMetabolism",repos='http://cran.us.r-project.org')
library(poweRlaw)
library(ggpubr)
library(maptools)
library(StreamMetabolism)
library(RColorBrewer)
library(reshape)
# install.packages("circular", repos="http://R-Forge.R-project.org")
library(circular)
library(wGribDat)
# devtools::install_github("keesmulder/circglmbayes")
library(circglmbayes)
library(lme4)
require(devtools)
# install_version("ncdf4", version = "1.16", repos = "http://cran.us.r-project.org")
# library(ncdf4)
# install.packages("gdalUtils")
# devtools::install_github("mdsumner/ncdf4")
library(gdalUtils)
library(CircStats)
library(Cairo)
options(timeout = 800)
# load in data 
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
    spTrav = c(D18$spTrav, D19$spTrav),
    recalSp = c(D18$recalSp, D19$recalSp),
    distFk = c(D18$distFromFk, D19$distFromFk),
    tripL = c(D18$tripL, D19$tripL),
    tkb = c(D18$tkb, D19$tkb),
    dv = c(D18$dv, D19$dv),
    UTME = c(D18$UTME, D19$UTME),
    UTMN = c(D18$UTMN, D19$UTMN))
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")
# functions to extract HYSPLIT model data
find_hull <- function(df) df[chull(df$lon, df$lat), ]
TrackDisp <- function(DT, lat, lon, hrs){
  time <- format(DT, "%Y-%m-%d %H:%M") - lubridate::hours(9)
  dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = lat, lon = lon, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm(time) + lubridate::hours(hrs),
    release_end = lubridate::ymd_hm(time) + lubridate::hours(hrs) + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm(time) + lubridate::hours(hrs),
    end_time = lubridate::ymd_hm(time),
    direction = "forward", 
    met_type = "reanalysis",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()
  dispersion_tbl <- dispersion_model %>% get_output_tbl()
  hulls <- ddply(dispersion_tbl, "hour", find_hull)
  outpt <- list("DispModel" = dispersion_model,"partDisp" = dispersion_tbl, "partPoly" = hulls)
}
TrackTraj <- function(DT, lat, lon, hrs){
  time <- format(DT - lubridate::hours(9), "%Y-%m-%d") # convert to UTC
  hr <- format(DT, "%H")
  trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = lat,
    lon = lon,
    height = 10,
    duration = hrs,
    days = time,
    daily_hours = hr,
    direction = "backward",
    met_type = "reanalysis",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()
  return(trajectory_model %>% get_output_tbl())
}
ListD <- c(Dat,Dat19)
StLDisp <- function(utm){ # find changes in sign
  chgs <- diff(utm) >= 0
  swtch <- which(diff(chgs)!=0) + 1
}
StLCalc <- function(DT,UTMN,UTME){ # split data into steps, including with consistent sampling intervals
  deltaT <- difftime(DT[2:length(DT)], DT[1:(length(DT)-1)], units = "secs")
  # find the cutoff points from the time difference
  cutoff <- median(deltaT, na.rm = T) + 10
  # find the start and end times based on time differences over the cutoff value
  low<-(deltaT <= cutoff)
  strts <- which(diff(low) == 1) + 1
  if(deltaT[1] <= cutoff){
    strts <- c(1,strts)
  }
  ends <- which(diff(low) == -1) + 1
  if(strts[length(strts)] > ends[length(ends)]){
    ends <- c(ends, max(which(deltaT <= cutoff)) + 1)
  }
  # find the changes in sign of the data
  stepLsN <- NA
  for(g in 1:length(strts)){
    iso <- UTMN[strts[g]:ends[g]]
    stepLsN[g] <- sum(abs(diff(iso)))
  }
  stepLsE <- NA
  for(g in 1:length(strts)){
    iso <- UTME[strts[g]:ends[g]]
    stepLsE[g] <- sum(abs(diff(iso)))
  }
  data.frame(start=DT[strts], end=DT[ends], dispN = stepLsN, dispE = stepLsE, strtInd = strts, endInd = ends)
}
levy <- function(x,mu){
  x^-mu
}
outTraj <- vector(mode = "list", length = length(ListD))
for(b in 1:length(ListD)){
  sel<-ListD[[b]]
  sel$ForPoint <- 0
  sel$ForPoint[sel$tkb==1 | sel$dv == 1] <- 1
  sel$tToFor <- NA
  sel$tToFor[sel$ForPoint == 1] = 0
  sel$dToFor <- NA
  sel$dToFor[sel$ForPoint == 1] = 0
  forStart <- which(diff(sel$ForPoint == 1) == 1) + 1
  if(sel$ForPoint[1] == 1){
    forStart <- c(1, forStart)
  }
  forEnd <- which(diff(sel$ForPoint == 1) == -1) + 1
  for(starts in which(sel$ForPoint == 0)){
    if(any(sel$ForPoint[starts:nrow(sel)] == 1)){
      sel$tToFor[starts] <- difftime(sel$DT[starts-1+min(which(sel$ForPoint[starts:nrow(sel)] == 1))], sel$DT[starts], units = "secs")
      sel$dToFor[starts] <- sqrt((sel$UTMN[starts-1+min(which(sel$ForPoint[starts:nrow(sel)] == 1))] - sel$UTMN[starts])^2 + (sel$UTME[starts-1+min(which(sel$ForPoint[starts:nrow(sel)] == 1))] - sel$UTME[starts])^2)
    }
  }
  # remove non-flight values
  sel <- sel[which(sel$spTrav > 15),]
  # colnames(sel)
  sel$tdiff <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel) - 1)], units = 'secs'))
  # decide outgoing/incoming
  sel$appSpd <- c(NA, diff(sel[,grepl("Fk", names(sel))]))/sel$tdiff # approach speed to Fk Island
  sel$rtChg <- NA
  for(c in 1:nrow(sel)){
    sel$rtChg[c] <- mean(sel$appSpd[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]/as.numeric(sel$tdiff[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]))
  }
  # sel <- sel[which(sel$rtChg < 0),]
  UDsts <- StLCalc(sel$DT,sel$UTMN,sel$UTME)
  dt <- difftime(UDsts$end,UDsts$start,units="secs")
  aveHead <- NA
  trajHead <- NA
  trajSpd <- NA
  for(ind in 1:nrow(UDsts)){
    aveHead[ind] <- atan2(mean(diff(sel$UTMN[UDsts$strtInd[ind]:UDsts$endInd[ind]])), mean(diff(sel$UTME[UDsts$strtInd[ind]:UDsts$endInd[ind]])))
    trajs <- tryCatch({
      TrackTraj(sel$DT[UDsts$strtInd[ind]], sel$Lat[UDsts$strtInd[ind]], sel$Lon[UDsts$strtInd[ind]], 6)
    }, error = function(e){NA})
    if(is.na(trajs)){
      trajHead[ind] <- NA
      trajSpd[ind] <- NA
    } else {
      trCord.dec <- SpatialPoints(cbind(rev(trajs$lon), rev(trajs$lat)), proj4string=CRS("+proj=longlat"))
      trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
      trutmDiff <- cbind(diff(trCord.utm$coords.x2), diff(trCord.utm$coords.x1))
      trajHead[ind] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
      trajSpd[ind] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)}
  }
  outTraj[[b]] <- data.frame(DT = sel$DT[UDsts$strtInd], aveHd = aveHead, trjHd = trajHead, trjSpd = trajSpd, lat = sel$Lat[UDsts$strtInd], lon = sel$Lon[UDsts$strtInd], timeTo = sel$tToFor[UDsts$strtInd], distTo = sel$dToFor[UDsts$strtInd], rtChg = sel$rtChg[UDsts$strtInd])
}
ind=205
plot(sel$Lon[which(sel$DT==UDsts$start[ind]):which(sel$DT==UDsts$end[ind])],sel$Lat[which(sel$DT==UDsts$start[ind]):which(sel$DT==UDsts$end[ind])],type='l')
points(trajs$lon,trajs$lat)

ind<-100
ggplot(ListD[[1]][UDsts$strtInd[ind]:UDsts$endInd[ind],], aes(x=Lon,y=Lat)) +
  geom_path() +
  geom_point(data = outTraj[[1]][ind,], aes(x = lon, y = lat))
  geom_spoke(data = outTraj[[1]][ind,], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = scales::rescale(outTraj[[1]]$trjSpd[ind], c(.2, .8)))

plot(trajs$lon,trajs$lat)
for(b in 1:length(outTraj)){
  sel <- ListD[[b]]
  sel <- sel[which(sel$spTrav > 15),]
  UDsts <- StLCalc(sel$DT,sel$UTMN,sel$UTME)
  dist <- NA
  for(g in 1:nrow(UDsts)){
    dist[g] <- sum(sel$distTrav[(UDsts$strtInd[g]+1):UDsts$endInd[g]])
  }
  outTraj[[b]]$cumDist <- dist
}
# save(outTraj,file="/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/StepsTrajTimeChg.RData")
sel <- allD[allD$forage == 1,]
ggplot(sel, aes(x = lon, y = lat)) + geom_point()
# LOAD IN THE STEP LENGTHS TRAJECTORIES
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/StepsTrajTimeChg.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/StepsTrajTimeChg.RData")
}
# save(ListD, file="/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ListD.RData")
# LOAD IN THE LISTED DATA
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ListD.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/ListD.RData")
}
summary(outTraj)

ind<-100
ggplot(ListD[[1]][UDsts$strtInd[ind]:UDsts$endInd[ind],], aes(x=Lon,y=Lat)) +
  geom_path() +
  geom_point(data = outTraj[[1]][ind,], aes(x = lon, y = lat))
  geom_spoke(data = outTraj[[1]][ind,], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = scales::rescale(outTraj[[1]]$trjSpd[ind], c(.2, .8)))
ggplot(outTraj[[1]][outTraj[[1]]$distTo<10*10^3 & outTraj[[1]]$rtChg <0,],aes(x = lon,y=lat)) +
  geom_spoke(aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")), radius = scales::rescale(outTraj[[1]]$trjSpd[outTraj[[1]]$distTo<10*10^3 & outTraj[[1]]$rtChg <0], c(.2, .8))) +
  scale_colour_gradient(low = "yellow", high= "red") +
  geom_spoke(aes(x = lon, y = lat, angle = aveHd), arrow = arrow(length = unit(0.05,"inches")), radius = scales::rescale(outTraj[[1]]$trjSpd[outTraj[[1]]$distTo<10*10^3 & outTraj[[1]]$rtChg <0], c(.2, .8)))
colnames(outTraj[[1]])


# add the duration of the proceeding foraging point
allEth <- c(Eth,Eth19)
colnames(allEth[[1]])
for(b in 1:length(outTraj)){
  outTraj[[b]]$nxtForDur <- NA
  for(g in 1:nrow(outTraj[[b]])){
    if(any(allEth[[b]]$DT > outTraj[[b]]$DT[g], na.rm=T)){
      loc <- min(which(allEth[[b]]$DT > outTraj[[b]]$DT[g]))
      if(any(allEth[[b]]$cat[loc:nrow(allEth[[b]])] == "Forage", na.rm=T)){
        nxtFor <- min(which(allEth[[b]]$cat[loc:nrow(allEth[[b]])] == "Forage")) + loc - 1
        strT <- allEth[[b]]$DT[nxtFor]
        endT <- allEth[[b]]$DT[min(which(allEth[[b]]$cat[nxtFor:nrow(allEth[[b]])] != "Forage")) + nxtFor - 1]
        outTraj[[b]]$nxtForDur[g] <- as.numeric(difftime(endT,strT,units='secs'))
      }
    }
  }
}
# repeat for windDat
WindDat$yr <- format(WindDat$DT, format = "%Y")
indivs <- unique(WindDat[c("ID","yr")])
ethInds <- data.frame(ID=NA,DT=NA)
for(b in 1:length(allEth)){
  ethInds[b,] <- cbind(allEth[[b]]$tagID[1],format(allEth[[b]]$DT[1],format="%Y"))
}
WindDat$nxtForDur <- NA
for(b in 1:nrow(indivs)){
  if(length(which(WindDat$ID == indivs$ID[b] & WindDat$yr == indivs$yr[b])) != 0){
    selInd <- which(WindDat$ID == indivs$ID[b] & WindDat$yr == indivs$yr[b])
    for(g in 1:length(selInd)){
      EthSel <- allEth[[which(ethInds$ID == indivs$ID[b] & ethInds$DT==indivs$yr[b])]]
      loc <- min(which(EthSel$DT > WindDat$DT[selInd[g]]))
      if(any(EthSel$cat[loc:nrow(EthSel)] == "Forage", na.rm=T)){
        nxtFor <- min(which(EthSel$cat[loc:nrow(EthSel)] == "Forage")) + loc - 1
        strT <- EthSel$DT[nxtFor]
        endT <- EthSel$DT[min(which(EthSel$cat[nxtFor:nrow(EthSel)] != "Forage")) + nxtFor - 1]
        WindDat$nxtForDur[selInd[g]] <- as.numeric(difftime(endT,strT,units='secs'))
      }
    }
  }
}

EthSel$cat[nxtFor:min(which(EthSel$cat[nxtFor:nrow(EthSel)] != "Forage")) + nxtFor - 1]
EthSel$cat[min(which(EthSel$cat[nxtFor:nrow(EthSel)] != "Forage")) + nxtFor]
EthSel$cat[34:43]

colnames(Eth[[1]])
summary(Dat[[1]]$Dur)
for(b in 1:length(ListD)){
  print(ListD[[b]]$tagID[1])
  print(allEth[[b]]$tagID[1])
}

for(show in 1:length(outTraj)){
  outTraj[[show]]$relH <- outTraj[[show]]$aveHd - outTraj[[show]]$trjHd
  outTraj[[show]]$relH[outTraj[[show]]$relH < -pi] <- outTraj[[show]]$relH[outTraj[[show]]$relH < -pi] + 2*pi
  outTraj[[show]]$relH[outTraj[[show]]$relH > pi] <- outTraj[[show]]$relH[outTraj[[show]]$relH > pi] - 2*pi
  outTraj[[show]]$relW <- NA
  outTraj[[show]]$relW[outTraj[[show]]$relH < pi/4 & outTraj[[show]]$relH > -pi/4] <- "Front"
  outTraj[[show]]$relW[outTraj[[show]]$relH < -pi/4 & outTraj[[show]]$relH > -3*pi/4] <- "Side"
  outTraj[[show]]$relW[outTraj[[show]]$relH > pi/4 & outTraj[[show]]$relH < 3*pi/4] <- "Side"
  outTraj[[show]]$relW[outTraj[[show]]$relH < -3*pi/4 | outTraj[[show]]$relH > 3*pi/4] <- "Behind"
  outTraj[[show]]$Phase <- NA
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo < 10*60] <- "<10"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 10*60 & outTraj[[show]]$timeTo < 20*60] <- "10-20"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 20*60 & outTraj[[show]]$timeTo < 30*60] <- "20-30"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 30*60 & outTraj[[show]]$timeTo < 40*60] <- "30-40"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 40*60 & outTraj[[show]]$timeTo < 50*60] <- "40-50"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 50*60 & outTraj[[show]]$timeTo < 60*60] <- "50-60"
  outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 60*60] <- "60+"
  # outTraj[[show]]$Phase[outTraj[[show]]$timeTo < 30*60] <- "Near"
  # outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 30*60] <- "Transit"
  outTraj[[show]]$proxim <- "Far"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo < 10^3] <- "<10"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 10^3 & outTraj[[show]]$distTo < 20^3] <- "10-20"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 20^3 & outTraj[[show]]$distTo < 30^3] <- "20-30"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 30^3 & outTraj[[show]]$distTo < 40^3] <- "30-40"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 40^3 & outTraj[[show]]$distTo < 50^3] <- "40-50"
  outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 50^3 & outTraj[[show]]$distTo < 60^3] <- "50-60"
  # fr <- ggplot(outTraj[[show]][outTraj[[show]]$proxim=="Far",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # nr<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="<10",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # Ten20<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="10-20",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # Twenty30<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="20-30",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # Thirty40<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="30-40",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # Forty50<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="40-50",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # Fifty60<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="50-60",], aes(x = relH*(180/pi))) +
  #   geom_histogram() + coord_polar()
  # png(paste("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/",as.character(ListD[[show]]$tagID[1]),format(outTraj[[show]]$DT[1],"%Y"),"OffTraj.png",sep=""))
  # print(ggarrange(nr, Ten20, Twenty30, Thirty40, Forty50, Fifty60, ncol= 3, nrow = 2, labels = c("<10km", "10-20km", "20-30km","30-40km","40-50km","50-60km")))
  # dev.off()
}
show <- 6
UDsts <- StLCalc(ListD[[show]]$DT,ListD[[show]]$UTMN,ListD[[show]]$UTME)
ggplot(ListD[[show]]) +
  geom_path(aes(x = Lon, y = Lat)) +
  xlim(c(143.5,144.8)) + ylim(c(42,43)) +
  geom_point(data = outTraj[[show]][UDsts$strtInd,], aes(x = lon, y = lat)) +
  geom_spoke(data = outTraj[[show]], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = scales::rescale(outTraj[[show]]$trjSpd, c(.2, .8))) +
  scale_colour_distiller(palette = "RdYlGn", name = "Approx. speed") +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black'))

# outTraj[[4]]$relH <- outTraj[[4]]$aveHd - outTraj[[4]]$trjHd
# outTraj[[4]]$relH[outTraj[[4]]$relH < -pi] <- outTraj[[4]]$relH[outTraj[[4]]$relH < -pi] + 2*pi
# outTraj[[4]]$relH[outTraj[[4]]$relH > pi] <- outTraj[[4]]$relH[outTraj[[4]]$relH > pi] - 2*pi
# ggplot(outTraj[[4]][outTraj[[4]]$timeTo > (5*60),], aes(x = distTo, y = relH)) +
#   geom_point() + xlim(c(0,25000))

outTraj[[show]]$relH <- outTraj[[show]]$aveHd - outTraj[[show]]$trjHd
outTraj[[show]]$relH[outTraj[[show]]$relH < -pi] <- outTraj[[show]]$relH[outTraj[[show]]$relH < -pi] + 2*pi
outTraj[[show]]$relH[outTraj[[show]]$relH > pi] <- outTraj[[show]]$relH[outTraj[[show]]$relH > pi] - 2*pi
outTraj[[show]]$relW <- NA
outTraj[[show]]$relW[outTraj[[show]]$relH < pi/4 & outTraj[[show]]$relH > -pi/4] <- "Front"
outTraj[[show]]$relW[outTraj[[show]]$relH < -pi/4 & outTraj[[show]]$relH > -3*pi/4] <- "Side"
outTraj[[show]]$relW[outTraj[[show]]$relH > pi/4 & outTraj[[show]]$relH < 3*pi/4] <- "Side"
outTraj[[show]]$relW[outTraj[[show]]$relH < -3*pi/4 | outTraj[[show]]$relH > 3*pi/4] <- "Behind"
outTraj[[show]]$Phase <- NA
outTraj[[show]]$Phase[outTraj[[show]]$timeTo < 10*60] <- "<10"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 10*60 & outTraj[[show]]$timeTo < 20*60] <- "10-20"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 20*60 & outTraj[[show]]$timeTo < 30*60] <- "20-30"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 30*60 & outTraj[[show]]$timeTo < 40*60] <- "30-40"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 40*60 & outTraj[[show]]$timeTo < 50*60] <- "40-50"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 50*60 & outTraj[[show]]$timeTo < 60*60] <- "50-60"
outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 60*60] <- "60+"
# outTraj[[show]]$Phase[outTraj[[show]]$timeTo < 30*60] <- "Near"
# outTraj[[show]]$Phase[outTraj[[show]]$timeTo >= 30*60] <- "Transit"
outTraj[[show]]$proxim <- "Far"
outTraj[[show]]$proxim[outTraj[[show]]$distTo < 10^3] <- "<10"
outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 10^3 & outTraj[[show]]$distTo < 20^3] <- "10-20"
outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 20^3 & outTraj[[show]]$distTo < 30^3] <- "20-30"
outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 30^3 & outTraj[[show]]$distTo < 40^3] <- "30-40"
outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 40^3 & outTraj[[show]]$distTo < 50^3] <- "40-50"
outTraj[[show]]$proxim[outTraj[[show]]$distTo >= 50^3 & outTraj[[show]]$distTo < 60^3] <- "50-60"
# ggplot(outTraj[[show]], aes(y = trjSpd, x = relH*(180/pi), colour = relW)) + coord_polar() + scale_x_continuous(limits = c(-180,180)) + geom_point()

# ggplot(outTraj[[show]], aes(x = relH*(180/pi), fill = proxim)) +
#   geom_histogram(alpha = .5, position = 'dodge')
fr <- ggplot(outTraj[[show]][outTraj[[show]]$proxim=="Far",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
nr<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="<10",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Ten20<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="10-20",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Twenty30<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="20-30",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Thirty40<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="30-40",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Forty50<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="40-50",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Fifty60<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="50-60",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()

ggplot(outTraj[[show]], aes(x = distTo*10^-3, y = relH*(180/pi))) + 
  geom_point() + xlim(c(0,60))

ggarrange(nr, Ten20, Twenty30, Thirty40, Forty50, Fifty60, ncol= 3, nrow = 2, labels = c("<10km", "10-20km", "20-30km","30-40km","40-50km","50-60km"))
ggplot()+
  geom_path(data=ListD[[show]],aes(x=Lon,y=Lat)) +geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(141, 147), ylim = c(39, 44)) + geom_point(data=FkOshi,aes(x =Long,y=Lat,colour = "green"))


Less10min<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="<10",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Tenmin20<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="10-20",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Twentymin30<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="20-30",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Thirtymin40<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="30-40",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Fortymin50<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="40-50",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Fiftymin60<-ggplot(outTraj[[show]][outTraj[[show]]$Phase=="50-60",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
hourPlus <- ggplot(outTraj[[show]][outTraj[[show]]$Phase=="60+",], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()

ggarrange(Less10min, Tenmin20, Twentymin30, Thirtymin40, Fortymin50, Fiftymin60, hourPlus, ncol= 4, nrow = 2, labels = c("<10", "10-20", "20-30","30-40","40-50","50-60","60+"))


hist(outTraj[[show]]$relH[outTraj[[show]]$proxim == "Inter"]*(180/pi))

ggplot(outTraj[[show]], aes(x = distTo*10^-3, y = relH, colour = relW)) + geom_point()
ggplot(outTraj[[show]][outTraj[[show]]$proxim == "Inter",], aes(x = timeTo, y = relH, colour = relW)) + geom_point()

ggplot(outTraj[[4]], aes(x = lon, y = lat)) +
  geom_point() + geom_spoke(data = outTraj[[4]], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd+pi), arrow = arrow(length = unit(0.05,"inches")),
  radius = scales::rescale(outTraj[[4]]$trjSpd, c(.2, .8))) +
  scale_colour_distiller(palette = "RdYlGn", name = "Approx. speed") +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black'))

# add in information about outgoing or returning
for(b in 1:length(ListD)){
  tdiff <- c(NA, difftime(ListD[[b]]$DT[2:nrow(ListD[[b]])], ListD[[b]]$DT[1:(nrow(ListD[[b]]) - 1)], units = 'secs'))
  # decide outgoing/incoming
  ListD[[b]]$appSpd <- c(NA, diff(ListD[[b]][,grepl("Fk", names(ListD[[b]]))]))/tdiff # approach speed to Fk Island
  ListD[[b]]$rtChg <- NA
  for(c in 1:nrow(ListD[[b]])){
    ListD[[b]]$rtChg[c] <- mean(ListD[[b]]$appSpd[ListD[[b]]$DT >= ListD[[b]]$DT[c] & ListD[[b]]$DT <= (ListD[[b]]$DT[c] + lubridate::hours(1))]/as.numeric(tdiff[ListD[[b]]$DT >= ListD[[b]]$DT[c] & ListD[[b]]$DT <= (ListD[[b]]$DT[c] + lubridate::hours(1))]))
  }
}

# FIND THE OFFSET TRAJECTORY AND ASSIGN DEPART/RETURN VALUES (FOR PROCEEDING HOUR AS PER SHIOMI 2012)
for(b in 1:length(outTraj)){
  # outTraj[[b]]$relH <- outTraj[[b]]$aveHd - outTraj[[b]]$trjHd
  # outTraj[[b]]$relH[outTraj[[b]]$relH < -pi] <- outTraj[[b]]$relH[outTraj[[b]]$relH < -pi] + 2*pi
  # outTraj[[b]]$relH[outTraj[[b]]$relH > pi] <- outTraj[[b]]$relH[outTraj[[b]]$relH > pi] - 2*pi
  # outTraj[[b]]$relW <- NA
  # outTraj[[b]]$relW[outTraj[[b]]$relH < pi/4 & outTraj[[b]]$relH > -pi/4] <- "Front"
  # outTraj[[b]]$relW[outTraj[[b]]$relH < -pi/4 & outTraj[[b]]$relH > -3*pi/4] <- "Side"
  # outTraj[[b]]$relW[outTraj[[b]]$relH > pi/4 & outTraj[[b]]$relH < 3*pi/4] <- "Side"
  # outTraj[[b]]$relW[outTraj[[b]]$relH < -3*pi/4 | outTraj[[b]]$relH > 3*pi/4] <- "Behind"
  # outTraj[[b]]$Phase <- NA
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo < 10*60] <- "<10"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 10*60 & outTraj[[b]]$timeTo < 20*60] <- "10-20"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 20*60 & outTraj[[b]]$timeTo < 30*60] <- "20-30"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 30*60 & outTraj[[b]]$timeTo < 40*60] <- "30-40"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 40*60 & outTraj[[b]]$timeTo < 50*60] <- "40-50"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 50*60 & outTraj[[b]]$timeTo < 60*60] <- "50-60"
  # outTraj[[b]]$Phase[outTraj[[b]]$timeTo >= 60*60] <- "60+"
  outTraj[[b]]$proxim <- "20+"
  outTraj[[b]]$proxim[outTraj[[b]]$distTo < 1*10^3] <- "<1"
  outTraj[[b]]$proxim[outTraj[[b]]$distTo >= 1*10^3 & outTraj[[b]]$distTo < 5*10^3] <- "1-5"
  outTraj[[b]]$proxim[outTraj[[b]]$distTo >= 5*10^3 & outTraj[[b]]$distTo < 10*10^3] <- "5-10"
  outTraj[[b]]$proxim[outTraj[[b]]$distTo >= 10*10^3 & outTraj[[b]]$distTo < 15*10^3] <- "10-15"
  outTraj[[b]]$proxim[outTraj[[b]]$distTo >= 15*10^3 & outTraj[[b]]$distTo < 20*10^3] <- "15-20"
  outTraj[[b]]$rtChg <- NA
  for(g in 1:nrow(outTraj[[b]])){
    outTraj[[b]]$rtChg[g] <- ListD[[b]]$rtChg[ListD[[b]]$DT == outTraj[[b]]$DT[g]]
  }
}

allTraj <- bind_rows(outTraj)
allTraj$relH <- allTraj$aveHd - allTraj$trjHd
# allTraj$relH <- allTraj$relH + pi
fr <- ggplot(allTraj[allTraj$proxim=="20+" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
nr<-ggplot(allTraj[allTraj$proxim=="<1" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Ten20<-ggplot(allTraj[allTraj$proxim=="1-5" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Twenty30<-ggplot(allTraj[allTraj$proxim=="5-10" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Thirty40<-ggplot(allTraj[allTraj$proxim=="10-15" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
Forty50<-ggplot(allTraj[allTraj$proxim=="15-20" & allTraj$rtChg < 0,], aes(x = relH*(180/pi))) +
  geom_histogram() + coord_polar()
ggarrange(nr, Ten20, Twenty30, Thirty40, Forty50, ncol= 3, nrow = 2, labels = c("<1km", "1-5km", "5-10km","10-15km","15-20km"))


p1 <- ggplot(allTraj[allTraj$rtChg<0,]) + geom_point(aes(x = aveHd, y = log10(distTo)), colour = "red")+ coord_polar(start = pi)
p2 <- ggplot(allTraj[allTraj$rtChg<0,]) + geom_point(aes(x = trjHd, y = log10(distTo)), colour = "green") + coord_polar(start=pi)
ggarrange(p1, p2, ncol=1,nrow=2)

for(show in 1:length(outTraj)){
  fr <- ggplot(outTraj[[show]][outTraj[[show]]$proxim=="20+" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  nr<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="<1" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  Ten20<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="1-5" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  Twenty30<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="5-10" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  Thirty40<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="10-15" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  Forty50<-ggplot(outTraj[[show]][outTraj[[show]]$proxim=="15-20" & outTraj[[show]]$rtChg < 0,], aes(x = relH*(180/pi))) +
    geom_density() + coord_polar(start=((-pi/2)*180/pi)) + xlim(c(-180,180)) + ylim(c(0,0.0045))
  png(paste("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/Figures/Outward/",as.character(ListD[[show]]$tagID[1]),format(outTraj[[show]]$DT[1],"%Y"),"OutWardOffTraj.png",sep=""),width=800,height=800)
  print(ggarrange(nr, Ten20, Twenty30, Thirty40, Forty50, ncol= 3, nrow = 2, labels = c("<1km", "1-5km", "5-10km","10-15km","15-20km")))
  dev.off()
}

ggplot(allTraj[allTraj$distTo<20*10^3 & allTraj$rtChg<0,], aes(x = relH)) + geom_density() + coord_polar(start=-pi/2)


ListD[[show]]$tagID[1]
png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/6_S1AirTraj.png", width = 700, height = 800)
dev.off()
ggplot(Sel[1:50,], aes(x = Lon, y = Lat, fill = FlSpeed, angle = FlHead, 
    radius = scales::rescale(FlSpeed, c(.2, .8)))) +
    geom_raster() +
    geom_spoke(arrow = arrow(length = unit(.05, 'inches'))) + 
    scale_fill_distiller(palette = "RdYlGn") + 
    coord_equal(expand = 0) + 
    theme(legend.position = 'bottom', 
          legend.direction = 'horizontal')


tags <- unique(allD$tagID)
years <- unique(allD$Year)
sel <- as.data.frame(allD[allD$tagID == tags[1] & allD$Year == years[1],])
comps <- vector(mode="list",length=length(ListD))
powLs <- vector(mode="list",length=length(ListD))
expLs <- vector(mode="list",length=length(ListD))
normLs <- vector(mode="list",length=length(ListD))
for(b in 1:length(ListD)){
  sel <- ListD[[b]][which(ListD[[b]]$spTrav > 15),]
  # colnames(sel)
  sel$tdiff <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel) - 1)], units = 'secs'))
  # decide outgoing/incoming
  sel$appSpd <- c(NA, diff(sel[,grepl("Fk", names(sel))]))/sel$tdiff # approach speed to Fk Island
  sel$rtChg <- NA
  for(c in 1:nrow(sel)){
    sel$rtChg[c] <- mean(sel$appSpd[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]/as.numeric(sel$tdiff[sel$DT >= sel$DT[c] & sel$DT <= (sel$DT[c] + lubridate::hours(1))]))
  }
  # take outgoing (-ve approach speeds) only
  out <- sel[which(sel$rtChg < 0),]
  ret <- sel[which(sel$rtChg > 0),]
  selTest <- StLCalc(sel$DT, sel$UTMN, sel$UTME)
  outTest <- StLCalc(out$DT, out$UTMN, out$UTME)
  retTest <- StLCalc(ret$DT, ret$UTMN, ret$UTME)

  powLs[[b]] <- conpl$new(Nls$disp[Nls$disp>0])
  est <- estimate_xmin(powLs[[b]])
  powLs[[b]]$setXmin(est)

  expLs[[b]] <- conexp$new(Nls$disp[Nls$disp>0])
  expLs[[b]]$setXmin(estimate_xmin(expLs[[b]]))

  normLs[[b]] <- conlnorm$new(Nls$disp[Nls$disp>0])
  est <- estimate_xmin(normLs[[b]])
  normLs[[b]]$setXmin(est)

  # powLs[[b]] <- conpl$new(Nls$disp[Nls$disp>0])
  # m_bl$setXmin(estimate_xmin(powLs[[b]]))
  normLs[[b]] <- conlnorm$new(Nls$disp[Nls$disp>0])
  normLs[[b]]$setXmin(powLs[[b]]$getXmin())
  normLs[[b]]$setPars(estimate_pars(normLs[[b]]))
  comps[[b]] <- compare_distributions(powLs[[b]], normLs[[b]])
  comp$p_two_sided
  comp$test_statistic
}
png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/1_S1StepLengths.png", width = 700, height = 800)
plot(powLs[[b]], xlab = expression('log'[10]*' Step length'), ylab = expression('log'[10]*' Rank'))
lines(powLs[[b]], col = 2, lwd = 2)
lines(expLs[[b]], col=4, lwd= 2)
lines(normLs[[b]], col = 3, lwd = 2)
dev.off()

plot(powLs[[b]])
lines(powLs[[b]], col=2, lwd=2)
lines(normLs[[b]], col=3,lwd=2)
plot(comps[[b]])
# remove non-flight values


rank(selTest$dispN)/max(rank(selTest$dispN))
Nls <- data.frame(disp = sort(selTest$dispE, decreasing = T)*10^-3, rank = rev(rank(sort(selTest$dispE, decreasing = T)*10^-3)/max(rank(sort(selTest$dispN, decreasing = T)*10^-3))))
plot(log10(Nls$disp*10^-3), log10(levy(Nls$disp*10^-3,2)))
plot(log10(Nls$disp),log10(levy(sort(selTest$dispN, decreasing = T),2)))
p1<-ggplot(Nls, aes(x = (disp), y = (rank))) +
  geom_point(pch=21) + scale_x_log10() + scale_y_log10()
p1 + geom_line(Nls, mapping=aes(x = (disp), y = (levy(disp,1.5))), colour = "red")

p1 + geom_smooth(method="lm", formula = y~x)

plot(rev(cumsum(log(rev(Nls$disp[Nls$disp>0])))))
plot(rev(seq_along(Nls$disp[Nls$disp>0])))
# working with poweRlaw
m_bl <- conpl$new(Nls$disp[Nls$disp>0])
est <- estimate_xmin(m_bl)
m_bl$setXmin(est)
plot(m_bl)
lines(m_bl, col = 2, lwd = 2)

m_exp <- conexp$new(Nls$disp[Nls$disp>0])
m_exp$setXmin(estimate_xmin(m_exp))
lines(m_exp, col=4, lwd= 2)

m_bl_ln <- conlnorm$new(Nls$disp[Nls$disp>0])
est <- estimate_xmin(m_bl_ln)
m_bl_ln$setXmin(est)
lines(m_bl_ln, col = 3, lwd = 2)

m_bl <- conpl$new(Nls$disp[Nls$disp>0])
m_bl$setXmin(estimate_xmin(m_bl))
m_bl_ln <- conlnorm$new(Nls$disp[Nls$disp>0])
m_bl_ln$setXmin(m_bl$getXmin())
m_bl_ln$setPars(estimate_pars(m_bl_ln))
comp <- compare_distributions(m_bl, m_bl_ln)
comp$p_two_sided
comp$test_statistic
plot(m_bl)
lines(m_bl, col=2, lwd=2)
lines(m_bl_ln, col=3,lwd=2)
plot(comp)

bs <- bootstrap(m_bl, xmins = seq(5,50,1))
plot(bs, trim = 0.1)

bs_p <- bootstrap_p(m_bl, no_of_sims = 1000, threads = 2)
bs_p$p

plot(bs_p)

plot(m_bl)
lines(m_bl, col = 2, lwd = 2)

plot(sort(levy(1:1000,2),decreasing=F),type='l')

plot(Nls$disp[Nls$disp>5],Nls$rank[Nls$disp>5],log="xy")
lines(Nls$disp[Nls$disp>5], Nls$disp[Nls$disp>5]^(-1.45))

plot(1:1000,levy(1:1000,1.5),log="xy")
  geom_line(aes(x = log10(disp), y = log10(levy((1:nrow(Nls)),-1))))
  
  
  +
  geom_line(aes(x = (levy(log10(disp),2)), y = log10(disp), colour = "red"))
  geom_line(aes(x = log10(as.numeric(Var1)), y = (levy(as.numeric(Var1),2)),colour = "red"))
levy(as.numeric(Nls$Var1),2)

plot(log10(Nls$lengths),log10(Nls$rank))
lines(log10(Nls$lengths),levy(log10(Nls$lengths),3))
ggplot(Nls) + 
  geom_point(aes(x = (log10(rank)), y = (log10(disp)))) +
  geom_line(aes(x = log10(1:nrow(Nls))), y = levy((1:nrow(Nls)),-2))
  scale_x_continuous() +
  geom_function(fun = ~ log10((1:nrow(Nls))^2), aes(colour = "red"))

plot(selTest$dispN[order(rank(selTest$dispN))])

plot(log10(order(selTest$dispN)),log10(selTest$dispN[order((selTest$dispN))]))
lines(log10(c(1:length(selTest$dispN))^(2)))
lines(log10(selTest$dispN[order(selTest$dispN)]^(1)))

# calculate step lengths as per Humphries et al. 2013
sel <- as.data.frame(allD[allD$tagID == allD$tagID[1] & allD$Year == allD$Year[1],])
sel$deltaT <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel)-1)], units = "secs"))
# point out where the time difference exceeds the median value
cutoff <- median(sel$deltaT, na.rm = T) + 10
stepAreas <- c(0,which(sel$deltaT > cutoff),length(sel))
stepLsN <- vector(mode = "list", length = length(stepAreas))
stepLsE <- vector(mode = "list", length = length(stepAreas))


tstr<-data.frame(dat=c(4,5,9,7,2,-4,-6,-2,3,6,5,2,4,6,9,3,0,-4,-8,-2,6,8,10,5,3,1,2,2),swtch=NA,
  dt=c(30,30,30,30,30,51,30,30,30,30,30,30,30,62,182,30,30,30,30,30,30,30,30,30,651,30,30,30))
plot(tstr,type="l")
points(which(diff(diff(tstr$dat) >= 0)!=0) + 1,tstr$dat[which(diff(diff(tstr$dat) >= 0)!=0) + 1])
tstr$swtch<-0
tstr$swtch[which(diff(diff(tstr$dat) >= 0)!=0) + 1]=1
swtch<-which(diff(diff(tstr) >= 0)!=0) + 1
strts <- which(tstr$dt > cutoff)+1
if(strts[length(strts)] > length(tstr$dt)){
  strts <- strts[1:(length(strts)-1)]
}
strts <- c(1,strts[!tstr$dt[strts] > cutoff])
ends <- NA
for(sts in strts){
  ends[sts] <- min(min(which(tstr$dt[sts:length(tstr$dt)] > cutoff))-1, nrow(tstr))
}
ends <- strts+sapply(strts, function(x) min(which(tstr$dt[x:length(tstr$dt)] > cutoff)))-1
tstr$DT <- c(0,cumsum(tstr$dt[2:nrow(tstr)]))
plot(tstr$DT, tstr$dt,type="l")
points(tstr$DT[strts],tstr$dt[strts])
points(tstr$DT[ends],tstr$dt[ends],pch=2)
points(which(tstr$dt > cutoff),tstr$dt[which(tstr$dt > cutoff)])
points(which(tstr$dt > cutoff)-1,tstr$dt[which(tstr$dt > cutoff)-1],pch=2)
sel$DT[81:120]

ggplot() +
  geom_point(aes(x=log10(levy(1:10000,2)), y = log10(1:10000)))
   +
  scale_y_log10() + scale_x_log10()
sel$deltaT[stepAreas]

df <- data.frame(raw=1:10000,lev=levy(1:10000,2))
ggplot(df) + geom_line(aes(y = log10(rank(desc(lev))), x =(lev)))
plot(log10(df$lev))
for(c in 1:(length(stepAreas) - 1)){
  isol <- diff(sel$UTMN[(stepAreas[c]0):stepAreas[c+1]])
  if(length(isol) == 0){
    next
  } else {
    zcs <- c(1, StLCalc(isol), length(isol))
    stepLsN[[c]] <- sapply(c(1:(length(zcs) - 1)), function(x) sum(abs(isol[zcs[x]:(zcs[x+1]-1)])))
  }
}


for(b in 1:(length(stepAreas)-1)){
diff(sel$UTMN[stepAreas[b]:stepAreas[b+1]]) > 0


  isol <- which(diff(sel$UTMN[stepAreas[b]:stepAreas[b + 1]]))
}

ggplot() + 
  geom_line(aes(x = sel$DT[2:nrow(sel)], y = diff(sel$UTMN)))

# CALCULATE RELATIVE WIND CONDITIONS FROM ESTIMATES WITH TIME/DISTANCE TO FORAGING
if(Sys.info()['sysname'] == "Darwin"){
    windLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
} else {
    windLoc = "F:/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/MinDat/"
}
windFiles <- dir(windLoc)
for(b in 1:length(windFiles)){
  if(b == 1){
    WindDat <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = F)
    colnames(WindDat) <- c("DT","Lat","Lon","Head","X","Y")
    WindDat$ID <- sub("*.csv", "", windFiles[b])
    Wind.dec <- SpatialPoints(cbind(WindDat$Lon,WindDat$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    WindDat$UTME <- coordinates(UTMdat)[, 1]
    WindDat$UTMN <- coordinates(UTMdat)[, 2] 
  } else {
    toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
    colnames(toAdd) <- c("DT","Lat","Lon","Head","X","Y")
    toAdd$ID <- sub("*.csv", "", windFiles[b])
    Add.dec <- SpatialPoints(cbind(toAdd$Lon,toAdd$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Add.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    toAdd$UTME <- coordinates(UTMdat)[, 1]
    toAdd$UTMN <- coordinates(UTMdat)[, 2] 
    WindDat <- rbind(WindDat, toAdd)
  }
}
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y-%m-%d %H:%M:%OS")
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
#REPEAT FOR 2019
if(Sys.info()['sysname'] == "Darwin"){
    windLoc19 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
} else {
    windLoc19 = "F:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/MinDat/"
}
windFiles19 <- dir(windLoc19,pattern=".csv")
for(b in 1:length(windFiles19)){
  if(b == 1){
    WindDat19 <- read.delim(paste(windLoc19, windFiles19[b], sep = ''), sep = ",", header = F)
    colnames(WindDat19) <- c("DT","Lat","Lon","Head","X","Y")
    WindDat19$ID <- sub("*.csv", "", windFiles19[b])
    Wind.dec <- SpatialPoints(cbind(WindDat19$Lon,WindDat19$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    WindDat19$UTME <- coordinates(UTMdat)[, 1]
    WindDat19$UTMN <- coordinates(UTMdat)[, 2] 
  } else {
    toAdd <- read.delim(paste(windLoc19, windFiles19[b], sep = ''), sep = ",", header = T)
    colnames(toAdd) <- c("DT","Lat","Lon","Head","X","Y")
    toAdd$ID <- sub("*.csv", "", windFiles19[b])
    Add.dec <- SpatialPoints(cbind(toAdd$Lon,toAdd$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Add.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    toAdd$UTME <- coordinates(UTMdat)[, 1]
    toAdd$UTMN <- coordinates(UTMdat)[, 2] 
    WindDat19 <- rbind(WindDat19, toAdd)
  }
}
WindDat19$DT <- as.POSIXct(WindDat19$DT, format = "%Y-%m-%d %H:%M:%OS")
WindDat19$timeTo <- NA
WindDat19$distTo <- NA

tags <- unique(WindDat19$ID)
for(g in 1:nrow(WindDat19)){
  if(any(allD$forage[allD$tagID == WindDat19$ID[g] & allD$Year == format(WindDat19$DT[g], "%Y") & allD$DT > WindDat19$DT[g]] == 1)){
    point <- which(allD$lat == WindDat19$Lat[g] & allD$lon == WindDat19$Lon[g] & allD$tagID == WindDat19$ID[g] & allD$Year == format(WindDat19$DT[b], "%Y"))
    forInd <- min(which(allD$forage[point:max(which(allD$tagID == WindDat19$ID[g] & allD$Year == format(WindDat19$DT[b], "%Y")))] == 1)) + point - 1
    WindDat19$timeTo[g] <- as.numeric(difftime(allD$DT[point+forInd],WindDat19$DT[g], units="secs"))
    WindDat19$distTo[g] <- sqrt((allD$UTMN[forInd] - allD$UTMN[point])^2 + (allD$UTME[forInd] - allD$UTME[point])^2)*10^-3
  } else {
    next
  }
}
# save(WindDat19,file="/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDat19.RData")

if(Sys.info()['sysname'] == "Darwin"){
  load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDat.RData")
  load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDat19.RData")
} else {
  load("F:/UTokyoDrive/PhD/Data/WindCalc/windDat.RData")
  load("F:/UTokyoDrive/PhD/Data/WindCalc/windDat19.RData")
}
# remove data points after the last foraging points
WindDat <- WindDat[!is.na(WindDat$distTo),]
WindDat19 <- WindDat19[!is.na(WindDat19$distTo),]
WindDat$WHd <- atan2(WindDat$X,WindDat$Y)
WindDat$RelHead <- WindDat$Head-WindDat$WHd
WindDat$RelHead[WindDat$RelHead < -pi] <- WindDat$RelHead[WindDat$RelHead < -pi] + 2*pi
WindDat$RelHead[WindDat$RelHead > pi] <- WindDat$RelHead[WindDat$RelHead > pi] - 2*pi
WindDat19$WHd <- atan2(WindDat19$X,WindDat19$Y)
WindDat19$RelHead <- WindDat19$Head-WindDat19$WHd
WindDat19$RelHead[WindDat19$RelHead < -pi] <- WindDat19$RelHead[WindDat19$RelHead < -pi] + 2*pi
WindDat19$RelHead[WindDat19$RelHead > pi] <- WindDat19$RelHead[WindDat19$RelHead > pi] - 2*pi

# combine '18, '19 data
WindDat<-rbind(WindDat,WindDat19)
WindDat$distFromFk <- NA
WindDat$rtChg <- NA
totalDat <- bind_rows(ListD)
for(b in 1:nrow(WindDat)){
  WindDat$distFromFk[b] <- totalDat$distFromFk[which(totalDat$DT > WindDat$DT[b] - lubridate::minutes(1) & totalDat$DT < WindDat$DT[b] + lubridate::minutes(1) & totalDat$Lat == WindDat$Lat[b] & totalDat$Lon == WindDat$Lon[b])]
  WindDat$rtChg[b] <- totalDat$rtChg[which(totalDat$DT > WindDat$DT[b] - lubridate::minutes(1) & totalDat$DT < WindDat$DT[b] + lubridate::minutes(1) & totalDat$Lat == WindDat$Lat[b] & totalDat$Lon == WindDat$Lon[b])]
}
WindDat$yrID <- NA
for(b in 1:nrow(WindDat)){
  WindDat$yrID[b] <- paste(totalDat$tagID[which(totalDat$DT > WindDat$DT[b] - lubridate::minutes(1) & totalDat$DT < WindDat$DT[b] + lubridate::minutes(1) & totalDat$Lat == WindDat$Lat[b] & totalDat$Lon == WindDat$Lon[b])], format(totalDat$DT[which(totalDat$DT > WindDat$DT[b] - lubridate::minutes(1) & totalDat$DT < WindDat$DT[b] + lubridate::minutes(1) & totalDat$Lat == WindDat$Lat[b] & totalDat$Lon == WindDat$Lon[b])],format="%Y"),sep="")
}
# save(WindDat,file="F:/UTokyoDrive/PhD/Data/WindCalc/windDatAll.RData")
# LOAD WIND DATA
if(Sys.info()['sysname'] == "Darwin"){
  load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDatAll.RData")
} else {
  load("F:/UTokyoDrive/PhD/Data/WindCalc/windDatAll.RData")
}


ggplot(WindDat, aes(y = distTo, x = Head)) + geom_point() + coord_polar()
ggplot(WindDat, aes(y = distTo, x = WHd)) + geom_point() + coord_polar()
ggplot(WindDat[WindDat$distTo < 10,], aes(x = RelHead)) + geom_histogram() + coord_polar(start=-pi/2) + xlim(c(-pi,pi))
Less1<- ggplot(WindDat[WindDat$distTo < 1,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) +
  labs(x="Relative wind heading",y="Count")
  rose.diag(WindDat$RelHead[WindDat$distTo <10], bins = 30)
  rose.diag(WindDat$RelHead[WindDat$distTo >=10 & WindDat$distTo <20])
One2<- ggplot(WindDat[WindDat$distTo < 2 & WindDat$distTo >= 1,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Two3<- ggplot(WindDat[WindDat$distTo < 3 & WindDat$distTo >= 2,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Three4<- ggplot(WindDat[WindDat$distTo < 4 & WindDat$distTo >= 3,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Four5<- ggplot(WindDat[WindDat$distTo < 5 & WindDat$distTo >= 4,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Five6<- ggplot(WindDat[WindDat$distTo < 6 & WindDat$distTo >= 5,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Six7<- ggplot(WindDat[WindDat$distTo < 7 & WindDat$distTo >= 6,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Seven8<- ggplot(WindDat[WindDat$distTo < 8 & WindDat$distTo >= 7,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Eight9<- ggplot(WindDat[WindDat$distTo < 9 & WindDat$distTo >= 8,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
Nine10<- ggplot(WindDat[WindDat$distTo < 10 & WindDat$distTo >= 9,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
png("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/1-5km.png",width=800,height=800)
ggarrange(Less1,One2,Two3,Three4,Four5,Five6, ncol=3,nrow=2, labels=c("<1km","1-2km","2-3km","3-4km","4-5km","5-6km"))
dev.off()
png("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/5-10km.png",width=800,height=800)
ggarrange(Four5,Five6,Six7,Seven8,Eight9,Nine10, ncol=3,nrow=2, labels=c("4-5km","5-6km","6-7km","7-8km","8-9km","9-10km"))
dev.off()


ggplot(WindDat,aes(y = nxtForDur, x = RelHead)) + geom_point() + coord_polar()
hist(WindDat$RelHead)

# align the data for plotting
WindDat$aligned <- WindDat$RelHead + pi
WindDat$aligned[WindDat$aligned > pi] <- WindDat$aligned[WindDat$aligned > pi] - 2*pi
# bin distances into kilometres 
breaks<- c(0,1,2,3,4,5,6,7,8,9,10,max(WindDat$distTo,na.omit=T)+1)
WindDat$bins<-cut(WindDat$distTo, breaks=breaks,include.lowest = T, right=F)
#calculate mean of each group
mnW <- ddply(WindDat, "bins", summarise, grp.mean=mean(aligned))
png()
ggplot(WindDat[WindDat$distTo < 10,], aes(x = aligned, colour = bins)) + geom_density(alpha=.2,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bins), geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("Density") + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black')) +
  scale_colour_discrete(name="Dist to foraging (km)")

breaks<-seq(from=0,to=round_any(max(WindDat$distTo),10,f=ceiling),by=10)
mnW <- ddply(WindDat, "bin10", summarise, grp.mean=mean(aligned))


ggplot(mnW, aes(y = 10*as.numeric(bin10), x = grp.mean)) + 
  geom_point() +
  scale_x_continuous(limits = c(-pi, pi),name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail"))

# bin data every 10 km
WindDat$bin10 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)
png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/DistRelDensity.png",width=800,height=800)
ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(WindDat$distTo),], aes(x = aligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.2,show.legend=FALSE)+stat_density(aes(x=aligned, colour=bin10), geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("Density") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 20,
        family = "Arial"), axis.text = element_text(size = 20, family = "Arial")) + 
  scale_colour_manual(name="Dist to foraging (km)", values = rev(brewer.pal(9,"YlOrRd")))
dev.off()

ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90,], aes(x = aligned, fill = bin10)) +#max(WindDat$distTo),], aes(x = aligned, colour = bin10)) +
  geom_density(alpha=.4) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("Density") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 20,
        family = "Arial"), axis.text = element_text(size = 20, family = "Arial"))
  
  geom_histogram(alpha=.2,aes(fill = bin10))
  geom_density(alpha=.2,show.legend=FALSE)+stat_density(aes(x=aligned, colour=bin10), geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("Density") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 20,
        family = "Arial"), axis.text = element_text(size = 20, family = "Arial")) + 
  scale_colour_manual(name="Dist to foraging (km)", values = rev(brewer.pal(9,"YlOrRd")))

AllTraj <- allTraj
TrLess1<- ggplot(AllTraj[AllTraj$distTo < 1*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) +
  labs(x="Relative wind heading",y="Count")
TrOne2<- ggplot(AllTraj[AllTraj$distTo < 2*10^3 & AllTraj$distTo >= 1*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrTwo3<- ggplot(AllTraj[AllTraj$distTo < 3*10^3 & AllTraj$distTo >= 2*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrThree4<- ggplot(AllTraj[AllTraj$distTo < 4*10^3 & AllTraj$distTo >= 3*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrFour5<- ggplot(AllTraj[AllTraj$distTo < 5*10^3 & AllTraj$distTo >= 4*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrFive6<- ggplot(AllTraj[AllTraj$distTo < 6*10^3 & AllTraj$distTo >= 5*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrSix7<- ggplot(AllTraj[AllTraj$distTo < 7*10^3 & AllTraj$distTo >= 6*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrSeven8<- ggplot(AllTraj[AllTraj$distTo < 8*10^3 & AllTraj$distTo >= 7*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrEight9<- ggplot(AllTraj[AllTraj$distTo < 9*10^3 & AllTraj$distTo >= 8*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
TrNine10<- ggplot(AllTraj[AllTraj$distTo < 10*10^3 & AllTraj$distTo >= 9*10^3,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + labs(x="Relative wind heading",y="Count")
# png("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/1-5km.png",width=800,height=800)
ggarrange(TrLess1,TrOne2,TrTwo3,TrThree4,TrFour5,TrFive6, ncol=3,nrow=2, labels=c("<1km","1-2km","2-3km","3-4km","4-5km","5-6km"))
# dev.off()
# png("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/5-10km.png",width=800,height=800)
ggarrange(TrFour5,TrFive6,TrSix7,TrSeven8,TrEight9,TrNine10, ncol=3,nrow=2, labels=c("4-5km","5-6km","6-7km","7-8km","8-9km","9-10km"))
# dev.off()
TrZero10 <- ggplot(AllTraj[AllTraj$distTo >= 0 & AllTraj$distTo < 10,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrTen20 <- ggplot(AllTraj[AllTraj$distTo >= 10 & AllTraj$distTo < 20,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrTwenty30 <- ggplot(AllTraj[AllTraj$distTo >= 20 & AllTraj$distTo < 30,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrThirty40 <- ggplot(AllTraj[AllTraj$distTo >= 30 & AllTraj$distTo < 40,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))

# png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/0-40kmTraj.png", width = 700, height = 700)
ggarrange(TrZero10,TrTen20,TrTwenty30,TrThirty40,ncol=2,nrow=2,labels=c("0-10km","10-20km","20-30km","30-40km"))
# dev.off()

for(b in 1:length(outTraj)){
  outTraj[[b]]$travDir <- NA
  outTraj[[b]]$travDir[outTraj[[b]]$relH > -pi/4 & outTraj[[b]]$relH < pi/4] <- "Away"
  outTraj[[b]]$travDir[outTraj[[b]]$relH > pi/4 & outTraj[[b]]$relH < 3*pi/4] <- "R side"
  outTraj[[b]]$travDir[outTraj[[b]]$relH > -3*pi/4 & outTraj[[b]]$relH < -pi/4] <- "L side"
  outTraj[[b]]$travDir[outTraj[[b]]$relH < -3*pi/4 | outTraj[[b]]$relH > 3*pi/4] <- "Toward"
}
b=1
ggplot(outTraj[[b]][outTraj[[b]]$rtChg<0 & outTraj[[b]]$distTo < 10*10^3 & outTraj[[b]]$distTo > 1*10^3,], aes(x =(distTo))) + geom_density(aes(fill = travDir), alpha = .5)

# find the relative directions for the wind data
WindDat$travDir <- NA
WindDat$travDir[WindDat$RelHead > -pi/4 & WindDat$RelHead < pi/4] <- "Away"
WindDat$travDir[WindDat$RelHead > pi/4 & WindDat$RelHead < 3*pi/4] <- "R side"
WindDat$travDir[WindDat$RelHead > -3*pi/4 & WindDat$RelHead < -pi/4] <- "L side"
WindDat$travDir[WindDat$RelHead < -3*pi/4 | WindDat$RelHead > 3*pi/4] <- "Toward"
WindDat19$travDir <- NA
WindDat19$travDir[WindDat19$RelHead > -pi/4 & WindDat19$RelHead < pi/4] <- "Away"
WindDat19$travDir[WindDat19$RelHead > pi/4 & WindDat19$RelHead < 3*pi/4] <- "R side"
WindDat19$travDir[WindDat19$RelHead > -3*pi/4 & WindDat19$RelHead < -pi/4] <- "L side"
WindDat19$travDir[WindDat19$RelHead < -3*pi/4 | WindDat19$RelHead > 3*pi/4] <- "Toward"
ggplot(data=WindDat[WindDat$distTo < 10,], aes(x = distTo, fill = as.factor(travDir))) +
  geom_density(alpha=.5)

plot(WindDat$distTo[WindDat$distTo < 10],WindDat$RelHead[WindDat$distTo < 10])
ggplot(data=WindDat[WindDat$distTo < 10,], aes(x = RelHead, y = distTo, colour = travDir)) +
  geom_point() + coord_polar()

allTraj <- bind_rows(outTraj)
allTraj <- allTraj[allTraj$distTo != 0,]
# remove where the trajectory model failed
allTraj <- allTraj[!is.na(allTraj$trjHd),]
allTraj$relH <- allTraj$aveHd - allTraj$trjHd
#calculate sunrise/sunset times
sriseset<-cbind(rep(NA,nrow(allTraj)),rep(NA,nrow(allTraj)))
for(b in 1:nrow(allTraj)){
  sriseset[b,] <- cbind(sunrise.set(allTraj$lat[b], allTraj$lon[b], format(allTraj$DT[b],format="%Y/%m/%d"), timezone="Asia/Tokyo")$sunrise,sunrise.set(allTraj$lat[b], allTraj$lon[b], format(allTraj$DT[b],format="%Y/%m/%d"), timezone="Asia/Tokyo")$sunset)
}


ggplot(allTraj[allTraj$rtChg<0,], aes(x = log10(cumDist))) + geom_density()
sum(is.na(allTraj$rtChg))
colnames(allTraj)

ggplot(allTraj[allTraj$rtChg<0,], aes(x = trjHd, y = distTo)) + geom_point() + coord_polar(start=-pi/2)
ggplot(allTraj[allTraj$rtChg<0,], aes(x = aveHd, y = distTo)) + geom_point() + coord_polar(start=-pi/2)

ggplot(allTraj[allTraj$rtChg<0 & allTraj$distTo <= 1*10^3,], aes(x = distTo)) + geom_density(aes(fill=travDir),alpha=.5)

# calculate high concentration chloro-A locations and test average air trajectories vs bird headings

dInfo <- info("erdMBchla3day")

format(range(ListD[[1]]$DT), format="%Y-%m-%d")
res <- griddap('erdMBchla8day',time = format(range(allD$DT), format="%Y-%m-%d"),latitude = c(min(D18$Lat), max(D18$Lat)),
    longitude = c(min(D18$Lon), max(D18$Lon)))
myFunc <- function(x) log(x)

format(allD$DT, format ="%Y-%m-%d" & min(abs(allD$Lat - res[[2]]$lat)) & min(abs(allD$Lon - res[[2]]$lon))


days <- unique(res[[2]]$time)
lapply(days, function(x) which(res[[2]]$chlorophyll[res[[2]]$time == x] == max(res[[2]]$chlorophyll[res[[2]]$time == x], na.rm = T)))
res[[2]][39533,c(2,3)]
summary(res[[2]])

ggplot(res[[2]][res[[2]]$time == days[1],], aes(x = lon, y = lat, fill = chlorophyll)) +
    geom_tile() + scale_fill_gradient(trans = 'log')


hist(log10(res[[2]]$chlorophyll[res[[2]]$time == days[1]]),1000)

lat <- ListD[[1]]$Lat
lon <- ListD[[1]]$Lon
trCord.dec <- SpatialPoints(cbind((lon), (lat)), proj4string=CRS("+proj=longlat"))
trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
norf <- trCord.utm$coords.x1
est <- trCord.utm$coords.x2
plot(lon,lat)
plot(norf,est)

AllTraj[1,]
a <- 130
traj <- TrackTraj(AllTraj$DT[a],AllTraj$lat[a],AllTraj$lon[a],6)
disp <- TrackDisp(AllTraj$DT[a],AllTraj$lat[a],AllTraj$lon[a],6)
gdtraj <- create_trajectory_model() %>%
  add_trajectory_params(
    lat = AllTraj$lat[4029],
    lon = AllTraj$lon[4029],
    height = 10,
    duration = 6,
    days = "2018-09-05",
    daily_hours = 12,
    direction = "backward",
    met_type = "gdas1",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model() %>% get_output_tbl()
traj %>% trajectory_plot()

a=4029
ggplot(AllTraj[a,]) +
  geom_spoke(data = AllTraj[a,], aes(x = lon, y = lat, colour = trjSpd, angle = aveHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = .5) +
  geom_point(data=traj, aes(x=lon,y=lat,shape=as.factor(traj_dt))) +
  geom_point(data=gdtraj, aes(x=lon,y=lat,shape=as.factor(traj_dt))) +
  scale_colour_distiller(palette="RdYlGn")

disp <- create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = 43, lon = 145, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm("2018-09-05 15:00", tz = "Asia/Tokyo") + lubridate::hours(6),
    release_end = lubridate::ymd_hm("2018-09-05 15:00", tz = "Asia/Tokyo") + lubridate::hours(6) + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm("2018-09-05 15:00", tz = "Asia/Tokyo") + lubridate::hours(6),
    end_time = lubridate::ymd_hm("2018-09-05 15:00", tz = "Asia/Tokyo"),
    direction = "backward", 
    met_type = "reanalysis",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model() %>% get_output_tbl()
disp %>% dispersion_plot()
disp
hist(disp$height)


wGribDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/gribSelected.csv", sep = ",")
wGribDat$DT <- as.POSIXct(wGribDat$DT, format="%Y-%m-%dT%H:%M:%S")
wGribDat$trajHead <- NA
wGribDat$trajSpd <- NA
for(ind in 1:nrow(wGribDat)){
  trajs <- tryCatch({TrackTraj(wGribDat$DT[ind], wGribDat$Lat[ind], wGribDat$Lon[ind], 6)
    }, error = function(e){c(NA)})
  if(!is.na(trajs)){
    trCord.dec <- SpatialPoints(cbind(rev(trajs$lon), rev(trajs$lat)), proj4string=CRS("+proj=longlat"))
    trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    trutmDiff <- cbind(diff(trCord.utm$coords.x2), diff(trCord.utm$coords.x1))
    wGribDat$trajHead[ind] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
    wGribDat$trajSpd[ind] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)
    } else {
      wGribDat$trajHead[ind] <- NA
      wGribDat$trajSpd[ind] <- NA
    }
}
wGribDat$estHead <- atan2(wGribDat$Y,wGribDat$X)
colnames(wGribDat)
ggplot(wGribDat, aes(x = WindHead, y = trajHead)) + geom_point() + geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y))
res<-cor.circular(wGribDat$estHead[!is.na(wGribDat$trajHead)], wGribDat$trajHead[!is.na(wGribDat$trajHead)], test = T)

res<-cor.circular(wGribDat$WindHead[!is.na(wGribDat$trajHead)], wGribDat$trajHead[!is.na(wGribDat$trajHead)], test = T)

res<-cor.circular(wGribDat$estHead, wGribDat$WindHead, test = T)

mdata <- melt(wGribDat[,2:ncol(wGribDat)], id=c("WindHead"))
png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/TrajEstcomparison.png", height = 800, width = 800)
ggplot(mdata[mdata$variable == "trajHead" | mdata$variable == "estHead",], aes(x = WindHead, y = value)) +
  geom_point(aes(fill = variable), shape = 21,size=3) + geom_line(data = data.frame("x" = -pi:pi, "y" = -pi:pi), aes(x = x, y = y)) +
  scale_y_continuous(name = "Estimated headings") + scale_x_continuous(name = "JMA wind headings") +
  scale_fill_manual(name = "", values=c("deepskyblue","red"), labels = c("Trajectory","Track estimate")) +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 20,
        family = "Arial"), axis.text = element_text(size = 20, family = "Arial"))
dev.off()


# TEST STRAIGHTNESS OF STEP LENGTHS
UDsts <- vector(mode="list", length=length(ListD))
for(b in 1:length(ListD)){
  sel<-ListD[[b]]
  # remove non-flight values
  sel <- sel[which(sel$spTrav > 15),]
  # sel <- sel[which(sel$rtChg < 0),]
  UDsts[[b]] <- StLCalc(sel$DT,sel$UTMN,sel$UTME)
}
colnames(UDsts[[1]])
for(b in 1:length(UDsts)){
  UDsts[[b]]$dist <- NA
  for(g in 1:nrow(UDsts[[b]])){
    UDsts[[b]]$dist[g] <- sqrt((ListD[[b]]$UTMN[UDsts[[b]]$strtInd[g]] - ListD[[b]]$UTMN[UDsts[[b]]$endInd[g]])^2 + (ListD[[b]]$UTME[UDsts[[b]]$strtInd[g]] - ListD[[b]]$UTME[UDsts[[b]]$endInd[g]])^2)
  }
  outTraj[[b]]$dist <- UDsts[[b]]$dist
}


ggplot(UDsts[[1]]) +
  geom_segment(aes(x = ListD[[1]]$Lon[strtInd], y = ListD[[1]]$Lat[strtInd], xend = ListD[[1]]$Lon[endInd], yend = ListD[[1]]$Lat[endInd]))

png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/TrackVTraj.png", width = 800, height = 900)
ggplot(wGribDat, aes(x = estHead, y = trajHead)) +
  geom_point(shape = "o",size=3) + geom_line(data = data.frame("x" = -pi:pi, "y" = -pi:pi), aes(x = x, y = y)) +
  scale_y_continuous(name = "Trajectory headings") + scale_x_continuous(name = "Track headings") +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 20,
        family = "Arial"), axis.text = element_text(size = 20, family = "Arial"))
dev.off()



res<-cor.circular(wGribDat$WindHead, wGribDat$estHead, test = T)
res<-cor.circular(atan2(wGribDat$X,wGribDat$Y), atan2(wGribDat$U,wGribDat$V),test=T)
ggplot(sel, aes(x = WindHead, y = EHead)) +
    geom_point() #+
    # geom_line(aes(x=-3:3,y=-3:3))
# res
spres <- cor.test(wGribDat$ESpd, wGribDat$WSpd)
summary(wGribDat)

# test the trajectory values for wind calculations
WindDat$trajHead <- NA
WindDat$trajSpd <- NA
for(ind in 1:nrow(WindDat)){
  trajs <- tryCatch({TrackTraj(WindDat$DT[ind], WindDat$Lat[ind], WindDat$Lon[ind], 6)
    }, error = function(e){c(NA)})
  if(!is.na(trajs)){
    trCord.dec <- SpatialPoints(cbind(rev(trajs$lon), rev(trajs$lat)), proj4string=CRS("+proj=longlat"))
    trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    trutmDiff <- cbind(diff(trCord.utm$coords.x2), diff(trCord.utm$coords.x1))
    WindDat$trajHead[ind] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
    WindDat$trajSpd[ind] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)
  } else {
    WindDat$trajHead[ind] <- NA
    WindDat$trajSpd[ind] <- NA
  }
}
colnames(WindDat)
ggplot(WindDat,aes(x=trajHead,y=WHd)) + geom_point()

Less1<- ggplot(WindDat[WindDat$distTo < 1,], aes(x = RelHead)) + geom_histogram(bins=50,aes(fill = RelHead)) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
One2<- ggplot(WindDat[WindDat$distTo < 2 & WindDat$distTo >= 1,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 12, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Two3<- ggplot(WindDat[WindDat$distTo < 3 & WindDat$distTo >= 2,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Three4<- ggplot(WindDat[WindDat$distTo < 4 & WindDat$distTo >= 3,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Four5<- ggplot(WindDat[WindDat$distTo < 5 & WindDat$distTo >= 4,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Five6<- ggplot(WindDat[WindDat$distTo < 6 & WindDat$distTo >= 5,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Six7<- ggplot(WindDat[WindDat$distTo < 7 & WindDat$distTo >= 6,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Seven8<- ggplot(WindDat[WindDat$distTo < 8 & WindDat$distTo >= 7,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Eight9<- ggplot(WindDat[WindDat$distTo < 9 & WindDat$distTo >= 8,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))
Nine10<- ggplot(WindDat[WindDat$distTo < 10 & WindDat$distTo >= 9,], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_y_continuous(limits=c(0,6))

# png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/1-6km20182019.png", width = 800, height = 900)
ggarrange(Less1,One2,Two3,Three4,Four5,Five6, ncol=2,nrow=3, labels=c("<1km","1-2km","2-3km","3-4km","4-5km","5-6km"))
# dev.off()

# png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/4-10km20182019.png", width = 800, height = 900)
ggarrange(Four5,Five6,Six7,Seven8,Eight9,Nine10, ncol=2,nrow=3, labels=c("4-5km","5-6km","6-7km","7-8km","8-9km","9-10km"))
# dev.off()
# ggplot(WindDat[WindDat$distTo > 5 & WindDat$distTo < 100,], aes(x = RelHead)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2)
Zero10 <- ggplot(WindDat[WindDat$distTo >= 0 & WindDat$distTo < 10,], aes(x = RelHead)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Ten20 <- ggplot(WindDat[WindDat$distTo >= 10 & WindDat$distTo < 20,], aes(x = RelHead)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Twenty30 <- ggplot(WindDat[WindDat$distTo >= 20 & WindDat$distTo < 30,], aes(x = RelHead)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Thirty40 <- ggplot(WindDat[WindDat$distTo >= 30 & WindDat$distTo < 40,], aes(x = RelHead)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))

# png("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/0-40km20182019.png", width = 700, height = 700)
ggarrange(Zero10,Ten20,Twenty30,Thirty40,ncol=2,nrow=2,labels=c("0-10km","10-20km","20-30km","30-40km"))
# dev.off()

AllTraj <- allTraj
AllTraj$relH[AllTraj$relH < -pi] <- AllTraj$relH[AllTraj$relH < -pi] + 2*pi
AllTraj$relH[AllTraj$relH > pi] <- AllTraj$relH[AllTraj$relH > pi] - 2*pi
TrLess1<- ggplot(AllTraj[AllTraj$distTo < 1*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrOne2<- ggplot(AllTraj[AllTraj$distTo < 2*10^3 & AllTraj$distTo >= 1*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrTwo3<- ggplot(AllTraj[AllTraj$distTo < 3*10^3 & AllTraj$distTo >= 2*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrThree4<- ggplot(AllTraj[AllTraj$distTo < 4*10^3 & AllTraj$distTo >= 3*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) +  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrFour5<- ggplot(AllTraj[AllTraj$distTo < 5*10^3 & AllTraj$distTo >= 4*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrFive6<- ggplot(AllTraj[AllTraj$distTo < 6*10^3 & AllTraj$distTo >= 5*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrSix7<- ggplot(AllTraj[AllTraj$distTo < 7*10^3 & AllTraj$distTo >= 6*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) +  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrSeven8<- ggplot(AllTraj[AllTraj$distTo < 8*10^3 & AllTraj$distTo >= 7*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrEight9<- ggplot(AllTraj[AllTraj$distTo < 9*10^3 & AllTraj$distTo >= 8*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrNine10<- ggplot(AllTraj[AllTraj$distTo < 10*10^3 & AllTraj$distTo >= 9*10^3 & AllTraj$rtChg < 0,], aes(x = relH)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
# png("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/1-5km.png",width=800,height=800)
ggarrange(TrLess1,TrOne2,TrTwo3,TrThree4,TrFour5,TrFive6, ncol=3,nrow=2, labels=c("<1km","1-2km","2-3km","3-4km","4-5km","5-6km"))
ggarrange(TrFour5,TrFive6,TrSix7,TrSeven8,TrEight9,TrNine10, ncol=3,nrow=2, labels=c("5-4km","5-6km","6-7km","7-8km","8-9km","9-10km"))

ggplot(AllTraj[AllTraj$distTo < 6*10^3 & AllTraj$distTo >= 5*10^3 & AllTraj$rtChg < 0,], aes(x = trjHd)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))


# ggplot(AllTraj[AllTraj$distTo > 5 & AllTraj$distTo < 100,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2)
TrZero10 <- ggplot(AllTraj[AllTraj$distTo >= 0 & AllTraj$distTo < 10,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrTen20 <- ggplot(AllTraj[AllTraj$distTo >= 10 & AllTraj$distTo < 20,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrTwenty30 <- ggplot(AllTraj[AllTraj$distTo >= 20 & AllTraj$distTo < 30,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
TrThirty40 <- ggplot(AllTraj[AllTraj$distTo >= 30 & AllTraj$distTo < 40,], aes(x = relH)) + geom_histogram(bins=100) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))

ggarrange(TrZero10,TrTen20,TrTwenty30,TrThirty40,ncol=2,nrow=2,labels=c("0-10km","10-20km","20-30km","30-40km"))

WindDat$Yr <- format(WindDat$DT, format="%Y")
windTags <- unique(WindDat[c("ID","Yr")])
plots <- vector(mode="list",length=nrow(windTags))
for(b in 1:nrow(windTags)){
  Less1<- ggplot(WindDat[WindDat$distTo < 1 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50,aes(fill = RelHead)) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
One2<- ggplot(WindDat[WindDat$distTo < 2 & WindDat$distTo >= 1 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 12, family = "Arial"))
Two3<- ggplot(WindDat[WindDat$distTo < 3 & WindDat$distTo >= 2 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Three4<- ggplot(WindDat[WindDat$distTo < 4 & WindDat$distTo >= 3 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Four5<- ggplot(WindDat[WindDat$distTo < 5 & WindDat$distTo >= 4 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Five6<- ggplot(WindDat[WindDat$distTo < 6 & WindDat$distTo >= 5 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Six7<- ggplot(WindDat[WindDat$distTo < 7 & WindDat$distTo >= 6 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Seven8<- ggplot(WindDat[WindDat$distTo < 8 & WindDat$distTo >= 7 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Eight9<- ggplot(WindDat[WindDat$distTo < 9 & WindDat$distTo >= 8 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
Nine10<- ggplot(WindDat[WindDat$distTo < 10 & WindDat$distTo >= 9 & WindDat$ID == windTags$ID[b] & WindDat$Yr == windTags$Yr[b],], aes(x = RelHead)) + geom_histogram(bins=50) + coord_polar(start=-2*pi/2) + xlim(c(-pi,pi)) + scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2), labels=c("Head","Side","Tail","Side")) +
  theme_bw() + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))
plots[[b]]<-ggarrange(Less1,One2,Two3,Three4,Four5,Five6,Six7,Seven8,Eight9,Nine10)
}
plots[[2]]
plots[[b]]<-ggarrange(Less1,One2,Two3,Three4,Four5,Five6,Six7,Seven8,Eight9,Nine10)
plots[[b]]<-ggarrange(Less1,One2)
plots[[b]]<-ggarrange(Nine10)
plots[[b]]<-ggarrange(Three4,Four5,Five6,Six7,Seven8,Eight9,Nine10)

plots[[b]]<-ggarrange(Two3,Three4,Four5,Five6,Six7,Seven8,Eight9,Nine10)




plots[[b]]<-ggarrange(Eight9,Nine10)


plots[[b]]<-NA

plots[[b]]<-ggarrange(Four5,Five6,Six7,Seven8,Eight9)

plots[[b]]<-ggarrange(Six7,Seven8,Eight9,Nine10)
plots[[b]]<-ggarrange(Two3,Three4,Four5,Five6,Six7,Eight9,Nine10)
plots[[b]]<-ggarrange(Less1,One2,Two3,Six7,Seven8,Eight9,Nine10,labels=c("<1","1-2","2-3","6-7","7-8","8-9","9-10"))
plots[[b]]<-ggarrange(Less1,Five6,Nine10,labels=c("<1","5-6","9-10"))
plots[[b]]<-ggarrange(Less1,One2,Seven8,Eight9,Nine10,labels=c("<1","1-2","7-8","8-9","9-10"))
plots[[b]]<-ggarrange(Six7,Seven8,Eight9,Nine10,labels=c("6-7","7-8","8-9","9-10"))

b=1
plots[[b]]
b=b+1

for(b in 1:length(plots)){
  png(paste("/Volumes/GoogleDrive/My Drive/PhD/Notes/WindNotes/IndivRelHead/",windTags[b,1],windTags[b,2],".png",sep=""),width=800,height=800)
  print(plots[[b]])
  dev.off()
}

# test the trajectories vs. wind headings
colnames(WindDat)
a<-80
lat <- WindDat$Lat[a]
lon <- WindDat$Lon[a]
days <- format(WindDat$DT[a],"%Y-%m-%d")
time <- format(WindDat$DT[a],"%Y-%m-%d %H:%M")
hr <- format(WindDat$DT[a],"%H")
hrs <- 6

for(b in 1:nrow(WindDat)){
  trajs <- 
    hysplit_trajectory(
      lat = lat,
      lon = lon,
      height = 10,
      duration = hrs,
      days = days,
      daily_hours = hr,
      direction = "forward",
      met_type = "gdas1",
      extended_met = TRUE,
      met_dir = paste(exec_loc, "met", sep = ""),
      exec_dir = paste(exec_loc, "exec", sep = "")
    ) 


  trCord.dec <- SpatialPoints(cbind(rev(trajs$lon), rev(trajs$lat)), proj4string=CRS("+proj=longlat"))
  trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
  trutmDiff <- cbind(diff(trCord.utm$coords.x2), diff(trCord.utm$coords.x1))
  WindDat$trajHead[b] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
  WindDat$trajSpd[b] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)
}
WindDat$TrRelHead <- WindDat$Head-WindDat$trajHead
WindDat$TrRelHead[WindDat$TrRelHead < -pi] <- WindDat$TrRelHead[WindDat$TrRelHead < -pi] + 2*pi
WindDat$TrRelHead[WindDat$TrRelHead > pi] <- WindDat$TrRelHead[WindDat$TrRelHead > pi] - 2*pi
ggplot(WindDat) + 
  geom_point(aes(x = RelHead, y = TrRelHead))

ggplot(WindDat[WindDat$distTo < 10,]) + 
  geom_point(aes(y = distTo, x = (RelHead))) + coord_polar() + xlim(c(-pi,pi))

ggplot(WindDat[WindDat$distTo < 10,]) + 
  geom_histogram(aes(x = (RelHead))) + coord_polar() + xlim(c(-pi,pi))


range(WindDat$RelHead)

hist(WindDat$aligned)

plot(WindDat$RelHead-WindDat$TrRelHead)

hist(WindDat$RelHead - WindDat$TrRelHead)

WindDat$TrRelHead*(180/pi)

WindDat$Head[a]
(WindDat$Head[a] - trajHead)*(180/pi)
colnames(WindDat)
trajectory_tbl <- trajectory_plot()

trajectory_model <- create_trajectory_model() %>%
  add_trajectory_params(
    lat = lat,
    lon = lon,
    height = 10,
    duration = hrs,
    days = time,
    daily_hours = hr,
    direction = "backward",
    met_type = "gdas1",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()
  return(trajectory_model %>% get_output_tbl())
a<-1:20
tst <- TrackTraj(WindDat$DT[a],WindDat$Lat[a],WindDat$Lon[a],6)


# MODELLING RELATIVE HEADINGS AS DISCRETE VARIABLE
# convert relative headings into discrete
colnames(WindDat)
summary(WindDat)
# or use Bayesian circular GLM package

circlm <- lm.circular(WindDat$RelHead, WindDat$distTo, type=c("c-l"),)

circglm <- circGLM(RelHead ~ distTo + nxtForDur, data = WindDat[WindDat$distTo < 20,])
circglm2 <- circGLM(RelHead ~ distTo, data = WindDat[WindDat$distTo < 20,])

hist(WindDat$distTo)
dat <- generateCircGLMData()
m <- circGLM(th ~ ., dat)
print(m, type="all")

testGLM <- circglmbayes(RelHead ~ distTo, WindDat)
plot(testGLM)

summary(testGLM)
print(testGLM)

# convert RelHead to single side
WindDat$oneSide <- WindDat$RelHead
WindDat$oneSide[WindDat$oneSide < 0] <- WindDat$oneSide[WindDat$oneSide < 0]*-1

ggplot(WindDat[WindDat$distTo < 40,]) + 
  geom_point(aes(x = distTo, y = oneSide))

hist(WindDat$oneSide)

glmtest <- lmer(oneSide ~ distTo | ID, WindDat)
plot(lmtest)


breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")
WindDat$windCond <- NA
WindDat$windCond[WindDat$aligned > -pi/2 & WindDat$aligned < pi/2] <- "Tail"
WindDat$windCond[WindDat$aligned < -pi/2 & WindDat$aligned > -3*pi/2] <- "Side"
WindDat$windCond[WindDat$aligned > pi/2 & WindDat$aligned < 3*pi/2] <- "Side"
WindDat$windCond[WindDat$aligned > 3*pi/2 | WindDat$aligned < -3*pi/2] <- "Tail"

# CALCULATE FORWARD DISPERSION MODEL FOR EACH FORAGING SPOT
# find foraging beginnings
options(timeout = 800)
allD$forage[is.na(allD$forage)] <- 0
forSt <- which(diff(allD$forage) == 1) + 1
if(allD$forage[1] == 1){
  forSt <- c(1, forSt)
}
forEd <- which(diff(allD$forage) == -1)
if(allD$forage[nrow(allD)] == 1){
  forEd <- c(forEd, nrow(allD))
}

TrackDisp <- function(DT, lat, lon, hrs){
  time <- format(DT - lubridate::hours(9), "%Y-%m-%d %H:%M")
  dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = lat, lon = lon, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm(time),
    release_end = lubridate::ymd_hm(time) + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm(time),
    end_time = lubridate::ymd_hm(time) + lubridate::hours(hrs),
    direction = "forward", 
    met_type = "reanalysis",
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()
  dispersion_tbl <- dispersion_model %>% get_output_tbl()
  hulls <- ddply(dispersion_tbl, "hour", find_hull)
  outpt <- list("DispModel" = dispersion_model,"partDisp" = dispersion_tbl, "partPoly" = hulls)
}

# disp1 <- vector(mode="list",length=length(forSt))
# disp2 <- vector(mode="list",length=length(forSt))
# disp3 <- vector(mode="list",length=length(forSt))
# forInfo <- data.frame("DT"=NA,"StInd"=NA,"Beh5"=NA)
# for(b in 1:length(forSt)){
#     disp1[[b]] <- tryCatch({TrackDisp(allD$DT[forSt[b]] - lubridate::hours(1),allD$lat[forSt[b]],allD$lon[forSt[b]], 3)
#     }, error = function(e){c(NA)})
#     disp2[[b]] <- tryCatch({TrackDisp(allD$DT[forSt[b]] - lubridate::hours(2),allD$lat[forSt[b]],allD$lon[forSt[b]], 3)
#     }, error = function(e){c(NA)})
#     disp3[[b]] <- tryCatch({TrackDisp(allD$DT[forSt[b]] - lubridate::hours(3),allD$lat[forSt[b]],allD$lon[forSt[b]], 3)
#     }, error = function(e){c(NA)})
# }
# bring in dispersal data
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ForageDisps.RData")
# } else {
#     load("F:/UTokyoDrive/PhD/Data/splitr/ForageDisps.RData")
# }
summary(disp1[[1]])
disp1[[1]]$partDisp

# calculate foraging starts/ends
allD$forage[is.na(allD$forage)] <- 0
forSt <- which(diff(allD$forage) == 1) + 1
if(allD$forage[1] == 1){
  forSt <- c(1, forSt)
}
forEd <- which(diff(allD$forage) == -1)
if(allD$forage[nrow(allD)] == 1){
  forEd <- c(forEd, nrow(allD))
}

# FORAGING SPOT DISPERSALS
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ForageDisps.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/ForageDisps.RData")
}
# test the average headings for birds as they approach foraging spots and compare with dispersals
allD$yrID <- paste(allD$Year,allD$tagID,sep="")
# calculate changes in location for northing and easting
allD$Ndiff <- c(NA, diff(allD$UTMN))
allD$Ediff <- c(NA, diff(allD$UTME))
# make sure the first change in location for each tag is NA (because no preceeding data)
for(b in unique(allD$yrID)){
  allD[min(which(allD$yrID == b)),c("Ndiff","Ediff")] <- NA
}
for(b in 1:length(forSt)){
  bf <- min(which(allD$DT >= (allD$DT[forSt[b]] - lubridate::hours(2)) & allD$yrID == allD$yrID[forSt[b]]))
  for(g in bf:forSt[b]){
    # calculate headings
    cbind(atan2(allD$Ndiff[g],allD$Ediff[g]), # bird heading
      atan2(allD$UTMN[forSt[b]] - allD$UTMN[g], allD$UTME[forSt[b]] - allD$UTME[g]), # heading to next foraging point
      allD$spTrav[g], # speed travelled
      allD$tripL[g],
      allD$Sex[g],
      over(SpatialPoints(cbind(allD$lon[g],allD$lat[g]),proj4string=CRS('+proj=longlat')),
        SpatialPoints(cbind(disp2[[b]]$partPoly$lon,disp2[[b]]$partPoly$lat),proj4string=CRS('+proj=longlat'))))
      )
  }
}

plot(allD$lon[g],allD$lat[g],xlim=c(142.5,142.8),ylim=c(39.94,39.96))
plot(SpatialPoints(cbind(disp2[[b]]$partPoly$lon,disp2[[b]]$partPoly$lat),proj4string=CRS('+proj=longlat')),add=T)

plot(disp2[[b]]$partPoly$lon,disp2[[b]]$partPoly$lat)
summary(disp2[[b]])
xlim(142.5,142.8)

plot(allD$lon[allD$DT > (allD$DT[forSt[b]] - lubridate::hours(2)) & allD$DT <= allD$DT[forSt[b]] & allD$yrID == allD$yrID[forSt[b]]],
  allD$lat[allD$DT > (allD$DT[forSt[b]] - lubridate::hours(2)) & allD$DT <= allD$DT[forSt[b]] & allD$yrID == allD$yrID[forSt[b]]])

tgYr <- allD[forSt,c("tagID","Year")]
for(b in 1:length(disp1)){
}

b=430
ggplot() + 
  geom_point(data=disp1[[b]]$partDisp, aes(x = lon, y = lat, colour = as.factor(hour))) +
  geom_path(data=allD[allD$DT > (allD$DT[forSt[b]] - lubridate::minutes(60)) & allD$DT <= allD$DT[forSt[b]]& allD$tagID == tgYr$tagID[b] & allD$Year == tgYr$Year[b],], aes(x = lon, y = lat))

# read in data from www.science.oregonstate.edu/ocean.productivity/standard.product.php
if(Sys.info()['sysname'] == "Darwin"){
    netcdLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/Oceanographic/"
} else {
    netcdLoc <- "F:/UTokyoDrive/PhD/Data/Oceanographic"
}
ncFiles <- list.files(netcdLoc, pattern = "*.hdf")
mcFiles <- ncFiles[!grepl(pattern="*.hdf.gz",list.files(netcdLoc))]
filename <- '/Volumes/GoogleDrive/My Drive/PhD/Data/Oceanographic/test.tif'
gdal_translate(ncFiles[1], dst_dataset = filename)
# Load the Geotiff created into R
r <- raster(filename)

  time <- format(allD$DT[forSt[start]] - lubridate::hours(3) - lubridate::hours(9), "%Y-%m-%d %H:%M")

lengths<-seq(from=0,to=max(WindDat$distTo)-1, by = 1)
lengthsL <- lengths+1 
binDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]]))),
  disp = unlist(lapply(1:length(lengths), function(x) circ.disp(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]])$var)))

binDat <- binDat[!is.na(binDat$disp),]
ggplot(binDat, aes(x = dist, y = disp)) + geom_line() + geom_vline(xintercept=10, linetype="dotted")
wGribDat$EHead <- atan2(wGribDat$X, wGribDat$Y)
wGribDat$WindHead <- atan2(wGribDat$U, wGribDat$V)
circ.cor(wGribDat$EHead, wGribDat$WindHead, test=T)
plot(wGribDat$EHead,wGribDat$WindHead)
lines(-3:3,-3:3)
# RAYLEIGH TEST OF UNIFORMITY (FROM VON MISES DIST)
# calculate means and find prob of uniformity
circ.mean(WindDat$RelHead[WindDat$distTo < 10])
rose.diag(WindDat$RelHead[WindDat$distTo < 10],50,prop=3)
r.test(WindDat$RelHead[WindDat$distTo < 10])
r.test(WindDat$RelHead[WindDat$distTo < 1])

circ.mean(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20])
r.test(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20])
circ.mean(WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30])
circ.mean(WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40])

circ.mean(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))
circ.mean(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))
r.test(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))

allTraj$relH[allTraj$relH < -pi] <- allTraj$relH[allTraj$relH < -pi] + 2*pi
allTraj$relH[allTraj$relH > pi] <- allTraj$relH[allTraj$relH > pi] - 2*pi
ggplot(allTraj[allTraj$distTo < 10*10^3,], aes(x = relH, y = distTo)) + 
 geom_point() + coord_polar()


# MODIFIED RAYLEIGHT TEST - TEST IF DISTRIBUTION MEAN IS WHAT IS EXPECTED

# HODGES-AJNE TEST FOR UNIFORMITY - ALTERNATIVE TO RAYLEIGHS W/OUT ASSUME SAMPLING FROM SPECIFIC DISTRIBUTION

# WILCOXON SIGNED-RANK - TEST SYMMETRY AROUND MEDIAN

# WATSON-WILLIAMS TWO-SAMPLE TEST - TEST IF 2 DISTS ARE DISTINCT (CAN BE EXTENDED TO MULTI-SAMPLE)
# OR WATSON TEST (FOR NON-UNIMODAL DATA)
# ALSO WHEELER-WATSON
watson.two(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2])
watson.wheeler.test(WindDat$RelHead[WindDat$distTo < 2], group = WindDat$distTo[WindDat$distTo < 2] < 1)
watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2]))
watson.two(WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2], WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3])
watson.two(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 4 & WindDat$distTo < 5]) 

watson.two(WindDat$RelHead[WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20])
watson.two(WindDat$RelHead[WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30])
watson.two(WindDat$RelHead[WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40],plot=T)
watson.two(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20], WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30])
watson.two(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40])




watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2]))
watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3]))
watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 3 & WindDat$distTo < 4]))



tst <- watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo < 4 & WindDat$distTo >= 3], WindDat$RelHead[WindDat$distTo >= 4 & WindDat$distTo < 5]))[[4]]
summary(tst)
tst[[4]]
# MANN-WHITNEY TEST FOR ANGULAR DISPERSION


# GENERATE KDE
################################################################################
############################ FORAGING BEHAVIOUR KDE ############################
################################################################################
colnames(ForDat)
colnames(ForDat19)
ForAll <- data.frame(DT = c(ForDat$DT, ForDat19$DT),
    Lat = c(ForDat[,2], ForDat19[,2]),
    Lon = c(ForDat[,3], ForDat19[,3]),
    tkb = c(ForDat[,4], ForDat19[,4]),
    dv = c(ForDat[,5], ForDat19[,5]),
    tagID = c(ForDat[,6], ForDat19[,6]),
    year = c(rep("2018", nrow(ForDat)), rep("2019", nrow(ForDat19))),
    levret = c(ForDat$levRet, ForDat19$levRet19),
    sex = c(ForDat$Sex, ForDat19$Sex),
    fBeh = c(ForDat$fBeh, ForDat19$fBeh),
    tripL = c(ForDat$tripL, ForDat19$tripL))

# create a SpatialPoints dataframe
spFD <- SpatialPointsDataFrame(cbind(ForAll$Lon,ForAll$Lat),proj4string=CRS('+proj=longlat'), data = ForAll[,c(1, 4:ncol(ForAll))])
# convert to UTM
UTMDat <- spTransform(spFD, CRS("+proj=utm +zone=54 +datum=WGS84"))
# split by year
spFD18 <- SpatialPointsDataFrame(cbind(ForDat$Lon,ForDat$Lat),proj4string=CRS('+proj=longlat'), data = ForDat[,c(1, 4:ncol(ForDat))])
UTMDat18 <- spTransform(spFD18, CRS("+proj=utm +zone=54 +datum=WGS84"))
spFD19 <-SpatialPointsDataFrame(cbind(ForDat19$Lon,ForDat19$Lat),proj4string=CRS('+proj=longlat'), data = ForDat19[,c(1, 4:ncol(ForDat19))])
UTMDat19 <- spTransform(spFD19, CRS("+proj=utm +zone=54 +datum=WGS84"))

# LDE based on sex (split by year)
#2018
sud18 <- kernelUD(UTMDat18[,18], h = 30*10^3)
image(sud18)
sud18.names <- names(sud18)
ud18 <- lapply(sud18, function(x) lapply(c(25,50,75,95), function(y) try(getverticeshr(x, y))))
sapply(1:length(ud18), function(i) {
    row.names(ud18[[i]]) <<- sud18.names[i]
})
sudF_poly18 <- Reduce(rbind,ud18[[1]])
# transform back to longlat
sudF_pLonLat18 <- spTransform(sudF_poly18, CRS("+proj=longlat +datum=WGS84"))
plot(sudF_poly18)
plot(sudF_pLonLat18)
sudM_poly18 <- Reduce(rbind,ud18[[2]])
sudM_pLonLat18 <- spTransform(sudM_poly18, CRS("+proj=longlat +datum=WGS84"))
plot(sudM_poly18)
sudFdf18 <- fortify(sudF_pLonLat18)
sudMdf18 <- fortify(sudM_pLonLat18)

# test for each foraging trip
# extract 8 day values across foraging trip dates

range(D18$DT)
for(g in 1:length(days)){
  isol = ListD[[1]][format(ListD[[1]]$DT,"%Y-%m-%d") == days[g],]
}
xcoord <- c(round(min(D18$Lon)/0.025)*0.025, round(max(D18$Lon)/0.025)*0.025)
ycoord <- c(round(min(D18$Lat)/0.025)*0.025, round(max(D18$Lat)/0.025)*0.025)
trange <- format(range(D18$DT),"%Y-%m-%d")

ext <- griddap('erdMBchla8day',
    time = trange, latitude = ycoord, longitude = xcoord)
chlorDat <- ext[[2]]
chlorMn <- chlorDat %>% dplyr::group_by(lat,lon) %>% dplyr::summarise(mean=mean(chlorophyll,na.rm=T))
ext19 <- griddap('erdMBchla8day',
  time=format(range(D19$DT),"%Y-%m-%d"), latitude = c(round(min(D19$Lat)/0.025)*0.025, round(max(D19$Lat)/0.025)*0.025),
  longitude=c(round(min(D19$Lon)/0.025)*0.025, round(max(D19$Lon)/0.025)*0.025))
chlorDat19 <- ext19[[2]]
chlorMn19 <- chlorDat19 %>% dplyr::group_by(lat,lon) %>% dplyr::summarise(mean=mean(chlorophyll,na.rm=T))
mycolor <- colors$temperature

ggplot(chlorMn19,aes(x = lon,y=lat,fill=log(mean))) +
  geom_raster(interpolate = F) +
  scale_fill_gradientn(colours = mycolor, na.value = NA) +
  theme_bw() + ylab("latitude") + xlab("longitude")

chlorSp <- SpatialPointsDataFrame(cbind(chlorMn$lon,chlorMn$lat), chlorMn, proj4string=CRS("+proj=longlat"))
D18sp <- SpatialPoints(cbind(D18$Lon, D18$Lat), proj4string=CRS("+proj=longlat"))
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
round(D18$Lat[1]/0.025)*0.025
rerddap::info('erdMBchla8day')

library(data.table)
D18ChlorLoc <- D18[c("Lat","Lon","Forage","Sex","tripL")]
D18ChlorLoc$latSeq <- round(D18ChlorLoc$Lat/0.025)*0.025
D18ChlorLoc$lonSeq <- round(D18ChlorLoc$Lon/0.025)*0.025
D18ChlorLoc$chlorophyll <- NA
for(b in 1:nrow(D18ChlorLoc)){
  D18ChlorLoc$chlorophyll[b] <- chlorMn$mean[dplyr::near(chlorMn$lat,D18ChlorLoc$latSeq[b]) & dplyr::near(chlorMn$lon,D18ChlorLoc$lonSeq[b])]
}

D19ChlorLoc <- D19[c("Lat","Lon","Forage","Sex","tripL")]
D19ChlorLoc$latSeq <- round(D19ChlorLoc$Lat/0.025)*0.025
D19ChlorLoc$lonSeq <- round(D19ChlorLoc$Lon/0.025)*0.025
D19ChlorLoc$chlorophyll <- NA
for(b in 1:nrow(D19ChlorLoc)){
  D19ChlorLoc$chlorophyll[b] <- chlorMn19$mean[dplyr::near(chlorMn19$lat,D19ChlorLoc$latSeq[b]) & dplyr::near(chlorMn19$lon,D19ChlorLoc$lonSeq[b])]
}

ggplot(D18ChlorLoc) + geom_histogram(aes(x = log(chlorophyll), fill = as.factor(Forage)),position="dodge")

ggplot(data = ext$data, aes(x = lon, y = lat, fill = log(chlorophyll))) +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradientn(colours = mycolor, na.value = NA) +
    theme_bw() + ylab("latitude") + xlab("longitude") + 
    geom_path(data = D18,aes(x = Lon,y=Lat))

dataInfo <- rerddap::info("erdMWchla1day")


## FIND DISTANCE FROM PRECEEDING FORAGING POINT
summary(WindDat)
