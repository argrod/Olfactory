library(rerddap)
library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
library(plotdap)
options(timeout = 200)
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
    met_dir = paste(exec_loc, "met", sep = ""),
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
    met_dir = paste(exec_loc, "met", sep = ""),
    exec_dir = paste(exec_loc, "exec", sep = "")) %>%
  run_model()
  return(trajectory_model %>% get_output_tbl())
}

min(D18$Lat)
max(D18$Lat)
min(D18$Lon)
max(D18$Lon)

StLCalc <- function(disp){
  # find the changes in sign of the data
  chgs <- diff(disp)
  cross <- which(diff(chgs))
}

# calculate step lengths as per Humphries et al. 2013
sel <- as.data.frame(allD[allD$tagID == allD$tagID[1] & allD$Year == allD$Year[1],])
sel$deltaT <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel)-1)], units = "secs"))
# point out where the time difference exceeds the median value
cutoff <- median(sel$deltaT, na.rm = T) + 10
stepAreas <- c(1,which(sel$deltaT > cutoff),length(sel))
stepLsN <- vector(mode = "list", length = length(stepAreas))
stepLsE <- vector(mode = "list", length = length(stepAreas))



for(b in 1:(length(stepAreas)-1)){
diff(sel$UTMN[stepAreas[b]:stepAreas[b+1]]) > 0


  isol <- which(diff(sel$UTMN[stepAreas[b]:stepAreas[b + 1]]))
}

ggplot() + 
  geom_line(aes(x = sel$DT[2:nrow(sel)], y = diff(sel$UTMN)))


# calculate high concentration chloro-A locations and test average air trajectories vs bird headings

dInfo <- info("erdMBchla3day")

res <- griddap('erdMBchla3day',time = c('2018-08-25','2018-08-31'),latitude = c(min(D18$Lat), max(D18$Lat)),
    longitude = c(min(D18$Lon), max(D18$Lon)))
myFunc <- function(x) log(x)

summary(res[[2]]$chlorophyll)


days <- unique(res[[2]]$time)
lapply(days, function(x) which(res[[2]]$chlorophyll[res[[2]]$time == x] == max(res[[2]]$chlorophyll[res[[2]]$time == x], na.rm = T)))
res[[2]][39533,c(2,3)]
summary(res[[2]])

ggplot(res[[2]][res[[2]]$time == days[1],], aes(x = lon, y = lat, fill = chlorophyll)) +
    geom_tile() + scale_fill_gradient(trans = 'log')
