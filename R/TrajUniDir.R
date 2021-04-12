library(rerddap)
library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
library(plotdap)
library(scales)
install.packages("poweRlaw")
library(poweRlaw)
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


################################################################################################################################
############################################################# FUNCTIONS ########################################################
################################################################################################################################
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
