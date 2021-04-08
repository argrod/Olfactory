library(rerddap)
library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
library(plotdap)
library(scales)
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

StLDisp <- function(utm){
  chgs <- diff(utm) >= 0
  swtch <- which(diff(chgs)!=0) + 1
}
StLCalc <- function(DT,UTMN,UTME){
  deltaT <- c(NA, difftime(DT[2:length(DT)], DT[1:(length(DT)-1)], units = "secs"))
  # find the cutoff points from the time difference
  cutoff <- median(deltaT, na.rm = T) + 10
  # find the start and end times based on time differences over the cutoff value
  strts <- which(deltaT > cutoff)+1
  strts <- c(1,strts[!deltaT[strts] > cutoff])
  ends <- c(which(deltaT > cutoff)-1,length(DT))
  ends <- ends[!deltaT[ends] > cutoff]
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
  data.frame(start=DT[strts], end=DT[ends], dispN = stepLsN, dispE = stepLsE)
}
levy <- function(x,mu){
  x^-mu
}
selTest <- StLCalc(sel$DT, sel$UTMN, sel$UTME)
rank(selTest$dispN)/max(rank(selTest$dispN))
Nls <- data.frame(disp = sort(selTest$dispN, decreasing = T)*10^-3, rank = rev(rank(sort(selTest$dispN, decreasing = T)*10^-3)/max(rank(sort(selTest$dispN, decreasing = T)*10^-3))))
plot(log10(Nls$disp*10), log10(levy(Nls$disp*10^-3,2)))
plot(log10(Nls$disp),log10(levy(sort(selTest$dispN, decreasing = T),2)))
p1<-ggplot(Nls, aes(x = (disp), y = (rank))) +
  geom_point(pch=21) + scale_x_log10() + scale_y_log10()
p1 + geom_line(Nls, mapping=aes(x = (disp), y = (levy(disp,3))), colour = "red")

plot(sort(levy(1:1000,2),decreasing=F),type='l')



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


tstr<-data.frame(dat=c(0,4,5,9,7,2,-4,-6,-2,3,6,5,2,4,6,9,3,0,-4,-8,-2,6,8,10,5,3,1,2,2),swtch=NA,
  dt=c(NA,30,30,30,30,30,51,30,30,30,30,30,30,30,62,182,30,30,30,30,30,30,30,30,30,651,30,30,30))
plot(tstr,type="l")
points(which(diff(diff(tstr$dat) >= 0)!=0) + 1,tstr$dat[which(diff(diff(tstr$dat) >= 0)!=0) + 1])
tstr$swtch<-0
tstr$swtch[which(diff(diff(tstr$dat) >= 0)!=0) + 1]=1
swtch<-which(diff(diff(tstr) >= 0)!=0) + 1
strts <- which(tstr$dt > cutoff)+1
strts <- c(1,strts[!tstr$dt[strts] > cutoff])
ends <- c(which(tstr$dt > cutoff)-1,nrow(tstr))
ends <- ends[!tstr$dt[ends] > cutoff]
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
