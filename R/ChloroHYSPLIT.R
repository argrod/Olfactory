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
# install.packages("plotdap",repos='http://cran.us.r-project.org')
library(poweRlaw)
library(ggpubr)
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

min(D18$Lat)
max(D18$Lat)
min(D18$Lon)
max(D18$Lon)

StLDisp <- function(utm){
  chgs <- diff(utm) >= 0
  swtch <- which(diff(chgs)!=0) + 1
}
StLCalc <- function(DT,UTMN,UTME){
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
    tryCatch({
      trajs <- TrackTraj(sel$DT[UDsts$strtInd[ind]], sel$Lat[UDsts$strtInd[ind]], sel$Lon[UDsts$strtInd[ind]], 6)
    }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
    tryCatch({
      trCord.dec <- SpatialPoints(cbind(rev(trajs$lon), rev(trajs$lat)), proj4string=CRS("+proj=longlat"))
      trCord.utm <- spTransform(trCord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
      trutmDiff <- cbind(diff(trCord.utm$coords.x1), diff(trCord.utm$coords.x2))
      trajHead[ind] <- atan2(mean(trutmDiff[,1]), mean(trutmDiff[,2]))
      trajSpd[ind] <- mean(sqrt(trutmDiff[,1]^2 + trutmDiff[,2]^2))/(length(trutmDiff[,1])*3600)
    }, error = function(e){c(NA)})
  }
  outTraj[[b]] <- data.frame(DT = sel$DT[UDsts$strtInd], aveHd = aveHead, trjHd = trajHead, trjSpd = trajSpd, lat = sel$Lat[UDsts$strtInd], lon = sel$Lon[UDsts$strtInd], timeTo = sel$tToFor[UDsts$strtInd], distTo = sel$dToFor[UDsts$strtInd])
}

sel <- allD[allD$forage == 1,]
ggplot(sel, aes(x = lon, y = lat)) + geom_point()
# LOAD IN THE STEP LENGTHS TRAJECTORIES
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/StepsTraj.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/StepsTraj.RData")
}
summary(outTraj)


show <- 3
# ggplot(ListD[[show]]) +
#   geom_path(aes(x = Lon, y = Lat)) +
#   xlim(c(143.5,144.8)) + ylim(c(42,43)) +
#   geom_point(data = outTraj[[show]][UDsts$strtInd,], aes(x = lon, y = lat)) +
#   geom_spoke(data = outTraj[[show]], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")),
#   radius = scales::rescale(outTraj[[show]]$trjSpd, c(.2, .8))) +
#   scale_colour_distiller(palette = "RdYlGn", name = "Approx. speed") +
#   theme(legend.position = 'bottom', legend.direction = 'horizontal') +
#   theme_bw() + theme(panel.grid = element_blank()) +
#   theme(panel.border = element_rect(colour = 'black'))

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


ggplot(outTraj[[show]], aes(x = lon, y = lat)) +
  geom_point() + geom_spoke(data = outTraj[[show]], aes(x = lon, y = lat, colour = trjSpd, angle = trjHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = scales::rescale(outTraj[[show]]$trjSpd, c(.2, .8))) +
  scale_colour_distiller(palette = "RdYlGn", name = "Approx. speed") +
  theme(legend.position = 'bottom', legend.direction = 'horizontal') +
  theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black'))


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
for(b in 1:nrow(WindDat)){
  if(any(allD$forage[allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y") & allD$DT > WindDat$DT[b]] == 1)){
    point <- which(allD$lat == WindDat$Lat[b] & allD$lon == WindDat$Lon[b] & allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y"))
    forInd <- min(which(allD$forage[point:max(which(allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y")))] == 1)) + point - 1
    WindDat$timeTo[b] <- as.numeric(difftime(allD$DT[point+forInd],WindDat$DT[b], units="secs"))
    WindDat$distTo[b] <- sqrt((allD$UTMN[forInd] - allD$UTMN[point])^2 + (allD$UTME[forInd] - allD$UTME[point])^2)*10^-3
  }
}

load("F:/UTokyoDrive/PhD/Data/WindCalc/windDat.RData")
WindDat$WHd <- atan2(WindDat$X,WindDat$Y)
WindDat$RelHead <- WindDat$Head-WindDat$WHd
WindDat$RelHead[WindDat$RelHead < -pi] <- WindDat$RelHead[WindDat$RelHead < -pi] + 2*pi
WindDat$RelHead[WindDat$RelHead > pi] <- WindDat$RelHead[WindDat$RelHead > pi] - 2*pi
ggplot(WindDat, aes(y = distTo, x = Head)) + geom_point() + coord_polar()
ggplot(WindDat, aes(y = distTo, x = WHd)) + geom_point() + coord_polar()







# calculate high concentration chloro-A locations and test average air trajectories vs bird headings

dInfo <- info("erdMBchla3day")

format(range(ListD[[1]]$DT), format="%Y-%m-%d")
res <- griddap('erdMBchla8day',time = format(range(allD$DT), format="%Y-%m-%d"),latitude = c(min(D18$Lat), max(D18$Lat)),
    longitude = c(min(D18$Lon), max(D18$Lon)))
myFunc <- function(x) log(x)

format(allD$DT, format ="%Y-%m-%d" & min(abs(allD$Lat - res[[2]]$lat)) & min(abs(allD$Lon - res[[2]]$lon))



summary(res[[2]]$chlorophyll)


days <- unique(res[[2]]$time)
lapply(days, function(x) which(res[[2]]$chlorophyll[res[[2]]$time == x] == max(res[[2]]$chlorophyll[res[[2]]$time == x], na.rm = T)))
res[[2]][39533,c(2,3)]
summary(res[[2]])

ggplot(res[[2]][res[[2]]$time == days[1],], aes(x = lon, y = lat, fill = chlorophyll)) +
    geom_tile() + scale_fill_gradient(trans = 'log')


hist(log10(res[[2]]$chlorophyll[res[[2]]$time == days[1]]),1000)
