# wind and foraging calculations

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
# install.packages("bpnreg")
library(bpnreg)
options(timeout = 800)

##############################################################################
######################### BRING IN THE FORAGING DATA #########################
##############################################################################

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
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")

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
        WindDat$yrID <- paste("2018_",sub("\\_S.*","",windFiles[b]),sep="")        
    } else {
        toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = F)
        toAdd$yrID <- paste("2018_",sub("\\_S.*","",windFiles[b]),sep="")
        WindDat <- rbind(WindDat, toAdd)
    }
}

# repeat for 2019
if(Sys.info()['sysname'] == "Darwin"){
    windLoc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
} else {
    windLoc <- 'E:/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/'
}
windFiles <- dir(windLoc)

for(b in 1:length(windFiles)){
    toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = F)
    toAdd$yrID <- paste("2019_",sub("\\_S.*","",windFiles[b]),sep="")
    WindDat <- rbind(WindDat, toAdd)
}
colnames(WindDat) <- c("DT","lat","lon","head","X","Y","yrID")
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y-%m-%d %H:%M:%OS", tz = "")

WindDat$WHead <- atan2(WindDat$Y, WindDat$X)
WindDat$WSpeed <- sqrt(WindDat$X^2 + WindDat$Y^2)
Wind.dec <- SpatialPoints(cbind(WindDat$lon,WindDat$lat), proj4string = CRS("+proj=longlat"))
UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
WindDat$UTME <- UTMdat$coords.x1
WindDat$UTMN <- UTMdat$coords.x2

WindDat$timeTo <- NA
WindDat$distTo <- NA
allD$forage[is.na(allD$forage)] <- 0
tags <- unique(WindDat$ID)
allD$yrID <- paste(format(allD$DT,'%Y'),"_",sub("\\_S.*","",allD$tagID),sep="")
for(b in 1:nrow(WindDat)){ # for each row in WindDat
  if(any(allD$forage[allD$yrID == WindDat$yrID[b] & allD$DT > WindDat$DT[b]] == 1)){ # if there are any foraging points after timepoint
    point <- min(which(allD$DT > WindDat$DT[b] & allD$forage == 1 & allD$yrID == WindDat$yrID[b]))
    
    # point <- which(allD$yrID == WindDat$yrID[b]) # find the timepoint of the next foraging points

    # forInd <- min(which(allD$forage[point:max(which(allD$tagID == WindDat$ID[b] & allD$Year == format(WindDat$DT[b], "%Y")))] == 1)) + point - 1 # select the minimum (i.e. the next one)
    
    # WindDat$timeTo[b] <- as.numeric(difftime(allD$DT[point+forInd],WindDat$DT[b], units="secs"))
    # WindDat$distTo[b] <- sqrt((allD$UTMN[forInd] - allD$UTMN[point])^2 + (allD$UTME[forInd] - allD$UTME[point])^2)*10^-3
    WindDat$timeTo[b] <- as.numeric(difftime(allD$DT[point],WindDat$DT[b], units="secs"))
    WindDat$distTo[b] <- sqrt((allD$UTMN[point] - WindDat$UTMN[b])^2 + (allD$UTME[point] - WindDat$UTME[b])^2)*10^-3
  }
}
WindDat$spTrav <- NA
for(b in 1:nrow(WindDat)){
  inds <- which(allD$yrID == WindDat$yrID[b] & allD$DT > (WindDat$DT[b] - lubridate::seconds(5)) & allD$DT < (WindDat$DT[b] + lubridate::seconds(5)))
  WindDat$spTrav[b] <- mean(allD$spTrav[inds])
}
WindDat$tripL <- NA
for(b in 1:nrow(WindDat)){
    ind = which(abs(allD$DT[allD$yrID == WindDat$yrID[b]] - WindDat$DT[b]) == min(abs(allD$DT[allD$yrID == WindDat$yrID[b]] - WindDat$DT[b])))
    WindDat$tripL[b] <- allD$tripL[which(allD$yrID == WindDat$yrID[b])[ind]]
}
# add foraging number for each indiv
WindDat$forNo <- NA
# add a column of transitions to foraging (foraging starts)
forChg = diff(allD$forage) == 1
if(allD$forage[1] == 1){
    forChg <- c(TRUE,forChg)
} else {
    forChg <- c(FALSE,forChg)
}
allD$forChg <- forChg
# add up all previous transitions to foraging, therefore giving the foraging number
for(b in 1:nrow(WindDat)){
    # find next foraging point
    if(any(allD$yrID == WindDat$yrID[b] & allD$DT > WindDat$DT[b] & allD$forage == 1)){
        nxtForInd <- min(which(allD$yrID == WindDat$yrID[b] & allD$DT > WindDat$DT[b] & allD$forage == 1),na.rm=T)
        # sum(allD$forage[1:nxtForInd][allD$yrID == WindDat$yrID[b]])
        WindDat$forNo[b] <- sum(allD$forChg[allD$yrID == WindDat$yrID[b] & allD$DT <= allD$DT[nxtForInd]])
    }
}

# add straightness values for each WindDat row (5 mins either side)
straightness <- function(wdt,wID,window,tType="mins") {
    return(distHaversine(cbind(allD$lon[min(which(allD$DT >= (wdt - as.difftime(window/2,units=tType)) & allD$yrID == wID))],
        allD$lat[min(which(allD$DT >= (wdt - as.difftime(window/2,units=tType)) & allD$yrID == wID))]),
        cbind(allD$lon[max(which(allD$DT <= (wdt + as.difftime(window/2,units=tType)) & allD$yrID == wID))],
        allD$lat[max(which(allD$DT <= (wdt + as.difftime(window/2,units=tType)) & allD$yrID == wID))]))/
    sum(distHaversine(cbind(allD$lon[allD$DT >= (wdt - as.difftime(window/2,units=tType)) & allD$DT <= (wdt + as.difftime(window/2,units=tType)) & allD$yrID == wID],
        allD$lat[allD$DT >= (wdt - as.difftime(window/2,units=tType)) & allD$DT <= (wdt + as.difftime(window/2,units=tType)) & allD$yrID == wID]))))
}

# NO IDEA WHY THIS ISN'T WORKING WITH MAPPLY OR APPLY. NEEDS WORK AS SLOW AS HELL
tst <- NA
for(b in 11757:nrow(WindDat)){
    tst[b] <- straightness(WindDat$DT[b],WindDat$yrID[b],10)
}
WindDat$straightness <- c(tst,rep(NA,nrow(WindDat)-length(tst)))

if(Sys.info()['sysname'] == "Darwin"){
    save(WindDat,file='/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalculations1819.RData')
} else {
    save(WindDat,file='E:/My Drive/PhD/Data/WindCalculations1819.RData')
}


WindDat$RelHead <- WindDat$head-WindDat$WHead
WindDat$RelHead[WindDat$RelHead < -pi] <- WindDat$RelHead[WindDat$RelHead < -pi] + 2*pi
WindDat$RelHead[WindDat$RelHead > pi] <- WindDat$RelHead[WindDat$RelHead > pi] - 2*pi
WindDat$aligned <- WindDat$RelHead + pi
WindDat$aligned[WindDat$aligned > pi] <- WindDat$aligned[WindDat$aligned > pi] - 2*pi


# TEST THE OFFSET FROM HEADWINDS (PI/2 RWH)
WindDat$HdOff <- abs(pi/2 - WindDat$RelHead)
# ggplot(WindDat[WindDat$distTo < 90,]) +
    # geom_point(aes(x = HdOff, y = distTo)) +
    # coord_polar()
WindDat$HdOff[WindDat$HdOff > pi] <- WindDat$HdOff[WindDat$HdOff > pi] - pi
plot(WindDat$distTo,WindDat$HdOff)

# READ IN RAW GPS DATA

if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2018Shearwater/AxyTrek/"
} else {
    fileloc <- "E:/My Drive/PhD/Data/2018Shearwater/AxyTrek/"
}

tags <- dir(fileloc)
GPSDat <- data.frame("GPS"=double(),"Lat"=numeric(),"Lon"=numeric(),"yrID"=character())
for(tag in tags){
    tmp <-  read.table(paste(fileloc,tag,'/',dir(paste(fileloc,tag,sep=""),pattern="*.txt"),sep=""),
        sep='\t',header=FALSE,colClasses=c("character","numeric","numeric",rep("NULL",5)))
    colnames(tmp) <- c("DT","Lat","Lon")
    tmp$yrID <- paste("2018_",sub("*_S.*","",dir(paste(fileloc,tag,sep=""),pattern="*.txt")),sep="")
    tmp$DT <- as.POSIXct(tmp$DT, format = "%d/%m/%Y,%H:%M:%S") + (9*3600)
    GPSDat <- rbind(GPSDat,tmp)
    rm(tmp)
}
# include 2019 data
if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/AxyTrek/"
} else {
    fileloc <- "E:/My Drive/PhD/Data/2019Shearwater/AxyTrek/"
}
tags <- dir(fileloc)
for(tag in tags){
    tmp <-  read.table(paste(fileloc,tag,'/',dir(paste(fileloc,tag,sep=""),pattern="*.txt"),sep=""),
        sep='\t',header=FALSE,colClasses=c("character","numeric","numeric",rep("NULL",5)))
    colnames(tmp) <- c("DT","Lat","Lon")
    tmp$yrID <- paste("2019_",sub("*_S.*","",dir(paste(fileloc,tag,sep=""),pattern="*.txt")),sep="")
    tmp$DT <- as.POSIXct(tmp$DT, format = "%Y/%m/%d,%H:%M:%S") + (9*3600)
    GPSDat <- rbind(GPSDat,tmp)
    rm(tmp)
}
# add UTM values and format datetime
GPS.dec <- SpatialPoints(cbind(GPSDat$Lon,GPSDat$Lat), proj4string = CRS("+proj=longlat"))
UTMdat <- spTransform(GPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
GPSDat$UTME <- UTMdat$coords.x1
GPSDat$UTMN <- UTMdat$coords.x2

# calculate the immediate bird heading for each WindDat value
WindDat$bHead <- NA
for(b in 1:nrow(WindDat)){
    # pointInd <- max(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT <= WindDat$DT[b]),na.rm=T)
    # pointAfter <- min(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT > WindDat$DT[b]),na.rm=T)
    WindDat$bHead[b] <- atan2(GPSDat$UTMN[min(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT > WindDat$DT[b]),na.rm=T)] - GPSDat$UTMN[max(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT <= WindDat$DT[b]),na.rm=T)],
        GPSDat$UTME[min(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT > WindDat$DT[b]),na.rm=T)] - GPSDat$UTME[max(which(GPSDat$yrID == WindDat$yrID[b] & GPSDat$DT <= WindDat$DT[b]),na.rm=T)])
}

# WindDat = na.omit(WindDat)
WindDat$rwh <- WindDat$bHead - WindDat$WHead
WindDat$rwh[WindDat$rwh < -pi] <- WindDat$rwh[WindDat$rwh < -pi] + 2*pi
WindDat$rwh[WindDat$rwh > pi] <- WindDat$rwh[WindDat$rwh > pi] - 2*pi
WindDat$rAligned <- WindDat$rwh + pi
WindDat$rAligned[WindDat$rAligned > pi] <- WindDat$rAligned[WindDat$rAligned > pi] - 2*pi

plot(WindDat$RelHead,WindDat$rwh)

# RETRY DISTRELHEAD WITH NEW BIRD HEADINGS
breaks<-seq(from=0,to=round_any(max(WindDat$distTo,na.rm=T),10,f=ceiling),by=10)
WindDat$bin10 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(WindDat, "bin10", summarise, grp.mean=mean(rAligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- WindDat[WindDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90,], aes(x = rAligned, colour = bin10)) +#max(WindDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.5,show.legend=FALSE)+stat_density(aes(x=rAligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_manual(name="Distance to next \nforaging spot (km)", values = rev(brewer.pal(9,"Blues")),
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(WindDat$bin10[WindDat$distTo < 90])))),", (", as.character(unlist(bin10ns[,2])),")",sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="Proportion across all birds (%)", breaks=seq(0,0.6,.1),labels=seq(0,60,10))


# split into long and short foraging trips
LwDat <- WindDat[WindDat$tripL >= 2,]
breaks<-seq(from=0,to=round_any(max(LwDat$distTo,na.rm=T),10,f=ceiling),by=10)
LwDat$bin10 <- cut(LwDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(LwDat, "bin10", summarise, grp.mean=mean(rAligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- LwDat[LwDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
lDists <- ggplot(LwDat[LwDat$distTo > 0 &LwDat$distTo < 90,], aes(x = rAligned, colour = bin10)) +#max(LwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.5,show.legend=FALSE)+stat_density(aes(x=rAligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_manual(name="Distance to next \nforaging spot (km)", values = rev(brewer.pal(9,"Blues")),
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(LwDat$bin10[LwDat$distTo < 90])))),", (", as.character(unlist(bin10ns[,2])),")",sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="Proportion across all birds (%)", breaks=seq(0,0.6,.1),labels=seq(0,60,10))

# split into long and short foraging trips
SwDat <- WindDat[WindDat$tripL < 2,]
breaks<-seq(from=0,to=round_any(max(SwDat$distTo,na.rm=T),10,f=ceiling),by=10)
SwDat$bin10 <- cut(SwDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(SwDat, "bin10", summarise, grp.mean=mean(rAligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- SwDat[SwDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
sDists <- ggplot(SwDat[SwDat$distTo > 0 &SwDat$distTo < 90,], aes(x = rAligned, colour = bin10)) +#max(SwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.5,show.legend=FALSE)+stat_density(aes(x=rAligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_manual(name="Distance to next \nforaging spot (km)", values = rev(brewer.pal(9,"Blues")),
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(SwDat$bin10[SwDat$distTo < 90])))),", (", as.character(unlist(bin10ns[,2])),")",sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="Proportion across all birds (%)", breaks=seq(0,0.6,.1),labels=seq(0,60,10))

ggarrange(lDists,sDists,nrow=2)

# circular bayesian mixed effects model
# explanatory variables: distance from (and/or time from or interaction term?) + foraging duration?

# example from Cremers et al. 2018
# fit.Maps <- bpnme(pred.I = Error.rad ~ Maze + Trial.type + L.c + (1|Subject), data = Maps, its = 10000, burn = 1000, n.lag = 3, seed = 101)
WindSel <- na.omit(WindDat)
WindSel$yrIDn <- as.numeric(unclass(as.factor(WindSel$yrID)))
WindSel$bin10n <- as.numeric(WindSel$bin10)
WindSel$forNo <- as.numeric(WindSel$forNo)
class(WindSel$bin10n)
fit.Frst <- bpnme(pred.I = rAligned ~ bin10n + (1|yrIDn) + (1|forNo),data=WindSel, its = 10000, burn = 1000, n.lag=3, seed = 101)
colnames(WindSel)
class(WindSel$bin10)
