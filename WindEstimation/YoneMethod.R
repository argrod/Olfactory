library(lubridate)
library(caTools)
library(sp)
if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/AxyTrek/"
} else {
    fileloc <- "F:/UTokyoDrive/PhD/Data/2016Shearwater/AxyTrek/"
}

detFl <- function(DT, lat, lon, fs){
    GPS.dec <- SpatialPoints(cbind(lon,lat) , proj4string = CRS("+proj=longlat"))
    UTMDat <- spTransform(GPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    UTME <- UTMDat$lon
    UTMN <- UTMDat$lat
    distTrav <- sqrt(diff(UTME)^2 + diff(UTMN)^2)
    spTrav <- c(NA, distTrav/as.numeric(difftime(DT[2:length(DT)],DT[1:(length(DT)-1)],units = "secs")))
    distTrav <- c(NA, distTrav)
    FLs <- rep(0, length(DT))
    FLs <- runmax(spTrav, ceiling(10/fs), endrule = "NA", align = "left")
    FLs <- FLs > 4
    # remove first and last min
    sts <- which(diff(FLs) == 1) + 1
    eds <- which(diff(FLs) == -1)
    if(sts[1] > eds[1]){
        sts <- c(1, sts)
    }
    if(eds[length(eds)] < sts[length(sts)]){
        eds <- c(eds, length(FLs))
    }
    # remove first and last minute
    for(b in 1:length(sts)){
        FLs[sts[b]:(sts[b] + max(which(DT[sts[b]:length(DT)] <= (DT[sts[b]] + lubridate::minutes(1)))))] <- FALSE
        FLs[min(which(DT[1:eds[b]] >= (DT[eds[b]] - lubridate::minutes(1)))):eds[b]] <- FALSE
    }
    return(data.frame(distTrav = distTrav, spTrav = spTrav, FLs = FLs))
}


files <- list.files(fileloc, pattern = ".txt", recursive = T)
files <- strsplit(files, "/")
tags <- sapply(files, function(x) x[1])
files <- sapply(files, function(x) x[2])
for(b in 1:length(tags)){
    dat <- read.delim(paste(fileloc, tags[b],"/", files[b], sep = ""), sep = "\t", header = F)
    dat <- data.frame(DT=dat[,1], lat = dat[,2], lon = dat[,3])
    dat$DT <- as.POSIXct(dat$DT, format = "%d/%m/%Y %H:%M:%S", tz = "GMT")
    dat$DT <- with_tz(dat$DT, tzone = "Japan")
    dat <- cbind(dat, detFl(dat$DT,dat$lat,dat$lon,1))
    for(b in max(which(dat$DT <= (dat$DT[1] + lubridate::seconds(150)))):(nrow(dat) - 5*60)){
        st <- min(which(dat$DT[1:b] >= dat$DT[b] - lubridate::seconds(15)))

    }
}

###################################################################################################
################################ COMPARE YONE AND GOTO METHOD RESULTS #############################
###################################################################################################

library(lubridate)
if(Sys.info()['sysname'] == "Darwin"){
    yoneloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
    gotoloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
} else {
    yoneloc <- "G:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
    gotoloc <- "G:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/MinDat/"
}

yonefiles <- list.files(yoneloc, pattern = '*.txt')
gotofiles <- list.files(gotoloc, pattern = '*.csv')

yoneDat <- vector(mode = "list", length = length(yonefiles))
gotoDat <- vector(mode = "list", length = length(gotofiles))
for(b in 1:length(yoneDat)){
    yoneDat[[b]] <- read.delim(paste(yoneloc,yonefiles[b],sep = ""),sep = ",", header = T)
    yoneDat[[b]]$DT <- sub(',',' ',yoneDat[[b]]$time)
    yoneDat[[b]]$DT <- as.POSIXct(yoneDat[[b]]$DT,tz = "")
    # remove non-calculated values
    yoneDat[[b]] <- yoneDat[[b]][yoneDat[[b]]$wDir != 0,]
    gotoDat[[b]] <- read.delim(paste(gotoloc,gotofiles[b],sep = ""), sep = ",", header = F)
    colnames(gotoDat[[b]]) <- c("DT","Lat","Lon","Head","X","Y")
    gotoDat[[b]]$DT <- as.POSIXct(gotoDat[[b]]$DT, tz = "")
}

overDF <- data.frame(gotoDT=POSIXct(),yoneDT = POSIXct(),gotoHead = numeric(), gotoX = numeric(), gotoY = numeric(), yoneY = numeric(), yoneX = numeric(),yoneResN = numeric())
overlDat <- vector(mode="list",length=length(yoneDat))
yoneDatrem <- yoneDat
for(b in 1:length(yoneDat)){
    overlDat[[b]] <- overDF
    yoneDatrem[[b]] <- yoneDatrem[[b]][yoneDatrem[[b]]$aveDir != 0,]
    for(g in 1:nrow(gotoDat[[b]])){
        if(any(which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(25)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(25))))){
            inds <- which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(25)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(25)))
            overlDat[[b]][g,] <- data.frame(gotoDT=gotoDat[[b]]$DT[g],
                yoneDT = mean(yoneDat[[b]]$DT[inds]),
                gotoHead = as.numeric(gotoDat[[b]]$Head[g]),
                gotoX = as.numeric(gotoDat[[b]]$X[g]),
                gotoY = as.numeric(gotoDat[[b]]$Y[g]),
                yoneY = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*sin(yoneDatrem[[b]]$wDir[inds]))),
                yoneX = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*cos(yoneDatrem[[b]]$wDir[inds]))),
                yoneResN = mean(yoneDat[[b]]$Resnorm[inds]))
        } else {
            overlDat[[b]][g,] <- cbind(NA,NA,NA,NA,NA,NA,NA,NA)
        }
    }
}

library(dplyr)
allOverL <- bind_rows(overlDat)

allOverL <- allOverL[!is.na(allOverL$gotoY),]
library(circular)
library(ggplot2)
allOverL$gotoDir <- atan2(allOverL$gotoY,allOverL$gotoX)
allOverL$yoneDir <- atan2(allOverL$yoneY,allOverL$yoneX)
allOverL$gotoSpd <- sqrt(allOverL$gotoX^2 + allOverL$gotoY^2)
allOverL$yoneSpd <- sqrt(allOverL$yoneX^2 + allOverL$yoneY^2)

res<-cor.circular(allOverL$yoneDir,allOverL$gotoDir, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
hdval <- ggplot(allOverL, aes(x = gotoDir, y = yoneDir)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=.5,y=3.1,label="corr = 0.27",hjust=0) +
    annotate("text",x=.5,y=2.6,label="p == 0",parse=T,hjust=0) +
    # annotate("text",x=-pi,y=0,label="n = 79",hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

splm <- lm(yoneSpd ~ gotoSpd, data = allOverL)
summary(splm)

ggplotRegressionNoDot <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  stat_smooth(method = "lm", col = "blue")
}

spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(allOverL$gotoSpd),y=0:max(allOverL$gotoSpd)),aes(x=x,y=y),colour="red",linetype='dashed') +
    # geom_point(pch = 21, fill = allOverL$yoneResN) +
    annotate("text",x=7.,y=10.5,label="y = 0.07x + 1.9",hjust=0) +
    annotate("text",x=7.,y=9.5,label="p < 2 %*% 10^{-3}",parse=T,hjust=0) +
    annotate("text",x=7.,y=8.5,label=expression(paste(R^2," = 0.007")),hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name=(("Estimated wind speed (m/s)"))) + scale_x_continuous(name=(("JMA wind speed (m/s)")),limits=c(min(allOverL$gotoSpd),max(allOverL$gotoSpd)))

ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"MethodCompare19.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 3.5, units = "in")


##################################################################################################
################################ COMPARE SUBSAMPLED 2014 YONE METHOD #############################
##################################################################################################

library(lubridate)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(circular)
library(sp)
library(rgdal)
if(Sys.info()['sysname'] == "Darwin"){
    oneloc14 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc14 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/5sFix/"
} else {
    oneloc14 <- "F:/UTokyoDrive/PhD/Data/2014Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc14 <- "F:/UTokyoDrive/PhD/Data/2014Shearwater/WindEst/YoneMet/5sFix/"
}

onefiles14 <- list.files(oneloc14, pattern = '*.txt')
fivefiles14 <- list.files(fiveloc14, pattern = '*.txt')

if(Sys.info()['sysname'] == "Darwin"){
    oneloc16 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc16 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/5sFix/"
} else {
    oneloc16 <- "F:/UTokyoDrive/PhD/Data/2016Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc16 <- "F:/UTokyoDrive/PhD/Data/2016Shearwater/WindEst/YoneMet/5sFix/"
}

onefiles16 <- list.files(oneloc16, pattern = '*.txt')
fivefiles16 <- list.files(fiveloc16, pattern = '*.txt')

if(Sys.info()['sysname'] == "Darwin"){
    oneloc17 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc17 <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/5sFix/"
} else {
    oneloc17 <- "F:/UTokyoDrive/PhD/Data/2017Shearwater/WindEst/YoneMet/1sFix/"
    fiveloc17 <- "F:/UTokyoDrive/PhD/Data/2017Shearwater/WindEst/YoneMet/5sFix/"
}

onefiles17 <- list.files(oneloc17, pattern = '*.txt')
fivefiles17 <- list.files(fiveloc17, pattern = '*.txt')

combD <- vector(mode = "list",length=length(c(onefiles14,onefiles16,onefiles17)))
for(b in 1:length(onefiles14)){
    ODat<-read.delim(paste(oneloc14,onefiles14[b],sep = ""),sep = ",",header=T)
    FDat<-read.delim(paste(fiveloc14,fivefiles14[b],sep = ""),sep = ",",header=T)
    ODat$DT <- as.POSIXct(ODat$time,format="%d-%b-%Y %H:%M:%S")
    FDat$DT <- as.POSIXct(FDat$timeSub,format="%d-%b-%Y %H:%M:%S")
    colnames(ODat) <- c("time","lon","lat","wSpd","wDir","vA","Resnorm","DT")
    colnames(FDat) <- c("timeSub","lonSub","latSub","wSpd","wDir","vA","ResnSub","DT")
    # remove NaN rows
    if(any(!is.na(ODat$wDir))){
        ODat <- ODat[!is.na(ODat$wDir),]
        FDat <- FDat[!is.na(FDat$wDir),]
        # go through each FDat row and find a nearby (in time) ODat row
        combD[[b]] <- data.frame(OwSpd = double(), FwSpd = double(), OwDir = double(), FwDir = double(), minDiff = double(), distDiff = double(), DT = double(), Oresnorm = double(),Fresnorm = double())
        for(g in 1:nrow(FDat)){
            minDiff <- min(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))))
            ind <- which(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))) == minDiff)
            # find the distance between measures (metres)
            OGGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]) , proj4string = CRS("+proj=longlat"))
            UTMDat <- spTransform(OGGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
            OGUTME <- UTMDat$coords.x1
            OGUTMN <- UTMDat$coords.x2
            SubGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]) , proj4string = CRS("+proj=longlat"))
            UTMDat <- spTransform(SubGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
            SubUTME <- UTMDat$coords.x1
            SubUTMN <- UTMDat$coords.x2    
            combD[[b]][g,] <- cbind(ODat$wSpd[ind],FDat$wSpd[g],ODat$wDir[ind], FDat$wDir[g],minDiff,as.numeric(sqrt((OGUTME - SubUTME)^2 + (OGUTMN - SubUTMN)^2)),FDat$DT[g],ODat$Resnorm[ind],FDat$ResnSub[g])
        }
    }
}
# add 2016 data
for(b in 1:length(onefiles16)){
    ODat<-read.delim(paste(oneloc16,onefiles16[b],sep = ""),sep = ",",header=T)
    FDat<-read.delim(paste(fiveloc16,fivefiles16[b],sep = ""),sep = ",",header=T)
    ODat$DT <- as.POSIXct(ODat$time,format="%d-%b-%Y %H:%M:%S")
    FDat$DT <- as.POSIXct(FDat$timeSub,format="%d-%b-%Y %H:%M:%S")
    # remove NaN rows
    ODat <- ODat[!is.na(ODat$wDir),]
    FDat <- FDat[!is.na(FDat$wDir),]
    # go through each FDat row and find a nearby (in time) ODat row
    combD[[b + length(onefiles14)]] <- data.frame(OwSpd = double(), FwSpd = double(), OwDir = double(), FwDir = double(), minDiff = double(), distDiff = double(), DT = double(), Oresnorm = double(),Fresnorm = double())
    for(g in 1:nrow(FDat)){
        minDiff <- min(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))))
        ind <- which(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))) == minDiff)
        OGGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]), proj4string = CRS("+proj=longlat"))
        UTMDat <- spTransform(OGGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
        OGUTME <- UTMDat$coords.x1
        OGUTMN <- UTMDat$coords.x2
        SubGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]) , proj4string = CRS("+proj=longlat"))
        UTMDat <- spTransform(SubGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
        SubUTME <- UTMDat$coords.x1
        SubUTMN <- UTMDat$coords.x2    
        combD[[b + length(onefiles14)]][g,] <- cbind(ODat$wSpd[ind],FDat$wSpd[g],ODat$wDir[ind], FDat$wDir[g],minDiff,as.numeric(sqrt((OGUTME - SubUTME)^2 + (OGUTMN - SubUTMN)^2)),FDat$DT[g],ODat$Resnorm[ind],FDat$ResnSub[g])
    }
}

for(b in 1:length(onefiles17)){
    ODat<-read.delim(paste(oneloc17,onefiles17[b],sep = ""),sep = ",",header=T)
    FDat<-read.delim(paste(fiveloc17,fivefiles17[b],sep = ""),sep = ",",header=T)
    ODat$DT <- as.POSIXct(ODat$time,format="%Y/%m/%d,%H:%M:%S")
    FDat$DT <- as.POSIXct(FDat$timeSub,format="%Y/%m/%d,%H:%M:%S")
    # remove NaN rows
    ODat <- ODat[!is.na(ODat$wDir),]
    FDat <- FDat[!is.na(FDat$wDir),]
    # go through each FDat row and find a nearby (in time) ODat row
    combD[[b + length(onefiles14) + length(onefiles16)]] <- data.frame(OwSpd = double(), FwSpd = double(), OwDir = double(), FwDir = double(), minDiff = double(), distDiff = double(), DT = double(), Oresnorm = double(),Fresnorm = double())
    for(g in 1:nrow(FDat)){
        minDiff <- min(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))))
        ind <- which(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))) == minDiff)
        OGGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]), proj4string = CRS("+proj=longlat"))
        UTMDat <- spTransform(OGGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
        OGUTME <- UTMDat$coords.x1
        OGUTMN <- UTMDat$coords.x2
        SubGPS.dec <- SpatialPoints(cbind(ODat$lon[ind],ODat$lat[ind]) , proj4string = CRS("+proj=longlat"))
        UTMDat <- spTransform(SubGPS.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
        SubUTME <- UTMDat$coords.x1
        SubUTMN <- UTMDat$coords.x2    
        combD[[b + length(onefiles14) + length(onefiles16)]][g,] <- cbind(ODat$wSpd[ind],FDat$wSpd[g],ODat$wDir[ind], FDat$wDir[g],minDiff,as.numeric(sqrt((OGUTME - SubUTME)^2 + (OGUTMN - SubUTMN)^2)),FDat$DT[g],ODat$Resnorm[ind],FDat$ResnSub[g])
    }
}

plot(combD[[1]]$OwSpd,combD[[1]]$FwSpd)
plot(combD[[1]]$OwDir,combD[[1]]$FwDir)

allSubD <- bind_rows(combD)
allSubD <- allSubD[allSubD$minDiff < 10,]

ggplot(allSubD) + geom_point(aes(x = OwSpd, y = FwSpd))

ggplot(allSubD) + geom_point(aes(x = OwDir, y = FwDir)) + scale_x_continuous(limits=c(-pi,pi)) + scale_y_continuous(limits=c(-pi,pi))

ggplot(allSubD) + geom_bar(aes(x = Oresnorm))

cor.circular(allSubD$OwDir, allSubD$FwDir, test = T)

allSubD$OwDir[allSubD$OwDir > pi] <- allSubD$OwDir[allSubD$OwDir > pi] - 2*pi
allSubD$FwDir[allSubD$FwDir > pi] <- allSubD$FwDir[allSubD$FwDir > pi] - 2*pi

ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(pch=21,fill="red") +
  stat_smooth(method = "lm", col = "blue")
}

res<-cor.circular(allSubD$OwDir, allSubD$FwDir, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
hdval <- ggplot(allSubD, aes(x = OwDir, y = FwDir)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=-pi,y=1,label="corr = 0.94",hjust=0) +
    annotate("text",x=-pi,y=0.5,label="p = 0",hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Subsampled wind heading (rad)') + scale_x_continuous(name='Original data wind headings (rad)')

splm <- lm(OwSpd ~ FwSpd, data = allSubD)

spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(allSubD$FwSpd),y=0:max(allSubD$FwSpd)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text", x = 0, y = 16, label = "y == 0.92*x + 0.40",parse = T,hjust=0) +
    annotate("text", x = 0, y = 14.5, label = "p < 2.2 %*% 10^{-16}",parse=T,hjust=0) + 
    annotate("text",x=0,y=13,label=expression(paste(R^2," = 0.90")),hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name=(("Subsampled data wind speed (m/s)"))) + scale_x_continuous(name=(("Original data wind speed (m/s)")))

if(Sys.info()['sysname'] == "Darwin"){
  figLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Figures/Olfactory/"
  # figLoc <- "/Documents/GitHub/PhD/Olfactory/"
} else {
  figLoc <- "F:/UTokyoDrive/PhD/Figures/Olfactory/"
  # figLoc <- "F:/Documents/GitHub/PhD/Olfactory/"
}
ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"SubSampVal.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 3.5, units = "in")





allComb <- bind_rows(combD)
ggplot(allComb) + geom_point(aes(x = as.numeric(OwDir), y = as.numeric(FwDir))) + scale_x_continuous(name = "Original direction (rad)",limits=c(-pi,pi)) + scale_y_continuous(name = "Subsampled direction (rad)",limits=c(-pi,pi)) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank())



subOvDF <- data.frame(DT=POSIXct(), fiveX = numeric(), fiveY = numeric(), oneY = numeric(), oneX = numeric())
subOvlDat <- vector(mode="list",length=length(oneDat))
oneDatrem <- oneDat
for(b in 1:length(oneDat)){
    subOvlDat[[b]] <- subOvDF
    oneDatrem[[b]] <- oneDatrem[[b]][oneDatrem[[b]]$aveDir != 0,]
    for(g in 1:nrow(fiveDat[[b]])){
        if(any(which(oneDatrem[[b]]$DT > (fiveDat[[b]]$DT[g] - lubridate::minutes(1)) & oneDatrem[[b]]$DT < (fiveDat[[b]]$DT[g] + lubridate::minutes(1))))){
            inds <- which(oneDatrem[[b]]$DT > (fiveDat[[b]]$DT[g] - lubridate::minutes(1)) & oneDatrem[[b]]$DT < (fiveDat[[b]]$DT[g] + lubridate::minutes(1)))
            subOvlDat[[b]][g,] <- data.frame(DT=fiveDat[[b]]$DT[g],fiveX = as.numeric(fiveDat[[b]]$X[g]),fiveY = as.numeric(fiveDat[[b]]$Y[g]),oneY = as.numeric(mean(oneDatrem[[b]]$wSp[inds]*sin(oneDatrem[[b]]$wDir[inds]))),oneX = as.numeric(mean(oneDatrem[[b]]$wSp[inds]*cos(oneDatrem[[b]]$wDir[inds]))))
        } else {
            subOvlDat[[b]][g,] <- cbind(NA,NA,NA,NA,NA)
        }
    }
}

library(dplyr)
allOverL <- bind_rows(subOvlDat)

allOverL <- allOverL[!is.na(allOverL$fiveY),]
plot(atan2(as.numeric(allOverL$fiveY),as.numeric(allOverL$fiveX)),atan2(as.numeric(allOverL$oneY),as.numeric(allOverL$oneX)))

plot(atan2(as.numeric(allOverL$fiveY),as.numeric(allOverL$fiveX)))

plot(atan2(as.numeric(allOverL$oneY),as.numeric(allOverL$oneX)))

library(ggplot2)
ggplot(allOverL) + geom_point(aes(y = atan2(as.numeric(oneY),as.numeric(oneX)),x = atan2(as.numeric(fiveY),as.numeric(fiveX))),data=allOverL[allOverL$oneX != 0,]) + scale_y_continuous(name = "one method wind direction (rad)") + scale_x_continuous(name = "five method wind direction (rad)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(circular)
cor.circular(atan2(allOverL$fiveY,allOverL$fiveX),atan2(allOverL$oneY,allOverL$oneX),test=T)

ggplot(allOverL) + geom_point(aes(y = oneX,x = fiveX),data=allOverL[allOverL$oneX != 0,])


sum(allOverL$X == 0)
plot(fiveDat[[b]]$DT,rep(1,nrow(fiveDat[[b]])))
points(oneDat[[b]]$DT,rep(5,nrow(oneDat[[b]])))

for(g in 1:nrow(fiveDat[[b]])){
print(min(abs(fiveDat[[b]]$DT[g] - oneDat[[b]]$DT)))
}


##################################################################################################
#################################### VALIDATION WITH JMA RESULTS #################################
##################################################################################################

library(circular)
library(ggplot2)
library(ggpubr)
# 2014 and 2016
valDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/Comb201420162017Validation.csv",sep=",")
valDat$Time <- sub("T"," ",valDat$Time)
valDat$Time <- as.POSIXct(valDat$Time,format="%Y-%m-%d %H:%M:%S",tz = "")
head(valDat)
cor.circular()
head(valDat)
    
ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(pch=21,fill="red") +
  stat_smooth(method = "lm", col = "blue")
}

res<-cor.circular(valDat$estHead,valDat$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
hdvalSub <- ggplot(valDat, aes(x = gribHead, y = estHead)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=-pi,y=.5,label="corr = 0.31",hjust=0) +
    annotate("text",x=-pi,y=0.,label="p = 0.002",hjust=0) +
    annotate("text",x=-pi,y=-0.5,label="n = 99",hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

splm <- lm(estSpeed ~ gribSpeed, data = valDat)
summary(splm)
spdvalSub <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(valDat$gribSpeed),y=0:max(valDat$gribSpeed)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=min(valDat$gribSpeed),y=12.5,label="y = 0.33x + 2.94",hjust=0) +
    annotate("text",x=min(valDat$gribSpeed),y=11.5,label="p = 0.07",hjust=0) +
    annotate("text",x=min(valDat$gribSpeed),y=10.5,label=expression(paste(R^2," = 0.012")),hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name=(("Estimated wind speed (m/s)"))) + scale_x_continuous(name=(("JMA wind speed (m/s)")),limits=c(min(valDat$gribSpeed),max(valDat$gribSpeed)))

ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"141617windVal.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 3.5, units = "in")

##################################################################################################
#################################### VALIDATION WITH JMA RESULTS #################################
##################################################################################################

# 2019
valDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMethodValidation.csv",sep=",")
valDat$Time <- sub("T"," ",valDat$Time)
valDat$Time <- as.POSIXct(valDat$Time,format="%Y-%m-%d %H:%M:%S",tz = "")
head(valDat)
cor.circular()
head(valDat)

ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(pch=21,fill="red") +
  stat_smooth(method = "lm", col = "blue")
}

res<-cor.circular(valDat$estHead,valDat$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
hdval <- ggplot(valDat, aes(x = gribHead, y = estHead)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=-pi,y=.5,label="corr = 0.19",hjust=0) +
    annotate("text",x=-pi,y=-0.,label="p < 0.01",hjust=0) +
    # annotate("text",x=-pi,y=0,label="n = 79",hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

splm <- lm(estSpeed ~ gribSpeed, data = valDat)
summary(splm)
spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(valDat$gribSpeed),y=0:max(valDat$gribSpeed)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=min(valDat$gribSpeed),y=12.5,label="y = 0.05x + 2",hjust=0) +
    annotate("text",x=min(valDat$gribSpeed),y=11.5,label="p = 0.4",hjust=0) +
    annotate("text",x=min(valDat$gribSpeed),y=10.5,label=expression(paste(R^2," = 0.003")),hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name=(("Estimated wind speed (m/s)"))) + scale_x_continuous(name=(("JMA wind speed (m/s)")),limits=c(min(valDat$gribSpeed),max(valDat$gribSpeed)))

ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"19windVal.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 3.5, units = "in")
nrow(valDat)

#####################################################################################################
################################ TEST METHOD ON ORIGINALLY SAMPLED DATA #############################
#####################################################################################################

OGdat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/OG201420162017Validation.csv", sep = ",", header = T)
OGdat$Time <- sub("T"," ",OGdat$Time)
OGdat$Time <- as.POSIXct(OGdat$Time,format="%Y-%m-%d %H:%M:%S",tz = "")
ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(pch=21,fill="red") +
  stat_smooth(method = "lm", col = "blue")
}

res<-cor.circular(OGdat$estHead,OGdat$gribHead, test = T)
res
splm <- lm(estSpeed ~ gribSpeed, data = OGdat)
summary(splm)

hdvalOG <- ggplot(OGdat, aes(x = gribHead, y = estHead)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=-pi,y=.5,label="corr = 0.21",hjust=0) +
    annotate("text",x=-pi,y=-0.,label="p < 0.03",hjust=0) +
    annotate("text",x=-pi,y=-0.5,label="n = 93",hjust=0) +
    # annotate("text",x=-pi,y=0,label="n = 79",hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

spdvalOG <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(OGdat$gribSpeed),y=0:max(OGdat$gribSpeed)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=min(OGdat$gribSpeed),y=13.5,label="y = 0.23x + 3.6",hjust=0) +
    annotate("text",x=min(OGdat$gribSpeed),y=12.5,label="p = 0.3",hjust=0) +
    annotate("text",x=min(OGdat$gribSpeed),y=11.5,label=expression(paste(R^2," = 0.013")),hjust=0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name=(("Estimated wind speed (m/s)"))) + scale_x_continuous(name=(("JMA wind speed (m/s)")),limits=c(min(OGdat$gribSpeed),max(OGdat$gribSpeed)))

# combo original and subsample
ggarrange(hdvalOG, hdvalSub, spdvalOG, spdvalSub, ncol=2,nrow=2, labels=c("a)","b)","c)","d)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"141617windVal.svg",sep=""), device="svg", dpi = 300, height = 6.5,
      width = 6.5, units = "in")

##########################################################################################################################
########################################### COMPARISON WITH CYBEROCEAN CALCULATIONS ######################################
##########################################################################################################################

# bring in circle calculations
if(Sys.info()['sysname'] == "Darwin"){
    cDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/WindEst/YoneMet/N18031WindYone.txt",sep=",",header=T)
} else {
    cDat <- read.delim("F:/UTokyoDrive/PhD/Data/2018Shearwater/GypsyTranslocation/WindEst/YoneMet/N18031WindYone.txt",sep=",",header=T)
}
cDat$time <- as.POSIXct(cDat$time,format="%d-%b-%Y %H:%M:%S",tz="")
# remove extra rows
cDat <- cDat[!cDat$aveDir == 0,]

# read CyberOcean data
if(Sys.info()['sysname'] == "Darwin"){
    oDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/GypsyTranslocation/CyberOceanwind18031.csv",sep=",",header=T)
} else {
    oDat <- read.delim("F:/UTokyoDrive/PhD/Data/2018Shearwater/GypsyTranslocation/CyberOceanwind18031.csv",sep=",",header=T)
}
oDat$datetime <- as.POSIXct(oDat$datetime, format= "%Y/%m/%d %H:%M:%S")

combD <- data.frame(tdiff = double(), oLat = double(), oLon = double(), cWdir = double(), cWspd = double(), oWdir = double(), oWspd = double())
# go through each row of data and find nearest values
for(b in 1:nrow(oDat)){
    ind <- which(abs(as.numeric(difftime(oDat$datetime[b],cDat$time,units="secs"))) == min(abs(as.numeric(difftime(oDat$datetime[b],cDat$time,units="secs")))))
    tdiff <- min(abs(as.numeric(difftime(oDat$datetime[b],cDat$time,units="secs"))))
    combD[b,] <- data.frame(tdiff, oDat$latitude[b], oDat$longitude[b], cDat$wDir[ind]*(180/pi), cDat$wSp[ind], oDat$wind_direction.deg[b], oDat$wind_speed.m.s[b])
}

plot(combD$oWspd,combD$cWspd)
plot(combD$oWdir,combD$cWdir)
