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
    yoneloc <- "F:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
    gotoloc <- "F:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/MinDat/"
}

yonefiles <- list.files(yoneloc, pattern = '*.txt')
gotofiles <- list.files(gotoloc, pattern = '*.csv')

yoneDat <- vector(mode = "list", length = length(yonefiles))
gotoDat <- vector(mode = "list", length = length(gotofiles))
for(b in 1:length(yoneDat)){
    yoneDat[[b]] <- read.delim(paste(yoneloc,yonefiles[b],sep = ""),sep = ",", header = T)
    yoneDat[[b]]$DT <- sub(',',' ',yoneDat[[b]]$time)
    yoneDat[[b]]$DT <- as.POSIXct(yoneDat[[b]]$DT,tz = "")
    gotoDat[[b]] <- read.delim(paste(gotoloc,gotofiles[b],sep = ""), sep = ",", header = F)
    colnames(gotoDat[[b]]) <- c("DT","Lat","Lon","Head","X","Y")
    gotoDat[[b]]$DT <- as.POSIXct(gotoDat[[b]]$DT, tz = "")
}

overDF <- data.frame(DT=POSIXct(),gotoHead = numeric(), gotoX = numeric(), gotoY = numeric(), yoneY = numeric(), yoneX = numeric())
overlDat <- vector(mode="list",length=length(yoneDat))
yoneDatrem <- yoneDat
for(b in 1:length(yoneDat)){
    overlDat[[b]] <- overDF
    yoneDatrem[[b]] <- yoneDatrem[[b]][yoneDatrem[[b]]$aveDir != 0,]
    for(g in 1:nrow(gotoDat[[b]])){
        if(any(which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(1)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(1))))){
            inds <- which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(1)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(1)))
            overlDat[[b]][g,] <- data.frame(DT=gotoDat[[b]]$DT[g],gotoHead = as.numeric(gotoDat[[b]]$Head[g]),gotoX = as.numeric(gotoDat[[b]]$X[g]),gotoY = as.numeric(gotoDat[[b]]$Y[g]),yoneY = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*sin(yoneDatrem[[b]]$wDir[inds]))),yoneX = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*cos(yoneDatrem[[b]]$wDir[inds]))))
        } else {
            overlDat[[b]][g,] <- cbind(NA,NA,NA,NA,NA,NA)
        }
    }
}

library(dplyr)
allOverL <- bind_rows(overlDat)

allOverL <- allOverL[!is.na(allOverL$gotoY),]
plot(atan2(as.numeric(allOverL$gotoY),as.numeric(allOverL$gotoX)),atan2(as.numeric(allOverL$yoneY),as.numeric(allOverL$yoneX)))

plot(atan2(as.numeric(allOverL$gotoY),as.numeric(allOverL$gotoX)))

plot(atan2(as.numeric(allOverL$yoneY),as.numeric(allOverL$yoneX)))

library(ggplot2)
ggplot(allOverL) + geom_point(aes(y = atan2(as.numeric(yoneY),as.numeric(yoneX)),x = atan2(as.numeric(gotoY),as.numeric(gotoX))),data=allOverL[allOverL$yoneX != 0,]) + scale_y_continuous(name = "Yone method wind direction (rad)") + scale_x_continuous(name = "Goto method wind direction (rad)") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(circular)
cor.circular(atan2(allOverL$gotoY,allOverL$gotoX),atan2(allOverL$yoneY,allOverL$yoneX),test=T)

ggplot(allOverL) + geom_point(aes(y = yoneX,x = gotoX),data=allOverL[allOverL$yoneX != 0,])


sum(allOverL$X == 0)
plot(gotoDat[[b]]$DT,rep(1,nrow(gotoDat[[b]])))
points(yoneDat[[b]]$DT,rep(5,nrow(yoneDat[[b]])))

for(g in 1:nrow(gotoDat[[b]])){
print(min(abs(gotoDat[[b]]$DT[g] - yoneDat[[b]]$DT)))
}
