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