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
    # remove non-calculated values
    yoneDat[[b]] <- yoneDat[[b]][yoneDat[[b]]$wDir != 0,]
    gotoDat[[b]] <- read.delim(paste(gotoloc,gotofiles[b],sep = ""), sep = ",", header = F)
    colnames(gotoDat[[b]]) <- c("DT","Lat","Lon","Head","X","Y")
    gotoDat[[b]]$DT <- as.POSIXct(gotoDat[[b]]$DT, tz = "")
}

overDF <- data.frame(DT=POSIXct(),gotoHead = numeric(), gotoX = numeric(), gotoY = numeric(), yoneY = numeric(), yoneX = numeric(),yoneResN = numeric())
overlDat <- vector(mode="list",length=length(yoneDat))
yoneDatrem <- yoneDat
for(b in 1:length(yoneDat)){
    overlDat[[b]] <- overDF
    yoneDatrem[[b]] <- yoneDatrem[[b]][yoneDatrem[[b]]$aveDir != 0,]
    for(g in 1:nrow(gotoDat[[b]])){
        if(any(which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(1)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(1))))){
            inds <- which(yoneDatrem[[b]]$DT > (gotoDat[[b]]$DT[g] - lubridate::minutes(1)) & yoneDatrem[[b]]$DT < (gotoDat[[b]]$DT[g] + lubridate::minutes(1)))
            overlDat[[b]][g,] <- data.frame(DT=gotoDat[[b]]$DT[g],gotoHead = as.numeric(gotoDat[[b]]$Head[g]),gotoX = as.numeric(gotoDat[[b]]$X[g]),gotoY = as.numeric(gotoDat[[b]]$Y[g]),yoneY = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*sin(yoneDatrem[[b]]$wDir[inds]))),yoneX = as.numeric(mean(yoneDatrem[[b]]$wSp[inds]*cos(yoneDatrem[[b]]$wDir[inds]))),yoneResN = mean(yoneDat[[b]]$Resnorm[inds]))
        } else {
            overlDat[[b]][g,] <- cbind(NA,NA,NA,NA,NA,NA,NA)
        }
    }
}

library(dplyr)
allOverL <- bind_rows(overlDat)

allOverL <- allOverL[!is.na(allOverL$gotoY),]
plot(atan2(as.numeric(allOverL$gotoY),as.numeric(allOverL$gotoX)),atan2(as.numeric(allOverL$yoneY),as.numeric(allOverL$yoneX)))

plot(atan2(as.numeric(allOverL$gotoY),as.numeric(allOverL$gotoX)))

plot(atan2(as.numeric(allOverL$yoneY),as.numeric(allOverL$yoneX)))

library(circular)
cor.circular(atan2(as.numeric(allOverL$gotoY),as.numeric(allOverL$gotoX)),atan2(as.numeric(allOverL$yoneY),as.numeric(allOverL$yoneX)),test=T)

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


##################################################################################################
################################ COMPARE SUBSAMPLED 2014 YONE METHOD #############################
##################################################################################################

library(lubridate)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(circular)
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

combD <- vector(mode = "list",length=length(c(onefiles14,onefiles16)))
for(b in 1:length(onefiles14)){
    ODat<-read.delim(paste(oneloc14,onefiles14[b],sep = ""),sep = ",",header=T)
    FDat<-read.delim(paste(fiveloc14,fivefiles14[b],sep = ""),sep = ",",header=T)
    ODat$DT <- as.POSIXct(ODat$time,format="%d-%b-%Y %H:%M:%S")
    FDat$DT <- as.POSIXct(FDat$timeSub,format="%d-%b-%Y %H:%M:%S")
    # remove NaN rows
    if(any(!is.na(ODat$wDir))){
    ODat <- ODat[!is.na(ODat$wDir),]
    FDat <- FDat[!is.na(FDat$wDir),]
    # go through each FDat row and find a nearby (in time) ODat row
    combD[[b]] <- data.frame(OwSpd = double(), FwSpd = double(), OwDir = double(), FwDir = double(), minDiff = double(), DT = double(), Oresnorm = double(),Fresnorm = double())
    for(g in 1:nrow(FDat)){
        minDiff <- min(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))))
        ind <- which(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))) == minDiff)
        combD[[b]][g,] <- cbind(ODat$wSpd[ind],FDat$wSpd[g],ODat$wDir[ind], FDat$wDir[g],minDiff,FDat$DT[g],ODat$Resnorm[ind],FDat$ResnSub[g])
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
    combD[[b + length(onefiles14)]] <- data.frame(OwSpd = double(), FwSpd = double(), OwDir = double(), FwDir = double(), minDiff = double(), DT = double(), Oresnorm = double(),Fresnorm = double())
    for(g in 1:nrow(FDat)){
        minDiff <- min(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))))
        ind <- which(abs(as.numeric(difftime(FDat$DT[g],ODat$DT, units = "secs"))) == minDiff)
        combD[[b + length(onefiles14)]][g,] <- cbind(ODat$wSpd[ind],FDat$wSpd[g],ODat$wDir[ind], FDat$wDir[g],minDiff,FDat$DT[g],ODat$Resnorm[ind],FDat$ResnSub[g])
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
    annotate("text",x=-2.,y=0.2,label="corr = 0.954 \np = 0") +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Subsampled wind heading (rad)') + scale_x_continuous(name='Original data wind headings (rad)')

splm <- lm(OwSpd ~ FwSpd, data = allSubD)

spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(allSubD$FwSpd),y=0:max(allSubD$FwSpd)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=4,y=15.5,label="y = 0.614x + 0.9\np = 1.361") +
    annotate("text",x=4,y=13.,label=expression(paste(R^2," = 0.88"))) +
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
ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-.25,vjust=2)
ggsave(paste(figLoc,"SubSampVal.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 5, units = "in")
 




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