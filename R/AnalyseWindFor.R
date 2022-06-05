install.packages("Gmisc")
install.packages("sf")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("lubridate")
install.packages("rgdal")
install.packages("adehabitatHR")
install.packages("ggsn")
install.packages("ggspatial")
install.packages("sp")
install.packages("plyr")
install.packages("dplyr")
install.packages("mapdata")
install.packages("rerddap")
install.packages("data.table")
install.packages("ggplot2")
install.packages("viridis")
install.packages("MASS")
install.packages("diagram")
install.packages("ggthemes")
install.packages("extrafont")
install.packages("bpnreg")

library(Gmisc)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(lubridate)
library(rgdal)
library(adehabitatHR)
library(ggsn)
library(ggspatial)
library(sp)
library(plyr)
library(dplyr)
library(mapdata)
library(rerddap)
library(data.table)
library(ggplot2)
library(viridis)
library(MASS)
library(diagram)
library(ggthemes)
library(extrafont)
library(CircMLE)
library(rgeos)
library(bpnreg)
library(circular)
library(geosphere)

#################################################################################
######################## BRING IN THE FORAGING ESTIMATES ########################
#################################################################################

if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/DatEth2018.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/DatEth2019.RData")
    outloc <- "/Volumes/GoogleDrive/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("E:/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("E:/My Drive/PhD/Data/DatEth2018.RData")
    load("E:/My Drive/PhD/Data/DatEth2019.RData")
    outloc <- "E:/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
}
# for(d in Dat){
#     d$distTrav <- c(NA,distHaversine(cbind(d$Lon[1:(nrow(d)-1)],d$Lat[1:(nrow(d)-1)]),
#         cbind(d$Lon[2:nrow(d)],d$Lat[2:nrow(d)])))
#     d$spTrav <- c(NA,(d$distTrav[2:nrow(d)]/1000)/as.numeric(difftime(d$DT[2:nrow(d)],d$DT[1:(nrow(d)-1)],units="mins")))
# }
D18 <- bind_rows(Dat)
# remove values with unrealistic speeds
D18 <- D18[-which(D18$spTrav > 100),]
# for(d in Dat19){
#     d$distTrav <- c(NA,distHaversine(cbind(d$Lon[1:(nrow(d)-1)],d$Lat[1:(nrow(d)-1)]),
#         cbind(d$Lon[2:nrow(d)],d$Lat[2:nrow(d)])))
#     d$spTrav <- c(NA,(d$distTrav[2:nrow(d)]/1000)/as.numeric(difftime(d$DT[2:nrow(d)],d$DT[1:(nrow(d)-1)],units="mins")))
# }
D19 <- bind_rows(Dat19)
D19 <- D19[-which(D19$spTrav > 100),]
allD <- data.frame(DT=c(D18$DT, D19$DT),
    lat = c(D18$Lat, D19$Lat),
    lon = c(D18$Lon, D19$Lon),
    forage = c(D18$Forage,D19$Forage),
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
allD$yrID <- paste(format(allD$DT,"%Y"),sub('\\_S.*','',allD$tagID),sep="_")

###############################################################################
######################## BRING IN THE WIND ESTIMATIONS ########################
###############################################################################

if(Sys.info()['sysname'] == "Darwin"){
    # windLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalculations1819.RData")
} else {
    # windLoc <- 'E:/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/'
    load('E:/My Drive/PhD/Data/WindCalculations1819.RData')
}

################################################################################################################
################################  HR AND RALEIGH TEST P VALUES AND AVE HEADINGS ################################
################################################################################################################

distGaps <- seq(0,90,10)
distGapsL <- distGaps+10
one2Ten <- vector(mode="list", length = length(distGaps))
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}
avRelHd <- NA
pvals<-vector(mode="list",length=length(distGaps))
wDat <- na.omit(WindDat)
for(b in 1:length(distGaps)){
    RaylT <- r.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])))
    tst<-HR_test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])),iter=999)
    pvals[[b]] <- data.frame(Distance=paste0(as.character(distGaps[b]),"-",as.character(distGapsL[b])),RlP = RaylT$p.value,RlR = RaylT$r.bar,HRp = tst[2])
}
save(pvals,file='E:/My Drive/PhD/Data/pvals.RData')

# repeat for 1 km bins from 10 downwards
distGaps <- seq(0,9,1)
distGapsL <- distGaps+1
pvals1km<-vector(mode="list",length=length(distGaps))
for(b in 1:length(distGaps)){
    RaylT <- r.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])))
    tst<-HR_test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])),iter=999)
    pvals1km[[b]] <- data.frame(Distance=paste0(as.character(distGaps[b]),"-",as.character(distGapsL[b])),RlP = RaylT$p.value,RlR = RaylT$r.bar,HRp = tst[2])
}
save(pvals1km,file='E:/My Drive/PhD/Data/pvals1km.RData')

distGaps <- seq(0,100,10)
distGapsL <- distGaps+10
wwLs <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wwLs[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$tripL > 2])),
        na.omit(WindDat$distTo[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$tripL > 2]) > distGapsL[b])
}

wwt <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wwt[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10)])),
        na.omit(WindDat$distTo[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10)]) > distGapsL[b])
}

WindDat$OutHome <- WindDat$OutHm > 0
WindDat$OutHome[WindDat$OutHome == TRUE] <- "Home"
WindDat$OutHome[WindDat$OutHome == FALSE] <- "Out"
wwOut <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wwOut[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$OutHome == "Out"])),
        na.omit(WindDat$distTo[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$OutHome == "Out"]) > distGapsL[b])
}


wwSs <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wwSs[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$tripL <= 2])),
        na.omit(WindDat$distTo[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 10) & WindDat$tripL <= 2]) > distGapsL[b])
}

distGaps <- seq(0,100,10)
distGapsL <- distGaps+10
wws <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wws[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])),
        na.omit(WindDat$tripL[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]] > 2))
}

ggplot(WindDat[WindDat$distTo > 10 & WindDat$distTo <= 40,]) +
    geom_point(aes(x = RelHead, y= distTo, colour = tripL > 2)) + coord_polar()

distGaps <- seq(0,10,1)
distGapsL <- distGaps+1
wws <- vector(mode="list",length=length(distGaps)-1)
for(b in 1:length(wws)){
    wws[[b]] <- watson.wheeler.test(circular(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 1)])),
        na.omit(WindDat$distTo[WindDat$distTo >= distGaps[b] & WindDat$distTo < (distGapsL[b] + 1)]) > distGapsL[b])
}

# # repeat for both long and short trips
# wLong <- wDat[wDat$tripL > 2,]
# wShort <- wDat$wDat$tripL <= 2,]
# pvalsLS<-vector(mode="list",length=length(distGaps))
# for(b in 1:length(distGaps)){
#     RaylTS <- r.test(wShort$rwh[wShort$distTo >= distGaps[b] & wShort$distTo < distGapsL[b]])
#     tstS<-HR_test(wShort$rwh[wShort$distTo >= distGaps[b] & wShort$distTo < distGapsL[b]])
#     RaylTL <- r.test(wLong$rwh[wLong$distTo >= distGaps[b] & wLong$distTo < distGapsL[b]])
#     tstL<-HR_test(wLong$rwh[wLong$distTo >= distGaps[b] & wLong$distTo < distGapsL[b]])
#     pvalsLS[[b]] <- data.frame(Distance=paste0(as.character(distGaps[b]),"-",as.character(distGapsL[b])),SRlP = RaylTS$p.value,SRlR = RaylTS$r.bar,SHRp = tstS[2],
#         LRlP = RaylTL$p.value,LRlR = RaylTL$r.bar,LHRp = tstL[2])
# }
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/pvalsLS.RData")
} else {
    load('E:/My Drive/PhD/Data/pvalsLS.RData')
}
for(df in pvalsLS){
    rownames(df) <- NULL
}
pvalsLS <- bind_rows(pvalsLS)
pvalsLS$LRlP
ggplot(pvalsLS) + 
    geom_line(aes(x = dist, y = LRlR, colour = "deepskyblue")) +
    geom_point(aes(x = dist, y = LRlR, fill = "deepskyblue"), pch=21, alpha = pvalsLS$LRlP < 0.01) + 
    geom_line(aes(x = dist, y = SRlR, colour = "red")) +
    geom_point(aes(x = dist, y = SRlR, fill = "red"),pch = 21, alpha = pvalsLS$SRlP < 0.01)

tstL <- data.frame(mean=sapply(1:10, function(x) mean.circular(circular(wDat$RelHead[wDat$distTo > distGaps[x] & wDat$distTo <= distGapsL[x] & wDat$tripL > 2]))), dev = sapply(1:10, function(x) angular.deviation(circular(wDat$RelHead[wDat$distTo > distGaps[x] & wDat$distTo <= distGapsL[x] & wDat$tripL > 2]))))

tstS <- data.frame(mean=sapply(1:10, function(x) mean.circular(circular(wDat$RelHead[wDat$distTo > distGaps[x] & wDat$distTo <= distGapsL[x] & wDat$tripL <= 2]))),dev=sapply(1:10, function(x) angular.deviation(circular(wDat$RelHead[wDat$distTo > distGaps[x] & wDat$distTo <= distGapsL[x] & wDat$tripL <= 2]))))

tstL$length<-"Long"
tstL$RlR <- pvalsLS$LRlR
tstL$Dist <- distGapsL
tstL$rP <- pvalsLS$LRlP
tstL$hP <- pvalsLS$LHRp
tstS$length<-"Short"
tstS$RlR <- pvalsLS$SRlR
tstS$Dist <- distGapsL
tstS$rP <- pvalsLS$SRlP
tstS$hP <- pvalsLS$SHRp
tstAll <- rbind(tstL,tstS)

plot(as.circular(0, zero = 0))

x <- as.circular(pi, control.circular=list(units="radians", zero=pi))
y <- conversion.circular(circular(pi), zero=pi)
res <- plot(x)
points(y, col=2, plot.info=res)

HR_test(as.circular(wDat$RelHead[wDat$distTo > 10 & wDat$distTo <= 20 & wDat$tripL > 2],units="radians"))


ggplot() + geom_point(data=tstAll,aes(x = mean, y = Dist, fill = RlR, pch = length),size=3) +
    coord_polar(start=pi) + scale_shape_manual(name = "Trip type", values=c(21,24)) +
    scale_fill_distiller(name=expression(bar(italic(r))),palette="YlOrRd") + theme_bw() +
    scale_x_continuous(name = "Relative wind ")

#################################################################################################################
##########################  WATSON-WHEELER TEST FOR 10KM BINS FOR LONG AND SHORT TRIPS###########################
#################################################################################################################

distGaps <- seq(0,90,10)
distGapsL <- distGaps+10
b <- 1

LS10kWW <- lapply(distGaps, function(x) watson.wheeler.test(WindDat$RelHead[WindDat$distTo > x & WindDat$distTo < (x + 10)],
    group = WindDat$tripL[WindDat$distTo > x & WindDat$distTo < (x + 10)]))


watson.wheeler.test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]],
    group = WindDat$tripL[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]] <= 2)

WindDat$DofY <- strftime(WindDat$DT, format= "%j")
WindDat$year <- year(WindDat$DT)
totst <- unique(WindDat[,c("DofY","year")])

wwTests <- lapply(distGaps, function(x) watson.two.test(WindDat$RelHead[WindDat$distTo > x & WindDat$distTo <= x + 10],
    WindDat$RelHead[WindDat$distTo > (x+10) & WindDat$distTo <= (x + 20) & WindDat$tripL <= 2]))


sapply(wwTests, function(x) print(x$statistic))
WWs <- data.frame(dGaps = distGaps, pvals = array(unlist(sapply(wwTests, '[',4))))

# run watson-wheeler tests for 10km bins for outgoing birds as they approach foraging spots
distGaps <- seq(0,90,10)
WindDat$OutHm <- WindDat$OutHm > 0      # OutHm false for outward legs
wwTestsL <- lapply(distGaps, function(x) watson.two.test(WindDat$RelHead[WindDat$distTo > x & WindDat$OutHm == F & WindDat$distTo <= (x + 10) & WindDat$tripL > 2],
    WindDat$RelHead[WindDat$distTo > (x+10) & WindDat$distTo <= (x + 20) & WindDat$OutHm == F & WindDat$tripL > 2]))
wwTestsS <- lapply(distGaps, function(x) watson.two.test(WindDat$RelHead[WindDat$distTo > x & WindDat$OutHm == F & WindDat$distTo <= (x + 10) & WindDat$tripL <= 2],
    WindDat$RelHead[WindDat$distTo > (x+10) & WindDat$distTo <= (x + 20) & WindDat$OutHm == F & WindDat$tripL <= 2]))
# warnings produced are simply about coercing RelHead to circular
WWs <- data.frame(dGaps = distGaps, pvals = array(unlist(sapply(wwTests, '[',1))))

wwTestsAll <- lapply(distGaps, function(x) watson.two.test(WindDat$RelHead[WindDat$distTo > x & WindDat$distTo <= (x + 10)],
    WindDat$RelHead[WindDat$distTo > (x+10) & WindDat$distTo <= (x + 20)]))

wwTestsall <- lapply(distGaps, function(x) watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo > x & WindDat$distTo <= (x + 10)],
    WindDat$RelHead[WindDat$distTo > (x+10) & WindDat$distTo <= (x + 20)])))


testdata = circular::rvonmises(20, mu = circular::circular(pi), kappa = 3)
testdata = na.omit(circular::circular(WindDat$RelHead[WindDat$distTo < 20 & WindDat$distTo >= 10]))
HR_test(testdata, iter = 999)

wwTests[[1]]

tripLengths <- allD %>% dplyr::group_by(yrID,tripL > 2) %>% dplyr::summarise(n = max(distFk))
dailyLengths <- as.data.frame(allD %>% dplyr::group_by(yrID,Day) %>% dplyr::summarise(lengths = sum(distTrav),
    speeds = mean(spTrav)))

mean(tripLengths$n[tripLengths[,2] == F],na.rm=T)
mean(tripLengths$n[tripLengths[,2] == T],na.rm=T)

unique(allD$yrID)
sum(allD$distTrav[allD$yrID == "2019_2018-05" & allD$Day == 4])/1000

mean(dailyLengths$lengths,na.rm=T)
hist(allD$distTrav)

mean(dailyLengths$speeds,na.rm=T)*24


tester <- allD[allD$yrID == unique(allD$yrID)[3] & allD$tripL == 1,]
tester$long <- tester$lon
tester$sumLength <- cumsum(tester$distTrav/1000)
ggplot(tester) +
    geom_point(aes(x = lon, y= lat, colour = sumLength))


ggplot() +
    geom_path(data = tester[tester$sumLength <= 400,], aes(x = lon, y = lat),colour="red")+
    geom_path(data = tester[tester$sumLength <= 300,], aes(x = lon, y = lat),colour="blue")+
    geom_path(data = tester[tester$sumLength <= 200,], aes(x = lon, y = lat),colour="purple")+
    ggsn::scalebar(data = tester[tester$sumLength <= 300,], dist = 50, model = 'WGS84',transform=T,dist_unit="km", height = .05,location="topleft")
    #   st.dist = .1,x.min = 142, x.max = 142.5, y.min = 39.15, y.max = 39.2, location = 'bottomleft',box.fill=c("black","white"))


sum(allD$distTrav[allD$yrID == unique(allD$yrID)[2] & allD$tripL == 1])/1000

watson.wheeler.test(WindDat$RelHead[WindDat$distTo > 40 & WindDat$distTo <= 60],
    group = WindDat$distTo[WindDat$distTo > 40 & WindDat$distTo <= 60] > 50)

watson.two.test(WindDat$RelHead[WindDat$distTo > 40 & WindDat$distTo <= 50],WindDat$RelHead[WindDat$distTo > 50 & WindDat$distTo <= 60])

plot(WindDat$RelHead[WindDat$distTo > 40 & WindDat$distTo <= 50],col="red")
points(WindDat$RelHead[WindDat$distTo > 50 & WindDat$distTo <= 60],col="blue")

b<-6
ggplot(na.omit(WindDat[which(WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]),])) + 
        geom_point(aes(x = RelHead, y = distTo, fill = tripL <= 2), pch = 21) + coord_polar(start=pi) + 
        scale_x_continuous(name = "Relative wind headings (rad)", limits = c(-pi,pi)) +
        scale_fill_manual(name = "Trip length", breaks = c(FALSE,TRUE), labels = c("Long","Short"), values=c("red","deepskyblue"))
for(b in 1:length(distGaps)){
    ggplot(na.omit(WindDat[which(WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]),])) + 
        geom_density(aes(x = RelHead, fill = tripL <= 2), alpha = .3) + coord_polar(start=pi) + 
        scale_x_continuous(name = "Relative wind headings (rad)", limits = c(-pi,pi)) +
        scale_fill_manual(name = "Trip length", breaks = c(FALSE,TRUE), labels = c("Long","Short"), values=c("red","deepskyblue"))
}

#################################################################################################################
############################  FINDING FORAGING WITH WIND CALCULATED BEFORE (30 MINS) ############################
#################################################################################################################

# deal with on tag-by-tag basis
# wind data (WindDat)
# foraging data (Dat[[tag#]])
tg <- 1
DatSel <- Dat[[tg]]
DatSel$forageLong[DatSel$tkb == 1 | DatSel$dv == 1] = 1
DatSel$forageLong[is.na(DatSel$forageLong)] = 0
WindSel <- WindDat[WindDat$ID == DatSel$tagID[tg],]
#find shared latlons of tag and wind data
dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
splits <- which(dt > 70)
splits <- c(0, splits, nrow(WindSel))
DatSel$windCal <- 0
for(b in 1:(length(splits) - 1)){
    DatSel$windCal[DatSel$DT > (WindSel$DT[splits[b] + 1] - 30) & DatSel$DT < (WindSel$DT[splits[b + 1] -1] + 30)] = 1
}
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$forageLong[b:nrow(DatSel)] != 1)){
        DatSel$tFromFor[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromFor[b] <- difftime(DatSel$DT[b + min(which(DatSel$forageLong[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# find points where there a wind calculations and <1 hour to foraging
DatSel$windNear <- DatSel$windCal == 1 & (DatSel$tFromFor < (60*60))

# # find turning points (0 xings of UTME and UTMN)
# pos <- which(diff(diff(DatSel$UTME) > 0) != 0) + 2
# DatSel$turn[pos] <- 1
# # timelag until turn
# b=1
# while(b<nrow(DatSel)){
#     if(all(DatSel$turn[b:nrow(DatSel)] != 1)){
#         DatSel$tFromTurn[b:nrow(DatSel)] <- NA
#         b <- nrow(DatSel)
#     } else {
#         DatSel$tFromTurn[b] <- difftime(DatSel$DT[b + min(which(DatSel$turn[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
#         b = b+1
#     }
# }
# where nearby wind calculations (5 mins) and nearby turns (<5 mins) overlap
# DatSel$turnNear <- DatSel$tFromTurn < (5*60) & DatSel$windNear == 1
# alternatively, find turning points of > 30 degrees
DatSel$dx = c(NA,diff(DatSel$UTME))
DatSel$dy = c(NA,diff(DatSel$UTMN))
DatSel$angle = atan2(DatSel$dy, DatSel$dx)
DatSel$lgTurn = DatSel$angle*(180/pi) > 30
# time to large turn
b=1
while(b<nrow(DatSel)){
    if(all(DatSel$lgTurn[b:nrow(DatSel)] != 1)){
        DatSel$tFromLgTurn[b:nrow(DatSel)] <- NA
        b <- nrow(DatSel)
    } else {
        DatSel$tFromLgTurn[b] <- difftime(DatSel$DT[b + min(which(DatSel$lgTurn[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
        b = b+1
    }
}
# where nearby wind calculations (5 mins) and nearby large turns (>30 deg, <5 mins) overlap
DatSel$lgTurnNear <- DatSel$tFromLgTurn < (5*60) & DatSel$windNear == 1

# create a mean vector for each selected period
st = which(diff(DatSel$turnNear) == 1) + 1
ed = which(diff(DatSel$turnNear) == -1) + 1
if(length(st) < length(ed)){
    st <- c(1, st)
}
if(length(ed) < length(st)){
    ed <- c(ed, nrow(DatSel))
}
for(b in 1:length(st)){
    DatSel$mnWindDx <- mean
}

######################################################################################################################
########################### TAKE NEAR FORAGING AND FIND AVE HEADING 30 MINS PRIOR ####################################
######################################################################################################################

for(tg in 1:length(Dat)){
    DatSel = Dat[[tg]]
    DatSel$forageLong[DatSel$tkb == 1 | DatSel$dv == 1] = 1
    DatSel$forageLong[is.na(DatSel$forageLong)] = 0
    WindSel = WindDat[WindDat$ID == unique(DatSel$tagID),]
    dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
    splits <- which(dt > 70)
    splits <- c(0, splits, nrow(WindSel))
    WindSel$nrFor = NA
    for(b in 1:nrow(WindSel)){
        WindSel$nrFor[b] = any(DatSel$forageLong[DatSel$DT >= WindSel$DT[b] & DatSel$DT < (WindSel$DT[b] + 3600)] == 1)
    }
    WindSel$BSpd = NA
    for(b in 1:(length(splits) - 1)){
        cord.dec <- SpatialPoints(cbind(WindSel$lon[(splits[b]+1):splits[b+1]], WindSel$lat[(splits[b]+1):splits[b+1]]),
            proj4string=CRS("+proj=longlat"))
        cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

        X <- cord.UTM$coords.x1
        Y <- cord.UTM$coords.x1
        vg_x_obs <- X[2:length(X)] - X[1:(length(X) - 1)]
        vg_y_obs <- Y[2:length(Y)] - Y[1:(length(Y) - 1)]
        times = WindSel$DT[(splits[b]+1):splits[b+1]]
        dt = as.numeric(difftime(times[2:length(times)],times[1:(length(times)-1)], units = 'secs'))
        g_speed <- sqrt(vg_y_obs^2 + vg_x_obs^2)/dt
        WindSel$BSpd[(splits[b]+1):splits[b+1]] = c(NA, g_speed)
    }
    if(tg == 1){
        WFor <- WindSel[WindSel$nrFor == 1,]
    } else {
        toAdd <- WindSel[WindSel$nrFor == 1,]
        WFor <- rbind(WFor, toAdd)
    }
}
WFor$dir <- NA
WFor$dir[(WFor$head - WFor$WHead) > (-pi/4) & (WFor$head - WFor$WHead) < (pi/4)] = "Same direction"
WFor$dir[(WFor$head - WFor$WHead) < (-pi/4) & (WFor$head - WFor$WHead) > (-3*pi/4)] = "L-Side"
WFor$dir[(WFor$head - WFor$WHead) < (-3*(pi/4)) | (WFor$head - WFor$WHead) > (3*(pi/4))] = "Toward"
WFor$dir[(WFor$head - WFor$WHead) < (3*pi/4) & (WFor$head - WFor$WHead) > (pi/4)] = "R-Side"
# plot the difference between the two
ggplot(WFor) +
    geom_point(aes(x = (head - WHead)*(180/pi), y = WSpeed)) +
    coord_polar() + scale_x_continuous(limits = c(0,180))

ggplot(WFor) +
    geom_point(aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed)) +
    coord_polar() + scale_x_continuous(limits = c(0,360))

ggplot() +
    geom_point(data = WFor[(WFor$head - WFor$WHead) < pi/2 & (WFor$head - WFor$WHead) >(-pi/2),],
        mapping=aes(x = (head - WHead), y = WSpeed), colour = "green") +
        geom_point(data = WFor[(WFor$head - WFor$WHead) < (-pi/2) | (WFor$head - WFor$WHead) >pi/2,],
        mapping=aes(x = (head - WHead), y = WSpeed), colour = "red") +
    coord_polar() + scale_x_continuous(limits = c(-pi,pi))



ggplot() +
    geom_point(data = WFor[((WFor$head+3) - (WFor$WHead + 3))*(180/pi) < 90,],
        mapping=aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed), colour = "green") +
        geom_point(data = WFor[((WFor$head+3) - (WFor$WHead + 3))*(180/pi) > 90,],
        mapping=aes(x = ((head+3) - (WHead+3))*(180/pi), y = WSpeed), colour = "red") +
    coord_polar() + scale_x_continuous(limits = c(0, 360))

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed, fill = dir), pch = 21) +
    coord_polar() + scale_x_continuous(limits = c(-pi, pi))

Wplot <- ggplot(WFor) +
    geom_point(aes(x = (WHead), y = WSpeed), pch = 2) +
    coord_polar() + scale_x_continuous(limits = c(-pi,pi))

Bplot <- ggplot(WFor) +
    geom_point(aes(x = head, y = BSpd), pch = 4)+
    coord_polar() + scale_x_continuous(limits = c(-3,3))

ggarrange(Wplot, Bplot)


dt <- difftime(WindDat$DT[2:nrow(WindDat)],WindDat$DT[1:(nrow(WindDat) - 1)], "units", "secs")
splits <- which(dt > 0 & dt > 70)
splits <- c(0, splits, nrow(WindDat))
WindDat$BSpd = NA
for(b in 1:(length(splits) - 1)){
    cord.dec <- SpatialPoints(cbind(WindDat$lon[(splits[b]+1):splits[b+1]], WindDat$lat[(splits[b]+1):splits[b+1]]),
        proj4string=CRS("+proj=longlat"))
    cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

    X <- cord.UTM$coords.x1
    Y <- cord.UTM$coords.x1
    vg_x_obs <- X[2:length(X)] - X[1:(length(X) - 1)]
    vg_y_obs <- Y[2:length(Y)] - Y[1:(length(Y) - 1)]
    times = WindDat$DT[(splits[b]+1):splits[b+1]]
    dt = as.numeric(difftime(times[2:length(times)],times[1:(length(times)-1)], units = 'secs'))
    g_speed <- sqrt(vg_y_obs^2 + vg_x_obs^2)/dt
    WindDat$BSpd[(splits[b]+1):splits[b+1]] = c(NA, g_speed)
}

# plot bird speed against wind heading
ggplot(WindDat, aes(x = WHead, y = BSpd)) +
    geom_point() + coord_polar()





# BWcomp <- data.frame(aWh=0,aWs=0,aBh=0,aBs=0,nrFor=0)
# for(b = 1:(length(splits)-1)){
#     BWcomp$aWh[b] = atan2(mean(WindSel$X[(splits[b] + 1):splits[b+1]]), mean(WindSel$Y[(splits[b] + 1):splits[b+1]]))
#     BWcomp$aWs[b] = mean(sqrt(WindSel$X[(splits[b] + 1):splits[b+1]]^2 + WindSel$Y[(splits[b] + 1):splits[b+1]]^2))
#     BWcomp$aBh[b] = 
# }



# dt <- difftime(WindSel$DT[2:nrow(WindSel)],WindSel$DT[1:(nrow(WindSel) - 1)], "units", "secs")
# splits <- which(dt > 70)
# splits <- c(0, splits, nrow(WindSel))
# DatSel$windCal <- 0
# for(b in 1:(length(splits) - 1)){
#     DatSel$windCal[DatSel$DT > (WindSel$DT[splits[b] + 1] - 30) & DatSel$DT < (WindSel$DT[splits[b + 1] -1] + 30)] = 1
# }
# b=1
# while(b<nrow(DatSel)){
#     if(all(DatSel$forageLong[b:nrow(DatSel)] != 1)){
#         DatSel$tFromFor[b:nrow(DatSel)] <- NA
#         b <- nrow(DatSel)
#     } else {
#         DatSel$tFromFor[b] <- difftime(DatSel$DT[b + min(which(DatSel$forageLong[b:nrow(DatSel)] == 1)) - 1], DatSel$DT[b], units = "secs")
#         b = b+1
#     }
# }
# windDiff <- as.numeric(difftime(tail(WindSel$DT, -1), head(WindSel$DT, -1), units = 'secs'))
# splits <- which(windDiff > (60*2))
# WindSel$aveX <- NA
# WindSel$aveY <- NA
# for(b in 1:length(splits)){
#     if(b == 1){
#         WindSel$aveX[1:splits[b]] <- sum((WindSel$X[1:splits[b]]), na.omit = T)
#         WindSel$aveY[1:splits[b]] <- sum((WindSel$Y[1:splits[b]]), na.omit = T)
#     } else if(b == length(splits)){
#         WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveX[splits[b]:nrow(WindSel)] <- sum((WindSel$X[splits[b]:nrow(WindSel)]), na.omit = T)
#         WindSel$aveY[splits[b]:nrow(WindSel)] <- sum((WindSel$Y[splits[b]:nrow(WindSel)]), na.omit = T)
#     } else {
#         WindSel$aveX[splits[b - 1]:splits[b]] <- sum((WindSel$X[splits[b - 1]:splits[b]]), na.omit = T)
#         WindSel$aveY[splits[b - 1]:splits[b]] <- sum((WindSel$Y[splits[b - 1]:splits[b]]), na.omit = T)
#     }
# }

# WindSel$headAve <- atan2(WindSel$aveY, WindSel$aveX)

# plot(WindSel$headAve)

# ggplot(data = DatSel, aes(x = Lon, y = Lat)) +
# geom_path() +
# geom_point(data = DatSel[DatSel$windCal == 1, ], aes(x = Lon, y = Lat), pch = 1) + 
#     geom_spoke(data = WindSel, arrow = arrow(length = unit(WindSel$FlSpeed/max(WindSel$FlSpeed)*0.15, 'inches')),
#         aes(x = lon, y = lat, angle = headAve, col = FlSpeed, radius = scales::rescale(FlSpeed, c(.1, .5)))) +
#     scale_colour_gradient("Wind speed", low = "yellow", high = "red")
#     #geom_point(data = DatSel[DatSel$windCal == T, ], aes(x = Lon, y = Lat), pch = 1)



# ggplot(data = DatSel) +
# geom_path(aes(x = Lon, y = Lat)) + 
# geom_point(data = DatSel[DatSel$turnNear == 1, ], aes(x = Lon, y = Lat))
# geom_point(aes(x = Lon, y = Lat, col = tFromFor)) +
# geom_point(data = DatSel[DatSel$Forage == 1,], aes(x = Lon, y = Lat), col = 'red')

sel <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/gribSelected.csv", sep = ",", header = T)
sel$EHead <- atan2(sel$X,sel$Y)
sel$WHead <- atan2(sel$U,sel$V)
ggplot() +
    geom_point(aes(y = EHead, x = WHead), data = sel, pch = 21, fill = "deepskyblue1") +
    geom_line(aes(x = (-pi:pi), y = (-pi:pi))) +
    theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated headings") +
    scale_x_continuous("JMA headings")
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/WindCorr.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()
library(circular)
res <- cor.circular(sel$WHead, sel$EHead, test = T)
res
ggplot() +
    geom_point(aes(x = WSpd, y = ESpd), data= sel, pch = 21, fill = "orangered2") +
    geom_line(aes(x = (0:12), y =(0:12))) +
    theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)") +
    scale_x_continuous("JMA wind speed (m/s)")
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/SpeedCorr.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed, fill = dir), pch = 21, size = 2) +
    coord_polar() + scale_x_continuous("Wind heading relative to bird (rad)",limits = c(-pi, pi)) +
    scale_fill_discrete("Wind heading") +
    theme_bw() + 
    theme(panel.border = element_blank(), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)")
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindCol.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

ggplot() +
    geom_point(data = WFor,
        mapping=aes(x = head - WHead, y = WSpeed), pch = 21, fill = "black", size = 1) +
    coord_polar() + scale_x_continuous("Wind heading relative to bird (rad)",limits = c(-pi, pi)) +
    theme_bw() + 
    theme(panel.border = element_blank(), text = element_text(size = 14,
        family = "Arial"), axis.text = element_text(size = 14, family = "Arial")) +
    scale_y_continuous("Estimated wind speed (m/s)")
ggsave("/Volumes/GoogleDrive/My Drive/PhD/Conferences/2021SeabirdSymposium/RelWindBW.png", , dpi = 300, height = 6,
    width = 6, units = "in", family = "Arial")
dev.off()

cor.test(sel$ESpd, sel$WSpd)


    scale_fill_viridis_d(option="magma")
    scale_fill_manual(values=c("#ffffbf","#ffffbf","#2c7bb6","#d7191c"))


###################################################################################################
####################################### HYPLIT ####################################################
###################################################################################################

library(devtools)
devtools::install_github("rich-iannone/splitr")
library(splitr)
install.packages("here",,"https://mac.R-project.org")
library(here)
library(dplyr)
library(ggplot2)
library(plyr)

setwd("~/Documents/SplitR_wd")

trajectory <- 
  hysplit_trajectory(
    lat = 42.83752,
    lon = -80.30364,
    height = 10,
    duration = 24,
    daily_hours = c(0, 6, 12, 18),
    direction = "forward",
    met_type = "gdas1",
    extended_met = TRUE)

if(Sys.info()['sysname'] == "Darwin"){
    exec_loc <- "/Users/aran/hysplit/"
} else {
    exec_loc <- "C:/hysplit/"
}

trajectory <- 
  hysplit_trajectory(
    lat = 39,
    lon = 143,
    height = 10,
    duration = 4,
    days = "2018-09-01",
    daily_hours = c(14),
    direction = "backward",
    met_type = "reanalysis",
    extended_met = FALSE,
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  )

trajectory_model <-
  create_trajectory_model() %>%
  add_trajectory_params(
    lat = 39,
    lon = 143,
    height = 10,
    duration = 6,
    days = "2018-09-01",
    daily_hours = c(14),
    direction = "backward",
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
  run_model()

trajectory_tbl <- trajectory_model %>% get_output_tbl()

trajectory_tbl
trajectory_tbl %>% trajectory_plot()
trajectory_model %>% trajectory_plot()


dispersion_model <-
  create_dispersion_model() %>%
  add_source(
    name = "particle",
    lat = 41.0, lon = 143.0, height = 10,
    rate = 5, pdiam = 15, density = 1.5, shape_factor = 0.8,
    release_start = lubridate::ymd_hm("2018-09-01 14:31"),
    release_end = lubridate::ymd_hm("2018-09-01 14:31") + lubridate::hours(2)
  ) %>%
  add_dispersion_params(
    start_time = lubridate::ymd_hm("2018-09-01 14:31"),
    end_time = lubridate::ymd_hm("2018-09-01 14:31") + lubridate::hours(6),
    direction = "backward", 
    met_type = "reanalysis",
    # met_dir = here::here("met"),
    exec_dir = paste(exec_loc, "exec", sep = "")
  ) %>%
    run_model()

dispersion_tbl <- dispersion_model %>% get_output_tbl()
dispersion_model %>% dispersion_plot()


find_hull <- function(df) df[chull(df$lon, df$lat), ]
hulls <- ddply(dispersion_tbl, "hour", find_hull)

plot <- ggplot(data = dispersion_tbl, aes(x = lon, y = lat, colour=as.factor(hour), fill = as.factor(hour), group = as.factor(hour))) +
# geom_point() + 
geom_polygon(data = hulls, alpha = 0.5) +
labs(x = "Lon", y = "Lat")
plot

#################################################################
#################### ANALYSE 2019 YONE DATA #####################
#################################################################

# bring in wind and foraging data
if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
} else {
    fileloc <- "E:/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
    forLoc <- "E:/My Drive/PhD/Data/2019Shearwater/TxtDat/AxyTrek/AlgorithmOutput/PredictedForage/"
}

files <- dir(fileloc)
forFiles <- dir(forLoc, pattern = "*ForageGPS.txt")

tags = unique(sub("_S.*", "",forFiles))

# read in the foraging data
# forD <- vector(mode="list", length=length(tags))
# for (t in 1:length(tags)) {
#     dayFiles = forFiles[grepl(paste("^",tags[t],"_S.*",sep=""),forFiles)]
#     ds <- data.frame(DT=character(),lat=numeric(),lon=numeric(),Forage=integer())
#     outpt <- ds
#     for (d in 1:length(dayFiles)){
#         outpt <- rbind(outpt,read.delim(paste(forLoc,dayFiles[d],sep=""),sep=",",header=T))
#     }
#     outpt$DT <- as.POSIXct(outpt$DT,"%d-%b-%Y %H:%M:%S",tz="")
#     outpt$ID <- tags[t]
#     w.dec <- SpatialPoints(cbind(outpt$Lon,outpt$Lat),proj4string = CRS("+proj=longlat"))
#     UTMdat <- spTransform(w.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
#     outpt$UTME <- UTMdat$coords.x1
#     outpt$UTMN <- UTMdat$coords.x2
#     forD[[t]] <- outpt
# }
# read in wind data
wind <- vector(mode="list",length=length(tags))
for(b in 1:length(files)){
    wind[[b]] <- read.delim(paste(fileloc,files[b],sep=""),sep=",",header=T)
    wind[[b]]$time <- as.POSIXct(wind[[b]]$time,"%Y/%m/%d,%H:%M:%S",tz="")
    w.dec <- SpatialPoints(cbind(wind[[b]]$lon,wind[[b]]$lat),proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(w.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    wind[[b]]$UTME <- UTMdat$coords.x1
    wind[[b]]$UTMN <- UTMdat$coords.x2
    wind[[b]] <- wind[[b]][!is.nan(wind[[b]]$Resnorm),]
    wind[[b]]$timeTo <- NA
    wind[[b]]$distTo <- NA
    for (g in 1:nrow(wind[[b]])){
        if (any(forD[[b]]$Forage[forD[[b]]$DT > wind[[b]]$time[g]] == 1)) {
            point <- which(abs(forD[[b]]$DT - wind[[b]]$time[g]) == min(abs(forD[[b]]$DT - wind[[b]]$time[g])))
            forInd <- min(which(forD[[b]]$Forage[point:nrow(forD[[b]])] == 1)) + point - 1 # select the minimum (i.e. the next one)
            wind[[b]]$timeTo[g] <- as.numeric(difftime(forD[[b]]$DT[forInd],wind[[b]]$time[g],units="secs"))
            wind[[b]]$distTo[g] <- sqrt( (forD[[b]]$UTMN[forInd] - forD[[b]]$UTMN[point])^2 + 
                (forD[[b]]$UTME[forInd] - forD[[b]]$UTME[point])^2 ) * 10^(-3)
        }
    }
    wind[[b]]$rwh <- wind[[b]]$aveDir - wind[[b]]$wDir
    wind[[b]]$ID <- tags[b]
}
dat$time <- as.POSIXct(dat$time, format = "%Y/%m/%d,%H:%M:%S", tz = "")
dat$rwh <- dat$aveDir - dat$wDir
dat$rwh[dat$aveDir == 0] <- NA
dat$U <- dat$wSp*cos(dat$wDir)
dat$V <- dat$wSp*sin(dat$wDir)

selectDat <- dat[!is.na(dat$rwh),]

colnames(dat)
unique()
ggplot() + 
    geom_point(data = dat[dat$ID == unique(dat$ID)[1],],
        aes(x=time,y=rwh))

ggplot(dat[dat$distTo < 10 & dat$forage != 1,]) + geom_point(aes(y = rwh, x = distTo)) + scale_x_continuous(limits=c(-180,180))
hist(dat$distTo)

ggplot(windAll) + 
    geom_point(aes(x = rwh*(pi/180), y = distTo)) + coord_polar() + scale_x_continuous(limits = c(-pi,pi))

FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
FkOshi.dec <- SpatialPoints(cbind(FkOshi$Long,FkOshi$Lat),proj4string=CRS('+proj=longlat'))
FkOshi.UTM <- spTransform(FkOshi.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
FkUTME <- FkOshi.UTM$coords.x1
FkUTMN <- FkOshi.UTM$coords.x2
for (b in 1:length(wind)){
    wind[[b]]$distFromFk <- sqrt((wind[[b]]$UTMN - FkUTMN)^2 + (wind[[b]]$UTME - FkUTME)^2) * 10^(-3)
    wind[[b]]$tSinceStart <- difftime(wind[[b]]$time, wind[[b]]$time[1], units="secs")
}

# go through each individuals tracks
ggplot(wind[[6]]) +
    geom_point(aes(x=rwh,y=timeTo/60)) + coord_polar() #+
    #geom_line(aes(x=rwh,y=distFromFk))
    # + facet_grid(~ ID) 
colnames(windAll)

ggplot(windAll[windAll$distTo < 30,]) +
    geom_point(aes(x = rwh*(pi/180), y = distTo, fill = ID),pch=21, alpha=.3) + coord_polar(start=pi)


dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
if(any(gribFls == "Z__C_RJTD_20190824150000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"))

###################################################################
#################### 2019 YONE WIND FORAGING ######################
###################################################################


if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "/Volumes/GoogleDrive/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("E:/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("E:/My Drive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "E:/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
}
D19 <- bind_rows(Dat19)
allD <- data.frame(DT=D19$DT,
    lat = D19$Lat,
    lon = D19$Lon,
    tagID = D19$tagID,
    Day = D19$Day,
    Sex = D19$Sex,
    distTrav = D19$recalDist,
    spTrav = D19$spTrav,
    recalSp = D19$recalSp,
    distFk = D19$distFromFk,
    tripN = D19$tripN,
    tripL = D19$tripL,
    tkb = D19$tkb,
    dv = D19$dv,
    UTME = D19$UTME,
    UTMN =  D19$UTMN)
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
# allD$forBeh <- NA
# allD$forBeh[allD$dv == 1] <- "Dive"
# allD$forBeh[allD$tkb == 1] <- "Surf"
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")

# CALCULATE RELATIVE WIND CONDITIONS FROM ESTIMATES WITH TIME/DISTANCE TO FORAGING
if(Sys.info()['sysname'] == "Darwin"){
    windLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
} else {
    windLoc = "E:/UTokyoDrive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
}
windFiles <- dir(windLoc)
for(b in 1:length(windFiles)){
  if(b == 1){
    WindDat <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
    colnames(WindDat) <- c("DT","Lat","Lon","Head","wDir","wSp","Resnorm")
    WindDat$ID <- sub("*WindYone.txt", "", windFiles[b])
    Wind.dec <- SpatialPoints(cbind(WindDat$Lon,WindDat$Lat), proj4string = CRS("+proj=longlat"))   
    UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    WindDat$UTME <- UTMdat$coords.x1
    WindDat$UTMN <- UTMdat$coords.x2
  } else {
    toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
    colnames(toAdd) <- c("DT","Lat","Lon","Head","wDir","wSp","Resnorm")
    toAdd$ID <- sub("*WindYone.txt", "", windFiles[b])
    Add.dec <- SpatialPoints(cbind(toAdd$Lon,toAdd$Lat), proj4string = CRS("+proj=longlat"))
    UTMdat <- spTransform(Add.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
    toAdd$UTME <- UTMdat$coords.x1
    toAdd$UTMN <- UTMdat$coords.x2 
    WindDat <- rbind(WindDat, toAdd)
  }
}
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y/%m/%d,%H:%M:%OS",tz="UTC")
# remove missing wind values
WindDat <- WindDat[!is.nan(WindDat$Resnorm),]
# remove points within 5km of FkOshi
FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
FkOshi.dec <- SpatialPoints(cbind(FkOshi$Long,FkOshi$Lat),proj4string=CRS('+proj=longlat'))
FkOshi.UTM <- spTransform(FkOshi.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
FkUTME <- FkOshi.UTM$coords.x1
FkUTMN <- FkOshi.UTM$coords.x2

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

#####################################################################################################
############################### BAYESIAN CIRCULAR MIXED EFFECTS MODEL ###############################
#####################################################################################################

fit.Motor <- bpnr(pred.I = Phaserad ~ 1 + Cond, data = Motor, its = 10000, burn = 100, n.lag = 3, seed = 101)


traceplot(fit.Motor,parameter="beta1")



#####################################################################################################
########################################### EXTRA FIGURES ###########################################
#####################################################################################################

ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    geom_path(data=ListD[[6]][ListD[[6]]$DT > as.POSIXct("2018/09/08 10:20:00") & ListD[[6]]$DT > as.POSIXct("2018/09/08 10:30:00"),], aes(x = Lon, y = Lat)) + geom_point(data=allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[6] & allD$forage == 1 & 
      allD$DT > as.POSIXct("2018/09/08 10:20:00") & allD$DT > as.POSIXct("2018/09/08 10:30:00"),],
    aes(x = lon, y = lat), pch = 21, fill = "deepskyblue") +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[6],], aes(x = Lon, y = Lat, colour = WSpd, angle = WHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = .5*(WindDat$WSpd[WindDat$yrID == indivWinds[6]]/max(WindDat$WSpd[WindDat$yrID == indivWinds[6]]))) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  geom_segment(aes(x=sbst$lon[seq(1,nrow(sbst)-1,200)],xend=sbst$lon[seq(2,nrow(sbst),200)],
    y=sbst$lat[seq(1,nrow(sbst)-1,200)],yend=sbst$lat[seq(2,nrow(sbst),200)]),arrow = arrow(length = unit(0.1,"inches"))) +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 8,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
    annotation_scale(location = 'br') +
    scale_y_continuous(breaks = c(39,40,41,42), labels = c("39","40","41","42"), name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(labels = c("140", "141", "142", "143", "144"), name = paste("Longitude (","\u00b0E",")", sep = ""))


# adding + 
#   geom_spoke(data=WindDat[WindDat$yrID == "2018_9",], aes(x = lon, y = lat, angle = WHead), colour = "black", arrow = arrow(length = unit(0.02,"inches")),
#   radius = .3*(WindDat$WSpeed[WindDat$yrID == "2018_9"]/max(WindDat$WSpeed[WindDat$yrID == "2018_9"])),size=1.3) +
#   geom_spoke(data=WindDat[WindDat$yrID == "2018_9",], aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.02,"inches")),
#   radius = .3*(WindDat$WSpeed[WindDat$yrID == "2018_9"]/max(WindDat$WSpeed[WindDat$yrID == "2018_9"]))) +
#   scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "GnBu") + 
#   theme_bw() + theme(panel.grid = element_blank()) +
#     theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8)) + 
#     annotation_scale(location = 'br') +
#     scale_y_continuous(breaks = c(39:44), labels = c("39","40","41","42","43","44"), name = paste("Latitude (","\u00b0N",")", sep = "")) +
#     scale_x_continuous(labels = c("140", "141", "142", "143", "144","145","146"), name = paste("Longitude (","\u00b0E",")", sep = ""))

# ggsave("/Volumes/GoogleDrive/My Drive/PhD/Admin/AORIPresentation/Animation/Finished.png", device = "png", dpi = 300, height = 5,
#     width = 5, units = "in")
# dev.off()
# colnames(WindDat)

#####################################################################################################
####################################### STRAIGHTNESS OF TRACK #######################################
#####################################################################################################

library(geosphere)
# create distances between positions using haversine method
HdistTrav <- lapply(unique(allD$yrID), function(x) c(NA,distHaversine(cbind(allD$lon[allD$yrID==x][1:(length(allD$lon[allD$yrID==x])-1)],allD$lat[allD$yrID==x][1:(length(allD$lat[allD$yrID==x])-1)]),
    cbind(allD$lon[allD$yrID==x][2:(length(allD$lon[allD$yrID==x]))],allD$lat[allD$yrID==x][2:(length(allD$lat[allD$yrID==x]))]))))
# finding minimum time (with try catch exception)
minT <- function(dt, x, window, units="mins") {
    if(!any(dt <= x - as.difftime(window,units=units))){
        out = 1
    } else {
        out = 
    }
}

ggplot(WindDat) +
    geom_point(aes(x=distTo,y=straightness))



colnames(WindDat)
unique(WindDat$bin10)[3]
x <- na.omit(WindDat$rwh[WindDat$bin10 == unique(WindDat$bin10)[3]])
y <- na.omit(WindDat$rwh[WindDat$bin10 == unique(WindDat$bin10)[1]])

cbind(sort(x %% (2*pi)),rep(1,length(x)))
2.4 %% (2*pi)

###############################################################################################
#################### CHECK THE GENERAL WIND DIRECTION FOR THE WHOLE PERIOD ####################
###############################################################################################

plot(WindDat$DT[year(WindDat$DT) == 2018],WindDat$WHead[year(WindDat$DT) == 2018])

ggplot(WindDat) + 
    stat_density(aes(x = WHead,fill=as.factor(year(DT))),alpha=.4)

ggplot(WindDat) +  
    geom_point(aes(x = lat, y = WHead))

###############################################################################################
################ COLLECT DATA FOR EACH FORAGING NUMBER, FIND AVERAGE WITH RBAR ################
###############################################################################################

breaks<-seq(from=0,to=round_any(max(WindDat$distTo,na.rm=T),5,f=ceiling),by=5)
WindDat$bin5 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)

uniqTr <- WindDat %>% dplyr::group_by(forNo,yrID,bin5,tripL > 2) %>% 
    dplyr::summarise(rbar = r.test(RelHead)$r.bar, pval = r.test(RelHead)$p.value,
    mnHead = circ.mean(RelHead))
uniqTr$dist <- as.numeric(sub(",.*","",sub("[][]","",as.character(uniqTr$bin5)))) + 5
uniqTr <- as.data.frame(uniqTr)
colnames(uniqTr) <- c("forNo","yrID","bin5","tripL","rbar","pval","mnHead","dist")
uniqTr$aligned <- uniqTr$mnHead + pi
uniqTr$aligned[uniqTr$aligned > pi] <- uniqTr$aligned[uniqTr$aligned > pi] - 2*pi
long <- ggplot(uniqTr[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == T,]) +
    stat_density(aes(x = aligned,colour = bin5),position = "identity", fill = NA, size = 1.1) +
    scale_colour_viridis(discrete = T)


ggplot(uniqTr[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == T,]) +
    geom_point(aes(x = aligned,colour = bin5),position = "identity", fill = NA, size = 1.1) +
    scale_colour_viridis(discrete = T)

short <- ggplot(uniqTr[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == F,]) +
    stat_density(aes(x = aligned,colour = bin5),position = "identity", fill = NA, size = 1.1) +
    scale_colour_viridis(discrete = T)

ggplot() +
    stat_density(aes(x = uniqTr$mnHead))

long <- ggplot(WindDat[WindDat$distTo < 50 & WindDat$tripL > 2,]) + stat_density(aes(x = aligned, colour = bin5),
    position = "identity", fill = NA, size = 1.1) +
    scale_colour_viridis(discrete=T)

short <- ggplot(WindDat[WindDat$distTo < 50 & WindDat$tripL <= 2,]) + stat_density(aes(x = aligned, colour = bin5),
    position = "identity", fill = NA, size = 1.1) +
    scale_colour_viridis(discrete=T)

ggarrange(long,short,nrow=1,common.legend = T)

###############################################################################################
############################ ATTEMPT AT CIRCULAR LINEAR REGRESSION ############################
###############################################################################################

# convert relative wind heading into cross, tail, or headwind category
WindDat$WindCat <- NA
WindDat$WindCat[abs(WindDat$RelHead) < (50*pi/180)] <- "Tailwind"
WindDat$WindCat[abs(WindDat$RelHead) > (50*pi/180) & abs(WindDat$RelHead) < (130*pi/180)] <- "SideWind"
WindDat$WindCat[abs(WindDat$RelHead) > (130*pi/180)] <- "HeadWind"
# WindDat$WindCat[WindDat$RelHead <= pi/4 & WindDat$RelHead >= -pi/4] <- "Tailwind"
# WindDat$WindCat[(WindDat$RelHead < -pi/4 & WindDat$RelHead > -3*pi/4) | (WindDat$RelHead > pi/4 & WindDat$RelHead < 3*pi/4)] <- "Sidewind"
# WindDat$WindCat[WindDat$RelHead >= 3*pi/4 | WindDat$RelHead <= -3*pi/4] <- "Headwind"
WindDat$WindCat <- as.factor(WindDat$WindCat)
WindDat$yrID <- as.factor(WindDat$yrID)
WindDat$seq <- as.factor(WindDat$seq)
WindDat$tripL <- as.factor(WindDat$tripL > 2)
WindDat$timeTo <- WindDat$timeTo/60
ggplot(WindDat) + 
    geom_point(aes(x = distTo, y = WSpeed, colour = WindCat))

library(lme4)
tst <- lmer(timeTo ~ WindCat + WSpeed + spTrav + tripL + (1 | yrID) + (1 | seq), data = WindDat)

summary(tst)
plot(tst)
lm.circular()

library(circglmbayes)
# Simulate data
x3 <- cbind(x = c(2001:2014, 2001:2014),
            g = c(rep(0, 14), rep(1, 14)))

y3 <- (1 + 2*atan(c(scale(x3) %*% c(.4, .1))) +
         rvonmises(28, mu = circular(0), kappa = 50)) %% (2*pi)

# Add interaction
y3 <- (y3 + c(14:1 / 6, rep(0, 14))) %% (2*pi)

# Collect and plot data
df <- data.frame(cbind(x3, y = y3))
plot(y ~ x, col = ifelse(df$g, "red", "blue"), df)

ggplot(df) + geom_point(aes(x = x, y = y, colour = g)) + coord_polar()

# Run model without interaction
cgm_noint <- circGLM(y ~ x + g, data = df)

# Run model with interaction
cgm_int <- circGLM(y ~ x * g, data = df)

cgm_noint
cgm_int

WindDat$rhd <- WindDat$RelHead + pi

bpnTest <- data.frame(rhd = as.numeric(WindDat$rhd),
    distTo = as.numeric(WindDat$distTo),
    timeTo = as.numeric(WindDat$timeTo),
    WSpeed = as.numeric(WindDat$WSpeed),
    yrID = as.factor(WindDat$yrID),
    year = as.factor(WindDat$year))

fit.wind <- bpnme(rhd ~ distTo*timeTo + WSpeed + (1|yrID) + (1|year),data=na.omit(bpnTest))
traceplot(fit.wind)


cgm_s <- circGLM(rhd ~ distTo + WSpeed, data = WindDat[!is.na(WindDat$distTo),])
df$y
WindDat$RelHead[1:20]
# basis of the model 
# relWindHead ~ distTo + wSpd + 

# remove non significant data
uniqTrSig <- uniqTr[uniqTr$pval < 0.05,]
distGaps <- seq(0,90,10)
pvalsUniq<-vector(mode="list",length=length(distGaps))
for(b in 1:length(distGaps)){
    RaylTS <- r.test(uniqTrSig$mnHead[uniqTrSig$tripL == F & uniqTrSig$dist > distGaps[b] & uniqTrSig$dist < (distGaps[b] + 10)])
    # tstS<-HR_test(wShort$rwh[wShort$distTo >= distGaps[b] & wShort$distTo < distGapsL[b]])
    RaylTL <- r.test(uniqTrSig$mnHead[uniqTrSig$tripL == T & uniqTrSig$dist > distGaps[b] & uniqTrSig$dist < (distGaps[b] + 10)])
    # tstL<-HR_test(wLong$rwh[wLong$distTo >= distGaps[b] & wLong$distTo < distGapsL[b]])
    pvalsUniq[[b]] <- data.frame(Distance=paste0(as.character(distGaps[b]),"-",as.character(distGaps[b]+10)),SRlP = RaylTS$p.value,SRlR = RaylTS$r.bar,
        LRlP = RaylTL$p.value,LRlR = RaylTL$r.bar)
}
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive/My Drive/PhD/Data/pvalsUniq.RData")
# } else {
#     load('E:/My Drive/PhD/Data/pvalsUniq.RData')
# }

# bayesian circular regression model based on projected normal distribution

# using a mixed-effects model. 
# first convert required columns to numeric
WindDat$RelHead <- as.numeric(WindDat$RelHead)
WindDat$distTo <- as.numeric(WindDat$distTo)
WindDat$WSpeed <- as.numeric(WindDat$WSpeed)
WindDat$tripL <- as.numeric(WindDat$tripL)
WindDat$OutHm <- as.numeric(WindDat$OutHm)
bcr <- bpnme(pred.I = RelHead ~ distTo + WSpeed + tripL + OutHm + (1|yrID),
    data = na.omit(WindDat), its = 1000, burn = 100, n.lag = 3, seed = 101)

WindDat$OutHm <- as.factor(WindDat$OutHm)

lmer(RelHead ~ distTo*timeTo + WSpeed + tripL > 2 + OutHm + (1|yrID), data=WindDat, REML=T)

library(mgcv)
library(lme4)
library(car)
library(circular)

WindDat$offset <- abs(WindDat$RelHead) # estimate a non-linear relative wind direction variable
WindDat$tripSL <- NA
WindDat$tripSL[WindDat$tripL > 2] = "L" # add a short/long trip factor
WindDat$tripSL[WindDat$tripL <= 2] = "S" # add a short/long trip factor
WindDat$yrID <- as.factor(WindDat$yrID)
WindDat$timeTo <- WindDat$timeTo/60 # convert timeTo to minutes
WindDat$year <- year(WindDat$DT)

ggplot(WindDat[WindDat$distTo < 100,], aes(x = distTo, y = RelHead)) + 
    geom_point()

# fitting a mixed model to the "offset" version of RelHead
distOff <- lmer(offset ~ distTo*timeTo + WSpeed + tripSL + (1 | yrID) + (1 | year), data = WindDat)
# residuals and qqplot show a lot of broken assumptions
plot(distOff)
summary(distOff)
qqPlot(residuals(distOff))
scatter.smooth(residuals(distOff) ~ fitted(distOff)) # residual plot

# add a random effect for yrID and foraging number
WindDat$yrIDForNo <- as.factor(paste(WindDat$yrID,as.character(WindDat$forNo),sep="_"))
# switch instead to using GAM with circular spline. In this case, RelHead becomes explanatory variable
distGam <- gam(distTo ~ s(RelHead, bs = "cc") + s(yrID, bs = "re") + s(year, bs = "re") + tripSL + s(WSpeed),
    data = WindDat, method = "REML")
gam.check(distGam)
summary(distGam)
par(mfrow=c(3,2))
plot(distGam,shade=T)





# attempt with lm.circular
comboPreds <- c()


ggplot(WindDat) + geom_point(aes(x = RelHead, y = distTo)) +
    geom_smooth(aes(x = RelHead, y = distTo))

plot(WindDat$RelHead, WindDat$spTrav)

wSpTrav <- lm.circular(y = circular(WindDat$RelHead+pi,units="radians"),
    x = )

lines(seq(min(residuals(distOff)),max(residuals(distOff)),length.out=nrow(WindDat)))



# max distances for each foraging point
# first, add foraging numbers to allD
allD$forNo <- NA
allD$forage[is.na(allD$forage)] <- 0
for(id in unique(allD$yrID)){
    forChg <- diff(allD$forage[allD$yrID == id])
    forSt <- which(forChg == 1) + 1
    forEd <- which(forChg == -1)
    if(forSt[1] > forEd[1]){
        forSt <- c(1,forSt)
    }
    if(forSt[length(forSt)] > forEd[length(forEd)]){
        forEd <- c(forEd,sum(allD$yrID == id))
    }
    for(ind in 1:length(forSt)){
        allD$forNo[allD$yrID == id][forSt[ind]:forEd[ind]] <- ind
    }
}
# now add foraging numbers for lead up as well 
allD$forNoLead <- NA
for(id in unique(allD$yrID)){
    forChg <- diff(allD$forage[allD$yrID == id])
    forSt <- which(forChg == 1) + 1
    forEd <- which(forChg == -1)
    if(forSt[1] > forEd[1]){
        forSt <- c(1,forSt)
    }
    if(forSt[length(forSt)] > forEd[length(forEd)]){
        forEd <- c(forEd,sum(allD$yrID == id))
    }
    for(ind in 1:length(forSt)){
        if(ind == 1){
            allD$forNoLead[allD$yrID == id][1:forEd[ind]] <- ind
        } else {
            allD$forNoLead[allD$yrID == id][(forEd[ind-1] + 1):forEd[ind]] <- ind
        }
    }
}
straightness <- function(dat, times, window){
    for(time in times){
        distHaversine(cbind(dat$lon[dat$DT == (time - as.difftime(window/2,units="mins"))],
                dat$lat[dat$DT == (time - as.difftime(window/2,units="mins"))]),
            cbind(dat$lon[dat$DT == (time + as.difftime(window/2,units="mins"))],
                dat$lat[dat$DT == (time + as.difftime(window/2,units="mins"))]))
        }
}
straightness(allD,allD$DT[which(allD$yrID == id & allD$forNoLead == b & allD$DT < 
            allD$DT[min(which(allD$yrID == id & allD$forNoLead == b & allD$forage == 1))-1])],
            5)

# add straightness index over a 5 minute period
# ratio of straight distance by distance travelled
allD$straightness <- NA
for(id in unique(allD$yrID)){
    for(b in unique(allD$forNoLead[allD$yrID == id])){
        # create 5 min sequence of time sequence for each value
        timeZ <- allD$DT[which(allD$yrID == id & allD$forNoLead == b & allD$DT < 
            allD$DT[min(which(allD$yrID == id & allD$forNoLead == b & allD$forage == 1))-1])]
        


        
        allD$forage[which(allD$yrID == id & allD$forNoLead == b)]
    }
}
allD[which(allD$yrID == id & allD$forNoLead == 2),]
# go through and calculate distances to foraging spot
for(id in unique(allD$yrID)){
    # go through each foraging number
    for(b in na.omit(unique(allD$forNo[allD$yrID == id]))){
        # find first foraging point
        nxtFor <- min(which(allD$forNo[allD$yrID == id] == b))
        if(nxtFor == 1){ # skip if first value is foraging
            next
        }
        if(b == 1){ 
            allD$distToNextFP[allD$yrID == id][1:(nxtFor-1)] <- distHaversine(cbind(allD$lon[allD$yrID == id][1:(nxtFor-1)],
                allD$lat[allD$yrID == id][1:(nxtFor-1)]),cbind(allD$lon[allD$yrID == id][nxtFor],allD$lat[allD$yrID == id][nxtFor]))/1000
        } else {
            allD$distToNextFP[allD$yrID == id][(max(which(allD$forNo[allD$yrID == id] == (b-1))) + 1):(nxtFor - 1)] <- distHaversine(cbind(allD$lon[allD$yrID == id][(max(which(allD$forNo[allD$yrID == id] == b-1))+1):(nxtFor-1)],
                allD$lat[allD$yrID == id][(max(which(allD$forNo[allD$yrID == id] == b-1))+1):(nxtFor-1)]),
                cbind(allD$lon[allD$yrID == id][nxtFor],allD$lat[allD$yrID == id][nxtFor]))
        }
    }
}

ggplot(allD,aes(x = timeD, y = distToNextFP, colour = yrID)) + geom_line()  


allD$distToNextFP <- allD$distToNextFP/1000

# add time from start (to normalise across individuals)
allD$timeD <- NA
for(id in unique(allD$yrID)){
    allD$timeD[allD$yrID == id] = difftime(allD$DT[allD$yrID == id],allD$DT[allD$yrID == id][1],units="mins")
}

spHd <- lm(spTrav ~ offset + WSpeed, data = WindDat)
summary(spHd)

spFun <- function(x) x - 
ggplot(WindDat, aes(x = offset, y = spTrav)) + geom_point() +

ggplot() + geom_line(colour = "red", data = spHdf, aes(x = predOut, y = spdPred))

ggplot() + 
    geom_point(data = WindDat, aes(x = offset, y= spTrav))
    
     +
     geom_line(data=spHdf, aes(x = predOut, y = spdPred), colour = "red")

WindDat$offset <- abs(WindDat$RelHead)

ggplot(allD, aes(x = timeD, y = distToNextFP, colour = yrID)) + geom_line()

################################################################################################
########################################## HMM FAILED ##########################################
################################################################################################

# HMM for movement modes in different winds
library(remotes)
install_github("bmcclintock/momentuHMM")
library(momentuHMM)
# separate into segments with consecutive data
sep <- which(abs(diff(WindDat$DT)) > 70)
WindDat$seq <- 0
for(b in 1:length(sep)){
    if(b == 1){
        WindDat$seq[1:sep[b]] <- b
    } else if(b == length(sep)){
        WindDat$seq[(sep[b-1]+1):sep[b]] <- b
        WindDat$seq[(sep[b]+1):nrow(WindDat)] <- b + 1
    } else {
        WindDat$seq[(sep[b-1]+1):sep[b]] <- b
    }
}
# format for momentuHMMData
colnames(WindDat)
# add distance travelled column
WindDat$distTrav <- NA
for(b in 1:nrow(WindDat)){
  inds <- which(allD$yrID == WindDat$yrID[b] & allD$DT > (WindDat$DT[b] - lubridate::seconds(5)) & allD$DT < (WindDat$DT[b] + lubridate::seconds(5)))
  WindDat$distTrav[b] <- mean(allD$distTrav[inds])
}
# add angle change
WindDat$angle <- NA
for(b in unique(WindDat$seq)){
    WindDat$angle[WindDat$seq == b] = c(atan2(diff(WindDat$UTME[WindDat$seq == b]),diff(WindDat$UTMN[WindDat$seq == b])),NA)
}
mformWD <- WindDat[c("lat","lon","RelHead","WSpeed","distTo","timeTo","OutHm","seq","tripL")]
colnames(mformWD) <- c("x","y","RelHead","WSpeed","distFP","timeFP","OutHm","ID","tripL")
# change trip length to category (1 = long -> 3+ days, 0 = short -> <=2 days)
mformWD$tripL <- as.factor(mformWD$tripL > 2)
# change OutHm to category (1 = outward, 0 = homeward)
mformWD$OutHm <- as.factor(mformWD$OutHm < 0)
# remove data where birds not approaching FPs
mformWD <- na.omit(mformWD)

# remove data with < 3 observations
mformWD <- as.data.frame(mformWD %>%
    group_by(ID) %>%
    filter(n() > 3)    )

lWD <- lapply(unique(mformWD$ID), function(x) prepData(mformWD[mformWD$ID == x,],
    type="LL", coordNames=c("y","x"), covNames=c("RelHead","WSpeed","distFP","tripL","OutHm")))

lWDSimple <- prepData(mformWD[mformWD$ID == unique(mformWD$ID)[15],],
    type="LL", coordNames=c("y","x"))

mformWD[mformWD$ID == 240,]

nbStates <- 3
stateNames <- c("transit","olfactory","visual")
dists <- list(step="gamma", angle = "wrpcauchy")
par0_m1 <- list()
plot(lWDSimple)

ggplot(na.omit(WindDat)) + 
    geom_point(aes(x = lon, y = lat, colour = yrID))

groupSum <- WindDat %>% group_by(seq) %>%
    summarise(n())
groupSum$seq[groupSum['n()'] == 107]



# create a list of distributions
nbStates <- 3 # number of behaviour states
stateNames <- c("transit","olfactory","visual")
dists = <- list(step = "gamma",
    angle = "wrpcauchy",
    RelHead = "vm", # relative wind heading
    WSpeed = "weibull", # wind speed
    distFP = "Poisson",
    OutHm = "Categorical")
Par0 = list(step = )

# create two models, one with ID (individual-level effects) and one without
MIfitHMM(lWD,nbStates=3,dist=list())
lWD[1]
prepData(mformWD[mformWD$ID == 1,],
    type="LL", coordNames=c("y","x"), covNames=c("RelHead","WSpeed","distFP","OutHm"))

install.packages("mitools")
library(mitools)
momentuHMM::MIfitHMM(mmntDat )

testdata <- as.data.frame(listWD[1],col.names=names(listWD[1]))


colnames(testdata) <- c("ID","x","y")
test <- prepData(listWD[1][c("ID","x","y")],type="LL")
# we have n series of consecutive data which can be used for momentuHMM data
momentuHMMData()

# Tutorial
### Load raw data
rawHaggis <- read.csv("C:/Users/arang/Downloads/rawHaggises.csv")
### Process data
processedHaggis<-prepData(data=rawHaggis,covNames=c("slope","temp"))
### Fit HMM
# initial step distribution natural scale parameters
stepPar0 <- c(1,5,0.5,3) # (mu_1,mu_2,sd_1,sd_2)
# initial angle distribution natural scale parameters
anglePar0 <- c(0,0,1,8) # (mean_1,mean_2,concentration_1,concentration_2)
fitHaggis <- fitHMM(data = processedHaggis, nbStates = 2,
dist = list(step = "gamma", angle = "vm"),
Par0 = list(step = stepPar0, angle = anglePar0),
formula = ~ slope + I(slope^2),
estAngleMean = list(angle=TRUE))
range(WindDat$RelHead)


processWind <- prepData(data=WindDat, covNames=c())


library(moveHMM)
set.seed(1122334455)

## Simulate covariate values 
# (slopes in degrees
# around 10 degrees - the latter won't affect the state switching in the model below).
# The model will be such that the haggises are most likely to be in the exploratory 
# state (state 2) when at slopes of around 20 degrees (the slope that perfectly matches 
# the difference in their leg lengths). For slopes close to 0 and slopes close to 40
# the animals become essentially immobile (due to the differences in their leg lengths) 
# and hence need to slowly crawl back (state 1) to slope levels better suited for them.
# We are using a quadratic predictor to achieve this setup.

slopes <- NULL
for (haggis in 1:15) {
  # for each of the 15 haggises
  arsim <- rep(NA,400)
  arsim[1] <- runif(1
  for (k in 2:400)
    arsim[k] <- 0.9*(arsim[k-1]-0.4)+0.4+rnorm(1
  
  slope <- 40*plogis(arsim)
  slopes <- c(slopes
}

slopes <- NULL		
for (haggis in 1:15) {		
  # for each of the 15 haggises	 simulate 400 slope values	
  arsim <- rep(NA,400)	
  arsim[1] <- runif(1,0.5,1)
  for (k in 2:400)		
    arsim[k] <- 0.9*(arsim[k-1]-0.4)+0.4+rnorm(1,0,0.7)
  		
  slope <- 40*plogis(arsim)		
  slopes <- c(slopes	slope)	
}		
		
temps <- NULL		


temps <- NULL
for (haggis in 1:15) {
  # for each of the 15 haggises
  arsim2 <- rep(NA
  arsim2[1] <- rnorm(1
  for (k in 2:400)
    arsim2[k] <- 0.9*(arsim2[k-1]-10)+10+rnorm(1
  
  temps <- c(temps
}

# data frame of covariate values
covs <- data.frame(slope=slopes

# specify parameters of the step length and turning angle distributions
stepPar <- c(1
anglePar <- c(pi
# (angle mean of -0.3 in state 2
# clockwise around hills when active)

# specify regression coefficients for the transition probabilities
beta <- matrix(c(-3.5
0.35
-0.01
0
               nrow=4

# simulate wild haggis movement data (15 haggises)
dataraw <- simData(nbAnimals=15
                   stepPar=stepPar
                   covs=covs

# only keep relevant columns (animals' ID
rawHaggis <- dataraw[
