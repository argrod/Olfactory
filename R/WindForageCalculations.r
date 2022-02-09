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

save(WindDat,file='E:/My Drive/PhD/Data/WindCalculations1819.RData')

WindDat$RelHead <- WindDat$head-WindDat$WHead
WindDat$RelHead[WindDat$RelHead < -pi] <- WindDat$RelHead[WindDat$RelHead < -pi] + 2*pi
WindDat$RelHead[WindDat$RelHead > pi] <- WindDat$RelHead[WindDat$RelHead > pi] - 2*pi
WindDat$aligned <- WindDat$RelHead + pi
WindDat$aligned[WindDat$aligned > pi] <- WindDat$aligned[WindDat$aligned > pi] - 2*pi


breaks<-seq(from=0,to=round_any(max(WindDat$distTo,na.rm=T),10,f=ceiling),by=10)
mnW <- ddply(WindDat, "bin10", summarise, grp.mean=mean(aligned))
WindDat$bin10 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- WindDat[WindDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))

ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(WindDat$distTo),], aes(x = aligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.2,show.legend=FALSE)+stat_density(aes(x=aligned, colour=bin10), size=1.5, geom="line",position="identity") +# coord_polar(start=pi) +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_colour_manual(name="Distance to next \nforaging spot (km)", values = rev(brewer.pal(9,"YlOrRd")),
    labels=paste(sort(unique(WindDat$bin10[WindDat$distTo < 90])),", n = ", as.character(unlist(bin10ns[,2])),sep="")) +
  scale_y_continuous(name="Proportion across all birds (%)", breaks=seq(0,0.3,.1),labels=seq(0,30,10))
ggsave(paste(figLoc,"DistRelDensity.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")