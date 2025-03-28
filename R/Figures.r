library(rerddap)
library(splitr)
library(plyr)
library(dplyr)
library(ggplot2)
library(devtools)
# devtools::install_github("ropenscilabs/rnaturalearth")
# devtools::install_github("ropenscilabs/rnaturalearthdata")
# install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")
# install.packages('sf')
library(sf)
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
library(circular)
library(wGribDat)
library(grid)
# devtools::install_github("keesmulder/circglmbayes")
library(circglmbayes)
library(lme4)
# library(ncdf4)
# install.packages("CircMLE")
library(CircStats)
library(Cairo)
library(extrafont)
library(CircMLE)
library(png)
library(ggspatial)
library(egg)
library(geosphere)
options(timeout = 800)

###################################################################################################################################
############################################################# READ IN #############################################################
###################################################################################################################################

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
# for(d in Dat){
#     d$distTrav <- c(NA,distHaversine(cbind(d$Lon[1:(nrow(d)-1)],d$Lat[1:(nrow(d)-1)]),
#         cbind(d$Lon[2:nrow(d)],d$Lat[2:nrow(d)])))
#     d$spTrav <- c(NA,d$distTrav[2:nrow(d)]/as.numeric(difftime(d$DT[2:nrow(d)],d$DT[1:(nrow(d)-1)],units="secs")))
# }
D18 <- bind_rows(Dat)
# for(d in Dat19){
#     d$distTrav <- c(NA,distHaversine(cbind(d$Lon[1:(nrow(d)-1)],d$Lat[1:(nrow(d)-1)]),
#         cbind(d$Lon[2:nrow(d)],d$Lat[2:nrow(d)])))
#     d$spTrav <- c(NA,d$distTrav[2:nrow(d)]/as.numeric(difftime(d$DT[2:nrow(d)],d$DT[1:(nrow(d)-1)],units="secs")))
# }
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
    distFk = c(D18$distFk, D19$distFromFk),
    tripN = c(D18$tripN, D19$tripNo),
    tripL = c(D18$tripL, D19$tripL))
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
allD$yrID <- paste(format(allD$DT,'%Y'),"_",sub("\\_S.*","",allD$tagID),sep="")
# allD$forBeh <- NA
# allD$forBeh[allD$dv == 1] <- "Dive"
# allD$forBeh[allD$tkb == 1] <- "Surf"
japan <- ne_countries(scale = "medium", country
 = "Japan", returnclass = "sf")
# distances <- allD %>% filter(forage==1) %>% group_by(tagID, forBeh, Year) %>% dplyr::summarise(mxdist=median(distFk))
# distances %>% group_by(forBeh,Year) %>% dplyr::summarise(mean(mxdist))
# FORAGING SPOT DISPERSALS
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/splitr/ForageDisps.RData")
# } else {
#     load("E:/My Drive/PhD/Data/splitr/ForageDisps.RData")
# }
# calculate foraging starts/ends
allD$forage[is.na(allD$forage)] <- 0
forSt <- which(diff(allD$forage) == 1) + 1
if(allD$forage[1] == 1){
  forSt <- c(1, forSt)
}
forEd <- which(diff(allD$forage) == -1)
if(allD$forage[nrow(allD)] == 1){
  forEd <- c(forEd, nrow(allD))
}

# # LOAD IN THE STEP LENGTHS TRAJECTORIES
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/splitr/StepsTrajTimeChgNew.RData")
# } else {
#     load("E:/My Drive/PhD/Data/splitr/StepsTrajTimeChgNew.RData")
# }
# # save(ListD, file="/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/splitr/ListD.RData")
# # LOAD IN THE LISTED DATA
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/splitr/ListD.RData")
# } else {
#     load("E:/My Drive/PhD/Data/splitr/ListD.RData")
# }
# LOAD WIND DATA
if(Sys.info()['sysname'] == "Darwin"){
  load("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalculations1819.RData")
} else {
  load("E:/My Drive/PhD/Data/WindCalculations1819.RData")
}
# allTraj <- bind_rows(outTraj)
# allTraj$relH <- allTraj$aveHd - allTraj$trjHd 

##################################################################################################################
################################################ OUTPUT LOCATIONS ################################################
##################################################################################################################

# figure locations
if(Sys.info()['sysname'] == "Darwin"){
  figLoc <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Figures/WindManuscript/"
  # figLoc <- "/Documents/GitHub/PhD/Olfactory/"
} else {
  figLoc <- "E:/My Drive/PhD/Figures/WindManuscript/"
  # figLoc <- "F:/Documents/GitHub/PhD/Olfactory/"
}
WindDat$WSpd <- sqrt(WindDat$X^2 + WindDat$Y^2)

# LESS THAN 10KM
ggplot(WindDat[WindDat$distTo < 10,]) + 
    geom_histogram(aes(x = RelHead), colour = "black", bins = 50) + 
    geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo < 10]), size = 1.3, linetype = 1) +
    coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("-180/180","-90","0","90"), limits = c(-pi,pi)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 12, family = "Arial"))
ggplot(WindDat) +
  geom_histogram(aes(x=RelHead,colour=yrID),position='dodge')

##################################################################################################################
################################### HEADINGS AND HR TEST FOR 10 KM BINS 100:10 ###################################
##################################################################################################################

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
pvals<-NA
for(b in 1:length(distGaps)){
    RaylT <- r.test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    tst<-HR_test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    # if(RaylT$p.value > 0.05){
    pvals[b] <- tst[2]
    if(tst[2] > 0.01){
      res <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) +
        geom_histogram(aes(x = RelHead), bins = round(nrow(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],])/2))
      ggplot_build(res)$data
        one2Ten[[b]] <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
          geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_sqrt(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p > 0.05", sep = ""))
    } else {
        roseplt <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 50,
              fill = "#d9d9d9") + scale_y_sqrt(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial"))
        avRelHd[b] <- circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
        one2Ten[[b]] <- roseplt + geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          xend = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          y = 0, yend = sqrt(max(ggplot_build(roseplt)$data[[1]]$count))*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = paste(as.character(distGapsL[b]),"  : ",as.character(distGaps[b]), "km, n = ", nrow(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]),
          ", ",length(unique(WindDat$yrID[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]))," individuals",sep="")) +
          annotate("text", x = pi/4,
            y = max(ggplot_build(roseplt)$data[[1]]$count),label = paste("p < ",as.character(signif.ceiling(RaylT$p.value,3)),sep=""))
    }
    # Cairo(width=4, height = 4, file = paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
    print(one2Ten[[b]])
    ggsave(paste(figLoc,"RelHeadRayleigh",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""), device="svg", dpi = 300, height = 3,
      width = 3, units = "in")
    dev.off()
    rm(tst)
}
summary(allD)
save(pvals,avRelHd,one2Ten, file="E:/My Drive/PhD/Data/WindCalc/headings.RData")

distGaps <- seq(0,90,10)
distGapsL <- distGaps+10
# avRelHd <- NA
outVals <- data.frame(dist=distGapsL,avRelHd=rep(NA,length(distGaps)),avDisp=rep(NA,length(distGaps)))
for(b in 1:length(distGaps)){
    RaylT <- CircStats::r.test(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]))
    # tst<-HR_test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    if(RaylT$p.value < 0.05){
    # if(tst[2] < .01){
      outVals$avRelHd[b] <- circ.mean(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]))
      outVals$avDisp[b] <- circ.disp((na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])))$var
    } else {
      outVals$avRelHd[b] <- NA
      outVals$avDisp[b] <- NA
    }
}

p1 <- ggplot(outVals) + geom_point(aes(x = avRelHd,y=dist)) + coord_polar(start=pi) + theme_bw() +
  scale_x_continuous(name="Average relative wind heading (rad)") +
  scale_y_continuous(name = "Distance to next foraging point (km)")+ 
  theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10),
    axis.text = element_text(size = 8))

p2 <- ggplot(outVals) + geom_point(aes(x = dist,y=avDisp)) + theme_bw() +
  scale_x_continuous(name = "Distance to next foraging point (km)") +
  scale_y_continuous(name = "Circular variance") + 
  theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10),
    axis.text = element_text(size = 8))

ggarrange(p1,p2,nrow=1,widths=c(1,1))
# apply(distGaps, function(x) watson.wheeler.test(WindDat$RelHead[WindDat$distTo >= x & WindDat$distTo < (x+1)]))

ggplot(WindDat, aes(x = spTrav, y = WSpeed)) + geom_point() +# coord_polar(start = pi) +
  # scale_colour_gradient(name = expression(paste("Wind speed (",ms^{-1},")")),low="blue",high="red") +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Ground speed (",ms^{-1},")",sep="")))

ggplot(WindDat[WindDat$rtChg < 0 & WindDat$distFk < 100,], aes(x = fminRelHd, y = WSpd)) + geom_point(pch =21, fill = "red") + coord_polar(start = pi) +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Estimated wind speeds (",ms^{-1},")",sep="")))
ggsave(paste(figLoc,"SpeedAngles.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")
unique(WindDat$yrID[WindDat$rtChg < 0 & WindDat$distFk < 100])
ggplot(WindDat[WindDat$distTo < 50,], aes(x = fminRelHd, y = distTo)) + geom_point(pch =21, fill = "red") + coord_polar(start = pi) +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Estimated wind speeds (",ms^{-1},")",sep="")))

ggplot(allListD[allListD$Forage == 1 & allListD$rtChg < 0 & allListD$tripL > 2,], aes(x = distFk)) + 
  geom_histogram(bins=50) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial")) +
  scale_x_continuous(name = "Distance from nest colony (km)") +
  scale_y_continuous("Count")

ggplot(WindDat[WindDat$distTo > 50,], aes(x = fminRelHd, y = WSpd)) + geom_point() + coord_polar(start=pi)
ggplot(WindDat[WindDat$DT < as.POSIXct("2019/01/01"),], aes(x = DT, y = WSpd)) + geom_point()


WindDat$fminRelHd <- WindDat$fminHd - WindDat$WHead
WindDat$fminRelHd[WindDat$fminRelHd < -pi] <- WindDat$fminRelHd[WindDat$fminRelHd < -pi] + 2*pi
WindDat$fminRelHd[WindDat$fminRelHd > pi] <- WindDat$fminRelHd[WindDat$fminRelHd > pi] - 2*pi

ggplot(data.frame(dist=distGaps[!is.na(avRelHd)],aveHd=avRelHd[!is.na(avRelHd)]), aes(x = aveHd, y = dist)) + geom_point(pch=21, colour="black",fill="deepskyblue") +
  coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name="Distance to next foraging (km)")
ggsave(paste(figLoc,"0-10Aves.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")
(watson.wheeler.test(WindDat$RelHead[WindDat$distTo < 10], group = WindDat$yrID[WindDat$distTo < 10]))
relNos <- WindDat[WindDat$distTo < 10 & WindDat$distTo > 2,] %>% group_by(yrID) %>% dplyr::summarise(n=n())
watson.wheeler.test(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID %in% relNos$yrID[relNos$n > 10]],
  group = WindDat$yrID[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID %in% relNos$yrID[relNos$n > 10]])

indivPlots = vector(mode="list",length=length(relNos$yrID[relNos$n > 10]))
for(b in 1:length(relNos$yrID[relNos$n > 10])){
  RaylT <- r.test(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b]])
    if(RaylT$p.value > 0.05){
        indivPlots[[b]] <- ggplot(WindDat[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b],]) + 
          geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
          labs(title = as.character(relNos$yrID[relNos$n > 10][b]))
    } else {
        roseplt <- ggplot(WindDat[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b],]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
        indivPlots[[b]] <- roseplt + geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b]]),
          xend = circ.mean(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b]]),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = as.character(relNos$yrID[relNos$n > 10][b]))
  }
}
indivHeads <- vector(mode="list",length=length(relNos$yrID[relNos$n > 10]))
for(b in 1:length(indivHeads)){
  RaylT <- r.test(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b]])
    if(RaylT$p.value < 0.05){
      indivHeads[[b]] = cbind(circ.mean(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 2 & WindDat$yrID == relNos$yrID[relNos$n > 10][b]]),
        RaylT$r.bar, relNos$yrID[relNos$n > 10][b])
    }
}
indivHeads[sapply(indivHeads, is.null)] <- NULL
indivHeads <- data.frame(matrix(unlist(indivHeads), nrow=length(indivHeads), byrow=TRUE))
colnames(indivHeads) <- c("AveHead","Rbar","yrID")
ggplot(indivHeads) +
  geom_segment(aes(x = as.numeric(AveHead), xend=as.numeric(AveHead),y=0,yend=as.numeric(Rbar),colour=yrID), lineend="round",
          arrow = arrow(length = unit(.25, "cm")), size=1.5) + coord_polar(start = pi)+ 
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + scale_y_continuous(name=expression(bar(r))) +
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA),
    text = element_text(size = 12,
    family = "Arial"), axis.text = element_text(size = 12, family = "Arial")) + scale_colour_discrete(name="Tag ID")
ggsave(paste(figLoc,"RelHeadIndivs.svg",sep=""), device="svg", dpi = 300, height = 10,
      width = 10, units = "in")


length(indivPlots)
ggarrange(indivPlots[[1]],indivPlots[[2]],indivPlots[[3]],indivPlots[[4]])
ggarrange(indivPlots[[5]],indivPlots[[6]],indivPlots[[7]],indivPlots[[8]])
ggarrange(indivPlots[[9]],indivPlots[[10]],indivPlots[[11]],indivPlots[[12]],indivPlots[[13]])

for(b in 1:length(distGaps)){
    RaylT <- r.test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    roseplt <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
        geom_point(aes(x = RelHead, y = WSpd)) + scale_y_continuous(name = "Count") +
        coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
        theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
            family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
    one2Ten[[b]] <- roseplt + geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
      xend = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
      y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
      arrow = arrow(length = unit(.25, "cm"))) +
      labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p < ",as.character(signif.ceiling(RaylT$p.value,3)),sep=""))
    print(one2Ten[[b]])
    ggsave(paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".eps",sep=""), device="eps", dpi = 300, height = 3.5,
      width = 3, units = "in")
}

ggplot(WindDat[WindDat$distTo < 10,]) + 
        geom_point(aes(x = RelHead, y = WSpd)) + scale_y_continuous(name = "Wind speed (m/s)") +
        coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
        theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
            family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
ggsave(paste(figLoc,"Less10Spd.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")

ggplot(WindDat[WindDat$rtChg < 0,]) + 
  geom_histogram(aes(x = RelHead), colour = "black", bins = 100, fill = "#d9d9d9") + scale_y_continuous(name = "Count") + coord_polar() +
  geom_segment(x = circ.mean(WindDat$RelHead[WindDat$rtChg < 0]),
          xend = circ.mean(WindDat$RelHead[WindDat$rtChg < 0]), y = 0, yend = max(ggplot_build(ggplot(WindDat[WindDat$rtChg < 0,]) + 
  geom_histogram(aes(x = RelHead), colour = "black", bins = 100, fill = "#d9d9d9"))$data[[1]]$count), colour = "#6a51a3", lineend="round",
          arrow = arrow(length = unit(.25, "cm")))
r.test(WindDat$RelHead[WindDat$rtChg < 0])

# ggplot(WindDat[WindDat$rtChg > 0,]) + 
#   geom_histogram(aes(x = RelHead), colour = "black", bins = 100, fill = "#d9d9d9") + scale_y_continuous(name = "Count") + coord_polar() +
#   geom_segment(x = circ.mean(WindDat$RelHead[WindDat$rtChg > 0]),
#           xend = circ.mean(WindDat$RelHead[WindDat$rtChg > 0]), y = 0, yend = max(ggplot_build(ggplot(WindDat[WindDat$rtChg > 0,]) + 
#   geom_histogram(aes(x = RelHead), colour = "black", bins = 100, fill = "#d9d9d9"))$data[[1]]$count), colour = "#6a51a3", lineend="round",
#           arrow = arrow(length = unit(.25, "cm")))
# r.test(WindDat$RelHead[WindDat$rtChg > 0])

geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          xend = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count), colour = "#6a51a3", lineend="round",
          arrow = arrow(length = unit(.25, "cm")))


ggarrange(one2Ten[[1]],one2Ten[[2]],one2Ten[[3]],one2Ten[[4]],one2Ten[[5]], nrow=2, ncol=3, labels=c("0-2km","2-4km","4-6km","6-8km","8-10km"))
circ.mean(WindDat$RelHead[WindDat$distTo < 10])
rose.diag(WindDat$RelHead[WindDat$distTo < 10],50,prop=3)
r.test(WindDat$RelHead[WindDat$distTo < 10])
r.test(WindDat$RelHead[WindDat$distTo < 1])

circ.mean(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 20])
r.test(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 20])
circ.mean(WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30])
circ.mean(WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40])

circ.mean(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))
circ.mean(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))
r.test(na.omit(allTraj$relH[allTraj$distTo < 10*10^3]))

allTraj$relH[allTraj$relH < -pi] <- allTraj$relH[allTraj$relH < -pi] + 2*pi
allTraj$relH[allTraj$relH > pi] <- allTraj$relH[allTraj$relH > pi] - 2*pi
ggplot(allTraj[allTraj$distTo < 10*10^3,], aes(x = relH, y = distTo)) + 
 geom_point() + coord_polar()

distGaps <- seq(0,40,10)
distGapsL <- distGaps+10
for(b in 1:(length(distGaps)-1)){
  print(watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]], WindDat$RelHead[WindDat$distTo >= distGaps[b+1] & WindDat$distTo < distGapsL[b+1]]))[[4]])
  # watson.two(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]], WindDat$RelHead[WindDat$distTo >= distGaps[b+1] & WindDat$distTo < distGapsL[b+1]], alpha = .05, plot = T)
}
install.packages("DescTools")
library(DescTools)

watson.two(WindDat$RelHead[WindDat$distTo >= 0 & WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40], alpha = .05, plot = T)
watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo >= 0 & WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40]))

watson.two(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2],plot=T,alpha=.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2], WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3],alpha=0.05,plot=T)
watson.two(WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3], WindDat$RelHead[WindDat$distTo >= 3 & WindDat$distTo < 4],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3], WindDat$RelHead[WindDat$distTo >= 4 & WindDat$distTo < 5],plot=T,alpha=0.05)

watson.two(WindDat$RelHead[WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 20],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo < 10 & WindDat$distTo > 20], WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40],plot=T,alpha=0.05)

ggplot(ListD[[6]]) + 
  geom_path(aes(x = Lon, y = Lat)) + geom_point(data=ListD[[6]][ListD[[6]]$Forage == 1,],aes(x=Lon,y=Lat), pch = 21, fill = "red") +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(140, 144), ylim = c(39, 42.5))

lengths<-seq(from=0,to=max(WindDat$distTo,na.rm=T)-1, by = 1)
lengthsL <- lengths+1
WindDat$trip <- NA
WindDat$trip[WindDat$tripL == 1] <- "Short"
WindDat$trip[WindDat$tripL > 1] <- "Long"
shrt <- WindDat[WindDat$trip == "Short",]
lng <- WindDat[WindDat$trip == "Long",]
shrtbinDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(shrt$RelHead[shrt$distTo >= lengths[x] & shrt$distTo < lengthsL[x]]))),
  disp = unlist(lapply(1:length(lengths), function(x) circ.disp(shrt$RelHead[shrt$distTo >= lengths[x] & shrt$distTo < lengthsL[x]])$var)))
shrtbinDat <- shrtbinDat[!is.na(shrtbinDat$disp),]
lngbinDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(lng$RelHead[lng$distTo >= lengths[x] & lng$distTo < lengthsL[x]]))),
  disp = unlist(lapply(1:length(lengths), function(x) circ.disp(lng$RelHead[lng$distTo >= lengths[x] & lng$distTo < lengthsL[x]])$var)))
lngbinDat <- lngbinDat[!is.na(lngbinDat$disp),]

# PLOT AVE HEADINGS IN 10KM BLOCKS AS THEY LEAVE FK ISLAND
leavelengths<-seq(from=0,to=max(WindDat$distFk)-1, by = 5)
leavelengthsL <- leavelengths+5
leavingshrt <- WindDat[((WindDat$OutHm < 0) & (WindDat$tripL == FALSE) & (WindDat$forNo == 1)),]
leaveShrtDat <- data.frame(dist = rep(NA,length(leavelengths)), aveHd = rep(NA,length(leavelengths)), disp = rep(NA,length(leavelengths)), hrP = rep(NA,length(leavelengths)))
for(b in 1:length(leavelengths)){
  if(sum(leavingshrt$distFk >= leavelengths[b] & leavingshrt$distFk < leavelengthsL[b]) > 20){
    leaveShrtDat$dist[b] = leavelengths[b]
    leaveShrtDat$aveHd[b] = circ.mean(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[b] & leavingshrt$distFk < leavelengthsL[b]])
    leaveShrtDat$disp[b] = circ.disp(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[b] & leavingshrt$distFk < leavelengthsL[b]])$rbar
    leaveShrtDat$hrP[b] = as.numeric(HR_test(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[b] & leavingshrt$distFk < leavelengthsL[b]])[2])
  }
}
save(leaveShrtDat,file="/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalc/shrtLeaveDat.RData")
leaveShrtDat <- leaveShrtDat[!is.na(leaveShrtDat$dist),]
# leaveShrtDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[x] & leavingshrt$distFk < leavelengthsL[x]]))),
  # disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[x] & leavingshrt$distFk < leavelengthsL[x]])$var)),
  # uniP = unlist(lapply(1:length(leavelengths), function(x) HR_test(leavingshrt$RelHead[leavingshrt$distFk >= leavelengths[x] & leavingshrt$distFk < leavelengthsL[x]])[2])))
leavinglng <- WindDat[WindDat$rtChg < 0 & WindDat$tripL > 2,]
leaveLngDat <- data.frame(dist = rep(NA,length(leavelengths)), aveHd = rep(NA,length(leavelengths)), disp = rep(NA,length(leavelengths)), hrP = rep(NA,length(leavelengths)))
for(b in 1:length(leavelengths)){
  if(sum(leavinglng$distFk >= leavelengths[b] & leavinglng$distFk < leavelengthsL[b]) > 20){
    leaveLngDat$dist[b] = leavelengths[b]
    leaveLngDat$aveHd[b] = circ.mean(leavinglng$RelHead[leavinglng$distFk >= leavelengths[b] & leavinglng$distFk < leavelengthsL[b]])
    leaveLngDat$disp[b] = circ.disp(leavinglng$RelHead[leavinglng$distFk >= leavelengths[b] & leavinglng$distFk < leavelengthsL[b]])$rbar
    leaveLngDat$hrP[b] = as.numeric(HR_test(leavinglng$RelHead[leavinglng$distFk >= leavelengths[b] & leavinglng$distFk < leavelengthsL[b]])[2])
  }
}
# leaveLngDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leavinglng$RelHead[leavinglng$distFk >= leavelengths[x] & leavinglng$distFk < leavelengthsL[x]]))),
#   disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leavinglng$RelHead[leavinglng$distFk >= leavelengths[x] & leavinglng$distFk < leavelengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(leavelengths), function(x) HR_test(leavinglng$RelHead[leavinglng$distFk >= leavelengths[x] & leavinglng$distFk < leavelengthsL[x]])[2])))
save(leaveLngDat,file="/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalc/lngLeaveDat.RData")
leaveLngDat <- leaveLngDat[!is.na(leaveLngDat$disp),]
# leaving <- WindDat[WindDat$rtChg < 0 ,]
# leaveDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leaving$RelHead[leaving$distFk >= leavelengths[x] & leaving$distFk < leavelengthsL[x]]))),
#   disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leaving$RelHead[leaving$distFk >= leavelengths[x] & leaving$distFk < leavelengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(leavelengths), function(x) r.test(leaving$RelHead[leaving$distFk >= leavelengths[x] & leaving$distFk < leavelengthsL[x]])$p.value)))
# leaveDat <- leaveDat[!is.na(leaveDat$disp),]

ggplot() + geom_point(data=leaveShrtDat[leaveShrtDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "red"), pch = 21) +
  geom_point(data=leaveLngDat[leaveLngDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "deepskyblue"), pch = 21) + 
  coord_polar(start=pi) + scale_x_continuous(name="Average heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
  scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  scale_fill_manual(name = "Trip length",values=c("red","deepskyblue"),labels=c("Long (2+ days)","Short (<2 day)"))
ggsave(paste(figLoc,"LeavingHeadings.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")

ggplot() + geom_point(data=leaveShrtDat[leaveShrtDat$hrP < 0.01 & leaveShrtDat$dist <= 200,], aes(x = aveHd, y = dist, fill = "red"), pch = 21) +
  geom_point(data=leaveLngDat[leaveLngDat$hrP < 0.01 & leaveLngDat$dist <= 200,], aes(x = aveHd, y = dist, fill = "deepskyblue"), pch = 21) + 
  coord_polar(start=pi) + scale_x_continuous(name="Average heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
  scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  scale_fill_manual(name = "Trip length",values=c("red","deepskyblue"),labels=c("Long (2+ days)","Short (1 day)"))
ggsave(paste(figLoc,"LeavingHeadings<200.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")

ggplot(WindDat[WindDat$rtChg < 0 & WindDat$distFk < 200 & WindDat$tripL > 2,], aes(x = RelHead)) + geom_density(alpha=.3) + 
  coord_polar(start = pi)

lm.circular(y=as.circular(WindDat$RelHead), x=WindDat$distFk, type="c-l",verbose=T, init=c(1))
x<-cbind(rnorm(10),rep(1,10))
y<-circular(2*atan(c(x%*%c(5,1))))+rvonmises(10, mu=circular(0), kappa=100)
lm.circular(y=y, x=x, init=c(5,1), type='c-l', verbose=TRUE)

speedFit <- lm.circular(y=as.circular(WindDat$RelHead),
  x = WindDat$spTrav, type = "c-l", verbose = T, init= c(1))


fit.Motor <- bpnr(pred.I = Phaserad ~ 1 + Cond, data = Motor, its = 10000, burn = 100, n.lag = 3, seed = 101)

fit.RelH <- bpnr(pred.I = (RelHead+pi) ~ 1 + distFk + (1|yrID) + WSpd, data = WindDat[WindDat$tripL > 2,], its = 10000, burn = 200, n.lag = 20, seed = 101)
traceplot(fit.RelH)

plot(WindDat$RelHead ~ WindDat$distFk)
WindDat <- WindDat[,-which(names(WindDat) %in% c("timeTo","forHd","spTrav","minHd"))]
WindDat$yrIDnm <- as.numeric(factor(WindDat$yrID))
testr <- bpnme(pred.I = RelHead ~ distFk + WSpd + (1|yrIDnm), data = WindDat[1:20,], its = 10000, burn = 1000, n.lag = 3)
WindDat$offset
tstr <- lm(offset ~ distTo + distFk + WSpd, data = WindDat)
summary(tstr)
plot(tstr)

fit.Maps <- bpnme(pred.I = Error.rad ~ Maze + Trial.type + L.c + (1|Subject), data = Maps, its = 10000, burn = 1000, n.lag = 3, seed = 101)

traceplot(testr)
summary(WindDat)

#####################################################################################
`################## DENSITY PLOTS FOR LONG AND SHORT ###`############################
#####################################################################################

geomplots <- vector(mode="list",length=10)
for(b in 1:length(distGaps)){
    geomplots[[b]] <- ggplot() + 
        stat_density(data = WindDat[WindDat$distTo > distGaps[b] & WindDat$distTo < (distGaps[b] + 10) & WindDat$tripL > 2,],
            aes(x = aligned, colour = "deepskyblue"),size=1.1, geom="line") +
        stat_density(data = WindDat[WindDat$distTo > distGaps[b] & WindDat$distTo < (distGaps[b] + 10) & WindDat$tripL <= 2,],
            aes(x = aligned, colour = "red"),size=1.1, geom="line") +
        scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) +
        ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
        scale_colour_manual(name = "Length", values = c("deepskyblue","red"),
          labels = c("Long","Short"))
}
geomplots[[3]]
ggarrange(geomplots[[1]],geomplots[[2]],geomplots[[3]],geomplots[[4]],
  geomplots[[5]],geomplots[[6]],geomplots[[7]],geomplots[[8]],geomplots[[9]],
  geomplots[[10]], nrow = 5, ncol = 2, labels = as.character(1:10))
####################################################################################
###################### ATTEMPT AT BUILDING A GAM ###################################
####################################################################################

# try building GAM
library(tidyverse)
set.seed(123)
trSize <- floor(nrow(WindDat) * .8)
trset <- sample(seq_len(nrow(WindDat)), size = trSize)
train.set <- WindDat[trset,]
test.set <- WindDat[-trset,]
library(mgcv)
gamWind <- gam(offset ~ distTo + WSpd + s(yrIDnm, bs = 're'), data = train.set, method = "REML")
summary(gamWind)
plot(gamWind)
summary(WindDat)

allLeavePlt <- ggplot() + geom_point(data=leaveShrtDat[leaveShrtDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "red"), pch = 21) +
  geom_point(data=leaveLngDat[leaveLngDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "deepskyblue"), pch = 21) + 
  coord_polar(start=pi) + scale_x_continuous(name="Average relative wind heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
  scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  scale_fill_manual(name = "Trip length",values=c("red","deepskyblue"),labels=c("Long (2+ days)","Short (<2 day)"))
library(grid)
vp <- viewport(width = 0.5,height=0.5,x=0.25,y=0.75)
smllPlt <- ggplot() + geom_point(data=leaveShrtDat[leaveShrtDat$hrP < 0.01 & leaveShrtDat$dist <= 100,], aes(x = aveHd, y = dist, fill = "red"), pch = 21) +
  geom_point(data=leaveLngDat[leaveLngDat$hrP < 0.01 & leaveLngDat$dist <= 100,], aes(x = aveHd, y = dist, fill = "deepskyblue"), pch = 21) + 
  coord_polar(start=pi) + scale_x_continuous(name="",breaks=c(-pi,-pi/2,0,pi/2), limits=c(-pi,pi)) +
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name="",breaks=c(50)) +
  theme(legend.position = "none", panel.border = element_blank(),axis.text.y = element_blank(),
         axis.ticks.y = element_blank(), axis.text.x = element_blank()) + annotate("text", x = pi/4, y = 50, label = "50", size = 3)

ggdraw() +
  draw_plot(allLeavePlt) +
  draw_plot(smllPlt, x = 0.067,y=.525,width=0.35,height=0.35)

# RETURNING
# PLOT AVE HEADINGS IN 10KM BLOCKS AS THEY ret FK ISLAND
retlengths<-seq(from=max(WindDat$distFk)-1, to = 5, by = -5)
retlengthsL <- retlengths+5
leavingshrt <- WindDat[WindDat$rtChg < 0 & WindDat$tripL <= 2,]
retShrtDat <- data.frame(dist = rep(NA,length(retlengths)), aveHd = rep(NA,length(retlengths)), disp = rep(NA,length(retlengths)), hrP = rep(NA,length(retlengths)))
for(b in 1:length(retlengths)){
  if(sum(leavingshrt$distFk >= retlengths[b] & leavingshrt$distFk < retlengthsL[b]) > 20){
    retShrtDat$dist[b] = retlengths[b]
    retShrtDat$aveHd[b] = circ.mean(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[b] & leavingshrt$distFk < retlengthsL[b]])
    retShrtDat$disp[b] = circ.disp(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[b] & leavingshrt$distFk < retlengthsL[b]])$rbar
    retShrtDat$hrP[b] = as.numeric(HR_test(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[b] & leavingshrt$distFk < retlengthsL[b]])[2])
  }
}
save(retShrtDat,file="/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalc/shrtretDat.RData")
retShrtDat <- retShrtDat[!is.na(retShrtDat$dist),]

# retShrtDat <- data.frame(dist = retlengths, aveHd = unlist(lapply(1:length(retlengths), function(x) circ.mean(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[x] & leavingshrt$distFk < retlengthsL[x]]))),
  # disp = unlist(lapply(1:length(retlengths), function(x) circ.disp(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[x] & leavingshrt$distFk < retlengthsL[x]])$var)),
  # uniP = unlist(lapply(1:length(retlengths), function(x) HR_test(leavingshrt$RelHead[leavingshrt$distFk >= retlengths[x] & leavingshrt$distFk < retlengthsL[x]])[2])))
leavinglng <- WindDat[WindDat$rtChg < 0 & WindDat$tripL > 2 & WindDat$distFk <= 200,]
retLngDat <- data.frame(dist = rep(NA,length(retlengths)), aveHd = rep(NA,length(retlengths)), disp = rep(NA,length(retlengths)), hrP = rep(NA,length(retlengths)))
for(b in 1:length(retlengths)){
  if(sum(leavinglng$distFk >= retlengths[b] & leavinglng$distFk < retlengthsL[b]) > 20){
    retLngDat$dist[b] = retlengths[b]
    retLngDat$aveHd[b] = circ.mean(leavinglng$RelHead[leavinglng$distFk >= retlengths[b] & leavinglng$distFk < retlengthsL[b]])
    retLngDat$disp[b] = circ.disp(leavinglng$RelHead[leavinglng$distFk >= retlengths[b] & leavinglng$distFk < retlengthsL[b]])$rbar
    retLngDat$hrP[b] = as.numeric(HR_test(leavinglng$RelHead[leavinglng$distFk >= retlengths[b] & leavinglng$distFk < retlengthsL[b]])[2])
  }
}
# retLngDat <- data.frame(dist = retlengths, aveHd = unlist(lapply(1:length(retlengths), function(x) circ.mean(leavinglng$RelHead[leavinglng$distFk >= retlengths[x] & leavinglng$distFk < retlengthsL[x]]))),
#   disp = unlist(lapply(1:length(retlengths), function(x) circ.disp(leavinglng$RelHead[leavinglng$distFk >= retlengths[x] & leavinglng$distFk < retlengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(retlengths), function(x) HR_test(leavinglng$RelHead[leavinglng$distFk >= retlengths[x] & leavinglng$distFk < retlengthsL[x]])[2])))
save(retLngDat,file="/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/WindCalc/lngretDat.RData")
retLngDat <- retLngDat[!is.na(retLngDat$disp),]

ggplot() + geom_point(data=retShrtDat[retShrtDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "red"), pch = 21) +
  geom_point(data=retLngDat[retLngDat$hrP < 0.01,], aes(x = aveHd, y = dist, fill = "deepskyblue"), pch = 21) + 
  coord_polar(start=pi) + scale_x_continuous(name="Average heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
  scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  scale_fill_manual(name = "Trip length",values=c("red","deepskyblue"),labels=c("Long (2+ days)","Short (<2 day)"))
ggsave(paste(figLoc,"ReturningHeadings.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")


# ggsave(paste(figLoc,"LeavingHeadingsInset.svg",sep=""), device="svg", dpi = 300, height = 5,
#       width = 5, units = "in")

# look at distribution of foraging distances from the nest colony

# ggplot(leaveDat[leaveDat$uniP < 0.05,], aes(x = aveHd, y = dist, fill = tripL > 1)) + geom_point(pch=21) + coord_polar(start=pi) + scale_x_continuous(name="Average heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
#   scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
#     family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
# colnames(WindDat)
# ggplot(WindDat[WindDat$rtChg < 0,], aes(x = RelHead, y=distFk,fill=trip)) + geom_point(pch=21,colour = 'black') + coord_polar(start=pi)
# ggplot(WindDat[WindDat$rtChg < 0,], aes(x = RelHead, fill = trip)) + geom_histogram(colour='black',bins=50,alpha=.6) + coord_polar(start=pi) +
#   scale_y_continuous(name="Count") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
#     family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
# ggsave(paste(figLoc,"LeavingHeadings.eps",sep=""), device="eps", dpi = 300, height = 4,
#       width = 4, units = "in")
ggplot(leaveDat[leaveDat$uniP < 0.05,], aes(y = disp, x = dist)) + geom_point() + scale_x_continuous(name="Distance from nest site (km)") +
  scale_y_continuous(name="Angular dispersal") + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
ggsave(paste(figLoc,"LeavingDispersals.eps",sep=""), device="eps", dpi = 300, height = 4,
      width = 4, units = "in")

lengths<-seq(from=0,to=max(WindDat$distTo)-1, by = 1)
lengthsL <- lengths+1 
binDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]]))),
  disp = unlist(lapply(1:length(lengths), function(x) circ.disp(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]])$var)))
# Cairo(width=8, height = 8, file = paste(figLoc,"DispersalOverDist.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
ggplot(binDat[binDat$dist < 150,], aes(x = dist, y = disp)) + geom_line() + geom_vline(xintercept=10, linetype="dotted") +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
            family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_x_continuous(name="Distance to next foraging (km)") +
    scale_y_continuous(name = "Angular dispersal")
ggsave(paste(figLoc,"DispersalOverDist.svg",sep=""), device="svg", dpi = 300, height = 4,
      width = 6, units = "in")
# GO THROUGH AND GET SAME PLOT FOR EACH INDIVIDUAL
dispDistPlots <- vector(mode="list",length=length(unique(WindDat$yrID)))
indivWinds <- unique(WindDat$yrID)
for(b in 1:length(dispDistPlots)){
  lengths<-seq(from=0,to=max(WindDat$distTo[WindDat$yrID == indivWinds[b]],na.rm=T)-1, by = 1)
  lengthsL <- lengths+1 
  binDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x] & WindDat$yrID == indivWinds[b]]))),
    disp = unlist(lapply(1:length(lengths), function(x) circ.disp(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x] & WindDat$yrID == indivWinds[b]])$var)))
  dispDistPlots[[b]] <- binDat[!is.na(binDat$disp),]
  dispDistPlots[[b]]$yrID <- rep(indivWinds[b], nrow(dispDistPlots[[b]]))
}
allDispDist <- bind_rows(dispDistPlots)
ggplot(allDispDist[allDispDist$dist < 10,], aes(x = dist, y = disp, colour = yrID)) + geom_line() +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_x_continuous(name="Distance to next foraging (km)") +
    scale_y_continuous(name = "Angular dispersal") + scale_colour_discrete(name = "Tag ID and year")
ggsave(paste(figLoc,"DispersalOverDistIndivs.svg",sep=""), device="svg", dpi = 300, height = 4,
      width = 6, units = "in")


breaks<-seq(from=0,to=round_any(max(WindDat$distTo,na.rm=T),10,f=ceiling),by=10)
WindDat$bin10 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(WindDat, "bin10", summarise, grp.mean=mean(rAligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- WindDat[WindDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
toFor <- ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90 & !is.na(WindDat$distTo),], aes(x = aligned, colour = bin10)) +#max(WindDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.5,show.legend=FALSE)+stat_density(aes(x=aligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(WindDat$bin10[WindDat$distTo < 90])))),", (", as.character(unlist(bin10ns[,2])),")",sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8)) + 
  scale_y_continuous(name="Proportion across all birds (%)", limits=c(0,.65),breaks=seq(0,0.6,.1),labels=seq(0,60,10))+
  annotate("text",x = -2.8, y = .6, label = "a)")
ggsave(paste(figLoc,"DistRelDensity.svg",sep=""), device="svg", dpi = 300, height = 8,
      width = 9, units = "cm")

ggplot(WindDat) + geom_point(aes(x = aligned, y = distFk)) + coord_polar()

ggplot(WindDat) + geom_point(aes(x = Fkbin10ns, y = aligned)) + 
    aggregate(aligned ~ as.factor(Fkbin10ns), WindDat, median)[,2]

ggplot(WindDat %>% group_by(Fkbin10) %>% summarise(n=n())) +
  geom_point(aes(x = Fkbin10,y=n))

plot(WindDat$distFk)

library(ggpubr)
breaks<-seq(from=0,to=200,by=50)
WindDat$Fkbin10 <- cut(WindDat$distFk, breaks = breaks, include.lowest=T,right=F)
WindDat$Fkbin10 <- sub("\\[","",as.character(WindDat$Fkbin10))
WindDat$Fkbin10 <- sub("\\]","",as.character(WindDat$Fkbin10))
WindDat$Fkbin10 <- sub(")","",as.character(WindDat$Fkbin10))
WindDat$Fkbin10 <- sub(",",":",as.character(WindDat$Fkbin10))
unique(WindDat$Fkbin10)
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
Fkbin10ns <- WindDat[WindDat$distFk < 200 & WindDat$distTo > 50 & WindDat$tripL,] %>% group_by(Fkbin10) %>% dplyr::summarise(length(unique(yrID)))
# FIX ORDER OF DISTANCES AND ONLY INCLUDE DATA PRIOR TO 50KM FROM FIRST VARIABLE
fromFk <- ggplot(WindDat[WindDat$distFk < 200 & WindDat$distTo > 50 & !is.na(WindDat$Fkbin10),]) + stat_density(aes(x = aligned, colour = Fkbin10),size=1.1,geom="line",position="identity") +  
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name = "Distance from \ncolony (km)", discrete=T,
  breaks=c("0:50","50:100","100:150","150:200"),
  labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(na.omit(WindDat$Fkbin10[WindDat$distFk < 500]))))),", (", as.character(unlist(na.omit(Fkbin10ns)[,2])),")",sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8)) + 
  scale_y_continuous(name="Proportion across all birds (%)", limits=c(0,.65), breaks=seq(0,0.6,.1),labels=seq(0,60,10)) +
  annotate("text",x = -2.8, y = .6, label = "b)")
ggsave(paste(figLoc,"FromFkDistRelDensity.svg",sep=""), device="svg", dpi = 300, height = 8,
      width = 9, units = "cm")

library(gridExtra)
gA <- ggplotGrob(toFor + rremove("ylab") + rremove("xlab"))
gB <- ggplotGrob(fromFk + rremove("ylab") + rremove("xlab"))
gA$widths <- gB$widths
combfig <- grid.arrange(gA,gB,nrow=2)

annotate_figure(combfig,left=textGrob("Proportion across all birds (%)", rot = 90,
  gp = gpar(fontsize = 10)),
  bottom = textGrob("Relative wind heading",hjust=.8,
    gp = gpar(fontsize = 10)),
  )
ggsave(paste(figLoc,"CombDistRelDensity.svg",sep=""), device="svg", dpi = 300, height = 10,
      width = 8.7, units = "cm")

# split into long and short foraging trips
LwDat <- WindDat[WindDat$tripL > 2,]
breaks<-seq(from=0,to=round_any(max(LwDat$distTo,na.rm=T),10,f=ceiling),by=10)
LwDat$bin10 <- cut(LwDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(LwDat, "bin10", summarise, grp.mean=mean(aligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- LwDat[LwDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
lDists <- ggplot(LwDat[LwDat$distTo > 0 &LwDat$distTo < 50,], aes(x = aligned, colour = bin5)) +#max(LwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin5), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(LwDat$bin5[LwDat$distTo < 50])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,0.9,.2),labels=seq(0,90,20),limits = c(0,.95)) +
  annotate("text",label="Long", x = -pi, y = .85, hjust = 0)

# split into long and short foraging trips
SwDat <- WindDat[WindDat$tripL <= 2,]
breaks<-seq(from=0,to=round_any(max(SwDat$distTo,na.rm=T),10,f=ceiling),by=10)
SwDat$bin10 <- cut(SwDat$distTo, breaks = breaks, include.lowest=T,right=F)
mnW <- ddply(SwDat, "bin10", summarise, grp.mean=mean(aligned))
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
bin10ns <- SwDat[SwDat$distTo < 90,] %>% group_by(bin10) %>% dplyr::summarise(length(unique(yrID)))
sDists <- ggplot(SwDat[SwDat$distTo > 0 &SwDat$distTo < 50,], aes(x = aligned, colour = bin5)) +#max(SwDat$distTo),], aes(x = aligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin5), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(SwDat$bin5[SwDat$distTo < 50])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,0.9,.2),labels=seq(0,90,20),limits=c(0,.95)) +
  annotate("text", label = "Short", x = -pi, y = .85, hjust = 0)

fig <- ggarrange(lDists,sDists,nrow=2,common.legend=T,legend = "right")
annotate_figure(fig,left=text_grob("Proportion of relative wind headings (%)",rot=90,size=10,vjust=1.5))
ggsave(paste(figLoc,"LSDistRelDensity.svg",sep=""), device="svg", dpi = 300, height = 10,
      width = 8.7, units = "cm")

ggplot(WindDat[WindDat$distTo < 50,]) + 
  geom_point(aes(x = RelHead, y = distTo,colour = tripL > 2)) + coord_polar()



###################################################################################
############################## WIND VS TRAVEL SPEED ###############################
###################################################################################

library(MASS)
ggplot(WindDat,aes(x = RelHead, y = spTrav)) + geom_point(size = .5) +
  scale_x_continuous(name = "Relative wind heading (rad)", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels = c(expression(pi),expression(-frac(pi,2)),0,expression(frac(pi,2)),expression(pi))) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(name = "Travel speed (kph)")
ggsave(paste0(figLoc,"SpeedRelHead.svg"), device = "svg", dpi = 300, height = 8,
  width = 8.7, units = "cm")

wDat <- na.omit(WindDat)
LoutwDat <- wDat[wDat$OutHm < 0 & wDat$tripL > 2,]
LhmwDat <- wDat[wDat$OutHm > 0 & wDat$tripL > 2,]
SoutwDat <- wDat[wDat$OutHm < 0 & wDat$tripL <= 2,]
ShmwDat <- wDat[wDat$OutHm > 0 & wDat$tripL <= 2,]
LoutDists <- ggplot(LoutwDat[LoutwDat$distTo > 0 &LoutwDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(LoutwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(LoutwDat$bin10[LoutwDat$distTo < 90])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,1,.1),labels=seq(0,100,10),limits=c(0,1.1))

SoutDists <- ggplot(SoutwDat[SoutwDat$distTo > 0 &SoutwDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(SoutwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(SoutwDat$bin10[SoutwDat$distTo < 90])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,1,.1),labels=seq(0,100,10),limits=c(0,1.1))

# split into long and short foraging trips
LhmDists <- ggplot(LhmwDat[LhmwDat$distTo > 0 &LhmwDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(LhmwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(LhmwDat$bin10[LhmwDat$distTo < 90])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,1,.1),labels=seq(0,100,10),limits=c(0,1.1))

ShmDists <- ggplot(ShmwDat[ShmwDat$distTo > 0 &ShmwDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(ShmwDat$distTo),], aes(x = rAligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")geom_density(alpha=.5,show.legend=FALSE)+
  stat_density(aes(x=aligned, colour=bin10), size=1.1, geom="line",position="identity") +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("") + theme_bw() + theme(panel.grid = element_blank()) +
  scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=paste(gsub(",",":",gsub('[[)]',"",sort(unique(ShmwDat$bin10[ShmwDat$distTo < 90])))),sep="")) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_y_continuous(name="", breaks=seq(0,1,.1),labels=seq(0,100,10),limits=c(0,1.1))

figAll <- ggpubr::ggarrange(LoutDists,SoutDists,LhmDists,ShmDists,nrow=2,ncol=2,labels=c("L Out","S Out","L Home","S Home"), common.legend = T, legend = "right")
annotate_figure(figAll,left=text_grob("Proportion of relative wind headings (%)",rot=90,size=10,vjust=1.5))

tst <- sapply(unique(allD$yrID), function(x) c(NA,distHaversine(cbind(allD$lon[which(allD$yrID==x)[1:(length(which(allD$yrID==x))-1)]],
  allD$lat[which(allD$yrID==x)[1:(length(which(allD$yrID==x))-1)]]),
  cbind(allD$lon[which(allD$yrID==x)[2:(length(which(allD$yrID==x)))]],
  allD$lat[which(allD$yrID==x)[2:(length(which(allD$yrID==x)))]]))))
tster <- NA
for(b in 1:length(tst)){
  if(b == 1){
    tster <- array(tst[[b]])
  } else {
    tster <- c(tster,array(tst[[b]]))
  }
}
allD$havDist <- tster

SLtripDists <- allD %>% dplyr::group_by(yrID,tripL>2,tripN) %>% dplyr::summarise(n = sum(havDist,na.rm=T))
colnames(SLtripDists) <- c("yrID","SL","tripN","n")
# mean short or long trip length
SLtripDists %>% dplyr::group_by(SL) %>% dplyr::summarise(mean(n))
# mean daily distance travelled
dailyDist <- allD %>% dplyr::group_by(yrID,Day) %>% dplyr::summarise(dist=sum(havDist,na.rm=T))
mean(dailyDist$dist/1000)

ggplot(allD[allD$tripL == 1,]) +
  geom_path(aes(x=lon,y=lat,colour=yrID)) +
  ggsn::scalebar(dist=100,model="WGS84",transform=T,dist_unit="km",x.min=min(allD$lon[allD$tripL == 1]),x.max=max(allD$lat[allD$tripL==1]),y.min=min(allD$lat[allD$tripL==1]),y.max=max(allD$lat[allD$tripL==1]))

sbst = allD[allD$yrID == indivWinds[6],]
WindDat$WSpd <- sqrt(WindDat$X^2 + WindDat$Y^2)
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    geom_path(data=ListD[[6]][ListD[[6]]$DT > as.POSIXct("2018/09/08 10:20:00") & ListD[[6]]$DT > as.POSIXct("2018/09/08 10:30:00"),], aes(x = Lon, y = Lat)) + geom_point(data=allD[allD$yrID == indivWinds[6] & allD$forage == 1 & 
      allD$DT > as.POSIXct("2018/09/08 10:20:00") & allD$DT > as.POSIXct("2018/09/08 10:30:00"),],
    aes(x = lon, y = lat), pch = 21, fill = "deepskyblue") +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[6],], aes(x = Lon, y = Lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.05,"inches")),
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
ggsave(paste(figLoc,"ExampleWind.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")

sbst = allD[allD$yrID == yrid[1] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),]
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(142.3, 142.7), ylim = c(39.5, 39.75)) +
    geom_path(data = allD[allD$yrID == yrid[1] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),],aes(x=lon,y=lat)) +
  geom_spoke(data = WindDat[WindDat$yrID == yrid[1] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00"),],
    aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpd[WindDat$yrID == yrid[1] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00")]/max(WindDat$WSpd[WindDat$yrID == yrid[1]]))) +  
  geom_point(data = allD[allD$yrID == yrid[1] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00") & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=2) +
  geom_segment(aes(x=sbst$lon[seq(1,nrow(sbst)-1,50)],xend=sbst$lon[seq(2,nrow(sbst),50)],
    y=sbst$lat[seq(1,nrow(sbst)-1,50)],yend=sbst$lat[seq(2,nrow(sbst),50)]),arrow = arrow(length = unit(0.08,"inches"))) +
  scale_x_continuous(name="Lon", breaks=seq(142.3,142.7,.2)) +
  scale_y_continuous(name="Lat", breaks=seq(39.5,39.7,.2)) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + 
    annotation_scale(location = 'br')
ggsave(paste(figLoc,"ExampleWindNearForage.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")

ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(140, 144), ylim = c(39, 42.5)) +
    geom_path(data = allD[allD$yrID == indivWinds[6],],aes(x=lon,y=lat)) +
  geom_point(data = allD[allD$yrID == indivWinds[6] & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=2) + 
  scale_x_continuous(name="Lon", breaks=seq(140,144,1)) +
  scale_y_continuous(name="Lat", breaks=seq(39,42.5,1)) +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + 
    annotation_scale(location = 'br')
ggsave(paste(figLoc,"ExampleForageAll.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")
    

ncCast <- st_cast(japan,"POLYGON")
# create plot for inset
nc2 = lapply(ncCast, function(x) as(ncCast,'Spatial'))
nc2@data$id = rownames(nc2@data)
nc2  = fortify(nc2) %>% inner_join(nc2@data)
#> Regions defined for each Polygons
#> Joining, by = "id"

# axis labels correctly moved to the right
ggplot() + geom_polygon(aes(x=long,y=lat),colour="grey",data = nc2[[2]]) + 
  scale_y_continuous(position='right') + coord_equal()
colnames(japan)

ggplot(allD[allD$yrID == yrid[1],],aes(x=lon,y=lat)) + geom_path() +
  geom_spoke(data = WindDat[WindDat$yrID == yrid[1] ,],
    aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.03,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpd[WindDat$yrID == yrid[1]]/max(WindDat$WSpd[WindDat$yrID == yrid[1]])))

yrid <- unique(WindDat$yrID)
sbst = allD[allD$yrID == yrid[1] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),]
inset <- ggplot() + geom_path(data = sbst,aes(x=lon,y=lat)) +
  # annotation_scale(location = 'br') +
  geom_spoke(data = WindDat[WindDat$yrID == yrid[1] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00"),],
    aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.03,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpd[WindDat$yrID == yrid[1] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00")]/max(WindDat$WSpd[WindDat$yrID == yrid[1]]))) +  
  geom_segment(aes(x=sbst$lon[seq(1,nrow(sbst)-1,50)],xend=sbst$lon[seq(2,nrow(sbst),50)],
    y=sbst$lat[seq(1,nrow(sbst)-1,50)],yend=sbst$lat[seq(2,nrow(sbst),50)]),arrow = arrow(length = unit(0.08,"inches"))) +
  geom_point(data = allD[allD$yrID == yrid[1] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00") & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=3) +
    scale_y_continuous(name="",breaks=c(39.5,39.6,39.7), position = "right",
      labels=as.character(c(39.5,39.6,39.7))) +
    scale_x_continuous(name="") +
  geom_text(data = data.frame(x=c(142.3,142.5,142.7),labs=as.character(c(142.3,142.5,142.7))),
            aes(x=x,label=labs, 
            y = 39.44), 
            hjust="middle", vjust="bottom",size=3.528*1.5) +
  geom_text(data = data.frame(y=c(39.5,39.6,39.7),labs=as.character(c(39.5,39.6,39.7))),
            aes(y=y,label=labs, 
            x = 142.695), 
            hjust="middle", vjust="middle",size=3.528*1.5) +
    scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
    theme(legend.position = "none",
      axis.ticks.length=unit(-0.25, "cm")) +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 10)) + 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title = element_blank(),
        axis.ticks.length=unit(-0.25, "cm"),
        axis.text = element_blank()) +
  theme(legend.position = "none") + 
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
  panel.background = element_rect(fill = "transparent",colour = NA)) +
  ggsn::scalebar(dist = 10, model = 'WGS84',transform=T,dist_unit="km", height = .05,
      st.dist = .1,x.min = min(sbst$lon), x.max = 142.65, y.min = min(sbst$lat), y.max = max(sbst$lat), location = 'topleft',box.fill=c("black","white"))
inset
# example wind with inset
full <- ggplot() + geom_path(data = allD[allD$yrID == yrid[1],], aes(x = lon, y = lat)) +
  geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  geom_point(data = allD[allD$yrID == yrid[1] & allD$forage == 1,], aes(x = lon, y = lat, fill = "deepskyblue"),pch = 21) + 
  geom_spoke(data = WindDat[WindDat$yrID == yrid[1],], aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.03,"inches")),
  radius = .3*(WindDat$WSpd[WindDat$yrID == yrid[1]]/max(WindDat$WSpd[WindDat$yrID == yrid[1]]))) + 
  coord_sf(xlim = c(139, 144), ylim = c(39, 42.5)) +
  scale_x_continuous(name="Longitude") +
  scale_y_continuous(name="Latitude") +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 10)) + 
    ggsn::scalebar(dist = 50, model = 'WGS84',transform=T,dist_unit="km", st.bottom = F,
      st.dist = .045,x.min = 139, x.max = 143.8, y.min = 39, y.max = 42.5, location = 'bottomright')
library(patchwork)
full + 
  inset_element(inset,0,0,0.5,.5,)
ggsave(paste(figLoc,"ExampleWindInset.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")    

a <- 5
WindDat[WindDat$yrID == indivWinds[a] & WindDat$distTo < 10 & (WindDat$RelHead < -2.6 | WindDat$RelHead > 2.6),c(1:3,26)]

a=8
sbst = allD[allD$yrID == indivWinds[a],]
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(141, 145), ylim = c(39, 43)) + geom_path(data=allD[allD$yrID == indivWinds[a],], aes(x = lon, y = lat)) + geom_point(data=allD[allD$yrID == indivWinds[a] & allD$forage == 1,], aes(x = lon, y = lat), pch = 21, fill = "deepskyblue") +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[a],], aes(x = Lon, y = Lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.05,"inches")),
  radius = .5*(WindDat$WSpd[WindDat$yrID == indivWinds[a]]/max(WindDat$WSpd[WindDat$yrID == indivWinds[a]]))) +
  geom_segment(aes(x=sbst$lon[seq(1,nrow(sbst)-1,300)],xend=sbst$lon[seq(2,nrow(sbst),300)],
    y=sbst$lat[seq(1,nrow(sbst)-1,300)],yend=sbst$lat[seq(2,nrow(sbst),300)]),arrow = arrow(length = unit(0.1,"inches"))) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 8,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
    annotation_scale(location = 'br') 
+
    scale_y_continuous(breaks = c(39,40,41,42), labels = c("39","40","41","42"), name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(labels = c("140", "141", "142", "143", "144"), name = paste("Longitude (","\u00b0E",")", sep = ""))

tmp50 <- vector(mode='list',length=20)
for(b in 1:length(yrid)){
  tmp50[[b]] <- sapply(unique(WindDat$forNo[WindDat$yrID == yrid[b]]), function(x)
    nrow(WindDat[WindDat$forNo == x & WindDat$yrID == yrid[b] & WindDat$distTo < 50,]))
}

## ANIMATION OF BIRD AND WIND TRAVELLING TOWARDS FORAGING POINT
# using tag 4
selWDat <- WindDat[WindDat$yrID == "2018_9" & WindDat$DT <= as.POSIXct('2018/09/04 5:55:09',format="%Y/%m/%d %H:%M:%S",tz="") & WindDat$DT > as.POSIXct('2018/09/04 4:0:00',tz=''),]
selD <- allD[allD$yrID == "2018_9" & allD$DT <= as.POSIXct('2018/09/04 5:55:09',format="%Y/%m/%d %H:%M:%S",tz="") & allD$DT > as.POSIXct('2018/09/04 4:0:00',format="%Y/%m/%d %H:%M:%S",tz=""),]
ggplot() + geom_point(aes(x = lon,y=lat),data=selD) + 
  geom_point(aes(x=lon,y=lat,fill=forage),pch=21,data=selD[selD$forage==1,]) +
  geom_spoke(data = selWDat, aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), radius = .03*selWDat$WSpeed/max(selWDat$WSpeed,na.omit=T))

library(ggsn)
library(png)
library(grid)
library(magick)

brd <- image_read("/Users/aran/Documents/GitHub/Olfactory/R/birdTop.png")
brd <- image_transparent(brd,'white')
raster <- as.raster(image)

p1 <- ggplot() +
    coord_sf(xlim = c(144.7, 145.01), ylim = c(42.7, 43)) + 
    geom_path(data=selD,aes(x = lon, y = lat), alpha = 0) +
    geom_spoke(data = selWDat, aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), radius = .03*selWDat$WSpeed/max(selWDat$WSpeed,na.omit=T), alpha = 0) +
    geom_point(aes(x=lon,y=lat,fill="#CC3300"),pch=21,data=selD[selD$forage==1,]) +
    scale_fill_manual(name = "Foraging", values="#CC3300",labels="", drop = F) +
    guides(fill = guide_legend(override.aes = list(fill="#CC3300",pch=21,alpha=1))) +
    ggsn::scalebar(dist=5,model="WGS84",transform=T,dist_unit="km",x.min=144.7,x.max=145.01,y.min=42.7,y.max=43,nudge_y=.025,location="bottomleft",border.size=.4) +
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10)) +
    scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
    scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p1
p1 + annotate(geom="raster",x=144.8,y=144.9,angle=14)
selD$birdAngle <- c(atan2(diff(selD$UTMN),diff(selD$UTME)),NA) / (pi/180)
# convert to 0 degrees as north coordinates
selD$birdAngle <- (-selD$birdAngle + 90) %% 360

# animation frames save location
if(Sys.info()['sysname'] == "Darwin"){
    dir_out <- "/Volumes/GoogleDrive-102199952889875375671/My Drive/PD/Figures/Visuals/WindToFP/"
} else {
    dir_out <- "I:/My Drive/PD/Figures/Visuals/WindToFP/"
}

for(b in selD$DT){
  adding <- p1 
  adding <- adding +
  geom_point(aes(x=lon,y=lat),fill="#CC3300",colour="#CC3300",pch=21,data=selD[selD$forage==1,],stroke=1) +
  annotation_custom(rasterGrob(magick::image_rotate(brd,selD$birdAngle[selD$DT == b]),height=.4,width=.4),
    xmin = selD$lon[selD$DT == b]-.05, ymin = selD$lat[selD$DT==b]-.05,
    xmax = selD$lon[selD$DT == b]+.05, ymax = selD$lat[selD$DT==b]+.05) + 
  geom_path(data=selD[selD$DT <= b,],aes(x = lon, y = lat)) +
  geom_spoke(data = selWDat[which(selWDat$DT <= b),], aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), radius = .03*selWDat[max(which(selWDat$DT <= b)),]$WSpeed/max(selWDat[selWDat$DT <= b,]$WSpeed,na.omit=T),lwd=1)
    
  fp <- file.path(dir_out, paste0(as.character(b),".png"))

  ggsave(adding, filename=fp, width = 5, height = 5, units = "in",
        device = "png")
}

p1 + annotation_custom(rasterGrob(magick::image_rotate(brd,90)))

imgs <- list.files(dir_out,full.names=T)
img_list <- lapply(imgs,image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 5)

## save to disk
image_write(image = img_animated,
            path = paste0(dir_out,"FPOlfactory.gif"))

WindDat[WindDat$yrID == yrid[1],],aes(x=lon,y=lat)
# RELATIVE HEADINGS AS BIRDS LEAVE COLONY
one2Ten <- vector(mode="list",length=10)
distGaps <- seq(0,9,1)
distGapsL <- distGaps+1
for(b in 1:length(distGaps)){
    RaylT <- r.test(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]))
    tst<-HR_test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    # if(RaylT$p.value > 0.05){
    if(tst[2] > 0.01){
        one2Ten[[b]] <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
          geom_histogram(aes(x = RelHead), colour = "black", bins = 30, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p > 0.05", sep = ""))
    } else {
        roseplt <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 30, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distFk >= distGaps[b] & WindDat$distFk < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
        one2Ten[[b]] <- roseplt + geom_segment(x = circ.mean(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])),
          xend = circ.mean(na.omit(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p < ",as.character(signif.ceiling(tst[2],3)),sep=""))
        # + geom_label(aes(x = pi/4, y = max(ggplot_build(roseplt)$data[[1]]$count)), label = paste("p value = ",as.character(signif(RaylT$p.value, 3)), sep = ""))
    }
    # Cairo(width=8, height = 8, file = paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
    print(one2Ten[[b]])
    # ggsave(paste(figLoc,"RelHeadLeavingHR",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""), device="svg", dpi = 300, height = 3.5,
      # width = 3, units = "in")
}
one2Ten[[3]]
# allList <- bind_rows(ListD)
# allList$levRet <- c(allList$levRet[!is.na(allList$levRet)],allList$levRet19[!is.na(allList$levRet19)])
# WindDat$levRet <- NA
# WindDat$tripL <- NA
# for(b in 1:nrow(WindDat)){
#   WindDat$levRet[b] <- allList$levRet[which(allList$DT > (WindDat$DT[b] - lubridate::seconds(30)) & allList$DT < (WindDat$DT[b] + lubridate::seconds(30)) & paste(allList$tagID,format(allList$DT,"%Y"),sep="") == WindDat$yrID[b])]
#   WindDat$tripL[b] <- allList$tripL[which(allList$DT > (WindDat$DT[b] - lubridate::seconds(30)) & allList$DT < (WindDat$DT[b] + lubridate::seconds(30)) & paste(allList$tagID,format(allList$DT,"%Y"),sep="") == WindDat$yrID[b])]
# }
# save(WindDat,file="E:/My Drive/PhD/Data/WindCalc/windDatAll.RData")

allTraj <- bind_rows(outTraj)
allTraj$distTo <- allTraj$distTo*10^-3
allTraj$relH <- allTraj$aveHd - allTraj$trjHd

allTraj$relH[allTraj$relH < pi] = allTraj$relH[allTraj$relH < pi] + 2*pi
allTraj$relH[allTraj$relH > pi] = allTraj$relH[allTraj$relH > pi] - 2*pi
allTraj <- na.omit(allTraj)
distGaps <- seq(0,9,1)
distGapsL <- distGaps+1
one2TenTr <- vector(mode="list",length=length(distGaps))
for(b in 1:length(distGaps)){
    RaylT <- r.test(allTraj$relH[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b]])
    # tst<-HR_test(allTraj$relH[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b]])
    if(RaylT$p.value > 0.05){
    # if(tst[2] > 0.05){
        one2TenTr[[b]] <- ggplot(allTraj[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b],]) + 
          geom_histogram(aes(x = relH), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p > 0.05", sep = ""))
    } else {
        roseplt <- ggplot(allTraj[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b],]) + 
            geom_histogram(aes(x = relH), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(allTraj$relH[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
        one2TenTr[[b]] <- roseplt + geom_segment(x = circ.mean(allTraj$relH[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b]]),
          xend = circ.mean(allTraj$relH[allTraj$distTo >= distGaps[b] & allTraj$distTo < distGapsL[b]]),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p < ",as.character(signif.ceiling(RaylT$p.value,3)),sep=""))
        # + geom_label(aes(x = pi/4, y = max(ggplot_build(roseplt)$data[[1]]$count)), label = paste("p value = ",as.character(signif(RaylT$p.value, 3)), sep = ""))
    }
    # Cairo(width=8, height = 8, file = paste(figLoc,"relH",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
    print(one2Ten[[b]])
    # ggsave(paste(figLoc,"relHHermansRasson",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""), device="svg", dpi = 300, height = 3.5,
    #   width = 3, units = "in")
}
one2TenTr[[4]]

ggplot(allTraj, aes(x = relH, y = distTo)) + geom_point() + coord_polar()

ggplot(WindDat[WindDat$distTo < 0,]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks =c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
ggsave(paste(figLoc,"EmptyWindPlot.svg",sep=""), device="svg", dpi = 300, height = 3.5,
      width = 3, units = "in")

# watson goodness of fit test for von mises dist. at 5km intervals
inters <- seq(52.5, 2.5, by=-1)
# create a circular class from RelHead
relHCirc <- circular(WindDat$RelHead, unit = "radians",rotation="counter")
aveHCirc <- circular(WindDat$head, units = "radians",rotation="counter")
lapply(inters, function(x) watson.test(relHCirc[WindDat$distTo > (x - 5) & WindDat$distTo <= x], alpha = 0.05, dist = "vonmises"))

# run for each individual
allOut = vector(mode="list",length=length(unique(WindDat$yrID)))
for(yrid in 1:length(unique(WindDat$yrID))){
  sel <- WindDat[WindDat$yrID == unique(WindDat$yrID)[yrid],]
  LcirDistEst <- bind_rows(lapply(inters, function(x) circ.disp(relHCirc[which(sel$distTo > (x - 2.5) & sel$distTo <= (x + 2.5) & sel$tripL > 2)])))
  # add the 90th and 10th %ile of est wind speed
  LcirDistEst[,c("ten","ninety")] <- bind_rows(lapply(inters, function(x) quantile(sel$WSpd[which(sel$distTo > (x - 2.5) & sel$distTo <= (x + 2.5) & sel$tripL > 2)], probs = c(.1,.9))))
  LcirDistEst$meanSpd <- unlist(lapply(inters,function(x) mean(sel$WSpd[which(sel$distTo > (x - 2.5) & sel$distTo <= (x + 2.5) & sel  $tripL > 2)])))
  LaveHDistEst <- bind_rows(lapply(inters, function(x) circ.disp(aveHCirc[which(sel$distTo > (x - 2.5) & sel$distTo <= (x + 2.5) & sel$tripL > 2)])))  
  allOut[[yrid]] <- c(data.frame("dis"=inters,"yrID"=unique(WindDat$yrID)[yrid]),LcirDistEst)
}
allRbar = bind_rows(allOut)
colnames(allRbar)
ggplot() + geom_line(data=allRbar,aes(x = dis, y = rbar,colour=yrID)) +
  # geom_line(data = data.frame(inters,aveHDistEst), aes(x = inters, y = rbar, linetype = "dotted")) +
  scale_x_reverse(name = "") + scale_y_continuous(name = expression(bar(italic(r))))

LcirDistEst <- bind_rows(lapply(inters, function(x) circ.disp(relHCirc[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL > 1)])))
# add the 90th and 10th %ile of est wind speed
LcirDistEst[,c("ten","ninety")] <- bind_rows(lapply(inters, function(x) quantile(WindDat$WSpd[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL > 1)], probs = c(.1,.9))))
LcirDistEst$meanSpd <- unlist(lapply(inters,function(x) mean(WindDat$WSpd[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL > 1)])))
LaveHDistEst <- bind_rows(lapply(inters, function(x) circ.disp(aveHCirc[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL > 1])))
Lnums <- sapply(inters, function(x) sum(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL > 1))

ScirDistEst <- bind_rows(lapply(inters, function(x) circ.disp(relHCirc[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 1)])))
# add the 90th and 10th %ile of est wind speed
ScirDistEst[,c("ten","ninety")] <- bind_rows(lapply(inters, function(x) quantile(WindDat$WSpd[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 1)], probs = c(.1,.9))))
ScirDistEst$meanSpd <- unlist(lapply(inters,function(x) mean(WindDat$WSpd[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 1)])))
SaveHDistEst <- bind_rows(lapply(inters, function(x) circ.disp(aveHCirc[which(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 1)])))
Snums <- sapply(inters, function(x) sum(WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 1))

sapply(inters, function(x) length(unique(WindDat$yrID[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5)])))

plot(rev(inters[cirDistEst$n>30]),rev(cirDistEst$rbar[cirDistEst$n>30]))

# dispersal of relHead vs dispersal of bird head
dispG <- ggplot() + geom_line(data=data.frame(inters, LcirDistEst),aes(x = inters, y = rbar,linetype="solid")) +
  geom_line(data=data.frame(inters, ScirDistEst),aes(x = inters, y = rbar,linetype="dotted")) +
  # geom_line(data = data.frame(inters,aveHDistEst), aes(x = inters, y = rbar, linetype = "dotted")) +
  scale_x_reverse(name = "Distance to next foraging point (km)") + scale_y_continuous(name = expression(bar(italic(r)))) +
  scale_linetype_manual(name = "", values=c("dotted","solid"), labels=c("Short trip\n(<3 days)", "Long trip\n(>2 days")) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) #+
  # annotate("text", x = 50, y = 0.8, label = "A)")

ggsave(paste(figLoc,"dispDist.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 8.7, units = "cm")



# dispG <- 
ggplot(data.frame(inters, cirDistEst), aes(x = inters, y = rbar)) + geom_line() +
  scale_x_reverse(name = "") + scale_y_continuous(name = expression(bar(italic(r)))) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  annotate("text", x = 50, y = 0.8, label = "A)")
ggsave(paste(figLoc,"dispDist.svg",sep=""), device="svg", dpi = 300, height = 3,
      width = 5, units = "in")
# spG <- 
spG <- ggplot(data.frame(inters, cirDistEst), aes(x = inters, y = rbar)) + 
  geom_ribbon(aes(ymin = ten, ymax = ninety), fill = 'grey70') +
  geom_line(aes(x = inters, y= meanSpd)) + scale_x_reverse(name = "Distance to next foraging spot (km)") + scale_y_continuous(name = expression(paste("Estimated wind speed (",ms^{-1},")")))+
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  annotate("text", x = 50, y = 8, label = "B)")
# ggarrange(dispG,spG, nrow=2,ncol=1, labels = c("a)","b)"),hjust = -1.25)
ggsave(paste(figLoc,"dispSpeed.svg",sep=""), device="svg", dpi = 300, height = 3,
      width = 5, units = "in")

# plot_grid(dispG,spG,ncol = 1, align="v",axis="lr")
# ggsave(paste(figLoc,"dispAll.svg",sep=""), device="svg", dpi = 300, height = 5,
      # width = 5, units = "in")

ggplot(data.frame(inters, cirDistEst), aes(x = inters, y = rbar)) + geom_line() + 
  geom_ribbon(aes(ymin = ten, ymax = ninety), fill = 'grey70') +
  scale_x_reverse(name = "Distance to next foraging spot (km)") + scale_y_continuous(name = expression(bar(italic(r))),
    sec.axis=sec_axis(~.*30,name="Wind speed")) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial"))
ggsave(paste(figLoc,"dispDist.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")

WindDat$fminRelH <- WindDat$fminHd-WindDat$WHead

plot(WindDat$fminRelH, WindDat$RelHead)

plot(WindDat$fminHd-WindDat$WHead,WindDat$Head-WindDat$WHead)


WindDat$fminRelH[WindDat$fminRelH < -pi] <- WindDat$fminRelH[WindDat$fminRelH < -pi] + 2*pi
WindDat$fminRelH[WindDat$fminRelH > pi] <- WindDat$fminRelH[WindDat$fminRelH > pi] - 2*pi
ggplot(WindDat[WindDat$distTo<50,],aes(y=distTo, x = fminRelH)) + geom_point() + coord_polar(start=pi)




plot(inters[aveHDistEst$n>30],aveHDistEst$rbar[aveHDistEst$n>30])

# model the change in dispersals
dispDat <- data.frame(disp = cirDistEst$rbar[cirDistEst$n>30], distTo = inters[cirDistEst$n>30])
dispMod <- lm(disp ~ poly(distTo,5,raw=T), data = dispDat)
plot(dispDat$distTo, dispDat$disp)
lines(seq(1,max(dispDat$distTo),1), predict(dispMod,data.frame(distTo=seq(1,max(dispDat$distTo),1))),col='red')

fit1<-fitdistr(dispDat$disp, 'exponential')
ks.test(dispDat$disp,fit1$estimate)

plot(aveHDistEst$rbar)

plot(RTestSpc)


plot(WindDat$Head[50:80],ylim=c(-pi,pi))
points(WindDat$forHd[50:80]+pi)

a <- 1
plot(rvonmises(100, vmEst[[a]]$mu, vmEst[[a]]$kappa))

allEst <- bind_rows(vmEst)
plot(allEst$kappa)
# search for change points in mean direction
mndrChg <- change.pt(WindDat$RelHead)

plotedf(WindDat$RelHead+pi)

circ.plot(WindDat$RelHead[WindDat$distTo < 5 & WindDat$distTo > 0])

ggplot(WindDat, aes(x = RelHead, y = WSpd, colour = distTo)) + geom_point() + coord_polar()

#calculate proportion of wind estimates vs gps durations
tagYrs <- unique(WindDat$yrID)
sumsts <- NA
windMins <- NA
totalDur <- NA
for(g in 1:length(tagYrs)){
  windMins[g] <- nrow(WindDat[WindDat$yrID == tagYrs[g],])
  AxDatFull <- allD[which(paste(allD$tagID,allD$Year,sep="") == tagYrs[g]),]
  timepoint <- AxDatFull$DT
	tsel <- seq(timepoint[1], timepoint[length(timepoint)], by = 60)
	dtFull <- as.numeric(diff(timepoint))
  # find data that line up to the new timepoints
  select <- NA
  for(b in 1:length(tsel)){
      choose <- which(timepoint >= (tsel[b] - (median(dtFull)/2)) & timepoint <= (tsel[b] + (median(dtFull)/2)))
      if(length(choose) != 0){
          if(length(choose) > 1){
              select[b] <- choose[which.min(abs(tsel[b] - timepoint[choose]))]
          } else {
              select[b] <- choose
          }
      }
  }
  AxDat <- AxDatFull$DT[na.omit(select)]
  totalDur[g] <- length(AxDat)
}
mean(windMins)
sd(windMins)
mean(totalDur)
sd(totalDur)

length(unique(WindDat))

mean(windMins[grepl("2018",tagYrs)]/totalDur[grepl("2018",tagYrs)])
sd(windMins[grepl("2018",tagYrs)]/totalDur[grepl("2018",tagYrs)])
mean(windMins[grepl("2019",tagYrs)]/totalDur[grepl("2019",tagYrs)])
sd(windMins[grepl("2019",tagYrs)]/totalDur[grepl("2019",tagYrs)])


allD %>% group_by(tagID,Sex,Year) %>% dplyr::summarise(max(distTrav))

foragePoints

(unique(WindDat$tripL[WindDat$distTo < 10]))
colnames(WindDat)
hist(WindDat$distTo)

# wind validation

outPt <- read.delim("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2019Shearwater/WindEst/Validation.csv", sep = ",", header=T)
outPt$Time <- as.POSIXct(outPt$Time,format="%Y-%m-%dT%H:%M:%S")
# cor.circular(outPt$est)

ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point(pch=21,fill="red") +
  stat_smooth(method = "lm", col = "blue")
}

res<-cor.circular(outPt$estHead,outPt$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
hdval <- ggplot(outPt, aes(x = gribHead, y = estHead)) +
    geom_point(pch=21,fill="deepskyblue") +
    geom_line(data=data.frame(x=-pi:pi,y=-pi:pi),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=-pi,y=2,label="corr = 0.4", hjust = 0) +
    annotate("text",x=-pi,y=1.5,label="p < 9 %*% 10^{-7}", parse = T, hjust = 0) +
    annotate("text",x=-pi,y=1,label="n == 152", parse = T, hjust = 0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

splm <- lm(estSpeed ~ gribSpeed, data = outPt)

spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(outPt$gribSpeed),y=0:max(outPt$gribSpeed)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=0,y=12.5,label="y = 0.52x + 0.65", hjust= 0) +
    annotate("text",x=0,y=11.5,label="p < 3 %*% 10^{-16}", parse = T, hjust= 0) +
    annotate("text",x=0,y=10.5,label=expression(paste(R^2," = 0.38")), hjust= 0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name="Estimated wind speed (m/s)") + scale_x_continuous(name=(("JMA wind speed (m/s)")))

ggarrange(hdval,spdval, ncol=1,nrow=2, labels=c("a)","b)"),hjust=-3,vjust=2)
ggsave(paste(figLoc,"windVal.svg",sep=""), device="svg", dpi = 300, height = 7,
      width = 3.5, units = "in")

# calculate the relative heading to next foraging spot
WindDat$forHd <- NA
for(b in 4008:nrow(WindDat)){
  nxtFor <- min(which(paste(allD$tagID,allD$Year,sep="") == WindDat$yrID[b] & allD$DT > WindDat$DT[b] & allD$forage == 1))
  colnames(WindDat)
  WindDat$forHd[b] <- atan2(allD$UTMN[nxtFor]-WindDat$UTMN[b],allD$UTME[nxtFor]-WindDat$UTME[b])
}

# WIND VALIDATION
valDat <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/gribs/gribSelectedProper.csv"

ggplot(WindDat[WindDat$distTo < 20,], aes(x = RelHead, y = spTrav, colour = WSpd)) +
  geom_point() + coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks =c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
      family = "Arial"), axis.text = element_text(size = 8, family = "Arial"))  + scale_y_continuous(name = expression(paste("Travel speed (",ms^{-1},")"))) +
  scale_colour_gradient(name = expression(paste("Wind speed (",ms^{-1},")")),low="blue",high="red")
ggsave(paste(figLoc,"groundWindSpeed.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")
# mean est wind speeds for quadrant

tstr <- lm(spTrav ~ WSpd, data = WindDat)
summary(tstr)

ggplot(WindDat[WindDat$distTo < 50,], aes(x = RelHead, y = distTo, colour = yrID)) +
  geom_point() + scale_x_continuous(name = "Relative wind heading", breaks =c(180,-90,0,90)) + coord_polar() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
      family = "Arial"), axis.text = element_text(size = 8, family = "Arial"))  + scale_y_continuous(name = expression(paste("Ground speed (",ms^{-1},")"))) +
  scale_colour_gradient(name = expression(paste("Wind speed (",ms^{-1},")")),low="blue",high="red")
# ggsave(paste(figLoc,"groundWindSpeed.svg",sep=""), device="svg", dpi = 300, height = 5,
      # width = 5, units = "in")

inters <- seq(52.5, 2.5, by=-1)
# create a circular class from RelHead
avHCirc <- circular(WindDat$fminHd, unit = "radians",rotation="counter")

avDistEst <- bind_rows(lapply(inters, function(x) circ.disp(avHCirc[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5)])))
avDistEst$p.value <- unlist(lapply(inters, function(x) r.test(avHCirc[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5)])$p.value))


tstr <- WindDat %>% group_by(yrID) %>% dplyr::summarise(sapply(inters, function(x) circ.disp(RelHead[distTo > (x - 2.5 & distTo <= (x + 2.5))])))

tstr[1,]


ggplot(data.frame(inters,avDistEst), aes(x = inters, y = rbar)) + geom_point()


ggplot(WindDat, aes(x=RelHead,y=WSpd,colour=distFk>10)) + geom_point()+coord_polar(start=pi)

WindDat$diffOf <- pi - abs(WindDat$RelHead)
WindDat$diffOf[WindDat$diffOf < ]
ggplot(WindDat[WindDat$distTo < 25,], aes(y = (pi-abs(RelHead))*(180/pi), x = distTo)) + geom_point()

# produce a figure to illustrate the wind estimation method
WindDat %>% group_by(yrID) %>% dplyr::summarise(n=n())
plot(WindDat$DT[WindDat$yrID == "7_S12018" & format(WindDat$DT,"%Y-%m-%d") == as.POSIXct("2018-09-01",tz="")],WindDat$tripL[WindDat$yrID == "7_S12018" & format(WindDat$DT,"%Y-%m-%d") == as.POSIXct("2018-09-01",tz="")])

WindDat$DT[WindDat$yrID == "7_S12018" & WindDat$DT > as.POSIXct("2018-09-01 11:00:00",tz="") & WindDat$DT < as.POSIXct("2018-09-01 12:00:00",tz="")]

# run the minute detection method for 7_S12018 for center at 3118
egDat <- data.frame(spd=r,hed=d)
ggplot(egDat, aes(x = cos(hed), y = sin(hed))) + geom_point() + coord_polar() + scale_x_continuous(limits=c(-pi,pi)) +
  geom_segment(aes(x = meangd, xend = meangd, y=0, yend = mean(spd)), colour = "black", arrow=arrow(length=unit(.5,'cm'))) +
  geom_segment(aes(x = WindDat$Head[6762], xend = WindDat$Head[6762], y = 0, yend = WindDat$WSpd[6762]), colour = "blue") +
  geom_segment(aes(x = atan2(Y[end]-Y[st],X[end]-X[st]), xend = atan2(Y[end]-Y[st],X[end]-X[st]),
    y = 0, yend = mean(spd)))


# Using raster
ggplot(egDat, aes(x=hed, y=spd) ) +
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(-pi,pi)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  ) + coord_polar()

ggplot() +
  geom_sf(data=japan,fill="#3849B5") + 
  theme(panel.background = element_rect("#F1F1EF"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank()) +
  scale_x_continuous("",limits = c(132,145)) +
  scale_y_continuous("",limits = c(30,45)) +
  geom_point(data=FkOshi, aes(x=Long,y=Lat),pch=21,fill="#E82237",size = 5)
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Conferences/2021Funder/FkOshimap.svg", dpi = 300,height=8.81,units="cm")

################################################################################
########################### FIGURE FOR MEXT SLIDE ##############################
################################################################################

FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
ggplot() + geom_sf(data = japan, fill = 'grey64', colour = 'grey64') +
    coord_sf(xlim = c(139, 147), ylim = c(39, 45)) +
    geom_path(data = allD,aes(x=lon,y=lat),size=.1,colour="grey40") +
  geom_point(data = allD[allD$forage==1,],aes(x=lon,y=lat,colour='#CC3300'), pch=16,size=2) + 
  scale_x_continuous(name="Longitude", breaks=seq(139,147,1)) +
  scale_y_continuous(name="Latitude", breaks=seq(39,45,1)) +
  geom_point(data = FkOshi, aes(x = Long, y = Lat, colour = "chartreuse3"), fill = "chartreuse3", pch = 24,size = 3) +
  scale_colour_manual(name = "", values = c("#CC3300","black"), labels=c("Nest colony","Foraging points")) +
  scale_shape_manual(name = "", values = c(24,16)) +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12),
    axis.text = element_text(size = 10),legend.position="top",
    legend.text = element_text(size = 12)) + 
    guides(colour = guide_legend(override.aes = list(size=3,shape=c(24,16),colour=c("black","#CC3300"))))
    #annotation_scale(location = 'br')
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Figures/Grants/ExampleForageAll.svg", device="svg", dpi = 300, height = 10,
      width = 10, units = "in")

###############################################################################
######################## INDIVIDUAL TRACK AND WINDS ###########################
###############################################################################

tmp50 <- WindDat[WindDat$distTo < 50,] %>% group_by(yrID,forNo) %>%
  summarise(wind_num = n()) %>% filter(wind_num == max(wind_num))

tmp50
# take example of tag b
b = 11
sub = WindDat[WindDat$forNo == tmp50$forNo[b] & WindDat$distTo < 50 & WindDat$yrID == tmp50$yrID[b],]
nxtForPt = min(which(allD$forage == 1 & allD$DT > sub$DT[1] & allD$yrID == tmp50$yrID[b]))

sub = na.omit(sub)
ggplot() +
  geom_spoke(data = sub,
    aes(x = UTME, y = UTMN, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), alpha = .6,
    radius = 2000*sub$WSpeed/max(sub$WSpeed,na.rm=T)) +
  geom_point(data=allD[nxtForPt,],aes(x=UTME,y=UTMN,fill="deepskyblue"),pch=21,size=3) +
  geom_path(data=allD[allD$DT > (sub$DT[1]-1800) & allD$DT <= allD$DT[nxtForPt] & allD$yrID == tmp50$yrID[b],],aes(x=UTME,y=UTMN)) +
  scale_x_continuous(name="Easting") +
  scale_y_continuous(name="Northing") +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial"))

ggplot() + 
  geom_quiver(data=sub,aes(x=lon,y=lat,u=WSpd*cos(WHead),v=WSpd*sin(WHead)))


ggplot() + 
  geom_point(data=sub,aes(x=lon,y=lat,colour=DT))

ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(142.3, 142.7), ylim = c(39.5, 39.75)) +
    geom_path(data = allD[allD$yrID == indivWinds[a] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),],aes(x=lon,y=lat)) +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[a] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00"),],
    aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpd[WindDat$yrID == indivWinds[a] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00")]/max(WindDat$WSpd[WindDat$yrID == indivWinds[a]]))) +  
  geom_point(data = allD[allD$yrID == indivWinds[a] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00") & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=2) +
  # geom_segment(data = allD[allD$yrID==indivWinds[a],],aes(x=lon[seq(1,nrow(sbst)-1,50)],xend=lon[seq(2,nrow(sbst),50)],
    # y=lat[seq(1,nrow(sbst)-1,50)],yend=lat[seq(2,nrow(sbst),50)]),arrow = arrow(length = unit(0.08,"inches"))) +
  scale_x_continuous(name="Lon", breaks=seq(142.3,142.7,.2)) +
  scale_y_continuous(name="Lat", breaks=seq(39.5,39.7,.2)) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + 
    annotation_scale(location = 'br')


ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    geom_path(data=allD[allD$yrID == tmp50$yrid[1],], aes(x = lon, y = lat)) + geom_point(data=allD[allD$yrID == tmp50$yrid[1] & allD$forage == 1 & allD$DT > as.POSIXct("2018/09/08 10:20:00") & allD$DT > as.POSIXct("2018/09/08 10:30:00"),],
    aes(x = lon, y = lat), pch = 21, fill = "deepskyblue") +
  geom_spoke(data = WindDat[WindDat$yrID == tmp50$yrid[1],], aes(x = lon, y = lat, colour = WSpd, angle = WHead), arrow = arrow(length = unit(0.05,"inches")),
  radius = .5*(WindDat$WSpd[WindDat$yrID == tmp50$yrid[1]]/max(WindDat$WSpd[WindDat$yrID == tmp50$yrid[1]]))) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  geom_segment(aes(x=allD[allD$yrID == tmp50$yrid[1],]$lon[seq(1,nrow(allD[allD$yrID == tmp50$yrid[1],])-1,200)],xend=allD[allD$yrID == tmp50$yrid[1],]$lon[seq(2,nrow(allD[allD$yrID == tmp50$yrid[1],]),200)],
    y=allD[allD$yrID == tmp50$yrid[1],]$lat[seq(1,nrow(allD[allD$yrID == tmp50$yrid[1],])-1,200)],yend=allD[allD$yrID == tmp50$yrid[1],]$lat[seq(2,nrow(allD[allD$yrID == tmp50$yrid[1],]),200)]),arrow = arrow(length = unit(0.1,"inches"))) +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 8,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
    annotation_scale(location = 'br') +
    coord_sf(xlim = c(139, 145), ylim = c(39, 43)) +
    scale_y_continuous(breaks = c(39,41,43), labels = c("39","41","43"), name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(breaks=c(139,141,143,145),labels = c("139", "141", "143", "145"), name = paste("Longitude (","\u00b0E",")", sep = ""))



# extra figures for AORI presentation
sbst = allD[allD$yrID == indivWinds[5] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),]
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(142.3, 142.7), ylim = c(39.5, 39.75)) +
    geom_path(data = allD[allD$yrID == indivWinds[5] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),],aes(x=lon,y=lat)) +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[5] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00"),],
    aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.05,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpeed[WindDat$yrID == indivWinds[5] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00")]/max(WindDat$WSpeed[WindDat$yrID == indivWinds[5]])),size=1.5) +  
  geom_point(data = allD[allD$yrID == indivWinds[5] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00") & allD$forage == 1,],aes(x=lon,y=lat, fill="deepskyblue"), pch=21,size=2) +
  geom_segment(aes(x=sbst$lon[seq(1,nrow(sbst)-1,50)],xend=sbst$lon[seq(2,nrow(sbst),50)],
    y=sbst$lat[seq(1,nrow(sbst)-1,50)],yend=sbst$lat[seq(2,nrow(sbst),50)]),arrow = arrow(length = unit(0.08,"inches"))) +
  scale_x_continuous(name="Lon", breaks=seq(142.3,142.7,.2)) +
  scale_y_continuous(name="Lat", breaks=seq(39.5,39.7,.2)) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "GnBu") +
  scale_fill_manual(name = "", values = c("deepskyblue",""), labels=c("Foraging points","")) +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10),
    axis.text = element_text(size = 10)) 

ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Admin/AORIPresentation/Animation/NearWind.png", device = "png", dpi = 300, height = 5,
    width = 5, units = "in")





WindDat$WSpd <- sqrt(WindDat$X^2 + WindDat$Y^2)
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
  geom_spoke(data = WindDat, aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.03,"inches")),
  radius = .4*(WindDat$WSpd/max(WindDat$WSpd)),size=.7) +
  coord_sf(xlim = c(141, 147), ylim = c(39, 44)) +
  scale_colour_viridis(name="Wind Speed (m/s)", option = "magma") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8)) + 
    annotation_scale(location = 'br') +
    scale_y_continuous(name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(name = paste("Longitude (","\u00b0E",")", sep = ""))
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Admin/AORIPresentation/Animation/ALlWind.png", device = "png", dpi = 300, height = 5,
    width = 5, units = "in")
  

allPts <- data.frame(x = allD$lon,y=allD$lat)
ch <- chull(allPts)
coords <- allPts[c(ch,ch[1]),]

FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
ggplot() + 
  geom_polygon(data=coords,aes(x=x,y=y,fill = "brown3"),alpha=.4) +
  geom_sf(data = japan, colour = '#969696',alpha=1) +
  geom_point(data = FkOshi, aes(x = Long, y = Lat,colour = "black"),fill="chartreuse3", pch = 24,size = 4) +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_fill_manual(name="",labels="Area usage",values="brown3") +
  scale_colour_manual(name="",labels='Colony',values="black") +
  guides(colour = guide_legend(override.aes = list(shape=24,colour="black"))) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8),legend.text = element_text(size = 12)) + 
  annotation_scale(location = 'br') +
  scale_y_continuous(name = paste("Latitude (","\u00b0N",")", sep = "")) +
  scale_x_continuous(name = paste("Longitude (","\u00b0E",")", sep = ""))
ggsave("/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Admin/AORIPresentation/Animation/Area.png", device = "png", dpi = 300, height = 5,
    width = 5, units = "in")

ggplot() + 
    geom_polygon(data=birdPath,aes(x=x,y=y),fill="") +
    geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    geom_spoke(data = WindDat, aes(x = lon, y = lat, colour = WSpeed, angle = WHead), arrow = arrow(length = unit(0.03,"inches")),
    radius = .4*(WindDat$WSpd/max(WindDat$WSpd)),size=.7) +
    coord_sf(xlim = c(141, 147), ylim = c(39, 44)) +
    scale_colour_viridis(name="Wind Speed (m/s)", option = "magma") +
    theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10), axis.text = element_text(size = 8)) + 
    annotation_scale(location = 'br') +
    scale_y_continuous(name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(name = paste("Longitude (","\u00b0E",")", sep = ""))

########################################################################
###################### LONG AND SHORT WIND USE #########################
########################################################################

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
    stat_density(aes(x = aligned,colour = bin5), geom = "line",
      position = "identity", size = 1.1) +
    scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=gsub(",",":",gsub('[[)]',"",sort(unique(uniqTr$bin5[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == T]))))) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10),
    axis.text = element_text(size = 8)) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_x_continuous(name = "", breaks = c(-pi,0,pi), labels = c("Tail","Head","Tail"), limits = c(-pi,pi))  +
  scale_y_continuous(name="Proportion across all birds (%)", breaks=seq(0,0.9,.1),labels=seq(0,90,10),limits=c(0,.85))+
  annotate("text",x = -3, y = .85, label = "Long", hjust=0)

# ggplot(uniqTr[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == T,]) +
#     geom_point(aes(x = mnHead,y = dist, fill = rbar), pch = 21, size = 3) +
#     coord_polar()

short <- ggplot(uniqTr[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == F,]) +
    stat_density(aes(x = aligned,colour = bin5),position = "identity", geom = "line", size = 1.1) +
    scale_colour_viridis(name="Distance to next \nforaging spot (km)", discrete = T,
    labels=gsub(",",":",gsub('[[)]',"",sort(unique(uniqTr$bin5[uniqTr$dist < 50 & uniqTr$pval < 0.05 & uniqTr$tripL == F]))))) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10),
    axis.text = element_text(size = 8)) + 
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_x_continuous(name = "", breaks = c(-pi,0,pi), labels = c("Head","Tail","Tail"), limits = c(-pi,pi))  +
  scale_y_continuous(name="", breaks=seq(0,0.9,.1),labels=seq(0,90,10),limits=c(0,.85))+
  annotate("text", x = -3, y = .85, label = "Short",hjust=0)

combdist5 <- ggarrange(long,short,nrow=2,common.legend = T,legend = "right")
annotate_figure(combdist5,bottom=textGrob("Relative wind heading"))
ggsave(paste(figLoc,"ShortLongGrouped.svg",sep=""), device="svg", dpi = 300, height = 9,
      width = 11.4, units = "cm")

ggarrange(fig, combdist5)
