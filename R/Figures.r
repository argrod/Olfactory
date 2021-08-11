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
library(circular)
library(wGribDat)
# devtools::install_github("keesmulder/circglmbayes")
library(circglmbayes)
library(lme4)
# library(ncdf4)
# install.packages("CircMLE")
library(CircStats)
library(Cairo)
library(extrafont)
library(CircMLE)
library(ggspatial)
library(png)
library(ggspatial)
options(timeout = 800)

if(Sys.info()['sysname'] == "Darwin"){
    # load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "/Volumes/GoogleDrive/My Drive/PhD/Manuscripts/BehaviourIdentification/Figures/"
} else {
    # load("F:/UTokyoDrive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("F:/UTokyoDrive/PhD/Data/20182019AnalysisDat.RData")
    outloc <- "F:/UTokyoDrive/PhD/Manuscripts/BehaviourIdentification/Figures/"
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
# allD$forBeh <- NA
# allD$forBeh[allD$dv == 1] <- "Dive"
# allD$forBeh[allD$tkb == 1] <- "Surf"
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")
# distances <- allD %>% filter(forage==1) %>% group_by(tagID, forBeh, Year) %>% dplyr::summarise(mxdist=median(distFk))
# distances %>% group_by(forBeh,Year) %>% dplyr::summarise(mean(mxdist))
# FORAGING SPOT DISPERSALS
# if(Sys.info()['sysname'] == "Darwin"){
#     load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ForageDisps.RData")
# } else {
#     load("F:/UTokyoDrive/PhD/Data/splitr/ForageDisps.RData")
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

# LOAD IN THE STEP LENGTHS TRAJECTORIES
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/StepsTrajTimeChgNew.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/StepsTrajTimeChgNew.RData")
}
# save(ListD, file="/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ListD.RData")
# LOAD IN THE LISTED DATA
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ListD.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/ListD.RData")
}
# LOAD WIND DATA
if(Sys.info()['sysname'] == "Darwin"){
  load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDatAll.RData")
} else {
  load("F:/UTokyoDrive/PhD/Data/WindCalc/windDatAll.RData")
}
allTraj <- bind_rows(outTraj)
allTraj$relH <- allTraj$aveHd - allTraj$trjHd

# figure locations
if(Sys.info()['sysname'] == "Darwin"){
  figLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Figures/Olfactory/"
  # figLoc <- "/Documents/GitHub/PhD/Olfactory/"
} else {
  figLoc <- "F:/UTokyoDrive/PhD/Figures/Olfactory/"
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

one2Ten <- vector(mode="list", length = length(distGaps))
distGaps <- seq(0,40,10)
distGapsL <- distGaps+10
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
save(pvals,avRelHd,one2Ten, file="F:/UTokyoDrive/PhD/Data/WindCalc/headings.RData")

distGaps <- seq(0,190,10)
distGapsL <- distGaps+10
avRelHd <- NA
for(b in 1:length(distGaps)){
    RaylT <- r.test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    # tst<-HR_test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    if(RaylT$p.value < 0.05){
      avRelHd[b] <- circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    } else {
      avRelHd[b] <- NA
    }
}
# apply(distGaps, function(x) watson.wheeler.test(WindDat$RelHead[WindDat$distTo >= x & WindDat$distTo < (x+1)]))

ggplot(WindDat, aes(x = spTrav, y = WSpd)) + geom_point() +# coord_polar(start = pi) +
  # scale_colour_gradient(name = expression(paste("Wind speed (",ms^{-1},")")),low="blue",high="red") +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Ground speed (",ms^{-1},")",sep="")))

ggplot(WindDat[WindDat$rtChg < 0 & WindDat$distFromFk < 100,], aes(x = fminRelHd, y = WSpd)) + geom_point(pch =21, fill = "red") + coord_polar(start = pi) +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Estimated wind speeds (",ms^{-1},")",sep="")))
ggsave(paste(figLoc,"SpeedAngles.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")
unique(WindDat$yrID[WindDat$rtChg < 0 & WindDat$distFromFk < 100])
ggplot(WindDat[WindDat$distTo < 50,], aes(x = fminRelHd, y = distTo)) + geom_point(pch =21, fill = "red") + coord_polar(start = pi) +
  scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,family = "Arial"),
    axis.text = element_text(size = 8, family = "Arial")) + scale_y_continuous(name=expression(paste("Estimated wind speeds (",ms^{-1},")",sep="")))

ggplot(allListD[allListD$Forage == 1 & allListD$rtChg < 0 & allListD$tripL > 2,], aes(x = distFromFk)) + 
  geom_histogram(bins=50) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial")) +
  scale_x_continuous(name = "Distance from nest colony (km)") +
  scale_y_continuous("Count")

ggplot(WindDat[WindDat$distTo > 50,], aes(x = fminRelHd, y = WSpd)) + geom_point() + coord_polar(start=pi)
ggplot(WindDat[WindDat$DT < as.POSIXct("2019/01/01"),], aes(x = DT, y = WSpd)) + geom_point()


WindDat$fminRelHd <- WindDat$fminHd - WindDat$WHd
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

lengths<-seq(from=0,to=max(WindDat$distTo)-1, by = 1)
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
leavelengths<-seq(from=0,to=max(WindDat$distFromFk)-1, by = 5)
leavelengthsL <- leavelengths+5
leavingshrt <- WindDat[WindDat$rtChg < 0 & WindDat$tripL <= 2,]
leaveShrtDat <- data.frame(dist = rep(NA,length(leavelengths)), aveHd = rep(NA,length(leavelengths)), disp = rep(NA,length(leavelengths)), hrP = rep(NA,length(leavelengths)))
for(b in 1:length(leavelengths)){
  if(sum(leavingshrt$distFromFk >= leavelengths[b] & leavingshrt$distFromFk < leavelengthsL[b]) > 20){
    leaveShrtDat$dist[b] = leavelengths[b]
    leaveShrtDat$aveHd[b] = circ.mean(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[b] & leavingshrt$distFromFk < leavelengthsL[b]])
    leaveShrtDat$disp[b] = circ.disp(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[b] & leavingshrt$distFromFk < leavelengthsL[b]])$rbar
    leaveShrtDat$hrP[b] = as.numeric(HR_test(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[b] & leavingshrt$distFromFk < leavelengthsL[b]])[2])
  }
}
save(leaveShrtDat,file="/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/shrtLeaveDat.RData")
leaveShrtDat <- leaveShrtDat[!is.na(leaveShrtDat$dist),]
# leaveShrtDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[x] & leavingshrt$distFromFk < leavelengthsL[x]]))),
  # disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[x] & leavingshrt$distFromFk < leavelengthsL[x]])$var)),
  # uniP = unlist(lapply(1:length(leavelengths), function(x) HR_test(leavingshrt$RelHead[leavingshrt$distFromFk >= leavelengths[x] & leavingshrt$distFromFk < leavelengthsL[x]])[2])))
leavinglng <- WindDat[WindDat$rtChg < 0 & WindDat$tripL > 2,]
leaveLngDat <- data.frame(dist = rep(NA,length(leavelengths)), aveHd = rep(NA,length(leavelengths)), disp = rep(NA,length(leavelengths)), hrP = rep(NA,length(leavelengths)))
for(b in 1:length(leavelengths)){
  if(sum(leavinglng$distFromFk >= leavelengths[b] & leavinglng$distFromFk < leavelengthsL[b]) > 20){
    leaveLngDat$dist[b] = leavelengths[b]
    leaveLngDat$aveHd[b] = circ.mean(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[b] & leavinglng$distFromFk < leavelengthsL[b]])
    leaveLngDat$disp[b] = circ.disp(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[b] & leavinglng$distFromFk < leavelengthsL[b]])$rbar
    leaveLngDat$hrP[b] = as.numeric(HR_test(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[b] & leavinglng$distFromFk < leavelengthsL[b]])[2])
  }
}
# leaveLngDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[x] & leavinglng$distFromFk < leavelengthsL[x]]))),
#   disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[x] & leavinglng$distFromFk < leavelengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(leavelengths), function(x) HR_test(leavinglng$RelHead[leavinglng$distFromFk >= leavelengths[x] & leavinglng$distFromFk < leavelengthsL[x]])[2])))
save(leaveLngDat,file="/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/lngLeaveDat.RData")
leaveLngDat <- leaveLngDat[!is.na(leaveLngDat$disp),]
# leaving <- WindDat[WindDat$rtChg < 0 ,]
# leaveDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(leaving$RelHead[leaving$distFromFk >= leavelengths[x] & leaving$distFromFk < leavelengthsL[x]]))),
#   disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(leaving$RelHead[leaving$distFromFk >= leavelengths[x] & leaving$distFromFk < leavelengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(leavelengths), function(x) r.test(leaving$RelHead[leaving$distFromFk >= leavelengths[x] & leaving$distFromFk < leavelengthsL[x]])$p.value)))
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

ggplot(WindDat[WindDat$rtChg < 0 & WindDat$distFromFk < 200 & WindDat$tripL > 2,], aes(x = RelHead)) + geom_density(alpha=.3) + 
  coord_polar(start = pi)

lm.circular(as.circular(WindDat$RelHead), x=WindDat$distFromFk, type="c-l",verbose=T, init=c(1))
x<-cbind(rnorm(10),rep(1,10))
y<-circular(2*atan(c(x%*%c(5,1))))+rvonmises(10, mu=circular(0), kappa=100)
lm.circular(y=y, x=x, init=c(5,1), type='c-l', verbose=TRUE)

fit.Motor <- bpnr(pred.I = Phaserad ~ 1 + Cond, data = Motor, its = 10000, burn = 100, n.lag = 3, seed = 101)

fit.RelH <- bpnr(pred.I = (RelHead+pi) ~ 1 + distFromFk + (1|yrID) + WSpd, data = WindDat[WindDat$tripL > 2,], its = 10000, burn = 200, n.lag = 20, seed = 101)
traceplot(fit.RelH)

plot(WindDat$RelHead ~ WindDat$distFromFk)
WindDat <- WindDat[,-which(names(WindDat) %in% c("timeTo","forHd","spTrav","minHd"))]
WindDat$yrIDnm <- as.numeric(factor(WindDat$yrID))
testr <- bpnme(pred.I = RelHead ~ distFromFk + WSpd + (1|yrIDnm), data = WindDat[1:20,], its = 10000, burn = 1000, n.lag = 3)
WindDat$offset
tstr <- lm(offset ~ distTo + distFromFk + WSpd, data = WindDat)
summary(tstr)
plot(tstr)

fit.Maps <- bpnme(pred.I = Error.rad ~ Maze + Trial.type + L.c + (1|Subject), data = Maps, its = 10000, burn = 1000, n.lag = 3, seed = 101)

traceplot(testr)
summary(WindDat)

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
retlengths<-seq(from=max(WindDat$distFromFk)-1, to = 5, by = -5)
retlengthsL <- retlengths+5
leavingshrt <- WindDat[WindDat$rtChg < 0 & WindDat$tripL <= 2,]
retShrtDat <- data.frame(dist = rep(NA,length(retlengths)), aveHd = rep(NA,length(retlengths)), disp = rep(NA,length(retlengths)), hrP = rep(NA,length(retlengths)))
for(b in 1:length(retlengths)){
  if(sum(leavingshrt$distFromFk >= retlengths[b] & leavingshrt$distFromFk < retlengthsL[b]) > 20){
    retShrtDat$dist[b] = retlengths[b]
    retShrtDat$aveHd[b] = circ.mean(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[b] & leavingshrt$distFromFk < retlengthsL[b]])
    retShrtDat$disp[b] = circ.disp(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[b] & leavingshrt$distFromFk < retlengthsL[b]])$rbar
    retShrtDat$hrP[b] = as.numeric(HR_test(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[b] & leavingshrt$distFromFk < retlengthsL[b]])[2])
  }
}
save(retShrtDat,file="/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/shrtretDat.RData")
retShrtDat <- retShrtDat[!is.na(retShrtDat$dist),]\

# retShrtDat <- data.frame(dist = retlengths, aveHd = unlist(lapply(1:length(retlengths), function(x) circ.mean(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[x] & leavingshrt$distFromFk < retlengthsL[x]]))),
  # disp = unlist(lapply(1:length(retlengths), function(x) circ.disp(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[x] & leavingshrt$distFromFk < retlengthsL[x]])$var)),
  # uniP = unlist(lapply(1:length(retlengths), function(x) HR_test(leavingshrt$RelHead[leavingshrt$distFromFk >= retlengths[x] & leavingshrt$distFromFk < retlengthsL[x]])[2])))
leavinglng <- WindDat[WindDat$rtChg < 0 & WindDat$tripL > 2 & WindDat$distFromFk <= 200,]
retLngDat <- data.frame(dist = rep(NA,length(retlengths)), aveHd = rep(NA,length(retlengths)), disp = rep(NA,length(retlengths)), hrP = rep(NA,length(retlengths)))
for(b in 1:length(retlengths)){
  if(sum(leavinglng$distFromFk >= retlengths[b] & leavinglng$distFromFk < retlengthsL[b]) > 20){
    retLngDat$dist[b] = retlengths[b]
    retLngDat$aveHd[b] = circ.mean(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[b] & leavinglng$distFromFk < retlengthsL[b]])
    retLngDat$disp[b] = circ.disp(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[b] & leavinglng$distFromFk < retlengthsL[b]])$rbar
    retLngDat$hrP[b] = as.numeric(HR_test(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[b] & leavinglng$distFromFk < retlengthsL[b]])[2])
  }
}
# retLngDat <- data.frame(dist = retlengths, aveHd = unlist(lapply(1:length(retlengths), function(x) circ.mean(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[x] & leavinglng$distFromFk < retlengthsL[x]]))),
#   disp = unlist(lapply(1:length(retlengths), function(x) circ.disp(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[x] & leavinglng$distFromFk < retlengthsL[x]])$var)),
#   uniP = unlist(lapply(1:length(retlengths), function(x) HR_test(leavinglng$RelHead[leavinglng$distFromFk >= retlengths[x] & leavinglng$distFromFk < retlengthsL[x]])[2])))
save(retLngDat,file="/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/lngretDat.RData")
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
# ggplot(WindDat[WindDat$rtChg < 0,], aes(x = RelHead, y=distFromFk,fill=trip)) + geom_point(pch=21,colour = 'black') + coord_polar(start=pi)
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
  lengths<-seq(from=0,to=max(WindDat$distTo[WindDat$yrID == indivWinds[b]])-1, by = 1)
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


breaks<-seq(from=0,to=round_any(max(WindDat$distTo),10,f=ceiling),by=10)
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



bin10s <- bin10ns[,2]
typeof(WindDat$bin10)


sbst = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[6],]
WindDat$WSpd <- sqrt(WindDat$X^2 + WindDat$Y^2)
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
ggsave(paste(figLoc,"ExampleWind.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")

sbst = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),]
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(142.3, 142.7), ylim = c(39.5, 39.75)) +
    geom_path(data = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00"),],aes(x=lon,y=lat)) +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[a] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00"),],
    aes(x = Lon, y = Lat, colour = WSpd, angle = WHd), arrow = arrow(length = unit(0.05,"inches")), alpha = .6,
    radius = .3*(WindDat$WSpd[WindDat$yrID == indivWinds[a] & WindDat$DT > as.POSIXct("2018-08-29 05:00:00") & WindDat$DT < as.POSIXct("2018-08-29 07:00:00")]/max(WindDat$WSpd[WindDat$yrID == indivWinds[a]]))) +  
  geom_point(data = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a] & allD$DT > as.POSIXct("2018-08-29 05:00:00") & allD$DT < as.POSIXct("2018-08-29 07:00:00") & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=2) +
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
    geom_path(data = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[6],],aes(x=lon,y=lat)) +
  geom_point(data = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[6] & allD$forage == 1,],aes(x=lon,y=lat, fill=factor(forage)), pch=21,size=2) + 
  scale_x_continuous(name="Lon", breaks=seq(140,144,1)) +
  scale_y_continuous(name="Lat", breaks=seq(39,42.5,1)) +
  scale_fill_manual(name = "Foraging points", values = "deepskyblue", labels="") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + 
    annotation_scale(location = 'br')
ggsave(paste(figLoc,"ExampleForageAll.svg",sep=""), device="svg", dpi = 300, height = 6,
      width = 6, units = "in")
    
    
a <- 5
WindDat[WindDat$yrID == indivWinds[a] & WindDat$distTo < 10 & (WindDat$RelHead < -2.6 | WindDat$RelHead > 2.6),c(1:3,26)]

a=8
sbst = allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a],]
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(141, 145), ylim = c(39, 43)) + geom_path(data=allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a],], aes(x = lon, y = lat)) + geom_point(data=allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[a] & allD$forage == 1,], aes(x = lon, y = lat), pch = 21, fill = "deepskyblue") +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[a],], aes(x = Lon, y = Lat, colour = WSpd, angle = WHd), arrow = arrow(length = unit(0.05,"inches")),
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


# RELATIVE HEADINGS AS BIRDS LEAVE COLONY
one2Ten <- vector(mode="list",length=5)
distGaps <- seq(0,40,10)
distGapsL <- distGaps+10
for(b in 1:length(distGaps)){
    RaylT <- r.test(WindDat$RelHead[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b]])
    tst<-HR_test(WindDat$RelHead[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b]])
    # if(RaylT$p.value > 0.05){
    if(tst[2] > 0.05){
        one2Ten[[b]] <- ggplot(WindDat[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b],]) + 
          geom_histogram(aes(x = RelHead), colour = "black", bins = 30, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p > 0.05", sep = ""))
    } else {
        roseplt <- ggplot(WindDat[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b],]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 30, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Tail","Side"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
        one2Ten[[b]] <- roseplt + geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b]]),
          xend = circ.mean(WindDat$RelHead[WindDat$distFromFk >= distGaps[b] & WindDat$distFromFk < distGapsL[b]]),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p < ",as.character(signif.ceiling(RaylT$p.value,3)),sep=""))
        # + geom_label(aes(x = pi/4, y = max(ggplot_build(roseplt)$data[[1]]$count)), label = paste("p value = ",as.character(signif(RaylT$p.value, 3)), sep = ""))
    }
    # Cairo(width=8, height = 8, file = paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
    print(one2Ten[[b]])
    ggsave(paste(figLoc,"RelHeadLeavingHR",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""), device="svg", dpi = 300, height = 3.5,
      width = 3, units = "in")
}

# allList <- bind_rows(ListD)
# allList$levRet <- c(allList$levRet[!is.na(allList$levRet)],allList$levRet19[!is.na(allList$levRet19)])
# WindDat$levRet <- NA
# WindDat$tripL <- NA
# for(b in 1:nrow(WindDat)){
#   WindDat$levRet[b] <- allList$levRet[which(allList$DT > (WindDat$DT[b] - lubridate::seconds(30)) & allList$DT < (WindDat$DT[b] + lubridate::seconds(30)) & paste(allList$tagID,format(allList$DT,"%Y"),sep="") == WindDat$yrID[b])]
#   WindDat$tripL[b] <- allList$tripL[which(allList$DT > (WindDat$DT[b] - lubridate::seconds(30)) & allList$DT < (WindDat$DT[b] + lubridate::seconds(30)) & paste(allList$tagID,format(allList$DT,"%Y"),sep="") == WindDat$yrID[b])]
# }
# save(WindDat,file="F:/UTokyoDrive/PhD/Data/WindCalc/windDatAll.RData")

allTraj <- bind_rows(outTraj)
allTraj$distTo <- allTraj$distTo*10^-3
allTraj$relH <- allTraj$aveHd - allTraj$trjHd

allTraj$relH[allTraj$relH < pi] = allTraj$relH[allTraj$relH < pi] + 2*pi
allTraj$relH[allTraj$relH > pi] = allTraj$relH[allTraj$relH > pi] - 2*pi
allTraj <- na.omit(allTraj)
one2TenTr <- vector(mode="list",length=5)
distGaps <- seq(0,8,2)
distGapsL <- distGaps+2
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
    # print(one2Ten[[b]])
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
aveHCirc <- circular(WindDat$Head, units = "radians",rotation="counter")
lapply(inters, function(x) watson.test(relHCirc[WindDat$distTo > (x - 5) & WindDat$distTo <= x], alpha = 0.05, dist = "vonmises"))

cirDistEst <- bind_rows(lapply(inters, function(x) circ.disp(relHCirc[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 2])))
# add the 90th and 10th %ile of est wind speed
cirDistEst[,c("ten","ninety")] <- bind_rows(lapply(inters, function(x) quantile(WindDat$WSpd[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 2], probs = c(.1,.9))))
cirDistEst$meanSpd <- unlist(lapply(inters,function(x) mean(WindDat$WSpd[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5) & WindDat$tripL <= 2])))
aveHDistEst <- bind_rows(lapply(inters, function(x) circ.disp(aveHCirc[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5)])))

plot()

sapply(inters, function(x) length(unique(WindDat$yrID[WindDat$distTo > (x - 2.5) & WindDat$distTo <= (x + 2.5)])))

plot(rev(inters[cirDistEst$n>30]),rev(cirDistEst$rbar[cirDistEst$n>30]))

# dispersal of relHead vs dispersal of bird head
dispG <- ggplot() + geom_line(data=data.frame(inters, cirDistEst),aes(x = inters, y = rbar,linetype="solid")) +
  geom_line(data = data.frame(inters,aveHDistEst), aes(x = inters, y = rbar, linetype = "dotted")) +
  scale_x_reverse(name = "") + scale_y_continuous(name = expression(bar(italic(r)))) +
  scale_linetype_manual(name = "", values=c("dotted","solid"), labels=c("Average Bird\nHeadings", "Relative Wind\nHeadings")) +
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
  annotate("text", x = 50, y = 0.8, label = "A)")

ggsave(paste(figLoc,"dispDist.svg",sep=""), device="svg", dpi = 300, height = 3,
      width = 5, units = "in")



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

WindDat$fminRelH <- WindDat$fminHd-WindDat$WHd

plot(WindDat$fminRelH, WindDat$RelHead)

plot(WindDat$fminHd-WindDat$WHd,WindDat$Head-WindDat$WHd)


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

outPt <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv", sep = ",", header=T)
outPt$Time <- as.POSIXct(outPt$Time,format="%Y-%m-%dT%H:%M:%S")
cor.circular(outPt$est)

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
    annotate("text",x=-pi,y=0.5,label="corr = 0.417", hjust = 0) +
    annotate("text",x=-pi,y=0.,label="p < 4 %*% 10^{-7}", parse = T, hjust = 0) +
    annotate("text",x=-pi,y=-0.5,label="n == 152", parse = T, hjust = 0) +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
    scale_y_continuous(name='Estimated headings (rad)') + scale_x_continuous(name='JMA headings (rad)')

splm <- lm(estSpeed ~ gribSpeed, data = outPt)

spdval <- ggplotRegression(splm) +
    geom_line(data=data.frame(x=0:max(outPt$gribSpeed),y=0:max(outPt$gribSpeed)),aes(x=x,y=y),colour="red",linetype='dashed') +
    annotate("text",x=0,y=12.5,label="y = 0.41x + 1.3", hjust= 0) +
    annotate("text",x=0,y=11.5,label="p < 2 %*% 10^{-10}", parse = T, hjust= 0) +
    annotate("text",x=0,y=10.5,label=expression(paste(R^2," = 0.24")), hjust= 0) +
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
valDat <- "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv"

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
  geom_point() + scale_x_continuous(name = "Relative wind heading", breaks =c(180,-90,0,90)) + 
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


ggplot(WindDat, aes(x=RelHead,y=WSpd,colour=distFromFk>10)) + geom_point()+coord_polar(start=pi)

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
