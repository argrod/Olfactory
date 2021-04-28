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
library(ncdf4)
# install.packages("CircStats")
library(CircStats)
library(Cairo)
library(extrafont)
options(timeout = 800)

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
# FORAGING SPOT DISPERSALS
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/ForageDisps.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/ForageDisps.RData")
}
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
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/splitr/StepsTrajTimeChg.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/splitr/StepsTrajTimeChg.RData")
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
  # figLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Figures/Olfactory/"
  figLoc <- "/Documents/GitHub/PhD/Olfactory/"
} else {
  # figLoc <- "F:/UTokyoDrive/PhD/Figures/Olfactory/"
  figLoc <- "F:/Documents/GitHub/PhD/Olfactory/"
}

# LESS THAN 10KM
ggplot(WindDat[WindDat$distTo < 10,]) + 
    geom_histogram(aes(x = RelHead), colour = "black", bins = 50) + 
    geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo < 10]), size = 1.3, linetype = 1) +
    coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("-180/180","-90","0","90"), limits = c(-pi,pi)) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 12,
        family = "Arial"), axis.text = element_text(size = 12, family = "Arial"))

one2Ten <- vector(mode="list",length=5)
distGaps <- seq(0,50,10)
distGapsL <- distGaps+10
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}
for(b in 1:length(distGaps)){
    RaylT <- r.test(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]])
    if(RaylT$p.value > 0.05){
        one2Ten[[b]] <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
          geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
          coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Side","Tail"), limits = c(-pi,pi)) + 
          theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
              family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p > 0.05", sep = ""))
    } else {
        roseplt <- ggplot(WindDat[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b],]) + 
            geom_histogram(aes(x = RelHead), colour = "black", bins = 50, fill = "#d9d9d9") + scale_y_continuous(name = "Count") +
            # geom_vline(xintercept = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]), linetype = 1, colour = "red") +
            coord_polar(start=pi) + scale_x_continuous(name = "Relative wind heading", breaks = c(pi,-pi/2,0,pi/2), labels = c("Head","Side","Side","Tail"), limits = c(-pi,pi)) + 
            theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
                family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
        one2Ten[[b]] <- roseplt + geom_segment(x = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          xend = circ.mean(WindDat$RelHead[WindDat$distTo >= distGaps[b] & WindDat$distTo < distGapsL[b]]),
          y = 0, yend = max(ggplot_build(roseplt)$data[[1]]$count)*RaylT$r.bar, colour = "#fd8d3c", lineend="round",
          arrow = arrow(length = unit(.25, "cm"))) +
          labs(title = paste(as.character(distGaps[b])," - ",as.character(distGapsL[b]), "km, p < ",as.character(signif.ceiling(RaylT$p.value,3)),sep=""))
        # + geom_label(aes(x = pi/4, y = max(ggplot_build(roseplt)$data[[1]]$count)), label = paste("p value = ",as.character(signif(RaylT$p.value, 3)), sep = ""))
    }
    # Cairo(width=8, height = 8, file = paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
    print(one2Ten[[b]])
    ggsave(paste(figLoc,"RelHead",as.character(distGaps[b]),"-",as.character(distGapsL[b]),".eps",sep=""), device="eps", dpi = 300, height = 3.5,
      width = 3, units = "in")
}


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

circ.mean(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20])
r.test(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20])
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

watson.two(WindDat$RelHead[WindDat$distTo >= 0 & WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40], alpha = .05, plot = T)
watson.wheeler.test(list(WindDat$RelHead[WindDat$distTo >= 0 & WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40]))

watson.two(WindDat$RelHead[WindDat$distTo < 1], WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2],plot=T,alpha=.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 1 & WindDat$distTo < 2], WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3],alpha=0.05,plot=T)
watson.two(WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3], WindDat$RelHead[WindDat$distTo >= 3 & WindDat$distTo < 4],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 2 & WindDat$distTo < 3], WindDat$RelHead[WindDat$distTo >= 4 & WindDat$distTo < 5],plot=T,alpha=0.05)

watson.two(WindDat$RelHead[WindDat$distTo < 10], WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 10 & WindDat$distTo < 20], WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30],plot=T,alpha=0.05)
watson.two(WindDat$RelHead[WindDat$distTo >= 20 & WindDat$distTo < 30], WindDat$RelHead[WindDat$distTo >= 30 & WindDat$distTo < 40],plot=T,alpha=0.05)



lengths<-seq(from=0,to=max(WindDat$distTo)-1, by = 1)
lengthsL <- lengths+1 
binDat <- data.frame(dist = lengthsL, aveHd = unlist(lapply(1:length(lengths), function(x) circ.mean(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]]))),
  disp = unlist(lapply(1:length(lengths), function(x) circ.disp(WindDat$RelHead[WindDat$distTo >= lengths[x] & WindDat$distTo < lengthsL[x]])$var)))
binDat <- binDat[!is.na(binDat$disp),]

# PLOT AVE HEADINGS IN 10KM BLOCKS AS THEY LEAVE FK ISLAND
leavelengths<-seq(from=0,to=max(WindDat$distFromFk)-1, by = 1)
leavelengthsL <- leavelengths+1
leaving <- WindDat[WindDat$RtChg < 0,]
leaveDat <- data.frame(dist = leavelengths, aveHd = unlist(lapply(1:length(leavelengths), function(x) circ.mean(WindDat$RelHead[WindDat$distFromFk >= leavelengths[x] & WindDat$distFromFk < leavelengthsL[x]]))),
  disp = unlist(lapply(1:length(leavelengths), function(x) circ.disp(WindDat$RelHead[WindDat$distFromFk >= leavelengths[x] & WindDat$distFromFk < leavelengthsL[x]])$var)),
  uniP = unlist(lapply(1:length(leavelengths), function(x) r.test(WindDat$RelHead[WindDat$distFromFk >= leavelengths[x] & WindDat$distFromFk < leavelengthsL[x]])$p.value)))
leaveDat <- leaveDat[!is.na(leaveDat$disp),]

ggplot(leaveDat[leaveDat$uniP < 0.05,], aes(x = aveHd, y = dist)) + geom_point() + coord_polar(start=pi) + scale_x_continuous(name="Average heading",breaks=c(-pi,-pi/2,0,pi/2),labels=c("Head","Side","Tail","Side"), limits=c(-pi,pi)) +
  scale_y_continuous(name="Distance from nest site (km)") + theme_bw() + theme(panel.grid.minor = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
    family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) 
ggsave(paste(figLoc,"LeavingHeadings.eps",sep=""), device="eps", dpi = 300, height = 4,
      width = 4, units = "in")
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
ggplot(allDispDist[allDispDist$dist < 40,], aes(x = dist, y = disp, colour = yrID)) + geom_point() + geom_vline(xintercept=10, linetype="dotted") +
    theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
            family = "Arial"), axis.text = element_text(size = 10, family = "Arial")) + scale_x_continuous(name="Distance to next foraging (km)") +
    scale_y_continuous(name = "Angular dispersal")



breaks<-seq(from=0,to=round_any(max(WindDat$distTo),10,f=ceiling),by=10)
mnW <- ddply(WindDat, "bin10", summarise, grp.mean=mean(aligned))
WindDat$bin10 <- cut(WindDat$distTo, breaks = breaks, include.lowest=T,right=F)
# Cairo(width=15, height = 15, file = paste(figLoc,"DistRelDensity.svg",sep=""),type="svg", bg = "transparent", dpi = 300, units="in")
ggplot(WindDat[WindDat$distTo > 0 &WindDat$distTo < 90,], aes(x = aligned, colour = bin10)) +#max(WindDat$distTo),], aes(x = aligned, colour = bin10)) +
  # geom_histogram(alpha=.2,fill=NA,position="dodge")
  geom_density(alpha=.2,show.legend=FALSE)+stat_density(aes(x=aligned, colour=bin10), geom="line",position="identity") +# coord_polar(start=pi) +
  scale_x_continuous(name = "Relative wind heading", breaks=c(-pi, -pi/2, 0, pi/2, pi), labels=c("Tail","Side","Head","Side","Tail")) + ylab("Density") + theme_bw() + theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
  scale_colour_manual(name="Distance to foraging (km)", values = rev(brewer.pal(9,"YlOrRd")), labels=c("0 - 10","10 - 20","20 - 30","30 - 40","40 - 50","50 - 60","60 - 70","70 - 80","80 - 90"))
ggsave(paste(figLoc,"DistRelDensity.eps",sep=""), device="eps", dpi = 300, height = 4,
      width = 6, units = "in")


WindDat$WSpd <- sqrt(WindDat$X^2 + WindDat$Y^2)
ggplot() + geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(140, 144), ylim = c(39, 42.5)) + geom_path(data=allD[paste(allD$tagID, allD$Year, sep = "") == indivWinds[6],], aes(x = lon, y = lat)) +
  geom_spoke(data = WindDat[WindDat$yrID == indivWinds[6],], aes(x = Lon, y = Lat, colour = WSpd, angle = WHd), arrow = arrow(length = unit(0.05,"inches")),
  radius = .5*(WindDat$WSpd[WindDat$yrID == indivWinds[6]]/max(WindDat$WSpd[WindDat$yrID == indivWinds[6]]))) +
  scale_colour_distiller(name="Wind Speed (m/s)", direction = 1, palette = "YlOrRd") +
  theme_bw() + theme(panel.grid = element_blank()) +
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 8,
        family = "Arial"), axis.text = element_text(size = 8, family = "Arial")) + 
    annotation_scale(location = 'br') +
    scale_y_continuous(breaks = c(39,40,41,42), labels = c("39","40","41","42"), name = paste("Latitude (","\u00b0N",")", sep = "")) +
    scale_x_continuous(labels = c("140", "141", "142", "143", "144"), name = paste("Longitude (","\u00b0E",")", sep = ""))
ggsave(paste(figLoc,"ExampleWind.svg",sep=""), device="svg", dpi = 300, height = 5,
      width = 5, units = "in")
