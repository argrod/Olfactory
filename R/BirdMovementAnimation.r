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
install.packages("magick")

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
library(rgeos)
library(magick)
library(RColorBrewer)

#################################################################################
######################## BRING IN THE FORAGING ESTIMATES ########################
#################################################################################

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
    recalSp = c(D18$recalSp, D19$recalSp),
    distFk = c(D18$distFromFk, D19$distFromFk),
    tripL = c(D18$tripL, D19$tripL),
    tkb = c(D18$tkb, D19$tkb),
    dv = c(D18$dv, D19$dv))
allD$Year <- format(allD$DT, format = "%Y")
allD$forage <- allD$dv == 1 | allD$tkb == 1
allD$yrID <- paste(format(allD$DT,'%Y'),"_",sub("\\_S.*","",allD$tagID),sep="")
japan <- ne_countries(scale = "medium", country = "Japan", returnclass = "sf")

# prepare data (2018)
selD <- allD[format(allD$DT,"%Y") == "2018",]
# reorder by time
selD <- selD[order(selD$DT),]

# prepare n colours (based on number of individuals in data)
n <- length(unique(selD$tagID))
IndCols <- brewer.pal(length(unique(selD$tagID)), name = "Paired")
names(IndCols) <- levels(as.factor(selD$tagID))
colScale <- scale_colour_manual(name = "tagID", values = IndCols)

p1 = ggplot(selD) +
    geom_point(aes(x = lon, y = lat, colour = tagID), alpha = 0) +
    scale_colour_manual(name = "Tag ID", values = IndCols, drop = F)+
    guides(colour = guide_legend(override.aes = list(alpha=1))) +
    geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(140, 146), ylim = c(39, 44)) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10)) +
    scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
    
if(Sys.info()['sysname'] == "Darwin"){
    dir_out <- "/Volumes/GoogleDrive/My Drive/PhD/Figures/Visuals/2018Tracks/"
} else {
    dir_out <- "E:/My Drive/PhD/Figures/Visuals/2018Tracks/"
}

tseq <- seq(selD$DT[1], selD$DT[nrow(selD)],by = "hour")
for (b in 1:length(tseq)) {
    # select data within 5 mins
    subD <- selD[selD$DT > (tseq[b] - minutes(60)) & selD$DT <= tseq[b],]
    # number of indivs
    Indivs <- unique(subD$tagID)
    adding <- p1
    for (ind in 1:length(Indivs)) {
        # select only that indiv
        indD <- subD[subD$tagID == Indivs[ind],]
        # apply alpha gradient by time difference from original time
        agrad <- 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/(60*60))
        adding <- adding + geom_point(data = indD, aes(x = lon, y = lat, group = tagID,
        colour = tagID), alpha = agrad, size = 3)+ colScale + 
        scale_colour_manual(name = "Tag ID", values = IndCols, drop = F)
    }
    adding <- adding + annotate(geom="text", x = 140, y = 44, label = as.character(tseq[b]), hjust = 0)
    fp = file.path(dir_out, paste0(as.character(b),".png"))

    ggsave(plot = adding,
        filename = fp,
        device = "png")
}

imgs <- list.files(dir_out, full.names = T)
imgList <- lapply(imgs, image_read)

imgJoined <- image_join(imgList)

imgAnim <- image_animate(imgJoined, fps = 10)

image_write(image = imgAnim, path = paste0(dir_out,"2018Move.gif"))

adding + geom_point(data = indD, aes(x = lon, y = lat,
            colour = as.character(IndCols$Cols[IndCols$Inds == ind]))),
       alpha = agrad), size = 2, pch = 1) +
       theme(legend.position = "none")


ggplot() + geom_point(data = indD, aes(x = lon, y = lat))
ggplot() + geom_point(data = indD, aes(x = lon, y = lat,
    colour = as.character(IndCols$Cols[IndCols$Inds == ind]),
    alpha = agrad))


ind=Indivs[1]
 indD <- subD[subD$tagID == ind,]
        # apply alpha gradient by time difference from original time
        agrad <- 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/(30*60))
tst = p1 + geom_point(data = indD, aes(x = lon, y = lat,
    colour = as.character(IndCols$Cols[IndCols$Inds == ind]),
    alpha = agrad)) + theme(legend.position="none")

ind=Indivs[2]
 indD <- subD[subD$tagID == ind,]
        # apply alpha gradient by time difference from original time
        agrad <- 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/(30*60))
tst = tst + geom_point(data = indD, aes(x = lon, y = lat,
    colour = as.character(IndCols$Cols[IndCols$Inds == ind]),
    alpha = agrad)) + theme(legend.position="none")

    tst +(geom_point(aes(x = indD$lon, y = indD$lat,
    colour = as.character(IndCols$Cols[IndCols$Inds == ind]),
    alpha = agrad)) + theme(legend.position="none"))

adding +(geom_point(data = indD, aes(x = lon, y = lat,
            colour = as.character(IndCols$Cols[IndCols$Inds == ind]),
       alpha = 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/30*60)), size = 2, pch = 1) +
       theme(legend.position = "none"))

################################################################################
######################### ANIMATE WITH FORAGING POINTS #########################
################################################################################


# prepare data (2018)
selD <- allD[format(allD$DT,"%Y") == "2018",]
# reorder by time
selD <- selD[order(selD$DT),]

# prepare n colours (based on number of individuals in data)
n <- length(unique(selD$tagID))
IndCols <- brewer.pal(length(unique(selD$tagID)), name = "Paired")
names(IndCols) <- levels(as.factor(selD$tagID))
colScale <- scale_colour_manual(name = "tagID", values = IndCols, drop = F)

FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
p1 = ggplot() +
    geom_sf(data = japan, fill = '#969696', colour = '#969696') +
    coord_sf(xlim = c(140, 146), ylim = c(39, 44)) + 
    geom_point(data=selD,aes(x = lon, y = lat,fill="deepskyblue"),colour="deepskyblue", alpha = 0) +
    geom_point(data = FkOshi, aes(x = Long, y = Lat, fill = "chartreuse3"), colour = "black", pch = 24,size = 3) +
    geom_point(data=selD[selD$forage==1,],aes(x=lon,y=lat,fill="#CC3300"),colour="#CC3300",alpha=0,stroke=2)+
    scale_fill_manual(name = "", values=c("deepskyblue","chartreuse3","#CC3300"),labels = c('Nest colony',"Kevin","Foraging"), drop = F)+
    guides(fill = guide_legend(override.aes = list(colour=c("black","deepskyblue","#CC3300"),fill=c("chartreuse3","deepskyblue","white"),stroke=c(1,1,2),pch=c(24,16,1)))) +
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(panel.border = element_rect(colour = 'black', fill = NA), text = element_text(size = 10)) +
    scale_x_continuous("Longitude") + scale_y_continuous("Latitude")
p1
if(Sys.info()['sysname'] == "Darwin"){
    dir_out <- "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Admin/AORIPresentation/Animation/WithPrev/"
} else {
    dir_out <- "E:/My Drive/PhD/Admin/AORIPresentation/Animation/WithPrev/"
}

selD <- allD[allD$yrID == "2018_6",]
tseq <- seq(selD$DT[1], selD$DT[nrow(selD)],by = "30 mins")
selD$forage[is.na(selD$forage)] = 0
for(b in 1:length(tseq)){
    adding <- p1
    # add previous foraging points and track, comment out if not needed
    adding <- adding + geom_path(data=selD[selD$DT <= tseq[b],],aes(x=lon,y=lat),size=.3,alpha=.5) +
        geom_point(data=selD[selD$DT <= tseq[b] & selD$forage == 1,],aes(x=lon,y=lat),pch=1,colour="#CC3300",size=1,stroke=1.2)
    agrad <- 1 + as.numeric(difftime(selD$DT[selD$DT > (tseq[b] - minutes(30)) & selD$DT <= tseq[b]],max(selD$DT[selD$DT <= tseq[b]]),units='secs'))/(3600)
        if(length(agrad) > 1){
            agrad[2:length(agrad)] <- agrad[2:length(agrad)] - .2
        }
        # agrad <- 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/(60*60))
        # select only that indiv
        # indD <- subD[subD$tagID == Indivs[ind],]
        # apply alpha gradient by time difference from original time
        # agrad <- 1-(as.numeric(difftime(selD$DT[b],indD$DT,units="secs"))/(60*60))
        adding <- adding + geom_point(data = selD[selD$DT > (tseq[b] - minutes(30)) & selD$DT <= tseq[b],],
            aes(x = lon, y = lat),colour = 'deepskyblue', alpha = agrad, size = 3)
        # add dot if foraging
    if(any(selD$forage[selD$DT > (tseq[b] - minutes(30)) & selD$DT <= tseq[b]]==1)){
        adding <- adding + geom_point(data=selD[selD$DT > (tseq[b] - minutes(30)) & selD$DT <= tseq[b] & selD$forage == 1,],
            aes(x = lon, y = lat),pch=1, colour = "#CC3300", size = 2,alpha = .7,stroke = 2)
    }
    if(b > 1){
        # add expanding dot for previous foraging
        if(any(selD$forage[selD$DT > (tseq[b-1] - minutes(30)) & selD$DT <= tseq[b-1]]==1)){
            adding <- adding + geom_point(data=selD[selD$DT > (tseq[b-1] - minutes(30)) & selD$DT <= tseq[b-1] & selD$forage == 1,],
                aes(x = lon, y = lat),pch=1,colour = '#CC3300', size = 3,alpha = .4,stroke = 2)
        }
    }
    if(b > 2){
        # add expanding dot for previous foraging
        if(any(selD$forage[selD$DT > (tseq[b-2] - minutes(30)) & selD$DT <= tseq[b-2]]==1)){
            adding <- adding + geom_point(data=selD[selD$DT > (tseq[b-2] - minutes(30)) & selD$DT <= tseq[b-2] & selD$forage == 1,],
                aes(x = lon, y = lat),pch=1, colour = "#CC3300", size = 4,alpha = .3,stroke = 2)
        }
    }
    if(b > 3){
        # add expanding dot for previous foraging
        if(any(selD$forage[selD$DT > (tseq[b-3] - minutes(30)) & selD$DT <= tseq[b-3]]==1)){
            adding <- adding + geom_point(data=selD[selD$DT > (tseq[b-3] - minutes(30)) & selD$DT <= tseq[b-3] & selD$forage == 1,],
                aes(x = lon, y = lat),pch=1, colour = "#CC3300", size = 5,alpha = .2,stroke = 2)
        }
    }
    if(b > 4){
        # add expanding dot for previous foraging
        if(any(selD$forage[selD$DT > (tseq[b-4] - minutes(30)) & selD$DT <= tseq[b-4]]==1)){
            adding <- adding + geom_point(data=selD[selD$DT > (tseq[b-4] - minutes(30)) & selD$DT <= tseq[b-4] & selD$forage == 1,],
                aes(x = lon, y = lat),pch=1, colour = "#CC3300", size = 6,alpha = .1,stroke =2)
        }
    }
    adding <- adding + annotate(geom="text", x = 140, y = 44, label = as.character(tseq[b]), hjust = 0)
    fp = file.path(dir_out, paste0(as.character(b),".png"))

    ggsave(plot = adding,
        filename = fp, width = 5, height = 5, units = "in",
        device = "png")

    rm(adding)
}
