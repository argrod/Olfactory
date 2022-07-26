# read in yoneEstimates
if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive-102199952889875375671/My Drive/PD/BiP/YoneMethodComparison/"
} else {
    fileloc <- "G:/PD/BiP/YoneMethodComparison/"
}
files = dir(fileloc)
compDat <- vector(mode="list",length=length(files))
for(b in 1:length(files)){
    compDat[[1]] <- read.delim(paste0(fileloc,files[b]),header=T,sep=',')
    compDat[[1]]$Time <- as.POSIXct(compDat[[1]]$Time,format='%d-%b-%Y %H:%M:%S')
    compDat[[1]]$Tag <- sub("\\Comparison.*","",files[b])
}

library(dplyr)
library(tibble)
newDat <- compDat[[1]][compDat[[1]]$fs == 1,c("Time","wSpeed","wDir","aSpeed","resnorm","fs","treatment")]
for(c in 1:nrow(newDat)){
    if(any((compDat[[1]]$Time > (newDat$Time[c] - 5)) & (compDat[[1]]$Time < (newDat$Time[c] + 5)) & (compDat[[1]]$fs == 5))){
        newDat <- rbind(newDat, data.frame(Time=newDat$Time[c],compDat[[1]][(compDat[[1]]$Time > (newDat$Time[c] - 5)) & (compDat[[1]]$Time < (newDat$Time[c] + 5)) & (compDat[[1]]$fs != 1),c("wSpeed","wDir","aSpeed","resnorm","fs","treatment")]))
    } else {
        newDat <- rbind(newDat,data.frame(Time=newDat$Time[c],
            wSpeed = NA, wDir = NA, aSpeed = NA, resnorm = NA, fs = 5, treatment = NA))
    }
    if(any((compDat[[1]]$Time > (newDat$Time[c] - 5)) & (compDat[[1]]$Time < (newDat$Time[c] + 5)) & (compDat[[1]]$fs == 10))){
        newDat <- rbind(newDat, data.frame(Time=newDat$Time[c],compDat[[1]][(compDat[[1]]$Time > (newDat$Time[c] - 5)) & (compDat[[1]]$Time < (newDat$Time[c] + 5)) & (compDat[[1]]$fs != 1),c("wSpeed","wDir","aSpeed","resnorm","fs","treatment")]))
    }
}


ggplot(newDat[newDat$treatment=="5min",]) + geom_point(aes(x = Time, y = wSpeed, colour = as.factor(fs)))

(compDat[[1]]$Time > (newDat$Time[c] - 5)) & (compDat[[1]]$Time < (newDat$Time[c] + 5)) & (compDat[[1]]$fs != 5)

cor.circular(newDat$wDir[(newDat$treatment=="5min") & (newDat$fs == 1) & ],newDat$wDir[newDat$treatment=="5min" & newDat$fs == 5], test = T)

comp5 <- vector(mode="list",length=4)
freqs <- c(5,10,30,60)
allComp <- bind_rows(compDat)
for(c in 1:4){
    # split groups in fs or window
    min5 <- allComp[allComp$treatment=="5min",]
    comp5[[c]] <- data.frame(wSpOG = NA, wDirOG = NA, aSpOG = NA, resnormOG = NA, wSpTest = NA, wDirTest = NA, aSpTest = NA, resnormTest = NA)
    min51s <- which(min5$fs == 1)
    for(b in min51s){
        if(any((min5$Time > (min5$Time[b]-5)) &
            (min5$Time < (min5$Time[b]+5)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))){
                comp5[[c]] <- rbind(comp5[[c]],
                    data.frame(wSpOG = min5$wSpeed[b],wDirOG = min5$wDir[b], aSpOG = min5$aSpeed[b],resnormOG = min5$resnorm[b],wSpTest = min5$wSpeed[which((min5$Time > (min5$Time[b]-5)) &
            (min5$Time < (min5$Time[b]+5)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            wDirTest = min5$wDir[which((min5$Time > (min5$Time[b]-5)) &
            (min5$Time < (min5$Time[b]+5)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            aSpTest = min5$aSpeed[which((min5$Time > (min5$Time[b]-5)) &
            (min5$Time < (min5$Time[b]+5)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            resnormTest = min5$resnorm[which((min5$Time > (min5$Time[b]-5)) &
            (min5$Time < (min5$Time[b]+5)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))]))
            }
    }
}

compWin <- vector(mode="list",length=4)
freqs <- c(5,10,30,60)
for(c in 1:4){
    # split groups in fs or window
    min5 <- allComp[allComp$treatment=="300points",]
    compWin[[c]] <- data.frame(wSpOG = NA, wDirOG = NA, aSpOG = NA, resnormOG = NA, wSpTest = NA, wDirTest = NA, aSpTest = NA, resnormTest = NA)
    min51s <- which(min5$fs == 1)
    for(b in min51s){
        if(any((min5$Time > (min5$Time[b]-60)) &
            (min5$Time < (min5$Time[b]+60)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))){
                compWin[[c]] <- rbind(compWin[[c]],
                    data.frame(wSpOG = min5$wSpeed[b],wDirOG = min5$wDir[b], aSpOG = min5$aSpeed[b],resnormOG = min5$resnorm[b],wSpTest = min5$wSpeed[which((min5$Time > (min5$Time[b]-60)) &
            (min5$Time < (min5$Time[b]+60)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            wDirTest = min5$wDir[which((min5$Time > (min5$Time[b]-60)) &
            (min5$Time < (min5$Time[b]+60)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            aSpTest = min5$aSpeed[which((min5$Time > (min5$Time[b]-60)) &
            (min5$Time < (min5$Time[b]+60)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))],
            resnormTest = min5$resnorm[which((min5$Time > (min5$Time[b]-60)) &
            (min5$Time < (min5$Time[b]+60)) &
            (min5$fs == freqs[c]) & 
            (min5$Tag == min5$Tag[b]))]))
            }
    }
}

plot(compWin[[1]]$wDirOG,compWin[[1]]$wDirTest)
library(ggpubr)

g5 <- ggplot(comp5[[1]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g10 <- ggplot(comp5[[2]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g30 <- ggplot(comp5[[3]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g60 <- ggplot(comp5[[4]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))

ggarrange(g5,g10,g30,g60,labels=c("5","10","30","60"))

library(circular)
circTests <- vector('list',length=4)

circ_eqn <- function(x,y){
    m <- suppressWarnings(cor.circular(x,y,test=T))
    eq <- substitute(italic(cor) == a *","~~italic(P)~"="~b, 
         list(a = format(unname(m$cor), digits = 2),
              b = format(unname(m$p.value), digits = 2)))
    as.character(as.expression(eq));
}

for(b in 1:4){
    test <- cor.circular(comp5[[b]]$wDirOG,comp5[[b]]$wDirTest,test=T)
    circTests[[b]] <- ggplot() +
        geom_point(aes(x=wDirOG,y=wDirTest),comp5[[b]]) +
        geom_line(aes(x =-pi:pi,y=-pi:pi),colour="red",size=1.3,
        linetype='dashed') +
        geom_text(aes(x=-3,y=2.3),hjust=0,label=circ_eqn(comp5[[b]]$wDirOG,comp5[[b]]$wDirTest),parse=T) +
        theme_bw() +
        scale_x_continuous("Original wind direction") +
        scale_y_continuous("Subsample wind direction")
}
ggarrange(circTests[[1]],circTests[[2]],circTests[[3]],circTests[[4]],
    labels=c("5s","10s","30s","60s"), hjust=0)
ggsave("/Volumes/GoogleDrive-102199952889875375671/My Drive/PD/BiP/YoneMethodComparison/Figures/WindDir.pdf",
    device="pdf",dpi=300,height=6,width=6)

g5 <- ggplot(comp5[[1]]) +
    geom_point(aes(x = wSpOG, y = wSpTest))
g10 <- ggplot(comp5[[2]]) +
    geom_point(aes(x = wSpOG, y = wSpTest))
g30 <- ggplot(comp5[[3]]) +
    geom_point(aes(x = wSpOG, y = wSpTest))
g60 <- ggplot(comp5[[4]]) +
    geom_point(aes(x = wSpOG, y = wSpTest))

ggarrange(g5,g10,g30,g60,labels=c("5","10","30","60"))

lm_eqn <- function(x,y){
    m <- lm(y ~ x);
    eq <- substitute(italic(y) == a + b *""* italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

lmFS <- vector('list',length=4)
for(b in 1:4){
    lmFS[[b]] <- ggplot() + geom_point(aes(x = wSpOG,y=wSpTest),data=comp5[[b]]) +
        geom_smooth(aes(x=wSpOG,y=wSpTest),data=comp5[[b]],method='lm') +
        geom_text(aes(x = .1, y = 4.3), hjust=0,
        label=lm_eqn(comp5[[b]]$wSpOG,comp5[[b]]$wSpTest), parse = T) + 
        theme_bw() +
        scale_x_continuous("Original wind speed") +
        scale_y_continuous("Subsample wind speed",limits=c(0,5))
}
ggarrange(lmFS[[1]],lmFS[[2]],lmFS[[3]],lmFS[[4]],
    labels=c("5s","10s","30s","60s"),hjust=0)
ggsave("/Volumes/GoogleDrive-102199952889875375671/My Drive/PD/BiP/YoneMethodComparison/Figures/WindSpeed.pdf",
    device="pdf",dpi=300,height=6,width=6)

g5 <- ggplot(compWin[[1]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g10 <- ggplot(compWin[[2]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g30 <- ggplot(compWin[[30]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))
g60 <- ggplot(compWin[[60]]) +
    geom_point(aes(x = wDirOG, y = wDirTest))

ggarrange(g5,g10,g30,g60)