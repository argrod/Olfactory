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

plot(comp5[[2]]$wDirOG,comp5[[2]]$wDirTest)
