if(Sys.info()['sysname'] == "Darwin"){
    fileloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/GPSDominantFrequencies/"
} else {
    fileloc <- "E:/My Drive/PhD/Data/GPSDominantFrequencies/"
}
freqFiles <- dir(fileloc)
freq <- vector(mode="list",length=length(freqFiles))
for(b in 1:length(freqFiles)){
	freq[[b]] <- data.frame(read.delim(paste0(fileloc,freqFiles[b]),sep=",",header=T))
	freq[[b]]$DT <- as.POSIXct(freq[[b]]$DT,format="%Y-%m-%dT%H:%M:%OS",tz="")
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

for(b in 1:length(freq)){

}
getmode(diff(freq[[b]]$DT))

library(dplyr)
allFreq <- bind_rows(freq)
colnames(allFreq) <- c("ID","Time","y","x","domFreq")

b <- 1
difftime(freq[[b]]$DT,units="secs")

crwOut <- crawlWrap

prepOut <- prepData(allFreq,type="LL"