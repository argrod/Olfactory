library(lubridate)
library(sp)
library(raster)
library(circular)

if(Sys.info()['sysname'] == "Darwin"){
    windLoc <- '/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/'
	estLoc <- '/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/Comparison/'
} else {
    windLoc <- 'F:/UTokyoDrive/PhD/Data/WindEstTest/2018/'
	estLoc <- 'F:/UTokyoDrive/PhD/Data/WindEstTest/Comparison/'
}
windFiles <- dir(windLoc)
estFiles <- dir(estLoc, pattern = 'Z.*.csv')

for(b in 1:length(windFiles)){
    if(b == 1){
        WindDat <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
        WindDat$ID <- sub("*LatLon.txt", "", windFiles[b])
		Wind.dec <- SpatialPoints(cbind(WindDat$Lon,WindDat$Lat), proj4string = CRS("+proj=longlat"))
		UTMdat <- spTransform(Wind.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
		WindDat$UTME <- coordinates(UTMdat)[, 1]
		WindDat$UTMN <- coordinates(UTMdat)[, 2] 
    } else {
        toAdd <- read.delim(paste(windLoc, windFiles[b], sep = ''), sep = ",", header = T)
        toAdd$ID <- sub("*LatLon.txt", "", windFiles[b])
		Add.dec <- SpatialPoints(cbind(toAdd$Lon,toAdd$Lat), proj4string = CRS("+proj=longlat"))
		UTMdat <- spTransform(Add.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
		toAdd$UTME <- coordinates(UTMdat)[, 1]
		toAdd$UTMN <- coordinates(UTMdat)[, 2] 
        WindDat <- rbind(WindDat, toAdd)
    }
}
WindDat$DT <- as.POSIXct(WindDat$DT, format = "%Y-%m-%d %H:%M:%OS")

EstDat <- list(NA)
for(b in 1:length(estFiles)){
	EstDat[[b]] <- read.delim(paste(estLoc, estFiles[b], sep = ""), sep = ",", header = F)
	colnames(EstDat[[b]]) <- c("Lat","Lon","U","V")
	Est.dec <- SpatialPoints(cbind(EstDat[[b]]$Lon,EstDat[[b]]$Lat), proj4string = CRS("+proj=longlat"))
	UTMdat <- spTransform(Est.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
	EstDat[[b]]$UTME <- coordinates(UTMdat)[, 1]
	EstDat[[b]]$UTMN <- coordinates(UTMdat)[, 2] 
}


EstTimes <- as.POSIXct(paste(substr(estFiles,11,14),substr(estFiles,15,16),
		substr(estFiles,17,18), substr(estFiles,19,20), sep = "/"), format = "%Y/%m/%d/%H")

EstCompare <- WindDat[format(WindDat$DT, format = "%H") %in% c(0,3,6,9,12,15,18,21) & format(WindDat$DT, format = "%M") == "0" |
	format(WindDat$DT, format = "%H") %in% c(2,5,8,11,14,17,20,23) & format(WindDat$DT, format = "%M") == "59",]

EstComp <- unique(EstCompare[c("DT","ID")])
EstComp$DT <- round(EstComp$DT, units = "hours")
EstCompDT <- unique(EstComp$DT)

# download relevant wind datafiles
for(b in 1:length(EstCompDT)){
	download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",
		format(EstCompDT[b], format = "%Y"), "/", format(EstCompDT[b], format = "%m"), "/", format(EstCompDT[b], format = "%d"), "/",
		"Z__C_RJTD_",format(EstCompDT[b], format = "%Y"), format(EstCompDT[b], format = "%m"), format(EstCompDT[b], format = "%d"),
		format(EstCompDT[b], format = "%H"),"0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin", sep = ""), 
		destfile = paste("F:/UTokyoDrive/PhD/Data/WindEstTest/Comparison/","Z__C_RJTD_",format(EstCompDT[b], format = "%Y"),
		format(EstCompDT[b], format = "%m"), format(EstCompDT[b], format = "%d"), format(EstCompDT[b], format = "%H"),
		"0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin", sep = ""))
}

#EstCompare <- WindDat[format(WindDat$DT, format = "%M") == c(59,0),]

EstTimes = as.POSIXct(substr(estFiles, 11,22), format = "%Y%m%d%H%M")
CompDat <- data.frame(estU = NA, estV = NA, compU = NA, compV = NA)
for(b in 1:nrow(EstComp)){
	selTime <- EstComp$DT[b]
	selTag <- EstComp$ID[b]
	if(any(EstTimes == selTime)){
		num = which(EstTimes == selTime)
		# select estimated winds <25mins apart
		Sel <- WindDat[WindDat$DT >= (selTime - (50*60)) & WindDat$DT <= (selTime + (50*60)) & WindDat$ID == selTag,]
		Selsp <- SpatialPoints(cbind(Sel$Lon,Sel$Lat), proj4string = CRS("+proj=longlat"))
		Selbuff <- buffer(Selsp, width = 5000)
		est <- EstDat[[num]]
		estsp <- SpatialPoints(cbind(est$Lon,est$Lat), proj4string = CRS("+proj=longlat"))
		est <- est[which(over(estsp, Selbuff) == 1),]
		CompDat[b,] <- cbind(mean(WindDat$X[WindDat$DT > (selTime - (25*60)) & WindDat$DT < (selTime + (25*60)) & WindDat$ID == selTag]),
			mean(WindDat$Y[WindDat$DT > (selTime - (25*60)) & WindDat$DT < (selTime + (25*60)) & WindDat$ID == selTag]),
			mean(est$U), mean(est$V))
	}	
}

CompDat$estHead <- atan2(CompDat$estV,CompDat$estU)
CompDat$estSpeed <- sqrt(CompDat$estU^2 + CompDat$estV^2)
CompDat$compHead <- atan2(CompDat$compV,CompDat$compU)
CompDat$compSpeed <- sqrt(CompDat$compU^2 + CompDat$compV^2)

plot(CompDat$compSpeed,CompDat$estSpeed)
plot(CompDat$compHead,CompDat$estHead)
res<-cor.test(CompDat$estSpeed, CompDat$compSpeed, method = "pearson")
res
res2 <- cor.circular(CompDat$estHead, CompDat$compHead, test = T)


buffer(sel)

latdiff <- NA
londiff <- NA
u <- NA
v <- NA
X <- NA
Y <- NA
for(b in 1:nrow(EstCompare)){
	if(sum((abs(as.numeric(difftime(EstCompare$DT[b], EstTimes, units = "secs"))) < 60)) == 0){
		latdiff[b] <- NA
		londiff[b] <- NA
		u[b] <- NA
		v[b] <- NA
		X[b] <- NA
		Y[b] <- NA
	} else {
		sel <- EstDat[[which(abs(as.numeric(difftime(EstCompare$DT[b], EstTimes, units = "secs"))) < 60)]]
		lat <- round(EstCompare$Lat[b]/.05)*.05
		lon <- round(EstCompare$Lon[b]/.05)*.05
		latdiff[b] <- (lat - EstCompare$Lat[b])
		londiff[b] <- abs(lon - EstCompare$Lon[b])
		u[b] <- sel$U[sel$Lat == lat & sel$Lon == lon]
		v[b] <- sel$V[sel$Lat == lat & sel$Lon == lon]
		X[b] <- EstCompare$X[b]
		Y[b] <- EstCompare$Y[b]
	}
}

round(EstCompare$Lat[b]/.05)*.05
diff(sel$Lon[1:500])



install.packages('rNOMADS', repos = 'https://cloud.r-project.org/')
library(rNOMADS)
install.packages('readxl', repos = 'https://cloud.r-project.org/')
library(readxl)
install.packages('data.table', repos = 'https://cloud.r-project.org/')
library(data.table)
install.packages('dplyr', repos = 'https://cloud.r-project.org/')
library(dplyr)

grbfileloc <- '/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/Comparison'
files <- list.files(grbfileloc, pattern = '.bin', recursive = T)

GribInfo(paste(grbfileloc,files[1], sep = '/'), file.type = "grib2")

links <- LinkExtractor('http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/2018/')

GribGrab("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/2018/08/29/")

OGLoc <- '/Users/arang/OneDrive - The University of Tokyo/PhD/Data/2018Shearwater/AxyTrek/'
OGfiles <- list.files(OGLoc, pattern = '.txt', recursive = T)
OGtags <- gsub("^.*?/","",OGfiles)
OGtags <- gsub(".txt","",OGtags)
PredLoc <- "/Users/arang/OneDrive - The University of Tokyo/PhD/Data/WindEstTest/"
Preds <- list.files(PredLoc, pattern = ".csv")
tags <- gsub("WindEst","",Preds)
tags <- gsub(".csv","",tags)

for(a in 1:length(tags)){
	Dat<-read.delim(paste(PredLoc,Preds[a],sep=""), sep=",", header = F)
	colnames(Dat) <- c("DT","Head","X","Y")
	Dat$DT <- as.POSIXct(as.character(Dat$DT), format = '%Y-%m-%d %H:%M:%OS')
	OGDat <- read.delim(paste(OGLoc,OGfiles[which(tags[a] == OGtags)], sep = ''), sep = '\t', header = F)
	colnames(OGDat) <- c('DT','Lat','Lon')
	OGDat$DT <- as.POSIXct(OGDat$DT, format = '%d/%m/%Y,%H:%M:%OS')

	for(b in 1:nrow(Dat)){
		Dat$Lat[b] <- OGDat$Lat[max(which(OGDat$DT <= Dat$DT[b]))]
		Dat$Lon[b] <- OGDat$Lon[max(which(OGDat$DT <= Dat$DT[b]))]
	}

	write.table(Dat,paste(PredLoc,tags[a],'LatLon.txt',sep=''),row.names=F,col.names=T, sep = ",")
}
GribInfo("/Users/arang/Documents/gribs/Z__C_RJTD_20180829000000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin", file.type="grib2")
variables <- c("UGRD","VGRD")
levels <- ("10 m above ground")
AllGrib <- ReadGrib("/Users/arang/Documents/gribs/Z__C_RJTD_20180829000000_MSM_GPV_Rjp_Lsurf_FH34-39_grib2.bin", levels, variables)
colnames(AllGrib)

concGrib <- SubsetNOMADS(AllGrib, lon = c(min(Dat$Lon):max(Dat$Lon)), lat = c(min(Dat$Lat):max(Dat$Lat)))

# read in the data files with the latlongs
LLPreds <- list.files(PredLoc, pattern = "LatLon.txt")

for(a in 1:length(tags)){
	TempDat<-read.delim(paste(PredLoc,LLPreds[a],sep=""), sep=",", header = T)
	TempDat$Tag <- rep(tags[a],nrow(TempDat))
	if(a == 1){
		Dat <- TempDat
	} else {
		Dat <- rbind(Dat,TempDat)
	}
}

minDat <- min(unique(as.Date(Dat$DT)))
maxDat <- max(unique(as.Date(Dat$DT)))
minLat <- min(Dat$Lat)
maxLat <- max(Dat$Lat)
minLon <- min(Dat$Lon)
maxLon <- max(Dat$Lon)

write.table(Dat,paste(PredLoc,'AllLatLon.txt',sep=''),row.names=F,col.names=T)
