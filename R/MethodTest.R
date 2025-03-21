if(Sys.info()['sysname'] == "Darwin"){
    testDat <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Papers/Goto2017/1700097_Real_Track_Data.csv", sep = ',', header = T)
} else {
    testDat <- read.delim("F:/UTokyoDrive/PhD/Papers/Goto2017/1700097_Real_Track_Data.csv", sep = ',', header = T)
}

cord.dec <- SpatialPoints(cbind(testDat$long, testDat$lat), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=55 +datum=WGS84"))

X <- cord.UTM$coords.x1
Y <- cord.UTM$coords.x2
ID <- as.factor(testDat$Bird_ID)
DT <- as.POSIXct(paste(testDat$Year, "/", testDat$Month, "/", testDat$Day, " ", testDat$Hour, ":", testDat$Min, ":", testDat$Sec, sep = ""),
    format = "%Y/%m/%d %H:%M:%S", tz = "")
newv <- data.frame(ID = ID, timepoint = DT, lat = testDat$lat, long = testDat$long, X = X, Y = Y)
newv$dist <- NA
newv$sp <- NA
newv$dt <- NA
for(b in unique(ID)){
    temp <- 
    newv$dist[newv$ID == b] <- c(NA, sqrt(diff(newv$X[newv$ID == b])^2 + diff(newv$Y[newv$ID == b])^2))
    newv$dt[newv$ID == b] <- c(NA, difftime(newv$timepoint[(min(which(newv$ID == b)) + 1):(max(which(newv$ID == b)))],
        newv$timepoint[min(which(newv$ID == b)):(max(which(newv$ID == b)) - 1)], units = "secs"))
}
newv$sp <- newv$dist/newv$dt

select <- testDat[testDat$Bird_ID == testDat$Bird_ID[1],]
select$DT <- as.POSIXct(paste(select$Year, "/", select$Month, "/", select$Day, " ", select$Hour, ":", select$Min, ":", select$Sec, sep = ""),
    format = "%Y/%m/%d %H:%M:%S", tz = "")
select <- select[order(select$DT),]
Selcord.dec <- SpatialPoints(cbind(select$long, select$lat), proj4string=CRS("+proj=longlat"))
Selcord.UTM <- spTransform(Selcord.dec, CRS("+proj=utm +zone=55 +datum=WGS84"))
SelX <- Selcord.UTM$coords.x1
SelY <- Selcord.UTM$coords.x2
select$dir <- c(NA, atan2(diff(SelY), diff(SelX)))
select$dist <- c(sqrt(diff(SelX)^2 + diff(SelY)^2), NA)
select$dt <- NA
select$dt[1:(nrow(select) - 1)] <- difftime(select$DT[2:nrow(select)], select$DT[1:(nrow(select) - 1)], units = "secs")
select$spd <- select$dist/select$dt

for(b in 1:length(unique(testDat$Bird_ID))){
    select <- testDat[testDat$Bird_ID == testDat$Bird_ID[1],]
    select$DT <- as.POSIXct(paste(select$Year, "/", select$Month, "/", select$Day, " ", select$Hour, ":", select$Min, ":", select$Sec, sep = ""),
        format = "%Y/%m/%d %H:%M:%S", tz = "")
    select <- select[order(select$DT),]
    Selcord.dec <- SpatialPoints(cbind(select$long, select$lat), proj4string=CRS("+proj=longlat"))
    Selcord.UTM <- spTransform(Selcord.dec, CRS("+proj=utm +zone=55 +datum=WGS84"))
    SelX <- Selcord.UTM$coords.x1
    SelY <- Selcord.UTM$coords.x2
    select$track_direction <- c(NA, atan2(diff(SelY), diff(SelX)))
    select$dist <- c(sqrt(diff(SelX)^2 + diff(SelY)^2), NA)
    select$dt <- NA
    select$dt[1:(nrow(select) - 1)] <- difftime(select$DT[2:nrow(select)], select$DT[1:(nrow(select) - 1)], units = "secs")
    select$track_speed <- select$dist/select$dt
    if(b == 1){
        allD <- select
    } else {
        allD <- rbind(allD, select)
    }
}
allD <- allD[!is.na(allD$track_direction),]
allD <- allD[!is.na(allD$track_speed),]

sum(is.na(allD$track_direction))
sum(is.na(allD$track_speed))

colnames(testDat)

plot(select$spd, select$track_speed)

plot(atan2(diff(SelY[1:10]), diff(SelX[1:10]))+2*pi, select$track_direction[2:10])

cbind(atan2(diff(SelY[1:10]), diff(SelX[1:10])), select$track_direction[2:10])

median(atan2(diff(SelY), diff(SelX)) + pi)
median(select$track_direction)

plot(SelX,SelY)

select

colnames(select)
median(atan2(diff()))




plot(testDat$track_speed[testDat$Bird_ID == ID[1]], newv$sp[newv$ID == ID[1]])

median(testDat$track_direction[testDat$Bird_ID == ID[1]])*(180/pi)
(median(atan2(diff(Y[newv$ID == ID[1]]),diff(X[newv$ID == ID[1]]))))*(180/pi)

x <- X[newv$ID == ID[1]]
y <- Y[newv$ID == ID[1]]

median(atan2(x[1:(length(x) - 1)] - x[2:length(x)], y[1:(length(y) - 1)] - y[2:length(y)]))*(180/pi)

median(atan2(diff(x), diff(y)))*(180/pi)
median(testDat$track_direction[testDat$Bird_ID == ID[1]])*(180/pi)




plot(newv$sp,ylim=c(0,100))

plot(testDat$long,testDat$lat)

plot(X,Y)

unlist(sapply(ID, function(x)
    c(NA, sqrt(diff(newv$X[newv$ID == x]))^2 + sqrt(diff(newv$Y[newv$ID == x]))^2)))


distances <- sqrt(diff(X)^2 + diff(Y)^2)
track_direction <- atan2(diff(Y),diff(X))
dt <- difftime()
track_speed <- (distances)/testDat$dt[2:(nrow(testDat))]
colnames(testDat)
testDat[1:5,3]
plot(track_speed*3.6, testDat$track_speed[2:(nrow(testDat))])
plot(testDat$track_speed)
plot(track_speed)
testDat[1,]
