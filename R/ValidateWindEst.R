
if(Sys.info()['sysname'] == "Darwin"){
    windest <- read.delim("/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/AllLatLon.txt", sep = ' ', header = T)
} else {
    windest <- read.delim("F:/UTokyoDrive/PhD/Data/WindEstTest/AllLatLon.txt", sep = ' ', header = T)
}
windest$DT <- as.POSIXct(windest$DT, format = "%Y-%m-%d %H:%M:%S")
colnames(windest)
