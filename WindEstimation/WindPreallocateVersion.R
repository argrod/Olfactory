
# code to format AxyTrek data for use in Goto2017 wind estimation method

install.packages('circular', repos = 'https://cloud.r-project.org/')
install.packages('readxl', repos = 'https://cloud.r-project.org/')
install.packages('data.table', repos = 'https://cloud.r-project.org/')
install.packages('dplyr', repos = 'https://cloud.r-project.org/')
install.packages('adehabitatLT', repos = 'https://cloud.r-project.org/')
install.packages('chron', repos = 'https://cloud.r-project.org/')

library(circular)
library(readxl)
library(data.table)
library(dplyr)
library(adehabitatLT)
library(chron)

### Functions

###############################################
# Log-likelihood of the model #################
###############################################
Likelihoodww<-function(data1,data2,cv){
	return(function(par){
		a  <-par[1]
		b  <-cv/gamma(1+1/a)
		mx <-par[2]
		my <-par[3]
		wx <-par[4]
		wy <-par[5]
		L    <- 0
		for(i in 1:length(data1)){ 
			rr<-sqrt((data1[i]*cos(data2[i])-wx)^2+(data1[i]*sin(data2[i])-wy)^2)
			rx<- (data1[i]*cos(data2[i])-wx)/rr
			ry<- (data1[i]*sin(data2[i])-wy)/rr
			lp<- (a-2)*log(rr)-(rr/b)^a+mx*rx+my*ry+log(a)-log(b)+(1-a)*log(b)-log(besselI(sqrt(mx^2+my^2),0,))
			L<- L+lp					} 
			return(L)
		}
		)
}		



###############################################
# Standard deviation of Weibull distribution ##     
###############################################                                 
Weibull_sd<-function(a,b){
	w_sd<-b*sqrt(gamma(1+2/a)-gamma(1+1/a)*gamma(1+1/a)) 
	return(w_sd)
}



###############################################
# Mean of Weibull distribution ################
###############################################
Weibull_mean<-function(a,b){
	w_m<-b*gamma(1+1/a) 
	return(w_m)
}



###############################################
# Standard deviation of von-Mises distribution (approximated #to Gaussian) 
###############################################
Von_Mises_sd<-function(kappa){
	vm_sd<-1/sqrt(kappa)
	return(vm_sd)
}

###############################################
# recalculate speeds using moving window
###############################################
movSpd <- function(DTGPS, UTMN, UTME, f){
    if(length(UTMN) != length(UTME)){
        stop("Lengths of UTM values not equal")
    }
    spTrav <- NA
    dst <- NA
    dt <- NA
    for(b in 1:length(UTMN)){
        if(b <= (f/2 + 1)){
            dst[b] <- sqrt((UTME[ceiling(f/2)] - UTME[1])^2 + (UTMN[ceiling(f/2)] - UTMN[1])^2)
            dt[b] <- as.numeric(difftime(DTGPS[ceiling(f/2)],DTGPS[1], units = 'secs'))
            spTrav[b] <- dst[b]/dt[b]
        } else if(b >= (length(UTMN) - ceiling(f/2 + 1))){
            dst[b] <- sqrt((UTME[length(UTME)] - UTME[length(UTME) - ceiling(f/2)])^2 + (UTMN[length(UTMN)] - UTMN[length(UTMN) - ceiling(f/2)])^2)
            dt[b] <- as.numeric(difftime(DTGPS[length(DTGPS)],DTGPS[length(DTGPS) - ceiling(f/2)], units = 'secs'))
            spTrav[b] <- dst[b]/dt[b]
        } else {
            dst[b] <- sqrt((UTME[b + ceiling(f/2)] - UTME[b - ceiling(f/2)])^2 + (UTMN[b + ceiling(f/2)] - UTMN[b - ceiling(f/2)])^2)
            dt[b] <- as.numeric(difftime(DTGPS[b + ceiling(f/2)], DTGPS[b - ceiling(f/2)], units = 'secs'))
            spTrav[b] <- dst[b]/dt[b]
        }
    }
    movSpeed <- list("distance" = dst, "timediff" = dt, "speed" = spTrav)
    return(movSpeed)
}

########################################################

# read in
if(Sys.info()['sysname'] == "Darwin"){
	outfile<-paste("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/AxyFS/",sep=",")
	fileloc <- "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/AxyTrek/"
} else {
	outfile<-paste("F:/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/AxyFS/",sep=",")
	fileloc <- "F:/UTokyoDrive/PhD/Data/2018Shearwater/AxyTrek/"
}
files <- list.files(fileloc, pattern = '.txt', recursive = T)

tags <- unique(sub('.*/','',files))
tags <- unique(sub('*.txt','',tags))

s_intervals <- c(5,30,30,30,30,30,30,30,30,30,5)

for(a in 1:length(files)){
	AxDat <- read.delim(paste(fileloc,files[a], sep = ''), sep = '\t', header = F)
	timepoint <- as.POSIXct(as.character(AxDat[,1]), format = '%d/%m/%Y,%H:%M:%OS') + (9*3600)
	tsel <- seq(timepoint[1], timepoint[length(timepoint)], by = 60)
	dt <- as.numeric(difftime(timepoint[2:length(timepoint)], timepoint[1:(length(timepoint) - 1)], units = "secs"))
	lat <- AxDat[,2]
	long <- AxDat[,3]

	cord.dec <- SpatialPoints(cbind(long, lat), proj4string=CRS("+proj=longlat"))
	cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

	X <- cord.UTM$long
	Y <- cord.UTM$lat

	FkOshi <- data.frame('Lat'=39.402289,'Long'=141.998165)
	FkOshi.dec <- SpatialPoints(cbind(FkOshi$Long,FkOshi$Lat),proj4string=CRS('+proj=longlat'))
	FkOshi.UTM <- spTransform(FkOshi.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
	FkUTME <- coordinates(FkOshi.UTM)[,1]
	FkUTMN <- coordinates(FkOshi.UTM)[,2]

	DistFromFk <- (sqrt((X-FkUTME)^2 + (Y-FkUTMN)^2))*10^-3

	close <- which(DistFromFk > 5)

	distances <- sqrt(diff(X)^2 + diff(Y)^2)
	track_direction <- atan2(diff(Y),diff(X))
	track_speed <- distances/dt
	
	temp <- movSpd(timepoint, X, Y, 5)
    recalDist <- temp[[1]]
    recaldt <- temp[[2]]
    recalSp <- temp[[3]]
	
	sampling_interval <-  s_intervals[a]   #  sampling interval [sec]
	time_window <- 51     #  time length of time window [min] *Give an odd number!
	cutlength <- ceiling((45/51)*(time_window*(60/sampling_interval)))     # minimum number of data points (track vectors) included in a time window [points]         
	cutv <- 4.1667 # minimum ground speed [m/sec]

	###         Condition 3: give mean air speed value      #####
	constv <- 34.7/3.6 # mean air speed [m/sec]
	#We gave the mean ground speed of streaked shearwater (Shiomi #et. al. 2012) as the mean air speed

	
	# sampling_interval <- median(dt)
	# time_window <- 51
	# cutlength <- round(.8*(time_window*(60/sampling_interval)))
	# cutv <- 4.1667
	# constv <- 34.7/3.6
	# error_of_sampling_interval <- 5
	cutt<- sampling_interval + error_of_sampling_interval  
# upper value of sampling interval [sec]
	windwidth <- time_window - 1                             
# length of time window(velocity)  [min]

	# FL <- kph>15
	# suitable_index <- data.frame(Start=double(),End=double(),Breaks=double())

	# for(b in 1:(length(FL)-(time_window*sampling_interval))){
	# 	if(any(which(difftime(timepoint[b:length(timepoint)],timepoint[b],units='mins') > 51))){
	# 		EndInd <- min(which(difftime(timepoint[b:length(timepoint)],timepoint[b],units='mins') > 51)) + b
	# 		if(sum(FL[b:EndInd]) > cutlength){
	# 			suitable_index[b,1] <- b
	# 			suitable_index[b,2] <- EndInd
	# 		}
	# 	}
	# }
	# suitable_index <- suitable_index[-which(is.na(suitable_index$End)),]
	
	# unik <- !duplicated(suitable_index$End)
	# USuitInd <- suitable_index[unik,]

	# # FitSeg <-  data.frame(Start=double(),End=double(),Dur=double(),StT=character(),stringsAsFactors=F)
	# # if(length(end) == length(st)){
	# # 	for(b in 1:length(st)){
	# # 		if(end[b]-st[b] > cutlength && ((as.numeric(timepoint[end[b]])-as.numeric(timepoint[st[b]]))/60) > time_window){
	# # 			FitSeg[b,1] <- st[b]
	# # 			FitSeg[b,2] <- end[b]
	# # 			FitSeg[b,3] <- (as.numeric(timepoint[end[b]])-as.numeric(timepoint[st[b]]))/60
	# # 			FitSeg[b,4] <- as.character(timepoint[st[b]])
	# # 		}
	# # 	}
	# # 	FitSeg <- na.omit(FitSeg)
	# # } else {
	# # 	for(b in 1:length(end)){
	# # 		if(b == length(end)){
	# # 			if(nrow(AxDat)-st[b+1] > cutlength && ((as.numeric(timepoint[nrow(AxDat)])-as.numeric(timepoint[st[b+1]]))/60) > 51){
	# # 				FitSeg[b,1] <- st[b+1]
	# # 				FitSeg[b,2] <- nrow(AxDat)
	# # 				FitSeg[b,3] <- (as.numeric(timepoint[nrow(AxDat)])-as.numeric(timepoint[st[b+1]]))/60
	# # 				FitSeg[b,4] <- as.character(timepoint[st[b+1]])
	# # 			}
	# # 		} 
	# # 		if(end[b]-st[b] > cutlength && ((as.numeric(timepoint[end[b]])-as.numeric(timepoint[st[b]]))/60) > time_window){
	# # 			FitSeg[b,1] <- st[b]
	# # 			FitSeg[b,2] <- end[b]
	# # 			FitSeg[b,3] <- (as.numeric(timepoint[nrow(AxDat)])-as.numeric(timepoint[st[b]]))/60
	# # 			FitSeg[b,4] <- as.character(timepoint[st[b]])
	# # 		}
	# # 	}
	# # }
	# # FitSeg <- na.omit(FitSeg)
	# # if(nrow(FitSeg) == 0){
	# # 	next
	# # }
	# for(b in 1:nrow(USuitInd)){
	# 	if(b == 1){
	# 		ToUse <- USuitInd$Start[b]:USuitInd$End[b]
	# 	} else {
	# 		ToUse <- c(ToUse,USuitInd$Start[b]:USuitInd$End[b])
	# 	}
	# }
	
	# TestSel <- ToUse
	rrow <- track_speed
	drow <- track_direction
	tp <- timepoint
	dt <- as.numeric(diff(tp))
	X <- X
	Y <- Y

	startpoint<-floor((windwidth*(60/sampling_interval))/2)
	endpoint<-length(rrow)-startpoint


	for(center in startpoint:endpoint){
########################################################


# select data points in time window

		windwidthsec<-(windwidth/2)*60+15 #half length of time window(velocity) [sec]

		inter<-0　# time from center to end of the window
		for(qf in 1:(length(dt)-center)){
			inter<-inter+dt[center+qf]
			if(inter>windwidthsec)
				break
		}

		inter<-0　# time from start to center of the window
		for(qb in 0:center-1){
			inter<-inter+dt[center-qb]
			if(inter>windwidthsec)
				break
		}

		end <- center+qf-1 # a data point at the end of a time window
		st <- center-qb   # a data point at the start of a time window
########################################################

		answ<-NULL
		r<-as.list(NULL)
		d<-as.list(NULL)
		index<-as.list(NULL)

#STEP 1
		for(k in st:end){	
			if((rrow[k]>cutv)&&(dt[k]<cutt)&&(drow[k]!=100)){
				r    <-append(r,rrow[k])
				d    <-append(d,drow[k])
				index<-append(index,k)
			} #else {print(paste(k,'Failed step 1',sep=' '))}
		}
		r<-unlist(r)
		d<-unlist(d)
		index<-unlist(index)
# r: speed of track vectors in the time window 
# d: direction of track vectors in the time window

#STEP 2
		if(length(r)>=cutlength){
			max_like <- "NaN"

			hdtry=3;
			for(id_hd in -hdtry:hdtry){ #Change heading direction 
			rr<-as.list(NULL)
			dd<-as.list(NULL)
			iindex<-as.list(NULL)

			for(k in 1:length(r)){
		#STEP 3
				if(r[k]>cutv){
					rr    <-append(rr,r[k])
					dd    <-append(dd,d[k])
					iindex<-append(iindex,index[k])
				}  else {print('Failed step 3')}
			}

			r<-unlist(rr)
			d<-unlist(dd)
			index<-unlist(iindex)

			seg_length<-length(r)
#if(seg_length<cutlength)
#break;

			inithd_first<-id_hd/hdtry*pi/2 #initial value of heading
			inita <- 0 #initial angle
			while(inita < 5 ){
				inita<-abs(rnorm(1,12.5,5))
			}

			meangd<-atan2(mean(sin(d)),mean(cos(d))) #mean track heading
			inithd<-meangd+inithd_first 
			initkappa<- A1inv(mean(cos(d-meangd))) #estimate concentration parameter for Von Mises distribution
			initmux<-initkappa*cos(inithd)
			initmuy<-initkappa*sin(inithd)
			initwx <-mean(r)*cos(meangd)-constv*cos(inithd)
			initwy <- mean(r)*sin(meangd)-constv*sin(inithd)

# MLE
			answ<-optim(par=c(inita,initmux,initmuy,initwx,initwy),fn=Likelihoodww(r,d,constv),control=list(fnscale=-1))

			yoko<-Von_Mises_sd(sqrt(answ$par[3]*answ$par[3]+answ$par[2]*answ$par[2]))*Weibull_mean(answ$par[1],constv/gamma(1+1/answ$par[1]))
#    yoko: SD of the heading vector distribution to the
#perpendicular direction relative to the mean direction	
#     
			tate<-Weibull_sd(answ$par[1],constv/gamma(1+1/answ$par[1]))
#    tate: SD of the heading vector distribution along the mean direction

#Repeat MLE to ensure the convergence of optimisation
	#STEP 4
			if(tate!="NaN"){
#if(yoko>tate){
				for(try_no in 1:100){
					inita  <-answ$par[1]
					initmux<-answ$par[2]
					initmuy<-answ$par[3]
					initwx <-answ$par[4]
					initwy <-answ$par[5]
					answ<-optim(par=c(inita,initmux,initmuy,initwx,initwy),fn=Likelihoodww(r,d,constv),control=list(fnscale=-1))
					if(answ$convergence==0)
						break;
					yoko<-Von_Mises_sd(sqrt(answ$par[3]*answ$par[3]+answ$par[2]*answ$par[2]))*Weibull_mean(answ$par[1],constv/gamma(1+1/answ$par[1]))#/2
					tate<-Weibull_sd(answ$par[1],constv/gamma(1+1/answ$par[1]))
					if(tate=="NaN")
						break;
				}
			}# else {print('Failed step 4')}

			yoko<-Von_Mises_sd(sqrt(answ$par[3]*answ$par[3]+answ$par[2]*answ$par[2]))*Weibull_mean(answ$par[1],constv/gamma(1+1/answ$par[1]))#/2
			tate<-Weibull_sd(answ$par[1],constv/gamma(1+1/answ$par[1]))

	#STEP 5
			if(tate!="NaN"){
		#STEP 6
				if(answ$convergence==0){

#calcurate speed (nr) and direction(nd) of heading
			#STEP 7
					if(seg_length>=cutlength){
						nr<-c(1:length(r))
						nd<-c(1:length(d))
						wx<-answ$par[4]
						wy<-answ$par[5]
						for(i in 1:length(r)){
							nvx<-r[i]*cos(d[i])-wx
							nvy<-r[i]*sin(d[i])-wy
							nr[i]<-sqrt(nvx*nvx+nvy*nvy)
							nd[i]<-atan2(nvy,nvx)
						}

#############################
### goodness of fit tests ###
#############################
#speed of heading vecotr
						nrp<-ks.test(nr,"pweibull",answ$par[1],constv/gamma(1+1/answ$par[1]))

#direction of heading vecotr
						mu<-atan2(answ$par[3],answ$par[2])
						kappa<-sqrt(answ$par[2]*answ$par[2]+answ$par[3]*answ$par[3])
						ndp<-ks.test(nd,"pvonmises",mu,kappa)

#correlation test between direction and speed of heaeding
#vector
						cnrnd<-cor.test(nr,nd,method="p")
#############################



						condition1    <- if((yoko/tate)>1){1}else{0}
#condition1: 1(condition1 satisfied) 0(not satisfied)

						condition2    <- if((cos(meangd)*cos(mu)+sin(meangd)*sin(mu))>0){1}else{0}
#condition2: 1(condition2 satisfied) 0(not satisfied)

						conditionfit  <- if(((nrp$p.value)>0.05)&&((ndp$p.value)>0.05)&&((cnrnd$p.value)>0.05)){1}else{0}
# conditionfit: 1(passed goodnes of fittest) 0(not passed)

						condition_all <- condition1*condition2*conditionfit
#condition_all: 1(all condition 1,2 and goodness of fit are satisfied) 0(not satisfied)

				#STEP 8
						if((max_like=="NaN")&&(condition_all==1)){
							max_like <- answ$value
							ans_best <- answ
						} #else {print('Failed step 8')}

				#STEP 9
						if(((answ$value)>max_like)&&(condition_all==1)){
							max_like <- answ$value
							ans_best <- answ
						} #else {print('Failed step 9')}
					} else {print('Failed step 7')}
				}# else {print('Failed step 6')}
			
			}# else {print('Failed step 5')}
			}

		# output estimation result
			#STEP 10
			if(max_like!="NaN"){
				outwa<-list(tp[center],atan2(ans_best$par[3],ans_best$par[2]),ans_best$par[4],ans_best$par[5])
				write.table(outwa,paste(outfile,tags[a],'.csv', sep = ''),quote=F,col.name=F,row.name=F,append=T,sep=", ")
			}#  else {print('Failed step 10')}
	}#  else {print('Failed step 2')}
}
}
