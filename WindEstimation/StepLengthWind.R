library(rerddap)
library(splitr)
library(dplyr)
library(plyr)
library(ggplot2)
library(rnaturalearth)
library(sp)
library(circular)
library(readxl)
library(data.table)
library(dplyr)
library(adehabitatLT)
library(chron)
if(Sys.info()['sysname'] == "Darwin"){
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("/Volumes/GoogleDrive/My Drive/PhD/Data/Temp2018.RData")
} else {
    load("F:/UTokyoDrive/PhD/Data/2019Shearwater/2019Dat.RData")
    load("F:/UTokyoDrive/PhD/Data/Temp2018.RData")
}
###############################################
# Zero crossings locations ####################
###############################################
StLDisp <- function(utm){
  chgs <- diff(utm) >= 0
  swtch <- which(diff(chgs)!=0) + 1
}
###############################################
# Step length calculator ######################
###############################################
StLCalc <- function(DT,UTMN,UTME){
  deltaT <- c(NA, difftime(DT[2:length(DT)], DT[1:(length(DT)-1)], units = "secs"))
  # find the cutoff points from the time difference
  cutoff <- median(deltaT, na.rm = T) + 10
  # find the start and end times based on time differences over the cutoff value
  strts <- which(deltaT > cutoff)+1
  strts <- c(1,strts[!deltaT[strts] > cutoff])
  ends <- c(which(deltaT > cutoff)-1,length(DT))
  ends <- ends[!deltaT[ends] > cutoff]
  # find the changes in sign of the data
  stepLsN <- NA
  for(g in 1:length(strts)){
    iso <- UTMN[strts[g]:ends[g]]
    stepLsN[g] <- sum(abs(diff(iso)))
  }
  stepLsE <- NA
  for(g in 1:length(strts)){
    iso <- UTME[strts[g]:ends[g]]
    stepLsE[g] <- sum(abs(diff(iso)))
  }
  data.frame(start=DT[strts], end=DT[ends], dispN = stepLsN, dispE = stepLsE)
}
###############################################
# Levy distribution estimation ################
############################################### 
levy <- function(x,mu){
  x^-mu
}
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
sel <- as.data.frame(allD[allD$tagID == allD$tagID[1] & allD$Year == allD$Year[1],])
sel$deltaT <- c(NA, difftime(sel$DT[2:nrow(sel)], sel$DT[1:(nrow(sel)-1)], units = "secs"))
sel$U <- NA
sel$V <- NA
selTest <- StLCalc(sel$DT, sel$UTMN, sel$UTME)
difftime(selTest$end,selTest$start, units= 'mins')
median(sel$deltaT,na.rm=T)

for(h in 1:nrow(selTest)){
    isolate <- sel[which(sel$DT == selTest$start[h]):which(sel$DT == selTest$end[h]),]
    lat <- isolate$lat
    long <- isolate$long
    dt <- isolate$deltaT[2:nrow(isolate)]
    vg_x_obs <- isolate$UTMN[2:nrow(isolate)] - isolate$UTMN[1:(nrow(isolate) - 1)]
	vg_y_obs <- isolate$UTME[2:nrow(isolate)] - isolate$UTME[1:(nrow(isolate) - 1)]

	g_speed <- sqrt(vg_y_obs^2 + vg_x_obs^2)/dt
	g_direction <- atan2(vg_y_obs, vg_x_obs)
    sampling_interval <- median(dt)
    time_window <- 51     #  time length of time window [min] *Give an odd number!
	cutlength <- ceiling((45/51)*(time_window*(60/sampling_interval)))     # minimum number of data 
    cutv <- 4.1667
	constv <- 34.7/3.6
    error_of_sampling_interval <- 5
	cutt<- sampling_interval + error_of_sampling_interval
    windwidth <- time_window - 1                             
# length of time window(velocity)  [min]
	# TestSel <- ToUse
	rrow <- g_speed
	drow <- g_direction
	tp <- isolate$DT
	# dt <- as.numeric(diff(tp))
	X <- isolate$UTMN
	Y <- isolate$UTME
    startpoint<-floor((windwidth*(60/sampling_interval))/2)
	endpoint<-length(rrow)-startpoint
    if(endpoint < 0 | is.na(endpoint)){
        next
    }
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
				break
				outwa<-list(tp[center],atan2(ans_best$par[3],ans_best$par[2]),ans_best$par[4],ans_best$par[5])
                sel$U[sel$DT == tp[center] & sel$lat == lat[center] & sel$lon == long[center]] = ans_best$par[4]
                sel$V[sel$DT == tp[center] & sel$lat == lat[center] & sel$lon == long[center]] = ans_best$par[5]
				# write.table(outwa,paste(outfile,tags[a],'.csv', sep = ''),quote=F,col.name=F,row.name=F,append=T,sep=", ")
			}#  else {print('Failed step 10')}
	}#  else {print('Failed step 2')}
}
}



