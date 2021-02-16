
# code to format AxyTrek data for use in Goto2017 wind estimation method

install.packages('circular', repos = 'https://cloud.r-project.org/')
library(circular)
install.packages('readxl', repos = 'https://cloud.r-project.org/')
library(readxl)
install.packages('data.table', repos = 'https://cloud.r-project.org/')
library(data.table)
install.packages('dplyr', repos = 'https://cloud.r-project.org/')
library(dplyr)
install.packages('adehabitatLT', repos = 'https://cloud.r-project.org/')
library(adehabitatLT)
install.packages('chron', repos = 'https://cloud.r-project.org/')
library(chron)

# read in
datin <- paste("/Users/arang/OneDrive - The University of Tokyo/PhD/Data/2018Shearwater/AxyTrek/2018-01/1_S2.txt")
outfile<-paste("~/Desktop/WindEstTest1.csv",sep=",")

AxDat <- read.delim(datin, sep = '\t', header = F)
timepoint <- as.POSIXct(as.character(AxDat[,1]), format = '%d/%m/%Y,%H:%M:%OS')
dt <- as.numeric(diff(timepoint))
#subsample to 1 Hz data
sampling_interval <- median(dt)
#sampled per minute
s_per_min <- 60/sampling_interval
#index series for extracting subsampled data
to_sub <- seq(from=1,to=nrow(AxDat),by=s_per_min)
#subsample
AxDat <- AxDat[to_sub,]
sampling_interval <- median(dt)

timepoint <- as.POSIXct(as.character(AxDat[,1]), format = '%d/%m/%Y,%H:%M:%OS')
dt <- as.numeric(diff(timepoint))
lat <- AxDat[,2]
long <- AxDat[,3]

cord.dec <- SpatialPoints(cbind(long, lat), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))

X <- coordinates(cord.UTM)[,1]
Y <- coordinates(cord.UTM)[,2]

distances <- sqrt(diff(X)^2 + diff(Y)^2)
track_direction <- atan2(diff(Y),diff(X))
track_speed <- distances/dt
kph <- distances*10^-3/dt*3600

# # remove points where speed is unrealistic
# while(max(kph) > 80){
# 	toRm <- which(kph > 80) + 1
# 	timepoint <- timepoint[-toRm]
# 	X <- X[-toRm]
# 	Y <- Y[-toRm]
# 	dt <- as.numeric(diff(timepoint))
# 	distances <- sqrt(diff(X)^2 + diff(Y)^2)
# 	track_direction <- atan2(diff(Y),diff(X))
# 	track_speed <- distances/dt
# 	kph <- distances*10^-3/dt*3600
# }

sampling_interval <- 60
time_window <- 51
cutlength <- 45
cutv <- 4.1667
constv <- 34.7/3.6

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
		})}		



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


########################################################

	title<-list("time_point_center","heading_est","flow_x_est","flow_y_est")
	write.table(title,outfile,quote=F,col.name=F,row.name=F,append=F,sep=", ")


	error_of_sampling_interval <- 5
	cutt<- sampling_interval + error_of_sampling_interval  
# upper value of sampling interval [sec]
	windwidth <- time_window - 1                             
# length of time window(velocity)  [min]

	rrow <- track_speed
	drow <- track_direction
	# set time (add 9 hours to convert from GMT)
	tp <- timepoint+(9*3600)

TestSel <- which(tp > as.POSIXct('02/09/2018 14:00:00', format = '%d/%m/%Y %H:%M:%OS') &
	tp < as.POSIXct('02/09/2018 16:00:00', format = '%d/%m/%Y %H:%M:%OS'))
rrow <- rrow[TestSel]
drow <- drow[TestSel]
tp <- tp[TestSel]
dt <- as.numeric(diff(tp))
X <- X[TestSel]
Y <- Y[TestSel]

	startpoint<-floor(windwidth*(60/sampling_interval)/2)
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
			} else {print(paste(k,'Failed step 1',sep=' '))}
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
			} else {print('Failed step 4')}

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
						} else {print('Failed step 8')}
						
						#STEP 9
						if(((answ$value)>max_like)&&(condition_all==1)){
							max_like <- answ$value
							ans_best <- answ
						} else {print('Failed step 9')}
					} else {print('Failed step 7')}
				} else {print('Failed step 6')}
			} else {print('Failed step 5')}
		}

# output estimation result
		#STEP 10
		if(max_like!="NaN"){
			outwa<-list(tp[center],atan2(ans_best$par[3],ans_best$par[2]),ans_best$par[4],ans_best$par[5])
			write.table(outwa,outfile,quote=F,col.name=F,row.name=F,append=T,sep=", ")
		}  else {print('Failed step 10')}
	}  else {print('Failed step 2')}
}
