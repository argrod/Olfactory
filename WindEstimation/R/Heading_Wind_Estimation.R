#Heading_Wind_Estimation.R

#Please install 'circular' package
install.packages("circular", dependencies = TRUE, repos = 'https://cloud.r-project.org/')
library(circular)
#Please install 'circular' package

#*******************************************************
# This porgram estimates heading direction and flow vecotr #from tracking data 
# When you input the data "Simulation_trck.csv" generated #from Sample_Track_Simulation.R to this program,
# this program also shows the histogram of the difference #between the value of the estimated and true value of heading #direction, x component of flow and y component of flow after #estimation. 
#
#
#*** input data format   ***
# input data should contain following headers and #corresponding data sets (sample data #"Simulation_Track.csv" can be generated using #"Sample_Track.R")
#   $timepoint       : time point of data
#   $dt              : elapsed time from  the fix recorded at previous time point [sec]
#   $track_speed     : speed of track vector [m/sec]
#   $track_direction : direction of track vector [rad]  
#                      (East:0 radian,  North:pi/2 radian, #West:pi radian, South:pi*3/2 radian, When direction can't #be calcurated for track speed=0[m/sec]: 100)
#
#*** output data format   ***
# this program outputs following estimations resluts
# $time_point_center : time point at the center of time window
# $heading_est       : estimated mean heading direction
# $flow_x_est        : estimated x component of flow
# $flow_y_est        : estimated y component of flow
#  
#*******************************************************


########################################################
### Please input following parameter values.  ########################################################

### Specify the place and name of track data 
infile <- paste(("~/Downloads/1700097_real_track_data.csv")
# sample data named "Simulation_Track.csv" can be generated #using "Sample_Track_Simulation.R"

### Specify the place and name of output data 
outfile<-paste("C:/Users/arang/Desktop/AmendedV.csv",sep=", ")


sampling_interval <- 60     #  sampling interval [sec]
time_window <- 51     #  time length of time window [min] *Give an odd number!
cutlength <- 45     # minimum number of data points (track vectors) included in a time window [points]         
cutv <- 4.1667 # minimum ground speed [m/sec]

###         Condition 3: give mean air speed value      #####
constv <- 34.7/3.6 # mean air speed [m/sec]
#We gave the mean ground speed of streaked shearwater (Shiomi #et. al. 2012) as the mean air speed


########################################################
### After checking (1) format of the input data and (2) values of above parameter, please conduct this program.    
########################################################











### Funcitons   

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

	error_of_sampling_interval <-5
	cutt<- sampling_interval+error_of_sampling_interval  
# upper value of sampling interval [sec]
	windwidth<-time_window-1                             
# length of time window(velocity)  [min]

	raw <- read.csv(infile)
	tp <- raw$timepoint        # time point of data
	dt <- raw$dt               # elapsed time from  the fix recorded at previous time point [sec]
	rrow <- raw$track_speed      # speed of track vector [m/sec]
	drow <- raw$track_direction  # direction of track vector [rad]  
 # (East:0 radian,  North:pi/2 radian, West:pi radian, #South:pi*3/2 radian,When direction can't be calcurated for #track speed =0[m/sec]: 100)

	startpoint<-floor(windwidth/2)
	endpoint<-length(rrow)-floor(windwidth/2)

	for(center in startpoint:endpoint){
########################################################


# select data points in time window

		windwidthsec<-(windwidth/2)*60+15 #half length of time window(velocity) [sec]

		inter<-0　# time from center to end of the window
		for(qf in 1:(length(dt)-center)){
			inter<-inter+dt[center+qf]
			if(inter>windwidthsec)break}

			inter<-0　# time from start to center of the window
			for(qb in 0:center-1){
				inter<-inter+dt[center-qb]
				if(inter>windwidthsec)break}

				end <- center+qf-1 # a data point at the end of a time window
				st <- center-qb   # a data point at the start of a time window
########################################################

				answ<-NULL
				r<-as.list(NULL)
				d<-as.list(NULL)
				index<-as.list(NULL)

				for(k in st:end){	
					if((rrow[k]>cutv)&&(dt[k]<cutt)&&(drow[k]!=100)){
						r    <-append(r,rrow[k])
						d    <-append(d,drow[k])
						index<-append(index,k)
					}
				}
				r<-unlist(r)
				d<-unlist(d)
				index<-unlist(index)
# r: speed of track vectors in the time window 
# d: direction of track vectors in the time window

				if(length(r)>=cutlength){
					max_like <- "NaN"

					hdtry=3;
					for(id_hd in -hdtry:hdtry){ #Change heading direction 
					rr<-as.list(NULL)
					dd<-as.list(NULL)
					iindex<-as.list(NULL)

					for(k in 1:length(r)){
						if(r[k]>cutv){
							rr    <-append(rr,r[k])
							dd    <-append(dd,d[k])
							iindex<-append(iindex,index[k])
						}
					}

					r<-unlist(rr)
					d<-unlist(dd)
					index<-unlist(iindex)

					seg_length<-length(r)
#if(seg_length<cutlength)
#break;

					inithd_first<-id_hd/hdtry*pi/2 #initial value of heading
					inita <- 0
					while(inita < 5 ){
						inita<-abs(rnorm(1,12.5,5))
					}

					meangd<-atan2(mean(sin(d)),mean(cos(d)))
					inithd<-meangd+inithd_first 
					initkappa<- A1inv(mean(cos(d-meangd)))
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

#Repeant MLE to ensure the convergence of optimisation
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
					}

					yoko<-Von_Mises_sd(sqrt(answ$par[3]*answ$par[3]+answ$par[2]*answ$par[2]))*Weibull_mean(answ$par[1],constv/gamma(1+1/answ$par[1]))#/2
					tate<-Weibull_sd(answ$par[1],constv/gamma(1+1/answ$par[1]))

					if(tate!="NaN"){
						if(answ$convergence==0){

#calcurate speed (nr) and direction(nd) of heading
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


								if((max_like=="NaN")&&(condition_all==1)){
									max_like <- answ$value
									ans_best <- answ
								}

								if(((answ$value)>max_like)&&(condition_all==1)){
									max_like <- answ$value
									ans_best <- answ
								}
							}
						}
					}
				}

# output estimation result
				if(max_like!="NaN"){
					outwa<-list(tp[center],atan2(ans_best$par[3],ans_best$par[2]),ans_best$par[4],ans_best$par[5])
					write.table(outwa,outfile,quote=F,col.name=F,row.name=F,append=T,sep=", ")
				}
			}
		}

		est<-read.csv(outfile)
		true<-read.csv(infile)
		dev.new()
		par(mfrow=c(3,1))
		if((length(true$mean_flow_x)!=0)&&(length(true$mean_flow_y)!=0)&&(length(true$mean_head_direction)!=0)){
			flow_x_est <- est$flow_x_est
			flow_y_est <- est$flow_y_est
			heading_est <- 180/pi*(est$heading)
			flow_x_true <- mean(true$mean_flow_x)
			flow_y_true <- mean(true$mean_flow_y)
			heading_true <- 180/pi*mean(true$mean_head_direction)

			hist(heading_est-heading_true, main="Difference between estimated and true headings [in degrees]")
			hist(flow_x_est-flow_x_true, main="Difference between estimated and true flows (x components) [m/sec]")
			hist(flow_y_est-flow_y_true, main="Difference between estimated and true flows (y components) [m/sec]")
		}





