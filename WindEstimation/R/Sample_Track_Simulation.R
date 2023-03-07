#Please install 'circular' package
install.packages("circular", dependencies = TRUE)
library(circular)
set.seed(1234)


######################################################
# This code simulates track of an animal which travels in moving fluid.  																					   
#																																							   
# After simulation, this code outputs following data.																										   
# The outputted data "Simulation_Track.csv" can be
#inputted to 																				                   
# "Heading_Wind_Estimation.R"(The R code for estimating #heading and wind vector from track data). 															   
#																																							   
#              timepoint: time point of data																											           
#             trackspeed: speed of track vector [m/sec]																										   
#         trackdirection: direction of track vector [rad]																									   
#                     dt: elapsed time from  the fix recorded at previous time point [sec]                  												       
#                      X: x component of fix [m]																										       
#                      Y: y component of fix [m]																							                   
#  mean_heading_speed: mean speed of heading vector [m/s]																								   
# mean_heading_direction: mean direction of heading #vector [m/s]																						           
#         mean_flow_x: x component of mean flow vecotr [m/s]																						            
#         mean_flow_y: y component of mean flow vecotr [m/s]																								   
#																																					           
# Also, this code plots the simulated track and #distribution of track vecotrs 																				   
#(green dots:track vector, black arrow: mean track vector, #red arrow: mean heading vecotr, blue arrow: mean flow #vector).					                   
				     #####################################################
### Please input following parameter values.  #####################################################
                                                                                                                                                              
# Please specify the place the data will be outputted																										    
  outfile <- paste("~/Desktop/Simulation_Track.csv")																										   				
																																							   
# Parameters
                                                                                                                                                                                                                                                                                                               
#              dt: sampling interval [sec]																													   
#               n: sample size of simulated data(number of #ground velocity. The number of fixes is n+1.)														   
#            gv_x: x component of mean track vector [m/sec]																									   
#            gv_y: y component of mean track vector [m/sec]																									   
#  mean_air_speed: mean air speed [m/sec]																												       
#              mu: mean heading direction [radian]  (East:0 #radian,  North:pi/2 radian, West:pi radian, #South:pi*3/2 radian)									   
#         sa_para: SD of heading vecotr along the mean #vector(see Supplementary Note 2)																	       
#         sa_orth: SD of heading vecotr perpendicular #direction relative to the mean vector(see Supplementary #Note 2)										   
#																																							   
#       windnoise: SD of wind vector  
# windnoise[[1]]:SD of #wind vector along the mean vector																	   
# windnoise[[2]]:SD of the wind vector perpendicular #direction relative to the mean vector (see Supplementary #Note 3)*     
																																							   
# obs_noise_ratio: Ratio between "SD of observation error" #and "mean distance between successive fixes"(see #Supplementary Note 3)**                            
#																																							   
#   *windnoise used in Supplementary Note 3                                                                                                                    
# ( c(0.25,0.25), c(0.5,0.5), c(1,1), c(0.25,0.5), #c(0.5,1), c(1,2), c(0.5,0.25), c(1,0.5), c(2,1) )   														   
# **obs_noise_ratio used in Supplementary Note 3 																											   
# (0, 0.05, 0.1, 0.15) 																																		   
																																							   
    samplingtime <- 60																												                           
              n <- 100 																																		   
	       gv_x <-  10																																		   
           gv_y <-   0  																																		   
 mean_air_speed <-  34.7/3.6																																   
		    mu <-  pi/4																																		   
        sa_para <-   0.5      																																   
        sa_orth <-   1.5																																   
      windnoise <- c(0.0,0.0)																																   
obs_noise_ratio <- 0.0                        																												   
											
																												   
######################################################
# After checking the above parameter values, please conduct
#this code.      #####################################################



### calculate shape parameter from mean and sd ###
init_a_numeric<-function(x,v_mean,v_sd){
	c(v_mean^2*(gamma(1+2/x)-gamma(1+1/x)^2)-v_sd*v_sd*gamma(1+1/x)^2)
}
##################################################

### calucuration of parameters

 #calucurate shape parameter of Weibull distribution
   fn<-function(x) init_a_numeric(x, mean_air_speed, sa_para)
  ans<-uniroot(fn,c(1,500))
  a<-ans$root
 #a: shape parameter of Weibull distribution (Air speed)

  b<- mean_air_speed/gamma(1+1/a)   
 #b: scale parameter of Weibull distribution (Air speed)

  kappa <- (mean_air_speed/sa_orth)^2
 #kappa: concentration parameter of von Mises distribution #(Heading) 
 
s_obs <- mean_air_speed*obs_noise_ratio
#s_obs:observation noise SD
s_ws_para  <-windnoise[[1]]
#s_ws_para: wind speed noise SD
s_ws_orth  <-windnoise[[2]]
#s_ws_orth: wind perpendicular noise SD

mean_wx<-gv_x-mean_air_speed*cos(mu)
mean_wy<-gv_y-mean_air_speed*sin(mu)
mean_wind_direction <-atan2( mean_wy,mean_wx)
mean_wind_speed      <- sqrt( mean_wx^2 + mean_wy^2 )

#### end of calcuration of parameters


### Wind vector (w_x, w_y) #####################################################

wx_s_list<-rnorm(n,0,s_ws_para)
wy_s_list<-rnorm(n,0,s_ws_orth)
wx_sr_list<-cos(mean_wind_direction)*wx_s_list-sin(mean_wind_direction)*wy_s_list
wy_sr_list<-sin(mean_wind_direction)*wx_s_list+cos(mean_wind_direction)*wy_s_list
#wx_sr_list: x component of noise of wind vector
#wy_sr_list: y component of noise of wind vector

w_x<- mean_wx+wx_sr_list
w_y<-mean_wy+wy_sr_list
# w_x: x component of wind vector
# w_y: y component of wind vector
######################################################


### Heading vecotr (va_x, va_y) #################################################
heading    <- as.numeric(rvonmises(n,mu,kappa))
air_speed <- rweibull(n,shape=a,scale=b)

va_x<- air_speed*cos(heading)
va_y<- air_speed*sin(heading)
# va_x: x component of heading vector
# va_y: y component of heading vector
######################################################
### "Observed" track vecotr ######################################################
###     g_speed_obs: speed of track vector [m/sec]     ###########################
### g_direction_obs: direction of track vecotr [m/sec] ###########################

vg_x<- air_speed*cos(heading)+ w_x
vg_y<- air_speed*sin(heading)+ w_y
#true ground velocity

x<-c(0:n)
y<-c(0:n)
for(t in 1:n){
	x[t+1] <- x[t]+vg_x[t]
	y[t+1] <- y[t]+vg_y[t]
	}
#true animal fixes

x_obs <-x + rnorm(length(x),0,s_obs)
y_obs<- y + rnorm(length(y),0,s_obs)
#observed fixes

vg_x_obs <- x_obs[2:(n+1)]- x_obs[1:n]
vg_y_obs <- y_obs[2:(n+1)]- y_obs[1:n]
#observed ground velocity

g_speed_obs      <- sqrt(vg_x_obs^2 + vg_y_obs^2)
#observed ground velocity speed
g_direction_obs <- atan2(vg_y_obs,vg_x_obs)
# observed ground velocity direction

######################################################

#################################################################   End of the simulation                                                   
######################################################


### plot the simulation results
edge <- max( abs(y_obs- mean(y_obs)),abs(x_obs- mean(x_obs)))+20
xmean<-mean(x_obs)
ymean<-mean(y_obs)
dev.new()
plot(x_obs,y_obs, xlim=c(xmean-edge,xmean+edge), ylim=c(ymean-edge,ymean+edge),asp=1, main="Simulated track")

scale <- 1.5*sqrt(gv_x^2+gv_y^2)
track_vector_x<-vg_x_obs
track_vector_y<-vg_y_obs
dev.new()
plot(track_vector_x,track_vector_y , xlim=c(-scale,scale), ylim=c(-scale,scale),asp=1, col="green", main="Track vector distribution of simulated track")
arrows(0, 0, mean_wx, mean_wy, col="blue")
#arrows(0, 0, mean_air_speed*cos(mu), mean_air_speed*sin(mu), col="red") 
arrows(mean_wx, mean_wy, mean_air_speed*cos(mu)+mean_wx, mean_air_speed*sin(mu)+mean_wy, col="red") 
arrows(0, 0, gv_x, gv_y, col="black")

### output the simulation results
outdata <- data.frame(timepoint=c(1:n), track_speed =g_speed_obs, track_direction = g_direction_obs,dt=rep(samplingtime,n),X=x_obs[2:(n+1)],Y=y_obs[2:(n+1)],mean_head_speed=rep(mean_air_speed,n),mean_head_direction=rep(mu,n), mean_flow_x=rep(mean_wx,n), mean_flow_y=rep(mean_wy,n))
write.table(outdata, outfile,quote =F, row.names=F,sep=", ")
