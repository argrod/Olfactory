library(ggplot2)
library(ggforce)
library(extrafont)

if(Sys.info()['sysname'] == "Darwin"){
  figLoc <- "/Volumes/GoogleDrive/My Drive/PhD/Figures/Olfactory/"
  # figLoc <- "/Documents/GitHub/PhD/Olfactory/"
} else {
  figLoc <- "G:/UTokyoDrive/PhD/Figures/Olfactory/"
  # figLoc <- "F:/Documents/GitHub/PhD/Olfactory/"
}

# run MinuteHeading_WindEstimation with these variables
a = 1
center = 1149
# which(as.character(tp) == "2018-09-03 18:33:34")

    yoko<-Von_Mises_sd(sqrt(ans_best$par[3]*ans_best$par[3]+ans_best$par[2]*ans_best$par[2]))*Weibull_mean(ans_best$par[1],constv/gamma(1+1/ans_best$par[1]))
#    yoko: SD of the heading vector distribution to the
#perpendicular direction relative to the mean direction	
#     
					tate<-Weibull_sd(ans_best$par[1],constv/gamma(1+1/ans_best$par[1]))
#    tate: SD of the heading vector distribution along the mean direction

    	nrp<-ks.test(nr,"pweibull",ans_best$par[1],constv/gamma(1+1/ans_best$par[1]))

#direction of heading vecotr
								mu<-atan2(ans_best$par[3],ans_best$par[2])
								kappa<-sqrt(ans_best$par[2]*ans_best$par[2]+ans_best$par[3]*ans_best$par[3])
								ndp<-ks.test(nd,"pvonmises",mu,kappa)

plot(long[st:end],lat[st:end])

# find mean heading
# split in x and y
xComp <- rrow[st:end] * cos(drow[st:end])
yComp <- rrow[st:end] * sin(drow[st:end])



plot(mean(xComp),mean(yComp))
ggplot() + 
    geom_circle(aes(x0 = 0, y0 = 0, r = sqrt((mean(xComp) + ans_best$par[4])^2 + (mean(yComp) + ans_best$par[5])^2)),size = 1.5,colour="#225ea8") +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      text = element_text(size = 10),axis.text = element_text(size = 8)) +
    scale_x_continuous(name = expression(paste("U component (m", s^{-1},")"))) +
    scale_y_continuous(name = expression(paste("V component (m", s^{-1},")"))) +
    geom_circle(aes(x0 = 0, y0 = 0, r = sqrt(mean(xComp)^2 + mean(yComp)^2)),size = 1.5,colour="#238443") +
    geom_segment(aes(x = 0, y = 0, xend = mean(xComp), yend = mean(yComp)),
    arrow = arrow(length = unit(0.03, "npc"),type="closed"),size=1.1,linejoin="round") +
    geom_segment(aes(x = mean(xComp), y = mean(yComp), xend = mean(xComp) + ans_best$par[4], yend = mean(yComp) + ans_best$par[5]),
    arrow = arrow(length = unit(0.03, "npc"),type="closed"), colour = "blue",size=1,linejoin="round") +
    geom_segment(aes(x = 0, y = 0, xend = mean(xComp) + ans_best$par[4], yend = mean(yComp) + ans_best$par[5]),
    arrow = arrow(length = unit(0.03, "npc"),type="closed"), colour= "#ef3b2c",size=1,linejoin="round")

ggsave(paste(figLoc,"DistMethod.svg",sep=""), device="svg", dpi = 300, height = 2,
      width = 2, units = "in")
dev.off()

ggplot() + geom_path(aes(x = long[st:end], y = lat[st:end])) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 10),axis.text = element_text(size = 8)) +
  scale_x_continuous(name = expression(paste("Longitude ("*degree*")")), breaks = c(142.55,142.600)) + 
  scale_y_continuous(name = expression(paste("Latitude ("*degree*")")),breaks = c(41.0,41.1,41.2)) + 
  geom_segment(aes(x=long[seq(st,end-1,round((end-st-1)/5))],
  y = lat[seq(st,end,round((end-st-1)/5))],
  xend = long[seq(st+1,end,round((end-st+1)/5))],
  yend = lat[seq(st+1,end,round((end-st+1)/5))]),
  arrow = arrow(length = unit(0.1,"inches")))

ggsave(paste(figLoc,"DistSegment.svg",sep=""), device="svg", dpi = 300, height = 2,
      width = 2, units = "in")
dev.off()

ggplot() + geom_point(aes(x = d, y = r)) + 
    scale_x_continuous(limits=c(-pi,pi))

ggplot() + geom_point(aes(x = drow[st:end], y = rrow[st:end])) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 10),axis.text = element_text(size = 8)) +
  scale_x_continuous(name = expression(paste("Headings (rad)"))) + 
  scale_y_continuous(name = expression(paste("Speed (m"*s^{-1}*")")))

ggsave(paste(figLoc,"DistSpdHd.svg",sep=""), device="svg", dpi = 300, height = 2,
      width = 2, units = "in")
dev.off()
