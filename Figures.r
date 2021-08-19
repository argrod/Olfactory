	
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
    geom_segment(aes(x = 0, y = 0, xend = mean(xComp), yend = mean(yComp)),
    arrow = arrow(length = unit(0.03, "npc"))) +
    geom_segment(aes(xend = mean(xComp), yend = mean(yComp), x = mean(xComp) + ans_best$par[4], y = mean(yComp) + ans_best$par[5]),
    arrow = arrow(length = unit(0.03, "npc")), colour = "blue") +
    geom_segment(aes(x = 0, y = 0, xend = mean(xComp) + ans_best$par[4], yend = mean(yComp) + ans_best$par[5]),
    arrow = arrow(length = unit(0.03, "npc")), colour = "red")

library(ggplot2)
ggplot() + geom_path(aes(x = long[st:end], y = lat[st:end]))

ggplot() + geom_point(aes(x = d, y = r)) + 
    scale_x_continuous(limits=c(-pi,pi))

ggplot() + geom_point(aes(x = drow[st:end], y = rrow[st:end]))
inithd_first


## m=matrix(data=sample(rnorm(100,mean=0,sd=2)), ncol=10)
## this function makes a graphically appealing heatmap (no dendrogram) using ggplot
## whilst it contains fewer options than gplots::heatmap.2 I prefer its style and flexibility
 
ggheat=function(m, rescaling='none', clustering='none', labCol=T, labRow=T, border=FALSE, 
heatscale= c(low='blue',high='red'))
{
  ## the function can be be viewed as a two step process
  ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
  ## using simple options or by a user supplied function
  ## 2. with the now resahped data the plot, the chosen labels and plot style are built
 
  require(reshape)
  require(ggplot2)
 
  ## you can either scale by row or column not both! 
  ## if you wish to scale by both or use a differen scale method then simply supply a scale
  ## function instead NB scale is a base funct
 
  if(is.function(rescaling))
  { 
    m=rescaling(m)
  } 
  else 
  {
    if(rescaling=='column') 
      m=scale(m, center=T)
    if(rescaling=='row') 
      m=t(scale(t(m),center=T))
  }
 
  ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
  ## if you want a different distance/cluster method-- or to cluster and then scale
  ## then you can supply a custom function 
 
  if(is.function(clustering)) 
  {
    m=clustering(m)
  }else
  {
  if(clustering=='row')
    m=m[hclust(dist(m))$order, ]
  if(clustering=='column')  
    m=m[,hclust(dist(t(m)))$order]
  if(clustering=='both')
    m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
  }
    ## this is just reshaping into a ggplot format matrix and making a ggplot layer
 
  rows=dim(m)[1]
  cols=dim(m)[2]
  melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
    g=ggplot(data=melt.m)
 
  ## add the heat tiles with or without a white border for clarity
 
  if(border==TRUE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
  if(border==FALSE)
    g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
 
  ## add axis labels either supplied or from the colnames rownames of the matrix
 
  if(labCol==T) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
  if(labCol==F) 
    g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
 
  if(labRow==T) 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))    
    if(labRow==F) 
    g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))    
 
  ## get rid of grey panel background and gridlines
 
  g2=g2+opts(panel.grid.minor=theme_line(colour=NA), panel.grid.major=theme_line(colour=NA),
  panel.background=theme_rect(fill=NA, colour=NA))
 
  ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
  return(g2+scale_fill_continuous("", heatscale[1], heatscale[2]))
 
}
 
  ## NB because ggheat returns an ordinary ggplot you can add ggplot tweaks post-production e.g. 
  ## data(mtcars)
  ## x= as.matrix(mtcars)
  ## ggheat(x, clustCol=T)+ opts(panel.background=theme_rect(fill='pink'))