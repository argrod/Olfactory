using RData
using DataFrames
using Plots

wEst = load("/Volumes/GoogleDrive/My Drive/PhD/Data/WindCalc/windDat.RData", convert=true)
wEst = DataFrame(wEst["WindDat"])

plot(wEst.RelHead.*(180/pi),wEst.distTo,seriestype=:scatter, proj=:polar, m=:red, bg=:black)