using Pkg
Pkg.add(["GRIB","DataFrames","Query","CSV","Geodesy","Plots","StructArrays","RecursiveArrayTools","Statistics","RCall","NetCDF"."DelimitedFiles","TimeZones","JLD","DataDeps"])
using DataFrames
using Query
using GRIB
using Dates
using CSV
using Geodesy
using Plots
using StructArrays
using RecursiveArrayTools
using Statistics
using RCall
using NetCDF
using DelimitedFiles
using TimeZones
using JLD
using DataDeps
using GMT

## FUNCTIONS
## FUNCTION TO FIND THE INDECES OF SUITABLE WIND ESTIMATIONS
function gribCond(dt, h, md)
    hr = Dates.hour(round(dt, Dates.Hour))
    min = Dates.minute(dt)
    iszero(hr%h) && (min > 60 - md || min < md)
end
## FUNCTION TO DEFINE THE NEAREST h HOUR VALUE
function gribCondGT(dt, h)
    hrs = (0:23)[iszero.((0:23).%h) .== 1]
    hr = Dates.hour(round(dt, Dates.Hour))
    chosen = maximum(hrs[hrs .<= hr])
    DateTime(Int(Dates.year(dt)),Int(Dates.month(dt)),Int(Dates.day(dt)),chosen,0,0)
end
## GENERATE GRIB FILENAME FROM DATETIME
function gribGen(dt)
    yr = string(Dates.year(round(dt, Dates.Hour)))
    mn = string(Dates.month(round(dt, Dates.Hour)))
    if length(mn) == 1
        mn = "0"*mn
    end
    dy = string(Dates.day(round(dt, Dates.Hour)))
    if length(dy) == 1
        dy = "0"*dy
    end
    hr = string(Dates.hour(round(dt, Dates.Hour)))
    if length(hr) == 1
        hr = "0"*hr
    end
    "Z__C_RJTD_"*yr*mn*dy*hr*"0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"
end
##
function roundF(vals, rnum)
    rnum*(round(vals)/rnum)
end
## GENERATE STRING OF GRIB TIMES
function gribT(dt)
    yr = string(Dates.year(round(dt, Dates.Hour)))
    mn = string(Dates.month(round(dt, Dates.Hour)))
    if length(mn) == 1
        mn = "0"*mn
    end
    dy = string(Dates.day(round(dt, Dates.Hour)))
    if length(dy) == 1
        dy = "0"*dy
    end
    hr = string(Dates.hour(round(dt, Dates.Hour)))
    if length(hr) == 1
        hr = "0"*hr
    end
    return yr*mn*dy*hr
end
## EXTRACT VALUE AND LAT LONS FROM GRIBFILE (MUST BE IN GRIB DIRECTORY)
function ConanThe(GribN)
    msg = GribFile(GribN) do f
        Message(f)
    end
    keylist = Vector{String}()
    for key in keys(msg)
        key != "values" && key != "latitudes" && key != "longitudes" && push!(keylist, key)
    end
    keylist = unique(keylist)
    df = DataFrame((Symbol(key) => Any[] for key in keylist)...)
    GribFile(GribN) do f
        for message in f
            push!(df, tuple((message[key] for key in keylist)...))
        end
    end
    df
end
## EXTRACT 10M EASTING WIND VECTOR
function xtractUWind00(GribN)
    GribFile(GribN) do f
        for message in f
            message["forecastTime"] == 0 && message["shortName"] == "10u" && return message
        end
    end
end
## EXTRACT 10M NORTHING WIND VECTOR
function xtractVWind00(GribN)
    GribFile(GribN) do f
        for message in f
            message["forecastTime"] == 0 && message["shortName"] == "10v" && return message
        end
    end
end
## EXTRACT GRIBDATA FOR SELECTED TIME
function gribSelect(dt, lat, lon)
    gribob = gribGen(dt)
    # fls = readdir()
    # if any(fls.==gribob)
    #     println("Failed on ", gribob)
    # end
    UWind = @NamedTuple{lon, lat, value}(data(xtractUWind00(gribob)))
    VWind = @NamedTuple{lon, lat, value}(data(xtractVWind00(gribob)))
    nrlon = unique(UWind.lon[(abs.(UWind.lon .- roundF(lon,0.0625)) .== min(abs.(UWind.lon .- roundF(lon,0.0625))...))])
    nrlat = unique(UWind.lat[(abs.(UWind.lat .- roundF(lat,0.05)) .== min(abs.(UWind.lat .- roundF(lat,0.05))...))])
    return UWind.value[@. (UWind.lon .== nrlon) & (UWind.lat .== nrlat)]..., VWind.value[@. (UWind.lon .== nrlon) & (UWind.lat .== nrlat)]..., nrlat, nrlon
end
## FIND LAT AND LON WITHIN 5 KM OF TRACK WITHIN tgap MINS OF t
function findGribPoints(t,tgap,buff)
    # isolate points along track to find lat lon ±50mins
    isol = EDat[(t - Minute(tgap)) .< EDat.DT .< (t + Minute(tgap)),:]
    @rput isol buff
    R"""
    library(sp)
    library(raster)
    library(plyr)
    for (group in isol$ID){
        df <- isol[isol$ID == group,]
        pnt <- SpatialPoints(cbind(df$Lon,df$Lat),proj4string=CRS("+proj=longlat"))
        out <- data.frame(lon=numeric(),lat=numeric())
        pntBuff <- buffer(pnt,buff)
        # create grid of grib lat lons
        bndris <- as.data.frame(bbox(pntBuff[1][1,1]))
        gribGrid <- SpatialPoints(expand.grid(seq(round_any(bndris$min[1],0.0625),round_any(bndris$max[1],0.0625),0.0625),
            seq(round_any(bndris$min[2],0.05),round_any(bndris$max[2],0.05),0.05)),proj4string=CRS("+proj=longlat"))
        out <- rbind(out,coordinates(gribGrid[which(over(gribGrid, pntBuff) == 1)]))
        }
    colnames(out) <- c("lon","lat")
    """
    @rget out df
    return out, df
end
# round to nearest n
function roundNearest(x,n)
    return round(x/n) * n
end
# find nearest 
# round lat lon to nearest 0.05 and 0.0625 respectively
lon = roundNearest(lon,0.0625)
lat = roundNearest(lat,0.05)

a=1
isol = EDat[(t[a] - Minute(51)) .< EDat.DT .< (t[a] + Minute(51)),:]
Plots.scatter()
for ids in unique(isol.ID)
    Plots.plot!(isol.Lon[isol.ID .== ids],isol.Lat[isol.ID .== ids])
end

isol.ID .== "10_S1"

## FIND LAT AND LON WITHIN buff M OF TRACK
function findGribPointsGT(ts,buff)
    # isolate points along track to find lat lon ±50mins
    isol = EDat[EDat.RoundDT .== ts,:]
    isol = EDat
    lon = isol.Lon
    lat = isol.Lat
    @rput lon lat buff
    R"""
    library(sp)
    library(raster)
    library(plyr)
    pnts <- SpatialPoints(cbind(lon,lat),proj4string=CRS("+proj=longlat"))
    pntsBuff <- buffer(pnts, buff)
    # plot(pntsBuff)
    # points(pnts,cex=.1)
    # create grid of grib lat lons
    bndris <- as.data.frame(bbox(pntsBuff[1][1,1]))
    gribGrid <- SpatialPoints(expand.grid(seq(round_any(bndris$min[1],0.0625),round_any(bndris$max[1],0.0625),0.0625),
        seq(round_any(bndris$min[2],0.05),round_any(bndris$max[2],0.05),0.05)),proj4string=CRS("+proj=longlat"))
    out <- as.data.frame(coordinates(gribGrid[which(over(gribGrid, pntsBuff) == 1)]))
    colnames(out) <- c("lon","lat")
    """
    @rget out
    return out
end
## FIND AVERAGE WIND HEADING WITHIN 5KM FOR GIVEN DATETIME LAT LON WITHIN tgap MINUTES
function gribAve(tPoint,tgap, buff)
    iso = EDat[(tPoint - Minute(tgap)) .< EDat.DT .< (tPoint + Minute(tgap)),:]
    gribPts = findGribPoints.(iso.Lon, iso.Lat, buff)
    # gribD = gribSelect.(fill(dt,nrow(gribPts)), gribPts.lat, gribPts.lon)
    # U = reduce(vcat,[x[1] for x in gribD])
    # V = reduce(vcat,[x[2] for x in gribD])
    # heads = atan(sin(mean(U))/cos(mean(V)))
    # # heads = atan.(sum(U),sum(V))
    # spds = mean(sqrt.(U.^2 + V.^2))
    return gribPts
end
## CALCULATE AVERAGE HEADING OF BIRDS
function aveBirdHd(head)
    atan(mean(sin.(head)),mean(cos.(head)))
end
## CALCULATE AVERAGE HEADINGS FROM X Y COMPONENTS
function aveWind(X, Y)
    atan(mean(X),mean(Y))
end
## CALCULATE AVERAGE ESTIMATED WIND WITHIN tgap MINUTES OF TIME tPoints
function birdWindow(tPoints,tgap)
    iso = EDat[(tPoints - Minute(tgap)) .< EDat.DT .< (tPoints + Minute(tgap)),:]
    isol = filter(:DT => x -> gribCond(x, 3, tgap), iso)
    aveWHd = aveWind(isol.X,isol.Y)
    # iso = EDat[(tPoints[b] - Minute(50)) .< EDat.DT .< (tPoints[b] + Minute(50)),:]
    # isol = filter(:DT => x -> gribCond(x, 3, 50), iso)
    # gribAve.(isol.DT, isol.Lon,isol.Lat)
    return aveWHd
end
## CALCULATE AVERAGE ESTIMATED WIND FOR ALL POINTS WITHIN 3 HOUR BLOCK
function birdWindowGT(tPoints)
    iso = EDat[EDat.RoundDT .== tPoints,:]
    aveWHd = aveWind(iso.X,iso.Y)
    # iso = EDat[(tPoints[b] - Minute(50)) .< EDat.DT .< (tPoints[b] + Minute(50)),:]
    # isol = filter(:DT => x -> gribCond(x, 3, 50), iso)
    # gribAve.(isol.DT, isol.Lon,isol.Lat)
    return aveWHd
end
##
function birdWSpeed(tPoints,tgap)
    iso = EDat[(tPoints - Minute(tgap)) .< EDat.DT .< (tPoints + Minute(tgap)),:]
    isol = filter(:DT => x -> gribCond(x, 3, tgap), iso)
    aveSpd = sqrt(mean(isol.X)^2 + mean(isol.Y)^2)
    return aveSpd
end
# ##
function birdWSpeedGT(tPoints)
    iso = EDat[EDat.RoundDT .== tPoints,:]
    aveSpd = sqrt(mean(iso.X)^2 + mean(iso.Y)^2)
    return aveSpd
end 

# function WindVal(t, lat, lon)
#     # calculate mean estimated heading for ±25 min
#     estHdAve = birdWindow(t)
#     # find all bird locations within ±50 mins and all grib locations within 5km
#     gribloc = findGribPoints(EDat, t)
# end

function readPyOut(filename)
    test = CSV.File(filename,header=true) |> DataFrame
    test.ID = fill(match(r".*PythonWinds/(.*).csv",filename).captures[1],nrow(test))
    return test
end

############################################################################################################
######################################## TEST PYTHON WIND ESTIMATES ########################################
############################################################################################################

# load in the wind calculations
estLoc = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/PythonWinds/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = filter(x->contains(x,".csv"), readdir(estLoc,join=true))
EDat = vcat(readPyOut.(EFiles)...)

## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
rename!(EDat, [:Time,:Lat,:Lon,:Head,:X,:Y,:ID])
df = dateformat"y-m-d H:M:S"
EDat.DT = DateTime.(EDat.Time, df)
sel = filter(:DT => x -> gribCond(x, 3, 25), EDat)

gribsels = unique(gribGen.(sel.DT))

# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
##
WLoc = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/gribs/"
cd(WLoc)
# @rput gribsels

# gribTimes = unique(gribT.(sel.DT))
# @rput gribTimes

# # surely at this stage there is a download link capacity within Julia - avoid the change to R where possible

# ## DOWNLOAD THE REQUIRED GRIBFILES
# R"""
# options(timeout=1000)
# dloadLoc  = "/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/gribs/"
# gribFls = list.files(dloadLoc,pattern="*grib2.bin")
# for(b in 1:length(gribsels)){
#     if(any(gribFls == gribsels[b])){
#         next
#     } else {
#         download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
#             destfile = paste(dloadLoc,gribsels[b], sep = ""))
#         }
# }
# """

# outfile = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/Temp.txt"
# open(outfile, "w") do f
#   writedlm(f, gribsels,',')
# end

outPt = DataFrame();

Plots.scatter(ylims=(-π,π),xlims=(-π,π))

compT = filter(:DT => x -> gribCond(x, 3, 25), EDat)
compT.wHead = atan.(compT.Y,compT.X)
compT.wSpd = sqrt.(compT.Y.^2 + compT.X.^2)
# if nrow(compT) > 0
    t = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd(WLoc)
    gribLs = findGribPoints.(t,50,5000)
    # aveHeds = zeros(length(t),1)
    # aveSpds = zeros(length(t),1)
    for tPoint = 1:length(t)
        gribD = gribSelect.(fill(t[tPoint],nrow(gribLs[tPoint][1])),gribLs[tPoint][1].lat,gribLs[tPoint][1].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds = atan(mean(U),mean(V))
        aveSpds = sqrt(mean(U)^2 + mean(V)^2)
        append!(outPt,DataFrame(hcat(t[tPoint], aveHeds, aveSpds, atan(mean(gribLs[tPoint][2].X),mean(gribLs[tPoint][2].Y)),sqrt(mean(gribLs[tPoint][2].X)^2 + mean(gribLs[tPoint][2].Y)^2)), :auto))
    end
    # rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed])
    # outPt = vcat(outPt, add)
# end
# end
rename!(outPt, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed])

Plots.scatter(outPt.gribSpeed,outPt.estSpeed,label=nothing,xaxis="JMA wind speed (m/s)",yaxis="Python-estimated wind speed (m/s)")
Plots.plot!(0:10,0:10,label=nothing)
Plots.scatter(outPt.gribHead,outPt.estHead,label=nothing,xaxis="JMA wind direction (rad)",yaxis="Python-estimated wind direction (rad)")
Plots.plot!(-π:π,-π:π,label=nothing)

CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/Validation.csv", outPt)

for x in 1:length(t)
    isol = EDat[(t[x] - Minute(50)) .< EDat.DT .< (t[x] + Minute(50)),:]
    lon = isol.Lon
    lat = isol.Lat
    CSV.write("/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/" * "csvOut" * string(x) * ".txt",DataFrame(lon=lon,lat=lat))
end

using Glob
inFiles = glob("*rOut.csv","/Users/aran/Library/CloudStorage/GoogleDrive-a-garrod@g.ecc.u-tokyo.ac.jp/My Drive/PD/Data/TestingData/")

gribLs = [CSV.File(inFiles[x],header=true) |> DataFrame for x in 1:length(inFiles)]


gribD = gribSelect.(sel.DT, sel.Lat, sel.Lon)
##
sel.U = reduce(vcat,[x[1] for x in gribD])
sel.V = reduce(vcat,[x[2] for x in gribD])
sel.WindHead = atan.(sel.V,sel.U)
sel.EHead = atan.(sel.Y,sel.X)
# plot(sel.WindHead, sel.Head, seriestype = :scatter)
# plot(sqrt.(sel.X.^2 + sel.Y.^2), sqrt.(sel.U.^2 + sel.V.^2), seriestype = :scatter)
sel.ESpd = sqrt.(sel.X.^2 + sel.Y.^2)
sel.WSpd = sqrt.(sel.U.^2 + sel.V.^2)
CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv", outPt)
##
@rput outPt
R"""
library(circular)
library(ggplot2)
res<-cor.circular(outPt$estHead,outPt$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
ggplot(outPt, aes(x = EHead, y = gribHead)) +
    geom_point() #+
    # geom_line(aes(x=-3:3,y=-3:3))
res
# spres <- circular::cor.test(sel$ESpd, sel$WSpd)
"""
## 
outPt
plot(sel.WindHead, sel.EHead, seriestype = scatter)


plot(sel.WindHead,sel.Head, seriestype = scatter)
plot(sel.ESpd,sel.WSpd, seriestype = scatter)

plot(outPt.estHead,outPt.gribHead,seriestype=scatter)
plot!(-3:3,-3:3)

plot(outPt.gribSpeed,outPt.estSpeed,seriestype=scatter)
plot!(0:15,0:15)

################################################################################################################
########################################### RUN FOR YONEHARA METHOD ############################################
################################################################################################################

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = DataFrame()
outPt = DataFrame()

for b = 1:length(EFiles)
    append!(EDat,[CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2019-"*EFiles[b][1:end-4],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
end
rename!(EDat, [:DT,:Lat,:Lon,:BirdHead,:WDir,:WSpeed,:Resnorm,:ID]) 
df = dateformat"y/m/d,H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
sel = filter(:DT => x -> gribCond(x, 3, 25), EDat)
gribsels = unique(gribGen.(sel.DT))

# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
# sel.GLat = roundF.(sel.Lat, .05)
# sel.GLon = roundF.(sel.Lon, .0625)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
cd(WLoc)
@rput gribsels

gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
## DOWNLOAD THE REQUIRED GRIBFILES
R"""
# save(list=c("gribsels","gribTimes"), file = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
# load("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
options(timeout=800)
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
for(b in 1:length(gribsels)){
    if(any(gribFls == gribsels[b])){
        next
    } else {
        download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
            destfile = paste(dloadLoc,gribsels[b], sep = ""))
        }
}
"""

outPt = DataFrame()
for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2019-"*EFiles[b][1:end-4],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:BirdHead,:WDir,:WSpeed,:Resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!(EDat.WSpeed .== 0),:]
    # convert to DateTime
    df = dateformat"y/m/d,H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        if length(gribD) > 0
            U = reduce(vcat,[x[1] for x in gribD])
            V = reduce(vcat,[x[2] for x in gribD])
            aveHeds[tPoint] = atan(mean(U),mean(V))
            aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
        else
            aveHeds[tPoint] = NaN
            aveSpds[tPoint] = NaN
        end
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)  
end
extra = outPt
## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
# rename!(EDat, [:DT,:Head,:X,:Y,:Lat,:Lon])
rename!(EDat, [:DT,:Lat,:Lon,:Forage,:DistTo,:BirdHead,:WDir,:WSpeed,:ID]) 
df = dateformat"d-m-y H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
Ll = LLA.(EDat.Lat,EDat.Lon)
UTMD = UTMZfromLLA(wgs84)
utmDat = map(UTMD, Ll); utmDat = StructArray(utmDat)
X = diff(utmDat.x)
Y = diff(utmDat.y)
head = atan.(X,Y)

CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMethodValidation.csv", outPt)

gribD = gribSelect.(sel.DT, sel.Lat, sel.Lon)
##
sel.U = reduce(vcat,[x[1] for x in gribD])
sel.V = reduce(vcat,[x[2] for x in gribD])
sel.WindHead = atan.(sel.V,sel.U)
sel.EHead = atan.(sel.Y,sel.X)
# plot(sel.WindHead, sel.Head, seriestype = :scatter)
# plot(sqrt.(sel.X.^2 + sel.Y.^2), sqrt.(sel.U.^2 + sel.V.^2), seriestype = :scatter)
sel.ESpd = sqrt.(sel.X.^2 + sel.Y.^2)
sel.WSpd = sqrt.(sel.U.^2 + sel.V.^2)
CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv", outPt)
##
@rput outPt
R"""
library(circular)
library(ggplot2)
res<-cor.circular(outPt$estHead,outPt$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
#ggplot(outPt, aes(x = EHead, y = gribHead)) +
 #   geom_point() #+
    # geom_line(aes(x=-3:3,y=-3:3))
summary(res)
# spres <- circular::cor.test(sel$ESpd, sel$WSpd)
"""
## 
outPt
plot(sel.WindHead, sel.EHead, seriestype = scatter)


plot(sel.WindHead,sel.Head, seriestype = scatter)
plot(sel.ESpd,sel.WSpd, seriestype = scatter)

plot(outPt.estHead,outPt.gribHead,seriestype=scatter,legend=false,xlabel="Estimated wind direction ( rad)", ylabel="Validation wind direction (rad)")
plot!(-pi:pi,-pi:pi)

plot(outPt.gribSpeed,outPt.estSpeed,seriestype=scatter,legend=false,xlabel="Estimated wind speed (m/s)", ylabel="Validation wind speed (m/s)",fmt=:png)
plot!(0:10,0:10)


# attempt reading in .nc file from sar-winds

#############################################################################################
######################## TESTING YONEHARA METHOD WITH 3 HOUR AVERAGE ########################
#############################################################################################

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = DataFrame()
outPt = DataFrame()

for b = 1:length(EFiles)
    append!(EDat,[CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2019-"*EFiles[b][1:end-4],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
end
rename!(EDat, [:DT,:Lat,:Lon,:BirdHead,:WDir,:WSpeed,:Resnorm,:ID]) 
df = dateformat"y/m/d,H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
# sel = filter(:DT => x -> gribCond(x, 3, 90), EDat)
gribsels = unique(gribGen.(gribCondGT.(EDat.DT,3)))

# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
# sel.GLat = roundF.(sel.Lat, .05)
# sel.GLon = roundF.(sel.Lon, .0625)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
cd(WLoc)
@rput gribsels

gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
## DOWNLOAD THE REQUIRED GRIBFILES
R"""
# save(list=c("gribsels","gribTimes"), file = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
# load("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
options(timeout=800)
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
for(b in 1:length(gribsels)){
    if(any(gribFls == gribsels[b])){
        next
    } else {
        download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
            destfile = paste(dloadLoc,gribsels[b], sep = ""))
        }
}
"""

outPt = DataFrame()
# save("/Volumes/GoogleDrive/My Drive/PhD/Data/3HourOutpt.jld","outPt",outPt)
load("/Volumes/GoogleDrive/My Drive/PhD/Data/3HourOutpt.jld")
# nrow(outPt)
# 26+27+27+25+31+25+32+23

# toAdd = DataFrame()
for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2019-"*EFiles[b][1:end-4],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:BirdHead,:WDir,:WSpeed,:Resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!(EDat.WSpeed .== 0),:]
    # convert to DateTime
    df = dateformat"y/m/d,H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    EDat.RoundDT = gribCondGT.(EDat.DT, 3)
    ts = unique(EDat.RoundDT)
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    # gribLs = findGribPointsGT.(ts,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        try
            gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
            if length(gribD) > 0
                U = reduce(vcat,[x[1] for x in gribD])
                V = reduce(vcat,[x[2] for x in gribD])
                aveHeds[tPoint] = atan(mean(U),mean(V))
                aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
            else
                aveHeds[tPoint] = NaN
                aveSpds[tPoint] = NaN
            end
        catch 
            print(tPoint)
            print(b)
        end
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    # add = DataFrame(hcat(birdWindowGT.(ts),birdWSpeedGT.(ts)), :auto)
    # append!(toAdd,add)
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindowGT.(ts),birdWSpeedGT.(ts),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end
extra = outPt
## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
# rename!(EDat, [:DT,:Head,:X,:Y,:Lat,:Lon])
rename!(EDat, [:DT,:Lat,:Lon,:Forage,:DistTo,:BirdHead,:WDir,:WSpeed,:ID]) 
df = dateformat"d-m-y H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
Ll = LLA.(EDat.Lat,EDat.Lon)
UTMD = UTMZfromLLA(wgs84)
utmDat = map(UTMD, Ll); utmDat = StructArray(utmDat)
X = diff(utmDat.x)
Y = diff(utmDat.y)
head = atan.(X,Y)

CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMethodValidation3HrAve.csv", outPt)

gribD = gribSelect.(sel.DT, sel.Lat, sel.Lon)
##
sel.U = reduce(vcat,[x[1] for x in gribD])
sel.V = reduce(vcat,[x[2] for x in gribD])
sel.WindHead = atan.(sel.V,sel.U)
sel.EHead = atan.(sel.Y,sel.X)
# plot(sel.WindHead, sel.Head, seriestype = :scatter)
# plot(sqrt.(sel.X.^2 + sel.Y.^2), sqrt.(sel.U.^2 + sel.V.^2), seriestype = :scatter)
sel.ESpd = sqrt.(sel.X.^2 + sel.Y.^2)
sel.WSpd = sqrt.(sel.U.^2 + sel.V.^2)
CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv", outPt)
##
@rput outPt
R"""
library(circular)
library(ggplot2)
res<-cor.circular(outPt$estHead,outPt$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
#ggplot(outPt, aes(x = EHead, y = gribHead)) +
 #   geom_point() #+
    # geom_line(aes(x=-3:3,y=-3:3))
summary(res)
# spres <- circular::cor.test(sel$ESpd, sel$WSpd)
"""
## 
outPt
plot(sel.WindHead, sel.EHead, seriestype = scatter)


plot(sel.WindHead,sel.Head, seriestype = scatter)
plot(sel.ESpd,sel.WSpd, seriestype = scatter)

plot(outPt.estHead,outPt.gribHead,seriestype=scatter,legend=false,xlabel="Estimated wind direction ( rad)", ylabel="Validation wind direction (rad)")
plot!(-pi:pi,-pi:pi)

plot(outPt.gribSpeed,outPt.estSpeed,seriestype=scatter,legend=false,xlabel="Estimated wind speed (m/s)", ylabel="Validation wind speed (m/s)",fmt=:png)
plot!(0:10,0:10)


# attempt reading in .nc file from sar-winds


#############################################################################
############### RUN FOR YONEHARA METHOD 2014, 2016, AND 2017 ################
#############################################################################

#############################################################################
################################# ORIGINAL ##################################
#############################################################################

estLoc14 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/1sFix/"
estLoc16 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/1sFix/"
estLoc17 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/1sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles14 = readdir(estLoc14,join=true)
EFiles16 = readdir(estLoc16,join=true)
EFiles17 = readdir(estLoc17,join=true)
EFiles = [EFiles14; EFiles16; EFiles17]
Tags14 = readdir(estLoc14)
Tags16 = readdir(estLoc16)
Tags17 = readdir(estLoc17)
Tags = [Tags14; Tags16; Tags17]
EDat = DataFrame()
for b = 1:length(EFiles)
    append!(EDat, [CSV.File(EFiles[b], header = true) |> DataFrame fill(EFiles[b][40:43]*"-"*Tags[b][1:end-12],nrow(CSV.File(EFiles[b], header = true) |> DataFrame))])
end
rename!(EDat, [:DT,:Lon,:Lat,:WSpeed,:WDir,:aSpeed,:resnorm,:ID])
# remove non-wind calculated rows
EDat = EDat[.!isnan.(EDat.WSpeed),:]
# convert to DateTime
df = dateformat"dd-uuu-yyyy H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
# add gribs
sel = filter(:DT => x -> gribCond(x, 3, 25), EDat)
gribsels = unique(gribGen.(sel.DT))
# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
cd(WLoc)
@rput gribsels

gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
## DOWNLOAD THE REQUIRED GRIBFILES
R"""
# save(list=c("gribsels","gribTimes"), file = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
# load("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
options(timeout=800)
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
for(b in 1:length(gribsels)){
    if(any(gribFls == gribsels[b])){
        next
    } else {
        download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
            destfile = paste(dloadLoc,gribsels[b], sep = ""))
        }
}
"""

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/1sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = DataFrame()
outPt = DataFrame()
for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2014-"*EFiles[b][1:end-4],nrow(EDat))
    rename!(EDat, [:DT,:Lon,:Lat,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"d-uuu-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

# for b = 1:length(EFiles)
#     append!(EDat, [CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2014-"*EFiles[b][1:end-12],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
# end

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/1sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)

for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2016-"*EFiles[b][1:end-12],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"d-uuu-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

# RUN FOR 2017
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/1sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)

for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2017-"*EFiles[b][1:end-12],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"d-uuu-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/OG201420162017Validation.csv", outPt)

#############################################################################
################################ SUBSAMPLED #################################
#############################################################################

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2014Shearwater/WindEst/YoneMet/5sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = DataFrame()
outPt = DataFrame()
for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2014-"*EFiles[b][1:end-4],nrow(EDat))
    rename!(EDat, [:DT,:Lon,:Lat,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"d-uuu-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

# for b = 1:length(EFiles)
#     append!(EDat, [CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2014-"*EFiles[b][1:end-12],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
# end

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2016Shearwater/WindEst/YoneMet/5sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)

for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2016-"*EFiles[b][1:end-12],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"d-uuu-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

# RUN FOR 2017
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/5sFix/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)

EDat = DataFrame()
for b = 1:length(EFiles)
    append!(EDat, [CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2017-"*EFiles[b][1:end-12],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
end
rename!(EDat, [:DT,:Lat,:Lon,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
df = dateformat"y/m/d,H:M:S"
EDat.DT = DateTime.(EDat.DT, df)

for b = 1:length(EFiles)
    EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
    EDat.ID = fill("2017-"*EFiles[b][1:end-12],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
    # remove non-wind calculated rows
    EDat = EDat[.!isnan.(EDat.WSpeed),:]
    # convert to DateTime
    df = dateformat"y/m/d,H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 5), EDat)
    ts = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
    gribLs = findGribPoints.(ts,5,5000)
    aveHeds = zeros(length(ts),1)
    aveSpds = zeros(length(ts),1)
    for tPoint = 1:length(ts)
        gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
    rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    append!(outPt,add)
end

CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2017Shearwater/WindEst/YoneMet/Comb201420162017Validation.csv", outPt)
# for b = 1:length(EFiles)
#     append!(EDat, [CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame fill("2016-"*EFiles[b][1:end-12],nrow(CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame))])
# end
rename!(EDat, [:DT,:Lon,:Lat,:WSpeed,:WDir,:aSpeed,:resnorm,:ID]) 
# remove non-wind calculated rows
EDat = EDat[.!isnan.(EDat.WSpeed),:]
# convert to DateTime
df = dateformat"d-uuu-y H:M:S"
EDat.DT = DateTime.(EDat.DT, df)

Ll = LLA.(EDat.Lat,EDat.Lon)
UTMD = UTMZfromLLA(wgs84)
utmDat = map(UTMD, Ll); utmDat = StructArray(utmDat)
X = diff(utmDat.x)
Y = diff(utmDat.y)
head = atan.(X,Y)
sel = filter(:DT => x -> gribCond(x, 3, 5), EDat)

gribsels = unique(gribGen.(sel.DT))

# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
# sel.GLat = roundF.(sel.Lat, .05)
# sel.GLon = roundF.(sel.Lon, .0625)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
cd(WLoc)
@rput gribsels

gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
## DOWNLOAD THE REQUIRED GRIBFILES
R"""
options(timeout=800)
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
for(b in 1:length(gribsels)){
    if(any(gribFls == gribsels[b])){
        next
    } else {
        download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
            destfile = paste(dloadLoc,gribsels[b], sep = ""))
        }
}
"""

ts = DateTime.(unique(gribT.(sel.DT)),"yyyymmddHH")
cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
gribLs = findGribPoints.(ts,5,5000)
aveHeds = zeros(length(ts),1)
aveSpds = zeros(length(ts),1)
for tPoint = 1:length(ts)
    gribD = gribSelect.(fill(ts[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
    U = reduce(vcat,[x[1] for x in gribD])
    V = reduce(vcat,[x[2] for x in gribD])
    aveHeds[tPoint] = atan(mean(U),mean(V))
    aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
end
# switch non-calculated wind speeds and directions to NaN
EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
add = DataFrame(hcat(ts, aveHeds, aveSpds, birdWindow.(ts,5),birdWSpeed.(ts,5),fill(EFiles[b][1:(end-4)],length(ts))), :auto)
rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
append!(outPt,add)

extra = outPt
## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
# rename!(EDat, [:DT,:Head,:X,:Y,:Lat,:Lon])
rename!(EDat, [:DT,:Lat,:Lon,:Forage,:DistTo,:BirdHead,:WDir,:WSpeed,:ID]) 
df = dateformat"d-m-y H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
Ll = LLA.(EDat.Lat,EDat.Lon)
UTMD = UTMZfromLLA(wgs84)
utmDat = map(UTMD, Ll); utmDat = StructArray(utmDat)
X = diff(utmDat.x)
Y = diff(utmDat.y)
head = atan.(X,Y)
sel = filter(:DT => x -> gribCond(x, 3, 25), EDat)

gribsels = unique(gribGen.(sel.DT))

# grib times to nearest hour
gTimes = round.(sel.DT, Dates.Hour)
# sel.GLat = roundF.(sel.Lat, .05)
# sel.GLon = roundF.(sel.Lon, .0625)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
cd(WLoc)
@rput gribsels

gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
## DOWNLOAD THE REQUIRED GRIBFILES
R"""
# save(list=c("gribsels","gribTimes"), file = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
load("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/2019YoneGribs.RData")
options(timeout=800)
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/"
gribFls = list.files(dloadLoc,pattern="*grib2.bin")
for(b in 1:length(gribsels)){
        download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
            destfile = paste(dloadLoc,gribsels[b], sep = ""))
}
"""

# load in the wind calculations
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
cd(estLoc)

    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 25), EDat)
    t = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd(WLoc)
    gribLs = findGribPoints.(t)
    aveHeds = zeros(length(t),1)
    aveSpds = zeros(length(t),1)
    for tPoint = 1:length(t)
        gribD = gribSelect.(fill(t[tPoint],nrow(gribLs[tPoint])),gribLs[tPoint].lat,gribLs[tPoint].lon)
        U = reduce(vcat,[x[1] for x in gribD])
        V = reduce(vcat,[x[2] for x in gribD])
        aveHeds[tPoint] = atan(mean(U),mean(V))
        aveSpds[tPoint] = sqrt(mean(U)^2 + mean(V)^2)
    end
    if g == 1
        outPt = DataFrame(hcat(t, aveHeds, aveSpds, birdWindow.(t),birdWSpeed.(t),fill(EFilesAll[g][1:(end-4)],length(t))), :auto)
        rename!(outPt, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    else
        add = DataFrame(hcat(t, aveHeds, aveSpds, birdWindow.(t),birdWSpeed.(t),fill(EFilesAll[g][1:(end-4)],length(t))), :auto)
        rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
        outPt = vcat(outPt, add)
    end
outPt

gribD = gribSelect.(sel.DT, sel.Lat, sel.Lon)
##
sel.U = reduce(vcat,[x[1] for x in gribD])
sel.V = reduce(vcat,[x[2] for x in gribD])
sel.WindHead = atan.(sel.V,sel.U)
sel.EHead = atan.(sel.Y,sel.X)
# plot(sel.WindHead, sel.Head, seriestype = :scatter)
# plot(sqrt.(sel.X.^2 + sel.Y.^2), sqrt.(sel.U.^2 + sel.V.^2), seriestype = :scatter)
sel.ESpd = sqrt.(sel.X.^2 + sel.Y.^2)
sel.WSpd = sqrt.(sel.U.^2 + sel.V.^2)
CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/gribSelectedProper.csv", outPt)
##
@rput outPt
R"""
library(circular)
library(ggplot2)
res<-cor.circular(outPt$estHead,outPt$gribHead, test = T)
res
# res<-cor.circular(atan2(sel$X,sel$Y), atan2(sel$U,sel$V))
#ggplot(outPt, aes(x = EHead, y = gribHead)) +
 #   geom_point() #+
    # geom_line(aes(x=-3:3,y=-3:3))
res
# spres <- circular::cor.test(sel$ESpd, sel$WSpd)
"""
## 
outPt
plot(sel.WindHead, sel.EHead, seriestype = scatter)


plot(sel.WindHead,sel.Head, seriestype = scatter)
plot(sel.ESpd,sel.WSpd, seriestype = scatter)

plot(outPt.estHead,outPt.gribHead,seriestype=scatter,legend=false,xlabel="Estimated wind direction ( rad)", ylabel="Validation wind direction (rad)")
plot!(-3:3,-3:3)

plot(outPt.gribSpeed,outPt.estSpeed,seriestype=scatter,legend=false,xlabel="Estimated wind speed (m/s)", ylabel="Validation wind speed (m/s)")
plot!(0:10,0:10)


