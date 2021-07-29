using Pkg
Pkg.add("GRIB")
Pkg.add("DataFrames")
Pkg.add("Query")
Pkg.add("CSV")
Pkg.add("Geodesy")
Pkg.add("Plots")
Pkg.add("StructArrays")
Pkg.add("RecursiveArrayTools")
Pkg.add("Statistics")
Pkg.add("RCall")
Pkg.add("NetCDF")
using GRIB
using DataFrames
using Query
using Dates
using CSV
using Geodesy
using Plots
using StructArrays
using RecursiveArrayTools
using Statistics
using RCall
using NetCDF

## FUNCTIONS
## FUNCTION TO FIND THE INDECES OF SUITABLE WIND ESTIMATIONS
function gribCond(dt, h, md)
    hr = Dates.hour(round(dt, Dates.Hour))
    min = Dates.minute(dt)
    iszero(hr%h) && (min > 60 - md || min < md)
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
## FIND LAT AND LON WITHIN 5 KM OF TRACK WITHIN 50 MINS OF t
function findGribPoints(t)
    # isolate points along track to find lat lon ±50mins
    isol = EDat[(t - Minute(50)) .< EDat.DT .< (t + Minute(50)),:]
    lon = isol.Lon
    lat = isol.Lat
    @rput lon lat
    R"""
    library(sp)
    library(raster)
    library(plyr)
    pnts <- SpatialPoints(cbind(lon,lat),proj4string=CRS("+proj=longlat"))
    pntsBuff <- buffer(pnts, 5000)
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
## FIND AVERAGE WIND HEADING WITHIN 5KM FOR GIVEN DATETIME LAT LON
function gribAve(tPoint)
    iso = EDat[(tPoint - Minute(50)) .< EDat.DT .< (tPoint + Minute(50)),:]
    gribPts = findGribPoints.(iso.Lon, iso.Lat)
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
## CALCULATE AVERAGE ESTIMATED WIND WITHIN 25 OF TIME tPoints
function birdWindow(tPoints)
    iso = EDat[(tPoints - Minute(25)) .< EDat.DT .< (tPoints + Minute(25)),:]
    isol = filter(:DT => x -> gribCond(x, 3, 25), iso)
    aveWHd = aveWind(isol.X,isol.Y)
    # iso = EDat[(tPoints[b] - Minute(50)) .< EDat.DT .< (tPoints[b] + Minute(50)),:]
    # isol = filter(:DT => x -> gribCond(x, 3, 50), iso)
    # gribAve.(isol.DT, isol.Lon,isol.Lat)
    return aveWHd
end
##
function birdWSpeed(tPoints)
    iso = EDat[(tPoints - Minute(25)) .< EDat.DT .< (tPoints + Minute(25)),:]
    isol = filter(:DT => x -> gribCond(x, 3, 25), iso)
    aveSpd = sqrt(mean(isol.X)^2 + mean(isol.Y)^2)
    return aveSpd
end
# ## 
# function WindVal(t, lat, lon)
#     # calculate mean estimated heading for ±25 min
#     estHdAve = birdWindow(t)
#     # find all bird locations within ±50 mins and all grib locations within 5km
#     gribloc = findGribPoints(EDat, t)
# end

# load in the wind calculations
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = []
for b = 1:length(EFiles)
    if b == 1
        EDat = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
        EDat.ID = fill("18"*EFiles[b][1:end-4],nrow(EDat))
    else Add = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
        Add.ID = fill("18"*EFiles[b][1:end-4],nrow(Add))
        EDat = vcat(EDat, Add)
    end
end
# repeat for 2019
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EFiles = EFiles[contains.(EFiles,".csv")]
EDat19 = []
for b = 1:length(EFiles)
    if b == 1
        EDat19 = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
        EDat19.ID = fill("19"*EFiles[b][1:end-4],nrow(EDat19))
    else Add19 = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
        Add19.ID = fill("19"*EFiles[b][1:end-4],nrow(Add19))
        EDat19 = vcat(EDat19, Add19)
    end
end
EDat = vcat(EDat,EDat19)

## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
# rename!(EDat, [:DT,:Head,:X,:Y,:Lat,:Lon])
rename!(EDat, [:DT,:Lat,:Lon,:Head,:X,:Y,:ID]) 
df = dateformat"y-m-d H:M:S"
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
EFiles = readdir(estLoc)
# repeat for 2019
estLoc19 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc19)
EFiles19 = readdir(estLoc19)
EFiles19 = EFiles19[contains.(EFiles19,".csv")]
EFilesAll = vcat(EFiles,EFiles19)

for g = 1:length(EFilesAll)
    if g > length(EFiles)
        cd(estLoc19)
        # read in and format data
        EDat = CSV.File(string(estLoc19,EFilesAll[g]), header = false) |> DataFrame
    else
        cd(estLoc)
        # read in and format data
        EDat = CSV.File(string(estLoc,EFilesAll[g]), header = false) |> DataFrame
    end
    EDat.ID = fill("18"*EFilesAll[g][1:end-4],nrow(EDat))
    rename!(EDat, [:DT,:Lat,:Lon,:Head,:X,:Y,:ID]) 
    df = dateformat"y-m-d H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
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
end

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

##############################################################################
########################## RUN FOR YONEHARA METHOD ###########################
##############################################################################

estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/YoneMet/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = DataFrame()
for b = 1:length(EFiles)
    if b == 1
        EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
        EDat.ID = fill("2019-"*EFiles[b][1:end-4],nrow(EDat))
        rename!(EDat, [:DT,:Lat,:Lon,:Forage,:DistTo,:BirdHead,:WDir,:WSpeed,:ID]) 
    else Add = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
        EDat = CSV.File(string(estLoc,EFiles[b]), header = true) |> DataFrame
        EDat.ID = fill("2019-"*EFiles[b][1:end-4],nrow(EDat))
        rename!(EDat, [:DT,:Lat,:Lon,:Forage,:DistTo,:BirdHead,:WDir,:WSpeed,:ID]) 
    end
    # remove non-wind calculated rows
    EDat = EDat[.!(EDat.WSpeed .== 0),:]
    # convert to DateTime
    df = dateformat"d-m-y H:M:S"
    EDat.DT = DateTime.(EDat.DT, df)
    # select time points
    compT = filter(:DT => x -> gribCond(x, 3, 25), EDat)
    t = DateTime.(unique(gribT.(compT.DT)),"yyyymmddHH")
    cd("/Volumes/GoogleDrive/My Drive/PhD/Data/gribs/")
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
    # switch non-calculated wind speeds and directions to NaN
    EDat.X = EDat.WSpeed.*(cos.(EDat.WDir))
    EDat.Y = EDat.WSpeed.*(sin.(EDat.WDir))
    if b == 1
        outPt = DataFrame(hcat(t, aveHeds, aveSpds, birdWindow.(t),birdWSpeed.(t),fill(EFilesAll[b][1:(end-4)],length(t))), :auto)
        rename!(outPt, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
    else
        add = DataFrame(hcat(t, aveHeds, aveSpds, birdWindow.(t),birdWSpeed.(t),fill(EFilesAll[b][1:(end-4)],length(t))), :auto)
        rename!(add, [:Time,:gribHead,:gribSpeed,:estHead,:estSpeed,:tagID])
        outPt = vcat(outPt, add)
    end
end

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

# attempt reading in .nc file from sar-winds
