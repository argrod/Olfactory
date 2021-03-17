using GRIB
using DataFrames
using Query
using Dates
using RCall
using CSV
using Geodesy
using Plots
using StructArrays
using RecursiveArrayTools
##

# load in the wind calculations
estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
# estLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/2018/"
cd(estLoc)
EFiles = readdir(estLoc)
EDat = []
for b = 1:length(EFiles)
    if b == 1
        EDat = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
    else Add = CSV.File(string(estLoc,EFiles[b]), header = false) |> DataFrame
        EDat = vcat(EDat, Add)
    end
end

## FORMAT DATA INTO UTM, FIND XY AND WIND HEADING
rename!(EDat, [:DT,:Lat,:Lon,:Head,:X,:Y])
df = dateformat"y-m-d H:M:S"
EDat.DT = DateTime.(EDat.DT, df)
Ll = LLA.(EDat.Lat,EDat.Lon)
UTMD = UTMZfromLLA(wgs84)
utmDat = map(UTMD, Ll); utmDat = StructArray(utmDat)
X = diff(utmDat.x)
Y = diff(utmDat.y)
head = atan.(X,Y)
## FUNCTION TO FIND THE INDECES OF SUITABLE WIND ESTIMATIONS
function gribCond(dt, h, md)
    hr = Dates.hour(round(dt, Dates.Hour))
    min = Dates.minute(dt)
    iszero(hr%h) && (min > 60 - md || min < md)
end
sel = filter(:DT => x -> gribCond(x, 3, 4), EDat)
##
function gribGen(dt)
    yr = string(Dates.year(dt))
    mn = string(Dates.month(dt))
    if length(mn) == 1
        mn = "0"*mn
    end
    dy = string(Dates.day(dt))
    if length(dy) == 1
        dy = "0"*dy
    end
    hr = string(Dates.hour(round(dt, Dates.Hour)))
    if length(hr) == 1
        hr = "0"*hr
    end
    "Z__C_RJTD_"*yr*mn*dy*hr*"0000_MSM_GPV_Rjp_Lsurf_FH00-15_grib2.bin"
end

gribsels = unique(gribGen.(sel.DT))
##
function roundF(vals, rnum)
    rnum*(round(vals)/rnum)
end
# sel.GLat = roundF.(sel.Lat, .05)
# sel.GLon = roundF.(sel.Lon, .0625)
##
WLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/gribs/"
cd(WLoc)
@rput gribsels
function gribT(dt)
    yr = string(Dates.year(dt))
    mn = string(Dates.month(dt))
    if length(mn) == 1
        mn = "0"*mn
    end
    dy = string(Dates.day(dt))
    if length(dy) == 1
        dy = "0"*dy
    end
    hr = string(Dates.hour(round(dt, Dates.Hour)))
    if length(hr) == 1
        hr = "0"*hr
    end
    return yr*mn*dy*hr
end
gribTimes = unique(gribT.(sel.DT))
@rput gribTimes
##
R"""
dloadLoc  = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/WindEst/MinDat/"
for(b in 1:length(gribsels)){
	download.file(paste("http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original/",substr(gribTimes[b],1,4),"/",substr(gribTimes[b],5,6),"/",substr(gribTimes[b],7,8),"/",gribsels[b],sep=""), 
		destfile = paste(dloadLoc,gribsels[b], sep = "")
		)
}
"""

##
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
##
function xtractUWind00(GribN)
    GribFile(GribN) do f
        for message in f
            message["forecastTime"] == 0 && message["shortName"] == "10u" && return message
        end
    end
end
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
    return UWind.value[@. (UWind.lon .== nrlon) & (UWind.lat .== nrlat)]..., VWind.value[@. (UWind.lon .== nrlon) & (UWind.lat .== nrlat)]...
end
##
gribD = gribSelect.(sel.DT, sel.Lat, sel.Lon)
sel.U = reduce(vcat,[x[1] for x in gribD])
sel.V = reduce(vcat,[x[2] for x in gribD])
sel.WindHead = atan.(sel.V,sel.U)
# plot(sel.WindHead, sel.Head, seriestype = :scatter)
# plot(sqrt.(sel.X.^2 + sel.Y.^2), sqrt.(sel.U.^2 + sel.V.^2), seriestype = :scatter)
sel.ESpd = sqrt.(sel.X.^2 + sel.Y.^2)
sel.WSpd = sqrt.(sel.U.^2 + sel.V.^2)
# CSV.write("/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/gribSelected.csv", sel)
##
@rput sel
R"""
library(circular)
res<-cor.circular(sel$WindHead, sel$Head, test = T)
# spres <- circular::cor.test(sel$ESpd, sel$WSpd)
"""
## 
plot(sel.WindHead,sel.Head, seriestype = scatter)
plot(sel.ESpd,sel.WSpd, seriestype = scatter)