using GRIB
using DataFrames
using Query
using Dates
using RCall
DatSum = DataFrame(ToL = String[], sNam = String[], nam = String[], levels = Int[], units = String[], dataTime = Int[])

##
msg = []
GribFile("flob.bin") do f
    for message in f
        msg = Message(f)
    end
end

##

GribFile("flob.bin") do f
    for message in f
        push!(DatSum, (message["typeOfLevel"],message["shortName"],message["name"],message["level"], message["units"], message["forecastTime"],))
    end
    # [(message["typeOfLevel"],message["shortName"],message["name"]) for message in f]
end
unique(DatSum[!,6])
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


##

uTest = @NamedTuple{lon, lat, value}(data(xtractUWind00("flob.bin")))
msg = xtractUWind00("flob.bin")
##

keylist = Vector{String}()
for key in keys(msg)
    push!(keylist, key)
end

##

UD = []
VD = []
GribFile("flob.bin") do f


    for message in f
        if message["name"] == "10 metre U wind component" && message["typeOfLevel"] == "heightAboveGround" && message["level"] == 10
            lon, lat, value = data(message)
            push!(UD, (; lon, lat, value))
        end
        if message["name"] == "10 metre V wind component" && message["typeOfLevel"] == "heightAboveGround" && message["level"] == 10
            lon, lat, value = data(message)
            push!(VD, (; lon, lat, value))
        end
    end
end

##



##

# sel = Query.@from i in DatSum begin
#     @where i.nam = "10 metre U wind component" | i.nam = "10 metre V wind component"
#     @collect DataFrame
# end
sel = DatSum |>
    @filter(_.nam == "10 metre U wind component") |>
    DataFrame



##

Index("flob.bin", "level") do index
    GRIB.select!(index, "level", 10)
    # GRIB.select!(index, "name", "10 metre U wind component")
    # GRIB.select!(index, "typeOfLevel", "heightAboveGround")
    for msg in index
        display(msg)
    end
end
# display(nms)

##

# j = [0.0,0,0,0,0]
# for i = 1:5
#     j = rand()
# end
##

# R VALIDATION SET, BRING IN DATA

##


if Sys.iswindows()
    windLoc = 'F:/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/MinDat/'
	estLoc = 'F:/UTokyoDrive/PhD/Data/WindEstTest/Comparison/'
	dloadLoc = 'F:/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/WindValidate/'
    gribLoc = "F:/Documents/UTokyoDrive/PhD/Data/2018Shearwater/WindEst/WindValidate/"
else
    windLoc = '/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/'
	estLoc = '/Volumes/GoogleDrive/My Drive/PhD/Data/WindEstTest/Comparison/'
	dloadLoc = '/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/'
    gribLoc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/WindValidate/"
end
gribFs = readdir(gribLoc)



##
# rcopy(b)
R"""
library(circular)
"""
rcopy(R"b")
