using CSV, DataFrames, PyPlot, Geodesy, Glob, Dates

# foraging locations
if Sys.iswindows()
    fileloc = "E:/My Drive/PhD/Data/2018Shearwater/TxtDat/AxyTrek/AlgorithmOutputUpdate/PredictedForage/"
else
    fileloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/TxtDat/AxyTrek/AlgorithmOutputUpdate/PredictedForage/"
end

files = glob("*ForageGPS.txt",fileloc)
forDat = DataFrame()
for b=1:length(files)
    append!(forDat,CSV.File(files[b],header=true,delim=",") |> DataFrame)
end
dt = dateformat"d/m/y H:M:S.s"
forDat.DT = DateTime.(forDat.DT,dt)

# wind data
if Sys.iswindows()
    windloc = "E:/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
else
    windloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
end
wfiles = glob("*.csv",windloc)
wDat = DataFrame()
for b=1:length(wfiles)
    append!(wDat,CSV.File(wfiles[b],header=false,delim=",") |> DataFrame)
end
rename!(wDat,[:DT,:Lat,:Lon,:BHead,:X,:Y])

clf()
q=quiver(wDat.Lon,wDat.Lat,wDat.X,wDat.Y)
ax = gca()
ax.quiverkey(q,X=141,Y = 40, U = 10,coordinates="figure", label="Quiver key, length = 10",labelpos = "E")
PyPlot.title("Quiver Plot Example")
gcf()


scatter(wDat.Lon,wDat.Lat)
gcf()