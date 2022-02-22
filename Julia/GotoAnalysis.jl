using CSV, DataFrames, PyPlot, Geodesy, Glob, Dates

# foraging locations
if Sys.iswindows()
    fileloc = "E:/My Drive/PhD/Data/2018Shearwater/TxtDat/AxyTrek/AlgorithmOutputUpdate/PredictedForage/"
else
    fileloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2018Shearwater/TxtDat/AxyTrek/AlgorithmOutputUpdate/PredictedForage/"
end

files = glob("*ForageGPS.txt",fileloc)
forDat = DataFrame()
for b=1:length(files)
    append!(forDat,CSV.File(files[b],header=true,delim=",") |> DataFrame)
end
dt = dateformat"d/m/y H:M:S.s"
forDat.DT = DateTime.(forDat.DT,dt)

occursin(r".*/[0-9]_S.*.csv",wfiles[1])

# wind data
if Sys.iswindows()
    windloc = "E:/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
else
    windloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/2018Shearwater/WindEst/MinDat/"
end
wfiles = glob("*.csv",windloc)
wDat = DataFrame()
for b=1:length(wfiles)
    append!(wDat,CSV.File(wfiles[b],header=false,delim=",") |> DataFrame)
end
rename!(wDat,[:DT,:Lat,:Lon,:BHead,:X,:Y])
wDat.wHd = atan.(wDat.Y,wDat.X)

clf()
q=quiver(wDat.Lon,wDat.Lat,wDat.X,wDat.Y)
ax = gca()
ax.quiverkey(q,X=141,Y = 40, U = 10,coordinates="figure", label="Quiver key, length = 10",labelpos = "E")
PyPlot.title("Quiver Plot Example")
gcf()

function WatsonTwoTestRad(x,y)
    a = [mod2pi(x) repeat(1,length(x))]
    sort!()
end

function watson.two.test(x,y,alpha=0)
end

mod2pi(-2.4)


names(wDat)