# FLAPPING RATES FOR AXYTREK DATA

using DataFrames, CSV, RCall, Geodesy, Dates, Distances, Statistics, Glob, CategoricalArrays, DSP, StatsPlots, Peaks
using Plots; theme(:dark)

# RECURSIVE FILE SEARCH
function rdir(dir::AbstractString, pat::AbstractString)
    result = String[]
    for (root, dirs, files) in walkdir(dir)
        append!(result, filter!(f -> occursin(Regex(pat), f), joinpath.(root, files)))
    end
    return result[occursin.(Regex(pat), result)]
end
# EXTRACT YEARIDs OF WIND DATA
function yrIDGather(dataLocation::String, IDpattern::String, year::Vector{Int64})

    #=
    dataLocation:           Placement search for recursive file search
    IDPattern:              Pattern to extract year ID patterns
    year:                   Select year for data search
    =#

    if length(year) == 1
        year = [year]
    end
    # find wind files
    wFiles = rdir.(dataLocation,"(?=".*join(string.(year), "Shearwater|").*"Shearwater).*(?=MinDat).*.csv")
    
    # isolate the unique datasets (year_tagID)
    yrIDs = unique(getindex.(match.(r"(\d+)Shearwater.*",wFiles),1) .* "_" .* getindex.(match.(Regex(IDpattern),wFiles),1))

    return yrIDs
end
# READ IN RAW AXYTREK ACCELERATION
function readinAxy(yrID::String)
    file=rdir(dataloc,"(?=" * yrID[1:4] * "Shearwater).*AxyTrek.*(?=" * yrID[6:end] * "_...csv)")
    dat = CSV.read(file[1],DataFrame,header=1,select=[:TagID,:Timestamp,:X,:Y,:Z]  
    )
    # remove missing rows
    dat = dropmissing(dat,:Timestamp)
    tst.Timestamp = DateTime.(tst.Timestamp,dateformat"d/m/y H:M:S.s") .+ Hour(9)
    return dat
end

# file locations for raw acceleration data
if Sys.iswindows()
    dataloc = "E:/My Drive/PhD/Data/"
else
    dataloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/"
end

# LOWPASS EQUIRIPPLE FILTER
function lowP(x,stop,pass,fs,pole)
    fil = remez(pole+1,[(0, stop) => 1, (pass, fs/2) => 0]; Hz = fs)
    return filtfilt(fil,x)
end
# extract dynamic and static acceleration
function dynstat(x,stop,pass,fs)
    statX = lowP(x,stop,pass,fs,100)
    dynX = x .- statX
    return statX, dynX
end

# find year and IDs of all wind data
yrIDs = yrIDGather(dataloc,".*MinDat[\\\\|/](.*)_S.*",[2018,2019])

yrid=yrIDs[1]
tst = readinAxy(yrid)

fs = 25 # sampling rate
# lowpass filter to separate dynamic and static acceleration
X,Y,Z = dynstat.([float.(tst.X),float.(tst.Y),float.(tst.Z)],Ref(1.0),Ref(1.5),Ref(fs))

# MATLAB process
winsize = Int(fs*4)
numoverlap=round(.9874*winsize)
win=Windows.hamming(winsize)
spectDynZ = spectrogram(Z[2],winsize,Int(numoverlap),fs=fs,window=win)

# plot(perDynZ.Freq)
# plot(perDynZ.freq,pow2db.(perDynZ.power))
time(spectDynZ)