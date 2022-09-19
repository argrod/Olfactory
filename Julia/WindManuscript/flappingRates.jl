# FLAPPING RATES FOR AXYTREK DATA
using DataFrames, CSV, RCall, Geodesy, Dates, Distances, Statistics, Glob, CategoricalArrays, DSP, StatsPlots, Peaks, StatsBase, KernelDensity
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
    dat.Timestamp = DateTime.(dat.Timestamp,dateformat"d/m/y H:M:S.s") .+ Hour(9)
    return dat
end
# READ IN GPS DATA
function readinGPS(yrID::String)
    file=rdir(dataloc,"(?=" * yrID[1:4] * "Shearwater).*AxyTrek.*(?=" * yrID[6:end] * "_...csv)")
    dat = CSV.read(file[1],DataFrame,header=1,select=[:TagID,:Timestamp,Symbol("location-lat"),Symbol("location-lon")])
    rename!(dat,[:TagID,:Timestamp,:lat,:lon])
    # remove missing rows
    dat = dropmissing(dat,:lat)
    dat.Timestamp = DateTime.(dat.Timestamp,dateformat"d/m/y H:M:S.s") .+ Hour(9)
    return dat
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
# find maximum frequencies above limit Hz
function maxFreqs(outLocation,yrID,fs)
    dat = readinAxy(yrID) # read in acceleration
    GPSdat = readinGPS(yrID) # read in GPS
    rename!(GPSdat,[:TagID,:Timestamp,:lat,:lon])
    GPSInds = findall(ismissing.(GPSdat.lat) .== false) # take indeces of GPS positions
    # select indeces from 2 seconds to end time - 2 seconds
    GPSInds = GPSInds[(fs*2) .<= GPSInds .<= (nrow(GPSdat) - (fs*2))]
    # lowpass filter to separate dynamic and static acceleration
    Z = dynstat.([float.(dat.X),float.(dat.Y),float.(dat.Z)],Ref(1.0),Ref(1.5),Ref(fs))[3]
    winsize = Int(fs*60)
    numoverlap=round(.5*winsize)
    win=Windows.hamming(winsize)
    spect = spectrogram(Z[2],winsize,Int(numoverlap),fs=fs,window=win) # spectrogram of dynamic dorsoventral acceleration
    maxix = findmax(spect.power[3 .<= spect.freq .<= 6,:];dims=1)[2] # find indeces of maximum frequencies
    maxix = getindex.(maxix,1) .+ findfirst(spect.freq .>= 3) .- 1
    maxFrec = [spect.freq[x] for x in maxix] # take the frequencies
    # take a range of acceleration indeces around the GPS positions
    mxFrcRanges = range.([maximum([Int(round(x/(30*fs))-(30/spect.time[1])),1]) for x in GPSInds],[minimum([Int(round(x/(30*fs))+(30/spect.time[1])),length(spect.time)]) for x in GPSInds],step=1);
    # take mean max frequencies for each GPS position
    GPSfrec = [mean(maxFrec[x]) for x in mxFrcRanges];
    # create an output
    namelist = [:yrID,:DT,:lat,:lon,:domFreq]
    df = DataFrame([repeat([yrID],length(GPSInds)),GPSdat[GPSInds,"Timestamp"],GPSdat[GPSInds,"lat"],GPSdat[GPSInds,"lon"],GPSfrec], namelist);
    CSV.write(outLocation * yrID * "domFreq.csv",df);
end
# find peaks and troughs of signal (all)
function pkstrghs(signal)
    peaks = Int[]
    if length(signal)>1
        for i=2:length(signal)-1
            if signal[i-1]<=signal[i]>signal[i+1]
                push!(peaks,i)
            end
        end
    end
    troughs = Int[]
    if length(signal)>1
        for i=2:length(signal)-1
            if signal[i-1]>signal[i]<=signal[i+1]
                push!(troughs,i)
            end
        end
    end
    # ensure starts with peaks and ends with trough
    if troughs[1] < peaks[1]
        troughs = troughs[2:end]
    end
    if peaks[end] > troughs[end]
        peaks = peaks[1:end-1]
    end
    peaks,troughs
end


# flapping vs gliding rates (minute average)
function flapglide(outLocation,yrID,fs)
    dat = readinAxy(yrID) # read in acceleration
    # lowpass filter to separate dynamic and static acceleration
    Z = dynstat(float.(dat.Z),1.0,1.5,fs)
    pks,trghs = pkstrghs(Z)
    den = kde(Z[pks] .- Z[trghs])
    plot(den)
    sep = den.x[findmin(den.density[findmax(den.density)[2]:findmax(den.density[.5 .<= den.x .<= 1])[2] + minimum(findall(den.x .>= .5))])[2]]
    flapMotion = sort([pks[findall((Z[pks] .- Z[trghs]) .> sep)];trghs[findall((Z[pks] .- Z[trghs]) .> sep)]])
    # group flaps by typically flapping frequency
    flaps = zeros(Int,length(Z))
    gaps = findall(diff(flapMotion) .> round(.5*fs))
    for b = eachindex(gaps)
        if b == 1
            flaps[flapMotion[1]:flapMotion[gaps[b]]] .= b
        elseif b == flapMotion[end]
            flaps[flapMotion[(gaps[b]+1):end]] .= b
        else
            flaps[flapMotion[(gaps[b-1]+1):gaps[b]]] .= b
        end
    end

    flSt = findall(diff(flaps) .> 0) .+ 1
    flEd = findall(diff(flaps) .< 0) .+ 1
    flBdur = (flEd .- flSt) / fs;

    GPSmdur = [mean(flBdur[(GPSdat.Timestamp[x] - Second(30)) .< dat.Timestamp[flSt] .> (GPSdat.Timestamp[x] + Second(30))]) for x in 1:nrow(GPSdat)]
    GPSstdur = [std(flBdur[(GPSdat.Timestamp[x] - Second(30)) .< dat.Timestamp[flSt] .> (GPSdat.Timestamp[x] + Second(30))]) for x in 1:nrow(GPSdat)]
    
    CSV.write(outLocation * yrID * "flapRatio.csv",DataFrame(DT=GPSdat.Timestamp, meanFlapRatio = GPSmdur, sdFlapRatio = GPSstdur)); 
end

# file locations for raw acceleration data
if Sys.iswindows()
    dataloc = "E:/My Drive/PhD/Data/"
else
    dataloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/"
end
# output locations for GPS dominant frequencies (dorsoventral)
if Sys.iswindows()
    outloc = "E:/My Drive/PhD/Data/GPSDominantFrequencies/"
else
    outloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/GPSDominantFrequencies/"
end
# find year and IDs of all wind data
yrIDs = yrIDGather(dataloc,".*MinDat[\\\\|/](.*)_S.*",[2018,2019])
fs = 25
for b in yrIDs
    maxFreqs(outloc,b,fs)
end

# output the flapping duration data
if Sys.iswindows()
    fgoutloc = "I:/My Drive/PD/Data/flappingDurations/"
else
    fgoutloc = "/Volumes/GoogleDrive-112399531131798335686/My Drive/PhD/Data/flappingDurations/"
end

[flapglide(fgoutloc,x,25) for x in yrIDs]

length(yrIDs)


yrid = yrIDs[14]
dat = readinAxy(yrid) # read in acceleration
GPSdat = readinGPS(yrid) # read in GPS
rename!(GPSdat,[:TagID,:Timestamp,:lat,:lon])
GPSInds = findall(ismissing.(GPSdat.lat) .== false) # take indeces of GPS positions
# select indeces from 2 seconds to end time - 2 seconds
GPSInds = GPSInds[(fs*2) .<= GPSInds .<= (nrow(GPSdat) - (fs*2))]
# lowpass filter to separate dynamic and static acceleration
Z = dynstat.([float.(dat.X),float.(dat.Y),float.(dat.Z)],Ref(1.0),Ref(1.5),Ref(fs))[3]
winsize = Int(fs*4)
numoverlap=round(.9874*winsize)
win=Windows.hamming(winsize)
spect = spectrogram(Z[2],winsize,Int(numoverlap),fs=fs,window=win) # spectrogram of dynamic dorsoventral acceleration
maxix = findmax(spect.power;dims=1)[2] # find indeces of maximum frequencies
maxix = getindex.(maxix,1)
maxFrec = [spect.freq[x] for x in maxix] # take the frequencies
GPSfs = Second(mode(diff(GPSdat.Timestamp[GPSInds]))).value # GPS sampling frequency (in seconds)
# take a range of acceleration indeces around the GPS positions
lower = Int.((GPSInds .- (fs*round(GPSfs/2))) .+ (2*fs))
lower[lower .< 1] .= 1
upper = Int.((GPSInds .+ (fs*round(GPSfs/2))) .+ fs*2)
upper[upper .> length(maxFrec)] .= length(maxFrec)
mxFrcRanges = range.(lower,upper);
# take mean max frequencies for each GPS position
GPSfrec = [mean(maxFrec[x]) for x in mxFrcRanges];
# create an output
namelist = [:yrID,:DT,:lat,:lon,:domFreq]
df = DataFrame([repeat([yrid],length(GPSInds)),GPSdat[GPSInds,"Timestamp"],GPSdat[GPSInds,"lat"],GPSdat[GPSInds,"lon"],GPSfrec], namelist);
CSV.write(outLocation * yrid * "domFreq.csv",df);

plot(GPSfrec,seriestype=:scatter)

GPSdat.Timestamp[GPSInds[50] + (fs*2)]
specTime[GPSInds[50]]

GPSdat.Timestamp[Int(GPSInds[50] + (fs*2) - (fs*round(GPSfs/2)))]
specTime[GPSInds[50]] - Second(round(GPSfs/2))

spect.time ./ .04

spect.time[end]

nrow(GPSdat) - fs*2

Second(dat.Timestamp[end] - dat.Timestamp[1])

maxFrec[(x - (round(GPSfs/2)*fs))]
spect.freq[maxix[1][2]]


maxix[14]

maxFr = maxFreqs(dat.Z,fs)

# MATLAB process
winsize = Int(fs*4)
numoverlap=round(.9874*winsize)
win=Windows.hamming(winsize)
spectDynZ = spectrogram(Z[2],winsize,Int(numoverlap),fs=fs,window=win)

# maximum frequencies of 3Hz+
threeHzmark = findfirst(spectDynZ.freq .== 3)
power3Up = spectDynZ.power[threeHzmark:end,:]
maxval,maxix = findmax(spectDynZ.power[threeHzmark:end,:];dims=1)
maxFreqs = [spectDynZ.freq[x] for x in getindex.(maxix,1) .+ threeHzmark]

plot(maxFreqs[1:100])

size(spectDynZ.power[13:end,:])

spectDynZ.power[threeHzmark:end,1]

range(13,52)

range(threeHzmark,size(spectDynZ.power,1))

