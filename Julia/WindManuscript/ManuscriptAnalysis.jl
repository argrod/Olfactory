using DataFrames, CSV, RCall, Geodesy, Dates, Statistics, Glob, CategoricalArrays, DSP, StatsPlots, Clustering
using Plots; theme(:dark)

# RECURSIVE FILE SEARCH
function rdir(dir::AbstractString, pat::AbstractString)
    result = String[]
    for (root, dirs, files) in walkdir(dir)
        append!(result, filter!(f -> occursin(Regex(pat), f), joinpath.(root, files)))
    end
    return result[occursin.(Regex(pat), result)]
end

# DATETIME SETTING FUNCTION FOR AVAILABLE DATEFORMATS
function robust_DateTime_parse(str)
    if any(occursin.(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],str))
        DateTime(str,dateformat"d-u-y H:M:S.s")
    elseif any(isequal('-'), str)
        DateTime(str, dateformat"y-m-d H:M:S")
    else
        DateTime(str, dateformat"d/m/y H:M:S.s")
    end
end

# FUNCTION FOR READING IN FORAGING AND WIND DATASETS
function readDat(dataLocation::String, IDpattern::String, AndWind::Bool, colnames::Vector{Symbol}, header::Vector{Int64}, year::Vector{Int64})

    #=
    dataLocation:           Placement search for recursive file searcg
    IDPattern:              Pattern to extract year ID patterns
    colnames:               Array of column name symbols (two arrays required if wind requested
    header:                 Number of header rows (two values required if wind requested)
    year:                   Select year for data search
    =#

    if length(year) == 1
        year = [year]
    end
    # find forageGPS files
    forFiles = rdir(dataLocation, "(?=" * join(string.(year), "Shearwater|") * "Shearwater)"* ".*(?=BehaviourDetection)" * ".*(?=PredictedForage)" * ".*(?=ForageGPS.txt)")
    # isolate the unique datasets (year_tagID)
    yrIDs = unique(getindex.(match.(r"(\d+)Shearwater.*",forFiles),1) .* "_" .* getindex.(match.(Regex(IDpattern),forFiles),1))

    # preallocate output
    ForRet = [DataFrame() for _ in 1:length(yrIDs)]
    # repeat for wind if desired
    if AndWind == true
        wFiles = rdir.(dataLocation,"(?=".*join(string.(year), "Shearwater|").*"Shearwater).*(?=MinDat).*.csv")
        WindRet = [DataFrame() for _ in 1:length(yrIDs)]
    end
    # read in data
    for tg in 1:length(ForRet)
        fFiles = forFiles[occursin.(yrIDs[tg][1:4] * "Shearwater",forFiles) .& occursin.(Regex("[\\\\|/]" * "(?=" * yrIDs[tg][6:end] * "_)"),forFiles)]
        for file in fFiles
            append!(ForRet[tg], hcat(CSV.read(file, DataFrame, header = header[1]), repeat([yrIDs[tg]], nrow(CSV.read(file, DataFrame, header = header[1])))), cols = :union)
        end
        rename!(ForRet[tg], colnames)
        if AndWind == true
            windFiles = wFiles[occursin.(yrIDs[tg][1:4] * "Shearwater",wFiles) .& occursin.(Regex("[\\\\|/]" * "(?=" * yrIDs[tg][6:end] * "_)"),wFiles)]
            for file in windFiles
                append!(WindRet[tg], hcat(CSV.read(file, DataFrame, header = header[2]), repeat([yrIDs[tg]], nrow(CSV.read(file, DataFrame, header = header[2])))), cols = :union)
            end
            rename!(WindRet[tg], [:DT,:lat,:lon,:head,:X,:Y,:yrID])
        end
        # assign datetime
        ForRet[tg].DT = robust_DateTime_parse.(ForRet[tg].DT)
        if AndWind == true
            WindRet[tg].DT = robust_DateTime_parse.(WindRet[tg].DT)
        end
    end
    return ForRet, WindRet
end

# add distance (m) and speed (kph) values
function dist(lat,lon)
    ll = LLA.(lat,lon)
    utmz = UTMZfromLLA(wgs84)
    return sqrt.(diff((p->p.x).(utmz.(ll))).^2 + diff((p->p.y).(utmz.(ll))).^2)
end
function distToOne(LL1,lats,lons)
    ll1 = LLA(LL1[1],LL1[2])
    ll2 = LLA.(lats,lons)
    utmz = UTMZfromLLA(wgs84)
    return sqrt.((utmz(ll1).x .- (p->p.x).(utmz.(ll2))).^2 .+ (utmz(ll1).y .- (p->p.y).(utmz.(ll2))).^2)
end

function speed(dt,lat,lon)
    tdiff = Dates.value.(Second.(diff(dt)))
    spTrav = (dist(lat,lon)./tdiff).*3.6
    return spTrav
end

# find the nearest time (index)
function findNearest(dt, time)
    argmin(abs.(dt .- time))
end

# calculate linearity
function linearity(dt,lat,lon,distance,twindow)
    out = fill(NaN, length(dt))
    trange = dt[dt .<= (dt[end] - Second(twindow*60/2))]
    pos = [findNearest(dt,x - (Second(60*twindow/2))):findNearest(dt,x+(Second(60*twindow/2))) for x in trange]
    dists = [dist([lat[pos[b][1]],lat[pos[b][end]]],[lon[pos[b][1]],lon[pos[b][end]]]) for b in 1:length(pos)]
    sumDists = [sum(distance[p]) for p in pos]
    out[1:length(dists)] = reduce(vcat,dists./sumDists)
    return out
end

# add foraging number to foraging data
function forNo(df)
    forageNo = repeat([NaN],nrow(df))
    forEd = findall(diff(df.forage) .== -1)
    if df.forage[1] == 1
        append!(forEd,nrow(df))
    end
    for b = 1:length(forEd)
        if b == 1
            forageNo[1:forEd[b]] .= Int(b)
        else
            forageNo[(forEd[b-1]+1):forEd[b]] .= Int(b)
        end
    end
    return forageNo
end

# ROUND TO NEAREST n
function roundNearest(x::Float64,n::Float64)
    round(x/n)*n
end

# BIN VALUES x INTO bins
function bin(x::Float64, bins::AbstractArray)
    pos = findlast(x .> bins)
    if !isnothing(pos)
        if bins[pos] == maximum(bins)
            out = string(maximum(bins)) * "+"
        else 
            out = string(bins[pos]) * " - " * string(bins[pos + 1])
        end
    else 
        out = missing
    end
    return out
end

# add foraging number to wind estimates
function importForNo(df,wf)
    forNo = repeat([NaN],nrow(wf))
    for b = 1:nrow(wf)
        forNo[b] = df.forNo[findNearest(df.DT,wf.DT[b])]
    end
    return forNo
end

# calculate relative wind heading
function rwh(head,x,y)
    rwh = head - atan(y,x)
    if rwh < -pi
        rwh = rwh + (2*pi)
    elseif rwh > pi
        rwh = rwh - (2*pi)
    end
    return rwh
end
# align relative wind heading so 0 is headwind
function align(rwh)
    aligned = rwh + pi
    if aligned > pi
        aligned = aligned - (2*pi)
    end
    return aligned
end
# ADD TRIP LENGTHS
function tripNL(df::DataFrame,trLengths::Dict)
    tripN = Vector{Int}(undef,nrow(df))
    tripL = Vector{Int}(undef,nrow(df))
    yrid = df.yrID[1]
    trips = trLengths[yrid]
    tripEnds = cumsum(trips)
    tripStarts = cat(dims=1,1,tripEnds[1:end-1] .+ 1)
    tagDates = unique(Dates.Date.(df.DT))
    tripEnds[end] = minimum([length(tagDates),tripEnds[end]])
    for x = 1:length(tripEnds)
        tripN[tagDates[tripStarts[x]] .<= Dates.Date.(df.DT) .<= tagDates[tripEnds[x]]] .= Int(x)
        tripL[tagDates[tripStarts[x]] .<= Dates.Date.(df.DT) .<= tagDates[tripEnds[x]]] .= Int(trips[x])
    end
    return tripN,tripL
end

# file locations for foraging and wind estimates
if Sys.iswindows()
    dataloc = "E:/My Drive/PhD/Data/"
else
    dataloc = "/Volumes/GoogleDrive/My Drive/PhD/Data/"
end

# # bring in FORAGING AND GPS DATA and WIND DATA (remove wind data by setting true to false)
fDat,wDat = readDat(dataloc,".*PredictedForage[\\\\|/](.*)_S.*-20.*",true,[:DT,:lat,:lon,:forage,:yrID],[1,0],[2018,2019]);

for x in fDat
    insertcols!(x, ncol(x) + 1, :distTrav => [dist(x.lat,x.lon);NaN])
    insertcols!(x, ncol(x) + 1, :spTrav => [speed(x.DT,x.lat,x.lon);NaN])
    insertcols!(x, ncol(x) + 1, :lin => linearity(x.DT,x.lat,x.lon,x.distTrav,51))
end

plot(fDat[6].lon,
fDat[6].lat,
color=:greens,
label="",
line_z=fDat[6].lin,
colorbar_title="Linearity (5 mins)",
xlabel="Lon",
ylabel="Lat")

# assign data to the wind data
df = vcat(fDat...)
for wf in wDat
    wf.timeFP = repeat([NaN],nrow(wf))
    wf.distFP = repeat([NaN],nrow(wf))
    wf.speed = repeat([NaN],nrow(wf))
    wf.lin = repeat([NaN],nrow(wf))
    for b = 1:nrow(wf)
        if any((df.DT .> wf.DT[b]) .& (df.forage .== 1) .& (df.yrID .== wf.yrID[b]))
            nxtFor = findfirst((df.DT .> wf.DT[b]) .& (df.forage .== 1) .& (df.yrID .== wf.yrID[b]))
            wf.timeFP[b] = Dates.value(Second(df.DT[nxtFor] - wf.DT[b]))
            wf.distFP[b] = dist([wf.lat[b],df.lat[nxtFor]],[wf.lon[b],df.lon[nxtFor]])[1]
            wf.speed[b] = mean(df.spTrav[(df.DT .> (wf.DT[b] - Second((5*60)/2))) .& (df.DT .< (wf.DT[b] + Second((5*60)/2))) .& (df.yrID .== wf.yrID[b])])
            wf.lin[b] = mean(df.lin[(df.DT .> (wf.DT[b] - Second((5*60)/2))) .& (df.DT .< (wf.DT[b] + Second((5*60)/2))) .& (df.yrID .== wf.yrID[b])])
        else
            break
        end
    end
end

# add foraging number to foraging data
[insertcols!(df,ncol(df)+1, :forNo => forNo(df)) for df in fDat]

# add foraging number to wind data
[insertcols!(wDat[b],ncol(wDat[b])+1, :forNo => importForNo(fDat[b],wDat[b])) for b in 1:length(fDat)];

# 10km distance bins to FP
bins = [0:10:roundNearest(maximum(vcat(wDat...).distFP[isnan.(vcat(wDat...).distFP) .== false]/1000),10.0);]
[insertcols!(wf,ncol(wf)+1, :distBin =>  bin.(wf.distFP./1000, Ref(bins))) for wf in wDat];

# calculate relative wind heading
[insertcols!(wf, ncol(wf)+1, :rwh => rwh.(wf.head,wf.X,wf.Y)) for wf in wDat];

# align relative wind heading so 0 is headwind
[insertcols!(wf, ncol(wf)+1, :aligned => align.(wf.rwh)) for wf in wDat];

tripLengths = Dict("2018_10" => [1,4,1],
"2018_11" => [1,3,1,1],
"2018_1" => [2,5,1,4],
"2018_2017-9" => [1,4],
"2018_3" => [3,1,4,1,1,1],
"2018_4" => [3,1,6],
"2018_5" => [8,1],
"2018_6" => [4,5],
"2018_7" => [2,5,1,1,2],
"2018_8" => [3,1,5,1],
"2018_9" => [2,1,4,2],
"2019_1" => [1,1,1,1],
"2019_2018-01" => [1,1,1,1,2],
"2019_2018-03" => [1,4],
"2019_2018-04" => [1,2,1,1],
"2019_2018-05" => [1,4],
"2019_2" => [1,3],
"2019_3" => [1,1,1,1,2],
"2019_4" => [1,1,1,2],
"2019_5" => [1,1,3])

for df in fDat
    tripNo,tripLe = tripNL(df,tripLengths)
    insertcols!(df,:tripN => tripNo, :tripL => tripLe)
end

# transfer to wind data
for wf in 1:length(wDat)
    wDat[wf].tripN = Vector{Int}(undef,nrow(wDat[wf]))
    wDat[wf].tripL = Vector{Int}(undef,nrow(wDat[wf]))
    for x in 1:nrow(wDat[wf])
        wDat[wf].tripN[x] = fDat[wf].tripN[findNearest(wDat[wf].DT[x],fDat[wf].DT)]
        wDat[wf].tripL[x] = fDat[wf].tripL[findNearest(wDat[wf].DT[x],fDat[wf].DT)]
    end
end

# find the values of all data within 1 hour of foraging
allfd = vcat(fDat...)
ftimes = allfd[findall(allfd.forage .== 1),[:DT,:yrID]]
# add distance sum to each foraging spot
allfd.distFP .= NaN
fPoints = findall(allfd.forage .== 1)
for b in 1:length(fPoints)
    allfd.distFP[((allfd.DT[fPoints[b]] - Minute(60)) .< allfd.DT .< allfd.DT[fPoints[b]]) .& (allfd.yrID .== allfd.yrID[fPoints[b]])] = distToOne([allfd.lat[fPoints[b]],allfd.lon[fPoints[b]]],
        allfd.lat[((allfd.DT[fPoints[b]] - Minute(60)) .< allfd.DT .< allfd.DT[fPoints[b]]) .& (allfd.yrID .== allfd.yrID[fPoints[b]])],
        allfd.lon[((allfd.DT[fPoints[b]] - Minute(60)) .< allfd.DT .< allfd.DT[fPoints[b]]) .& (allfd.yrID .== allfd.yrID[fPoints[b]])])
end
a = DataFrame()
for (j, k) in zip(ftimes.DT,ftimes.yrID)
    append!(a, allfd[((j - Minute(60)) .< allfd.DT .< j ) .& (allfd.yrID .== k),:],cols=:union)
end

scatter(a.distFP[a.distFP .< 1000],a.lin[a.distFP .< 1000])

allwf = vcat(wDat...)
@df allwf[(allwf.distFP .< 100*1000) .& (allwf.tripL .> 2),:] density(:aligned, group = (:distBin))


findall(abs.(diff(allwf.DT)) .> Second(70))
function consec(dt,gap)
    out = fill(0,length(dt))
    sep = findall(abs.(diff(allwf.DT)) .> Second(gap))
    for b in 1:length(sep)
        if b == 1
            out[1:sep[b]] .= b
        elseif b == length(sep)
            out[(sep[b-1]+1):sep[b]] .= b
            out[(sep[b]+1):length(dt)] .= b
        else
            out[(sep[b-1]+1):sep[b]] .= b
        end
    end
    return out
end
allwf.seq = consec(allwf.DT,70)    

# remove all points where bird not approaching FP
allwf = allwf[all(!isnan, (allwf.distFP); dims = 2),:]


# @rput wDat
# nbStates <- 3 

