using Pkg
Pkg.add(url="https://github.com/anowacki/CircStats.jl")
Pkg.add.(["Plots",
"CSV",
"GRIB",
"DataFrames",
"Query",
"Dates",
"CSV",
"Geodesy",
"StructArrays",
"RecursiveArrayTools",
"Statistics",
"RCall",
"NetCDF",
"DelimitedFiles",
"Glob",
"DataFrames",
"PlotlyJS",
"LsqFit"])

using Plots, CSV, GRIB, DataFrames, Query, Dates, CSV, Geodesy, StructArrays, RecursiveArrayTools, Statistics, RCall, NetCDF, DelimitedFiles, Glob, DataFrames, LsqFit

if Sys.iswindows()
    fileloc19 = "E:/My Drive/PhD/Data/2019Shearwater/AxyTrek/"
else
    fileloc19 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/AxyTrek/"
end
files = glob("*/*.txt",fileloc19)
EDat = DataFrame()
for b = 1:length(files)
    append!(EDat, hcat(CSV.File(files[b],header=false,delim="\t") |> DataFrame,fill(readdir(fileloc19)[b],nrow(CSV.File(files[b],header=false,delim="\t") |> DataFrame))))
end

rename!(EDat,[:DT,:lat,:lon,:c4,:c5,:c6,:c7,:c8,:ID])
dt = dateformat"y/m/d,H:M:S"
EDat.DT = DateTime.(EDat.DT,dt);

tags = unique(EDat.ID);
sel = EDat[findall(x->x==tags[1],EDat.ID),:];

# data = [(-2.0,1.0), (0.0, 3.0), (2.0, 1.0)]
# ptrs = 1:length(data)
#y = zeros(length(data))
# ea(c) = wy([c[1] c[2]]) * sind(gd) + wx([c[1] c[2]]) * cosd(gd) + 
    sqrt((wy([c[1] c[2]]) * sind(gd) + wx([c[1] c[2]]) * cos(gd))^2) - wx([c[1] c[2]])^2 - wx([c[1] c[2]])^2 + c[3]^2 - vg
# wy * sin(gd) + wx*cos(gd) + sqrt((wy*sin(gd) + wx*cos(gd))^2) - wx^2 - wx^2
# MATLAB VERSION
# wind x and y component
# @. wx(v) = (v[1] .* cos(v[2]));
# @. wy(v) = (v[1] .* sin(v[2]));
# # ground speed as a function of track direction (gd) (polar axis)
# ea = @(c)(wy([c(1) c(2)]).*sin(gd) + wx([c(1) c(2)]).*cos(gd) + ...
#     sqrt((wy([c(1) c(2)]).*sin(gd)+wx([c(1) c(2)]).*cos(gd)).^2 ...
#     - wy([c(1) c(2)]).^2 - wx([c(1) c(2)]).^2 + c(3)^2) - vg);
# @. model(vg, vd, c) = wy([c[1] c[2]])
# THE CIRCLE function
# circfun(vg,gd,c) = wy()

# write it as a function
# function windCircFit(vg, gd, c0, asDeg = false)
    # if asDeg == true
        wx(v) = v[1] * cosd(v[2])
        wy(v) = v[1] * sind(v[2])
        

# function scircwind(x,c)
#     wx(v) = v[1] * cos(v[2])
#     wy(v) = v[1] * sin(v[2])
#     vg = x[1]
#     gd = x[2]
#     α = c[1]
#     β = c[2]
#     γ = c[3]
#     wy([α,β])*sin(gd) + wx([α,β])*cos(gd) + sqrt((wy([α,β])*sin(gd) + wx([α,β])*cos(gd))^2 - wy([α,β])^2 - wx([α,β])^2 + γ^2) - vg
# end

# scircwind([spd[ss[300]] dir[ss[300]]],[3,0,9])

function flightMask(time,lat,lon,threshold,consecutive)
    # get speeds, directions, and a flight mask
    # requires datetime of data (in datetime format), latitudel, longitude, speed threshold for flight (m/s), and how many slow speed values constitute rest
    # first find UTM values
    Ll = LLA.(lat,lon)
    UTMD = UTMZfromLLA(wgs84)
    utmDat = map(UTMD, Ll); utmDat = DataFrame(utmDat)
    X = diff(utmDat.x)
    Y = diff(utmDat.y)
    distTrav = sqrt.(X.^2 + Y.^2)
    tdiff = Dates.value.(Second.(diff(time)))
    dir = atan.(Y,X)
    spd = distTrav./tdiff
    rest = spd .< threshold
    a = findall(x -> x == true, [false; rest])
    tvs = findall(x -> x == 1, diff([0;rest]))
    tve = findall(x -> x == -1, diff([rest;0]))
    rmv = findall(x -> x .<= consecutive, (tve .- tvs .+ 1))
    for x in rmv
        rest[tvs[x]:tve[x]] .= false
    end
    flight = rest .== false
    vs = findall(x -> x .== 1, diff([0;flight]))
    ve = findall(x -> x .== -1, diff([flight;0]))
    return spd, dir, vs, ve, rest
end
spd, dir, vs, ve, rest = flightMask(sel.DT,sel.lat,sel.lon,5,2);
# separate flight data into section for wind estimation
function getsection(F,Win,delta,fs,fe,time)
    # F - sampling frequency (Hz)
    # Win - sampling window (seconds)
    # delta - window step length (seconds)
    # fs and fe - output from flightMask
    window = Win * F
    del = delta * F
    
    # use only segments that exceed 7 minutes in length
    fsa = fs[(fe .- fs) .>= (7 * 60 * F)] .+ (60 * F)
    fea = fe[(fe .- fs) .>= (7 * 60 * F)] .- (60 * F)
    
    secnum = sum(floor.((fea .- fsa .- window)./del) .+ 1)
    
    ss = zeros(Int,trunc(Int,secnum),1)
    se = zeros(Int,trunc(Int,secnum),1)
    k = 1
    while k <= secnum
        for i = 1:length(fsa)
            for c = 0:floor((fea[i] - fsa[i] - window)./del)
                # ensure time difference is less than 6 minutes
                if time[trunc(Int,(fsa[i] + c * del + window))] - time[trunc(Int,(fsa[i] + c * del))] > Minute(6)
                    ss[k] = 0
                    se[k] = 0
                else
                    ss[k] = fsa[i] + c * del
                    se[k] = fsa[i] + c * del + window
                end
                k = k + 1
            end
        end
    end
    ss = filter(x->x≠0,ss)
    se = filter(x->x≠0,se)
    return ss,se
end
ss,se = getsection(.2,300,60,vs,ve,sel.DT);
function scircwind(gd,c)
    wx(v) = v[1] * cos(v[2])
    wy(v) = v[1] * sin(v[2])
    α,β,γ = c
    wy([α,β]).*sin.(gd) .+ wx([α,β]).*cos.(gd) .+ sqrt.((wy([α,β]).*sin.(gd) .+ wx([α,β]).*cos.(gd)).^2 .- wy([α,β]).^2 .- wx([α,β]).^2 .+ γ.^2)
end

function circwind(gd,c)
    wx(v) = v[1] * cos(v[2])
    wy(v) = v[1] * sin(v[2])
    α,β,γ = c
    wy([α,β])*sin(gd) + wx([α,β])*cos(gd) + sqrt((wy([α,β])*sin(gd) + wx([α,β])*cos(gd))^2 - wy([α,β])^2 .- wx([α,β]).^2 + γ^2)
end

circwind.(dir[ss[13]:se[13]],Ref([3,0,9]))

curve_fit(circwind, dir[ss[13]:se[13]], spd[ss[13]:se[13]],[3,0,9],lower = [0.0, -pi, 0.0], upper = [20.0, pi, 20.0])

scircwind.(dir[ss[13]:se[13]],[3,0,9])

fit=curve_fit(scircwind,dir[ss[g]:se[g]],spd[ss[g]:se[g]],[3.0,0.0,9.0],lower = [0.0, -pi, 0.0], upper = [20.0, pi, 20.0])
fit.param

Plots.plot(sel.lon,sel.lat)

fit = curve_fit(scircwind, dir[ss[13]:se[13]], spd[ss[13]:se[13]],[3.0,0.0,9.0],lower = [0.0, -pi, 0.0], upper = [20.0, pi, 20.0])
fit.param

Plots.scatter(dir[ss[13]:se[13]],spd[ss[13]:se[13]])
Plots.plot!([-pi:.01:pi],scircwind(-pi:.01:pi,fit.param))

Plots.plot(sel.lon,sel.lat,legend=false)
Plots.scatter!(sel.lon[reduce(vcat, ranges)],sel.lat[reduce(vcat, ranges)],markersize=3,legend=false)

findall(ss.==1318)
se

sel.lon[ranges]

[collect(ss[g]:se[g]) for g in 1:length(ss)]

ranges = [collect(ss[g]:se[g]) for g = 1:length(ss)]


identity(ranges)
ss:se


function calcWind(Dat::DataFrame,ID::String)
    sel = Dat[Dat.ID .== ID,:]
    spd, dir, vs, ve, rest = flightMask(sel.DT,sel.lat,sel.lon,5,2)
    ss,se = getsection(.2,300,60,vs,ve,sel.DT)
    fits = [curve_fit(scircwind, dir[ss[g]:se[g]], spd[ss[g]:se[g]],[3.0,0.0,9.0],lower = [0.0, -pi, 0.0], upper = [20.0, pi, 20.0]).param for g = 1:length(ss)].param
        fits = fit.param
 
fits = [0.0 0.0 0.0]
for b = 1:length(unique(EDat.ID))
    sel = EDat[EDat.ID .== tags[b],:]
    spd, dir, vs, ve, rest = flightMask(sel.DT,sel.lat,sel.lon,5,2)
    ss,se = getsection(.2,300,60,vs,ve,sel.DT)
    