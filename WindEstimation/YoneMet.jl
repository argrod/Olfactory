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

data = [(-2.0,1.0), (0.0, 3.0), (2.0, 1.0)]
	
ptrs = 1:length(data)
	
y = zeros(length(data))



ea(c) = wy([c[1] c[2]]) * sind(gd) + wx([c[1] c[2]]) * cosd(gd) + 
    sqrt((wy([c[1] c[2]]) * sind(gd) + wx([c[1] c[2]]) * cos(gd))^2) - wx([c[1] c[2]])^2 - wx([c[1] c[2]])^2 + c[3]^2 - vg


# wy * sin(gd) + wx*cos(gd) + sqrt((wy*sin(gd) + wx*cos(gd))^2) - wx^2 - wx^2

# MATLAB VERSION
# wind x and y component
@. wx(v) = (v[1] .* cos(v[2]));
@. wy(v) = (v[1] .* sin(v[2]));
# ground speed as a function of track direction (gd) (polar axis)
ea = @(c)(wy([c(1) c(2)]).*sin(gd) + wx([c(1) c(2)]).*cos(gd) + ...
    sqrt((wy([c(1) c(2)]).*sin(gd)+wx([c(1) c(2)]).*cos(gd)).^2 ...
    - wy([c(1) c(2)]).^2 - wx([c(1) c(2)]).^2 + c(3)^2) - vg);

@. model(vg, vd, c) = wy([c[1] c[2]])

wy([15,13])

# THE CIRCLE function
circfun(vg,gd,c) = wy()


xp = x .+ k.*x.*(x.^2+y.^2)
yp = y .+ k.*y.*(x.^2+y.^2)
x = [x for x=-2.:2. for y=-2.:2.]
y = [y for x=-2.:2. for y=-2.:2.]
in_data = [x y]
out_data = [xp yp]


# write it as a function
# function windCircFit(vg, gd, c0, asDeg = false)
    # if asDeg == true
        wx(v) = v[1] * cosd(v[2])
        wy(v) = v[1] * sind(v[2])
        


circwind(vg,gd,c) = wy([c[1],c[2]])*sin(gd) + wx([c[1],c[2]])*cos(gd) + sqrt((wy([c[1],c[2]])*sin(gd) + wx([c[1],c[2]])*cos(gd))^2 - wy([c[1],c[2]])^2 - wx([c[1],c[2]])^2 + c[3]^2) - vg


function scircwind(x,c)
    wx(v) = v[1] * cos(v[2])
    wy(v) = v[1] * sin(v[2])
    vg = x[1]
    gd = x[2]
    α = c[1]
    β = c[2]
    γ = c[3]
    wy([α,β])*sin(gd) + wx([α,β])*cos(gd) + sqrt((wy([α,β])*sin(gd) + wx([α,β])*cos(gd))^2 - wy([α,β])^2 - wx([α,β])^2 + γ^2) - vg
end

[spd[31] dir[31]][2]



scircwind([spd[ss[300]] dir[ss[300]]],[3,0,9])

[spd[ss[300]] dir[ss[300]]][1,1]

scircwind.([spd dir],Ref([3,0,9]))

[spd dir]

tst=[spd dir]
tst = Ref([3,0,9])



tst = tst[2]

@view tst[:,1]
tst[:,2]

[spd dir][2]


#         # circwind(vg,gd,c) = wy([c[1],c[2]]).*sind.(gd) + wx([c[1],c[2]]).cosd.(gd) + sqrt.((wy([c[1],c[2]]).*sind.(gd) + wx([c[1],c[2]]).*cosd.(gd)).^2 - wy([c[1],c[2]])^2 - wx([c[1],c[2]])^2 + c[3]^2) - vg
#     # else
#         wx(v) = v[1] * cosd(v[2])
#         wy(v) = v[1] * sind(v[2])
#         @. circwind(vg,gd,c) = wy([c[1],c[2]]).*sin.(gd) .+ wx([c[1],c[2]]).*cos.(gd) + sqrt.((wy([c[1],c[2]]).*sin.(gd) .+ wx([c[1],c[2]]).*cos.(gd)).^2 - wy([c[1],c[2]])^2 - wx([c[1],c[2]])^2 + c[3]^2) - vg
#     end
#     lb = [0, -pi, 0];
#     ub = [20, pi, 20];
#     vecx = spd[ss[300]:se[300]].*cos.(dir[ss[300]:se[300]]);
#     vecy = spd[ss[300]:se[300]].*sin.(dir[ss[300]:se[300]]);
#     fit = curve_fit(circwind,vecx,vecy,[3,0,9],lower=lb, upper = ub)
# circwind(spd[ss[300]:se[300]],dir[ss[300]:se[300]],[3,0,9])
# end


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

# spd2,dir2,vsShow,veShow = flightMask(sel.DT,sel.lat,sel.lon,tagDur,1)
# amount = zeros(Int,nrow(sel))
# for b = 1:length(vsShow)
#     amount[vsShow[b]:veShow[b]] .= 1
# end
plot(sel.lon,sel.lat,label="")
scatter!(sel.lon[amount .== 1],sel.lat[amount.==1],markersize = 2,markerstrokewidth=0,label="Flight")

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

function scircwind(x,c)
    wx(v) = v[1] * cos(v[2])
    wy(v) = v[1] * sin(v[2])
    vg =x[:,1][1]
    gd =x[:,2][1]
    α,β,γ = c
    wy([α,β])*sin(gd) + wx([α,β])*cos(gd) + sqrt((wy([α,β])*sin(gd) + wx([α,β])*cos(gd))^2 - wy([α,β])^2 - wx([α,β])^2 + γ^2) - vg
end

scircwind([spd[ss[300]] dir[ss[300]]],[3,0,9])

[spd[ss[300]] dir[ss[300]]][2]
[spd[ss[300]] dir[ss[300]]][:,2][1]
scircwind.([spd  dir],Ref([ 3, 0, 9]))

@view [spd[300] dir[300]][:,1]

scircwind([12 1],[3, 0, 9])

sin([spd[300] dir[300]][:,2][1])


[spd[ss[300]] dir[ss[300]]][1,1]

[scircwind([spd[x],dir[x]],[3,0,9]) for x = 1:length(spd)]

@view [spd[300],dir[300]][2]

scircwind([spd[300],dir[300]],[3,0,9])
scircwind.([spd,dir],Ref.([3,0,9]))

tst = Ref([3,0,9])

tst

f(x,y) = sum(x) + y
f([1,2,3], 4) #-> 10
f([1,2,3], 5) #-> 11

f.([1,2,3], [4,5])    # -> error

f.(([1,2,3],), [4,5])

scircwind.([spd,dir],([3,0,9],))

([3,0,9],)

eg = 300
egSel = sel[ss[eg]:se[eg],:]
plot(sel.lon[ss[eg]:se[eg]],sel.lat[ss[eg]:se[eg]])
xlabel!("Lon")
ylabel!("Lat")

plot(1:61,spd[ss[eg]:se[eg]])
xlabel!("Time (s)")
ylabel!("Speed (m/s)")

plot(sel.lon,sel.lat)
for x = 1:length(ss)
    scatter!(sel.lon[ss[x]:se[x]])
end

flRange = [ss[x]:se[x] for x = 1:length(ss)]

function windestimates(spd, dir, ss, se)
    vw = zeros(Float64,length(ss),1)
end

selDir = dir[ss[1]:se[2]]; selSpd = spd[ss[1]:se[2]]
sortInd = sortperm(selDir)
sortedDir = selDir[sortInd]
sortedSpd = selSpd[sortInd]
sinRel(θ,c) = c[1] .+ c[2].*sin.(θ) .+ c[3].*cos.(θ)
ft = curve_fit(sinRel,sortedDir,sortedSpd,[1.0,1.0,1.0]);

scatter(dir[ss[1]:se[2]],spd[ss[1]:se[2]])
plot!(sortedDir,sinRel(sortedDir,ft.param),ylim=(0,15))

length(sinRel(dir[ss[1]:se[2]].*(180/pi).+180,ft.param))
plot(sinRel(1:360,[5.7,0.7,2.6]))
ft.param
gd = dir[ss[1]:se[1]];
vg = spd[ss[1]:se[1]];
scatter(gd,vg,xlim=(-pi,pi))
xlabel!("Track direction")
ylabel!("Track speed (m/s)")
wx(v) = v[1] .* cos(v[2])
wy(v) = v[1] .* sin(v[2])
ea(gd,c) = wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd) .+ sqrt.((wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd)).^2 .- wy([c[1] c[2]])^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
lb = [0.0, -pi, 0.0]
ub = [20.0, pi, 20.0]
curve_fit(ea,vg,gd,[3.0, 0.0, 9.0],lower = lb, upper = ub)
function wind2dveclsq(vg,gd,c0)
    wx(v) = v[1] .* cos(v[2])
    wy(v) = v[1] .* sin(v[2])
    ea(gd,c) = wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd) .+ sqrt.((wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
    lb = [0 -pi 0]
    ub = [20 pi 20]
    ydata = ea(c0)
    curve_fit(ea,vg,gd,c0,lower = lb, upper = ub)
end
wind2dveclsq(vg,gd,[3,0,9])
f(p, c) = wy([c[1] c[2]]).*sin.(p) .+ wx([c[1] c[2]]).*cos.(p) .+ sqrt.((wy([c[1] c[2]]).*sin.(p) .+ wx([c[1] c[2]]).*cos.(p)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
curve_fit(f,gd,vg,[3.0 0.0 9.0])
wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd) .+ sqrt.((wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)

wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd) .+ sqrt.((wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd)).^2 .- wy([c[1] c[2]])^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)

data = [(-2.0,1.0), (0.0, 3.0), (2.0, 1.0)]
	
ptrs = 1:length(data)
	
y = zeros(length(data))
	
function circle_model(t,p)
	out = []
	x0,y0,r = p
	for ptr ∈ ptrs
		x,y = data[ptr]
		push!(out,(x-x0)^2+(y-y0)^2-r^2)
	end
	out
end
	
p0 = [0.0, 0.0, 1.0] # using 0.0 as initial radius fails
	
fit = curve_fit(circle_model, ptrs, y, p0)
	
fit.param