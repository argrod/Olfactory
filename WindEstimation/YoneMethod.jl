### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b2b504be-f728-42d5-bb28-04c7a566cf5d
begin
	using Plots
	using CSV
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
	using DelimitedFiles
	using Glob
	using DataFrames
	using LsqFit
	using PlutoUI
	using LaTeXStrings
end

# ╔═╡ 24dd02f0-f348-11eb-02fe-3309b1be5acc
md"""
# Fitting circle to track vector plot

Following discussions with 米原さん, after the [PNAS publication](https://www.pnas.org/content/113/32/9039 "Flight paths of seabirds soaring over the ocean surface enable measurement of fine-scale wind speed and direction") he tested another method of wind estimation, instead fitting a circle to a vector plot of the ground speed. Deviation of the fitted circle from $(0,0)$ indicates the wind vector.
"""

# ╔═╡ 0e34f851-4e25-4b6a-b76f-672b311a6819
# READ IN AND FORMAT
begin
	if Sys.iswindows()
		fileloc19 = "F:/UTokyoDrive/PhD/Data/2019Shearwater/AxyTrek/"
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
end

# ╔═╡ 19ad02e2-a725-42e4-8d3e-ecaa93a559a7
begin
	tags = unique(EDat.ID);
	sel = EDat[findall(x->x==tags[1],EDat.ID),:];
end

# ╔═╡ 0c8048da-df98-4e4a-ba2e-195ede0d6e5c
md"""
## Example tag

Plot of the tracks of inidiviual $(EDat.ID[1]) 
"""

# ╔═╡ 6f21012b-39d6-47d4-99b3-0dde5068c357
plot(sel.lon,sel.lat)

# ╔═╡ 62ee1055-dbd9-4f93-990d-ce40092e2b1b
begin
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
end

# ╔═╡ 99f6bc95-5032-46ef-bac5-5b01bd55846c
spd, dir, vs, ve, rest = flightMask(sel.DT,sel.lat,sel.lon,5,5);

# ╔═╡ 08cd4275-db5a-4553-bf0a-7b71a8927eb7
md"""
## Separating the track into flight segments

Using the criteria of a ground speed of $> 4 \textrm{ms}^{-1}$, the track is split into flight or rest. You can test other thresholds and their results below.
"""

# ╔═╡ cf3bd1d2-53ed-4b4c-ad8a-c3716b89b569
md"""
Flight speed threshold (ms$$^{-1}$$) $(@bind tagDur Slider(1:50, show_value = true))
"""

# ╔═╡ 1ad6e21c-05d1-4fa6-b21e-abe6858f1ed6
begin
	spd2,dir2,vsShow,veShow = flightMask(sel.DT,sel.lat,sel.lon,tagDur,5)
	amount = zeros(Int,nrow(sel))
	for b = 1:length(vsShow)
		amount[vsShow[b]:veShow[b]] .= 1
	end
	plot(sel.lon,sel.lat,label="")
	scatter!(sel.lon[amount .== 1],sel.lat[amount.==1],markersize = 2,markerstrokewidth=0,label="Flight")
end

# ╔═╡ 7c23b6d9-e659-430f-bc5e-ca24ef58196d
begin
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
end

# ╔═╡ 87a66166-31c6-460a-b534-a219d8d52a79
ss,se = getsection(.2,300,60,vs,ve,sel.DT);

# ╔═╡ a7cf59c7-74bd-4c55-a1bb-74fde41c7a6a
begin
	eg = 300
	egSel = sel[ss[eg]:se[eg],:]
	plot(sel.lon[ss[eg]:se[eg]],sel.lat[ss[eg]:se[eg]])
	xlabel!("Lon")
	ylabel!("Lat")
end

# ╔═╡ 80e537c2-76a1-4507-a83b-67e76ad873a5
begin
	plot(1:61,spd[ss[eg]:se[eg]])
	xlabel!("Time (s)")
	ylabel!("Speed (m/s)")
end

# ╔═╡ 9e7a192a-c823-467b-b743-36daf0337aa8


# ╔═╡ b012fbab-ddf1-4b58-840b-b5ab278f7150


# ╔═╡ 7020789b-aaad-431e-8ee9-d929852f70e1


# ╔═╡ 76b5a853-5e04-497c-a3c2-c3af800f1909
[3 0 9]

# ╔═╡ 2cf89d73-4442-4084-a1e9-2c85b38095de


# ╔═╡ 8526fc1f-c3da-4676-b804-bc98abd30290


# ╔═╡ 9f07e1c6-5afb-4ac9-b571-e22161ac84d1


# ╔═╡ 0b0bac61-86c9-4ce1-9f2d-a8a7c609ad76


# ╔═╡ 2045c3bb-9b50-44b5-8e9d-477f4717dd8a


# ╔═╡ 8e0b0f75-536c-4531-b346-39821ef10fec


# ╔═╡ ebab33b1-3aff-41d6-b04e-176f02342b35


# ╔═╡ 8f20f3bb-df17-4a62-a20a-857e9268abed


# ╔═╡ 89a89bb0-acc7-46e1-bcbd-39402e075c11


# ╔═╡ 205cd474-d5d0-4470-b56d-2ee9d37ed0ed
wy([3 0]*sin.(gd)

# ╔═╡ 77351942-34e9-44c9-95c5-94e210288bbd


# ╔═╡ 47fb8578-4c00-41c0-8586-515d8e8a6cfe


# ╔═╡ 9236a1f8-0e22-42bc-9a07-ef8ef0bf25bf


# ╔═╡ 3b69098f-3777-409e-a384-34b54b3c85e2


# ╔═╡ b1e7696b-9d09-44eb-86c8-2cfc9d4d9671


# ╔═╡ 729f983a-2794-48ea-9484-447c2982137b


# ╔═╡ 83f89135-fb60-4004-858b-955156c36a47


# ╔═╡ e0a93033-6792-46f0-ba01-a8e187f88955


# ╔═╡ dad257f4-afe4-4272-8a71-5d5cdb1c1afd


# ╔═╡ 4e0f5953-0d75-4b3b-ba06-2a4a562e56ec


# ╔═╡ 569163c3-e9d4-4edc-89cc-14ae2da12c92


# ╔═╡ 3fdd9c83-e2a4-4c05-8c05-3b8a66e1a1f7


# ╔═╡ 03fe0635-4e1d-4f87-8edc-080ca4935e8f


# ╔═╡ 8b9b9c15-eb5b-448b-835d-af7aa2a6448c


# ╔═╡ fb255db9-3fdd-48eb-a307-36e064e61ccd


# ╔═╡ 04ce9b62-abcb-4ed7-a35d-85565357932c
# wind2dveclsq(vg,gd,[3 0 9])

# ╔═╡ 11cdeaca-f0c1-47d2-9f4d-1ca0bbbc9d5e
# foist(v) = v[1] .* cos(v[2])
# 		secnd(v) = v[1] .* sin(v[2])
# 		full(c) = secnd([c[1] c[2]]).*sin(gd) + foist([c[1] c[2]]).*cos(gd) + sqrt((secnd([c[1] c[2]]).*sin(gd) + foist([c[1] c[2]]).*cos(gd)).^2 - secnd([c[1] c[2]]).^2 - foist([c[1] c[2]]).^2 + c[3]^2) - vg
# 		lb = [0 -pi 0]
# 		ub = [20 pi 20]
# 		ydata = ea(c0)

# ╔═╡ dc38cae8-9786-45fd-9e49-730c7b24ab90
begin 
	function windestimates(spd, dir, ss, se)
		vw = zeros(Float64,length(ss),1)
	end
end

# ╔═╡ bf268156-adae-4688-a94b-f2ed2e97d320
# begin
# 		wx(v) = v[1] .* cos(v[2])
# 		wy(v) = v[1] .* sin(v[2])
# 		ea(t,c) = wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd) .+ sqrt.((wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2) .- vg
# 		lb = [0 -pi 0]
# 		ub = [20 pi 20]
# 		# ydata = ea([0 0 0],gd,vg)
# 		# curve_fit(ea,[3 0 9],lower = lb, upper = ub)
# end

# ╔═╡ b5ea6fd2-a58c-4f78-a046-26081615170e
# curve_fit(ea,gx,gy,[3, 0, 9])

# ╔═╡ 9bc1d10b-5314-485b-86aa-3294d1c7f4c3
begin
	selDir = dir[ss[1]:se[2]]; selSpd = spd[ss[1]:se[2]]
	sortInd = sortperm(selDir)
	sortedDir = selDir[sortInd]
	sortedSpd = selSpd[sortInd]
	sinRel(θ,c) = c[1] .+ c[2].*sin.(θ) .+ c[3].*cos.(θ)
	ft = curve_fit(sinRel,sortedDir,sortedSpd,[1.0,1.0,1.0])
end

# ╔═╡ 9d5c901d-c3f9-4c56-bade-8a25b2be2fe7
ft.param

# ╔═╡ eaac9617-cdd5-49fb-80a2-d36b3ca15a68
begin
	scatter(dir[ss[1]:se[2]],spd[ss[1]:se[2]])
	plot!(sortedDir,sinRel(sortedDir,ft.param),ylim=(0,15))
end

# ╔═╡ 66f5865d-9236-4562-b7ce-cdc7440c4c0e


# ╔═╡ 800e1bad-8917-4d98-b0fb-73ef0469e1d6
length(sinRel(dir[ss[1]:se[2]].*(180/pi).+180,ft.param))

# ╔═╡ cacfb01e-31fc-4b0c-8b75-4a55fffdf4c8


# ╔═╡ 92ef42d9-ee38-4485-bbca-d5b4303dbb2a
plot(sinRel(1:360,[5.7,0.7,2.6]))

# ╔═╡ 2381a54b-ba61-468f-9287-bf3e2024b694


# ╔═╡ 493fba7a-e07f-496e-b40a-1f214dd4dfb9
ft.param

# ╔═╡ 115295f1-047a-4da4-a348-1b89e8c795bb
c = [3 0 9]

# ╔═╡ eed56293-6b89-416c-96a8-685b5957b4f2
begin
	gd = dir[ss[1]:se[1]];
	vg = spd[ss[1]:se[1]];
	
end

# ╔═╡ bb33aaed-8f48-4895-ae32-835d6d92508d
begin
	scatter(gd,vg,xlim=(-pi,pi))
	xlabel!("Track direction")
	ylabel!("Track speed (m/s)")
end

# ╔═╡ 82cf3981-1648-4430-baca-2fce80e40a6b
begin
	gx = vg .* cos.(gd)
	gy = vg .* sin.(gd)
	scatter(gx,gy)
	xlabel!("Track U component (m/s)")
	ylabel!("Track V component (m/s)")
end

# ╔═╡ ba6daf12-03d7-4714-800f-3ab12d02d5a2
begin
		wx(v) = v[1] .* cos(v[2])
		wy(v) = v[1] .* sin(v[2])
		ea(gd,c) = wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd) .+ sqrt.(wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd).^2 .- wy([c[1] c[2]])^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
		lb = [0, -pi, 0]
		ub = [20, pi, 20]
		curve_fit(ea,vg,gd,[3.0, 0.0, 9.0],lower = lb, upper = ub)
end

# ╔═╡ cabe59d8-b22e-4b5b-ad6e-7c01f709c4ea
begin
	function wind2dveclsq(vg,gd,c0)
		wx(v) = v[1] .* cos(v[2])
		wy(v) = v[1] .* sin(v[2])
		ea(gd,c) = wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd) .+ sqrt.((wy([c[1] c[2]]).*sin.(gd) .+ wx([c[1] c[2]]).*cos.(gd)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
		lb = [0 -pi 0]
		ub = [20 pi 20]
		ydata = ea(c0)
		curve_fit(ea,vg,gd,c0,lower = lb, upper = ub)
	end
end

# ╔═╡ 9658da94-80ea-4454-96ee-ecddd59ca2a6
wind2dveclsq(vg,gd,[3,0,9])

# ╔═╡ 48f99cd6-752f-4f99-80b2-6695369710f3
begin
	wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd) .+ sqrt.(wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd).^2 .- wy([c[1] c[2]])^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)
end

# ╔═╡ ac03761e-9d66-4f7f-874d-9c18de325336
wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd) .+ sqrt.(wy([c[1] c[2]])*sin.(gd) .+ wx([c[1] c[2]])*cos.(gd).^2 .- wy([c[1] c[2]])^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2)

# ╔═╡ cffa7c7b-b5b8-4636-93be-ff7bd155f961
length(wy([3 0])*sin.(gd) .+ wx([3 0])*cos.(gd) .+ sqrt.(wy([3 0])*sin.(gd) .+ wx([3 0])*cos.(gd).^2 .- wy([3 0])^2 .- wx([3 0]).^2 .+ c[3]^2))

# ╔═╡ 901d3194-b035-4154-9ed3-709320e114c1
vg

# ╔═╡ d015f470-943d-42f4-ad57-c32255e8e4e6
# begin 
# 	# a two-parameter exponential model
# # x: array of independent variables
# # p: array of model parameters

# model(c) = wy([c[1] c[2]]).*sin.(gd) + wx([c[1] c[2]]).*cos.(gd) .+ sqrt.((wy([c[1] c[2]]).*sin.(gd) + wx([c[1] c[2]]).*cos.(gd)).^2 .- wy([c[1] c[2]]).^2 .- wx([c[1] c[2]]).^2 .+ c[3]^2) .- vg

# # some example data
# # xdata: independent variables
# # ydata: dependent variable
# # xdata = range(0,stop=10,length=20)
# ydata = model([3 0 9])
# # p0 = [0.5, 0.5]

# fit = curve_fit(model, [3 0 9], ydata,[3 0 9])
# end

# ╔═╡ bdcadc88-cfda-48b3-b678-3d6e1db9effc
md"""
## Theory

With no wind effect on the bird as it flies, it should, theoretically, fly at the same speed regardless of the direction. Therefore, the track speed ($t_s$) and direction ($t_d$) should generate a relationship such that $t_s^2 + t_d^2 = r^2$ where $r$ is the radius of a circle centered around the point ($0,0$).

The addition of winds should move the origin point such that the speed of travel should be dictated by the travel direction (i.e. a westerly wind would increase travel speed to the east and decrease it to the west, shifting the center further to the right).

The wind vector made up of wind speed ($w_v$) and wind direction ($w_d$) can be converted into $u$ and $v$ components via:

$u = w_v \times \cos(w_d)$
$v = w_v \times \sin(w_d)$

Then a function can be written of the ground speed ($v_d$) as a function of track direction ($g_d$):

$f(w_d,w_v,a_v) = v\sin(g_d)+u\cos(g_d)+\sqrt{\big(v\sin(g_d)+u\cos(g_d)\big)^2-v^2-u^2}-v_g$

where $a_v$ is the air speed.
push!(out,(x-x0)^2+(y-y0)^2-r^2)
"""

# ╔═╡ 61470dff-2413-4c6f-97d5-2ec739774b1c
begin
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
end

# ╔═╡ 772339c6-6cd8-4936-be49-6b9f3a782eb8
plot(data); scatter!(fit.param)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
GRIB = "b16dfd50-4035-11e9-28d4-9dfe17e6779b"
Geodesy = "0ef565a4-170c-5f04-8de2-149903a85f3d"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
NetCDF = "30363a11-5582-574a-97bb-aa9a979735b9"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
CSV = "~0.8.5"
DataFrames = "~1.2.2"
GRIB = "~0.3.0"
Geodesy = "~1.0.1"
Glob = "~1.3.0"
LaTeXStrings = "~1.2.1"
LsqFit = "~0.12.1"
NetCDF = "~0.11.3"
Plots = "~1.20.0"
PlutoUI = "~0.7.9"
Query = "~1.0.0"
RCall = "~0.13.12"
RecursiveArrayTools = "~2.16.1"
StructArrays = "~0.6.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "2e004e61f76874d153979effc832ae53b56c20ee"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.22"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CSV]]
deps = ["Dates", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode"]
git-tree-sha1 = "b83aa3f513be680454437a0eee21001607e5d983"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.8.5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[CategoricalArrays]]
deps = ["DataAPI", "Future", "JSON", "Missings", "Printf", "RecipesBase", "Statistics", "StructTypes", "Unicode"]
git-tree-sha1 = "1562002780515d2573a4fb0c3715e4e57481075e"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f53ca8d41e4753c41cdafa6ec5f7ce914b34be54"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "0.10.13"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random", "StaticArrays"]
git-tree-sha1 = "ed268efe58512df8c7e224d2e170afd76dd6a417"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.13.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "344f143fa0ec67e47917848795ab19c6a455f32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.32.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6d1c23e740a586955645500bbec662476204a52c"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.1"

[[Crayons]]
git-tree-sha1 = "3f71217b538d7aaee0b69ab47d9b7724ca8afa0d"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.0.4"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d785f42445b63fc86caa08bb9a9351008be9b765"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.2.2"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4437b64df1e0adccc3e5d1adbc3ac741095e4677"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.9"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "85d2d9e2524da988bffaf2a381864e20d2dae08d"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.2.1"

[[DiskArrays]]
git-tree-sha1 = "599dc32bae654fa78056b15fed9b2af36f04ee44"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.2.11"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3889f646423ce91dd1055a76317e9a1d3a23fff1"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.11"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8c8eac2af06ce35973c3eadb4ab3243076a408e7"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.1"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "b5e930ac60b613ef3406da6d4f42c35d8dc51419"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.19"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GRIB]]
deps = ["eccodes_jll"]
git-tree-sha1 = "bc1b389c9ee2821c714dfd24763fce0bcf38666c"
uuid = "b16dfd50-4035-11e9-28d4-9dfe17e6779b"
version = "0.3.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[Geodesy]]
deps = ["CoordinateTransformations", "Dates", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "2341a0d40d1f96db72275d7ae97ad5600778f137"
uuid = "0ef565a4-170c-5f04-8de2-149903a85f3d"
version = "1.0.1"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fd83fa0bde42e01952757f01149dd968c06c4dba"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.0+1"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "44e3b40da000eab4ccb1aecdc4801c040026aeb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.13"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "81690084b6198a2e1da36fcfda16eeca9f9f24e4"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.1"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "7bd5f6565d80b6bf753738d2bc40a5dfea072070"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.2.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f8c673ccc215eb50fcadb285f522420e29e69e1c"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "0.4.5"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "144bab5b1443545bc4e791536c9f1eacb4eed06a"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.1"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetCDF]]
deps = ["DiskArrays", "Formatting", "NetCDF_jll"]
git-tree-sha1 = "23b0e32fde256a4e2e497e678abcf956ed26204b"
uuid = "30363a11-5582-574a-97bb-aa9a979735b9"
version = "0.11.3"

[[NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "0cf4d1bf2ef45156aed85c9ac5f8c7e697d9288c"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.702.400+0"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "bfd7d8c7fd87f04543810d9cbd3995972236ba1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "e39bea10478c6aff5495ab522517fae5134b40e3"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.20.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "0d1245a357cc61c8cd61934c07447aa569ff22e6"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.1.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "12fbe86da16df6679be7521dfb39fbc861e1dc7b"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.1"

[[Query]]
deps = ["DataValues", "IterableTables", "MacroTools", "QueryOperators", "Statistics"]
git-tree-sha1 = "a66aa7ca6f5c29f0e303ccef5c8bd55067df9bbe"
uuid = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
version = "1.0.0"

[[QueryOperators]]
deps = ["DataStructures", "DataValues", "IteratorInterfaceExtensions", "TableShowUtils"]
git-tree-sha1 = "911c64c204e7ecabfd1872eb93c49b4e7c701f02"
uuid = "2aef5ad7-51ca-5a8f-8e88-e75cf067b44b"
version = "0.9.3"

[[RCall]]
deps = ["CategoricalArrays", "Conda", "DataFrames", "DataStructures", "Dates", "Libdl", "Missings", "REPL", "Random", "Requires", "StatsModels", "WinReg"]
git-tree-sha1 = "80a056277142a340e646beea0e213f9aecb99caa"
uuid = "6f49c342-dc21-5d91-9882-a32aef131414"
version = "0.13.12"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "0426474f50756b3b47b08075604a41b460c45d17"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.16.1"

[[Reexport]]
git-tree-sha1 = "5f6c21241f0f655da3952fd60aa18477cf96c220"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.1.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "a3a337914a035b2d59c9cbe7f1a38aaba1265b02"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.6"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "508822dca004bf62e210609148511ad03ce8f1d8"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.0"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "885838778bb6f0136f8317757d7803e0d81201e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.9"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[StatsFuns]]
deps = ["LogExpFunctions", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "30cd8c360c54081f806b1ee14d2eecbef3c04c49"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.8"

[[StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "a209a68f72601f8aa0d3a7c4e50ba3f67e32e6f8"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.24"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "e36adc471280e8b346ea24c5c87ba0571204be7a"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.7.2"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "8fc12ae66deac83e44454e61b02c37b326493233"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[WinReg]]
deps = ["Test"]
git-tree-sha1 = "808380e0a0483e134081cc54150be4177959b5f4"
uuid = "1b915085-20d7-51cf-bf83-8f477d6f5128"
version = "0.3.1"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"

[[eccodes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "OpenJpeg_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "bab4294c7f62efa25346bcc77c3e899fd20ee3a7"
uuid = "b04048ba-5ccd-5610-b3f6-85129a548705"
version = "2.21.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─24dd02f0-f348-11eb-02fe-3309b1be5acc
# ╠═b2b504be-f728-42d5-bb28-04c7a566cf5d
# ╟─0e34f851-4e25-4b6a-b76f-672b311a6819
# ╟─19ad02e2-a725-42e4-8d3e-ecaa93a559a7
# ╟─0c8048da-df98-4e4a-ba2e-195ede0d6e5c
# ╟─6f21012b-39d6-47d4-99b3-0dde5068c357
# ╟─62ee1055-dbd9-4f93-990d-ce40092e2b1b
# ╟─99f6bc95-5032-46ef-bac5-5b01bd55846c
# ╟─08cd4275-db5a-4553-bf0a-7b71a8927eb7
# ╟─cf3bd1d2-53ed-4b4c-ad8a-c3716b89b569
# ╟─1ad6e21c-05d1-4fa6-b21e-abe6858f1ed6
# ╟─7c23b6d9-e659-430f-bc5e-ca24ef58196d
# ╠═87a66166-31c6-460a-b534-a219d8d52a79
# ╟─a7cf59c7-74bd-4c55-a1bb-74fde41c7a6a
# ╟─80e537c2-76a1-4507-a83b-67e76ad873a5
# ╟─bb33aaed-8f48-4895-ae32-835d6d92508d
# ╟─82cf3981-1648-4430-baca-2fce80e40a6b
# ╠═cabe59d8-b22e-4b5b-ad6e-7c01f709c4ea
# ╠═9658da94-80ea-4454-96ee-ecddd59ca2a6
# ╠═ba6daf12-03d7-4714-800f-3ab12d02d5a2
# ╠═9e7a192a-c823-467b-b743-36daf0337aa8
# ╠═b012fbab-ddf1-4b58-840b-b5ab278f7150
# ╠═48f99cd6-752f-4f99-80b2-6695369710f3
# ╠═7020789b-aaad-431e-8ee9-d929852f70e1
# ╠═76b5a853-5e04-497c-a3c2-c3af800f1909
# ╠═2cf89d73-4442-4084-a1e9-2c85b38095de
# ╠═8526fc1f-c3da-4676-b804-bc98abd30290
# ╠═9f07e1c6-5afb-4ac9-b571-e22161ac84d1
# ╠═0b0bac61-86c9-4ce1-9f2d-a8a7c609ad76
# ╠═2045c3bb-9b50-44b5-8e9d-477f4717dd8a
# ╠═ac03761e-9d66-4f7f-874d-9c18de325336
# ╠═8e0b0f75-536c-4531-b346-39821ef10fec
# ╠═cffa7c7b-b5b8-4636-93be-ff7bd155f961
# ╠═ebab33b1-3aff-41d6-b04e-176f02342b35
# ╠═8f20f3bb-df17-4a62-a20a-857e9268abed
# ╠═89a89bb0-acc7-46e1-bcbd-39402e075c11
# ╠═205cd474-d5d0-4470-b56d-2ee9d37ed0ed
# ╠═77351942-34e9-44c9-95c5-94e210288bbd
# ╠═47fb8578-4c00-41c0-8586-515d8e8a6cfe
# ╠═9236a1f8-0e22-42bc-9a07-ef8ef0bf25bf
# ╠═3b69098f-3777-409e-a384-34b54b3c85e2
# ╠═b1e7696b-9d09-44eb-86c8-2cfc9d4d9671
# ╠═729f983a-2794-48ea-9484-447c2982137b
# ╠═83f89135-fb60-4004-858b-955156c36a47
# ╠═e0a93033-6792-46f0-ba01-a8e187f88955
# ╠═dad257f4-afe4-4272-8a71-5d5cdb1c1afd
# ╠═4e0f5953-0d75-4b3b-ba06-2a4a562e56ec
# ╠═569163c3-e9d4-4edc-89cc-14ae2da12c92
# ╠═3fdd9c83-e2a4-4c05-8c05-3b8a66e1a1f7
# ╠═03fe0635-4e1d-4f87-8edc-080ca4935e8f
# ╠═8b9b9c15-eb5b-448b-835d-af7aa2a6448c
# ╠═901d3194-b035-4154-9ed3-709320e114c1
# ╠═fb255db9-3fdd-48eb-a307-36e064e61ccd
# ╠═04ce9b62-abcb-4ed7-a35d-85565357932c
# ╠═11cdeaca-f0c1-47d2-9f4d-1ca0bbbc9d5e
# ╠═dc38cae8-9786-45fd-9e49-730c7b24ab90
# ╠═bf268156-adae-4688-a94b-f2ed2e97d320
# ╠═b5ea6fd2-a58c-4f78-a046-26081615170e
# ╠═9bc1d10b-5314-485b-86aa-3294d1c7f4c3
# ╠═9d5c901d-c3f9-4c56-bade-8a25b2be2fe7
# ╠═eaac9617-cdd5-49fb-80a2-d36b3ca15a68
# ╠═66f5865d-9236-4562-b7ce-cdc7440c4c0e
# ╠═800e1bad-8917-4d98-b0fb-73ef0469e1d6
# ╠═cacfb01e-31fc-4b0c-8b75-4a55fffdf4c8
# ╠═92ef42d9-ee38-4485-bbca-d5b4303dbb2a
# ╠═2381a54b-ba61-468f-9287-bf3e2024b694
# ╠═493fba7a-e07f-496e-b40a-1f214dd4dfb9
# ╠═115295f1-047a-4da4-a348-1b89e8c795bb
# ╟─eed56293-6b89-416c-96a8-685b5957b4f2
# ╠═d015f470-943d-42f4-ad57-c32255e8e4e6
# ╟─bdcadc88-cfda-48b3-b678-3d6e1db9effc
# ╠═61470dff-2413-4c6f-97d5-2ec739774b1c
# ╠═772339c6-6cd8-4936-be49-6b9f3a782eb8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
