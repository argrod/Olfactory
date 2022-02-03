### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

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
end

# ╔═╡ 24dd02f0-f348-11eb-02fe-3309b1be5acc
md"""
# Fitting circle to track vector plot

Following discussions with 米原さん, after the [PNAS publication](https://www.pnas.org/content/113/32/9039 "Flight paths of seabirds soaring over the ocean surface enable measurement of fine-scale wind speed and direction") he tested another method of wind estimation, instead fitting a circle to a vector plot of the ground speed. Deviation of the fitted circle from $(0,0)$ indicates the wind vector.
"""

# ╔═╡ 044a9589-e372-432a-bdc7-b1003bdc4aa6
Base.load.path()

# ╔═╡ 432b5aea-9520-4dcf-8605-5a4c93e2c11c


# ╔═╡ 33c61589-def4-4fb8-ab34-a683a11205e8


# ╔═╡ 0e34f851-4e25-4b6a-b76f-672b311a6819
# READ IN AND FORMAT
begin
	fileloc19 = "/Volumes/GoogleDrive/My Drive/PhD/Data/2019Shearwater/AxyTrek/"
	files = glob("*/*.txt",fileloc19)
	EDat = DataFrame()
	for b = 1:length(files)
		append!(EDat, hcat(CSV.File(files[b],header=false,delim="\t") |> DataFrame,fill(readdir(fileloc19)[b],nrow(CSV.File(files[b],header=false,delim="\t") |> DataFrame))))
	end
	rename!(EDat,[:DT,:lat,:lon,:c4,:c5,:c6,:c7,:c8,:ID])
	dt = dateformat"y/m/d,H:M:S"
	EDat.DT = DateTime.(EDat.DT,dt)
end

# ╔═╡ 19ad02e2-a725-42e4-8d3e-ecaa93a559a7
begin
	tags = unique(EDat.ID)
	sel = EDat[findall(x->x==tags[1],EDat.ID),:]
end

# ╔═╡ 6f21012b-39d6-47d4-99b3-0dde5068c357


# ╔═╡ 39232859-9a15-450a-a79e-6632ecbed677


# ╔═╡ 12e5dc37-400d-4406-8de9-26a20df98eb3


# ╔═╡ d9e65c24-6cf0-47ac-8039-7d47186613cd


# ╔═╡ 49e45c19-8e5d-4f8d-8bf8-ecaa2e65707d


# ╔═╡ 16b94084-9283-4d4b-adb8-2b08e52462a8


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
		rmv = findall(x -> x .<= 5, (tve .- tvs .+ 1))
		for x in rmv
			rest[tvs[x]:tve[x]] .= false
		end
		vs = findall(x -> x == -1, diff([0;rest]))
		ve = findall(x -> x == 1, diff([rest;0]))
		return spd, dir, vs, ve
	end
end

# ╔═╡ 44183ed0-b125-42e4-bb34-d8b60d95d272
spd,dir,vs,ve=flightMask(sel.DT,sel.lat,sel.lon,4,5)

# ╔═╡ 69dc708a-2ee8-41f9-86d6-81c88f898069
rest = spd .< 4

# ╔═╡ 506ff3d6-bd06-4045-a418-e1888b1afee5
diff(rest)

# ╔═╡ 5f2b2319-ee04-46f9-8b18-390435b30504


# ╔═╡ 7c23b6d9-e659-430f-bc5e-ca24ef58196d
begin
	# separate flight data into section for wind estimation
	function getsection(F,Win,delta,fs,fe)
		# F - sampling frequency (Hz)
		# Win - sampling window (seconds)
		# delta - window step length (seconds)
		# fs and fe - output from flightMask
		window = Win * F
		del = delta * F
		
		# use only segments that exceed 7 minutes in length
		fsa = fs[(fe .- fs) .>= (7 * 60 * F)] + (60 * F)
		fea = fe[(fe .- fs) .>= (7 * 60 * F)] - (60 * F)
		
		secnum = sum(floor((fea .- fsa .- window)./del) + 1)
		
		ss = zeros(Int8,secnum,1); se = zeros(Int8,secnum,1)
		k = 1
		while k <= secnum
			for i = 1:length(fsa)
				for c = 0:floor((fea[i] - fsa[i] - window)./del)
					ss[k] = fsa[i] + c * del;
					se[k] = fsa[i] + c * del + window
					k = k + 1
				end
			end
		end
		return ss,se
	end
end

# ╔═╡ 85b28651-7cb0-4d62-8fb9-6926bbc0fc0d
begin
	window = 300 * .2
	del = 60*.2
	fsa = vs[(ve .- vs) .>= (7 * 60 * .2)] .+ (60 * .2)
	fea = ve[(ve .- vs) .>= (7 * 60 * .2)] .- (60 * .2)
end

# ╔═╡ 95e4ab43-3be8-47db-9cf7-44590a99ba99


# ╔═╡ cccf3dee-0ed2-4a1f-a828-2ab04331361a


# ╔═╡ 4b0b2af9-0760-477d-9394-12728641364d
sum((ve .- vs) .>= (7 * 60 * .2))

# ╔═╡ 49dbe98a-9f6c-400a-aa3f-1aa7f641d55d


# ╔═╡ 4c2b26a5-1b58-44b9-b4d8-347d44308ceb
vs[1]

# ╔═╡ ce81b81b-ad03-452c-8a46-30a21975dcce
ve[1]

# ╔═╡ 863e2c7d-b4d4-42c1-bd97-4d3b6b7f7c8f


# ╔═╡ 5b19e0ac-4318-4061-8344-03310a01a6ed


# ╔═╡ 3174b784-1bde-4fb2-b014-4ad049efa5bc


# ╔═╡ 704285fc-692d-43e1-a3dc-60f5a524991e


# ╔═╡ d09e7671-8bbe-4bd1-b1b5-aaddfcfe9cc8
(ve .- vs) .>= (7*60*9.2)

# ╔═╡ 1d0f649d-9de1-402e-8f44-5896c5b1d3dc


# ╔═╡ ec65b7ef-3016-4cb7-a414-b498bdfdf18f


# ╔═╡ 2acd5143-ad20-4049-964b-2c968e083740
length(vs)

# ╔═╡ f5ff84ff-8930-4dd8-92fe-ea72b62e5176


# ╔═╡ 7de594be-d8f0-412d-970a-380f3fdec814


# ╔═╡ 361d07d3-b7be-4eb0-bfe5-6ed46f94b063


# ╔═╡ effc3ba5-4dd7-4848-a8dd-c5f299c5e718
window = Win * F
		del = delta * F
		
		# use only segments that exceed 7 minutes in length
		fsa = fs[(fe .- fs) >= (7 * 60 * F)] + (60 * F)
		fea = fe[(fe .- fs) >= (7 * 60 * F)] - (60 * F)
		
		secnum = sum(floor((fea .- fsa .- window)./del) + 1)
		
		ss = zeros(Int8,secnum,1); se = zeros(Int8,secnum,1)
		k = 1
		while k <= secnum
			for i = 1:length(fsa)
				for c = 0:floor((fea[i] - fsa[i] - window)./del)
					ss[k] = fsa[i] + c * del;
					se[k] = fsa[i] + c * del + window
					k = k + 1
				end
			end
		end

# ╔═╡ a72ec19d-3c6f-43ba-9b62-8b7eaa5aadfd
ss,se = getsection(.2,300,60,vs,ve)

# ╔═╡ 2e531133-6ed9-4038-b5f7-7059c3198742


# ╔═╡ cabe59d8-b22e-4b5b-ad6e-7c01f709c4ea
begin
	function wind2dveclsq(vg,gd,c0)
		wx(v) = v[1] .* cos(v[2])
		wy(v) = v[1] .* sin(v[2])
		ea(c) = wv([c[1] c[2]]).*sin(gd) + wx([c[1] c[2]]).*cos(gd) + sqrt((wy([c[1] c[2]]).*sin(gd) + wx([c[1] c[2]]).*cos(gd)).^2 - wy([c[1] c[2]]).^2 - wx([c[1] c[2]]).^2 + c[3]^2) - vg
		lb = [0 -pi 0]
		ub = [20 pi 20]
		ydata = ea(c0)
		curve_fit(ea,c0,lower = lb, upper = ub)
	end
end

# ╔═╡ 11cdeaca-f0c1-47d2-9f4d-1ca0bbbc9d5e
foist(v) = v[1] .* cos(v[2])
		secnd(v) = v[1] .* sin(v[2])
		full(c) = secnd([c[1] c[2]]).*sin(gd) + foist([c[1] c[2]]).*cos(gd) + sqrt((secnd([c[1] c[2]]).*sin(gd) + foist([c[1] c[2]]).*cos(gd)).^2 - secnd([c[1] c[2]]).^2 - foist([c[1] c[2]]).^2 + c[3]^2) - vg
		lb = [0 -pi 0]
		ub = [20 pi 20]
		ydata = ea(c0)

# ╔═╡ b540a82e-95d1-47d4-bb8a-56196282eb73


# ╔═╡ fa8a55ad-1721-4bb2-a6fb-16bc81eaa87c


# ╔═╡ c6658800-236d-485a-a743-a5adc0266ef1


# ╔═╡ 946d9015-c75c-47f9-ae25-c09d002e6eda
begin
	m(t, p) = p[1] * exp.(p[2] * t)
	p0 = [0.5, 0.5]
	fit = curve_fit(m, tdata, ydata, p0)
end

# ╔═╡ Cell order:
# ╟─24dd02f0-f348-11eb-02fe-3309b1be5acc
# ╠═b2b504be-f728-42d5-bb28-04c7a566cf5d
# ╠═53a548a9-5e87-47bd-b5ee-f9acdb96ae5a
# ╠═044a9589-e372-432a-bdc7-b1003bdc4aa6
# ╠═432b5aea-9520-4dcf-8605-5a4c93e2c11c
# ╠═33c61589-def4-4fb8-ab34-a683a11205e8
# ╠═0e34f851-4e25-4b6a-b76f-672b311a6819
# ╠═19ad02e2-a725-42e4-8d3e-ecaa93a559a7
# ╠═6f21012b-39d6-47d4-99b3-0dde5068c357
# ╠═39232859-9a15-450a-a79e-6632ecbed677
# ╠═12e5dc37-400d-4406-8de9-26a20df98eb3
# ╠═d9e65c24-6cf0-47ac-8039-7d47186613cd
# ╠═49e45c19-8e5d-4f8d-8bf8-ecaa2e65707d
# ╠═16b94084-9283-4d4b-adb8-2b08e52462a8
# ╠═62ee1055-dbd9-4f93-990d-ce40092e2b1b
# ╠═44183ed0-b125-42e4-bb34-d8b60d95d272
# ╠═69dc708a-2ee8-41f9-86d6-81c88f898069
# ╠═506ff3d6-bd06-4045-a418-e1888b1afee5
# ╠═5f2b2319-ee04-46f9-8b18-390435b30504
# ╠═7c23b6d9-e659-430f-bc5e-ca24ef58196d
# ╠═85b28651-7cb0-4d62-8fb9-6926bbc0fc0d
# ╠═95e4ab43-3be8-47db-9cf7-44590a99ba99
# ╠═cccf3dee-0ed2-4a1f-a828-2ab04331361a
# ╠═4b0b2af9-0760-477d-9394-12728641364d
# ╠═49dbe98a-9f6c-400a-aa3f-1aa7f641d55d
# ╠═4c2b26a5-1b58-44b9-b4d8-347d44308ceb
# ╠═ce81b81b-ad03-452c-8a46-30a21975dcce
# ╠═863e2c7d-b4d4-42c1-bd97-4d3b6b7f7c8f
# ╠═5b19e0ac-4318-4061-8344-03310a01a6ed
# ╠═3174b784-1bde-4fb2-b014-4ad049efa5bc
# ╠═704285fc-692d-43e1-a3dc-60f5a524991e
# ╠═d09e7671-8bbe-4bd1-b1b5-aaddfcfe9cc8
# ╠═1d0f649d-9de1-402e-8f44-5896c5b1d3dc
# ╠═ec65b7ef-3016-4cb7-a414-b498bdfdf18f
# ╠═2acd5143-ad20-4049-964b-2c968e083740
# ╠═f5ff84ff-8930-4dd8-92fe-ea72b62e5176
# ╠═7de594be-d8f0-412d-970a-380f3fdec814
# ╠═361d07d3-b7be-4eb0-bfe5-6ed46f94b063
# ╠═effc3ba5-4dd7-4848-a8dd-c5f299c5e718
# ╠═a72ec19d-3c6f-43ba-9b62-8b7eaa5aadfd
# ╠═2e531133-6ed9-4038-b5f7-7059c3198742
# ╠═cabe59d8-b22e-4b5b-ad6e-7c01f709c4ea
# ╠═11cdeaca-f0c1-47d2-9f4d-1ca0bbbc9d5e
# ╠═b540a82e-95d1-47d4-bb8a-56196282eb73
# ╠═fa8a55ad-1721-4bb2-a6fb-16bc81eaa87c
# ╠═c6658800-236d-485a-a743-a5adc0266ef1
# ╠═946d9015-c75c-47f9-ae25-c09d002e6eda
