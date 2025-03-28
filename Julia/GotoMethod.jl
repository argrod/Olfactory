### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 6732848d-ce1c-4f64-8607-263124cc7745
begin
	using SpecialFunctions
	using Optim
	using Glob
	using CSV
	using DataFrames
	using Dates
	using Statistics
end

# ╔═╡ 6eb67569-9ce7-470a-94a4-26e092926cf4
begin
	function Likelihoodww(data1,data2,cv)
		return function(par)
			a = par[1]
			b = cv/gamma(1+1/a)
			mx = par[2]
			my = par[3]
			wx = par[4]
			xy = par[5]
			L = 0
			for i = 1:length(data1)
				rr = sqrt((data1[i]*cos(data2[i]) - wx)^2 + (data1[i] + sin(data2[i]) - wy)^2)
				rx = (data1[i]*cos(data2[i]) - wx)/rr
				ry = (data1[i]*sin(data2[i]) - wy)/rr
				lp = (a - 2)*log(rr) - (rr/b)^a + mx + rx + my + ry + log(a) - log(b) + (1 - a)*log(b) - log(besseli(sqrt(mx^2 + my^2),0))
				L = L + lp
			end
		end
	end
end

# ╔═╡ 2547ed04-0522-47e5-9a2a-80af02772a82
begin
	function weibullSD(a,b)
		b * sqrt(gamma(1+2/a) - gamma(1+1/a) * gamma(1+1/a))
	end
end	

# ╔═╡ 1dd1bf6b-2ae1-4749-923e-d6dd09e3c5b2
begin
	function weibull_mean(a,b)
		b+gamma(1+1/a)
	end
end

# ╔═╡ b1baff6d-69b6-46d0-9d8d-52a8e578f266
begin
	function vM_sd(kappa)
		1/sqrt(kappa)
	end
end

# ╔═╡ e220fdcd-74e5-4653-89bc-5cfbd58b163d
begin
	
	# read in
	if Sys.iswindows()
		infile = "E:/My Drive/PhD/Data/2018Shearwater/AxyTrek/"
		outfile = "E:/My Drive/PhD/Data/2018Shearwater/WindEst/MinuteGotoVersion/"
	else
		outfile = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/WindEst/MinuteGotoVersion/"
		infile = "/Volumes/GoogleDrive/My Drive/PhD/Data/2018Shearwater/AxyTrek/"
	end
	files = glob("*/*.txt",infile)
	# EDat = DataFrame()
	# for b = 1:length(files)
	# 	append!(EDat, hcat(CSV.File(files[b],header=false,delim="\t") |> DataFrame,fill(readdir(fileloc19)[b],nrow(CSV.File(files[b],header=false,delim="\t") |> DataFrame))))
	# end
	# rename!(EDat,[:DT,:lat,:lon,:c4,:c5,:c6,:c7,:c8,:ID])
	# dt = dateformat"y/m/d,H:M:S"
	# EDat.DT = DateTime.(EDat.DT,dt);
end

# ╔═╡ 6e34b25d-030c-4095-bacb-2cff4ae66caa
begin
	md"""
	### Defining variables
	The sampling interval, window for analysis, minimum number of points, minimum ground speed, and sampling interval error all need defining. They are as follows:
	
	| Variable | Value |
	|----------|-------|
	| Sampling Interval | 60 seconds |
	| Analysis window | 51 minutes |
	| Minimum points | 45 |
	| Minimum ground speed | 4 m$s^{-1}$ |
	| Sampling interval error | 5 |
	
	"""
end

# ╔═╡ ec89669f-b61f-4c5e-9cba-81b77ca465ad
begin
	sampling_interval = 60;     #  sampling interval [sec]
	time_window = 51;     #  time length of time window [min] *Give an odd number!
	cutlength = 45;     # minimum number of data points (track vectors) included in a time window [points]         
	cutv = 4.1667; # minimum ground speed [m/sec]

	###         Condition 3: give mean air speed value      #####
	constv = 34.7/3.6; # mean air speed [m/sec]
	#We gave the mean ground speed of streaked shearwater (Shiomi #et. al. 2012) as the mean air speed
	error_of_sampling_interval = 5;
	cutt = sampling_interval+error_of_sampling_interval; 
	# upper value of sampling interval [sec]
	windwidth = time_window-1;                           
	# length of time window(velocity)  [min]
end

# ╔═╡ 8777e91b-900b-457a-9d43-188f0f0cc70b
begin
	function findInGap(DT,time,gap)
		(findall(x -> (time - Dates.Millisecond(round(gap/2))) <= x <= (time + Dates.Millisecond(round(gap/2))),DT))
	end
end

# ╔═╡ f2908902-c0fe-4efb-8d99-59196927f140
begin
	function subSamp(DT::Vector{Dates.DateTime},gap::Int64=5)
		newT = (DT[1]:Dates.Minute(gap):DT[end])
		Scheck = median(Dates.value.(Dates.Millisecond.(diff(DT))))
		ind = findInGap.((DT,),newT,(Scheck,))
		# for b in 0:(length(newT)-1)
		# 	try
		# 		push!(findInGap(DT,newT[b+1]),ind)
		# 	catch
		# 		continue
		# 	end
		# end
		# filter!(!isempty,ind)
		return(ind)
	end
end

# ╔═╡ 6fe19134-6d35-4ddf-bc90-858c05027056
a = 18

# ╔═╡ 6ef4a57e-1f96-43d0-b21f-19a814405bcd
begin
	b = 1
	EDat = CSV.File(files[b],header=false,delim="\t") |> DataFrame
	rename!(EDat,[:DT,:lat,:lon,:altitude,:grSpd,:satCount,:hdop,:maxSigStr])
	df = DateFormat("d/m/y,H:M:S")
	EDat.DT = DateTime.(EDat.DT,df)
	tsel = (EDat.DT[1]:Dates.Minute(1):EDat.DT[end])
	select = Int[];
	AxDat = EDat[subSamp(EDat.DT,5),:]
	# for b = 1:length(tsel)
	# 	choose = which(EDat.DT .>= (tsel[b] - Dates.Second(dtFull)) & EDat.DT .<= (tsel[b] + (median(dtFull)/2)))
	# end
end

# ╔═╡ d264811a-3eda-47d4-aac4-b4e5b2456daf
newT = (EDat.DT[1]:Dates.Minute(5):EDat.DT[end])

# ╔═╡ 3b02c1bd-55c5-4d69-933a-557470c98b77
min(findall(x -> (newT[a] - Dates.Millisecond(round(5000/2))) <= x <= (newT[a] + Dates.Millisecond(round(5000/2))),EDat.DT)...)

# ╔═╡ c0763d32-5da6-4c86-acd3-2902b554513d
findInGap.((EDat.DT,),newT[a:a+10],(5000,))

# ╔═╡ 5f4a208c-0096-41f8-be17-64ffb9c86e9e
(findall(x -> (newT[a] - Dates.Millisecond(round(5000/2))) <= x <= (newT[a] + Dates.Millisecond(round(5000/2))),EDat.DT))

# ╔═╡ d717299f-41f5-49fe-8445-fb7cfb89502f
subSamp(EDat.DT,5)

# ╔═╡ dcb64cc3-540d-41ec-8b88-00299df45e83
typeof(subSamp(EDat.DT,5))

# ╔═╡ f7cc9490-b543-423c-bde9-7460dc4aabb9


# ╔═╡ cadedcef-9cb2-4acd-9662-d09b7d6bc69a


# ╔═╡ ad2582b6-cc8c-44ea-b749-10b60aaac586
  AxDat <- read.delim(paste(infile,files[a], sep = ''), sep = '\t', header = F)
    AxDat[, 1] <- as.POSIXct(as.character(AxDat[, 1]), format = '%d/%m/%Y,%H:%M:%OS') + (9*3600)
    tsel <- seq(AxDat[1, 1], AxDat[nrow(AxDat), 1], by = 60)
    dtFull <- as.numeric(diff(AxDat[, 1]))
    # find data that line up to the new timepoints
    select <- NA
    for(b in 1:length(tsel)){
        choose <- which(AxDat[, 1] >= (tsel[b] - (median(dtFull)/2)) & AxDat[, 1] <= (tsel[b] + (median(dtFull)/2)))
        if(length(choose) != 0){
            if(length(choose) > 1){
                select[b] <- choose[which.min(abs(tsel[b] - AxDat[choose, 1]))]
            } else {
                select[b] <- choose
            }
        }
    }
    AxDat <- AxDat[na.omit(select),]
    # remove superfluous values
    AxDat <- AxDat[, 1:3]
    colnames(AxDat) <- c("DT", "lat", "lon")
    AxDat$DT <- as.POSIXct(as.character(AxDat$DT), format = '%Y-%m-%d %H:%M:%OS') + (9*3600)
    dt <- as.numeric(difftime(AxDat$DT[2:nrow(AxDat)], AxDat$DT[1:(nrow(AxDat) - 1)]), units = "secs")
    lat <- AxDat[,2]
	long <- AxDat[,3]
	cord.dec <- SpatialPoints(cbind(long, lat), proj4string=CRS("+proj=longlat"))
	cord.UTM <- spTransform(cord.dec, CRS("+proj=utm +zone=54 +datum=WGS84"))
	X <- coordinates(cord.UTM)[,1]
	Y <- coordinates(cord.UTM)[,2]
    n <- nrow(AxDat) - 1
    vg_x_obs <- X[2:(n+1)] - X[1:n]
    vg_y_obs <- Y[2:(n+1)] - Y[1:n]
    g_speed_obs <- sqrt(vg_x_obs^2 + vg_y_obs^2)/dt
    g_direction_obs <- atan2(vg_y_obs, vg_x_obs)
    tp <- AxDat$DT
    rrow <- g_speed_obs
    drow <- g_direction_obs

    startpoint<-floor(windwidth/2)
	endpoint<-length(rrow)-floor(windwidth/2)


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.8.5"
DataFrames = "~1.2.2"
Glob = "~1.3.0"
Optim = "~1.4.1"
SpecialFunctions = "~1.6.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b6dec2ed4f10840e2cf836508525656450d4d289"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.26"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CSV]]
deps = ["Dates", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode"]
git-tree-sha1 = "b83aa3f513be680454437a0eee21001607e5d983"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.8.5"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bdc0937269321858ab2a4f288486cb258b9a0af7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "727e463cfebd0c7b999bbf3e9e7e16f254b94193"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.34.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

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
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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
git-tree-sha1 = "3ed8fa7178a10d1cd0f1ca524f249ba6937490c0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7c365bdef6380b29cfc5caaf99688cd7489f9b87"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.2"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

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

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InvertedIndices]]
deps = ["Test"]
git-tree-sha1 = "15732c475062348b0165684ffe28e85ea8396afc"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.0.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

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

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3d682c07e6dd250ed082f883dc88aee7996bf2cc"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.0"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "2ca267b08821e86c5ef4376cffed98a46c2cb205"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.1"

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

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Optim]]
deps = ["Compat", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7863df65dbb2a0fa8f85fcaf0a41167640d2ebed"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.4.1"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "bfd7d8c7fd87f04543810d9cbd3995972236ba1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "1.1.2"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "cde4ce9d6f33219465b55162811d8de8139c0414"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.2.1"

[[PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "54f37736d8934a12a200edea2f9206b03bdf3159"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.7"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

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

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
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

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═6732848d-ce1c-4f64-8607-263124cc7745
# ╠═6eb67569-9ce7-470a-94a4-26e092926cf4
# ╠═2547ed04-0522-47e5-9a2a-80af02772a82
# ╠═1dd1bf6b-2ae1-4749-923e-d6dd09e3c5b2
# ╠═b1baff6d-69b6-46d0-9d8d-52a8e578f266
# ╠═e220fdcd-74e5-4653-89bc-5cfbd58b163d
# ╠═6e34b25d-030c-4095-bacb-2cff4ae66caa
# ╟─ec89669f-b61f-4c5e-9cba-81b77ca465ad
# ╠═8777e91b-900b-457a-9d43-188f0f0cc70b
# ╠═f2908902-c0fe-4efb-8d99-59196927f140
# ╠═d264811a-3eda-47d4-aac4-b4e5b2456daf
# ╠═6fe19134-6d35-4ddf-bc90-858c05027056
# ╠═3b02c1bd-55c5-4d69-933a-557470c98b77
# ╠═c0763d32-5da6-4c86-acd3-2902b554513d
# ╠═5f4a208c-0096-41f8-be17-64ffb9c86e9e
# ╠═d717299f-41f5-49fe-8445-fb7cfb89502f
# ╠═dcb64cc3-540d-41ec-8b88-00299df45e83
# ╠═6ef4a57e-1f96-43d0-b21f-19a814405bcd
# ╠═f7cc9490-b543-423c-bde9-7460dc4aabb9
# ╠═cadedcef-9cb2-4acd-9662-d09b7d6bc69a
# ╠═ad2582b6-cc8c-44ea-b749-10b60aaac586
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
