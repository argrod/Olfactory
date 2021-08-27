### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 6732848d-ce1c-4f64-8607-263124cc7745
begin
	using Pkg
	Pkg.add(url="https://github.com/anowacki/CircStats.jl")
	Pkg.add("SpecialFunctions")
	using CircStats
	using SpecialFunctions
	Pkg.add("Optim")
	using Optim
end

# ╔═╡ 6eb67569-9ce7-470a-94a4-26e092926cf4
begin
	function Likelihoodww(data1,data2,cv)
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
		return function(par)
	end
end

# ╔═╡ Cell order:
# ╠═6732848d-ce1c-4f64-8607-263124cc7745
# ╠═6eb67569-9ce7-470a-94a4-26e092926cf4
