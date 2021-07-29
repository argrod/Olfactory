### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 2bf25444-1ad8-44dd-be96-c7cd9e230297
using Plots

# ╔═╡ fec68519-f89a-4c66-a722-b4d709e3ff08
pos1 = (0,0)

# ╔═╡ ec3ca6a4-1b42-407c-a3b5-d84434d16f46
# add a random addition or subtraction
step1() = rand( (-1, +1) )

# ╔═╡ bd41a18f-4a4d-4c86-aaab-6bec8a38c9c6
# random walk function
function walk1D(N)
	x = 0
	xs = [x]
	
	for i in 1:N
		x += step1() # x now equals the result of the RHS
		push!(xs,x) # add the x value to the end of the xs array
	end
	
	return xs
end

# ╔═╡ 99bc319d-c16b-4175-a749-e436a4099fff
xs = walk1D(100)

# ╔═╡ 47c586f9-67d4-4439-93b8-a9a60a8d99fe
begin
	plot()
	
	for i in 1:100
		plot!(walk1D(1000), leg=false, size=(500, 300), lw=2, alpha=0.5)
	end
	
	xlabel!("t")
	ylabel!("x position in space")
	plot!()
end

# ╔═╡ 35ba3a46-629d-47ff-a973-985a4074ea49
md"""
Adding in considerations for other numbers of dimensions
"""

# ╔═╡ 1f6acc1d-1094-4d02-9771-970d3385c3ca
function walk(initialise, step, N)
	x = initialise()
	xs = [x]
	
	for i in 1:N
		x += step()
		push!(xs, x)
	end
	
	return xs
end

# ╔═╡ 530a0d4a-9841-49bb-b70a-e0b3563f9692
begin
	initialise() = 0
	step() = rand( (-1, +1) )
end

# ╔═╡ be018c73-c7d3-46dc-b012-1ccef02aa837
directions = [ [1, 0], [0, 1], [-1, 0], [0, -1] ]

# ╔═╡ 53d85407-6a3a-4056-b3b3-a48d2eeceedf
begin
	initialise2D() = [0,0]
	step2D() = rand( directions )
end

# ╔═╡ cebc84f9-f427-4bbc-bbbe-cbef9b2f8713
walk(initialise2D,step2D,10)

# ╔═╡ ad9fb1b0-bf41-4e12-9063-93a516a83db7
abstract type Walker end

# ╔═╡ b4507aa4-6141-4669-901d-a7dfa851977b
begin
	struct Walker1D <: Walker
		pos::Int64
	end
	
	Walker1D() = Walker1D(0)
end

# ╔═╡ 988584da-181c-48bd-98da-b562a43c0015
position(w::Walker) = w.pos

# ╔═╡ 1f2f1768-0755-4f3d-91af-f6e4bed03c9c
step(w::Walker1D) = rand( (-1, -1) )

# ╔═╡ 0fc50e42-a41c-4d16-bb06-2cad71864251
w = Walker1D()

# ╔═╡ 52a204f7-3770-430f-a94f-23c6cd81ccf3
md"""
How to update w in this case?
"""

# ╔═╡ a6a3ea2e-6b80-45b9-a580-3b7cd996b0b9
w.pos = -1

# ╔═╡ 67a001fc-0fbc-4444-8ea0-4e36a5e624e6
md"""
The **struct** definition generates an immutable value, one that cannot be changed internally.
"""

# ╔═╡ 8c718100-e03f-4324-85aa-f2619bfbc4d5
mutable struct MutableWalker <: Walker
	pos::Int64
end

# ╔═╡ 6076bb5f-2c3f-46a6-a7ea-16842094380b
w2 = MutableWalker(0)

# ╔═╡ b2da1d18-b1b7-49ca-baa4-ced2270a3dd3
w2.pos = -1

# ╔═╡ 0df37a43-2ebd-41ef-8383-3dadef3240e8
w2

# ╔═╡ 001013af-19c4-4cdc-8b37-b721359ab275
md"""
Update the immutable walker - immutable objects may not need storing and so can be more efficient
"""

# ╔═╡ 9241a807-54cd-4efb-8ce6-02f6398f3278
w3 = Walker1D(1)

# ╔═╡ ec6f4592-716f-47ee-83f1-82941315a619
md"""
Now need a new function to **update** the immutable walker
"""

# ╔═╡ 519c3ba0-c4c8-47d1-961d-650d925f4748
md"""
So this creates a new object of type `W` with a new value `(position(w) + step)`. To do so, the function requires access to the type `W`. `where` allows this. `::` annotation refers to type, and `<:` requires the LHS to be of subtype RHS, in this case, `Walker`. So `W` is either the `Walker1D` or `Walker2D` type that we defined above.
"""

# ╔═╡ b4ab3dce-c8fc-4971-9f37-2370c72dc365
w

# ╔═╡ a9199172-c412-4dfd-943f-ed51f00faf3d
w # w is NOT being changed, a new walker is being produced

# ╔═╡ 8ae13b81-4bae-4ec3-8cda-8d06f3af86f1
struct Walker2D <: Walker
	x::Int
	y::Int
end

# ╔═╡ 9b55e449-8798-4df6-bbb6-a59f9191e4b8
position(w::Walker2D) = (w.x, w.y)

# ╔═╡ 195dc3d5-8143-474b-8155-74402b1bb5d2
update(w::W, step) where {W <: Walker} = W(position(w) + step)

# ╔═╡ d81f8213-fe0e-4fd6-8086-39bc9d2f130c
step(w::Walker2D) = rand( [ [1, 0], [0, 1], [-1 0], [0, -1] ] )

# ╔═╡ 1920a20d-9e16-4bcc-8707-6d8be09d1379
step(w)

# ╔═╡ d5069791-409d-4ac7-b7fd-f07544110499
update(w::Walker2D, step::Vector) = Walker2D(w.x + step[1], w.y + step[2])

# ╔═╡ e062bd31-4182-4ff1-96e3-113f7ad6272f
update(w, - 1)

# ╔═╡ 77ccb79f-80fb-4a85-a41a-4153eb0d83d5
w2D = Walker2D(0, 0)

# ╔═╡ cbe38e97-1589-476e-9ef4-cb20c186525d
position(w2D)

# ╔═╡ 7b4341d4-e26a-43db-aaf6-a034403a1439
update(w2D, [-1,-1])  

# ╔═╡ 96fcf139-866f-4edc-919e-14bfebd9ecc4
position(w2D)

# ╔═╡ e949175d-245f-4f27-a107-b8ef3aab7778
function trajectory(w::W, N) where {W <: Walker}
	ws = [position(w)]
	
	for i in 1:N
		pos = position(w)
		w = update(w, step(w))
		
		push!(ws, position(w))
	end
	
	return ws
end

# ╔═╡ 240b7fb7-c5eb-47e1-877a-4da04e991560
trajectory(Walker1D(0), 10)

# ╔═╡ 2eaf1c65-7391-4b8d-99ee-e54b0831434e
trajectory(Walker2D(0, 0), 10)

# ╔═╡ 8dd48f13-fe15-4b9e-b5e1-f1cb26652119
step([-1,-1])

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Plots = "~1.19.0"
"""

# ╔═╡ Cell order:
# ╠═2bf25444-1ad8-44dd-be96-c7cd9e230297
# ╠═fec68519-f89a-4c66-a722-b4d709e3ff08
# ╠═ec3ca6a4-1b42-407c-a3b5-d84434d16f46
# ╠═bd41a18f-4a4d-4c86-aaab-6bec8a38c9c6
# ╠═99bc319d-c16b-4175-a749-e436a4099fff
# ╠═47c586f9-67d4-4439-93b8-a9a60a8d99fe
# ╠═35ba3a46-629d-47ff-a973-985a4074ea49
# ╠═1f6acc1d-1094-4d02-9771-970d3385c3ca
# ╠═cebc84f9-f427-4bbc-bbbe-cbef9b2f8713
# ╠═530a0d4a-9841-49bb-b70a-e0b3563f9692
# ╠═be018c73-c7d3-46dc-b012-1ccef02aa837
# ╠═53d85407-6a3a-4056-b3b3-a48d2eeceedf
# ╠═ad9fb1b0-bf41-4e12-9063-93a516a83db7
# ╠═b4507aa4-6141-4669-901d-a7dfa851977b
# ╠═988584da-181c-48bd-98da-b562a43c0015
# ╠═1f2f1768-0755-4f3d-91af-f6e4bed03c9c
# ╠═0fc50e42-a41c-4d16-bb06-2cad71864251
# ╠═1920a20d-9e16-4bcc-8707-6d8be09d1379
# ╠═52a204f7-3770-430f-a94f-23c6cd81ccf3
# ╠═a6a3ea2e-6b80-45b9-a580-3b7cd996b0b9
# ╠═67a001fc-0fbc-4444-8ea0-4e36a5e624e6
# ╠═8c718100-e03f-4324-85aa-f2619bfbc4d5
# ╠═6076bb5f-2c3f-46a6-a7ea-16842094380b
# ╠═b2da1d18-b1b7-49ca-baa4-ced2270a3dd3
# ╠═0df37a43-2ebd-41ef-8383-3dadef3240e8
# ╟─001013af-19c4-4cdc-8b37-b721359ab275
# ╠═9241a807-54cd-4efb-8ce6-02f6398f3278
# ╟─ec6f4592-716f-47ee-83f1-82941315a619
# ╠═195dc3d5-8143-474b-8155-74402b1bb5d2
# ╟─519c3ba0-c4c8-47d1-961d-650d925f4748
# ╠═b4ab3dce-c8fc-4971-9f37-2370c72dc365
# ╠═e062bd31-4182-4ff1-96e3-113f7ad6272f
# ╠═a9199172-c412-4dfd-943f-ed51f00faf3d
# ╠═8ae13b81-4bae-4ec3-8cda-8d06f3af86f1
# ╠═9b55e449-8798-4df6-bbb6-a59f9191e4b8
# ╠═d81f8213-fe0e-4fd6-8086-39bc9d2f130c
# ╠═d5069791-409d-4ac7-b7fd-f07544110499
# ╠═77ccb79f-80fb-4a85-a41a-4153eb0d83d5
# ╠═cbe38e97-1589-476e-9ef4-cb20c186525d
# ╠═7b4341d4-e26a-43db-aaf6-a034403a1439
# ╠═96fcf139-866f-4edc-919e-14bfebd9ecc4
# ╠═e949175d-245f-4f27-a107-b8ef3aab7778
# ╠═240b7fb7-c5eb-47e1-877a-4da04e991560
# ╠═2eaf1c65-7391-4b8d-99ee-e54b0831434e
# ╠═8dd48f13-fe15-4b9e-b5e1-f1cb26652119
# ╟─00000000-0000-0000-0000-000000000001
