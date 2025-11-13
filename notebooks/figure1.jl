### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 2e84ad80-cffa-407f-8817-091f6e821d04
using Gnuplot

# ╔═╡ 84a77c02-6ea5-11ef-0ccc-29525292ecad
# Cahill, K. (2012). Models of membrane electrostatics. Physical Review E, 85(5), 051921. https://doi.org/10.1103/PhysRevE.85.051921

# we start by trying to reproduce some of the figures from this paper

# ╔═╡ af0d5461-a81f-4676-8ed1-7e9fbfcb3377
begin
	# Physical constants; done dirt cheap
	"Permittivity of free space, (C²N⁻¹m⁻²)."
	const ε_0 = ϵ_0 = 8.85418682e-12
	"Electron charge, (kgm²s⁻²)."
	const eV = q = ElectronVolt = 1.602176634e-19
	
end

# ╔═╡ c7de94fd-8792-4c7b-987f-1dedf9aa3fa8
# "Water absolute dielectric."
const ϵ_w = 80 * ϵ_0

# ╔═╡ 1a32993f-4aa1-4b96-90ec-e6aac423a2bb
# "Lipid absolute dielectric."
const ϵ_l = 2 * ϵ_0

# ╔═╡ f4da6d89-4a66-47b9-a05b-95b44ac45eb2
# "Cytosol absolutely dielectric"
const ϵ_c = 80 * ϵ_0

# ╔═╡ 630dd62a-c312-461f-a97b-b7b2b04c15a3
# nm, or at least 1E-9 scaling
const nm = 1E-9

# ╔═╡ c9b29a75-daa0-4c06-966a-4e59690abe41
p=(ϵ_w - ϵ_l)/(ϵ_w + ϵ_l)

# ╔═╡ 31a8eaba-0feb-4692-84ba-f75561db3729
p′=(ϵ_c - ϵ_l)/(ϵ_c + ϵ_l)

# ╔═╡ 12966c74-8de9-43cc-8cf9-cb816b1bfd23
ϵ_w

# ╔═╡ be91df5f-75fa-4096-a89f-995fd7c836d2
q/(4π*ϵ_0)

# ╔═╡ 589c65a3-58ea-489a-8f41-a55870b88181


# ╔═╡ 1cea7fcb-e445-45af-a0b8-18d67f009850


# ╔═╡ 589cf1f7-56a6-4dd3-a2c2-1bb16a1edbe3
# Eqn (10)
ϵ_wl=(ϵ_w + ϵ_l)/2

# ╔═╡ 42f3b835-c8bc-499d-93ce-c23e01be3748
ϵ_cl = (ϵ_c + ϵ_l)/2

# ╔═╡ 8286961c-6c21-4bc6-97aa-e584da2a81ed


# ╔═╡ da083ad4-c569-4e9e-847a-3f9c01d35906


# ╔═╡ a1ce501e-2822-42ce-9947-937fd7400bd1
t=5nm

# ╔═╡ 32037f8f-1bad-459a-afdd-312615a1eef5
Zs=collect(-(t+0.9nm):t/100:5nm)

# ╔═╡ 5630d8fe-d7ef-4b97-8bca-25cacd287fd2
h=-1nm

# ╔═╡ a3e0eb0f-0bf2-4e6c-9718-436c5e2c50e0
summand9(n,z,ρ)=p*p′^(n-1)/(√(ρ^2 + (z + 2n*t + h)^2))

# ╔═╡ ddf44a07-63d8-40b2-bc60-362eb87da3d4
# Figure 1 - Phospholipid bilayer, blue curve - charge in the water at (0,1)

# Eqn (9) - potential in water
function V_ww(z; ρ=sqrt(0.707nm^2+0.707nm^2)) 
	q/(4π*ϵ_w) * 
	(
		  1/√(ρ^2+(z-h)^2) 
		+ p/√(ρ^2+(z+h)^2)
		- p′*(1-p^2) * sum([ summand9(n,z,ρ) for n in 1:1000 ] ) 
	)

end

# ╔═╡ 6aa0eb1b-3831-478c-881e-dbdbb441938c
# Eqn (10)
summand10(n,z,ρ) = (p*p′)^n * 
( 
	1/√(ρ^2 + (z - 2n*t - h)^2)  
	- p′/√(ρ^2 + (z + 2*(n+1)*t + h)^2) )

# ╔═╡ 16c0c851-419e-4fd7-a696-6f9d8bed48fc
# Eqn (10) - lipid bilayer potential
function V_wl(z; ρ=sqrt(0.707nm^2+0.707nm^2))
	q/(4π*ϵ_wl) *
	sum([summand10(n,z,ρ) for n in 0:1000])
end

# ╔═╡ ee8a56ff-5626-4d1f-a323-7fd55b417b1e
summand11(n,z,ρ) = (p*p′)^n / √(ρ^2 + (z - 2n*t - h)^2)

# ╔═╡ bd2841c7-4cc5-43ce-9c0b-e6ca9080594e
# Eqn (11) - cytosol potential
function V_wc(z;ρ=sqrt(0.707nm^2+0.707nm^2))
	q*ϵ_l/(4π*ϵ_wl*ϵ_cl) *
	sum([summand11(n,z,ρ) for n in 0:1000])
end

# ╔═╡ e1edefc5-73cd-4bb6-8fea-d76691946c79
function V(z;ρ=sqrt(0.707nm^2+0.707nm^2))
	if z>=0 # in water z>=0
		V=V_ww(z,ρ=ρ)
	elseif z>-t # in lipid -t<z<0
		V=V_wl(z,ρ=ρ)
	else # in cytosol, z<-t
		V=V_wc(z,ρ=ρ)
	end
	V
end

# ╔═╡ 4e9fd964-70e3-4c1b-88e3-c430db7046b8


# ╔═╡ c5ccccb3-66c6-4d3e-82c6-4fa7e73c0965
length(Zs)

# ╔═╡ 07bc6449-5cb5-44f7-add9-dd13573776c1
Vs=[V(z, ρ=0.3nm) for z in Zs ].+[ V(z, ρ=1.0nm) for z in Zs ] 

# ╔═╡ c90f5dda-fadf-4d9a-9719-7156f8a40bd4
maximum(Vs)

# ╔═╡ 7d4b997a-2845-41b4-9834-62a7ae5c2062
minimum(Vs) 

# ╔═╡ e59df517-6aad-45ef-baad-baf4eae5df20
length(Vs)

# ╔═╡ e22c6183-13c1-451c-ae68-fe255cbdc7a9
begin
	@gp "set xlabel 'Distance from membrane (nm)'"
	@gp :- "set ylabel 'Potential (V)'"
	@gp :- "set arrow from 0,0 to 0,0.025 nohead lc rgb 'red'"
	@gp :- "set arrow from -5E-9,0 to -5E-9,0.025 nohead lc rgb 'red'"
	@gp :- Zs [V(z, ρ=1.0nm) for z in Zs ] "w l"
	@gp :- Zs [V(z, ρ=0.5nm) for z in Zs ] "w l"
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Gnuplot = "dc211083-a33a-5b79-959f-2ff34033469d"

[compat]
Gnuplot = "~1.6.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.4"
manifest_format = "2.0"
project_hash = "b2dea39cecb56c72baabcd1299cc6019ed5a1596"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Gnuplot]]
deps = ["ColorSchemes", "ColorTypes", "Colors", "DataStructures", "PrecompileTools", "REPL", "ReplMaker", "StatsBase", "StructC14N", "Test"]
git-tree-sha1 = "72b7242dccedbe153dadbf1e1412f9bff3d81bad"
uuid = "dc211083-a33a-5b79-959f-2ff34033469d"
version = "1.6.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.ReplMaker]]
deps = ["REPL", "Unicode"]
git-tree-sha1 = "f8bb680b97ee232c4c6591e213adc9c1e4ba0349"
uuid = "b873ce64-0db9-51f5-a568-4457d8e49576"
version = "0.2.7"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StructC14N]]
deps = ["DataStructures", "Test"]
git-tree-sha1 = "a3d153488e0fe30715835e66585532c0bcf460e9"
uuid = "d2514e9c-36c4-5b8e-97e2-51e7675c221c"
version = "0.3.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═84a77c02-6ea5-11ef-0ccc-29525292ecad
# ╠═2e84ad80-cffa-407f-8817-091f6e821d04
# ╠═af0d5461-a81f-4676-8ed1-7e9fbfcb3377
# ╠═c7de94fd-8792-4c7b-987f-1dedf9aa3fa8
# ╠═1a32993f-4aa1-4b96-90ec-e6aac423a2bb
# ╠═f4da6d89-4a66-47b9-a05b-95b44ac45eb2
# ╠═630dd62a-c312-461f-a97b-b7b2b04c15a3
# ╠═c9b29a75-daa0-4c06-966a-4e59690abe41
# ╠═31a8eaba-0feb-4692-84ba-f75561db3729
# ╠═12966c74-8de9-43cc-8cf9-cb816b1bfd23
# ╠═be91df5f-75fa-4096-a89f-995fd7c836d2
# ╠═ddf44a07-63d8-40b2-bc60-362eb87da3d4
# ╠═589c65a3-58ea-489a-8f41-a55870b88181
# ╠═1cea7fcb-e445-45af-a0b8-18d67f009850
# ╠═a3e0eb0f-0bf2-4e6c-9718-436c5e2c50e0
# ╠═589cf1f7-56a6-4dd3-a2c2-1bb16a1edbe3
# ╠═16c0c851-419e-4fd7-a696-6f9d8bed48fc
# ╠═6aa0eb1b-3831-478c-881e-dbdbb441938c
# ╠═42f3b835-c8bc-499d-93ce-c23e01be3748
# ╠═bd2841c7-4cc5-43ce-9c0b-e6ca9080594e
# ╠═ee8a56ff-5626-4d1f-a323-7fd55b417b1e
# ╠═32037f8f-1bad-459a-afdd-312615a1eef5
# ╠═8286961c-6c21-4bc6-97aa-e584da2a81ed
# ╠═da083ad4-c569-4e9e-847a-3f9c01d35906
# ╠═a1ce501e-2822-42ce-9947-937fd7400bd1
# ╠═5630d8fe-d7ef-4b97-8bca-25cacd287fd2
# ╠═e1edefc5-73cd-4bb6-8fea-d76691946c79
# ╠═4e9fd964-70e3-4c1b-88e3-c430db7046b8
# ╠═c90f5dda-fadf-4d9a-9719-7156f8a40bd4
# ╠═7d4b997a-2845-41b4-9834-62a7ae5c2062
# ╠═c5ccccb3-66c6-4d3e-82c6-4fa7e73c0965
# ╠═e59df517-6aad-45ef-baad-baf4eae5df20
# ╠═07bc6449-5cb5-44f7-add9-dd13573776c1
# ╠═e22c6183-13c1-451c-ae68-fe255cbdc7a9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
