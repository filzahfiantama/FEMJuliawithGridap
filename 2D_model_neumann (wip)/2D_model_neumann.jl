### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d73852a0-3b80-11ef-2ed1-f108b2d1ceeb
begin
	using PlutoUI, Gridap, GridapGmsh
	import Gmsh:gmsh
	import Gridap:solve
end

# ╔═╡ 70cfc41c-9845-4d83-96d2-dd31b96e8983
md"""
### Importing necessary packages
"""

# ╔═╡ d8be6042-cb6c-473a-8193-0845c334601e
TableOfContents()

# ╔═╡ 2e1d62ba-0de3-46b3-a680-b15f4c27d9ca
md"""
# 2D model of a Vertical Section of the Subsurface
"""

# ╔═╡ 7f996195-6a95-4556-a9f7-b573a97964b7
md"""
This notebook touches on a more complex scenario, where the model doesn't only have one boundary condition accross the boundary. We will go into further details as the notebook goes on.
"""

# ╔═╡ da73ce21-1023-4309-becc-b5c58cbfdff3
md"""
Following assumptions are made for this model to simplify but still accounting for boundary conditions in reality:
- The subsurface consists of only one rock type and is assigned one thermal conductivity value
- A constant temperature at the top of the model
- A constant heat flow at the bottom of the model
"""

# ╔═╡ 432fabbb-9d9f-46f5-a497-693f247af2eb
md"""
## Defining the Domain with a Mesh
"""

# ╔═╡ e9caa95e-ec8c-411f-82db-1d6ab163a690
md"""
The following file name will be used across the mesh generation processes
"""

# ╔═╡ 4b7ee5b5-2edc-414b-b695-5b2c009565b4
Filename = "2dmodel"

# ╔═╡ 46ecf0b8-06bf-47ce-b661-74a6a7a6b12a
md"""
The geometry basis for the mesh will be generated with the `gmsh` method from GMSH package as a `.geo` file, on top of which the mesh will be generated with the help of `GridapGmsh` package. At the end we will have a `.msh` file, which will be visualized and exported as an image with Paraview.
"""

# ╔═╡ c1469e5f-385b-470a-8898-8f4a8594442b
md"""
### Generating the Geometry
"""

# ╔═╡ 72da26fd-35f3-4957-a43a-a0139bf2616b
md"""
### Meshing on the Geometry
"""

# ╔═╡ 254b17e4-dc7d-4bf6-b71e-b7fad3f5990d
md"""
### Reading and Displaying the Mesh Data
"""

# ╔═╡ 30733946-1cee-4ff0-911b-6d2eb943a134
mesh = GmshDiscreteModel(Filename*".msh")

# ╔═╡ 7cfea374-92e4-4c43-bfbc-6ad9d1a8b66c
writevtk(mesh,Filename)

# ╔═╡ 994c101c-98ec-478a-9b72-7389483d142f
md"""
Here's the resulting mesh image exported from Paraview
$((PlutoUI.LocalResource("2dmodel.png")))
"""

# ╔═╡ c67f30f9-9101-462a-939e-fcde858f2bd3
md"""
## Problem Overview
"""

# ╔═╡ e23f5681-fb08-431a-b8ac-7d33686b7fad
md"""
Because the focus of this problem solving is the heat flow, introducing an additional (vector) variable of (negative) flux $\sigma = - \nabla u$ is advantageous.

Thus the partial differential equations we want to solve are:

$$\sigma + \nabla u=0\quad \text{in } \Omega,$$
$$\qquad \nabla\sigma =f\quad \text{in }\Omega.$$
$$\qquad\quad\ \  u=u_0\ \  \text{on } \Gamma_D,$$
$$\ \ \ \ -\sigma\cdot n=g\quad \text{on }\Gamma_N.$$


The same equations arise in connection with flow in porous media, where thery are referred to as Darcy flow.

After multiplying by test functions $\tau$ and $v$, integrating over the domain, and integrating term $\nabla\cdot\sigma v$ by parts, one obtains the following variational formulation: Find $\sigma\in\Sigma$ and $v\in V$ satisfying

$$\int_{\Omega}(\sigma\cdot\tau + \nabla u \cdot\tau)\text{ d}x = 0 \quad \forall \, \tau \in \Sigma$$

$$\int_{\Omega}\sigma\cdot \nabla v\text{ d}x = -\int_{\Omega}fv\text{ d}x - \int_{\Gamma_N}gv\text{ d}s\quad \forall \, v \in V$$

Compared to classical mixed formulation used in demo Mixed formulation for Poisson equation, the Dirichlet condition is here essential one and Neumann condition is natural.

To discretize the above formulation, two discrete function spaces $\Sigma_h \subset \Sigma$ and $V_h\subset V$ are needed to form a mixed function space $\Sigma_h \times V_h$. A stable choice of finite element spaces is to let $\Sigma_h$ be the discontinuous Raviart-Thomas elements of polynomial order $k$ and let $V_h$ be Lagrange elements of polynomial order $k+1$.
"""

# ╔═╡ ee7dd0e5-8faf-455d-abda-99bb5c4168bd
md"""
$((PlutoUI.LocalResource("2dmodel.png")))
"""

# ╔═╡ a2dda0d9-0e65-4292-a572-b17b710400ef
md"""
## Multi-field FE Spaces
"""

# ╔═╡ 34b37552-1fe9-4973-8a82-e4fb6ab6e869
md"""
Defining the constants
"""

# ╔═╡ 04230dac-dca0-41a8-8adf-e41e87e9fa7c
begin
	T_top = 11.0 # °C    -> Surface temperature
	𝜆 = 3.0 # W/mk  -> Thermal conductivity of the soil
	q = 0.08   # W/m^2 -> Heat flow (at the bottom)
	k = 1 #order
end

# ╔═╡ 8c1b37a6-9e64-4cff-916d-68cf5c48ad82
md"""
1. Lagrange element for temperature, which has a polynomial order of k+1 and dirichlet boundary condition on the top part of the boundary
"""

# ╔═╡ aa56c123-75b8-4ff8-86da-ef1414d4aac4
CGref = ReferenceFE(lagrangian, Float64, k+1)

# ╔═╡ d141053a-8c55-46c4-a5c6-816e0ce2279a
V = TestFESpace(mesh, CGref, conformity=:H1, dirichlet_tags="top")

# ╔═╡ b33f024c-4c47-4226-b5ed-0f1001594968
u_0(x) = T_top #Constant Temperature at the top (Diriclet Boundary)

# ╔═╡ ad904295-027a-4e5d-8db6-9f31987621c2
U = TrialFESpace(V,u_0)

# ╔═╡ 1b468ad2-18b0-45dd-9cfd-3f205471accd
md"""
2. Discontinuous Raviart-Thomas element for the heat flow with neumann boundary condition for the bottom part and a polynomial order of k
"""

# ╔═╡ bc0d9808-c2a4-4e20-a343-0e7cec922f1a
DRTref = ReferenceFE(raviart_thomas, Float64, k)

# ╔═╡ 9c66ab8e-6147-460e-8067-bab9cf08e979
τ = TestFESpace(mesh, DRTref, conformity=:Hdiv)

# ╔═╡ cffc9157-d432-4844-aaba-32037a23b3c8
σ = TrialFESpace(τ)

# ╔═╡ 44755203-0bc9-4386-ac87-8a18a51b1a06
md"""
3. Combining the single-field Test spaces into `VV` and Trial spaces into `UU` to get the mixed function spaces
"""

# ╔═╡ 8a72fa7f-72ea-426d-b9eb-96dd2adf35b6
VV = MultiFieldFESpace([τ,V])

# ╔═╡ 97db3417-29cb-471a-8f9c-47a7af33d157
UU = MultiFieldFESpace([σ,U])

# ╔═╡ bc0bfdd7-248e-4418-b248-9129e85100eb
md"""
## Numerical Integrations

For the computation of the integral on the domain $\Omega$ we define the integration mesh with `Triangulation`.
"""

# ╔═╡ 04440c4f-8c36-476e-b129-2d65e0fc3e4f
Ω = Triangulation(mesh)

# ╔═╡ c76c6d42-61da-4cd1-9e88-f5780ad5f7b6
dΩ = Measure(Ω,2)

# ╔═╡ ce77dda2-d1ef-45b8-9951-cb3e1e4a7fae
neumann_tags = ["bottom"]

# ╔═╡ f20b72ee-a187-4d31-8fc2-e4604ea50e1d
Γ = BoundaryTriangulation(mesh,tags=neumann_tags)

# ╔═╡ ccd2fa33-45e3-4802-b8d2-ff19bed4b381
dΓ = Measure(Γ,2)

# ╔═╡ 27b9583f-c315-4930-8c5e-9cfbf7de6865
md"""
## Defining the variational problem

To define the weak form, we need to define several functions.
"""

# ╔═╡ 6f102c98-dbca-4c9c-9093-9a190cb3142d
md"""
followed by the bi-linear form $a((σ,u),(τ,v)) = $$\int_{\Omega}(\sigma\cdot\tau + \nabla u \cdot\tau)\text{ d}\Omega$
"""

# ╔═╡ 1de8daa4-00f8-4364-b678-020fc311486f
a((σ,u),(τ,v)) = ∫( (σ⋅τ) + (𝜆*∇(u)⋅τ) )dΩ

# ╔═╡ f065b940-1f22-4ef2-ba97-e2414b40d481
md"""
and the linear form $L((τ,v)) = -\int_{\Omega}fv\text{ d}\Omega - \int_{\Gamma_N}gv\text{ d}\Gamma$.
"""

# ╔═╡ a60db35d-f8ac-4dc8-af93-eae05587e2e6
f(x) = 0.0

# ╔═╡ fbc127ad-23e3-4935-af53-9a6936e6423a
g(x) = q #constant heat flow at the bottom

# ╔═╡ b540fdb0-032f-4422-ba7d-e416c8a4a1f9
l((τ,v)) = -1*(∫(v*f)dΩ) + -1*(∫( v*q )dΓ)

# ╔═╡ 43a933b8-a613-422a-9c8e-2465029994e6
md"""
## Building the FE Problem

For the `AffineFEOperator` function in this case, we pass the multi-field test and trial spaces
"""

# ╔═╡ b195fcc6-35ba-4ec5-bca0-65113946bf93
op = AffineFEOperator(a, l, UU, VV)

# ╔═╡ 98016615-9fd2-446c-bf08-0d70dfe6fb33
md"""
## Solving the FE Problem
"""

# ╔═╡ 96e0610d-a17b-4970-8c9d-9c68c9b066d2
md"""
The solution is multi-field as we have a multi-field function space. The solution can be unpacked directly from the `solve` function. `σh` will capture heat flow and `uh` corresponds to temperature.
"""

# ╔═╡ eec8b332-e087-4654-8c2d-5d7be7d511df
σh,uh = solve(op)

# ╔═╡ b1cb51d3-9415-441e-836c-da3cba7066ae
md"""
**Unfortunately**, within the course of this project, I haven't been able to solve the problem, where more than one simple boundary conditions and/or function spaces are involved.

The generated error `"LinearAlgebra.SingularException(0)"` possibly indicates that:
- my definition of the multi-field system is still incorrect that causes lack of a unique solution
- the tolerance at the Dirichlet boundary to counteract numerical minimum deviations that was done in the python version is still necessary
"""

# ╔═╡ ee28167e-de89-467c-b40c-cef5dafb75a4
md"""
## Save solution as a vtk file for visualization

In case the problem were solved, the following cell will *react* and execute the generation of the vtk file automatically.
"""

# ╔═╡ 5ea7a834-5263-4d52-b5a7-97ad3522bfbc
writevtk(Ω,"$(Filename)_solution",cellfields=["heat flow"=>σh,"temperature"=>uh])

# ╔═╡ 4b702fc3-ead7-40dc-b543-9140b00c4064
md"""
# Next Step

There is a tutorial from Gridap on Darcy Equation that served as a good base for this notebook, but further research and experiment has to be conducted to implement the example in the mentioned tutorial for our use case in Geothermal Energy System. The potential of Gridap as an open-source tool for modelling using Finite Element Methods is still to be discovvered through a thorough research .
"""

# ╔═╡ 677fcfea-0eb3-429a-8558-54d39244afe8
md"""
# Function Definitions
"""

# ╔═╡ 1f942382-208b-4d24-89a9-10effad9ccbe
function write_geo(filename::String,density::Integer,x_end,y_end,x_0=0,y_top=0)	
	
	# Writing the content of the .geo file
	geo_text = """
	// Coordinates for rectangle corners
	Point(1) = {$x_0, -$y_end, 0, 1.0};
	Point(2) = {$x_0, $y_top, 0, 1.0};
	Point(3) = {$x_end, $y_top, 0, 1.0};
	Point(4) = {$x_end, -$y_end, 0, 1.0};
	
	// Lines connecting two points given on the right side
	Line(1) = {1, 2};
	Line(2) = {2, 3};
	Line(3) = {3, 4};
	Line(4) = {4, 1};

	// Connecting lines, creating a full rectangle
	Curve Loop(1) = {2, 3, 4, 1};

	// Defining the plane, later labeled as domain
	Plane Surface(1) = {1};
	Physical Surface("domain", 5) = {1};

	// Defining the boundaries and labeling it
	Physical Curve("top", 6) = {2};
	Physical Curve("sidewalls", 7) = {1, 3};
	Physical Curve("bottom", 8) = {4};
	
	// Preparing the geometry for meshing
	/// Defining the density of the mesh
	//// --> Dividing each side into n points, not how many rectangles!!!
	Transfinite Curve {1, 2, 3, 4} = $(density) Using Progression 1;
	Transfinite Surface {1};
	"""

	# Generating the .geo file
	cleanname = (split.(filename,["."]))[1][1]*".geo" #generally will accept all filename
	open(cleanname,"w") do file
		write(file, geo_text)
	end
end

# ╔═╡ c733e3e7-711d-46b3-9341-ced00e355324
begin
	x_min = 0     # m
	x_max = 7000  # m
	y_min = 0     # m
	y_max = 3600  # m

	n = 10

	write_geo(Filename,n,x_max,y_max)
end

# ╔═╡ cfef7f88-f6c8-48ca-8aa1-50d90f3646e8
function write_msh(filenamegeo::String,dim=2::Integer)
	filename = (split.(filenamegeo,["."]))[1][1] #removing any extension if given
	
	gmsh.initialize()              
	gmsh.model.add("Intro_Model")	# Defining a new model
	gmsh.open(filename*".geo")	# Reading the .geo file
	gmsh.model.mesh.generate(dim)		# Geenerating Mesh in 2D
	gmsh.write(filename*".msh")	# Writing the .msh file
	gmsh.finalize()
end

# ╔═╡ 3cb73eb9-67f5-42e8-b1ec-82666f57af52
write_msh(Filename)

# ╔═╡ f3a90024-9ae6-48a4-9f8e-41fbfecc35c4
md"""
# References

[1] https://gridap.github.io/Tutorials/stable/pages/t007_darcy/#Tutorial-7:-Darcy-equation-(with-RT)-1
[2] https://github.com/JuliaFEM/Gmsh.jl\
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Gmsh = "705231aa-382f-11e9-3f0c-b7cb4346fdeb"
Gridap = "56d4f2e9-7ea1-5844-9cf6-b9c51ca7ce8e"
GridapGmsh = "3025c34a-b394-11e9-2a55-3fee550c04c8"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Gmsh = "~0.3.1"
Gridap = "~0.18.2"
GridapGmsh = "~0.7.1"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "39d8c1a3f787a3ca0f817b23fdeedc792c781a1f"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "5c9b74c973181571deb6442d41e5c902e6b9f38e"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.12.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "8556500c18fcad8b4c44058e23fbc4a36143f6be"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.10.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
deps = ["Pkg"]
git-tree-sha1 = "daaeb6f7f77b88c072a83a2451801818acb5c63b"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "2.1.0"

[[deps.BSON]]
git-tree-sha1 = "4c3e506685c527ac6a54ccc0c8c76fd6f91b42fb"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "9a9610fbe5779636f75229e423e367124034af41"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.16.43"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "e24a6f390e5563583bb4315c73035b5b3f3e7ab4"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.4.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "260fd2400ed2dab602a7c15cf10c1933c59930a2"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.5"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

    [deps.Distances.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FLTK_jll]]
deps = ["Artifacts", "Fontconfig_jll", "FreeType2_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXft_jll", "Xorg_libXinerama_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "72a4842f93e734f378cf381dae2ca4542f019d23"
uuid = "4fce6fc7-ba6a-5f4c-898f-77e99806d6f8"
version = "1.3.8+0"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

    [deps.FillArrays.weakdeps]
    PDMats = "90014a1f-27ba-587c-ab20-58faa44d9150"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "2de436b72c3422940cbe1367611d137008af7ec3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.23.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg"]
git-tree-sha1 = "65af046f4221e27fb79b28b6ca89dd1d12bc5ec7"
uuid = "bd17208b-e95e-5925-bf81-e2f59b3e5c61"
version = "9.0.1+0"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.2.1+6"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Gmsh]]
deps = ["gmsh_jll"]
git-tree-sha1 = "6d815101e62722f4e323514c9fc704007d4da2e3"
uuid = "705231aa-382f-11e9-3f0c-b7cb4346fdeb"
version = "0.3.1"

[[deps.Gridap]]
deps = ["AbstractTrees", "BSON", "BlockArrays", "Combinatorics", "DataStructures", "DocStringExtensions", "FastGaussQuadrature", "FileIO", "FillArrays", "ForwardDiff", "JLD2", "JSON", "LineSearches", "LinearAlgebra", "NLsolve", "NearestNeighbors", "PolynomialBases", "QuadGK", "Random", "SparseArrays", "SparseMatricesCSR", "StaticArrays", "Statistics", "Test", "WriteVTK"]
git-tree-sha1 = "4918159105057659c25508ea2729fa7b5e4990ae"
uuid = "56d4f2e9-7ea1-5844-9cf6-b9c51ca7ce8e"
version = "0.18.2"

[[deps.GridapDistributed]]
deps = ["BlockArrays", "FillArrays", "Gridap", "LinearAlgebra", "MPI", "PartitionedArrays", "SparseArrays", "SparseMatricesCSR", "WriteVTK"]
git-tree-sha1 = "f5cb74ac35d39fedfb454d4614e25bc4d8d31569"
uuid = "f9701e48-63b3-45aa-9a63-9bc6c271f355"
version = "0.4.1"

[[deps.GridapGmsh]]
deps = ["Gridap", "GridapDistributed", "Libdl", "Metis", "PartitionedArrays", "gmsh_jll"]
git-tree-sha1 = "d57e69bba40c1e77bcf3a781aadbd846199ea251"
uuid = "3025c34a-b394-11e9-2a55-3fee550c04c8"
version = "0.7.1"

[[deps.HDF5_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "LibCURL_jll", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "OpenSSL_jll", "TOML", "Zlib_jll", "libaec_jll"]
git-tree-sha1 = "82a471768b513dc39e471540fdadc84ff80ff997"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.14.3+3"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1d334207121865ac8c1c97eb7f42d0339e4635bf"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.11.0+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be50fe8df3acbffa0274a744f1a99d29c45a57f4"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "59545b0a2b27208b0650df0a46b8e3019f85055b"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.4"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "PrecompileTools", "Reexport", "Requires", "TranscodingStreams", "UUIDs", "Unicode"]
git-tree-sha1 = "84642bc18a79d715b39d3724b03cbdd2e7d48c62"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.49"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearElasticity_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "71e8ee0f9fe0e86a8f8c7f28361e5118eab2f93f"
uuid = "18c40d15-f7cd-5a6d-bc92-87468d86c5db"
version = "5.0.0+0"

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

[[deps.METIS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1fd0a97409e418b78c53fac671cf4622efdf0f21"
uuid = "d00139f3-1899-568f-a2f0-47f597d42d70"
version = "5.1.2+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "80b2833b56d466b3858d565adcd16a4a05f2089b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.1.0+0"

[[deps.MMG_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "LinearElasticity_jll", "Pkg", "SCOTCH_jll"]
git-tree-sha1 = "70a59df96945782bb0d43b56d0fbfdf1ce2e4729"
uuid = "86086c02-e288-5929-a127-40944b0018b7"
version = "5.6.0+0"

[[deps.MPI]]
deps = ["Distributed", "DocStringExtensions", "Libdl", "MPICH_jll", "MPIPreferences", "MPItrampoline_jll", "MicrosoftMPI_jll", "OpenMPI_jll", "PkgVersion", "PrecompileTools", "Requires", "Serialization", "Sockets"]
git-tree-sha1 = "14cef41baf5b675b192b02a22c710f725ab333a7"
uuid = "da04e1cc-30fd-572f-bb4f-1f8673147195"
version = "0.20.20"

    [deps.MPI.extensions]
    AMDGPUExt = "AMDGPU"
    CUDAExt = "CUDA"

    [deps.MPI.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"

[[deps.MPICH_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Hwloc_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "4099bb6809ac109bfc17d521dad33763bcf026b7"
uuid = "7cb0a576-ebde-5e09-9194-50597f1243b4"
version = "4.2.1+1"

[[deps.MPIPreferences]]
deps = ["Libdl", "Preferences"]
git-tree-sha1 = "c105fe467859e7f6e9a852cb15cb4301126fac07"
uuid = "3da0fdf6-3ccc-4f1b-acd9-58baa6c99267"
version = "0.1.11"

[[deps.MPItrampoline_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "8c35d5420193841b2f367e658540e8d9e0601ed0"
uuid = "f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"
version = "5.4.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Metis]]
deps = ["CEnum", "LinearAlgebra", "METIS_jll", "SparseArrays"]
git-tree-sha1 = "5582d3b0d794280c9b818ba56ce2b35b108aca41"
uuid = "2679e427-3c69-5b7f-982b-ece356f1e94b"
version = "1.4.1"

    [deps.Metis.extensions]
    MetisGraphs = "Graphs"
    MetisLightGraphs = "LightGraphs"
    MetisSimpleWeightedGraphs = ["SimpleWeightedGraphs", "Graphs"]

    [deps.Metis.weakdeps]
    Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
    LightGraphs = "093fc24a-ae57-5d10-9952-331d41423f4d"
    SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"

[[deps.MicrosoftMPI_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f12a29c4400ba812841c6ace3f4efbb6dbb3ba01"
uuid = "9237b28f-5490-5468-be7b-bb81f5f5e6cf"
version = "10.1.4+2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "91a67b4d73842da90b526011fa85c5c4c9343fe0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.18"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OCCT_jll]]
deps = ["Artifacts", "FreeType2_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXft_jll", "Xorg_libXinerama_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "acc8099ae8ed10226dc8424fb256ec9fe367a1f0"
uuid = "baad4e97-8daa-5946-aac2-2edac59d34e1"
version = "7.6.2+2"

[[deps.OffsetArrays]]
git-tree-sha1 = "e64b4f5ea6b7389f6f046d13d4896a8f9c1ba71e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenMPI_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "MPIPreferences", "TOML"]
git-tree-sha1 = "e25c1778a98e34219a00455d6e4384e017ea9762"
uuid = "fe0851c0-eecd-5654-98d4-656369965a5c"
version = "4.1.6+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.PartitionedArrays]]
deps = ["CircularArrays", "Distances", "FillArrays", "IterativeSolvers", "LinearAlgebra", "MPI", "Printf", "Random", "SparseArrays", "SparseMatricesCSR"]
git-tree-sha1 = "149d2287770c6a533507d74beaa73d76c0727922"
uuid = "5a9dfac6-5c52-46f7-8278-5e2210713be9"
version = "0.3.4"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PolynomialBases]]
deps = ["ArgCheck", "AutoHashEquals", "FFTW", "FastGaussQuadrature", "LinearAlgebra", "Requires", "SimpleUnPack", "SpecialFunctions"]
git-tree-sha1 = "aa1877430a7e8b0c7a35ea095c415d462af0870f"
uuid = "c74db56a-226d-5e98-8bb0-a6049094aeea"
version = "0.4.21"

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

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SCOTCH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7110b749766853054ce8a2afaa73325d72d32129"
uuid = "a8d0f55d-b80e-548d-aff6-1a04c175f0f9"
version = "6.1.3+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseMatricesCSR]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38677ca58e80b5cad2382e5a1848f93b054ad28d"
uuid = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"
version = "0.6.7"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "20833c5b7f7edf0e5026f23db7f268e4f23ec577"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.6"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "d73336d81cafdc277ff45558bb7eaa2b04a8e472"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.10"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VTKBase]]
git-tree-sha1 = "c2d0db3ef09f1942d08ea455a9e252594be5f3b6"
uuid = "4004b06d-e244-455f-a6ce-a5f9919cc534"
version = "1.0.1"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams", "VTKBase"]
git-tree-sha1 = "48b9e8e9c83865e99e57f027d4edfa94e0acddae"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.19.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXft_jll]]
deps = ["Fontconfig_jll", "Libdl", "Pkg", "Xorg_libXrender_jll"]
git-tree-sha1 = "754b542cdc1057e0a2f1888ec5414ee17a4ca2a1"
uuid = "2c808117-e144-5220-80d1-69d4eaa9352c"
version = "2.3.3+1"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.gmsh_jll]]
deps = ["Artifacts", "Cairo_jll", "CompilerSupportLibraries_jll", "FLTK_jll", "FreeType2_jll", "GLU_jll", "GMP_jll", "HDF5_jll", "JLLWrappers", "JpegTurbo_jll", "LLVMOpenMP_jll", "Libdl", "Libglvnd_jll", "METIS_jll", "MMG_jll", "OCCT_jll", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "Xorg_libXft_jll", "Xorg_libXinerama_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bdc2fa0a123008ad941cabb0ad88c571e696af2e"
uuid = "630162c2-fc9b-58b3-9910-8442a8a132e6"
version = "4.13.0+1"

[[deps.libaec_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "46bf7be2917b59b761247be3f317ddf75e50e997"
uuid = "477f73a3-ac25-53e9-8cc3-50b2fa2566f0"
version = "1.1.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─70cfc41c-9845-4d83-96d2-dd31b96e8983
# ╠═d73852a0-3b80-11ef-2ed1-f108b2d1ceeb
# ╟─d8be6042-cb6c-473a-8193-0845c334601e
# ╟─2e1d62ba-0de3-46b3-a680-b15f4c27d9ca
# ╟─7f996195-6a95-4556-a9f7-b573a97964b7
# ╟─da73ce21-1023-4309-becc-b5c58cbfdff3
# ╟─432fabbb-9d9f-46f5-a497-693f247af2eb
# ╟─e9caa95e-ec8c-411f-82db-1d6ab163a690
# ╠═4b7ee5b5-2edc-414b-b695-5b2c009565b4
# ╟─46ecf0b8-06bf-47ce-b661-74a6a7a6b12a
# ╟─c1469e5f-385b-470a-8898-8f4a8594442b
# ╠═c733e3e7-711d-46b3-9341-ced00e355324
# ╟─72da26fd-35f3-4957-a43a-a0139bf2616b
# ╠═3cb73eb9-67f5-42e8-b1ec-82666f57af52
# ╟─254b17e4-dc7d-4bf6-b71e-b7fad3f5990d
# ╠═30733946-1cee-4ff0-911b-6d2eb943a134
# ╠═7cfea374-92e4-4c43-bfbc-6ad9d1a8b66c
# ╟─994c101c-98ec-478a-9b72-7389483d142f
# ╟─c67f30f9-9101-462a-939e-fcde858f2bd3
# ╟─e23f5681-fb08-431a-b8ac-7d33686b7fad
# ╟─ee7dd0e5-8faf-455d-abda-99bb5c4168bd
# ╟─a2dda0d9-0e65-4292-a572-b17b710400ef
# ╟─34b37552-1fe9-4973-8a82-e4fb6ab6e869
# ╠═04230dac-dca0-41a8-8adf-e41e87e9fa7c
# ╟─8c1b37a6-9e64-4cff-916d-68cf5c48ad82
# ╠═aa56c123-75b8-4ff8-86da-ef1414d4aac4
# ╠═d141053a-8c55-46c4-a5c6-816e0ce2279a
# ╠═b33f024c-4c47-4226-b5ed-0f1001594968
# ╠═ad904295-027a-4e5d-8db6-9f31987621c2
# ╟─1b468ad2-18b0-45dd-9cfd-3f205471accd
# ╠═bc0d9808-c2a4-4e20-a343-0e7cec922f1a
# ╠═9c66ab8e-6147-460e-8067-bab9cf08e979
# ╠═cffc9157-d432-4844-aaba-32037a23b3c8
# ╟─44755203-0bc9-4386-ac87-8a18a51b1a06
# ╠═8a72fa7f-72ea-426d-b9eb-96dd2adf35b6
# ╠═97db3417-29cb-471a-8f9c-47a7af33d157
# ╟─bc0bfdd7-248e-4418-b248-9129e85100eb
# ╠═04440c4f-8c36-476e-b129-2d65e0fc3e4f
# ╠═c76c6d42-61da-4cd1-9e88-f5780ad5f7b6
# ╠═ce77dda2-d1ef-45b8-9951-cb3e1e4a7fae
# ╠═f20b72ee-a187-4d31-8fc2-e4604ea50e1d
# ╠═ccd2fa33-45e3-4802-b8d2-ff19bed4b381
# ╟─27b9583f-c315-4930-8c5e-9cfbf7de6865
# ╟─6f102c98-dbca-4c9c-9093-9a190cb3142d
# ╠═1de8daa4-00f8-4364-b678-020fc311486f
# ╟─f065b940-1f22-4ef2-ba97-e2414b40d481
# ╠═a60db35d-f8ac-4dc8-af93-eae05587e2e6
# ╠═fbc127ad-23e3-4935-af53-9a6936e6423a
# ╠═b540fdb0-032f-4422-ba7d-e416c8a4a1f9
# ╟─43a933b8-a613-422a-9c8e-2465029994e6
# ╠═b195fcc6-35ba-4ec5-bca0-65113946bf93
# ╟─98016615-9fd2-446c-bf08-0d70dfe6fb33
# ╟─96e0610d-a17b-4970-8c9d-9c68c9b066d2
# ╠═eec8b332-e087-4654-8c2d-5d7be7d511df
# ╟─b1cb51d3-9415-441e-836c-da3cba7066ae
# ╟─ee28167e-de89-467c-b40c-cef5dafb75a4
# ╠═5ea7a834-5263-4d52-b5a7-97ad3522bfbc
# ╟─4b702fc3-ead7-40dc-b543-9140b00c4064
# ╟─677fcfea-0eb3-429a-8558-54d39244afe8
# ╟─1f942382-208b-4d24-89a9-10effad9ccbe
# ╟─cfef7f88-f6c8-48ca-8aa1-50d90f3646e8
# ╟─f3a90024-9ae6-48a4-9f8e-41fbfecc35c4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
