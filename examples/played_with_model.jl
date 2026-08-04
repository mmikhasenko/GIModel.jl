### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ce2ee43c-82a7-4438-a428-de2798083c90
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.offline()
	# GIPaper (comparison layer) dev-depends on GIModel (pure computation).
	Pkg.activate(joinpath(@__DIR__, "..", "GIPaper"))
	Pkg.instantiate()
	#
	using GIModel
	using GIPaper
	using TOML
end

# ╔═╡ 89136d22-4856-41e1-9c84-603f6a6a8193
md"""
## Parameters
"""

# ╔═╡ 69a25a08-e3dd-46e1-b9fa-3d7177d7bb66
params, all_masses = GIModel.load_parameters_and_quark_masses(GIPaper.model_parameters_path());

# ╔═╡ ef068e70-bf99-4b66-b058-4eb2d834d204
params

# ╔═╡ 73ff8dab-2349-4615-b485-d506059a0638
GIParameters

# ╔═╡ da043069-41bd-4a5c-a7ec-f49459255fa9
all_masses

# ╔═╡ 756a1ce4-3d68-4ae9-9e0e-a672b78245c1
md"""
### Bottomonium

The model no longer needs reference data to run: specify the meson by quark
flavors and the levels to compute, and `compute_spectrum` returns an organized
`Spectrum`. The comparison with the paper (GIPaper) is a separate step.
"""

# ╔═╡ c55cbbcf-13d0-426e-acb8-90f6a03cf762
label = "bottomonium"

# ╔═╡ e113ec80-2eba-4da1-b3be-0e74569449d7
reference_bbbar = GIPaper.load_reference_spectrum(GIPaper.reference_spectrum_path(label))

# ╔═╡ 3d960667-667e-4961-b5c3-18ecd7706b49
bbbar = Meson(all_masses, :b, :b)

# ╔═╡ 3df46f57-c312-479a-b6ce-7ebf1e078e4b
levels = spectrum_levels(3; L_labels = ("S", "P", "D"))

# ╔═╡ edcb9a64-f0c5-464f-9c0a-6cf853187f1e
central = central_spectrum(params, bbbar; levels)

# ╔═╡ 7b0ce4f4-5985-49bf-80c2-afce8ec54de8
spectrum = compute_spectrum(
	params,
	bbbar;
	levels = spectrum_levels(2; L_labels = ("S", "P", "D")),
	solver = RadialSolver(; kinetic = :relativistic),
)

# ╔═╡ 35378712-12cf-4322-9683-f58b5048416a
spectrum.states

# ╔═╡ 6ffbc3ab-3d46-4212-b943-19deb71fc71c
[(s.label, round(s.mass_GeV; digits = 4)) for s in spectrum.states]

# ╔═╡ fad8f37d-8938-47f9-a978-6b943a37b345
rows = compare_reference(
	params,
	all_masses,
	reference_bbbar;
	contact_hyperfine = true,
	use_fine_structure = params.fine_structure.enabled,
)

# ╔═╡ fddc66c5-3a9c-45da-b792-0d3fa7a84ff7
# GIPaper.write_residual_report("bottomonium_residuals.md", "Residuals: bottomonium", rows)

# ╔═╡ Cell order:
# ╠═ce2ee43c-82a7-4438-a428-de2798083c90
# ╟─89136d22-4856-41e1-9c84-603f6a6a8193
# ╠═69a25a08-e3dd-46e1-b9fa-3d7177d7bb66
# ╠═ef068e70-bf99-4b66-b058-4eb2d834d204
# ╠═73ff8dab-2349-4615-b485-d506059a0638
# ╠═da043069-41bd-4a5c-a7ec-f49459255fa9
# ╟─756a1ce4-3d68-4ae9-9e0e-a672b78245c1
# ╠═c55cbbcf-13d0-426e-acb8-90f6a03cf762
# ╠═e113ec80-2eba-4da1-b3be-0e74569449d7
# ╠═3d960667-667e-4961-b5c3-18ecd7706b49
# ╠═3df46f57-c312-479a-b6ce-7ebf1e078e4b
# ╠═edcb9a64-f0c5-464f-9c0a-6cf853187f1e
# ╠═7b0ce4f4-5985-49bf-80c2-afce8ec54de8
# ╠═35378712-12cf-4322-9683-f58b5048416a
# ╠═6ffbc3ab-3d46-4212-b943-19deb71fc71c
# ╠═fad8f37d-8938-47f9-a978-6b943a37b345
# ╠═fddc66c5-3a9c-45da-b792-0d3fa7a84ff7
