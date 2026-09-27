### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ df79ef8e-ba55-11f1-2806-4f8c87fd148f
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(@__DIR__)
	# 
	using GIModel.QuarkModelTransitions
	using GIModel
end

# ╔═╡ f1f369b1-954a-4ba9-91f9-87d0032d2f80
params, mq = load_parameters_and_quark_masses(default_parameters_path())

# ╔═╡ dc5a703d-ce66-49c4-b9e5-01827c449df2
spec = compute_spectrum(params, Meson(mq, :c, :c); levels = spectrum_levels(1))

# ╔═╡ e3431728-1c51-481f-b31f-f961540563c3
md"""
## Transitions
"""

# ╔═╡ 017f9de4-931d-49ec-b9e9-df2a3d9321dc
QuarkModelTransitions

# ╔═╡ Cell order:
# ╠═df79ef8e-ba55-11f1-2806-4f8c87fd148f
# ╠═f1f369b1-954a-4ba9-91f9-87d0032d2f80
# ╠═dc5a703d-ce66-49c4-b9e5-01827c449df2
# ╟─e3431728-1c51-481f-b31f-f961540563c3
# ╠═017f9de4-931d-49ec-b9e9-df2a3d9321dc
