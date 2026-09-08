### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ ce7a5438-ab59-11f1-3133-6917dafe3591
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	# 
	using GIModel
end

# ╔═╡ 780e10b9-5d0d-48e7-b1ba-ce16ba5623bf
Meson(HeavyQuark{:up}(1.628, :c), HeavyQuark{:up}(1.628, :c))

# ╔═╡ f16c704a-08e3-4582-9aca-450972a648e7
params, mq = load_parameters_and_quark_masses("../data/parameters.provisional.toml")

# ╔═╡ 5f3ba5ff-c58c-4af0-a24b-12d23cd9240d
Meson(mq, :c, :b)

# ╔═╡ c68e5a8a-2618-4acc-87ee-4cd69fdbb7f8
params

# ╔═╡ 142cc7fc-4e71-490b-a888-a3ddca1810e5
GIParameters

# ╔═╡ fe159a8e-6d93-4e7f-9080-011532210ac5
mq

# ╔═╡ d2c39c69-1dd2-4317-b915-3fbb4a1db0b7
my_meson = Meson(mq, :c, :b)

# ╔═╡ 90a315cc-fe43-4612-a816-ce17503ebf66
spectrum_levels(1)

# ╔═╡ f5167fb7-d6e3-468b-bbcb-9c9fbb19c7c9
spec = compute_spectrum(params, my_meson; levels = spectrum_levels(2))

# ╔═╡ a9733ab5-f534-49f6-81e3-490fd72aa995
spec.states[1]

# ╔═╡ 547ee88d-3494-49bc-bdda-f201a257696a
spec.states[6].corrected

# ╔═╡ 7157d6be-e1bb-4cbb-96e2-38ea5a5956ea
BasisState

# ╔═╡ 17453553-bca3-4b2a-b3f5-92fded6e55a4
MixedState

# ╔═╡ Cell order:
# ╠═ce7a5438-ab59-11f1-3133-6917dafe3591
# ╠═5f3ba5ff-c58c-4af0-a24b-12d23cd9240d
# ╠═780e10b9-5d0d-48e7-b1ba-ce16ba5623bf
# ╠═f16c704a-08e3-4582-9aca-450972a648e7
# ╠═c68e5a8a-2618-4acc-87ee-4cd69fdbb7f8
# ╠═142cc7fc-4e71-490b-a888-a3ddca1810e5
# ╠═fe159a8e-6d93-4e7f-9080-011532210ac5
# ╠═d2c39c69-1dd2-4317-b915-3fbb4a1db0b7
# ╠═90a315cc-fe43-4612-a816-ce17503ebf66
# ╠═f5167fb7-d6e3-468b-bbcb-9c9fbb19c7c9
# ╠═a9733ab5-f534-49f6-81e3-490fd72aa995
# ╠═547ee88d-3494-49bc-bdda-f201a257696a
# ╠═7157d6be-e1bb-4cbb-96e2-38ea5a5956ea
# ╠═17453553-bca3-4b2a-b3f5-92fded6e55a4
