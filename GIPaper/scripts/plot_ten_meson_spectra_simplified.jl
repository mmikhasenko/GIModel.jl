#!/usr/bin/env julia
# Low-lying, column-compressed companion to plot_ten_meson_spectra.jl.

"--simplified" in ARGS || push!(ARGS, "--simplified")
include("plot_ten_meson_spectra.jl")
