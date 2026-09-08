# Explicit diagnostic/plot adapter. Production code consumes `radial_wave`
# directly; tests that validate sampled arrays request that view deliberately.
function sampled_arrays(sol::ChannelRadialSolution; grid = nothing)
    sampled = if all(w -> w isa MeshWave, sol.waves)
        sol.waves
    else
        isnothing(grid) &&
            throw(ArgumentError("sampling an oscillator solution requires `grid`"))
        [sample_wave(w, grid) for w in sol.waves]
    end
    r = isnothing(grid) ? first(sampled).r : collect(Float64, grid)
    return sol.eigenvalues_GeV, hcat((w.u for w in sampled)...), r
end
