using Aqua

# Package-quality checks: method ambiguities, type piracy, unbound type
# parameters, stale dependencies and compat bounds. The persistent-task check
# precompiles a wrapper package and is left to local release checks.
@testset "Aqua" begin
    Aqua.test_all(GIModel; persistent_tasks = false)
end
