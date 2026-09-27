"""
Vacuum final state for a reduced meson-current matrix element.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = LeptonicCurrent(:electromagnetic, masses)
@assert isfinite(matrix_element(Vacuum(), operator, initial).value)
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct Vacuum end
"""
Two real photons; the implemented result is GI's integrated width amplitude.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = TwoPhotonAnnihilation(masses, AnnihilationTerm((:c, :c), 4/9))
@assert decay_width(TwoPhotonChannel(), operator, final) > 0
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct TwoPhotonChannel end
"""
Two gluons in the published lowest-order integrated-rate prescription.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = GluonicAnnihilation(masses, 0.3) # Explicit fixed coupling.
@assert decay_width(TwoGluonChannel(), operator, final) > 0
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct TwoGluonChannel end
"""
Three gluons in the published lowest-order integrated-rate prescription.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = GluonicAnnihilation(masses, 0.3)
@assert decay_width(ThreeGluonChannel(), operator, initial) > 0
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct ThreeGluonChannel end

"""
Charged lepton plus massless neutrino; CKM is supplied explicitly.

## Example

```julia
using GIModel, QuarkModelTransitions
params, masses = load_parameters_and_quark_masses(default_parameters_path())
basis = BasisState(1, "S", 1, 0)
spectrum = compute_spectrum(params, Meson(masses, :u, :s); levels=[basis],
    solver=FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1))
kaon = physical_state(spectrum, basis)
operator = LeptonicCurrent(:P_P, masses, AnnihilationTerm((:u, :s), 2sqrt(3)))
channel = LeptonNeutrinoChannel(0.10566; ckm=0.225)
@assert decay_width(channel, operator, kaon) > 0
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct LeptonNeutrinoChannel
    lepton_mass_GeV::Float64
    ckm::Float64
    function LeptonNeutrinoChannel(m::Real; ckm::Real)
        isfinite(m) && m >= 0 && isfinite(ckm) && 0 <= ckm <= 1 ||
            throw(ArgumentError("invalid lepton mass or CKM magnitude"))
        new(m, ckm)
    end
end
"""
Lepton pair in the massless approximation of GI Eq. D8.

## Example

```julia
using GIModel, QuarkModelTransitions
import QuarkModelTransitions as QMT
params, masses = load_parameters_and_quark_masses(default_parameters_path())
levels = [BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)]
# A modest grid for this example; check convergence for quantitative widths.
solver = FiniteDifferenceSolver(ngrid=240, rmax=24.0, nlevels_per_channel=1)
spectrum = compute_spectrum(params, Meson(masses, :c, :c); levels, solver)
initial = physical_state(spectrum, levels[1])
final = physical_state(spectrum, levels[2])
operator = LeptonicCurrent(:electromagnetic, masses)
@assert decay_width(MasslessLeptonPair(), operator, initial) > 0
```

## Related

[`matrix_element`](@ref), [`decay_width`](@ref).
"""
struct MasslessLeptonPair end
