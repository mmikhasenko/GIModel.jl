# Strong decays beyond the paper

Godfrey and Isgur computed strong decays ``A \to B\,\pi`` with the
elementary-emission operator of Eq. (19). For their Table V they used one
simplification: every meson got the same Gaussian wavefunction, with
``\beta = 0.40`` GeV. GIModel can instead evaluate the same operator between
the solved wavefunctions, each with its own size, radial nodes and mixing.
This tutorial builds the states, fixes the two couplings, and predicts decays
that were not used in the fit.

```@example strong
using GIModel
using GIModel.QuarkModelTransitions
params, mq = load_parameters_and_quark_masses(default_parameters_path())
solver = OscillatorSolver()
nothing # hide
```

## The operator

```math
H_\text{int} \propto g\,\boldsymbol\sigma\cdot\boldsymbol q
                 + h\,\boldsymbol\sigma\cdot\boldsymbol p'
```

A quark or antiquark in the parent emits a pseudoscalar with momentum
``\boldsymbol q``; ``\boldsymbol p'`` is the momentum of the emitting quark.
The direct coupling ``g`` and the recoil coupling ``h`` (both in
``\mathrm{GeV}^{-1}``) are phenomenological. The amplitude is linear in them,
so two numbers fitted to data fix every decay.

## Building the states

Solve the flavor channels we need. Each `BasisState` is a level to keep:

```@example strong
S1, S0 = BasisState(1, "S", 3, 1), BasisState(1, "S", 1, 0)
P1_triplet, P1_singlet = BasisState(1, "P", 3, 1), BasisState(1, "P", 1, 1)
solve(f1, f2, levels) = compute_spectrum(params, Meson(mq, f1, f2); levels, solver)
ud = solve(:u, :d, [S1, S0, P1_triplet, P1_singlet])
uu, dd = solve(:u, :u, [S1, S0]), solve(:d, :d, [S1, S0])
us, ds = solve(:u, :s, [S1, S0]), solve(:d, :s, [S0])
nothing # hide
```

Charged states come straight from their spectra:

```@example strong
rho, pi_plus, a1, b1 = (physical_state(ud, b) for b in (S1, S0, P1_triplet, P1_singlet))
K_star, K_plus, K_zero = physical_state(us, S1), physical_state(us, S0), physical_state(ds, S0)
(ρ = rho.mass_GeV, a₁ = a1.mass_GeV, b₁ = b1.mass_GeV, K⁺ = K_plus.mass_GeV)
```

Neutral isospin states are superpositions of ``u\bar u`` and ``d\bar d``. A
[`PhysicalState`](@ref) can be built from any signed components, so a small
helper combines solved states:

```@example strong
function combine(label, parts)
    components = [(basis = c.basis, coefficient = w * c.coefficient, wave = c.wave)
                  for (state, w) in parts for c in state.components]
    mass = sum(w^2 * state.mass_GeV for (state, w) in parts)
    PhysicalState(label, mass, components; provenance = (source = :isospin_combination,))
end
pi_zero = combine("π⁰", [(physical_state(uu, S0), 1 / √2), (physical_state(dd, S0), -1 / √2)])
omega = combine("ω", [(physical_state(uu, S1), 1 / √2), (physical_state(dd, S1), 1 / √2)])
pi_zero
```

The flavor components must be explicit (`:u`, `:d`): the averaged `:q` cannot
say whether a state is ``u\bar u + d\bar d`` or ``u\bar u - d\bar d``, so the
operator rejects it.

## Channels

A [`TwoMesonChannel`](@ref) is ordered as `(surviving, emitted)`, and the
emitted meson must be a pseudoscalar:

```@example strong
channels = [
    "ρ⁺ → π⁺ π⁰" => (TwoMesonChannel(pi_plus, pi_zero), rho),
    "K*⁺ → K⁰ π⁺" => (TwoMesonChannel(K_zero, pi_plus), K_star),
    "K*⁺ → K⁺ π⁰" => (TwoMesonChannel(K_plus, pi_zero), K_star),
    "a₁⁺ → ρ⁺ π⁰" => (TwoMesonChannel(rho, pi_zero), a1),
    "b₁⁺ → ω π⁺" => (TwoMesonChannel(omega, pi_plus), b1),
]
amplitude(name, g, h) = matrix_element(first(Dict(channels)[name]),
    PseudoscalarEmission(g, h, mq), last(Dict(channels)[name]))
nothing # hide
```

For two pseudoscalars in the final state, as in ``\rho \to \pi\pi``, either one
could be called the emitted meson. The operator evaluates the assignment you
give and never adds the other one silently.

## Fixing the couplings

**The ratio ``h/g``.** In the paper's single-Gaussian limit the Table IV
amplitudes reduce to two combinations, ``A = (g + h/4)\beta`` and
``S_0 = 3h\beta``. The paper's fit gives ``S_0/A = 3.291/1.644``, which fixes
the ratio:

```@example strong
r = 3.291 / 1.644
h_over_g = r / (3 - r / 4)
```

**The scale.** Fix ``g`` so that ``\Gamma(\rho \to \pi\pi) = 149.1`` MeV.
Since the width is quadratic in the couplings, one evaluation at ``g = 1``
is enough:

```@example strong
unit_width = decay_width(amplitude("ρ⁺ → π⁺ π⁰", 1.0, h_over_g))
g = sqrt(149.1 / unit_width)
h = h_over_g * g
(g = g, h = h)
```

## Predictions

```@example strong
for (name, _) in channels
    a = amplitude(name, g, h)
    waves = partial_waves(a)
    line = rpad(name, 14) * "q = $(round(Int, 1000a.momentum_GeV)) MeV   Γ = $(round(decay_width(a); digits = 1)) MeV"
    if length(waves) == 2
        line *= "   D/S = $(round(real(a[waves[2]] / a[waves[1]]); digits = 3))"
    end
    println(line)
end
```

Compare with measured values:

| decay | model | measured |
|---|---|---|
| ``\Gamma(K^{*+}\to K\pi)`` | sum of the two lines above | ``50.8 \pm 0.9`` MeV |
| ``K^0\pi^+ : K^+\pi^0`` | ``2 : 1`` | ``2 : 1`` (isospin) |
| ``b_1 \to \omega\pi``, D/S | see above | ``0.277 \pm 0.027`` |
| ``a_1 \to \rho\pi``, D/S | see above | ``-0.062 \pm 0.020`` |
| ``\Gamma(b_1)`` | see above | ``142 \pm 9`` MeV (``\omega\pi`` dominant) |

The D/S ratios test the internal structure of the P-wave mesons and do not
depend on the overall scale ``g``. The ``b_1`` ratio is reproduced well, and
the ``a_1`` ratio has the measured sign. The ``K^*`` width is overestimated:
the ratio ``h/g`` came from the paper's Gaussian-wavefunction fit, and the
decay momentum uses the model masses, which are a few MeV off.

## Kinematics use model masses

Momenta and widths above use the masses of the solved states: the model pion
is 149 MeV rather than 138 MeV, and the model ``K^*`` is 902 MeV rather than
892 MeV. [`mass_correction_factor`](@ref) shows how an amplitude would change at
the measured momentum, keeping the wavefunctions fixed:

```@example strong
a = amplitude("K*⁺ → K⁰ π⁺", g, h)
wave = only(partial_waves(a))
factor = mass_correction_factor(first(Dict(channels)["K*⁺ → K⁰ π⁺"]),
    PseudoscalarEmission(g, h, mq), K_star; target_momentum = 0.291, partial_wave = wave)
```

The factor multiplies that partial-wave amplitude; the phase-space factor
``q`` in the width changes separately.

## How the couplings enter

Because the amplitude is linear in ``g`` and ``h``, the D/S ratio of an axial
meson depends only on ``h/g``. Scanning it shows how sensitive each ratio is:

```@example strong
using CairoMakie
ratios = range(0.0, 1.6; length = 41)
d_over_s(name, x) = (a = amplitude(name, 1.0, x); w = partial_waves(a); real(a[w[2]] / a[w[1]]))
fig = Figure(size = (500, 320))
ax = Axis(fig[1, 1]; xlabel = "h / g", ylabel = "D / S")
lines!(ax, ratios, [d_over_s("b₁⁺ → ω π⁺", x) for x in ratios]; label = "b₁ → ωπ")
lines!(ax, ratios, [d_over_s("a₁⁺ → ρ⁺ π⁰", x) for x in ratios]; label = "a₁ → ρπ")
vlines!(ax, [h_over_g]; color = :gray, linestyle = :dash)
axislegend(ax; position = :lt)
fig
```

The dashed line is the ratio from the paper's fit. A global fit of ``g`` and
``h`` to several widths and ratios would be the next step for a quantitative
study.

## What this operator is, and is not

The emitted pion is treated as an elementary field: its mass and flavor
enter, but its wavefunction is not integrated. This is the Godfrey–Isgur
emission model, not a quark-pair-creation (``{}^3P_0``) model. It cannot
describe decays that need a new quark pair, such as ``a_1 \to K^*\bar K``
through ``s\bar s`` creation, or decays into two non-pseudoscalar mesons.
See [Scope and limitations](@ref).
