"""Published phenomenological defaults used by electromagnetic observables.

These are additional transition-model inputs, not spectrum parameters or
numerical tolerances. Individual functions accept keyword overrides.
"""
const ELECTROMAGNETIC_DEFAULTS = (
    magnetic_exponent = 0.7,
    electric_exponent = 0.5,
    charge_radius_exponent = 0.2,
    recoil_beta_GeV = 0.40,
)

"""Calibration defaults for the frozen Table IV/V model (amplitudes in MeV^1/2)."""
const STRONG_DECAY_DEFAULTS = (
    rho_amplitude = 12.4,
    B_amplitude = -11.0,
    beta_GeV = 0.40,
)
