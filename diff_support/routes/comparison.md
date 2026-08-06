# FD and HO route comparison

The two routes share physics parameters and the final symmetric spectral rule,
but not their numerical representation or preparation.

| Aspect | Finite difference | Harmonic oscillator |
|---|---|---|
| Representation | radial interior grid | finite analytic HO basis |
| Typical central matrix size | `ngrid`, default 450 | `nbasis`, default 24 |
| `p²` | FD tridiagonal from `h,r,L` | exact tridiagonal from `beta,L` |
| Position potential | diagonal samples on grid | Gauss--Laguerre/DVR matrix |
| Continuous physics inputs | `b,c,sigma0,s,m1,m2` | same |
| Extra discrete context | `ngrid,rmax` | `beta,nbasis,nq` |
| Adaptive/discrete choice | optional eigensolver only | quadrature order and beta-grid selection |
| Spatial mesh in central H | yes | no on active Appendix-A path |
| Constant work to prepare | grid, `p²` eigensystem | `p²` eigensystem, quadrature nodes/projector |
| Primary convergence check | grid refinement | basis/quadrature convergence and beta margin |
| Central spectral derivative | `v' H_dot v` | `v' H_dot v` on fixed-beta branch |
| Main nonsmoothness | eigenvalue/state crossings | crossings plus beta-branch switches |

## Shared API, separate dispatch

Production ownership is now `FiniteDifferenceSolver` versus
`OscillatorSolver`, with no basis marker in `GIParameters`. A future
user-facing name may be shared, but the prepared implementations should still
dispatch separately:

```julia
central_masses(theta, prepared::PreparedFDProblem)
central_masses(theta, prepared::PreparedHOBranch)
```

The HO orchestration layer additionally selects a branch and reports its margin.
There should be no boolean `basis=:fd/:ho` inside one large differentiated
function.

## Shared validation, route-specific tolerances

Both routes should use the same physics points and derivative directions where
possible, but their convergence tolerances cannot be assumed identical.

For each route, compare:

1. current primal implementation vs prepared implementation;
2. analytic checks such as `dE/dc=1`;
3. direct finite differences vs explicit spectral derivatives;
4. backend JVP/VJP vs the explicit spectral oracle;
5. coarse vs refined numerical representation.

Cross-route FD/HO agreement is a useful physics diagnostic, but it is not the
definition of derivative correctness. Each is an approximation with its own
convergence controls.
