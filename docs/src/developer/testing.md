# Testing and CI

## Test suites

```bash
julia --project=. -e 'using Pkg; Pkg.test()'          # GIModel + QuarkModelTransitions
julia --project=GIPaper -e 'using Pkg; Pkg.test()'    # GIPaper
```

The default GIModel suite is the quick one that CI runs. It includes the
transition tests and [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality
checks. Tests check invariants and independent limits (closed-form integrals,
sum rules, agreement of the two solvers) rather than frozen numbers where
possible.

Solver-convergence sweeps take a few minutes and run only on request:

```bash
GI_HEAVY_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

## Full verification

```bash
bash scripts/verify_packages.sh    # both test suites
bash scripts/verify_project.sh     # heavy tests, GIPaper tests, and every comparison report
```

`verify_project.sh` regenerates all reports in `GIPaper/reports/`. Review the
diff of those files: a physics change shows up there first.

## Documentation

The documentation is built with
[Documenter.jl](https://documenter.juliadocs.org) and
[DocumenterVitepress.jl](https://luxdl.github.io/DocumenterVitepress.jl), from
the repository root:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build runs every example on every page, so a successful build means the
documented workflows work. Node.js is needed for the VitePress stage;
`DOCS_MD_ONLY=true` stops after the Markdown stage, and `DOCS_DRAFT=true`
skips running examples while you edit prose. After a full build, preview the
site with `npm run docs:dev` from `docs/`.

Pages under **Recorded reports** are copied from `GIPaper/reports/` at build
time; regenerate the reports first if the physics changed.

## Continuous integration

| workflow | runs | on |
|---|---|---|
| `Test.yml`, `TestOnPRs.yml` | both test suites | pushes to `main`, pull requests |
| `Documentation.yml` | docs build and deployment to GitHub Pages | pushes to `main`, tags, pull requests (preview) |
| `Lint.yml` | offline link check of Markdown files | pushes, pull requests |
| `CompatHelper.yml`, `TagBot.yml` | dependency bounds, release tags | scheduled |

Deployment uses the `DOCUMENTER_KEY` secret. Documentation for pull requests
is deployed as a preview under `previews/PR<number>`.
