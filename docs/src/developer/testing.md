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

The site is built in two stages.

1. **Executed pages** (Getting started, most of the Manual, all Tutorials) are
   written as Quarto notebooks in `docs/quarto/*.qmd`. `julia docs/render.jl`
   runs them with [Quarto](https://quarto.org) in the environment
   `docs/quarto/Project.toml` and writes Markdown with outputs and figures to
   `docs/src/`. These generated files are committed. Rendering everything
   takes about 20–30 minutes; pass file names to render only some pages:

   ```bash
   julia docs/render.jl tutorials/charmonium.qmd
   ```

2. **Site assembly** runs no calculations:

   ```bash
   julia --project=docs -e 'using Pkg; Pkg.instantiate()'
   julia --project=docs docs/make.jl
   ```

   Documenter resolves cross-references and API docstrings, and
   DocumenterVitepress builds the static site with its bundled Node.js.
   `DOCS_MD_ONLY=true` stops after the Markdown stage. Preview a full build
   with `npm run docs:dev` from `docs/`.

Edit the `.qmd` file, never the generated `.md` next to it. Each generated page
records the hash of its source, and `make.jl` refuses to build when a source
changed without re-rendering. **Before a release, re-render all pages** so the
published outputs match the released code.

Pages under **Recorded reports** are copied from `GIPaper/reports/` at build
time; regenerate the reports first if the physics changed.

## Continuous integration

| workflow | runs | on |
|---|---|---|
| `Test.yml`, `TestOnPRs.yml` | both test suites | pushes to `main`, pull requests |
| `Documentation.yml` | site assembly from the committed pages and deployment to GitHub Pages (no calculations) | pushes to `main`, tags, pull requests (preview) |
| `Lint.yml` | offline link check of Markdown files | pushes, pull requests |
| `LearningSheets.yml` | builds the course PDFs and publishes them to the `learning-sheets` release | changes under `LearningTrack/` |
| `CompatHelper.yml`, `TagBot.yml` | dependency bounds, release tags | scheduled |

Deployment uses the `DOCUMENTER_KEY` secret. Documentation for pull requests
is deployed as a preview under `previews/PR<number>`.
