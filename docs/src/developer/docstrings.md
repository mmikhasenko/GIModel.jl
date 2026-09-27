# Writing docstrings

Docstrings are the first documentation most users read, in the REPL or in a
Pluto notebook. A good GIModel docstring answers three questions: **what is
this, how do I get it, and what can I do with it next?**

## Conventions

- Start with the signature, the purpose, and the units. Show real property
  access, such as `state.mass_GeV`.
- Keep example blocks to three or four lines, each showing one action. A
  longer list of independent inspection calls is fine.
- State prerequisites and link to the functions that produce them. For the
  bundled parameters use [`default_parameters_path`](@ref), so examples work
  from any directory.
- End with one `## Related` section: relevant links with a short description
  each. Link with ``[`name`](@ref)``; do not repeat the main explanation.
- Show the next operation, including in constructor overloads. Extended
  workflows belong in the Manual or Tutorials, not in every docstring.

## Links between docstrings

The ``[`name`](@ref)`` links become hyperlinks in this documentation and are
checked on every build: a link to a name without a docstring fails the build.
In the REPL they are shown as text, so write them where they help the reader
move from inputs to a meson, a spectrum, a state and its wavefunction.

`scripts/audit_doc_links.jl` prints the link graph between the core entries
of that workflow. Entries with no incoming links are hard to discover; entries
with no outgoing links are dead ends.

```bash
julia scripts/audit_doc_links.jl
```

## Exported names

Every exported name needs a docstring; the documentation build lists missing
ones. If a name is exported but not meant for users, reconsider the export.
