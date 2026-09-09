# Discoverability

GIModel help should answer: **what is this, how do I get it, and what can I do
with it next?** In Pluto, related links should let a reader move naturally from
inputs to a meson, a spectrum, a state, and its wavefunction components.

- Start with purpose, signature, and units. Show real property access.
- Keep example blocks to three or four lines, each showing one action.
  Longer lists of independent inspection methods are fine.
- State prerequisites and link to their producers. Use `pkgdir(GIModel)` for
  bundled data paths so examples work from any notebook directory.
- Use one `## Related` section with relevant links and short descriptions.
  Use ``[`name`](@ref)`` for links; avoid repeating the main explanation.
- Show the next operation, including in constructor overloads. Keep extended
  workflows in examples rather than repeating them in every docstring.

[The documentation graph](discoverability_graph.md) maps explicit links among
key public entries. Missing incoming links reveal hard-to-find entries; missing
outgoing links reveal dead ends. Link relevance matters more than link count.
The graph measures source links, not Pluto rendering or example correctness.

Regenerate the graph from the repository root:

```sh
python3 scripts/audit_doc_links.py
```

Run edited examples and check live help when changing documentation syntax.
