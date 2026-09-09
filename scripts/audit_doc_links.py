#!/usr/bin/env python3
"""Audit explicit @ref links in core API docstrings (not Julia call dependencies).

Run from any directory: python3 scripts/audit_doc_links.py
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
CORE = '''load_parameters load_quark_masses load_parameters_and_quark_masses
GIParameters QuarkMassTable ConstituentMasses Meson reduced_mass flavor_label
is_equal_flavor BasisState spectrum_levels
compute_spectrum MixedSpectrum MixedState StateMixing spectrum_state
radial_wave physical_components RadialWave sample_wave wave_norm
physical_state_amplitude compute_isoscalar_spectrum'''.split()
DOC = re.compile(r'"""(.*?)"""', re.S)
DECL = re.compile(r'\s*(?:#[^\n]*\n\s*)*(?:(?:mutable\s+)?struct\s+|abstract\s+type\s+|function\s+|const\s+)?([A-Za-z_]\w*)')
LINK = re.compile(r'\[`([^`]+)`\]\(@ref\)')


def collect():
    docs = {}
    for path in sorted((ROOT / 'src').glob('*.jl')):
        text = path.read_text()
        for match in DOC.finditer(text):
            declaration = DECL.match(text, match.end())
            if declaration:
                docs.setdefault(declaration[1], []).append(match[1])
    return {name: set(LINK.findall('\n'.join(bodies))) - {name} for name, bodies in docs.items()}


def metrics(graph):
    edges = {(a, b) for a in CORE for b in graph.get(a, set()) if b in CORE}
    incoming = {n: sum(b == n for a, b in edges) for n in CORE}
    outgoing = {n: sum(a == n for a, b in edges) for n in CORE}
    return edges, incoming, outgoing


def reachable(graph, start):
    seen, todo = {start}, [start]
    while todo:
        for target in graph.get(todo.pop(), set()) & set(CORE) - seen:
            seen.add(target)
            todo.append(target)
    return seen


def main():
    graph = collect()
    edges, incoming, outgoing = metrics(graph)
    lines = ['# Public API documentation graph', '',
             f'Explore {len(CORE)} public entries along the input → spectrum → wavefunction workflow.',
             'Arrows lead from a help entry to a related entry linked in its docstring.',
             'Start with the suggested reading paths, or expand the complete graph.', '']
    route = [
        ('QuarkMassTable', 'load_quark_masses'),
        ('load_quark_masses', 'load_parameters_and_quark_masses'),
        ('load_parameters', 'load_parameters_and_quark_masses'),
        ('load_parameters_and_quark_masses', 'Meson'),
        ('ConstituentMasses', 'Meson'), ('Meson', 'reduced_mass'),
        ('Meson', 'flavor_label'), ('Meson', 'is_equal_flavor'),
        ('Meson', 'spectrum_levels'),
        ('spectrum_levels', 'BasisState'), ('Meson', 'compute_spectrum'),
        ('compute_spectrum', 'MixedSpectrum'), ('MixedSpectrum', 'MixedState'),
        ('MixedSpectrum', 'spectrum_state'), ('spectrum_state', 'physical_components'),
        ('MixedState', 'StateMixing'), ('physical_components', 'RadialWave'),
        ('physical_components', 'physical_state_amplitude'),
        ('spectrum_state', 'radial_wave'), ('radial_wave', 'sample_wave'),
        ('radial_wave', 'wave_norm'),
        ('BasisState', 'compute_isoscalar_spectrum'),
    ]
    lines += ['', '## Suggested reading paths', '',
              'A selected subset of actual links makes the main workflow easier to read.',
              '', '```mermaid', 'flowchart LR']
    lines += [f'    {a} --> {b}' for a, b in route if (a, b) in edges]
    lines += ['```', '', '## Complete graph', '',
              '<details>', '<summary>Show all core documentation links</summary>', '',
              '```mermaid', 'flowchart LR']
    lines += [f'    {n}["{n}"]' for n in CORE]
    lines += [f'    {a} --> {b}' for a, b in sorted(edges)]
    lines += ['```', '', '</details>', '', '## Link coverage', '',
              f'The selected entries contain **{len(edges)}** links. From `Meson`,',
              f'**{len(reachable(graph, "Meson"))} of {len(CORE)}** entries are reachable (including itself).', '',
              '| Entry | Incoming links | Outgoing links |', '|---|---:|---:|']
    lines += [f'| `{n}` | {incoming[n]} | {outgoing[n]} |' for n in CORE]
    lines += ['',
              'Counts combine constructor and method docstrings and exclude self-links.',
              'This graph covers explicit docstring links among the selected entries;',
              'it does not represent function calls or every exported method.', '',
              'See [Discoverability](discoverability.md) for help conventions and regeneration.', '']
    output = ROOT / 'docs' / 'discoverability_graph.md'
    output.write_text('\n'.join(lines))
    print(f'{output}: {len(edges)} core links')

if __name__ == '__main__':
    main()
