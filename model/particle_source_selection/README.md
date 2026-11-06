# Family phases from observer records

Exact native event rates give the numbers \(\sqrt2\) and \(2/9\), but the
complete classical record law does not determine a coherent phase. The
[proof](PROOF.md) gives an exact family of indistinguishable classical
instruments whose coherent record erasures have different spectra. Cyclic
averaging of the native real family kernel also forces a repeated eigenvalue.
A physical family spectrum requires a source-selected phase-sensitive
preparation, observation or dynamics beyond these inputs.

The executable reconstructs native weights, checks all twenty inherited
three-cycles and verifies both counterfamily members on the entire matrix
space. No particle comparison data are imported.

```sh
python3 -m model.particle_source_selection.selection --check
python3 -m pytest -q model/particle_source_selection
```

[`results.json`](results.json) binds the source and contains the exact replay.
Omit `--check` to regenerate it after a reviewed source change.
