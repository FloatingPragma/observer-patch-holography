# Boundary interactions and source correlations

The selected polynomial state bounds its conditional regional correlations
and boundary-sector coherence by the actual interactions crossing the cut.
The comparison preserves the source's two regional marginals and its Gauss
pairing. It does not replace the preparation by a thermal state or discard
boundary charges. A finite polynomial correction is necessary even for a
noninteracting control. The resulting boundary bound has no universal
continuum area coefficient.

See [the proof](PROOF.md) and [source receipt](results.json).
Run `python -m model.particle_primary_boundary_correlation --check` and
`python -m pytest -q model/particle_primary_boundary_correlation -W error`.
