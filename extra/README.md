# OPH Focused Papers And Supplements

The main reading route lives in [`paper/`](../paper/). This directory contains focused papers that develop one mathematical, physical, computational, or interpretive branch in depth.

## Mathematical Foundations

- [Observation-Determined Normal Forms](observable_normal_forms.pdf) ([source](observable_normal_forms.tex), [bibliography](observable_normal_forms.bib)) develops stability, obstructions, and refinement for constraint and rewrite systems.
- [Finite Quantum Event Algebras](machine_checked_finite_event_algebras.pdf) ([source](machine_checked_finite_event_algebras.tex), [bibliography](machine_checked_finite_event_algebras.bib)) develops the finite event surface.
- [Explaining the Yang–Mills Mass Gap with Observer-Patch Repair Dynamics](yang_mills_gap_clay_problem.pdf) ([source](yang_mills_gap_clay_problem.tex)) develops the fixed-cutoff gap and continuum-transfer program.

## Quantitative And Physical Branches

- [The Positive-Chamber Koide Identity for Icosahedral Face Circulants](koide_identity_from_positive_c3_face_circulants.pdf) ([source](koide_identity_from_positive_c3_face_circulants.tex)) proves the exact positive-eigenvalue identity \(Q=1/3+(2/3)(|b|/a)^2\) and the conditional finite tracial Gelfand–Naimark–Segal balance. Physical charged-family attachment, phase, and numerical mass ratios are open; the target-informed numerical near-match is diagnostic.
- [The de Sitter Time-Advance Sign from a Finite Screen with Fixed Capacity](de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf) ([source](de_sitter_time_advance_sign_from_fixed_screen_capacity.tex)) proves the pure-de-Sitter shock normalization, finite entropy maximum, uniform capacity-transfer law, analytic curvature, and line-graph spectrum identity. The time-advance interpretation is conditional on an explicit physical dictionary.
- [The Fine-Structure Constant as an OPH Pixel Fixed Point](fine_structure_constant_derivation.pdf) ([source](fine_structure_constant_derivation.tex)) develops the local closure calculation and its certificates.
- [Photonic Fixed-Point Consensus for SHA-256d Proof of Work](Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf) ([source](Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.tex)) develops the optical constraint-and-repair architecture.
- [Observer-Patch Holography as a String-Vacuum Selector](observer_patch_holography_as_string_vacuum_selector.pdf) ([source](observer_patch_holography_as_string_vacuum_selector.tex)) develops the conditional string-sector selection program.

## Observer And Engineering Interpretations

- [Thinking as Patch-Net Fixed-Point Search](thinking_as_patch_net_fixed_point_search.pdf) ([source](thinking_as_patch_net_fixed_point_search.tex)) applies the observer-patch architecture to cognition and learning.

All root-level TeX papers in this directory can be rebuilt from the repository root with:

```bash
python3 tools/build_tex_papers.py --extra-only
```
