# OPH Focused Papers And Supplements

The main reading route lives in [`paper/`](../paper/). This directory contains the short informal case and focused papers that develop one mathematical, physical, computational, or interpretive branch in depth.

## Start Here

- [A Compact Case for OPH](compact_proof_of_oph.pdf) ([source](compact_proof_of_oph.tex)) is the seven-page, reader-friendly summary of the theory’s strongest converging results.

## Mathematical Foundations

- [Observation-Determined Normal Forms](observable_normal_forms.pdf) ([source](observable_normal_forms.tex), [bibliography](observable_normal_forms.bib)) develops stability, obstructions, and refinement for constraint and rewrite systems.
- [Machine-Checked Finite Quantum Event Algebras](machine_checked_finite_event_algebras.pdf) ([source](machine_checked_finite_event_algebras.tex), [bibliography](machine_checked_finite_event_algebras.bib)) connects the finite event surface to its Lean artifact.
- [Explaining the Yang–Mills Mass Gap with Observer-Patch Repair Dynamics](yang_mills_gap_clay_problem.pdf) ([source](yang_mills_gap_clay_problem.tex)) develops the fixed-cutoff gap and continuum-transfer program.

## Quantitative And Physical Branches

- [The Fine-Structure Constant as an OPH Pixel Fixed Point](fine_structure_constant_derivation.pdf) ([source](fine_structure_constant_derivation.tex)) develops the local closure calculation and its certificates.
- [Theoretical Bounds on $\chi_\nu$](chi_nu_susceptibility_bounds.pdf) ([source](chi_nu_susceptibility_bounds.tex), [collar-survival correction](chi_nu_collar_survival_presence_correction.md)) develops the coherent-matter susceptibility bounds and evidence conditions.
- [Photonic Fixed-Point Consensus for SHA-256d Proof of Work](Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.pdf) ([source](Photonic_fixed-point_consensus_for_SHA-256d_proof_of_work.tex)) develops the optical constraint-and-repair architecture.
- [Observer-Patch Holography as a String-Vacuum Selector](observer_patch_holography_as_string_vacuum_selector.pdf) ([source](observer_patch_holography_as_string_vacuum_selector.tex)) develops the conditional string-sector selection program.

## Observer And Engineering Interpretations

- [Thinking as Patch-Net Fixed-Point Search](thinking_as_patch_net_fixed_point_search.pdf) ([source](thinking_as_patch_net_fixed_point_search.tex)) applies the observer-patch architecture to cognition and learning.
- [Hacking the Simulation: The Anti-Gravity Exploit](hacking-the-simulation-anti-gravity-exploit.pdf) ([book source](hacking-the-simulation-anti-gravity-exploit/)) presents the dark-gravity engineering route and its experimental evidence protocol.

All root-level TeX papers in this directory can be rebuilt from the repository root with:

```bash
python3 tools/build_tex_papers.py --extra-only
```
