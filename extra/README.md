# OPH Focused Papers

The main reading route lives in [`paper/`](../paper/). This directory contains focused TeX papers that develop one mathematical, physical, computational, or interpretive branch in depth. Research notes written in Markdown live in [`docs/research/`](../docs/research/).

The papers sit in one tower over the [executable-model hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md). The executable model runs the OPH law, specified in the [simulator specification](../simulator/SPEC.md): each constituent moves by reflected routing, two constituents that meet bind into a pair through a contact law whose two rails form a clock, a capture vertex forms that pair and sends away a marker, and a newborn pair opens its hop over 128 updates. Observers who inhabit the simulation can describe their experiences as relativistic 3+1D continuum spacetimes. Each layer proves finite theorems under displayed premises, and a conclusion passes to the next layer when that premise is realized on the same model process. The [model page](../model/README.md) lists what is proven and measured on the law; the other premises are work in progress.

## Layer 1: repair mathematics

These two papers are standalone mathematics. They are tools of the layer in which local agreement repair reaches global consistency, and they assume no physical interpretation of their objects.

- [Observation-Determined Normal Forms](observable_normal_forms.pdf) ([source](observable_normal_forms.tex), [bibliography](observable_normal_forms.bib)) proves when a settled state is fixed by protected observations, how stably under refinement, which schedulers and attempt budgets reach it, and when no local repair can construct it.
- [Verified Projection-Event Calculus in Lean 4](machine_checked_finite_event_algebras.pdf) ([source](machine_checked_finite_event_algebras.tex), [bibliography](machine_checked_finite_event_algebras.bib)) is a Lean library of finite projection events, partition expectations, Lüders updates and CHSH bounds used by the record and instrument interfaces.

The [working terms](../docs/INTERPRETATION_AND_TERMS.md) distinguish the proposed complete-history interpretation from repair fixed points, stationary laws and local operational time. An operational observer need not be conscious.

## Layers 5 to 7: conditional physical supplements

These papers carry conditional branches, diagnostics and retrospective comparisons. Each states its hypothesis, what it proves unconditionally, and what the executable model must realize before the physical reading applies. None is a same-process derivation from the executable model.

- [A Conditional Identification of the Yang–Mills Gap with a Finite Repair Gap](yang_mills_gap_clay_problem.pdf) ([source](yang_mills_gap_clay_problem.tex)) supplements the gauge layer. It proves a certified relaxation gap for a finite collar-repair generator under a Dobrushin influence bound and tests the transfer receipt that would identify it with the pure Yang–Mills gap: on \(\mathbb Z_2\) lattice gauge theory the receipt holds in the free limit and fails at every tested interacting coupling. A positive repair gap or finite level spacing is not a physical mass gap, the full theory with a massless photon has no positive global vacuum gap, and the continuum premises a Clay-level construction requires stay conditional.
- [The Positive-Chamber Koide Identity for Icosahedral Face Circulants](koide_identity_from_positive_c3_face_circulants.pdf) ([source](koide_identity_from_positive_c3_face_circulants.tex)) supplements the matter layer. It proves the exact positive-eigenvalue identity \(Q=1/3+(2/3)(|b|/a)^2\) for Hermitian \(C_3\) face circulants and the conditional finite tracial Gelfand–Naimark–Segal balance.
- [The de Sitter Time-Advance Sign from a Finite Screen with Fixed Capacity](de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf) ([source](de_sitter_time_advance_sign_from_fixed_screen_capacity.tex)) supplements the gravity and cosmology layers. It proves the pure de Sitter shock normalization, the finite entropy maximum, the uniform capacity-transfer law, the analytic curvature and the line-graph spectrum identity. The time-advance sign and the fixed-capacity equation of state are conditional on fixed screen capacity and an explicit horizon-area dictionary that the executable model must realize.
- [The Fine-Structure Constant as an OPH Pixel Fixed Point](fine_structure_constant_derivation.pdf) ([source](fine_structure_constant_derivation.tex)) is a conditional layer-seven diagnostic. Two declared closure maps have interval-certified fixed points; they use declared closure coordinates and a gauge-running convention, and with the measured hadronic spectral input in place of the screening ansatz the endpoint excludes the measured value. The paper also reconstructs [A. D. Wissner-Gross's entangled-qudit candidate](https://alexwg.org/publications/ChainTheoryAlpha.pdf). The original spectral construction is credited to him; its physical electromagnetic identity and proposed same-cell attachment are conditional.

Archived manuscripts, including the string-vacuum selector and the patch-net model of cognition, are kept with their last built PDFs in [`paper/archive/`](../paper/archive/) and listed in the [paper set](../docs/PAPER_SET.md).

All root-level TeX papers in this directory can be rebuilt from the repository root with:

```bash
python3 tools/build_tex_papers.py --extra-only
```
