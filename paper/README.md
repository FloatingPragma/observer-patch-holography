# OPH Paper Index

This directory is the canonical publication surface for the main OPH papers. Each paper is kept as a TeX source beside its built PDF.

The papers distinguish a proposed complete-history interpretation from finite
repair theorems and conditional physical reconstruction. Internal computation
and timeless closure describe the same proposed structure; they do not equate
a solver counter, causal precedence, modular parameter or local clock reading.
The [working terms](../docs/INTERPRETATION_AND_TERMS.md) summarize these
definitions. Mathematical background supplied by a benchmark remains an
explicit premise even when the programme aims to reconstruct its physical role.

The papers form one tower over the executable model, and the
[paper set](../docs/PAPER_SET.md) maps each paper to its layer. Every paper
opens with the same statement of the
[executable-model hypothesis](../docs/OPH_EXECUTABLE_HYPOTHESIS.md) and then
its own hypothesis, what it proves, the open join it needs from the model and
what would refute it. The [basis](../docs/BASIS.md) states the closure
hypothesis behind the whole set.

For a first reading, begin with the flagship synthesis,
[*Observers Who Read Their Universe into Existence*](../flagship/from_observer_consensus_to_standard_physics.pdf),
which walks through the tower layer by layer. The construction paper,
[*An Executable Observer for Observer Patch Holography*](executable_observer_construction.pdf),
states the specification, the executed components and the ordered joins.

The current-source previews also include the
[28 October finite C/E campaign](../model/completion/README.md): continuing
record-buffer operation with accounted exports, an exact stationary reference
lifted to the actual complete records, and scoped branch-charge and modular-
access obstructions. This finite campaign alone grants no integrated gate acceptance.
The construction paper gives the detailed observer/reference join; the shared
summary in every primary paper preserves the same scope. The current previews
also include the [same-action bound-pair construction](../model/bound_pair/README.md),
its complete native clock corner, regional CAR maps and finite energy/entropy
identities. One composite has verified motion, spectral protection,
charge and complete spectral records. Its band is topologically trivial,
physical species and units remain unassigned, and the same-speed relativistic
mass relation fails. The interaction also prevents exact isolation of one
Weyl momentum pocket by deleting the rest. The papers state these results
alongside the remaining field and Einstein premises.

The [particle paper](deriving_the_particle_zoo_from_observer_consistency.pdf)
gives certified finite interacting spectra, their conserved radial currents,
a memory realization of the pixel-dependent charged-lepton operator, its
localized native-current inclusion, and a lattice texture energy bound.
The [geometry paper](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf)
states the corresponding conditional stress and generalized-volume entropy
relation. The [source particle map](../model/completion/SOURCE_PARTICLE_MAP.md)
connects each calculation to its readable inputs and mathematical domain.
A common physical defect assignment, interacting continuum and Einstein
realization require the stated identification and limit theorems.

| Paper | Role |
| --- | --- |
| [Observers Who Read Their Universe into Existence](../flagship/from_observer_consensus_to_standard_physics.pdf) ([source](../flagship/from_observer_consensus_to_standard_physics.tex)) | Flagship synthesis: the closure and executable-model hypotheses, adopted principles and tower layer by layer, with the status of each layer; released on [PhilPapers](https://philpapers.org/rec/MUEFOC) |
| [An Executable Observer for Observer Patch Holography](executable_observer_construction.pdf) ([source](executable_observer_construction.tex)) | Layers 2–4 as executed: the specification, executed components, premise map, and the ordered joins with their stopping tests |
| [Reality as a Consensus Protocol](reality_as_consensus_protocol.pdf) ([source](reality_as_consensus_protocol.tex)) | Layer 1: local agreement repair, quotient normal forms, provenance order, and cycle holonomy as the residue repair cannot remove |
| [Observers Are All You Need](observers_are_all_you_need.pdf) ([source](observers_are_all_you_need.tex)) | Layer 3: observers, records, quantum and thermodynamic identities, and the closure equations |
| [Recovering Observer Spacetime and Einstein Dynamics](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) ([source](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)) | Layers 4 and 6: clocks, causal order and geometry from the same process, and the conditional Einstein implication with its realization premises |
| [Standard Model Gauge Structure under Complete-Response and Observer-Transport Premises](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) ([source](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex)) | Layer 5, gauge half: gauge structure as overlap transport and the compact Lie type under complete-response premises |
| [Particles as Localized Topological Defects of Native Overlap Transport](deriving_the_particle_zoo_from_observer_consistency.pdf) ([source](deriving_the_particle_zoo_from_observer_consistency.tex)) | Layer 5, matter half: the localized-defect hypothesis, the native decision test and the mass obstructions |
| [Twelve-Port Icosahedral Carriers and Native Seam Repair](screen_microphysics_and_observer_synchronization.pdf) ([source](screen_microphysics_and_observer_synchronization.tex)) | Layer 2: the twelve-port carrier, seam rule, record laws and federation |

The flavor theorem is developed in [The Positive-Chamber Koide Identity for
Icosahedral Face Circulants](../extra/koide_identity_from_positive_c3_face_circulants.pdf),
and the finite de Sitter theorem in [The de Sitter Time-Advance Sign from a
Finite Screen with Fixed Capacity](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf).
Focused papers are listed in the [supplement index](../extra/). Papers that
left the active set keep their sources and last built PDFs in
[`archive/`](archive/); the [paper set](../docs/PAPER_SET.md) lists them.

Shared TeX lives in [`tex_fragments/`](tex_fragments/). The [BFT/QECC appendix](appendix_B_bft_qecc_extensions.tex) is included by *Reality as a Consensus Protocol*.

The [fine-structure supplement](../extra/fine_structure_constant_derivation.pdf)
retains the certified trial-map roots and separately reconstructs the
[entangled-qudit candidate of A. D. Wissner-Gross](https://alexwg.org/publications/ChainTheoryAlpha.pdf).
His exchange construction and exact scalar receive direct attribution.
The OPH complete-outcome, support-alignment and normalized-readout controls
retain the independent physical coupling and same-cell premises.

## Reproducibility

The shared release identifier lives in [`release_info.tex`](release_info.tex).
[`paper_release_manifest.json`](paper_release_manifest.json) records the paper
artifacts and the builder-stamped canonical book receipt.

From the repository root:

```bash
python3 tools/refresh_paper_release.py --preview
```

This rebuilds review PDFs and their manifest without requiring a version bump,
including when the visible identifier belongs to an existing release. After
corrections, main-branch PDFs can therefore be current-source previews with
the previous visible identifier; tagged or deposited release artifacts retain
their original bytes. After
an explicit maintainer decision under the release policy and manual review,
run `python3 tools/bump_paper_release.py` followed by
`python3 tools/refresh_paper_release.py --publication` to prepare a candidate
under a new release identifier. Publication mode rebuilds the canonical book
and binds its bytes to that identifier in the same manifest.
