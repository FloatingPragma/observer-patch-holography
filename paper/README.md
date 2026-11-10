# OPH Paper Index

This directory is the canonical publication surface for the main OPH papers. Each paper is kept as a TeX source beside its built PDF.

The papers distinguish a proposed complete-history interpretation from finite
repair theorems and conditional physical reconstruction. Internal computation
and timeless closure describe the same proposed structure; they do not equate
a solver counter, causal precedence, modular parameter or local clock reading.
The [working terms](../docs/INTERPRETATION_AND_TERMS.md) summarize these
definitions. Mathematical background supplied by a benchmark remains an
explicit premise even when the programme aims to reconstruct its physical role.

The [paper set](../docs/PAPER_SET.md) maps the logical dependencies. Each
paper states its own construction, assumptions, proofs and physical scope.
Shared project summaries are kept in the index; specialist papers develop
their mathematical subject without repeating the full program.

Begin with the flagship synthesis,
[*Observers Who Read Their Universe into Existence*](../flagship/from_observer_consensus_to_standard_physics.pdf).
The construction paper,
[*An Executable Observer for Observer Patch Holography*](executable_observer_construction.pdf),
states the repair and memory specification, complete observer records and
controlled field limits. The [basis](../docs/BASIS.md) distinguishes closure,
operational agreement, and the adopted state and response principles.

The [particle paper](deriving_the_particle_zoo_from_observer_consistency.pdf)
gives the native transfer and charged spectra, the complete conserved
quadratics, the pixel-dependent mass law and its memory realization, the common
localized current/mass operator, the compact finite-energy defect pair and its
quartic interaction, the twelve mass comparison values and the exact exclusions.
It takes the gauge type and matter module from the gauge paper.
The selected Wilson source also supplies an isolated finite gauge transition,
regional charged spin and pair-energy records, a controlled energy-current
response and a capacity-dependent bound on measured energy moments. The gauge
line loses original weight under refinement; these results do not identify
physical particle masses. The [Markdown dependency map](../model/particle_joint_reconstruction/DEPENDENCIES.md)
and [research account](../docs/research/OPH_NATIVE_PARTICLE_MASS_JOIN.md)
state their source domains and remaining physical requirements.
The [geometry paper](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf)
states the corresponding conditional stress and generalized-volume entropy
relation. The [source particle map](../model/completion/SOURCE_PARTICLE_MAP.md)
connects each calculation to its readable inputs and mathematical domain.
A common physical defect assignment, interacting continuum and Einstein
realization require the stated identification and limit theorems.

The flagship summarizes this construction alongside the strongest controlled
field limits and conditional Einstein theorem, with specialist references
for their proofs. The core's Gibbs response uses the same Hamiltonian as its
interacting poles. Distinct texture carriers and a supplied response kernel
do not establish a common physical spectrum or geometric entropy law.

| Paper | Role |
| --- | --- |
| [Observers Who Read Their Universe into Existence](../flagship/from_observer_consensus_to_standard_physics.pdf) ([source](../flagship/from_observer_consensus_to_standard_physics.tex)) | Flagship synthesis: the OPH proposal, finite theorem stack, source specification, physical reconstruction and status of each implication; [PhilPapers entry](https://philpapers.org/rec/MUEFOC) |
| [An Executable Observer for Observer Patch Holography](executable_observer_construction.pdf) ([source](executable_observer_construction.tex)) | Layers 2–4: repair and memory, complete observer records, clocks, field limits and physical identification conditions |
| [Reality as a Consensus Protocol](reality_as_consensus_protocol.pdf) ([source](reality_as_consensus_protocol.tex)) | Layer 1: local agreement repair, quotient normal forms, provenance order, and cycle holonomy as the residue repair cannot remove |
| [Observers Are All You Need](observers_are_all_you_need.pdf) ([source](observers_are_all_you_need.tex)) | Layer 3: observers, records, quantum and thermodynamic identities, and the closure equations |
| [Recovering Observer Spacetime and Einstein Dynamics](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) ([source](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)) | Layers 4 and 6: clocks, causal order and geometry from the same process, and the conditional Einstein implication with its realization premises |
| [Standard Model Gauge Structure under Complete-Response and Observer-Transport Premises](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) ([source](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex)) | Layer 5, gauge half: gauge structure as overlap transport and the compact Lie type under complete-response premises |
| [Particles from Native Overlap Transport: Common Currents, Defect Modes and a Pixel-Dependent Mass Operator](deriving_the_particle_zoo_from_observer_consistency.pdf) ([source](deriving_the_particle_zoo_from_observer_consistency.tex)) | Layer 5, matter half: native charged spectra, common gauge currents, defect modes, the pixel-dependent mass operator, the mass comparison and its exclusions |
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
