# OPH Paper Index

This directory is the canonical publication surface for the main OPH papers. Each paper is kept as a TeX source beside its built PDF.

For a first reading, begin with the short [compact
case](../extra/compact_proof_of_oph.pdf). Continue with *Observers Are All You
Need* for the observer interpretation and the two quantitative closures, then
use the two technical derivation papers for spacetime and gravity on one side
and Standard Model gauge structure on the other. The consensus, particle, and screen papers provide the finite
repair theorem, numerical continuations, and physical carrier architecture.
The focused [positive-chamber Koide paper](../extra/koide_identity_from_positive_c3_face_circulants.pdf)
isolates the exact face-circulant identity and conditional finite tracial
Gelfand--Naimark--Segal balance from the open physical family attachment,
phase, and numerical ratios.
The focused [finite de Sitter capacity paper](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf)
proves the smooth shock normalization, finite entropy and transfer laws, and
graph-spectrum identity, with the physical time-advance dictionary stated as
a separate conditional attachment.

| Paper | Role |
| --- | --- |
| [Recovering Observer Spacetime and Einstein Dynamics](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf) ([source](recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)) | Full celestial-screen, modular, event-manifold, stress, generalized-entropy, and Einstein chain |
| [Deriving Standard Model Gauge Structure](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf) ([source](deriving_standard_model_gauge_structure_from_observer_overlap_consistency.tex)) | Full categorical and finite-carrier gauge routes, conditional matter image, and physical global-form boundary |
| [Observers Are All You Need](observers_are_all_you_need.pdf) ([source](observers_are_all_you_need.tex)) | Observer-first synthesis and main account of the local $P$ and global $N$ closures |
| [Reality as a Consensus Protocol](reality_as_consensus_protocol.pdf) ([source](reality_as_consensus_protocol.tex)) | Finite repair, protected records, and quotient normal forms |
| [Deriving the Particle Zoo](deriving_the_particle_zoo_from_observer_consistency.pdf) ([source](deriving_the_particle_zoo_from_observer_consistency.tex)) | Particle carriers, hierarchy coordinates, declared flavor tests, and executable pipeline |
| [Federated Echosahedral Screen Microphysics](screen_microphysics_and_observer_synchronization.pdf) ([source](screen_microphysics_and_observer_synchronization.tex)) | Twelve-port source selector, central records, and observer synchronization |
| [Paradise as Fixed-Point Consensus](paradise_as_fixed_point_consensus.pdf) ([source](paradise_as_fixed_point_consensus.tex)) | Observer continuation and interpretation |

The shortest informal introduction is [A Compact Case for OPH](../extra/compact_proof_of_oph.pdf).
The flavor theorem is developed in [The Positive-Chamber Koide Identity for
Icosahedral Face Circulants](../extra/koide_identity_from_positive_c3_face_circulants.pdf).
The finite de Sitter theorem is developed in [The de Sitter Time-Advance Sign
from a Finite Screen with Fixed Capacity](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.pdf).
Focused papers are listed in the [supplement index](../extra/).

Shared TeX lives in [`tex_fragments/`](tex_fragments/). The [BFT/QECC appendix](appendix_B_bft_qecc_extensions.tex) is included by *Reality as a Consensus Protocol*.

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
manual review, run `python3 tools/bump_paper_release.py` followed by
`python3 tools/refresh_paper_release.py --publication` to prepare a candidate
under a new release identifier. Publication mode rebuilds the canonical book
and binds its bytes to that identifier in the same manifest.
