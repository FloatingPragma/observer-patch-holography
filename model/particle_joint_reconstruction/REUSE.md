# Reusable construction results and their domains

The governing acceptance conditions are in [dependencies.json](dependencies.json).
Each result is used on its actual source, state, reference and clock. A
theorem about a different program supplies a method or a comparison;
transferring its physical conclusions needs a proved map.

| Material | Reusable statement | Boundary of the reuse |
| --- | --- | --- |
| [Primary contract](../particle_common_contract/PROOF.md), [physical boundary state](../particle_physical_boundary_state/PROOF.md), [locality](../particle_primary_locality/PROOF.md) | Exact finite preparation and complete local records on one Wilson/Gauss family, with a conservative uniform influence cone. | No physical Lorentz transformations or species dictionary. |
| [Finite ground tail](../../Lean/QFT/GroundRecordTail.lean) and [positive-energy consumer](POSITIVE_ENERGY.md) | Original negative-transfer weights vanish for uniformly bounded physical probes along the selected cooling schedule. Time filtering gives uniformly continuous response families. | A single response limit can be zero and does not define the full field algebra. |
| [Weyl/contact continuum](../continuum_field_net/README.md) | Full moment hierarchy, GNS construction and finite propagation on a separately declared CAR family. | Its state and action do not replace the primary gauge state. The [centered branch](../continuum_free_limit/README.md) is free and massless. |
| [Selected-state GNS](../../Lean/QFT/ColimitSelectedStateGNS.lean) and [representation](../../Lean/QFT/ColimitGNSRepresentation.lean) | Mathematical representation interfaces after a compatible algebra and state have been supplied. | Finite type correctness does not supply physical vacuum convergence, noncollapse or Lorentz covariance. |
| [Spectral window access](../source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md) | Complete Ramsey scores, original spectral weights, finite coherent services and explicit approximation errors. | Its charged CAR experiment has its own source and ground. A window containing zero does not certify an excitation above zero. |
| [Energy/work identities](../particle_gravity_response/PROOF.md), [regional modular support](../particle_spatial_modular_response/PROOF.md) | Operator continuity, physical support and finite marginal entropy response. | Marginal entropy differs from the entropy of the complete apparatus/record output; geometric modular kernels are separate. |
| [Matched Einstein limit](../matched_einstein_limit/PROOF.md) | Ten-probe tensor inversion and an explicit Einstein residual from matched regional entropy, area and stress estimates. | Proper geometry and common errors small relative to the fourth power of the radius are hypotheses. |
| [Actual finite survival](../particle_primary_survival/PROOF.md) | A nonzero primary plaquette response and finite-filter positive band, using the actual polynomial passivity. | Times and bands depend on level and the lower weight vanishes; rare pure-gauge sectors cannot be normalized into a surviving source spectrum. |
| [Regular spatial algebra](../particle_primary_spatial_algebra/PROOF.md) | Continuous spatial translations commuting with source time, with actual finite-state bulk word transfer. | The chosen state's spatial invariance and canonical GNS momentum are unproved. |
| [Excitation transport](../particle_primary_mass_readout/PROOF.md) | Original transfer and final excitation laws have transport cost at most C² epsilon; declared SI conversion is immediate. | Physical rest shells, particle identities and an independent mass prediction require their own inputs. |
| [Primary capacity](../particle_primary_capacity/PROOF.md) | Exact degree-one ambient capacity and physical dimension bounds from local singlets. | Prepared entropy, correctable records, apparatus resources and cosmological capacity are different quantities. |
| [Conditional particle action](../conditional_standard_model/PROOF.md), [candidate](CANDIDATE.md) | Exact representations, invariant couplings and calibrated quadratic spectrum. | Chiral source attachment, physical poles, neutrino matching and independent predictions have their stated premises. |
| [Dense-correlation gap criterion](../particle_yang_mills_gap_contract/PROOF.md) | A continuum spectral exclusion from dense vacuum correlation decay and finite-energy noncollapse. | The photon-containing full theory and a pure compact-simple Yang-Mills source are different domains. |

The external-workspace `pro/Particles/REUSE_MAP.md` indexes the private
particle calculations and their public counterparts. Its recorded scan
covers 3,833 text files and identifies 77 selected candidates; this is a
file/search inventory, not a mathematical verification of every file.
Round 05's interacting Higgs and chiral-electric controls are represented
in the public packages. The separate private copies supply no independent
proof of a common interacting relativistic limit. Archived inverse mass
and Koide calculations retain their observable and calibration premises.

The standard [scaling-algebra framework of Buchholz and Verch](https://arxiv.org/abs/hep-th/9501063)
provides a useful precedent for separating existence of limit states from
their physical classification. Its hypotheses start from a local quantum
field theory; they do not establish those hypotheses for the finite source.
[Lüscher's lattice transfer-matrix construction](https://bib-pubdb1.desy.de/record/396349)
proves positivity for its specified Wilson lattice theory. Applying it to
a changed regulator or declaring a continuum mass gap requires the
corresponding identification and limiting estimates.

The productive order follows the shared algebra and its actual apparatus:
construct complete filtered records; take one state on all their words;
test nonzero positive-energy response with original weights; then attach
spatial covariance, interactions and species. Conserved source currents retain their finite theorem; gravity is deferred
while the active work resolves physical mass outputs and capacity dependence.
