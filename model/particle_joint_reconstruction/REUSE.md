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
| [Finite spatial Fourier instrument](../particle_primary_momentum/PROOF.md) | Complete physical aperture records, original spatial Gram weights and a reference-uniform error bound. | Open-edge terms and a wavevector label do not establish conserved particle momentum. |
| [Native channel labels](../particle_primary_channel_labels/PROOF.md) | Five conserved occupations with a physical hypercharge relation; a distinct cubic A1 observable has positive finite weight. | An invariant mixed state is not a scalar vacuum; cubic A1 is not a unique continuum-spin or species label. |
| [Source spectral scale](../particle_primary_spectral_scale/PROOF.md) | Exact decoupled gap and local perturbative packets; conditional fixed-band exclusion from a true uniform interacting gap. | The interacting endpoint gap, isolated branch, perturbation remainder and mass spectrum are unproved. Finite comparison and continuum claims have different error domains. |
| [Endpoint lower bound](../particle_primary_endpoint_bounds/PROOF.md), [Ritz upper bound](../particle_primary_variational_readout/PROOF.md) | Full-strength interacting ground enclosure, with a full-space residual and evaluated branch-isolation test. | The residual is broad and the retained-penalty isolation test fails; no isolated excitation or original probe residue is inferred. |
| [Joint spectroscopy](../particle_primary_joint_spectroscopy/PROOF.md) | Complete two-energy/address records, original positive broad band and same-member ground interval adapters. | Initial energy measurement disturbs the retained purification; a broad wavevector band is not a resolved relativistic mass branch. |
| [Occupation-sensitive electric bound](../particle_primary_occupation_bound/PROOF.md) | Member-zero ground floor and selected mean beta-defect count below5; averaged neutral-pair success at least11/16. | The reference enters an operator inequality, not the preparation. The retained penalty used for this weight fails its next-member admission. |
| [Exact local residual moments](../particle_primary_residual_refinement/PROOF.md) | Parent Haar/CAR moments bound the unchanged trial residual by5.231016 on member zero, retaining the full omitted space. | A residual enclosure does not isolate a level or assign the trial's spectral weights to the selected state. |
| [Original neutral response](../particle_primary_response_moments/PROOF.md) | Local non-Hermitian pair probes have frequency-six residual squared at most3 on member zero; original transfer weight at least1/2 in [2,10] and final excitation weight at least3/8 in [2,21/2], before clock errors. | Broad gauge-neutral, even response; no selected physical momentum, dispersion or particle name. |
| [Cubic/parity pair records](../particle_primary_pair_channels/PROOF.md) | Exact full-Fock pair components and a complete coherent apparatus have original positive-band weight in every member-zero cubic/parity isotypic channel and all32 spatial/internal records. | These are operator representations and aperture labels, without continuum spin, physical momentum, an isolated branch or a rest mass. |
| [Scalable collective response](../particle_primary_refinement_response/PROOF.md) | A physical tiled trial and block-sensitive electric estimate give selected defect density below97/10, original transfer weight at least3/100 and final excitation weight at least1/50 throughout the finite family. | The lower energy endpoint shrinks and the upper grows; fixed physical-energy noncollapse and capacity-independent particle masses are unproved. |
| [Polynomial modular work](../particle_primary_gravity_premises/PROOF.md) | The actual preparation's modular slope sandwiches source work and bounds regional marginal entropy response; an actual neutral flip has work above18 at member zero. | Proper geometry and area are separate inputs. Retained pure-reference relative entropy differs, and direct periodic-coframe transfer changes the open flat action. |
| [Concentrated pair and odd trial](../particle_primary_branch_resolution/PROOF.md) | Original member-zero pair weight at least3/20 in [4,8], negative real correlation at time1/2, and a physical twelve-packet scalar-odd Rayleigh calculation. | The band permits several frequencies; the full-complement estimate does not isolate a pole or threshold. Trial upper-bound differences are not masses. |
| [Retained observation origin](../particle_primary_origin_momentum/PROOF.md) | Keeping one origin throughout each history gives a stationary translation-invariant relative limit and strongly commuting canonical H>=0 and P, with quantitative finite bulk and collar bounds. | The finite preparation and reference are unchanged. Fixed-origin invariance, nonzero momentum weight and Lorentz dispersion are not inferred; the regular Fourier transfer has its own proved scope below. |
| [Exact collective CAR moment](../particle_primary_energy_tightness/PROOF.md) | Original weight at least3/100 in [1/(100 beta_n),46/a_n] throughout the selected family; a direct selected-occupation estimate also populates the retained-origin protocol. The selected negative hopping density supplies an exact number-density f-sum. | The actual bare physical operator-norm residual grows as1/a. Fixed-energy survival needs controlled selected weighted response or a different probe, not merely a sharper constant. The tested volume-normalized density instrument loses every fixed positive-energy band. |
| [Open boundary strain](../particle_primary_open_coframe/PROOF.md) | A declared coframe reproduces the literal flat primary Hamiltonian, with six complete mean readouts and K_I=-H-H_W. | A declared metric coupling is not proper geometry, a Lorentz stress tensor, a diamond modular kernel or an area law. |
| [Odd spectral residual and coherent response](../particle_primary_odd_residual/PROOF.md) | Full-space residual at most 5.365163597 certifies some actual A1u excitation in [0.54794769,15.00864937]; a separate selected-state coherent probe has original positive-band weight. | The eigenvalue is not isolated; the probe residue is not assigned to it. Operator parity is not continuum spin or a rest frame. |
| [Charged endpoint pair](../particle_primary_charged_string/PROOF.md) | A quantum-link character joins local Q6=+6 and -6 increments, with original transfer weight above 3/10 in [39/4,135/4] and complete failure records. | The physical operator is even and globally neutral; no electromagnetic-charge, single-fermion, spin or Standard Model identification follows. |
| [Complete momentum transfer](../particle_primary_momentum_transfer/PROOF.md) | Regular coherent-offset instruments converge to sinc-squared kernels of the same relative H and P, with original weights and explicit clock/quadrature/bulk errors. | A finite bin is not a momentum projection. Nonzero joint residue and an effective state-subsequence rate are unproved. |
| [Bounded local Weyl response](../particle_primary_local_fluctuation/PROOF.md) | Actual phase operations have original probability mass one, a uniform absolute first energy moment and positive mean injected energy; retained bulk-origin records have a controlled complete distance. | Energy probability tightness allows a purely elastic limit. Fixed positive-band survival and this probe's momentum attachment require additional estimates. |
| [Prepared A1u spectral window](../particle_primary_spectral_resolution/PROOF.md) | Full residual below 4.947619 encloses some actual A1u excitation in [1.06662597,14.69223697]. Explicit preparation retains the original source/reference state in physical memory; the complete prepared-state energy clock reports an excited window with probability above 1/4. | Ready memory and state-changing controlled gates are explicit resources. No isolated level, original mixed-state residue, autonomous controller or rest mass is inferred. |
| [Coherent checkerboard response](../particle_primary_staggered_response/PROOF.md) | Actual Wilson cancellation gives original selected-state weight above 17/1280 in [3,9] and complete reported weight above 539/40960 in [95/32,289/32]. | A2g is the operator representation; the mixed input need not have that type. The finite pattern supplies no conserved physical momentum or continuum spin. |
| [Coherent charge-transfer meter](../particle_primary_charged_readout/PROOF.md) | Two finite couplings record admitted opposite native increments while preserving each original 49-outcome energy/probe history and its reference; radix-five counters are alias free on {-6,0,+6}. | Arbitrary charge differences, absolute initial-charge measurement, electromagnetic identification and particle species are different claims. |
| [Mixed-history energy forms](../particle_primary_current_tail/PROOF.md) | Uniform energy-root drift for fixed mixed phase/filtered words gives one subsequence with continuous real-time histories, nonnegative GNS energy and complete finite pointer limits. | No nonelastic survival, effective subsequence rate, growing-reference trace-norm limit or spatial momentum for the unfiltered fields. A physical anchor excludes a uniform absolute commutator-norm shortcut. |
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

[Yarotsky's relative-perturbation theorem](https://arxiv.org/abs/math-ph/0412040)
provides a possible route to a uniform lattice gap after a local dressing:
the perturbation must split into a relative part with bound below one and
a bounded part satisfying the theorem's quantitative smallness condition.
The primary source has no verified constants admitting that theorem at
interaction coefficient one. Its Gauss subspace and cutoff family also
require an explicit application of the theorem's tensor-product hypotheses.
The companion [quasiparticle theorem](https://arxiv.org/abs/math-ph/0411042)
requires an isolated reference level separated from multiparticle sums,
sufficiently small perturbation and translation invariance. A decoupled
gap, local electric positivity or a finite Rayleigh quotient alone does
not verify these premises or establish relativistic dispersion.

The productive order follows the shared algebra and its actual apparatus:
use the constructed complete phase/filtered histories and their common
positive-energy state; test nonzero positive-energy response with original
weights; attach same-channel spatial covariance, interactions and species. Conserved source currents retain their finite theorem; gravity is deferred
while the active work resolves physical mass outputs and capacity dependence.
