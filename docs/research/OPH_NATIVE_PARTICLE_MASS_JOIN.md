# Native particles, currents and mass operators

The [particle-mass response audit](OPH_PARTICLE_MASS_RESPONSE.md) retains the
source-bound evidence and its tested mass prescriptions. The
[source particle map](../../model/completion/SOURCE_PARTICLE_MAP.md) gives the
complete finite coefficient map and its physical boundaries. The
[completion program](../COMPLETION_PROGRAM.md) separates these results from
acceptance of a physical interacting theory.

## Common localized states, native currents and pixel-dependent masses

The [memory realization](../../model/mcpr_native_memory_bridge/PROOF.md) and
[defect/current construction](../../model/defect_current_bridge/PROOF.md)
provide an explicit common operator model. The carrier, independent probe
durations, feedback gain and coherent field service are chosen simulator
design. Their consequences are derived from that specification; uniqueness
under bare overlap consistency is not a requirement for this construction.

The retained pixel-dependent charged-lepton law determines three coefficients
`m_i(P)/v(P)`. The complete C-duty reset controller realizes a positive family
matrix `M(P)` with exactly those eigenvalues. Its coefficient has trace `p`;
the gain is `sum_i(m_i/v)/p`, so the calculation retains every disabled record.
Both pixel-source branches admit two independent-duration solutions within
the native gate domain. The forward branch reads no charged-lepton masses at
runtime. Its mass-law architecture has target-informed provenance, so this
is a conditional realization rather than an independent spectrum prediction.

On the odd radial cube, let `hW=0`, `W*W=I_2`, and let `C` be the staggered
spin antisymmetry satisfying `hC=-C conjugate(h)` and
`C conjugate(W)=W i sigma_2`. Extend the physical CAR carrier by the actual
native fifteen-state matter representation and three declared family copies.
For its symmetric Higgs kernel `K`, the local pairing
`Delta=C tensor conjugate(K)` has the exact localized Nambu inclusion

```math
E=\operatorname{diag}(W\otimes I_{45},\overline W\otimes I_{45}),\qquad
BE=E\begin{pmatrix}0&i\sigma_2\otimes\overline K\\
 (i\sigma_2\otimes\overline K)^\dagger&0\end{pmatrix}.
```

All twelve actual native current operators intertwine through this same
inclusion. The charged-lepton block uses `Ye=sqrt(2)M(P)` with neutral Higgs
`H=(0,v(P)/sqrt(2))`; its localized singular energies therefore have the
retained three mass coefficients. A complete displayed kernel also uses
the declared Stage5 quark powers and three zero neutrino coefficients.
Its rank is 42, with three zero directions. The 180 localized Nambu rows
contain 84 positive, 84 negative and 12 zero eigenvalues. They count spin and
charge-conjugate coordinates, not 180 independent particle species.

These are spatially uniform internal current transformations. Dynamical gauge
links and local Gauss-law constraints are not supplied by this inclusion.

The full cube retains 2,430 physical annihilator modes and 4,860 Nambu rows.
In the chosen `v(P)` unit, the source radial complement has gap at least 62
and `||K||<1`, so the localized single-quasiparticle bands are separated from
the full complement. Bulk energies are `lambda(h) +/- s_j(K)`. This additive
identity supplies no relativistic dispersion relation.

Each complete memory record can drive the corresponding local pairing pulse
with fixed kinetic and quark terms. Every pulse preserves the same localized
inclusion. A proved `O(1/N)` weak-pulse bound, uniform over finite untouched
references at fixed volume, gives the stated quadratic Hamiltonian as the
limit of this correlated record process. Current probes use fresh auxiliary
legs. This is an explicit coherent service, not a claim that original
site-resolved field repair preserves the kernel. The finite unread channel
need not be unitary.

The CAR state inclusion has exact partial-trace recovery with a tracial
complementary state; that state is not a selected physical vacuum. The modes
are bound to a prescribed winding texture, and their occupation does not
create its winding. Three family copies are supplied, not counted by the
two-dimensional radial index. Autonomous defect formation, interacting pole
masses, physical spin and clock identification, and a common relativistic
continuum with Einstein dynamics are not established by this finite join.

The owning calculations can be reproduced from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.mcpr_native_memory_bridge --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.defect_current_bridge --check
```

## Finite interacting poles and texture constraints

[Complete winding certificates](../../model/certified_winding_core/PROOF.md)
retain all 28 source matrices and their multiplicities: 16 have exactly two
zero modes, six are invertible, and six are unresolved. Their complete-Fock
contact estimates certify 34 of 84 supplied cases with an isolated ground
pair and equal-gap addition/removal states. These interacting poles belong
to the two-spin source; they are not silently transferred into the extended
quadratic model above.

[The radial current classification](../../model/radial_current_map/PROOF.md)
gives exactly eta-pairing SU(2), modulo constants, for conserved sums of
onsite even CAR operators. This explains the need for an explicit internal
carrier extension. [The ordered-spin texture bound](../../model/lattice_texture_selection/PROOF.md)
forces finite-energy rotation textures to have eventually admissible,
degree-zero shells. Its crossing-energy floor does not prove a particle
lifetime or identify that spin texture with the prescribed fermion mass.

This bound also constrains the proposed identification: a radial hedgehog
has degree one on every surrounding shell. If it is identified with this
unscreened ordered-spin field, its total energy diverges with radius. The
finite bound-state energy then cannot represent the rest energy of an
isolated finite-energy defect. Gauge screening, a different localized
topology, or a different dynamical texture law requires an explicit source
construction and a recomputed spectrum.


## Declared field candidates and measured-mass comparison

The [particle-response construction](../../model/particle_response/README.md)
resolves two charge-one excitations of the same interacting dimer about its
supplied half-filled ground of `K=H−N/2`. The lower `N=1` hole and `N=3`
particle are odd-parity internal SU(2) doublets, related by an explicit
particle–hole symmetry. Both have model gap
`Delta=sqrt(17)/8−1/16 ≈ 0.4528882032` and local single-mode residue
`0.3106339063`. Complete source Ramsey records reproduce their addition and
removal correlations; all preparations, pointer outcomes and failures count.

The [comparison packet](../../model/particle_comparison/README.md) declares
conditional electron-like and positron-like assignments. One CODATA electron
mass calibration sets `kappa=(m_e c²)/Delta`. Both mass energies are then
`0.51099895069 MeV`, with fully correlated calibration uncertainty. Their
predicted normalized difference is zero and satisfies the explicitly
conditional PDG positronium bound. **No absolute particle mass is independently
predicted.** CAR number charge and an internal doublet have not been identified
with electromagnetic charge and Lorentz spin. A physical vacuum, continuum
dispersion, autonomous native protection and a same-excitation topological-defect map
are unproved. The two-site frequencies alone cannot supply those maps.
The [executable dictionary](../../model/particle_comparison/dictionary.py)
supplies their conditional SI operators and a common time/length conversion.
Its complete current/dipole response passes the finite external-electric-field
checks. The magnetic moment and QED photon interaction are undetermined.
A separate periodic four-site check fails the fixed-mass/bare-speed lattice
Dirac benchmark, including a diagnostic using its own rest gap. Hopping,
fermion wrapping signs and momenta agree with independent reconstruction.
The boundary/volume change and possible velocity renormalization prevent a
continuum no-go; two momentum-squared values cannot validate a fitted law.

An exact audit of all signed charge-one transition frequencies excludes every
assignment to the electron, muon and tau mass-ratio windows, even after an
arbitrary common chemical-potential shift and positive energy scale. This
rejects that finite mass attachment, not the full OPH proposal. The fixed
alpha/ACT data replay is also retrospective and supplies no universe-membership
test. The complete physical comparison gate remains open.

## Same-candidate transport and number-conserving reads

The [charged-particle transport construction](../../model/particle_transport/README.md)
uses these same dimer ground and charged doublets. A complete graded block
swap transports an unknown internal state, including its reference, between
two identical separated dimers and restores the vacated ground exterior.
Source words prepare the states from counted fixed-occupation inputs and
execute both transport and an exchange interferometer with complete pointer
and export-failure records. Ordinary swaps and which-path records supply
distinct controls. This is coherent supplied transport; the generator is
the sum of the two dimer generators, not the uniform four-site lattice.

An explicit isometry places two identical odd candidate columns among twelve
even ground blocks. Graded block permutations reproduce every movement,
rotation and complete C/E record branch of the existing negative H66 lift
on this 66-dimensional image. Multiplication gives the corresponding
finite-history intertwiner. This is a formal 48-mode embedding, not a
48-mode source execution or a native protecting invariant.

Neutral inter-dimer transfer has a complete three-band particle-hole response
with analytically fixed weights. A supplied occupied or empty reservoir and
number-conserving reflection probes recover the same addition/removal Green
functions through eight complete source readouts. The subtraction combines
record probabilities without postselection. Preparation, source controls,
reservoir charge, neutral pointers and Hamiltonian model time remain declared
resources. A local number read disturbs the charged lower band.
The failed dispersion benchmark and physical-particle identification gates
therefore remain unchanged.

## Native memory and energy conditions

The [mobile-memory construction](../../model/particle_screening/mobile.py)
carries the same unknown qubit between all 60 admissible core locations.
It keeps naturally moved cores and intervenes only at the six impossible
locations. The optimal long-run core-reset rate is `1/150` per native
attempt, compared with `1/15` for the fixed edge-core service. Coherent
recovery follows **every** attempt. The construction below establishes
controlled transport and its exact maintenance requirement; it supplies no
physical species assignment.

The original [exchange-sector packet](../../model/exchange_sector/README.md) gives a
conserved *source/core* overlap `q` on binary H66, an exact two-dimensional
`q=1` code and a source-executed one-seam recovery with explicit core reset
and archive cost. Its all-seam correctability is a finite code theorem; the
executed recovery uses the declared `(0,5)` seam. Every sharp `q=1` state
differs by at least `4/9` from the transported reference on the full ten-port
exterior. An uninterrupted positive odd closed exchange flag has a rank-one
source input domain, so it cannot carry an unknown coherent qubit through a
final decoder. Intermediate active recovery changes that grammar. These
results leave physical localization open; coherent statistics requires the
separate complete operations specified below.

The [screening successor](../../model/particle_screening/README.md) reuses the
finite all-seam code theorem to improve the candidate. Exact correction
forces its exterior support onto ports nonadjacent to both core ports.
There are four such ports for each of the 30 edge cores, two for each of
the 30 distance-two cores, and none for the six opposite cores. A qubit
distributed over four safe exterior ports has worst one-port trace distance
`61/450`, compared with `82/225` for the original code. This is optimal
among exact all-seam `q=1` qubits with an edge core. Without the code
constraint, the sharp-sector optimum is `8/225`; exact hypergeometric
averaging gives the corresponding bounds for every exterior cap size.

The screened code is algebraically correctable under all complete C/E branches,
and its [successor recovery service](../../model/screened_recovery/README.md)
source-executes all 208 positive C/E recovery branches. Exact composition
preserves arbitrary logical/reference states when immediate correction and
explicit core-label reset follow each native attempt. Its native core changes
with probability `1/15` per uniform-seam attempt and leaves the edge set
with probability `1/30`. Separately, a sharp `q=2` pair initially on an edge
has first-exit survival `(14/15)^t` and mean first exit 15 attempts under
uniform-seam, fair-coin scheduling. This
tests an invariant edge bond; it does not determine a particle lifetime in
an observer's calibrated clock.

For both codes, complete uniform-seam C/E records give movement probability
`1/6` per *attempt*, versus `5/33` for the uniform-source reference: an
excess `1/66`. Relative to the separately supplied fixed-core three-charge
reference, the excess is `17/1125`. Neither statistic includes coherent
recovery work or a calibrated duration. The
[co-moving balance result](../../model/clock_defect_energy/REFERENCE_BALANCE.md)
has a stationary joint source/core measure, but the native gradient form of
its conserved charge is exactly zero. Its sector probabilities do not select
a rest-energy splitting.

The [same-code energy packet](../../model/active_code_energy/README.md)
imports the screened code's canonical definition and makes that freedom explicit.
Each connected charge sector has a uniform stationary measure. Multiplying
the `q=1` sector weight by two preserves the native law and its conditional
`q=1` records, while shifting a separately declared Gibbs modular contrast
by `-log(2)`. Stationarity therefore cannot select the physical energy gap.
The stationary reference averages over core labels, whereas the codes start
at fixed core; their current excess is a diagnostic comparison.

The original code's executed recovery on seam `(0,5)` costs mean repair plus
recovery pulses `439409` for C and `863483/2` for E per attempt, plus its
preparation and stated archive services. These counts belong to that fixed
seam and original code. Combining them with uniform-seam activity would
require a common protocol. They measure compiler work, not rest energy.

The [clock-count inference](../../model/clock_defect_energy/CLOCK_COUNT_INFERENCE.md)
provides a useful clock join: after hiding the full native history, the joint
modular detector and field readout estimate their shared incident count more
accurately than the detector alone. A control that keeps both separate
marginals removes this gain. The parameter is a dimensionless incident count;
physical duration and its common action normalization remain open.

The [positive current form and product-mode comparison](../../model/clock_defect_energy/MODE_ENERGY.md)
remain finite checks. The product mode's rest gap and independent curvature
approach the same mass-shell value under source-word subdivision. That mode
belongs to a separately declared field source, with no identified map to the
active `q=1` code. Using `E=mc²` for a particle here requires one localized
excitation to supply a reference-subtracted rest action, moving kinetic
response with a common signal speed and calibrated clock, and conserved
stress. Defining `m=E/c²` before those independent checks would rename an
assumed energy. No particle species or mass in SI units is identified.

## Exchange topology and observable statistics

The [exchange-statistics construction](../../model/exchange_sector/STATISTICS.md)
completes the path classification on the declared native two-particle
configuration complex: its fundamental group is `Z/2`, not merely its first
homology. Its two irreducible flat unitary characters are `+1` and `-1`.
This classifies transport phases on that fixed complex; it does not enumerate
all possible defects on arbitrary carrier refinements.

The [finite presentation certificate](../../model/exchange_sector/PRESENTATION.md)
establishes the same classification directly from all 235 generators and
495 native relators. Its 234 verified noncommutative Tietze eliminations
leave `<g|g^2>`; the generator and native exchange both have odd ordered-cover
parity. The trace uses 167 identity substitutions and 67 single-letter
substitutions. This supplies a compact exact proof without adding higher
cells or relying on continuous surface classification.

The negative lift's moving blocks agree with number-conserving fermionic
bilinears on the two-occupation sector. Bare C/E uses the positive lift.
They have the same classical trajectories, but need not have the same complete quantum
outputs. A record-selected, source-compiled diagonal operation implements
the negative lift, and a coherent read distinguishes it from bare C/E.
The coherent exchange-loop experiment separately makes the loop sign
observable. Committing which-path information destroys that contrast.

The [signed conservation theorem](../../model/native_fixed_algebra/SIGNED.md)
shows that both lifts have the same faithful diagonal reference, full source
matrix algebra and scalar universally conserved source observables under
positive exploration. An exact operator-space sign involution intertwines
every complete branch, so their channel spectra and GNS contraction bounds
coincide. Its `-1/3` positivity counterexample prevents treating it as a
quantum channel or physical gauge. The negative lift includes its supplied
phase service, and the common contraction bound is per complete lift
application. It is neither a common physical duration nor a mass gap.

Thus the exchange sign has an explicit operation and measurement interface,
rather than being a name assigned to a retained classical path. The coherent
controls and sign choice remain supplied parts of the observer program.
Neither choice produces physical spin, a protected reference-compatible
localized excitation, a rest gap or a Standard Model species by itself.
Creation and annihilation are outside this fixed total-two carrier. CP-03,
CP-04 and CP-05 require those additional same-process identifications before
an observed particle can be assigned.

## Local field and reference attachments

The separate supplied [complete-Fock Weyl/contact field](../../model/full_fock_field/LOCAL_SIMULATION.md)
has a volume-linear split coefficient and complete local-history error
bounds using nested buffers. On the `2x2x2` field at model duration `1/32`,
one step and 10,560 pulses meet operator error `10^-4`; the norm-sum bound
requires ten steps. This retains every number sector and untouched reference.
The local schedules are independent of ambient volume, but their initial
restricted states and observation services remain supplied. These bounds
do not identify a field excitation with the actively recovered native qubit.

An [executed Gibbs preparation](../../model/field_thermal/README.md) supplies
the same interacting action's reference on one or two sites. All number
sectors, purifier, read outcomes and aborts are retained. The faithful field
marginal supports signed CAR restrictions and the noncommuting finite
entropy identity; the pure joint purification is not a faithful full-algebra
reference. This bounded preparation does not provide the separate many-site
field reference, select a vacuum or calibrate physical energy and time.

The [rational support retract](../../model/support_refinement/README.md)
preserves homology equivariantly across the 12-, 42-, 162- and 642-vertex
support levels. Its midpoint averages and chain homotopy supply a compatible
finite refinement interface. They are not vertex-valued port attachments
or quantum channels, so they do not transport a native particle into
the spatial field or generate responsive geometry.

## Completed recovery, transport and unit attachment

The screened code's own uniform-seam service costs `37499681/60` C or
`36870437/60` E source pulses per attempt on average, plus initial encoding
and declared apparatus services. Its [transport calculation](../../model/screened_transport/README.md)
has rank-eleven transient current covariance with trace `1/3`, but corrected
logical transfer is identity and the core returns to `(0,1)`. The earlier
mode profile gives current/variance ratios `13/60` and `13/300` on the two
logical basis states. It therefore does not select one wave rate for this
code. A globally invisible phase is not a measured zero rest energy.

The [quantified sector theorem](../../model/conserved_sector_action/README.md)
proves the stationary-reference freedom for every positive sector factor,
with exact native support checks. The [four shared-event experiment](../../model/clock_defect_energy/shared_event_clock/README.md)
executes ordered conditional ensemble clock and field readouts through
three FIFO exports. Raw probabilities fail a two-anchor affine comparison
at a reserved event; model-conditioned count inversion is a separate step.

The [SI utility](../../model/clock_defect_energy/SI_UNITS.md) takes explicit
seconds per model clock unit and a model signal speed. Under a declared
energy-frequency identification it converts a supplied coherent rest gap
to joules and kilograms. This candidate supplies no such gap or independent
massive dispersion, so no particle mass or GeV prediction is reported.

## Mobile quantum memory and optimal core resets

The [contract](../../model/particle_screening/mobile_contract.json) fixes
binary H66, an independent uniform choice among thirty seams, the complete
fair-coin C or E instrument, and an unknown source-only logical qubit with
any untouched reference. At every ready boundary the code must be sharp
`q=1` and exactly correctable for **all** seam instruments. The core is
classical metadata; logical information is not stored in a persistent
auxiliary. Codes and core-reset decisions may depend on the complete past.

This covers more than the equal-amplitude examples. For core `{a,b}`, write
an arbitrary encoding as `V = sum_j |j>_exterior tensor B_j`, where the two
rows of `B_j` represent occupation at `a` and `b`. The existing all-seam
flagged-correction theorem forces nonedge support and
`B_j* B_j = w_j I_2`. Only safe exterior ports survive. A nonzero square
`B_j` is `sqrt(w_j) U_j`, with `U_j` unitary, and `sum_j w_j=1`.
Scalar compression of the two core occupations then forces
`V* n_a V = V* n_b V = I_2/2`. Conversely, such matrices give an all-seam
code when `sum_j w_j U_j* Z U_j=0`: the complete binary C/E effects are
linear combinations of `I,n_u,n_v,n_u n_v`, and all native-edge pair
occupations vanish on this support. Thus this is a classification of the
admitted exact codes, including complex encodings, rather than an assumption
that the displayed examples exhaust them.

For every such code, an exterior neighbor of a core endpoint is unoccupied.
Selecting that seam, finding its core endpoint occupied, and obtaining the
moving coin has probability `(1/30)(1/2)(1/2)=1/120`. It follows that
the natural core walk is independent of the unknown qubit, its reference,
and every adaptive choice of admissible code. Its complete type census is:

| Core | Count | Safe exterior ports | Mean attempts to first impossible core |
| --- | ---: | ---: | ---: |
| Native edge | 30 | 4 | 150 |
| Distance two | 30 | 2 | 120 |
| Antipodal | 6 | 0 | No exact all-seam sharp-`q=1` qubit |

Before the first antipodal arrival, the killed transition matrix in
edge/distance-two order is

```text
Q = [[29/30, 1/30], [1/30, 19/20]].
```

The exact survival after `n` attempts is `e_initial Q^n 1`, with eigenvalues
`(115 +/- sqrt(17))/120`. Both are below one. A code that follows only the
natural core cannot continue indefinitely under this contract, even if its
encoding and recovery are redesigned after every record.

The continuing policy retains every legal arrival. At an antipodal arrival
it sets the core to the attempted seam, which is an edge, and recovers into
the code there. This decision uses only the committed native record and
previous core. It changes the declared relational boundary and the recovery
target; it is an active intervention. The corrected class matrix is
`[[29/30,1/30],[1/20,19/20]]`, with stationary weights `(3/5,2/5)`.
Starting from an edge, the exact expected reset count is
`(n-12+12(11/12)^n)/150`.

The rate is optimal over the whole admitted class. Let
`h(edge)=150`, `h(distance_two)=120`, `h(antipodal)=0`.
Every native attempt decreases `h` by one in conditional expectation;
each counted core intervention can increase it by at most 150. Telescoping
gives `E[N_n] >= max(0,(n-h(initial))/150)` for arbitrary history-dependent,
randomized and preemptive resets. The stated policy attains rate `1/150`.
It also minimizes the expected count for every fixed finite horizon. The
edge and distance-two value functions differ by
`(1-(11/12)^n)/5`, strictly less than the unit cost of an earlier reset.
The exact Bellman recurrence therefore keeps a legal core and, only at a
forbidden arrival, resets to an edge. The certificate checks the recurrence
symbolically and the tests solve the independent finite-horizon optimization.
This minimizes core-label interventions, not coherent pulse count or energy.
Changing the seam schedule, permitting only chosen-seam correctability,
storing the qubit in an auxiliary, or enlarging the carrier changes the
optimization problem.

The [compiled service](../../model/particle_screening/mobile_service.py)
uses actual proper port rotations of the two existing code prototypes.
Their stabilizers fix both logical columns, so all 60 frames are consistent
without another frame register. For each complete positive flag it implements

```text
R_h(c) K_h V_c = sqrt(p_h(c)) V_target(c,h).
```

The inverse branch encoder followed by the target encoder is a full
512-address unitary extension. Complete zero/abort outputs remain present:
valid zero-code-weight slots use identity recovery, and aborts terminate.
On encoded ideal inputs the zero slots have zero probability; their
implemented amplitudes are checked. The identity composes along every
finite history, so the complete tape contains no unknown logical information
and the same logical/reference state survives at the record-selected core.
This supplies a continuing **classically located** quantum memory, not a
coherent superposition of positions. Every closed corrected core loop has
identity logical transport in these frames; no fermionic phase is generated.

Source checks independently bind the primitive address chart, check all
emitted permutation layouts, and count every native/recovery pulse even when
the numerical macro is cached. Fresh full-register macro and dilation errors
give product composition diagnostics for all cores and raw flags through
100 attempts. Complete stopped histories actually execute safe motion,
antipodal reset and a closed exchange loop, retaining all complementary
quantum outputs without restart or renormalization. Fixed numerical gates
are replayed, not interpreted as hardware accuracy or infinite-time bounds.
The [receipt](../../model/particle_screening/mobile_results.json) gives
the per-core, worst-case and stationary resource budgets. Operational storage
is the seven-bit current core plus the existing complete native tape
(eight C or nine E bits per attempt); deterministic reset decisions and
targets can be reconstructed from it. Mixed-law selection needs its tag.
Source blanks, native apparatus services and growing export capacity remain
explicit resources.

The full 60-core stationary law assigns probability `1/50` to each edge and
`1/75` to each distance-two core. The reset controller introduces sixty
directed links with no reverse. Therefore native reversible-action theorems
do not automatically apply to this corrected process. No physical entropy
production or rest energy is obtained by assigning a price to a reset.

The bare-law and active-memory catalogue is:

| Object | Established invariant or dynamics | Particle assignment |
| --- | --- | --- |
| Double-load collision | Transient under complete native repair | None |
| Collision-free exchange history | Exact 234-step Tietze certificate gives fundamental group `Z/2`; an uninterrupted odd closed flag resolves its input | Bare C/E has the positive lift |
| Source-operated exchange statistics | Both flat characters have complete instruments, coherent readouts and equal channel spectra; negative moving blocks are CAR bilinears | Similarity is nonpositive; sign and controls supplied, no physical spin, protection or species assignment |
| Relational source/core overlap `q` | Conserved under simultaneous current transport | A relational sector, not a source-only charge |
| Sharp `q=2` edge pair | Mean native-edge first exit 15 attempts | No protected edge bond or measured decay law |
| Mobile sharp `q=1` qubit | Exact active correction, 60 admissible cores, optimal resets | A transported quantum memory; rotations act trivially on its logical qubit |
| Complete binary occupation carrier | Fixed algebra `C^11 direct_sum M_2(empty,full)` for either C/E law and exchange character | Global number charges and a distance-one global qubit; no local protected source charge |

The [complete occupation classification](../../model/native_fixed_algebra/OCCUPATION.md)
exhausts passive source-only conservation on all 4096 binary configurations.
Actual native effects separate every configuration except empty/full, and
signed native transport generates every matrix unit within each intermediate
number sector. The only nontrivial quantum multiplicity is therefore the
global empty/full span. Every proper-port-local fixed observable is scalar;
a one-port phase acts as a logical error, while an unread occupation
measurement dephases the memory. For a
fixed positive all-seam selector, the complete stationary-state family and
its limiting channel are also determined, without assigning a physical
clock or convergence rate.

The adopted V2 program fixes total number and forbids cross-number coherence,
so the extremal qubit is a separate component probe, not an admitted V2
observer. A complete twelve-occupation-bit source execution preserves the
extremal code under native repair, while a complementary six-occupation code
loses coherence to the retained native records. All preparation, readout,
complementary outputs and resource budgets are included. The classification
does not exclude moving particles carrying a global charge, active recovery,
relational source/controller invariants, or refinement protection. It
separates the exchange-path invariant from passive local state protection;
neither alone establishes physical particle identity.

This classification concerns the bare twelve-port, one-bit-per-port all-seam
law. The [charged internal factor](../../model/charged_defect/README.md),
[common spin factor](../../model/common_spin_protection/README.md) and
[protected cube sector](../../model/cube_eta_protection/README.md) use distinct
field-block instruments, with internal states or a retained outside projector.
Their positive protection results do not contradict the bare-law catalogue.
Even the spin field's matching 4096-state dimension does not identify its
factorization or instrument with that law; no general exclusion of passive
field protection follows.

None supplies the missing joint localization, coherent statistics, matched
rest/kinetic action and detector response required to name a physical species.
The mobile-memory objects above have no established Standard Model current
map. The separate extended radial construction at the start of this document
does have that operator inclusion, under its stated coherent-service design.
Neither result establishes all physical requirements of CP-03, CP-04 and CP-05.

```sh
python -m model.particle_screening.mobile_packet --check
python -m pytest -q -W error model/particle_screening
```
