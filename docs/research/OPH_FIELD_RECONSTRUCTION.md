# Field response and particle identification from observer records

## Chiral fields on an actual native source

The [sitewise chiral source](../../model/coherent_defect_wave/chiral_source/PROOF.md)
encodes the full `F(C2 tensor R15)` at every cell, with thirty fermion modes
and all `2^30` occupations. Each occupied object contains the actual paired
flux, two-component routing memory and a native dark cloud. A cell-local
antisymmetrization intertwines every actual subgroup current and preserves
the support of even gates in the declared graded source algebra. It avoids
a global particle-slot sorting operation. Ordinary qubit compilation with
bounded-range gates is a separate implementation question.

The internal module is `R15=Lambda^2(C5)+Lambda^4(C5)`. The twelve currents
of the trace-free subgroup act with weights

`(3bar,1,-2/3)+(3,2,1/6)+(1,1,1)+(1,2,-1/2)+(3bar,1,1/3)`.

This is one all-left-Weyl family, with right singlets represented by their
left-handed conjugates. The extra central native degree current is separate.
The native idle law fixes the even cloud, so `K=0` supplies no binding or
mass coefficients. Grading, cloud binding and the kinetic law are explicit
architecture choices.

### Propagation, reference and the regulator boundary

The [finite spectral regulator](../../model/coherent_defect_wave/chiral_regulator/PROOF.md)
uses `h_L(p)=-sigma.p tensor I15` on the complete antiperiodic cubic momentum
menu. All cutoff modes are retained; no mirror state is removed by
postselection. Filling the negative band gives a unique finite free ground
and positive particle/hole energies. The finite gap is `sqrt(3)*pi/L`,
which vanishes with physical volume and is not a rest mass. The particles
have helicity `-1/2`; their holes have helicity `+1/2` and conjugate currents.
They are antiparticles rather than independent same-charge mirror fields.

The finite position derivative is long-range. Its Fourier symbol equals
the free continuum symbol exactly on retained windows, and the bound
`||(I-Pi_N)f|| <= [1+(pi(N+1)/L)^2]^(-s/2) ||f||_Hs`
controls cutoff tails and all finite free-field words. Causal locality and
Lorentz covariance hold in that free effective continuum. There is no
strict finite-source light cone or proved interacting gauge continuum.
The free filled state is not substituted for a fluctuating-gauge Gauss
reference.

### Gauge, scalar and massive field interfaces

The source supplies finite sixty-state quantum links with exact endpoint
Ward identities on the same full matter occupation algebra. Their charged
raising coefficients are nilpotent and cannot act as identity links in a
vacuum. Supplied-background propagation and quantum-link propagation must
therefore retain distinct state and operator domains.

The [quantum Higgs construction](../../model/coherent_defect_wave/chiral_higgs/PROOF.md)
uses the actual invariant Weyl pair tensors, an explicit five-state charged
scalar memory and a Hermitian gauge-invariant interaction. Its complete
lepton/scalar calculation supplies a finite quantum response. A nonzero
neutral classical Higgs background is a separate benchmark, not the
expectation of a gauge-invariant finite reference.

For that benchmark, the [canonical particle dictionary](../../model/coherent_defect_wave/chiral_particle_map/PROOF.md)
is an exact change of the thirty original fermion variables. It combines
weak-doublet annihilators with conjugate-singlet creators to give seven
Dirac fields counting color, plus one neutral left-Weyl field. No independent
mirror family is added. Actual electric charge acts equally on both halves
of each massive field, color is vectorlike, and the weak currents act only
on the left components. The source-to-field map therefore includes currents
and canonical operators, rather than a dimension count alone.

| Conditional field role | Number counting color | Charge | Source rest energy |
| --- | ---: | ---: | --- |
| Up-type quark | 3 Dirac fields | `2/3` | `abs(v*y_u)` |
| Down-type quark | 3 Dirac fields | `-1/3` | `abs(v*y_d)` |
| Charged lepton | 1 Dirac field | `-1` | `abs(v*y_e)` |
| Neutral lepton | 1 left-Weyl field | `0` | `0` |

The paired branches have `E^2=|p|^2+abs(v*y)^2`. On the antiperiodic finite
menu no zero momentum is present; rest energy means the parameter of this
controlled field relation as physical volume grows. The scalar background
and three Yukawa values are inputs. No generation or measured species is
selected, and the massless one-family neutral field cannot reproduce neutrino
mass splittings. Unit conversion gives `m c^2=u_E*abs(v*y) GeV`; `u_E` is
an independent energy standard.

### One quantum Hamiltonian, physical subspace and clock

The [common quantum gauge/Higgs construction](../../model/coherent_defect_wave/chiral_joint/PROOF.md)
places the spectral coefficients, actual quantum-link path products, local
Yukawa terms and invariant scalar/link energies in one Hermitian H. Every
local total Gauss current commutes with H. For `beta>0`, the physical Gibbs reference
`rho_beta=P_G exp(-beta H) P_G / Tr(P_G exp(-beta H))` is faithful on its
nonempty physical subspace and stationary for the same `exp(-it H)`.
A complete success/failure filter prepares it from the maximally mixed
physical input, using a specified controller that preserves total Gauss
charge and total fermion parity, with its graded ancilla retained. The
explicit side-two norm bound is `||H||<=1884` at unit parameters. This is
a finite spectral-controller construction. The full
cubic-family eigensystem and efficient preparation are not supplied; the
complete two-cell physical spectrum below has a separate certified calculation.

A separate control on the same operator grammar retains all 122880 states
of two lepton cells, a quantum link and a quantum scalar. Starting from
empty matter, its first action creates a gauge-compatible pair/scalar;
its second propagates a created lepton to the other cell. The exact
arrival component has `||R H^2 Omega||^2=3`, with a strictly positive
complete finite-time arrival bound and the no-arrival alternative retained.
This input is explicitly distinct from the Gibbs reference. The common
finite construction supplies no interacting Lorentz vacuum, physical Higgs
phase, mass renormalization or GeV prediction. No primitive magnetic plaquette
term is inserted in the declared Hamiltonian; the matter-generated loop
calculation below uses its actual hopping terms.

### A complete two-cell physical ground and boundary sectors

The [uniform two-cell calculation](../../model/coherent_defect_wave/chiral_joint/two_cell/SPECTRUM_PROOF.md)
retains both quantum Yukawa terms, both scalar memories and the actual
quantum link on one 614400-state active source. Exact Gauss projection
leaves 1053 neutral physical states, 718 with right boundary charge `-1`
and 783 with right boundary charge `+1`. Every conserved occupation block
and every physical energy is enclosed; the calculation is not restricted
to states reached from an empty input. The neutral physical ground is
unique, has `(L,2Sz)=(2,0)` and energy approximately `-4.554826635` in the
specified clock. Its certified excitation gap exceeds `0.394804412`.
The exported energy intervals have radius `2^-28` source units. A complete
energy instrument prepares this ground from the physical maximally mixed
input with original probability `1/1053`, retaining the other outcomes and
a compensating preparation reference carrying the opposite L and twice-spin projection.

The right-negative-boundary ground is an odd-L doublet near `-3.502782787`.
The right-positive-boundary absolute ground near `-3.712362234` has L=4,
so a single removal from the L=2 neutral reference cannot reach it. The
complete right spin-up e^c instrument instead certifies the first visible
addition transfer in `[1.0520438402,1.0520438552]` and removal transfer in
`[0.9097986244,0.9097986394]`, in source energy units. Their positive
joint record weights are enclosed by `[0.2183892428,0.2183892504]` and
`[0.2308929339,0.2308929414]`, respectively. These probabilities retain the
reservoir setting factor one half and condition only on successful ground
preparation. The unconditional protocol multiplies them by `1/1053` and
retains preparation failure `1052/1053`. Removal reaches L=3; the L=4
absolute charged minimum is a zero response of this probe. These finite
energy transfers are visible detector lines, without a rest-mass or GeV
identification.

The [exact symmetries](../../model/coherent_defect_wave/chiral_joint/two_cell/SYMMETRIES.md)
prove the odd-sector Kramers degeneracy and map right addition to left
removal with the charge-compensating boundary transformed. They do not
require equality of addition and removal at the same pinned boundary.
The one open quantum link admits no literal source-field translation.
A separately declared two-oriented-link periodic source has an exact
antiperiodic translation and a complete momentum-phase PVM. A physical
local pair/scalar preparation gives all four phase weights `(1/2,0,1/2,0)`.
Its state is not the open-edge ground, and that ground's energies are not
transferred to the periodic Hamiltonian. Three-dimensional dispersion,
physical spin, interacting rest masses and an independent GeV scale remain
separate physical requirements.

### Controlled motion and the finite particle dependency path

The [propagation calculation](../../model/coherent_defect_wave/chiral_joint/two_cell/PROPAGATION_PROOF.md) prepares an actual
right-cell matter pair correlated with its quantum scalar, using the local
Yukawa pulse `exp(-i*pi*H_Y/4)`. All Gauss constraints hold, every preparation
branch is retained, and the left cell starts empty. At source time `1/20`,
a complete receiver-parity detector has arrival probability in
`[63853775/34359738368,128634893/68719476736]` on the open source. The separate
periodic source gives `[31793921/68719476736,32424155/68719476736]` and equal
opposite-direction laws. Exact degree-eight evolution and full-H remainder
bounds certify both intervals. Removing only hopping makes arrival exactly
zero at all times. This measures motion on a nonstationary physical input;
it does not identify a ground-state particle or a scattering amplitude.

The [periodic comparison](../../model/coherent_defect_wave/chiral_joint/two_cell/REFINEMENT_PROOF.md) retains a charge `-1`
compensator with a two-position register and translates it with the source.
A colocated e^c occupation and compensator satisfy every Gauss constraint;
the complete joint translation-phase probabilities are `(0,1/2,0,1/2)`.
These are source-plus-reference phases. Their reduction to source-particle
momentum requires its own reference argument. Two antiperiodic cells have
only one odd-sector momentum magnitude, so they cannot determine both a
mass and a speed. Four cells give two magnitudes; six give three, permitting
an additional shell consistency test once energies are calculated.

Inserting a blank second quantum link does not transfer the open spectrum
to the periodic source. On an explicitly physical input the compressed
Hamiltonian residual vanishes, while the full residual has squared norm
`3/16`, entirely in the omitted link states. No constant clock rescaling
or energy offset removes that leakage. Each refined source therefore needs
its own physical reference and complete energy/momentum records, or an
explicit full-residual comparison that bounds every omitted component.

| Finite step | Result on the declared source | Requirement for a particle mass |
| --- | --- | --- |
| 1. Hamiltonian and constrained carrier | Complete active Gauss bases and actual quantum scalar/link H | Preserve that operator/state ancestry under refinement |
| 2. Reference and preparation | Independently certified open and periodic grounds, with complete preparation laws `1/1053` and `1/3424`; full-family local ground doublet and all-sector exclusion | Control the physical reference on a growing three-dimensional family |
| 3. Energy instrument | [2627 original-attempt periodic records](../../model/coherent_defect_wave/chiral_joint/periodic_spectrum/PROOF.md), including equal visible charged thresholds and quantum Higgs overlap | Establish the moving branch and its physical mass observable |
| 4. Motion and momentum | Actual Peter–Weyl neutral and all-five-family arrival; periodic energy and joint-reference phases on one ground; full three-axis operator family | Source-relative dispersion across independent momentum magnitudes and sizes |
| 5. Candidate discrimination | [Actual native currents, CAR, Kramers and conjugation tests](../../model/coherent_defect_wave/chiral_joint/two_cell/CANDIDATES.md) | Separate elementary, composite and boundary-bound responses; no empirical species is assigned |
| 6. Size, cutoff and field transfer | A finite dressed interaction has a full error bound; raw long quantum paths fail the blank-reference uniform bound | Control the chosen link architecture, physical reference and interaction under refinement |

The finite completion of steps 1–3 and the motion controls in step 4 do not
supply the physical conclusions in the final column. Scalar/flavor response,
interacting relativistic control and an independent energy reference retain
their separate roles in a GeV comparison.

### Periodic quantum ground, charged fields and controlled interaction

The [periodic quantum spectrum](../../model/coherent_defect_wave/chiral_joint/periodic_spectrum/PROOF.md)
uses two oriented quantum links and both scalar memories on 3686400 active
states. Its full neutral Gauss space has dimension 3424, with a unique
`(L,2Sz,T)=(2,0,+1)` ground at `-4.33817023295` and gap greater than
`0.25770063535`, in the declared source clock. Each compensated charged
space has dimension 2456. Complete rational eigenbasis bounds include every
conserved sector. These values belong to this periodic Hamiltonian and are
calculated independently of the open-edge values above.

Translation composed with the actual conjugation maps addition to removal
at the same boundary. Both probes have their first transfer in
`[2083265509/2147483648,2083265637/2147483648]`, approximately
`0.97009612853`, with original-swap probability in
`[29863037/134217728,14931519/67108864]`, approximately `0.22249697`.
Ground preparation succeeds with probability `1/3424`; failure `3423/3424`
is retained. The complete adaptive energy and joint reference-phase
instrument has 2626 post-preparation records and 2627 original-attempt
records. The charged phase weights `(0,1/2,0,1/2)` hold at each energy
window because the compensator's two recorded positions are orthogonal.
They are joint source/reference phases, with only one squared momentum;
they do not determine a source-only rest mass.

The [quantum Higgs operator identities](../../model/coherent_defect_wave/chiral_joint/higgs_composites/PROOF.md)
attach the native up, down and charged-lepton pair tensors to their actual
Higgs-screened operators. All twelve current identities and the creation
equations hold on the entire thirty-mode local occupation space. In the
periodic lepton control, the elementary charged operator and its quantum
Higgs composite excite the same first energy line. The latter's squared
matrix element lies in `[56785389/134217728,28392695/67108864]`, approximately
`0.42308412`. This is an operator residue, not an implemented preparation
probability. The exact same-H equation retains the hopping contribution
and its interference. No scalar expectation value is substituted.

The neutral contraction `phi0*c_L11-phi1*c_L10` is a gauge singlet and is
fermion odd. It is noncanonical and has a nonzero interacting equation of
motion. Local conserved star parities make its separated two-point
correlator zero in a unique-ground or Gibbs reference. A physical neutral
propagation claim therefore requires an appropriate link-dressed operator;
the absence of a native neutral singlet partner proves no massless neutrino.

The [three-axis periodic family](../../model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md)
retains all native quark/lepton occupations and supplies three commuting
translations with exact Gauss-compatible path transport. A physical
pair/scalar preparation has 64 complete joint phase outcomes, eight of
weight `1/8`. Its full interacting ground is distinct from the calculated
reduced lepton ground. The raw raising-link paths amplify a blank-reference
internal state by `d^(r/2)` on a path of length r in a d-dimensional block.
A [contractive-link alternative and finite effective-interaction proof](../../model/coherent_defect_wave/chiral_joint/field_limit/PROOF.md)
make the required architecture and comparison choices explicit. On an
actual fifteen-state Gauss block at link penalty 128, a six-state dressed
effective generator has complete-record error below `10^-12*abs(t)`.
A separate member at penalty 65536 has a rigorously positive connected
energy. Neither estimate transfers the unit-coupling ground or proves a
nonzero interacting relativistic limit.

The [source spectral-edge attachment](../../model/coherent_defect_wave/chiral_joint/spectral_edges/PROOF.md)
feeds the actual original preparation weight into the existing threshold
readout theorem. Particle identification uses the appropriate observable:
a stable gapped pole, an infrared charged threshold with inclusive records,
an unstable resonance, or a specified renormalized confined-quark mass.
An isolated atom is not required for every physical charged sector.

The [complete neutral-composite instrument](../../model/coherent_defect_wave/chiral_joint/neutral_composite/PROOF.md)
and [scalar/link-loop instruments](../../model/coherent_defect_wave/chiral_joint/periodic_boson_records/PROOF.md)
use this same periodic ground. Neutral creation and removal each retain
all 298 target states, both translation phases and their preparation and
apparatus failures. Their first transfer is in
`[2681604835/4294967296,2681605091/4294967296]`, about `0.62435981` source
units. The conserved local-star symmetry pairs these phases and forbids
the separated bare two-point response; this line does not describe a
propagating neutrino. A local scalar-number probe and two closed-link
probes retain all 502 target states, elastic records and failures. All
three have positive weight at the same inelastic line in
`[228185229/134217728,228185237/134217728]`, about `1.70011247` source units.
This is an exact common-reference response, without a Higgs, photon or
weak-vector identification.

### Spatial symmetries, neutral motion and the full-family cell

The [orientation-closed link construction](../../model/coherent_defect_wave/chiral_joint/oriented_links/PROOF.md)
uses 119 quantum states per link, with both flux orientations. Its actual
native-current matrices satisfy every endpoint Ward identity, admit unitary
orientation reversal and have contractive full internal transporters.
All 24 proper cubic rotations, including antiperiodic seam phases, preserve
the full periodic Hamiltonian. A two-pi rotation is fermion parity. This
finite spatial symmetry is distinct from continuous Lorentz covariance.
These minimal link operators omit nonvacuum Clebsch transitions of literal
Peter-Weyl multiplication; that difference is proved with actual currents.

The [literal Peter–Weyl source](../../model/coherent_defect_wave/chiral_joint/peter_weyl_links/PROOF.md)
keeps the same eleven representation sectors and derives all 37 nonzero
coordinate channels, including 27 transitions between excited sectors.
The complete Haar-normalized internal transporters are contractions. All
native endpoint currents, proper cubic rotations and flux-orientation
identities hold on the full link spaces. The additional coordinate
transitions break the minimal link's conserved local-star tag. This is an
explicit link-architecture choice with all three Yukawa coefficients and
both penalties equal to one; its spectrum is not inherited from the
six-state periodic lepton model.

On this actual full-family Hamiltonian, a physical neutral scalar–lepton
packet moves without an added neutral hopping term. Its first nonzero
neighboring opposite-spin moment is `-151i/32` at fifth order. Every native
matter block contributes to this coefficient. The complete preparation,
fermionic reference and both spin branches are retained. A full-H remainder
proves positive arrival at `t=151/4641746545968`, with both the exclusive
event and a local packet detector specified. The small certified time
establishes a nonzero response; it supplies no practical propagation rate,
stationary particle branch or relativistic mass.

The [compensated family-transport calculation](../../model/coherent_defect_wave/chiral_joint/family_transport/PROOF.md)
uses the same source and a single static dual-R15 reference for all five
native matter blocks. It calculates the full current weights and Casimirs,
including the nonorthogonal color-current metric. Each physical spin/gauge
singlet preparation develops two orthogonal one-link arrival strings, with
total first-derivative norm squared `1/8`. All thirty source components are
evolved twice with every Yukawa and link term retained. At the common time
`1/1508`, a detector on the neighboring cell and its incident links has
prepared-minus-unprepared probability at least `1/81866304`. The reference
and both detector outcomes are retained. This distinguishes transported
compensated gauge excitations from vacuum preparation noise; it does not
establish a free charged particle or stationary mass shell.

The [contractible-square calculation](../../model/coherent_defect_wave/chiral_joint/matter_induced_plaquette/PROOF.md)
uses this same Peter–Weyl Hamiltonian and a static dual-native reference.
Four actual hops return a prepared matter excitation to its origin while
leaving a gauge-invariant four-link trace. The fourth moments are `-1/256`
for the e^c representation and `-1/768` for the colored up-conjugate
representation. All chronological orders, spin signs and native currents
are included. Full-state bounds certify nonzero energy-dependent eliminated
loop kernels and complete finite-time loop records. No magnetic coupling is
inserted to obtain these records. A pure-gauge generator, Coulomb/confining
phase and physical photon or gluon identification require separate results.

An [exact rational finite clock](../../model/coherent_defect_wave/chiral_joint/finite_clock/PROOF.md)
preserves these detector distinctions with an explicit implementation bound.
For the unchanged H, the unitary diagonal Padé step has operator error at
most `abs(t*B)^7/100800`; n subdivisions reduce it by `n^6`. One step
preserves the fifth-order neutral arrival, and two steps preserve the
five-family local contrast with lower bound `t^2/48`. A separate exact
72-state execution retains both loop probabilities and their complement.
The service is finite algebraic functional calculus; an efficient local
circuit on the huge cubic carrier is not proved.

A comparison [declared neutral composite hopping term](../../model/coherent_defect_wave/chiral_joint/neutral_mobility/PROOF.md)
removes the local-star obstruction for one physical odd packet. It preserves
Gauss constraints, global parity and proper cubic rotations. All base quark,
lepton, scalar and quantum-link interactions remain active. An even,
gauge-neutral preparation retains its fermionic reference and both spin
branches. At source time `1/353`, a complete neighboring packet detector
has arrival probability in `[1/498436,9/498436]`; its complement is retained.
The full second-action norm `15948017/128` bounds all omitted Taylor terms.
With this extra hopping coefficient zero the specified arrival is exactly
zero. This is a verified design change and finite physical transport;
it supplies no stationary vacuum, relativistic mass or GeV identification.

The [full-family quantum cell](../../model/coherent_defect_wave/chiral_joint/family_cell/PROOF.md)
activates all three unit Yukawa channels on `5*2^30` local states. Its complete
vacuum-connected graph has 97878 states. An exact equitable compression and
Krylov closure produce a 36-dimensional Gauss-singlet cyclic space, with all
36 positive energy probabilities certified by rational polynomial root
isolation. The complementary code outcome is retained with zero weight.
The graph's unique ground is near `-11.14102996` source units and its ground
filter succeeds with probability near `0.00047723947`. These are a
preparation-sector reference and energy records, without a global minimum
claim across other sectors. Independent complete short-time number records
resolve up-pair, down-pair and charged-lepton-pair production. The labels
specify source representations, not measured quarks or lepton generations.

The [physical neutral local response](../../model/coherent_defect_wave/chiral_joint/neutral_local_spectrum/PROOF.md)
starts from the gauge-singlet two-component preparation `D_up^dagger Omega`.
Its complete connected graph has 113948 states; exact compression gives
49 positive spectral weights and a retained zero-weight complementary
outcome. The component ground is a physical local spin-one-half state near
`-11.23058529425`, strictly below the vacuum-connected minimum by about
`0.08955533`. This comparison rules out the empty-register cyclic ground
as the global local reference when these two sectors are both admitted.
It is not a negative rest mass. The raw conserved charge has conjugation
law `L -> 2-L`, so its centered value is `L-1`; this odd ground has centered
charge zero. Its conjugation partner has opposite spin in the same raw
charge sector. A two-cell singlet of these local grounds is an explicit
even physical reference for uncoupled cells. The interacting many-cell
ground requires its own calculation, including sectors with link-compensated
local gauge charge.

The [complete physical local reference](../../model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md)
compares every zero-boundary-charge local Gauss sector. An exact particle–hole
transformation makes the Yukawa operator a collective `sp(4)+sl(3)` action.
The full `2^24` quark occupation character excludes every physical quark
module except one 140-dimensional module below energy `-11.13`. Its complete
physical coupling to all lepton and scalar states has dimension 104,
including 56 dark states at energy one. Exact rational ordering proves the
global physical ground is the spin-one-half doublet near `-11.23058529425`.
The first excitation is the pair of raw-L=0 and L=2 spin singlets; the
complete physical gap lies in
`[25207584560067/281474976710656,25207584560069/281474976710656]`, about
`0.08955533047607` source units. This strengthens the cyclic-sector comparison
without identifying that gap with a particle mass.

The local density `P_ground/2` is stationary and rotation invariant. Its
native mean occupations are exactly `(3,6,1,2,3)` in the up-conjugate,
quark-doublet, charged-singlet, lepton-doublet and down-conjugate blocks,
respectively. Total mean occupation is fifteen. These are interacting
reference occupations, not fifteen detected particles. The reference is
specified by its ground projector; this result supplies no arbitrary-input
preparation protocol or many-cell vacuum.

### A controlled interacting reference and its collective modes

The [complete atomic lattice bound](../../model/coherent_defect_wave/chiral_joint/atomic_reference_control/PROOF.md)
includes charged local matter compensated by every allowed quantum-link
sector. On all `5*2^30` unconstrained local states, the actual Hamiltonian
is greater than `-11.71`. Exact positivity checks cover all 1686 decisive
weight blocks. At unit link penalty, a nonblank link can lower the energy
of at most its two endpoint cells. The resulting full physical atomic gap
is greater than `0.0411705884953`, uniformly in the finite volume. The atomic
ground space is the tensor product of local doublets, of rank `2^V`.

The same full hopping operator has norm at most 360 on the eight-cell
member and zero first-order compression onto that ground space. Its
complete charged-state resolvent gives a second-order spin interaction
with the [exact coefficient](../../model/coherent_defect_wave/chiral_joint/charged_ground_exchange/PROOF.md)
`J in [41579260299095/281474976710656,5197407537387/35184372088832]`,
about `0.147719206819003`. Complete charged response polynomials have
degrees nine and seventeen. Separate quark and lepton spin conservation
makes the quark contribution to this spin coefficient zero at every
virtual time; its scalar-energy contribution is retained. No charged
intermediate state or native family is discarded.
An [exact sublattice-parity argument](../../model/coherent_defect_wave/chiral_joint/atomic_even_remainder/PROOF.md)
bounds the difference between the canonical same-H generator and its
second-order expression by `25*(abs(lambda)*360)^4/g^3`, where g is the
atomic gap bound. The side-four source contains same-sublattice paths;
this improved error bound is specific to side two.

The [complete spin geometry](../../model/coherent_defect_wave/chiral_joint/atomic_spin_geometry/PROOF.md)
uses an exact local spin-coordinate rotation to obtain
`K2=-24*A+2*J*H_cube`. All 256 cube states are retained. Its unique singlet
ground has energy about `-4.82008937437` in H_cube units; its first
triplet is at exactly `-4`. Combining the full spectrum with the quartic
remainder proves a unique **full physical interacting ground** for the
declared family `H_lambda=H0+lambda*T`, throughout
`0<lambda<=2^-29`. At the reported coefficient `lambda=2^-29`, the
ground gap lies in the outward source-unit interval
`[6.95741157484e-19,9.85455023443e-19]`. This is a rigorous
error-and-coefficient enclosure, not a measured mass estimate.

The actual source translation is retained under the spin-coordinate
rotation. The interacting ground has zero momentum. Its first isolated
cluster has one mode at each of `(pi,0,0)`, `(0,pi,0)` and `(0,0,pi)`;
actual cubic symmetry makes their energies equal. These are collective
even-parity source excitations. The rotated atomic-spin labels do not
establish Lorentz spin, and the three momenta do not identify three
generations. Their lowest gap is not a zero-momentum rest mass.

Small hopping is a declared architecture choice; its uniqueness and
equivalence to unit hopping are not required. Physical identification
requires propagation, detector records and a controlled refinement on
that chosen family. The ground theorem does not transfer to unit hopping,
and the finite gap does not prove a continuum mass gap. A formal complete
ground filter on the maximally mixed full carrier has success `1/D` and
failure `1-1/D`, with `D=(5*2^30)^8*119^24`; an efficient native implementation
and its full histogram are not supplied.

The [complete interacting spin detector](../../model/coherent_defect_wave/chiral_joint/collective_spin_response/PROOF.md)
uses this same source and its actual stationary ground at `lambda=2^-34`.
The full spin current `O_a=sum_x (-1)^x_a S_x,a` has norm 60. Its physical
pump retains a four-state scalar/vector reference, which compensates each
source momentum. The exact finite dilation includes pump failure. Full
operator and state bounds give a first-line probability in
`[0.00204894419,0.00211090162]` and an outside-line probability in
`[0.00001386725,0.00001938783]`, conditioned on ground preparation.
The excitation energy is in `[8.20758426381,8.21034719056]e-22` source units.
Initial failure `1-1/D`, pump failure, and line/outside for each of three
axes form eight positive coarse records. Their exact probability sum is
one. The energy and momentum filters are declared finite spectral services;
the result proves a visible collective line, including a positive spectral
complement, without an empirical particle or rest-mass assignment.

The [centre/current catalogue](../../model/particle_catalogue_correspondence/PROOF.md)
recomputes the six classes of the 32-state exterior register and the actual
fifteen-state internal current representation. Class labels retain colour
triality and weak centre parity. They lose integer hypercharge, electric
charge, full non-Abelian representation and physical fermion parity. In
particular, the colour singlet and octet share the same centre, and an even
two-quark channel can have the same gauge representation as an odd `d^c`
field. These exact witnesses constrain transfers from a defect register to
the full interacting source without invalidating its internal candidates.

### Complete detector records on the physical local ground

The [actual scalar–lepton composite instrument](../../model/coherent_defect_wave/chiral_joint/ground_composite_response/PROOF.md)
acts on the full physical local ground density. Four spin/addition/removal
settings retain all 26 target states in each charge sector. Exact projector
moments give nine positive energy weights per setting, including the
lowest physical local gap. The joint first-line probability is about
`0.0177861244421` per setting. The complete law has 44 records, including
failures and zero-weight complementary outcomes. Its parity-preserving
apparatus retains its fermionic reference; no successful branch is
renormalized.

The [scalar and total-Yukawa instrument](../../model/coherent_defect_wave/chiral_joint/ground_scalar_response/PROOF.md)
uses the full source identities `N^2=N`, `{N,Y}=Y` and `Y^2=H^2-H`.
From the ground energy E0, the only inelastic scalar frequency is
`1-2*E0`, about `23.4611705884953` source units. Its raw retarded residue
is about `0.249545806772616`. The actual current obeys
`(i[H,N])^2=Y^2`. All 22 instrument records are retained: six have positive
weight and sixteen are exactly zero. Elastic records and both full-carrier
failure operators are included. These are local stationary responses;
their transfer to the interacting many-cell reference requires an explicit
operator-and-state comparison.

| Field role | Exact finite attachment | Required physical identification |
| --- | --- | --- |
| Charged leptons | Native L and e^c currents; quantum Higgs composite; certified periodic charged spectral line | Propagating relativistic charged sector, generations, mass thresholds and decay laws |
| Up and down quarks | Native colored chiral weights, all-occupation Higgs equations, joint full-family production | Colored dynamics, confinement and a stated renormalized mass convention |
| Neutral leptons | Gauge-singlet odd composite, complete local spectral response, neutral motion from full Peter–Weyl transitions | A common stationary reference and moving relativistic branch; flavor and mass/mixing law |
| Scalar/Higgs channel | Five-state quantum memory, actual invariant pair tensors and complete same-ground scalar energy records | Propagating scalar response and its relation to fermion and vector masses |
| Photon and weak vectors | Native electroweak currents, gauge-compatible links and complete same-ground closed-loop records | Vector-channel attachment, mixing, massless/massive dispersion and physical polarizations |
| Gluons | Native color currents, full compensated colored transport and matter-generated contractible loops | Non-Abelian phase and color-singlet observable consequences; no free-gluon pole is required |

No gauge-weight table by itself assigns observed particle identities. The
computed finite records and the conditional background field dictionary
have separate state and parameter ancestry.

### Other source branches and the physical join

The [neutral graded Wilson field](../../model/coherent_defect_wave/many_body_route/PROOF.md)
and its [complete mass instrument](../../model/coherent_defect_wave/mass_instrument/PROOF.md)
retain their positive free gap two, particle/antiparticle vacuum and explicit
work-battery records. Their five face channels are gauge singlets. The
[finite native encounter](../../model/coherent_defect_wave/propagating_native_join/PROOF.md)
retains two propagated objects and a positive joint interaction response;
its bounded-contact limit applies only to the stated one-flux/one-native
window. Neither branch becomes the chiral interacting field by matching a
free dispersion formula.

The [source-specific dependency path](../../model/coherent_defect_wave/PARTICLE_PATH.md)
requires a controlled physical refinement of its interacting reference,
nontrivial covariant dynamics, scalar/flavor response, physical mass/scattering
records and independent experimental identification. Finite current and
charge-conjugation theorems supply their exact interfaces; they do not
select those physical inputs. The [atlas](OPH_PARTICLE_ATLAS_V3.svg) separates
the proved source constructions, conditional field roles and measured masses.

## Local coherent seam operations

The [coherent seam observer](../../model/coherent_seam_repair/README.md) is an
explicit executable architecture. A seam correction exchanges
`|->|0_record>` and `|+>|1_record>`, where
`|+/->=(|10>+/-|01>)/sqrt(2)`. Removed disagreement enters coherent memory.
Committing that memory after each event defines a different, dissipative law.
Exchange antisymmetry of these amplitudes does not establish fermionic statistics.

The seam operation and a reflection of its address register give a
1440-dimensional unitary `W=Q(2AA*-I)`. Its blank-memory discriminant is the
exact native one-load kernel `I-L/60`; compressed amplitudes obey a finite wave
recurrence. The clock counts W steps. Its two-step return probabilities differ
from the committed repair-attempt chain.

The same W implements a complete preparation with
`K_success=(2I+W)/3` and `K_failure=sqrt(2)(I-W)/3`. Its success state on the
walk-register algebra, from the maximally mixed input, is
`rho_W=[(5/4)I+(W+W*)/2]/524`, faithful and stationary.
Success probability is `131/810`; k independent restarts have failure
probability `(679/810)^k`. Failures are retained. Conditioning on a preparation
archive changes the reference. Its modular generator is a function of
`W+W*`, rather than the walk's quasienergy generator.

Optional register projectors and declared quarter-turn probes give complete
complex tomography on this finite space. Apparatus interventions are separated
from idle propagation. Their tomographic linear span does not implement every
positive operator-valued measurement.

The [encoded transport coupling](../../model/coherent_seam_repair/PROOF.md#gauge-coupled-walk)
attaches a declared five-level fibre to the moving source. The resulting
7200-dimensional walk preserves `Q5=diag(2,2,2,-3,-3)` and has its own reference
`[(5/4)I+(W_G+W_G*)/2]/2620`, with the same preparation success probability.
A central binary-cover link produces two nontrivial triangular holonomies
and changes an actual three-step return amplitude. This is a finite gauge
background, propagation and observable response in one source and clock.
The background has no autonomous gauge-field evolution or backreaction.
Its flux is invariant under recharting and removable by a link change.

The [finite agreement theorem](../../model/coherent_seam_repair/PREMISES.md)
proves equality of complete adaptive record laws under implemented chart
changes, with the preparation, links and apparatus transformed together.
Preparation failures and untouched quantum references are included. It also
lists the exact A1, A2, A3 and operational Born interfaces that require further
witnesses; a coherent execution alone does not supply every finite premise.

A separate [boundary-flux construction](../../model/coherent_seam_repair/PROTECTION.md)
stores thirty central link bits. Its interior repair moves a face defect while
preserving the region's odd flux when its boundary is fixed. The defect reaches
all five faces of the region. A permitted boundary flip can erase the original
pair. This component proves protection under a specified move grammar. Its
coupling to the wave is the joint construction above; physical particle
identification requires the remaining spatial, dynamical and readout interfaces.

The [checked result](../../model/coherent_seam_repair/results.json) supplies
a finite kinematic coupling. A controlled spatial family, protected excitation
dynamics, operational energy/momentum, physical spin/statistics and the
detector dictionary are requirements for a particle identification. No
interacting 3+1-dimensional limit or physical rest mass follows from the
finite carrier. The dependency map below applies separately to each program.

## Particle candidates on the Wilson and comparison sources

The [labelled catalogue](../../model/particle_candidate_confirmation/CATALOGUE.md)
groups eleven protected-object entries and sixteen source excitation IDs into
families. Conserved charges, protected memory, spatial winding and kinetic
winding retain their distinct meanings and construction domains. The bare
source entry is an exclusion control. No protected-object entry is identified
with a named Standard Model particle.

The [per-object particle and mass paths](../../model/particle_candidate_confirmation/MASS_PATHS.md)
state the first missing attachment for each protected entry and identify the
existing finite results that supply part of its calculation. A construction
can be continued directly; combining its records with the primary Wilson
source requires a transfer proof. No architectural uniqueness proof is needed.
For the protected hidden-spin factor, `H=h tensor I` makes spin-only probes
energy-elastic. The eta isometry obeys `H_6 V=V(H_4+U)`, an exact two-occupation
energy shift. Neither statement supplies an elementary rest mass.

The [candidate confirmation plan](../../model/particle_candidate_confirmation/CONFIRMATION_STEPS.md)
attaches six tests to every elementary target: source operator and charges;
populated energy-momentum response; physical transformations and statistics;
interactions and flavor; the appropriate mass observable; and independent
experimental comparison. These specialize the existing dependency map.
Charged thresholds, unstable-particle resonances and renormalized quark masses
have distinct definitions. Unit conversion supplies no species identification.

| Source family | Verified calculation on the declared source | Physical attachment required |
| --- | --- | --- |
| Dressed charged endpoints with compensation | [Separated complete records](../../model/particle_candidate_endpoint_threshold/PROOF.md) at shortest cube path lengths `1,2,3` have original recorded probabilities above `4/25,2/25,7/250` in `[55/32,1337/32]`, `[239/32,2161/32]`, `[263/32,3145/32]` source units, respectively. Each is a fixed member-zero apparatus setting. | A charged-sector threshold or controlled compensator limit, together with a populated momentum-dependent branch. The overlapping pair-transfer windows include flux and the opposite endpoint; dividing them by two gives no endpoint mass. |
| Six signed neutral composite families | [Forty component tests](../../model/particle_candidate_composite_spectrum/PROOF.md) certify original neighbor-versus-nonneighbor records. For the weak/singlet trion's G1 and H multiplets, [original two-point sum rules](../../model/particle_candidate_neutral_dispersion/PROOF.md) bound each four-centroid spread below `5/9`, every centroid between `16/3` and `38/3`, and their ratio within one multiplet below `53/48`. | Resolve the signed neighbor coefficient and an original spectral branch, then attach operational momentum, physical currents and Lorentz spin. Spectral centroids and the separate four-point final-energy records are not pole energies or rest masses. |
| Neutral radial scalar | [The original scalar spectrum](../../model/particle_candidate_scalar_spectrum/PROOF.md) has positive transfer weight above the fixed source frequency `1/40`. Complete clocks retain weight at least `3cq/4` in `[3/160,D+1/160]`, with explicit finite resources. The same detector has exact elastic weight at least `3481/160000` across the family and greater than `7/10` at member zero. | A populated scalar branch in a controlled compact energy window, with momentum/dispersion and the appropriate interaction vertices. The inelastic floor is extremely small and the upper endpoint grows with capacity. Neither total radial detection nor the positive event identifies a Higgs mass. |

For the scalar row, `D=U-L` is the authenticated full-source energy enclosure
width, `q=9/(128 a^2 D^2)`, and `c=2^(-N)` is the original density floor;
`N` is a density-bound exponent, not capacity in nats. The
[total radial detector](../../model/particle_candidate_scalar_gauge_response/PROOF.md)
and the inelastic frequency cut use the same preparation and source clock.
The positive inelastic estimate coexists with the displayed elastic floor.

The [local charged filter](../../model/particle_candidate_charged_propagation/PROOF.md)
also retains its distinct `3/80` original event in `[311/32,1081/32]`.
Admission of its local seed to the canonical energy/momentum algebra does
not transfer this weight to the inactive member-zero bulk aperture.

All failures, compensators and quantum references stay in their respective
complete records. The neutral neighbor-moment estimator has deterministic
bias at most `1/2048`, using approximately `9.44e21` bins per energy clock;
its value and sign are uncomputed. Scalar clock sizes contain the enormous
density-floor exponent. These are finite analytic and executable resource
certificates, without sampled full histograms or an established compact
native-observer apparatus implementation. Physical species and source-derived
GeV outputs are empty in the
[checked confirmation receipt](../../model/particle_candidate_confirmation/results.json).

The [particle explainer](OPH_PARTICLE_ATLAS_V3.svg) ([PNG](OPH_PARTICLE_ATLAS_V3.png))
and [candidate atlas](OPH_PARTICLE_CANDIDATE_ATLAS.svg) show the family labels,
finite results and physical confirmation requirements. Their maintained
renderer is [the atlas package](../../model/particle_atlas_explainer/README.md).

## Particle formation and observer time

The [catalogue's cosmological comparison](../../model/particle_candidate_confirmation/CATALOGUE.md#early-universe-comparison-and-source-time)
separates a hot thermal population from the vacuum particle mass table.
The source preparation is stationary under its Hamiltonian: for a fixed
source effect `E`, `[rho,H]=0` implies
`Tr(rho exp(iHt) E exp(-iHt)) = Tr(rho E)`. Conditional histories and
retained-reference records can change. A first click does not establish
a first production event, and a refinement index is not elapsed cosmic time.
An early-universe comparison needs an initial preparation, temperature and
observer-time dictionary, and complete production/annihilation records.
Neither mass ordering nor a defect's catalogue position supplies that history.

The [unified field-and-spectrum dependency graph](../../model/particle_joint_reconstruction/DEPENDENCIES.md)
and its [machine-readable acceptance conditions](../../model/particle_joint_reconstruction/dependencies.json)
are the single ten-node construction map. The active target is source-derived
particle masses in declared SI/GeV comparison units, with capacity dependence
explicit. Nodes 1-3 establish the finite program, physical preparation and
complete local apparatus. Nodes 4-8 include a common positive-energy
algebra of filtered probes and local phase operations, finite complete
energy records and transfer-to-excitation transport. A retained observation
origin gives a translation-invariant
relative-observation limit with commuting energy and momentum generators.
Nonzero momentum response, dispersion, species and mass identification
require further spectral control. The primary physical-mass table in GeV
is empty. Pure Yang-Mills
(node 9) and geometric gravity (node 10) are deferred branches. Comparison
results retain their own sources and do not add canonical dependency lists.

The [physical mass readout criterion](OPH_NATIVE_PARTICLE_MASS_JOIN.md#physical-mass-readout-criterion)
states the original joint energy-momentum weight test, charged-threshold
and neutral-multiplet alternatives, and the exact refinement obstructions.
Internal source labels and finite transition energies alone do not identify
particles or their relative physical masses.

The objects produced by the executable are finite simulation-generated
observables: events, apparatus records and observer comparisons of their
interpretations. Agreed records support effective spacetime and field
models. A smooth continuum is a mathematical limit or approximation with
stated error bounds, not an additional objective universe emitted by the
program. Agreement alone does not prove the dimension, relativistic
symmetry or dynamics of that effective description.

The [finite source/output construction](../../model/particle_observer_outputs/README.md)
connects complete spectral readouts, the faithful modular clock, source
coupling recovery, energy-accounted holonomy transport and comparison of
observer meanings. The [coefficient inverse](../../model/particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md)
removes an additive sector offset from observable family-coefficient
ratios in the prepared Q*=0 experiment. Its two-line inverse is not the
multi-line response of the full Q*=6 ground. Physical rest-mass
identification requires the interacting
spatial response of those same operators; no measured particle assignment
is inferred from the finite coefficient reconstruction.

## Mass outputs and the immediate physical target

The immediate target is a source-derived mass output in GeV. Geometric
gravity and the pure Yang-Mills gap branch are separate targets; they are
not prerequisites for a first particle comparison. The fixed source,
physical preparation, complete finite instruments and common bounded
filtered algebra are available. The following work order expands nodes
4-8 of the canonical map; it does not add another set of acceptance gates.

| Order | Required mass-readout connection | Existing result and remaining condition |
| --- | --- | --- |
| 1 | Resolve an original-weight excitation in a declared physical measurement domain. | Three conjugate native-occupation sector pairs have full positive thresholds on member zero. The original spin-resolved instrument reports weights above 121/512 (G1) and 61/128 (H) in [11/4,61/4]; the source-centered filter separately reports above 3/5 there. The sector resolvent and low-band weight concern a prepared trion; source faithfulness also gives a capacity-deteriorating positive original final-energy event. These are distinct probability bounds. An isolated branch or threshold onset with a retained residue and same-channel dispersion requires further calculation; a continuum claim also needs fixed-energy survival. |
| 2 | Attach momentum to the same source, state and clock. | One retained observation origin supplies commuting canonical H>=0 and P, with complete regular-probe detector kernels. The complete signed odd analyzer populates every member-zero Walsh row in both G1 and H on the original density. Actual phase-dependent energy laws, a resolved branch and a same-channel connection to canonical P or an error-controlled finite momentum readout require proof. Prepared cubic transport retains its separate domain; fixed-origin invariance is not assumed. |
| 3 | Extract a relativistic rest threshold, stable branch or resonance. | Original transfer-to-final-excitation transport and conditional rest intervals are proved. The same momentum-resolved excitation needs a controlled energy-momentum relation and retained spectral weight over the claimed measurement domain. A finite comparison may bound deviations directly; an exact continuum claim requires its limiting theorem. |
| 4 | Identify the propagating particle channel. | Each odd physical gauge-singlet space has six source spin states, decomposing as binary-cube G1 plus H; both components have exact nonzero cubic transport, and the weak creation apparatus has populated original G1/H records. Native-current commutators, graded exchange at distinct sites and opposite-sector charge conjugation are explicit. Lorentz spin, asymptotic particle statistics, physical currents and interactions must attach to the same spectral branch before assigning a Standard Model name. |
| 5 | Bound the complete prediction error and capacity dependence. | Finite preparation, filtering, clock and detector bounds are constructed. Joint volume, spacing, link-cutoff and observable-matching errors require control. The chosen sequence changes these together, so capacity alone does not determine a mass law. |
| 6 | Apply one explicit SI/GeV comparison scale. | Available immediately through a declared absolute length or clock unit. Keep that scale fixed or explicitly tracked across refinement; deriving a cesium atom is not a prerequisite for this convention. A scale fitted from a measured mass retains that calibration ancestry. |
| 7 | Compare independent observable masses and ratios. | The primary physical-mass table in GeV is empty. Stable masses, charged infrared thresholds and unstable resonances require their respective experimental definitions. Inputs used to fit a scale or coefficient are not independent predictions. |

Spectral response (1), operational momentum (2) and channel labels (4)
admit independent calculations on the fixed construction. Each carries its
own finite error bounds into (5). Their common branch is the input to the
rest-energy calculation (3). Unit conversion (6) is available independently;
the comparison (7) consumes the identified branch, its errors and that
declared conversion. No particle name or measured mass is used to certify
the preceding connections.

The required physical joins are a resolved primary response, same-state
momentum and spectral branch identification, with errors on the stated
domain. A first finite-capacity comparison does not require an exact
infinite-volume continuum theorem. The broader continuum claim requires
the corresponding uniform response and refinement proofs. Unit declarations
can accompany either calculation directly. Mass identification joins their
results on the same source. A failure to
match a measured ratio must first be localized to that chain; it does not
justify modifying coefficients using the target and relabelling the fit
as a derivation.

## Positive sector thresholds and conjugate source excitations

The [conserved-sector construction](../../model/particle_primary_sector_threshold/PROOF.md)
uses the exact occupation increments relative to the filled-beta anchor,
in native charge order `(-4,-3,1,2,6)`. The full Hamiltonian preserves
each increment. Full-occupation electric inequalities and exact physical
Gauss-singlet trials bound the bottom of each indicated sector above the
full source ground:

| Sector pair | Physical local gauge contraction | Source CAR parity | Lowest excitation interval, source units |
| --- | --- | --- | --- |
| `±(0,2,0,0,1)` | Weak-doublet epsilon pair and charge-six singlet | Odd | `[1.54674,12.39802]` |
| `±(1,0,0,2,0)` | Charge-minus-four and two charge-two color antifundamentals, contracted by epsilon | Odd | `[0.56642,12.72435]` |
| `±(0,0,0,3,-1)` | Charge-two color determinant and charge-six hole | Even | `[3.56682,15.71559]` |

Every state in a listed sector obeys its lower bound. The weak-sector
upper uses the [43-coordinate virtual-pair certificate](../../model/particle_primary_branch_certificate/PROOF.md);
the other uppers use exact 13-packet interacting Rayleigh quotients.
No invariance of these trial spans is assumed. These are member-zero threshold enclosures,
not rest masses or a positive global gap. Their six sectors are physical
reducing subspaces; their names describe source contractions and do not
assign measured particles.

A complete preparation makes each sector operational. An even joint
Householder prepares opposite sectors in two idle physical memories.
A graded source/memory exchange with its required parity correction moves
the original state and all reference correlations into retained memory;
the second memory keeps the compensating increment. Combined native
occupations and parity are preserved. One unchanged source-energy clock
then has original probability at least 255/512 in an interval wholly above
the ground upper. Ready-memory capacity and source work are explicit.
This changes the source marginal and does not assign that probability to
its original near-ground sector population.

The [source charge-conjugation theorem](../../model/particle_primary_charge_conjugation/PROOF.md)
gives more than equality of trial expectations: at every finite member,
combined full-CAR particle/hole transformation and gauge-coordinate
conjugation are a linear unitary preserving Gauss support, cutoff and H.
It maps `Delta N` to `-Delta N`, so the full opposite-sector spectra,
multiplicities and conjugate complete energy laws agree with the same
energy zero. The polynomial preparation is invariant; the canonical
reference transforms only in the comparison identity. Proper binary cube
rotations commute with this symmetry, while odd conjugate inversion eigenstates have
opposite inversion signs. This supplies source transformation information,
without a Lorentz or observed-antiparticle identification.
For the centered local density, `C(N_x-30)C†=-(N_x-30)`, so the two
single-probe phase signs have equal original energy laws at every origin
on the unchanged selected state. Their sign-averaged bounds therefore
also hold for either deterministic sign after the same origin average.
This does not permit flipping one sign inside an arbitrary conditional
history while leaving its other operations unchanged.

The same symmetry preserves each Hermitian source-action term separately.
For the fixed [regional energy allocation](../../model/particle_primary_stress/PROOF.md),
it follows that `C h_x C†=h_x` and
`C i[h_x,h_y] C†=i[h_x,h_y]`; C is linear and does not conjugate i.
The specified square-root energy/current instruments therefore give equal
complete classical histories on conjugate input states, with unitarily
related retained quantum outputs. Arbitrary extra noncovariant measurement
back-action is outside this comparison. This supplies a finite regional
energy property; proper Lorentz stress, geometry and Einstein dynamics
retain their independent premises.

A defensible finite GeV mass requires a populated branch or threshold
with physical momentum and controlled dispersion in one of these same
apparatus domains, followed by species/current matching and the declared
unit conversion. The positive source thresholds do not complete those
connections. The primary physical-mass table in GeV is empty.

## Finite composite spin, propagation and original-state creation

The [full local composite construction](../../model/particle_primary_composite_spin/PROOF.md)
expands both beta-spin rows in the actual gauge kernels. Each odd sector
has six states with source spin-Casimir eigenvalues `3/4` and `15/4`, of
multiplicities two and four. Exact physical binary-cube action gives
`G1 + H`; the conjugate sectors have opposite intrinsic inversion signs.
The even four-defect composite has eight states, with `T1 + E + T2`.
Native-current commutators and graded exchange between distinct sites are
explicit. Uniform orbital preparation makes these full-source cube
components, each with a positive sector threshold and complete reported
energy-window weight at least `255/512`. These are finite symmetry labels,
without a Lorentz little-group or Standard Model identification.

The [interacting composite transport](../../model/particle_primary_composite_transport/PROOF.md)
uses that same source at interaction coefficient one. For neighboring
localized odd triples, the first two transport moments vanish and the
third is exactly `-3/4`. On the entire six-state local gauge kernel it is
`-3 I6/4`, so both G1 and H components propagate. The degree-one cutoff
does not change this cubic moment because both endpoint first actions
are admitted. For the one-spin-row witnesses, the full-source Taylor remainder
proves prepared transition probability at least `t^6/256` for `0<|t|` up to about `4.945e-8` source
time units. All outcomes of the specified global bare-vector detector,
the compensating sector and the displaced source/reference are retained.
This is a conservative propagation proof, not a measured transport rate.

Eight coherent Walsh preparations have common mean `-711` and variance
`132301/735`, but centered third moments
`369869271/42875-(3/4)(3-2|k|)`. Their complete energy laws are therefore
different. Neither the Walsh label nor this moment is an energy eigenvalue
or conserved momentum. Exact electric-denominator third-order transport
coefficients about the decoupled action are also available; their
certified weak-interaction domain excludes the primary coefficient one.
No extrapolated band is used in the mass calculation.

The [spin-resolved original-state apparatus](../../model/particle_primary_selected_spin_response/PROOF.md)
uses all six weak-pair/singlet creation polynomials on the unchanged
selected density. Its complete seven-outcome instrument at each address
retains the compensating memory, failure and reference. At member zero,
reported original weights exceed `121/512` for G1 and `61/128` for H in
`[11/4,61/4]`. Throughout the source family, their reported weights
exceed `1/200` and `3/200` in `[1/(4a),71/(4a)]`. Five coherent
native-number meters record the same
`Delta=(0,2,0,0,1)` without measuring the unknown initial numbers.
These are probe representations, not total-spin labels of the mixed
output. Cube symmetry fixes equal component weights after the address
average.

The [original-state branch apparatus](../../model/particle_primary_original_branch/PROOF.md)
uses the exact selected polynomial slope to bracket creation and adjoint
reverse probabilities by `exp(-beta omega)` and
`exp(-beta_low omega)`, with `beta=22272` and `beta_low>22271` at
member zero. No Gibbs state is substituted. In `[1/4,71/4]`, each
chosen flat spatial phase has reported creation weight above `1/3`,
including G1 above `5/48` and H above `11/48`. The complete Walsh
probe retains 96 successes and eight failures before clocks and meters;
every spatial row has original weight above `1/24`, with G1 above
`5/384` and H above `11/384`. Adjoint annihilation and source-conjugate
hole creation are distinct experiments.

Complete native-number and total-cube projectors refine these energy
records without changing their classical energy marginal. They retain
all initial and final labels, with number increment `(0,2,0,0,1)`.
The uniform row permits 28 G1g and 44 Hg total-cube transitions. Specific
original branches from `(n0,A1g)` to `(n0+Delta,G1g)` and
`(n0+Delta,Hg)` have reported final-energy weights above `85c/4096`
and `85c/2048`, respectively, where `c=2^(-582407133)` at member zero.
Here `n0=(48,32,96,48,16)`. These are source-marginal positive-effect
bounds using an anchor witness,
without replacing the preparation or its reference. They do not assert
an initial-to-final energy-transfer interval. The same conserved labels
contain orthogonal source vectors with means `-711` and `-699`, proving
spectral diameter at least 12. Labels alone therefore do not select a
mass. Actual inversion parity is retained; it is not Lorentz spin. For every
normalized coherent component, the actual anchor residual gives original
normalized transfer variance at least `9c/4`; each full measure therefore
contains at least two frequencies. This tiny finite bound gives no line
positions or clock-resolved separation and permits an isolated pole with
additional spectral weight.

The [native odd-field sum rules](../../model/particle_primary_odd_sum_rules/PROOF.md)
retain creation and reflected adjoint-annihilation histories. Their
combined success effect is identical in every full Walsh row, and the
first spectral shell moments reduce exactly to `(m0,m1,0,0)`. The
[occupied-pair constraint](../../model/particle_candidate_neutral_dispersion/PROOF.md)
on the original weak/singlet trion gives `abs(m1) <= d delta/64`, where
`d=2,4` for G1,H and the actual defect bound has `delta<5`. Its value and
sign are uncomputed. The four normalized two-point centroids form an
arithmetic progression with spread below `5/9`; this is not a resolved
dispersion branch or a bound on particle mass ratios. A separate complete
four-point graded experiment distinguishes neighbors from nonneighbors
with original probability difference at least `c t^2/128`, at
`t=c/(64 R_H^2)`. This is a tiny finite propagation record, not an odd
two-point dispersion or a particle mass.

On every selected capacity, exact spatial inversion exchanges the native
`Gamma5` projectors while preserving their gauge representation and energy
laws. The proof supplies complete Wilson-dressed bilinear and clock
records, with all failures retained. A parity-equivariant identification
with physical handedness cannot discard one reflected low-energy response;
parity breaking or a different physical current dictionary needs an
explicit construction. Native chirality is not proved Lorentz chirality,
and paired records need not be distinct particles. The local weak and
color trion creators are gauge singlets, so their occupation labels do
not identify charged leptons or quarks. Charged regional fermions need a
compensating dressing or boundary/memory construction.

The [complete odd filter](../../model/particle_primary_odd_low_energy/PROOF.md)
uses the original four Slater probes with a triangular duration `a/6`
and analyzer frequency `9/a`. A finite coherent selector and the same
source clocks retain original reported blank weight above `3/5` in
`[11/4,61/4]` at member zero and at least `1/80` in
`[1/(4a),71/(4a)]` throughout the family. The analyzer frequency does
not shift the source energy zero. Summing all retained selector outputs
recovers the raw probe's energy law exactly. A useful successful filtered
field is permitted, but its own low-energy weight requires a lower bound.
The proved spin-resolved and source-centered bands scale as `1/a`.

The [full odd-sector resolvent](../../model/particle_primary_odd_resolvent/PROOF.md)
is bounded at the actual interaction coefficient one on member zero.
For the prepared local weak trion, the complete complement within its
reducing occupation sector gives two-sided self-energy and resolvent
enclosures below the rigorous sector floor `L=-722.22108`. Exact moments also give prepared spectral weight
above `1/26` in `[L,-714]`, and one complete source clock retains above
`1/27` in `[L-1/64,-714+1/64]`. These are absolute source energies with
the unchanged zero. The preparation moves the original state/reference
into retained memory. The actual source-marginal faithfulness floor
`rho>=c P_G` also transfers its implemented positive final-clock event:
the original fixed local odd probe has weight above `c/27`, or above
`c/864` for one retained row in the member-zero `4V` probe menu.
This event measures final energy without initial energy recording. It
neither transfers the prepared probability wholesale nor identifies an
original transfer-energy law. The authenticated member-zero floor is
`c=2^(-582407133)`. Its positivity supplies neither a practical event rate
nor a uniform continuum residue. Neither this interval
nor the resolvent bounds isolate a particle branch. The same parent fourth
moment improves the specified one-spin-row weak-composite transport guarantee to
`Pr(transfer)>=t^6/256` for `0<|t|<=2.77366e-6`.

The [full-sector virtual-pair certificate](../../model/particle_primary_branch_certificate/PROOF.md)
uses 43 admitted real trial coordinates, including two distinct quantum
links, to lower the full weak odd-sector upper below `-715.10018`.
The bound uses the complete source operator, not an invariant compression.
An energy-local polynomial floor
`c_low=2^(-b_D-ceil(2 beta(U_Delta-L_0)))=2^(-552953)` gives original
probability at least `2 c_low(1-delta)` per signed occupation sector, or twice that
for their union, in the widened final-energy band around
`[L_Delta,U_Delta]`. The finite clock loss `delta<1/512` is applied
before the tiny density factor. There is no creation probe, initial
energy clock or success normalization in this event.

The [candidate correspondence](../../model/standard_model_source_map/candidates/CANDIDATES.md)
compares 16 source-excitation rows and 11 existing protected-catalogue rows
with 18 named elementary targets. Exact gauge weights motivate family
hypotheses; they do not assign flavors, Lorentz spin or particle masses.
The actual `Q6=+6`, `G1` endpoint is a charged-lepton-conjugate-singlet-type
lead, shared by the electron, muon and tau hypotheses. The protected
catalogue and source-excitation identifiers remain distinct. The
[common-ground ratio calculation](../../model/standard_model_source_map/candidates/FINITE_SECTOR_RATIOS.md)
bounds every ordered ratio of the three member-zero neutral sector
thresholds below21.16, and the two odd-family ratios below15.31.
One common unit conversion cannot turn those particular thresholds into
the measured charged-lepton hierarchy. This excludes that direct mapping,
not other source excitations, capacities or dictionaries. The physical
particle and GeV-mass maps remain empty.

The [isolated gauge-sector transition](../../model/particle_primary_reduced_sector/PROOF.md)
is an actual source atom at member zero. In the empty-matter reducing
sector, the twelve-level first gauge cluster contains a simple `A1u`
level within that symmetry. Its transition from the sector gauge ground
lies in `[63-59/1351350,63+1/45675]` source units, with squared Wilson-loop
matrix element above `4/49`. The original selected transfer atom retains
its rare-sector weight; the complete clock event exceeds
`2^(-32406376)/13`. This is CAR even, and the initial gauge ground is
not the global primary vacuum. The extreme-sector weight vanishes along
refinement, so this exact finite line is not a surviving particle or a GeV
mass. It demonstrates full-complement isolation on an actual reducing
carrier without replacing the source state.

The [loop-fusion calculation](../../model/particle_primary_loop_fusion/PROOF.md)
retains the full empty-sector complement. Adjacent outward charge-six
loops cancel their shared edge and populate an isolated 44-level cluster
in `[1311/14,1335/14]`; disjoint opposite loops have a complete window
`[124,128]`. Their static energy deficit relative to twice the isolated
line lies between `30.642769` and `32.357165` source units. Each row of
the complete 30-pair menu retains original clock probability above
`2^(-32406376)/40`. These are static sequential loop-creation records,
not moving collisions or a physical binding mass; their rare-sector
weight vanishes on refinement.

The [separate collision controls](../../model/contact_scattering/README.md#a-complete-source-operated-interaction-experiment)
have 17 replayed regression tests. Their full four-mode, two-site
Weyl/contact experiment executes 22,531 native pulses and 45,062
load/export events, with every failure and the quantum reference retained.
It uses declared `c=1/8`, `U=1` and duration4 on a different action and
preparation from the primary Wilson source. Its nonzero interaction echo
is an executed interaction control, not a CERN cross section or a
collision of identified primary particles. The
[moving-packet family](../../model/scattering_packet_readout/PROOF.md) has a
controlled nonrelativistic collision detector; its fine-spacing packet
experiment is distinct from its executed finite five-qubit control.

The [regional fermion instrument](../../model/particle_primary_regional_fermion/PROOF.md)
retains a charged `G1` endpoint, its opposite-charge endpoint, their Wilson
link and the original reference. The complete Pauli pair analyzer has
original recorded weights above `5/64` and `15/64` in
`[69/8,279/8]` for scalar and vector coefficients. An alternative complete
96-row raw-spin menu records actual endpoint spin increments `+/-1/2`
and charge increments `+/-6` without initial spin or charge measurement.
Each of its four endpoint-spin entries has recorded weight above `1/16`
in `[77/8,271/8]`. These are source transformation and probe-time increment
labels. The pair intervention's energy is not an isolated endpoint mass;
Lorentz spin, a charged one-particle dressing and dispersion need proof.

Every energy in the weak-trion native occupation sector has binary-cube
multiplicity at least
two. A rank-one complement in the entire sector therefore retains a sector
ground vector and cannot have a raised bottom. Exact `G1`, `G2` and `H`
matrix units remove these compulsory partners while preserving the source
clock. The remaining resolvent estimate belongs to one diagonal irrep
multiplicity block. An isolated multiplet, its original creation weight,
dispersion and physical mass require that estimate and their own attachment.

The retained scalar onsite/electric-penalty comparison has a proved
complement-floor ceiling below `-715.4069`, below the authenticated weak
trial upper below `-715.10018`. That comparison cannot isolate the trial
band. This bounds the method, without excluding an actual isolated
sector; a connected or dressed complement requires its own estimate.

The [selected odd residual](../../model/particle_primary_selected_odd_response/PROOF.md)
places more than `1/50` original weight at energies of order `1/a` and
forces the positive-frequency first moment and unnormalized output
excitation mean above `1/(100a)`. A coexisting low-energy component is
permitted. The complete filter calculation bounds the raw odd fixed-band
limsup by `1/9`; its ideal zero-centered duration-one blank has weight
below `1/4` at every member. Neither upper bound proves a surviving residue. Static
number-phase conjugation only changes each retained local odd row by a
scalar phase. Even time-separated critical density-Weyl dressings of
depth `o(a^-1/2)` cannot make the complete odd output energy tight.

The [original density-shell calculation](../../model/particle_primary_original_shell_moments/PROOF.md)
computes an exact off-diagonal energy moment on this selected state in the
separate even sector `n_x=N_x-30`. At member zero,
`M^(1)=(tau/2)(3I-Adj_cube)`, with
`tau=-Tr(rho H_hop)/12>0.2373`. The complete sixteen-record density
analyzer has Walsh first moments `s tau/7200`, for `s=0,1,2,3` with
multiplicities `1,3,3,1`. Actual finite clocks retain strictly ordered
nonzero original means. The constant row is elastic and has positive
weight; the other rows have positive-frequency weight, so their normalized
energy shapes differ from the constant row. This evaluates even density
moments, not the four odd G1/H measures. The latter's normalized spectral
shapes and a common resolved branch require their own calculation.

The same calculation resolves each native internal density: its bond
coefficient obeys `tau_r>=2 d_r c beta_low>0`, with the actual faithful
source floor `c`; this proves native inelastic density response, not
creation of a named fermion. All total-density successes have combined
original probability below `0.0208`. Conditioning on those retained
successes is explicit and gives signed transfer mean above `0.019`, or
above `0.015` with the finite clocks. These conditional means are not
spectral lines or rest masses.

The [exact cube-kernel reduction](../../model/particle_primary_coherent_odd_response/KERNEL_PROOF.md)
reduces each component-summed G1/H response to four real spectral measures,
indexed by cell-pair Hamming distance. Off-diagonal measures may be signed.
Their four Walsh combinations have multiplicities `1,3,3,1`; this is a
spectral Gram decomposition, without an H-eigenstate or momentum assignment.
The next calculation is those measures on the original selected density.
Their actual values, normalized energy-shape differences and a common
branch or threshold are not fixed by symmetry or by prepared-vector
transport. Intensity contrast alone is insufficient for dispersion.
Operational momentum and controlled dispersion are the physical mass
connections. A continuum claim also requires fixed-energy original weight bounded away from zero under
refinement; the finite faithfulness floor and filtered upper bounds do
not supply that estimate.

For the six odd polynomials, the original kernel is

\[
 G_{xy}^{jk}(t)=\operatorname{Tr}\rho\,
 F_{x,j}^{\dagger}e^{itH}F_{y,k}e^{-itH}.
\]

After summing j within one multiplet lambda, proper cube rotations and
pair reversals give

\[
 \sum_{j\in\lambda}G_{xy}^{jj}(t)
 =\int e^{it\omega}\,dm_{d(x,y)}^\lambda(\omega),\qquad
 \begin{pmatrix}\Lambda_0\\\Lambda_1\\\Lambda_2\\\Lambda_3\end{pmatrix}
 =\underbrace{\begin{pmatrix}
 1&3&3&1\\1&1&-1&-1\\1&-1&-1&1\\1&-3&3&-1
 \end{pmatrix}}_{K}
 \begin{pmatrix}m_0\\m_1\\m_2\\m_3\end{pmatrix},\qquad K^2=8I.
\]

Here d is the cell-pair Hamming distance. Each Lambda_s is a positive
spectral measure, with multiplicities `1,3,3,1`. The chosen-phase creation
law is `Lambda_s/6`; one row of the complete Walsh menu is `Lambda_s/48`.
Their complete sum is `m_0/6`, and the four representative menu-row
measures p reconstruct `m=6Kp`. The same covariance reduction holds for
implemented clock events, retaining their approximation error. The
full component-resolved and cross-irrep kernels are outside this scalar
reduction, and time transforms of the real signed measures can be complex.

Nonzero off-diagonal measures are equivalent to nonidentical unnormalized
Walsh measures, but can change only their intensities while leaving every
normalized energy shape identical. Their source-specific values must
therefore establish an energy-shape or branch dependence before a
dispersion claim. An explicit separate faithful, stationary polynomial
example has nonzero prepared cubic transport and vanishing original
inter-site measures even under arbitrarily strong cooling. It proves
why prepared transport cannot replace this original-state calculation;
it is not a computation of the primary source. Branch identification,
physical momentum and the rest-energy relation consume this common kernel
with their stated finite or refinement errors.

The [critical-profile work calculation](../../model/particle_primary_critical_refinement/PROOF.md)
keeps the actual source and bounds its selected hopping coefficient by
`h_n>=h_*+3a^2/7-3a^4/28-a^8/128`, with `h_*>0.0424`.
For fixed physical tents, `W_n/h_n` approaches the three-dimensional
Dirichlet moment `(pi^2 s^2/6) integral |grad f|^2`, with explicit
`O(a^2)` lattice and clock errors and `O(a^4)` phase errors. Work ratios
cancel the unknown coefficient. This is an energy-weighted response;
a purely elastic probability limit remains possible. The fixed-anchor
relative-perturbation decomposition has bounded remainder norm at least
`sqrt(15)` per Wilson bond and supplies no uniform interacting gap.

The [original work-spectrum calculation](../../model/particle_primary_work_spectrum/PROOF.md)
instantiates the finite fluctuation theorems on the actual energy-index
marginal, with a symmetric doubly stochastic kernel. Its entropy variable
`log p_i-log q_j` differs from modular work `Delta K`;
`E exp(-Delta K)=1` and the source affinity is `E exp(-Delta K/2)`.
Paired source transitions sharpen the negative first moment to
`1/(3 beta_low)`. The [selected spectral-tail bound](../../model/particle_primary_weighted_spectral_tail/PROOF.md)
uses the exact stationary third-moment identity and the actual electric
expectation `<H_E><91 V/a` to prove `M2<=4,000,000/a`, independently
of link-cutoff degree. The positive first moment and negative-moment
bounds close a quadratic inequality for the original second moment;
no reference-state domination or trial-state expectation is substituted.
For the positive work floor `w_n` and certified second-moment cap `C_n`,
the original band `[w_n/4,4 C_n/w_n]` has probability at least
`w_n^2/(4 C_n)`. The unchanged complete clocks retain more than half
this floor in a slightly widened positive band. The guaranteed floor is
of order at least `a`, while the upper endpoint is `O(1/a)`.
Uniform fixed-energy survival requires further control of the selected
current double commutators; this bound permits rare cutoff-scale weight
carrying finite work. Source affinity and complete record distinguishability
also have positive finite bounds. The [equal-work parity control](../../model/particle_primary_equal_work_control/PROOF.md)
uses this same source, state, reference and H with a different complete
apparatus: apply local parity `(-1)^N_x` with probability `p_n=kbar/4`,
and otherwise idle, retaining the origin and coin. Its selected work
matches the critical probe exactly. Yet `p_n<=45a/32` for `n>=1`, its
original law is `(1-p_n)delta_0+p_n nu_parity`, and its second moment
is bounded on both sides by constants divided by `a`. The same clocks
retain positive mean work while every fixed positive-energy band loses
probability. Fixed-profile spatial work ratios are also reproduced.
This is an actual-source insufficiency result, not a collapse proof for
the original critical phase apparatus. Its selected current correlations
or original fixed-band detector probabilities require evaluation.

The [original local-energy response](../../model/particle_primary_energy_response/PROOF.md)
uses `A=h'_0` and the actual incoming current `B=i[H,A]`.
A full-CAR physical witness and the selected polynomial slope give
`Tr rho[A,[H,A]]>=c beta_low 27/(2a^4)>0`.
An explicit finite Cayley control pulse and both binary current outcomes
therefore have strictly positive original probability contrast and cyclic
work, with the reference retained. The finite source-relative derivative
is generated by the integrated same-source energy, and its coarse Ward
law retains the applied-control work term. These bounds deteriorate with
the source faithfulness floor. Identifying the control as a proper lapse,
or deriving a geometric modular kernel and Einstein response, requires
additional physical attachment.

The [cube-covariant energy allocation](../../model/particle_primary_covariant_stress/PROOF.md)
shares complete edge terms between endpoints and plaquettes among their
four vertices. Its density sums to the unchanged H; the difference from
the owner allocation is an explicit nearest-edge divergence. The current
`i[h'_x,h'_y]`, complete instruments and applied scalar-control work are
cube covariant and instantiate the finite energy-Ward sums. Charge
conjugation preserves them and annihilates the first regional entropy
response of admissible odd tangents. Proper geometry, Lorentz stress,
the diamond modular kernel, area response and matched continuum errors
are separate gravity premises, outside the first GeV mass comparison.

The [induced regional cut](../../model/particle_primary_geometric_cut/PROOF.md)
uses regional matter and wholly internal links, with crossing links in
the complement. Its actual marginal, boundary sectors and modular flow
are cube covariant. The representation entropy satisfies
`0<=S_edge<=ell b_R log(6)`, where `b_R` counts crossing links; this is
not a total entropy or finite area-coefficient theorem. A contained
regional energy differs from the covariant allocation by an explicit
boundary operator. Its modular generator is not identified with that
energy. Proper geometry and the matched Einstein estimates are separate
physical premises.

The [source boundary-information identity](../../model/particle_primary_boundary_information/PROOF.md)
keeps the induced cut, its central sectors and their global coherences.
For source marginals it proves
`I_ambient=I_alg+2 S_edge+C_center`, where
`C_center=D(rho||pinch(rho))`. The actual magnetic action forces
`C_center>=10(c beta_low/(210a))^2>0`; the complete center record costs
source work at least `20 c beta_low/(210a)^2>0`. The exact polynomial
slope also gives
`beta_low W<=C_center+D(pinch(rho)||rho)<=beta W`.
The regional-algebra information has the finite floor
`I_alg>=10 N c log(2)`, with the actual `c=2^(-N)`.
These floors deteriorate with capacity. The reference-inclusive pure-state
quantity is different, and none of these entropies is identified with
proper area, a geometric modular kernel or an Einstein coefficient.

The [Gauss-compatible boundary comparison](../../model/particle_primary_boundary_correlation/PROOF.md)
retains the physical representation singlets and both source regional
marginals while removing conditional multiplicity correlations and center
coherence. Exact relative entropy gives
`sum p_lambda I(tau_lambda)+C_center <=2 beta ||V_cut||+r`, where
`r=beta delta(U-L)` is the actual polynomial modular defect. Only crossing
hopping and mixed plaquettes enter `V_cut`, with
`||V_cut||<=(30 b+p_mix/7)/a` and `p_mix<=2b`.
Together with `H(p)<=ell b log(451)` and `S_edge<=ell b log(6)`, this
bounds full source mutual information by a boundary count plus `r`.
A noninteracting polynomial control proves that `r` cannot be omitted.
The coefficient diverges along refinement; proper area, a geometric
modular kernel and Einstein dynamics remain separate premises.

## Populated channels across the finite source family

The [pair-channel construction](../../model/particle_primary_pair_channels/PROOF.md)
uses four Pauli components in the native dimension-one charge-six block.
Each local contraction obeys `[C†,C]=I-D6/2`. The physical cube action and
beta inversion give local operator types `A1u + T1u`; the 32 site/component
coefficients decompose in each parity as `A1 + A2 + E + 2T1 + 2T2`.
All 32 coherent Walsh/component records have positive original transfer
weight. In the alternate cubic/parity apparatus, each isotypic component
of rank d has original weight at least `11d/1024` in
`[11/(2048 beta),40]` source units. The complete instrument retains 32
successes, 32 failures and both finite energy registers. The source marginal
is stationary; the initial energy readout changes its purification coherences.
The source matrix and complete spatial histogram are not executed.

These are actual source operator channels. The rank-one A1u mode is the
zero-Walsh scalar, while zero-Walsh vector components occupy only part of
the rank-six T1u space. Cubic scalar and vector labels do not fix continuum
spin, and a zero Walsh label does not establish a rest frame. The source's
exact polynomial population slope bounds low frequencies; no Gibbs state
or measured mass is substituted to obtain the positive windows.

The [refinement response](../../model/particle_primary_refinement_response/PROOF.md)
uses a physical eight-cell trial tiled through larger sources and an exact
full-occupation local electric estimate. The trial gives a ground-energy
upper bound without changing the selected polynomial preparation.
Cross-tile hopping vanishes by fixed tile occupations; cross-tile
plaquettes vanish by Haar integration over crossing links. On every
primary member the selected total beta-defect density is below 97/10.
In the native dimension-six block, a normalized collective pair has
`[C†,C]=(12-D1)/42`, so its averaged commutator is at least 1/20.
The [exact collective moment calculation](../../model/particle_primary_energy_tightness/PROOF.md)
bounds the dimensionless centered second moment by `3141/224` uniformly.
Original transfer weight is at least3/100 in `[1/(100 beta_n),46/a_n]`,
and final excitation weight is at least1/50 in
`[1/(100 beta_n),46/a_n+100 epsilon_n]`. Complete finite clocks retain
at least1/50 in their widened reported interval. Here D1 denotes defects
in the Q6=1 block, not a physical quark count. At member one, the separate
charge-six scalar has original weight above1/20 in a band with upper
endpoint48 source units. The source's physical product anchor rules out
removing the inverse-spacing scale by an operator-norm bound on these
unsmeared probes; it supplies no lower bound on their selected-state
weighted residual.

The [selected-state density sum rule](../../model/particle_primary_energy_tightness/PROOF.md)
also gives `-a <H_hop>/V>0.03459` on every primary member, using the actual
prepared state. Its exact number-density double commutator determines the
signed first energy moment of a complete cosine/sine/failure instrument.
For its declared normalization `Lambda=30V`, that moment is at most
`a/(6V)`; cooling bounds the original downward absolute first moment by
`epsilon`. Consequently its original weight above any fixed positive
energy delta is at most `(a/(6V)+epsilon)/delta`, which tends to zero.
This is a limitation of this volume-normalized density apparatus, not of
the collective-pair channel or every source observable. It gives a concrete
test for the mass route: a candidate's state-weighted moments and original
response must support its probe normalization and energy window together.

The collective-pair bounds supply populated finite responses at every level.
A continuum field needs
weight in a fixed nonzero physical energy interval, which the moving
endpoints do not establish. A finite-capacity mass can be tested first:
resolve an energy branch in these channels, attach the corresponding
physical momentum and controlled dispersion, then apply the admitted
SI/GeV scale. The physical-mass table remains empty until that spectral
identification is proved.

## Finite spectral and spatial inputs

The [complete spatial Fourier detector](../../model/particle_primary_momentum/PROOF.md)
coherently addresses actual physical plaquette probes and records all eight
outcomes for a four-position aperture. Its unfiltered original probability
at each blank Fourier label is at least `2^(-b_n)/360`, using the selected
state's faithful floor. Filtered records have their own complete error
bound. The finite aperture has an explicit edge term; these labels are
spatial wavevectors without a proved conserved-particle momentum assignment.
Energy cooling alone does not establish momentum concentration.

The [channel calculation](../../model/particle_primary_channel_labels/PROOF.md)
proves conservation of the five native occupation labels with
`Q6=(-4,-3,1,2,6)` and block dimensions `(3,2,6,3,1)`. Full physical Gauss
invariance imposes `sum_q q N_q=0`. The plaquette and its filter preserve
every count. A distinct average over its 24 cube rotations is a cubic A1
observable with a positive original finite response. A mixed rotation-invariant
state is not a scalar vacuum, and A1 does not determine continuum spin or
a particle name.

The [actual-source spectral expansion](../../model/particle_primary_spectral_scale/PROOF.md)
uses `Hhat(lambda)=H0+lambda Vpert`, with only `lambda=1` the primary
action. The decoupled physical ground has energy `-90V` and gap six.
One hopping packet in native block r has squared norm `d_r` and cost
`Delta_r=6+(105/2)C_r`. Exact local coefficients give
`Ehat(lambda)=-90V-lambda²(B b1+P p1)+O(lambda³)`, where
`b1=117266225/340570896` and `p1=463/233377200`. Here B and P count
actual bonds and plaquettes. The remainder and interacting endpoint gap
are unproved; the reference gap is not a particle mass. This makes
endpoint spectral control a concrete calculation on the specified action.

## Full-strength finite energy bounds and joint records

The [local electric bound](../../model/particle_primary_endpoint_bounds/PROOF.md)
controls the actual interacting Hamiltonian, with every occupation and
quantum link retained. Its one-link Schur estimate is summed before the
physical Gauss restriction. The independent
[eleven-vector Ritz calculation](../../model/particle_primary_variational_readout/PROOF.md)
uses exact hopping and loop packets on the same source. It supplies a
variational upper bound at interaction coefficient one, rather than an
extrapolated perturbative coefficient. The
[occupation-sensitive electric bound](../../model/particle_primary_occupation_bound/PROOF.md)
retains the filled-reference block and its coupling to four beta defects.
For member zero it gives the full physical ground enclosure

`-727.498195 <= E0 <= -723.7678199` in source energy units, of width below3.731.

These are absolute ground energies. The
[exact parent Haar and CAR moments](../../model/particle_primary_residual_refinement/PROOF.md)
bound the same trial's full-space residual by5.231016 source units. All
omitted states are retained; degree-one link compression only decreases
the computed parent norm. Neither the ground enclosure nor this trial
residual isolates a particle branch. The reference vector is an operator
calculation device and does not replace the selected polynomial state.

## A detected finite neutral excitation band

The [concentrated pair calculation](../../model/particle_primary_branch_resolution/PROOF.md)
retains four fifths of the native charge-six onsite penalty in the actual
member-zero Hamiltonian. The full-occupation electric estimate bounds the
selected charge-six beta-defect count below1.840. A beta defect is an
occupied positive-beta mode or an empty negative-beta mode; it is not a
protected topological object or a particle species.

For the uniform average of sixteen local matched probes
`B=c_positive† c_negative`, the original success weight exceeds113/128.
The all-occupation source commutator obeys
`||[H,B]-6B||² <=3/2+sqrt(2)<102/35`. These exact local bounds place at
least3/20 original transfer weight in `[4,8]`. Cooling by at most1/16
then gives final ground-relative excitation weight at least7/80 in
`[4,9]`. A complete seventeen-outcome probe and two finite clocks report
`[63/16,129/16]` with original probability at least9/64. All failures,
energy records and the idle reference are retained. Additional compilation,
preparation or sampling errors need their own allowances; the large source
histogram is not executed. The broader [local response bound](../../model/particle_primary_response_moments/PROOF.md)
also retains at least1/2 in `[2,10]`.

The same original correlation satisfies `Re F(1/2)<-3/8` at a fixed source
time. A complete34-outcome Ramsey apparatus retains a negative signed
score with its stated finite evolution error. This is a resolved finite
time contrast, not proof of a single spectral line or a threshold onset.
Its experiment averages local addresses; the weight is not assigned to
an individual coherent A1u component.

The [physical A1u spectral construction](../../model/particle_primary_spectral_resolution/PROOF.md)
uses twelve source packets with exact rational coefficients. Its center is
approximately -717.753576 and its full-complement residual is below 4.947619.
The ground enclosure therefore places some actual A1u excitation in
`[1.06662597,14.69223697]`. The residual includes every omitted direction;
it is not the residual of the small compressed matrix alone. Neither this
enclosure nor a difference of variational energy uppers is a mass.

The trial has a complete physical preparation and spectral-window readout.
An identical memory starts in the declared filled-negative-beta,
constant-link anchor. An admitted Householder prepares the trial there;
a graded exchange places it on the source and retains the original source
state and its reference correlations in that memory. The source Hamiltonian,
energy reference and clock stay fixed, while the source marginal changes
by this apparatus intervention. The exchange preserves the combined source
and memory occupations. The extra memory allocation uses at most 696 qubits, and mean
source work lies in `[5.95174423,9.74461872]` source units. Ready memory and
controlled gates are explicit resources; no autonomous energy-conserving
controller or efficient gate compilation is inferred.

In that prepared state the spectral interval of radius 23/4 around the
trial center has weight above 137/529. One complete finite energy clock,
with all outcomes retained, reports its radius-widened interval with
probability above 1/4; even its lower endpoint is more than 0.2486 above
the authenticated ground upper. This is a detected excited A1u window on
the prepared state. It does not assign that weight to one isolated line.

On the original selected state the
[coherent scalar response](../../model/particle_primary_odd_residual/PROOF.md)
`K=T/(16 sqrt(2))` instead has original transfer weight above 5/1024 in
`[3/4,45/4]`, final excitation weight above 1/2048 in `[3/4,53/4]`, and
complete reported weight at least 39/8192 in `[23/32,361/32]`. Its mixed
initial and final states need not have the A1u operator's type. These
probabilities and the prepared trial's probabilities belong to distinct
experiments; overlap of their intervals does not identify a common line.

The [checkerboard response](../../model/particle_primary_staggered_response/PROOF.md)
uses `w_x=(-1)^(x1+x2+x3)` in the native charge-six pair operator. Opposite
weights on adjacent sites cancel the pair-creation part of its Wilson
commutator. The actual selected-state defect bound gives original transfer
weight above 17/1280 in `[3,9]` and final excitation weight above 23/5120 in
`[3,10]`. Complete two-clock records retain weight above 539/40960 in
`[95/32,289/32]`. Its A2g label is an operator representation; the finite
checkerboard pattern is not an attached conserved momentum or a fixed
physical wavevector under refinement. The Hermitian partner's signed first
moment exceeds 531/1600, without identifying a particle or a mass.

The [retained-penalty bound](../../model/particle_primary_branch_resolution/PROOF.md)
excludes isolation certification by optimizing the retained onsite/electric
fractions in its Schur-bound family against the authenticated
ground upper: its complement floor is always below -724 and that upper is
above -724; other resolvents, sharper ground uppers and the physical gap are
outside this bound-method limitation.

## Charged endpoints, momentum kernels and local fluctuations

The [charged-string apparatus](../../model/particle_primary_charged_string/PROOF.md)
uses an actual native charge-six character to join opposite endpoint
matter increments Q6=+6 and -6. A constant-link projection is part of
the intervention, not a replacement of the source state. All 48 successes
and the failure are retained. Original member-zero transfer weight
exceeds 3/10 in [39/4,135/4]. The instrument is gauge invariant and
fermion-parity even, so this result identifies a charged endpoint pair,
not a single charged fermion, continuum spin or measured electromagnetic
charge. The [coherent charge-transfer meter](../../model/particle_primary_charged_readout/PROOF.md)
uses two controlled couplings to record the admitted increments without
first measuring absolute charges. A radix-five counter at each cell
separates `{-6,0,+6}`; the meter preserves each original energy/probe
history, including failures and quantum reference correlations. The
charged-pair success event keeps original weight above 3/10, and reported
weight above 1/4 in `[77/8,271/8]`. Arbitrary charge changes are outside
this alphabet. A terminal absolute-charge measurement has its own
back-action, and an initial charge measurement has a disturbed-state law.

The [phase-scanned spatial apparatus](../../model/particle_primary_populated_momentum/PROOF.md)
has 32 coherent address/internal rows and 64 complete probe outcomes.
For `K(q)=sum exp(-iq.x) B_xalpha/(16 sqrt(2))`, its exact Wilson
creation residual has squared coefficient norm `16 sum_i cos²(q_i/2)`.
This gives populated original energy bands for all phases and a tighter
band in a neighborhood of the checkerboard cancellation. With
`kappa=(pi/16,pi/16,pi/16)`, eight scalar analyzer records have a joint
band-union weight above 0.06646; complete clocks retain more than 1/16.
The corresponding original final-excitation union has weight above 15/256
under the same initial-energy cut. Cube and inversion covariance are exact.
These finite aperture settings are not conserved momentum eigenvalues or
a measured dispersion curve, and their source occupations differ from
those of the prepared odd-sector experiment.

The [momentum-transfer construction](../../model/particle_primary_momentum_transfer/PROOF.md)
attaches complete coherent aperture records to the canonical joint spectral
measure of H and P in the retained-origin regular representation. A
distinguished original success has kernel
\(\operatorname{sinc}^2(S(E-\nu))\prod_i\operatorname{sinc}^2(Q(p_i-q_i))\).
The full finite analyzer keeps every spatial/time Fourier label, inner
outcome and local failure. Its quadrature, source-clock and boundary
errors are explicit. The source subsequence has no supplied effective
convergence rate, and the kernel's weight may vanish. Finite Fourier
labels are filters, rather than exact finite-box momentum projections.
A populated joint spectral ridge is the remaining input to a dispersion
calculation on this domain.

The [local Weyl construction](../../model/particle_primary_local_fluctuation/PROOF.md)
provides a distinct bounded probe with a fixed physical tent radius:
\(\exp(\pm i\pi a\sum_y f_R(y-X)(N_y-30))\), where \(R=1/a\).
It retains the origin and sign. The source-state hopping bound forces
mean injected energy above 0.06149 uniformly, while the original absolute
first moment is at most 14179/32. Its unfiltered probability responses
are therefore tight and have continuous positive-energy subsequential
limits of total mass one. Boundary-wrapped apparatus profiles are removed
on the same retained-origin menu with explicit full-reference error.
The limit may be entirely elastic: a first-moment floor does not supply
uniform integrability or a nonzero fixed positive-energy residue.

The [mixed-history construction](../../model/particle_primary_current_tail/PROOF.md)
places these unfiltered phase operations and the fixed-time-filtered
physical probes in one common positive-energy representation. Its exact
source double-commutator bound is 2025/2. For any fixed mixed word W with
k phase factors and filter times T_j, the energy-root drift is bounded by
`D_W=45k/2+sum_j 2 sqrt(2/T_j)`. With the selected mean excess epsilon,
`Tr rho W* (H-E0) W <= (sqrt(epsilon)+D_W)^2`. This controls original
absolute energy first moments and state-weighted time continuity on the same
source, state and clock.

One subsequence consequently supplies all finite mixed real-time words,
a strongly continuous canonical time group with nonnegative energy, and
complete fixed finite pointer histories. The proof retains the source
marginal/reference distinction and does not assert a trace-norm limit of
the growing reference. Spatial continuity for these unfiltered fields,
an effective subsequence rate and nonelastic survival need their own
estimates. A physical anchor gives an origin-averaged squared commutator norm lower
bound of 160/a, so a uniform absolute operator-norm estimate cannot replace the
selected energy-form method.

A sufficient **unproved** survival estimate uses that same source marginal
and complete origin/sign apparatus. Put \(K_n=H_n-E_{0,n}\), let \(\mu_n\)
be its original energy-transfer measure, and set \(m=(16/9)h_*>0\).
For one fixed \(E>m/4\), require, for every sufficiently large member,
\[
 R_n(E):=\frac1{2V_n}\sum_{X,\sigma}
 \operatorname{Tr}\!\left[\rho_n U_{X,\sigma}^{\dagger}
 K_n\mathbf1_{(E,\infty)}(K_n)U_{X,\sigma}\right]
 \le\frac m4.
\]
Ground subtraction is only in this estimate; the source clock and energy
records are unchanged, and the apparatus retains the idle reference.
Initial ground-relative energy is nonnegative, so
\(\int_{\omega>E}\omega\,d\mu_n\le R_n(E)\). Splitting the proved
signed mean floor across low energies, the fixed band and this tail gives
\[
 m\le\frac m4+E\mu_n([m/4,E])+R_n(E),
 \qquad \mu_n([m/4,E])\ge\frac{m}{2E}.
\]
A fixed triangular time filter with \(TE\le1\) retains at least half
this original band weight, without success normalization. This sufficient
criterion assumes no uniform second moment and is weaker than proving
one. The energy-form bound does not supply the required weighted tail;
the [elastic-limit control](../../model/particle_primary_current_tail/PROOF.md#6-the-physical-boundary)
shows why positive mean and probability tightness do not suffice.

The mass calculation therefore needs a resolved populated branch and
momentum/dispersion for the same apparatus domain. The explicit prepared
A1u window supplies a finite excited-state experiment; the coherent
checkerboard and charged-transfer records keep their original selected-state
laws. Their physical spin/current and species assignments are separate
conditions. For a continuum particle field the local phase route needs
nonelastic survival in its common history representation. None of these
conditions requires geometric gravity for a finite-capacity comparison.

For the Weyl route the missing estimate can be stated as a selected-state
current correlation. Write the dimensionless hopping current on bond b as
\(\widehat j_b=i(W_b-W_b^\dagger)\), with hopping term
\(T_b=(W_b+W_b^\dagger)/a\), and put
\(\delta_{Xb}=\pi a\,\Delta_b f_R\). The two signs of the same unitary
give energy changes \(A_X\pm B_X\), where
\[
 A_X=\sum_b(\cos\delta_{Xb}-1)T_b,\qquad
 B_X=a^{-1}\sum_b\sin\delta_{Xb}\,\widehat j_b.
\]
Their original second moment is exactly
\(V^{-1}\sum_X\langle A_X^2+B_X^2\rangle_{\rho_n}\).
The even term has uniform norm at most 2025/4. The required additional
control concerns the quadratic form
\[
 \frac1{a^2V}\sum_X v_X^T\Gamma_n v_X,\qquad
 (v_X)_b=\sin\delta_{Xb},\quad
 (\Gamma_n)_{bb'}=\operatorname{Re}
       \langle\widehat j_b\widehat j_{b'}\rangle_{\rho_n}.
\]
The selected occupation and hopping-mean bounds do not bound these
cross-bond correlations. Since \(\|v_X\|^2=O(a)\), even a uniform bound
on \(\|\Gamma_n\|\) gives only an \(O(1/a)\) estimate. A uniform second
moment needs stronger suppression on these smooth-gradient modes;
vacuum-vector Haar moments cannot replace the selected expectation.

The source-marginal tail criterion above avoids requiring this stronger
second-moment bound. A certified
[fixed-time correlation contrast](../../model/particle_positive_window/PROOF.md)
is a separate sufficient witness; a fixed filter alone cannot prevent the elastic-limit counterexample.

The [complete-history strength bound](../../model/particle_primary_weighted_survival/PROOF.md)
excludes a normalization shortcut. For the bounded local density
`Q_X=sum f_R(N-30)/(30R³)`, the exact source double commutator is at
most `9a⁴/80`. A complete coherent circuit of phase calls
`exp(i pi r_j Q_X tensor A_j)`, with `||A_j||<=1` and total strength
`S=sum |r_j|`, has final source excitation mean at most
`(sqrt(epsilon)+3a²S/4)²`. Its original probability above any fixed
positive energy therefore vanishes when `a²S -> 0`, even with coherent
filtering, all failure records and arbitrary retained references.
A nonzero event floor p at energy eta requires
`S >= 4 (sqrt(p eta)-sqrt(epsilon))_+/(3a²)`.

The original local Weyl probe has `r=30/a²`, outside that vanishing
strength domain. Its nonelastic survival is neither proved nor excluded
by this restriction. The selected weighted current or energy-tail estimate
is required on that actual stronger probe. Uncontrolled waiting-time
comparisons have a full-reference bound in this circuit family; controlled
clock queries keep their observable ground-energy apparatus phase.


The [complete joint spectrometer](../../model/particle_primary_joint_spectroscopy/PROOF.md)
records an initial source energy, all spatial Fourier outcomes, and a final
source energy, retaining `2 M J²` fine outcomes and the reference. The
source marginal is unchanged by its initial energy readout; the joint
purification is disturbed and that record is retained. Explicit bounds
cover aliasing, both clock tails, clock approximation and the conversion
from original transfers to final excitation. An unfiltered zero-wavevector
channel has a proved positive broad band. That Fourier label is not an
established rest frame, and the full histogram and a resolved branch are
uncomputed. The ground adapter recomputes its lower and upper bounds at
the aperture's actual member, so member-zero bounds are not assigned to
member-two records. Two energy registers and their precision costs are
additional resources beyond the compact Ramsey apparatus.

## One specified source, state, reference and clock

The [primary construction contract](../../model/particle_common_contract/PROOF.md)
resolves the specification requirement for one field family. For integer
`n>=0`, it selects an open Wilson box of side `2^(2n+1)`, spacing `2^-n`,
per-link native-word degree `n+1`, all sixty CAR modes per cell, and
`e0=105, kappa0=1/7`. Each cell owns three outgoing links; the local
Gauss constraints include their exterior endpoints with zero external
charge. The exact operator, basis order and energy zero are fixed by
[contract.py](../../model/particle_common_contract/contract.py).
These are declared architecture choices without measured mass inputs.

The prepared density is `rho=Q²/Tr(Q²)` on the source-only Gauss carrier,
where `Q=P_G[T_m(beta(U-H)/(2s))]^s`. Its temperature, least dyadic scale
and least certified Taylor order follow a deterministic plan with mean
energy excess at most `2^(-n-4)`. The exact canonical Householder
preparation retains a full conjugate reference with idle Hamiltonian
zero. Both factors are separately Gauss invariant. Source stationarity
does not imply joint-purification stationarity. A compiled approximate
preparation needs its own deterministic tape and complete error bound.

The clock is `U(t)=exp(-itH)`, with `hbar=1` in source units. On the
physical support, `K=log(Tr Q²)I-2 log Q` is its modular Hamiltonian;
the polynomial state does not have the exact Gibbs identity
`tau=-beta s`. That identity belongs to the separately stated exact
Gibbs controls. Reference counter-evolution is an apparatus gate,
not an alternative idle Hamiltonian. Declared coframe responses vary
`H` at fixed bare spacing and fixed input state.

Typed consequence receipts recompute the actual source bounds, Gauss
preparation, quadratic gauge/matter characteristic match, and
[regional support floor](../../model/particle_spatial_modular_response/PROOF.md).
They reject substitutions of the source, state, reference, clock or result.
The inventory keeps scalar, directional-rotor and periodic detector controls
in their own domains. It does not transfer their distributions to this
primary open-box state. The finite sequence is fully specified. Nontrivial
physical joint refinement beyond the bounded filtered subsequence,
relativistic particle sectors and proper geometric readouts are output
obligations.

## Physical preparation with retained boundary charges

The [finite boundary-state construction](../../model/particle_physical_boundary_state/PROOF.md)
resolves preparation including boundary charges on every member of the
selected program. Its full-support polynomial density is prepared by the
specified deterministic source/reference map. Zero external charge at the
outer box does not erase charges at a regional cut: inside and outside
representations pair to a singlet.

For the origin cell and its three owned links, exact integer hypercharge
projectors give the complete boundary-charge menu. At the first member,
it has 1,331 outcomes, retaining vacuum and all zero-weight outcomes.
An explicit neutral pointer unitary records those projectors while
preserving the full Gauss constraints and the outside reference. It is
not assumed to commute with H or to perform work-free measurement.
Actual native plaquette-character vectors and their conjugates exhibit
ten distinct nonzero cut-charge tuples; the blank gives an eleventh.
The selected state's faithful support proves positive probability for
each, bounded below by `2^(-582407133)` at member zero. This very small
bound proves positivity, not an abundance prediction. Their regional
Schmidt ranks `d²` are representation dimensions, not particle counts.

The executed one-link control retains all 119 native coefficient,
conjugate-coefficient and constant states. The spatial preparation and
complete instrument are finite matrix constructions; their enormous full
matrices and gate tapes are not claimed as executed. Full nonabelian
multiplicities are retained by the carrier but are not inferred from the
coarse hypercharge labels. Compatibility under simultaneous volume,
spacing and cutoff refinement is a separate output requirement.

## Local records on the selected physical source

The [primary locality theorem](../../model/particle_primary_locality/PROOF.md)
and its [typed consumer](../../model/particle_primary_locality/locality.py)
instantiate the finite Wilson cone directly on the selected program.
A product per-link projection preserves each bond and plaquette support
and cannot increase its norm. The electric coefficient `105/(2a)` acts
within the owning cell and is removed by the onsite interaction picture.
For spatial weight two the source has weighted interaction rate at most
`5232/(7a)`, independent of the per-link degree and exterior volume.
The proof uses the local ambient Hamiltonian and then restricts its
all-state bound to Gauss-invariant instruments and the physical state.
It does not treat a globally projected Hamiltonian as a local interaction.

The actual observer apparatus is the complete origin-cell boundary-charge
pointer. A remote interior plaquette phase is compared with identity,
using the same selected canonical purification, Hamiltonian duration and
idle reference. Every local outcome and the quantum reference are retained;
the remote field and any private remote record are traced out. At source
duration `1/10000`, selected levels 0-3 have conservative complete-output
half-trace bounds `2^-2, 2^-7, 2^-30, 2^-124` for their stated far-face
plaquettes. A separate level-six comparison keeps assigned separation one
and the same duration, giving `2^-21`. These are analytic bounds with
validated finite source and apparatus inputs; the full spatial histories
are not executed by the certificate.

The bound preserves exact zero-time non-signalling. Its conservative cone
is uniform in the declared spatial scaling, without identifying it with
a measured light speed or Lorentz action. This is a selected-state
consumer of unified node 3.

## Common filtered observables and complete records

The [primary observable-algebra construction](../../model/particle_primary_field_algebra/PROOF.md)
uses bounded physical charge and plaquette probes, fixed source-time
filters and all finite products and time translates. Their separable
sequence algebra, modulo norm-null sections, admits one subsequential
stationary state for every finite multitime word. Its canonical GNS time
group is strongly continuous and fixes the cyclic vector. The
[selected cooling bound](../../model/particle_joint_reconstruction/POSITIVE_ENERGY.md)
removes negative energy-transfer weight, proving a nonnegative generator
on the whole cyclic representation. The filtered source cone gives scoped
spatial commutation. No inter-level state embeddings or change of source
clock enter this construction.

Continuous functional calculus supplies complete finite effect instruments
on this algebra. Every fixed finite history and its finite pointer density
converge on the same subsequence. This does not prove convergence of the
joint state with the growing canonical reference or completeness of an
infinite fixed-charge menu. The [energy-form extension](../../model/particle_primary_current_tail/PROOF.md)
includes all fixed mixed histories of the local phase operations and these
filtered probes on one common continuous positive-energy time representation.
The algebra can represent its chosen probes as scalars or zero. Nonzero
finite-energy fields, Lorentz covariance,
surviving interactions and controlled physical renormalization are
additional conditions of nodes 4-5.

The [filtered-probe apparatus](../../model/particle_filtered_probe/PROOF.md)
constructs a finite coherent implementation of a triangularly filtered
boundary-charge parity. Exact positive quadrature weights, source-clock
Cayley echoes and a neutral selector retain every failure outcome and the
idle reference. Complete Ramsey records give the original filtered norm
and real or imaginary correlation, without success normalization. The
finite precision planner includes state, clock and compiler allowances;
its exact apparatus-service route and conditional source-code compilation
route are distinct. The enormous primary matrices and tape are not executed.
An authenticated Hermitian-contraction dilation extends this apparatus to
the centered charge and plaquette menu of the common algebra. Its LCU
success effect is `B*B`. A separate square-root instrument measures the
filtered charge effect `E=B` or plaquette effect `E=(I+B)/2`. If the effect
error is at most `delta`, its complete Stinespring error is at most
`sqrt(2 delta)`; compiler, preparation and interleaved clock errors then
compose over each fixed finite history, including arbitrary references.

The [positive-window criterion](../../model/particle_positive_window/PROOF.md)
uses a lower original contrast and explicit low-frequency, negative-transfer
and ultraviolet allowances. A uniform positive bound in a fixed band would
establish finite-energy survival. The finite primary witness below proves
positivity in a different, level-dependent experiment; it does not discharge
that uniform condition or identify a particle.

## A finite primary response and its nonuniform weight

The [primary survival theorem](../../model/particle_primary_survival/PROOF.md)
uses the actual normalized plaquette and the selected state's faithful
floor `rho>=2^-N_floor P_G`. A physical blank/loop matrix element gives
`||[H,A]b||²=4949/(20a²)` on the full sixty-mode Wilson source. The matter
commutator annihilates this witness, and the magnetic cancellation uses
the common final compression; compressed multiplication operators are not
assumed to commute generally. The blank is a matrix-element witness,
not a replacement preparation.

The proved channel has the following exact source labels. Gauge-matter
evolution can spread the filtered operator into matter factors, while
preserving these commutators.

| Property of the actual plaquette channel | Established identification |
| --- | --- |
| Unsmeared operator | Normalized real plaquette trace, acting on link factors |
| Local gauge transformations | Gauge invariant and neutral |
| Matter number and fermion parity | `Delta N_F=0`; even, for both the original operator and its time filter |
| Physical Lorentz spin and particle content | Unidentified; no single charged-lepton or other Standard Model species assignment |

Number preservation and even parity permit neutral composite or pair
responses; they do not select that content or a physical boson species.
The [survival proof](../../model/particle_primary_survival/PROOF.md) retains
the full matter carrier and the number-preserving Hamiltonian in this
classification.

For certified source width W, choose `T_n=t_n=1/W_n` and
`r=(23/24)^4`. The original contrast obeys
`D_n>=(11/24)r(4949/20)2^-N_floor/(a²W²)>0`. An explicit finite
quadrature/Cayley apparatus retains at least half that contrast. Polynomial
passivity also gives positive original band weight in
`[(8/a)2^-ceil(N_floor/2),W]`, for both the ideal and implemented filter.
At level 1, `T=t=7/31358976` and `N_floor=33469453062583`. The full spatial
matrix is not executed. The weight floor and lower endpoint tend to zero,
the upper endpoint diverges, and the times shrink; fixed-filter noncollapse
does not follow. `N_floor` is distinct from capacity `N=ln D`.

Empty and fully occupied Fock sectors reduce the actual Hamiltonian to its
full product-group gauge action. Their combined original selected weight
is at most `7/(627*2^(8n+7))`. Conditioning their spectra discards this cost
and does not produce surviving primary masses. Selecting a compact-simple
group factor requires its own source or a proved reducing map.

## Spatial action, excitation transport and declared units

The [regular spatial algebra](../../model/particle_primary_spatial_algebra/PROOF.md)
adds bounded L1 spatial smears on the same source. An invariant exterior
blank extends the actual finite Gauss state; a bulk locality comparison
controls its infinite-lattice time evolution. Rounded lattice shifts define
an exact point-norm continuous R3 action on the regular asymptotic quotient,
commuting with source time. Equality of limiting word moments does not
assert an isometric finite-carrier embedding. The fixed-origin limiting
state is not proved spatially invariant.

The [retained-origin observation domain](../../model/particle_primary_origin_momentum/PROOF.md)
adds an independent recorded origin, uniform over a bulk cube. The same
origin controls every step of each complete history. This block-diagonal
translation preserves products, whereas an arithmetic observable average
or a fresh origin per step does not. Its collar has vanishing volume
fraction and diverging physical width. The actual locality bound makes
all fixed-history boundary errors vanish, including origin-conditioned
quantum references, without changing the source clock or preparation.

A common subsequence gives a translation-invariant relative-observation
state. Its canonical generators strongly commute, with `H>=0` and
`H Omega=P Omega=0`. They can be trivial; no nonzero field or dispersion
is implied. The fixed-origin state is not replaced or proved invariant.
Only fixed pointer histories and relative-algebra moments take limits,
not the growing complete origin/reference density. The same recorded-origin
collective probe has positive finite original weight from member seven
by its selected occupation bound. Nonzero fixed spatial/time-smear response
and the convergence of finite Fourier readouts to these spectral momenta
require separate estimates.

The [mass-readout theorem](../../model/particle_primary_mass_readout/PROOF.md)
couples original transfer `omega=E_j-E_i` to final excitation `z=E_j-E0`
with mean discrepancy at most `C² epsilon_n` for probe norm at most C.
A transfer window `[l,u]` of original weight w gives final weight at least
`max(0,w-C² epsilon_n/r0)` in `[l,u+r0]`. The full-width primary band above
transfers without additional loss, including its finite implemented filter.
These are excitation weights; a rest-mass interval additionally requires
the same source's momentum, initial/vacuum momentum control and spectral
branch with controlled shell error.

An absolute source length or tick duration may be declared immediately for
comparison. For T seconds per source time unit,
`GeV/source-energy=hbar_GeV_s/T`; matching a declared source speed v to c
gives `metres/source-coordinate=cT/v`. Specifying c alone leaves the
absolute scale free. The declaration and any calibration retain their
ancestry. Neither the declaration nor a broad excitation band supplies a
missing mass pole, and the primary physical-mass table in GeV is empty.

## Source capacity and execution cost

The [capacity calculation](../../model/particle_primary_capacity/PROOF.md)
uses `N=ln D` in nats. The smallest ambient source has
`N_ambient,0=8(60 ln2+3 ln119)=447.4096105...` nats. Physical Gauss
constraints give `4V ln2<=N_phys<=N_ambient`; the lower bound counts actual
local singlet occupation vectors, without asserting an invariant hopping
subspace. Selected-state entropy is a different quantity.

The level-five ambient upper bound is approximately `1.302193221e12` nats
for more than eight billion cells. This is neither an exact physical
capacity nor the capacity of the complete experiment: the full conjugate
reference doubles the allocated ambient log dimension. Apparatus, program
storage, controller precision, record archives and workspace add resources.
A generic carrier of capacity `10^12` nats has
`exp(10^12)` complex amplitudes. These bounds do not claim tractable dense
simulation or a cosmological capacity identification. They do not make
`10^12` nats necessary or sufficient for every implementation.

The [compact response apparatus](../../model/particle_compact_mass_probe/PROOF.md)
directly dilates the identical finite Hermitian filter `B_tilde`. Two neutral
apparatus qubits give four complete Ramsey outcomes, retaining the same
original successful scores and band weights. Its complete failure
instrument differs from the selector construction. Source, full reference
and these qubits have ambient capacity `2N_ambient+2 ln2`; this is a
one-invocation carrier count for the declared direct service. Its algebraic
root and controlled operations have no supplied bound on computation,
controller precision, compilation, extra compiler flags, workspace or
repeated record archives. Autonomous repair implementation and efficiency
are unproved. A smaller detector register supplies no rest-mass
identification or total machine budget.

Refinement changes volume, spacing, link cutoff and cooling together.
`H_n=Hhat_n/a_n` gives `E_n=2^n Ehat_n`; at a fixed declared time unit,
a finite nonzero rest-energy limit needs the relevant dimensionless
excitation to scale like a_n and retain original weight. Neither a
mass-capacity law nor an isolated branch follows from the capacity count.

## Primary energy transport and the gravity branch

The [same-source modular-work comparison](../../model/particle_primary_gravity_premises/PROOF.md)
bounds source-marginal relative entropy for physical source unitaries
between `beta(1-delta) W` and `beta W`. For Gauss-preserving tensor-local
unitaries it bounds the regional marginal response by the same source work.
Its global modular clock has an explicit error relative to the source
clock. The complete retained-reference entropy is a different functional.
The [boundary-completed coframe](../../model/particle_primary_open_coframe/PROOF.md)
fixes a declared spatial strain coupling with `H(I)=H_primary` exactly.
Six strain derivatives have complete two-record mean detectors and
full-carrier bounds, using the unchanged flat input. The response can
disturb that input. This source coupling supplies no proper metric,
Lorentz stress, diamond kernel or area law. Gravity is separate from the
finite mass calculation.


Gravity node 10 is deferred from the active mass target. Its finite
source-energy and entropy results retain the following domain.

The [primary energy consumer](../../model/particle_primary_stress/PROOF.md)
partitions the unchanged Wilson action into local cell energies, including
owned electric links and allocated bond and plaquette energies. The exact
currents `J_xy=i[h_x,h_y]` give complete regional boundary continuity and
finite detector/flux-integration error bounds on the selected source clock.
The physical density is stationary; its retained idle-reference purification
need not be. Readouts used to estimate an undisturbed energy balance are
alternative experiments on copies of that same preparation.

The same consumer uses the actual supported regional marginal and the
existing finite entropy theorem for a complete identity/unitary mixture.
Its source-marginal first law retains physical boundary sectors. The
classical coin and quantum reference remain in the observer output;
the coin-plus-region entropy is `H2(s)+S(sigma)`, distinct from the smooth
marginal functional. Node 10 therefore has same-source finite energy and
marginal-entropy premises, without a proper metric, Lorentz stress tensor,
geometric diamond kernel or area response. The matched Einstein theorem
requires those objects and common shrinking-ball error estimates.

## Conditional chiral candidate and observable spectrum

The [joint reconstruction candidate](../../model/particle_joint_reconstruction/CANDIDATE.md)
and [recomputed action receipt](../../model/particle_joint_reconstruction/candidate_results.json)
make the additional hypotheses explicit. The primary source has sixty CAR
modes per cell, from fifteen internal rows and four Dirac components.
The candidate extension has three all-left families and two Weyl spin
components, giving ninety static CAR coordinates, together with a Higgs
doublet. Equality of internal matrices does not give a 60-to-90 state,
Hamiltonian or clock map. A chiral regulator, mirror control, physical
vacuum and interacting correlation limit are extension premises.

Its two independent action producers agree on the supplied charged and
bosonic quadratic coefficients. Running and tree coefficients are not
measured poles, widths or complete detector probabilities. The minimal
renormalizable arm has three massless neutrinos and cannot reproduce the
nonzero oscillation splittings. A separate calibrated Weinberg-operator
arm supplies a consistent neutrino tensor and retains its lightest-mass,
Majorana-phase, energy-window and running/matching assumptions. Its input
splittings remain calibrations. Neither arm proves the physical matching
or identifies the complete spectrum of the primary field. The pure gauge
gap question is node 9; a photon-containing continuum has no positive
global vacuum gap.

## Interacting four-dimensional Yang–Mills and the vacuum gap

The existing continuum results have distinct conclusions and domains.
The [original Weyl/contact field net](../../model/continuum_field_net/README.md)
supplies a positive-energy, local 3+1-dimensional continuum hierarchy on a
subsequence, without a proof of Lorentz covariance or surviving interactions.
The [centered Weyl limit](../../model/continuum_free_limit/README.md) is
proved free and massless. The [joint spin construction](../../model/spin_joint_field_limit/PROOF.md)
has nonzero charged response at arbitrarily small positive energies in its
specified scaling. Supplied massive Dirac symbols, a fixed-spacing interacting
ground, and a nonrelativistic moving pair do not supply an interacting
relativistic continuum with a positive vacuum gap.

The [Yang–Mills paper](../../extra/yang_mills_gap_clay_problem.pdf) states its
continuum and spectral-transfer premises explicitly. In particular,
[RepairGapChain](../../Lean/ObserverPatchHolography/RepairGapChain.lean)
assumes the repair gap and equality of the relevant spectra. Its conclusion
does not establish either premise. The
[variable-rate comparison](../../Lean/Screen/KogutSusskindFiberRateComparison.lean)
allows fiber-dependent repair rates with a positive floor; it avoids an
unnecessary constant-rate identity. Applying it requires the actual
interacting quotient measure's uniform variance bound and a transfer to
its physical clock. A product-measure calculation supplies neither.

A theory containing a massless photon cannot have a positive global vacuum
gap. The [Clay target](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)
is pure Yang–Mills for every compact simple gauge group, with a nontrivial
four-dimensional continuum theory and the stated axiomatic and short-distance
properties. Selecting a Lie-algebra factor of the native representation does
not select an invariant subspace of the coupled field Hamiltonian. A pure-gauge
comparison needs a separately specified program or a proved decoupling map.
It does not replace the primary observer-field program silently.

The [local-correlation gap criterion](../../model/particle_yang_mills_gap_contract/PROOF.md)
provides a direct alternative to assumed spectral equality. If the same
reconstructed positive-energy time group has a dense family of centered
local vectors and their Euclidean correlations obey one exponential decay
rate, its spectrum has the corresponding gap. It is enough to prove the
correlation estimate with vanishing complete errors at each fixed time,
take the regulator limit, and then take long time. A retained finite-energy
vector with nonzero variance excludes an empty or infinite-energy limit.
These correlation bounds are not established for the primary field merely
by cooling it. A stationary approximate ground can contain negative-energy
transfers; the proof includes an explicit two-level counterexample.

These Yang-Mills conditions instantiate nodes 1-5 for a pure-gauge source
and node 9 of the unified graph. The prepared coupled matter source does
not supply that source or its vacuum by sharing a gauge Lie-algebra factor.
The exact decay and finite-energy noncollapse conditions are recorded in
[the gap contract](../../model/particle_yang_mills_gap_contract/PROOF.md);
they refine this branch rather than create another observer-preparation
checklist. Finite physical preparation and joint physical refinement have
separate acceptance conditions.
The selected action has `H_n=(1/a_n) Hhat_n` at fixed dimensionless couplings.
For a proof using uniform finite-sector gaps and a finite-energy
noncollapse witness, the dimensionless thresholds must be of order `a_n`.
A fixed positive lattice gap would instead diverge in that clock. A slower
prelimit mode whose original local weight vanishes need not survive in the
limiting observable sector. Changing the time unit changes every propagation speed and
requires a fresh common-clock comparison.

## Low-energy states, conjugate records and common propagation

Three constructions connect the field action to reusable finite consumers.
Their [source calculations](../../model/particle_spatial_ground_state/PROOF.md),
[spectral records](../../model/particle_ground_spectral_records/PROOF.md) and
[gauge characteristics](../../model/particle_gauge_characteristics/PROOF.md)
use no measured mass or coupling. The full spatial family and its
208-state directional control have distinct Hilbert spaces; the control
does not become a spatial vacuum by sharing the same operator formulas.

**A low-energy family without a gap oracle.** On each finite spatial box,
retain all sixty CAR modes per cell and a gauge-invariant tensor product
of native coordinate-word spaces, with degree at most `ell` on each link.
Compress the declared parent Hamiltonian into this space. This product
cutoff differs from the total-degree cutoff used by the separate finite
readout certificates. For `V` cells its dimension obeys
`D <= 2^(60V) 451^(3V ell)`. Explicit electric, magnetic and kinetic norm
bounds give a spectral interval `[L,U]` without diagonalizing this matrix.
For the same Hamiltonian's Gibbs target,

`Tr(rho_beta H) - E0 <= log(D)/beta`.

A positive, scaled Taylor polynomial in `U-H` prepares a stationary
approximation through its canonical purification. The joint unitary
retains a conjugate gauge reference; there is no selected success branch.
If its half-trace error is at most `eta`, its mean ground-energy excess
is at most `log(D)/beta + (U-L) eta`. Choosing temperature and polynomial
order from this formula reaches any positive finite energy tolerance.
This is an explicit finite construction with very large resources; no
efficient many-body ground algorithm or autonomous cooling law is inferred.
The [physical Gauss preparation](../../model/particle_gauss_ground_state/PROOF.md)
restricts this construction to `P_G`, the product of commuting local Haar
projectors. Its source and conjugate reference are separately invariant;
invariance of an unrestricted mixed density alone would not establish
this source-only support. A nonconstant Wilson-loop singlet has nonzero
magnetic coupling to the blank, so the physical interaction is nontrivial.

For a stationary state with excess at most `epsilon` above the restricted
ground, every bounded probe preserving the physical carrier satisfies
`omega(A* [H,A]) >= -||A||^2 epsilon`. A charged apparatus needs its full
compensated domain and the corresponding accessible-sector ground energy.
This supplies an approximate ground-state premise directly from a preparation
budget.
For the finite-range Wilson member, local interaction bounds and a
diagonal compactness argument give a stationary algebraic ground-state
subsequence at fixed spacing and fixed per-link cutoff. Spatial averaging
can make it translation invariant. The source-only Gauss version supplies
a ground state on the gauge-invariant observable algebra, with all CAR
occupations retained; it gives no computable selection, uniqueness or
spacing-refinement rate.
Small total energy excess alone does not prove convergence as spacing
tends to zero. The separate [Weyl/contact ground construction](../../model/thermodynamic_ground/README.md)
supplies a reusable compactness argument on its own source, without
transferring that source's state to the gauge-coordinate field.

The [fixed-volume link-cutoff theorem](../../model/particle_link_cutoff_limit/PROOF.md)
removes the electric cutoff on the same parent Wilson source. Native word
spaces form an electric form core; ground energies converge, and cooled
canonical preparations have trace-norm subsequences. Complete bounded
detector histories converge with the reference and every failure label
retained. This is a fixed-spacing, fixed-volume result. Combining it with
a separate fixed-cutoff volume subsequence does not prove the primary
contract's simultaneous refinement.

The spatial restriction uses the actual CAR ordering and outgoing-link
partition. For the primary Gauss state,
`rho_R >= 2^(-N) Tr_out(P_G)`, where `N` is computed from its exact
preparation and source bounds. The first member gives `N=582407133`.
This proves faithfulness on the supported regional algebra, with boundary
charge center weights retained; it does not assert full ambient rank or
identify the region with a proper geometric ball. The
[regional entropy construction](../../model/particle_spatial_modular_response/PROOF.md)
proves the finite first law with a relative-entropy remainder on that
support. Its four exact rational spinor/flux controls execute the complete
apparatus on separate full208 and physical78 states. Those internal cuts
are not spatial regions, and their modular energy differs from bare local
energy. A geometric modular kernel and an area response are additional
premises of the Einstein consumer.

**Original charged records and energy signs.** The unchanged 208-state
charge-six quantum-flux control admits complete compensated CAR
addition/removal probes. Charge conjugation exchanges opposite charge
sectors while preserving this Hamiltonian. Rational spectral enclosures
give its ground and charged thresholds. The
[complete charged-window apparatus](../../model/particle_charged_residue/PROOF.md)
certifies a nonzero original lowest-window weight near `0.498558` for
each charge sign, and excitation energy near `1.107713183`. Its rational
ground trial, forty complete detector records and bounded source-clock
Fourier plan belong to that finite control. These are energies in its
declared clock, with no GeV conversion or physical species assignment.
For a complete contraction instrument, the total original weight of
transfers `Delta E <= -delta` is at most `epsilon/delta` for every
positive threshold `delta`. The
[Lean theorem](../../Lean/QFT/GroundRecordTail.lean)
uses nonnegative raw spectral weights whose row sums do not exceed the
preparation weights. It needs no spectral gap and does not condition on
accepted detector events. This is the finite positive-energy input;
a physical spectrum condition also needs convergence and translations
on the interacting spatial state.

**A common quadratic gauge and matter cone.** The native fifteen-state
representation has trace indices `2` in each weak/color direction and
`10/3` for hypercharge `Y=q6/6`. The electric Casimir uses their inverse
trace metric. Expanding the actual Wilson plaquette action therefore
gives the same speed squared `e0*kappa0/15` in all twelve directions.
The lattice curl has two transverse modes and one longitudinal null at
nonzero momentum. The Wilson matter symbol has principal speed one.
The declared relation `e0*kappa0=15` matches their quadratic cones on
one clock. It does not fit measured data or force that architecture:
an unmatched member retains its calculated speed difference. The
commonly used choice `e0=1,kappa0=1/7` has gauge speed squared `1/105`.
Charge-coordinate rescaling leaves these speeds unchanged. Interacting
Lorentz covariance, physical photon/gluon states, running couplings and
chiral matter are separate consequences to establish on the prepared
spatial family.

<a id="required-outputs-of-a-declared-construction"></a>
<a id="six-output-obligations-evidence-and-exact-boundaries"></a>

## Construction outputs and comparison evidence

The [unified graph](../../model/particle_joint_reconstruction/DEPENDENCIES.md)
contains the physical acceptance conditions and source domains. Carrier,
controller, scheduler, regulator, preparation and apparatus may be declared
architecture choices. Settings selected using comparison data retain that
calibration ancestry. Self-reference is a consistency condition on actual
records; it supplies no substitute for their calculated laws, and requires
no proof that the architecture is unique.

The primary construction supplies one physical state, complete boundary
charge readout and finite locality comparison. Separate scalar and
rotor controls supply their own spectral, thermal and detector laws.
Their use as supporting calculations requires an explicit state, operation
and clock map. A same-clock coefficient inverse is not automatically a
physical particle-mass inverse. Likewise, whole-field crystal momentum,
proper time and single-particle momentum have different operator domains.
The detailed examples below record those distinctions within nodes 3-8,
rather than proposing a parallel construction checklist.

The [Einstein consumer](../../model/matched_einstein_limit/PROOF.md) is a
separate geometric target: proper stress and clocks, regional modular/area
response and shrinking-ball errors on one state. Its tensor decoder needs
ten timelike readouts; no repeated proof of its tensor algebra discharges
the physical realization of those readouts.

For the declared attachment `H_Q=(3I-Q)V/2`, the scalar vacuum at one tick
`t=pi/(6 sqrt(56))` has occupied probability `1/4` in sector `Q=+1` and
`3/4` in sector `Q=-1`. The geometric birth experiment retains no-pair/pair
weights `9/25,16/25`; its four original birth/occupation weights are
`27/100,9/100,16/100,48/100`. Native memory and face parity are alternative
controllers, not identified registers. Protected transport commutes with
the coupled Hamiltonian. Birth preserves the defect/work cost but changes
the coupled field's energy distribution; a separate two-time energy audit
records work `0,+sqrt(56),-sqrt(56)` with weights `9/25,8/25,8/25`.
That audit disturbs the scalar preparation and is a different experiment
from the occupation protocol.

The [fixed-spacing causal bound](../../model/particle_field_causality/PROOF.md)
is derived from local interaction norms, independent of
the unobserved exterior volume. Its declared control bounds a remote
intervention 32 cells away by total variation below `1e-7` at duration
`1/64`. It is an analytic bound on complete records, not an executed
large-Fock simulation or an exact relativistic light cone. Total-degree
representation compression has a nonzero distant commutator: apply the
locality estimate to the parent field and add cutoff errors separately.

Energy conservation follows from `J_xy=i[h_x,h_y]` for the actual partition
`H=sum_x h_x`; coarse regional fluxes cancel exactly. An energy response
to a supplied action coefficient does not by itself define a metric.
The declared coframe construction below supplies action-derived spatial
strain. The gravity consumer requires its identification with proper
stress, together with curvature, entropy and area on one state and
proper-clock family. These distinctions
prevent finite energy or modular identities from silently becoming the
Einstein equations.

### Prepared states and actual regional response

The [full scalar ground certificate](../../model/particle_scalar_ground/PROOF.md)
uses `H=W(diag(1,2,3))+W†+N_H`, with fermion offset zero and all
`3*2^18=786432` states. A signed three-mode decomposition reduces the
complete Gram comparison to blocks of dimension at most 729. Exact
integer Perron bounds and Sturm root counts isolate a unique gauge/spin
singlet with `E0≈-7.94630530665`, `Q*=6` and full isolation gap greater
than `0.04690056`. Its scalar occupation transition is instead
`Delta≈16.89261061331`. The complete two-time laws retain the same action
and clock. A six-fermion singlet seed and controlled preparation have
certified infidelity below `2.6e-11`; autonomous evolution from empty
Q*=0 cannot reach this different conserved sector. Preparation work is
assigned to the controller. The general block theorem transports by family
SVD, but these numerical values certify only the stated diagonal action.

The [full-ground charged instrument](../../model/particle_ground_charged/PROOF.md)
acts on that Q*=6 singlet, with compensated addition into Q*=7 and removal
into Q*=5. It calculates twelve complete family/multiplet/direction
channels. Literal Hamiltonian closures have dimensions 53 through 240;
their orthogonal complements and all spectral groups are retained.
Integer Gram/residual checks certify complete spectral projectors and
their original weights, with the ground-preparation error included.
For example, a first-family left-removal line lies in
`[0.0710179490,0.0710358526]`, with single-mode weight in
`[0.6615763742,0.6616186281]`. Its complete-instrument weight includes
the factor `4/18`. A corresponding addition line is near `0.48728`.
These finite lines are not mass shells or species labels.

For the same prepared reference, adding `epsilon Q*` translates addition
energies by `+epsilon` and removal energies by `-epsilon`, while preserving
their weights and conditional variances. Certified variance intervals
exclude equal translated addition/removal shapes in every multiplet.
The Q*=6 state remains the unique ground for the sufficient window
`|epsilon|<ground_gap/12`; outside it a global-ground comparison is
required. This obstruction does not exclude a different vacuum or assert
physical CPT violation. The tested diagonal action and the coherent source
matrix require separate ground and charged calculations.

The [gauge-word state](../../model/particle_field_preparation/PROOF.md)
is the normalized gauge average of
`product_e ((1+Re Tr_M15 U_e/15)/2)^L` on a periodic cubic carrier.
It belongs to a finite word code, has explicit normalization and electric
energy bounds, and can be prepared by a complete unitary. Its conservative
one-cell 20-bit gap planner uses `L=67108864`; the resulting certified
duration is extremely short and has no volume-uniform lower bound.
An executed L=1 character preparation supplies 1728 original electric
records, including 397 zero complements, and mean electric energy
`99/9100`. That small seed does not satisfy a useful gap guarantee.
The unbounded coordinate parent is an analytic comparison; executable
preparation and evolution use finite codes. Periodic preparation and
open-box locality comparisons retain their distinct boundary and
initial-marginal hypotheses.

The [energy-confinement construction](../../model/particle_field_preparation/REFINEMENT.md)
reuses the plaquette admissibility theorem with
`S=sum_p(1-Re Tr_M15 U_p/15)`. A Wilson gap below `1/2` implies
`S>1/48000`. For `H=H_E+K_m+kappa S` up to a constant, positive electric
energy and conserved total energy give the all-time bound
`Prob(Bad) <= 48000 [<S>initial+(<H_E>initial+2||K_m||)/kappa]`.
It also holds for the compressed bad-gap effect on every finite word
cutoff containing the preparation; this effect is not the gap of the
compressed Wilson matrix. The explicit finite-time gap preparation and the
all-time configuration are different choices. The latter is a one-cell
periodic construction with
`L=2^86`, word degree `3*2^86`, `m=128` and
`kappa≈5.4973e69`, giving all-time bad weight below `6.43e-13` and
polynomial/exact-sign final-record TV below `0.000208` on the periodic
parent at time one. Finite-cutoff, preparation, compiler and volume errors
are additional.
Those values have no physical gauge-coupling identification. Enlarging
kappa can freeze an interaction; uniform interacting refinement requires
its own common-state correlation and parameter estimates.
The Wilson-kernel gap controls a regulator. It is not a physical mass
gap: the flat free source has Wilson gap at least one while its overlap
field has a zero at momentum zero. The one-cell, order-128, time-one
configuration is distinct from the order-zero, side-21 experiment below.

The [prepared finite field experiment](../../model/particle_prepared_field_readout/PROOF.md)
composes preparation, locality, representation cutoff and a complete
occupation apparatus on one specified spatial program. It uses `m=0`,
`a=1`, periodic side 21, an interior radius-five box and the actual
whole-box marginal of its gauge-word preparation. All sixty origin modes
are filled and the neighboring cells are empty; total fermion number 60
is conserved without filtering other source occupations. The same large
plaquette coefficient is used in every comparison. Boundary charge is
retained by a conjugate preparation reference. A degree cutoff
`K=K0+12` gives a finite Hamiltonian and one repeated-unitary clock.

At its declared tick `T≈1.23896e-152`, all 61 origin-occupation records
have rational bounds. The event with occupation 59 has a positive finite
probability of order `1e-303`, from the actual hopping norm
`2880/289`. The composed periodic exact-sign to finite-record total
variation error is at most `1.05065e-149`, before implementation errors.
That absolute error exceeds the rare-event probability, so the exact-sign
parent interval has lower bound zero. The finite event is proved; its
positivity in the parent is not resolved by this comparison. The enormous
state and propagator are specified by finite words and bounds rather than
materialized. This small-tick construction supplies neither a uniform
interacting refinement nor an identified physical sector.

The obstacle to a fixed-duration volume limit is explicit in these same
bounds. The [preparation estimate](../../model/particle_field_preparation/REFINEMENT.md#all-time-bound-from-the-actual-conserved-hamiltonian)
chooses `kappa=96000(E_el+2B)/eta`, where
`E_el=(29/10) E L²`, `E=3s³` and `s` is the periodic side. For fixed
target `eta` and `L>=1`, this coefficient grows at least with the volume.
The [composed locality estimate](../../model/particle_prepared_field_readout/PROOF.md#3-periodic-locality-open-dynamics-and-representation-cutoff)
has `L_rate=2880/17+192 kappa` at `m=0` and contains
`2^(-R+4)[exp(2 L_rate T)-1]`. Along `s` proportional to `R` at fixed
positive duration `T`, its factor `exp(384 kappa T)` defeats the spatial
decay. This estimate therefore supplies no fixed-action thermodynamic
or refinement control. It is not a no-go theorem for the dynamics.
This projected-state estimate requires a different local-energy or
admissibility argument, or controlled isolated low modes, together with
local correlation bounds at nonvanishing physical duration that preserve
the interactions. Enlarging the box or repeating the tiny tick alone
does not provide that input. The mixed preparation below supplies a
separate family with fixed local parameters.

The [conditional local-density bound](../../model/particle_field_preparation/REFINEMENT.md#conditional-local-roughness-from-energy-density)
is one precise alternative. For a translation-invariant source marginal,
unchanged evolution, finite initial electric and plaquette densities
`e_E,s_0`, and `||K_m||<=b V`, each plaquette's roughness probability is
bounded by `min(1,48000[s_0+(e_E+2b)/kappa])` at every time.
There is no volume factor, provided those densities have uniform bounds
on the actual state family. The finite-code effect is the compression
of the full roughness detector, with its complement retained. The
filled-origin experiment is not translation invariant and does not
instantiate this lemma. A translation average is a different preparation.
The distinction between local spectral control and complete dynamical
response is specified by the
[operator interface](../../model/particle_field_preparation/REFINEMENT.md#the-missing-local-operator-estimate).

The [referenced mixed gauge state](../../model/particle_local_gauge_state/PROOF.md)
realizes those uniform initial-density hypotheses. It independently
normalizes each finite link polynomial and twirls the full product
density over the vertex gauge group, retaining the maximally mixed
sixty-mode Fock state at every site. This is a unitary mixture, without
amplitude projection or smooth-outcome conditioning. A conjugate charged
reference makes the joint preparation gauge invariant; the source density
itself includes nonsinglet sectors and is not a vacuum. The reference
Hamiltonian is zero during observation.

At spacing one and polynomial order 128, one choice uses `L=2^64` on
every link and `kappa≈2.98009318237e50` for every periodic volume.
The same energy argument bounds each specified rough-plaquette
probability at every model time by
`11567/17179869184≈6.73288013786e-7`. Finite word compression preserves
this bound for the compressed full roughness effect and its complement.
The exact one-open-link control retains 119 positive-support states,
including total nontrivial charge probability `1/91`; it is not substituted
for the correlated periodic site twirl. These results do not change the
distinct filled-origin experiment.

The [local Wilson spectral bound](../../model/particle_local_wilson_density/PROOF.md)
uses the full native kernel and a positive square decomposition. For
`Pi_low=1_(|H_W|<1/2)` and a set X of cells, it gives
`<Tr(P_X Pi_low P_X)> <=400 |X| sqrt(6 <S>/V)`.
On the declared state this is at most `2 sqrt(5 eta)|X|` at every time,
with target `eta=2^-20`, independently of volume. The normalized trace defines a complete binary
finite-code effect. The result bounds the local one-body
Hilbert–Schmidt error between the sign and polynomial kernels; it neither
proves a global gap nor spatial decay of each low mode.

The same mixed family has
[exactly matching local marginals](../../model/particle_local_history_family/PROOF.md).
Fixed-parameter polynomial locality therefore gives a thermodynamic
limit for every specified finite local instrument history at any fixed
duration. An open-box word code approximates its complete law. The
order-zero, time-one example allocates `2^-24` to each of the volume
and cutoff errors, totaling `2^-23` before implementation errors.
Its radius and code are enormous and are not executed as a simulation.
All 61 occupation labels and intervening records are retained. The exact
initial curvature of the two-time squared occupation increment is
`2880/289`; this is a nonzero dynamical diagnostic, without a finite-time
particle or scattering interpretation.

Full-volume canonical purifiers are not identified across volumes.
The default local-record comparison traces them. A fixed independent
reference, or the same local preparation channel acting on it, gives
compatible local reference outputs; arbitrary extensions require that
compatibility explicitly. Measurements can break translation symmetry,
so uninterrupted-evolution roughness bounds are not reused after them.

### Regulator choice and complete local quantum response

The [native sign-sensitivity calculation](../../model/particle_sign_response_boundary/PROOF.md)
tests an inference that the local spectral bound does not justify. A
legitimate sixty-column weak-link family has squared Wilson eigenvalues
`delta^2,1,3+2delta`, with ranks `16,28,16`. Its unit-electric derivative
of the sign has squared Hilbert–Schmidt norm `(1-delta)^2/delta^2`.
Giving that branch original weight delta makes its low-mode trace tend
to zero while this derivative cost diverges. Normalizable packet controls
retain their preparation flags; no uniform electric-energy bound is
asserted. This refutes a shortcut from unweighted low-mode density to
inverse-gap electric response. It is not a failure theorem for the
configured mixed preparation or for its fixed-time interacting records.

The [smooth field construction](../../model/particle_smooth_sign_field/PROOF.md)
chooses a different, explicit action:
`h_delta=[beta+sqrt(1+delta^2) H_W/sqrt(H_W^2+delta^2)]/a`, with fixed
positive delta. Its uniformly convergent finite-path expansion gives an
exponentially summable interaction on the full quantum electric/link and
fermion carrier. Local Duhamel bounds compare every specified complete
instrument history with polynomial approximants on the same input and
clock, uniformly in unobserved volume and untouched reference dimension.
No Wilson gap, low-mode estimate or link-derivative assumption enters
that comparison. Compatible local preparations supply its thermodynamic
record limit; full-volume purifiers are not identified.

For delta=1 and duration one, the certificate assigns `2^-20` each to
smooth-to-polynomial, finite-volume and representation-cutoff errors,
giving total `3/2^20` before preparation, apparatus and numerical errors.
It specifies finite resources without executing the enormous code.
The mixed input and fixed plaquette coefficient are retained; energy
estimates use this action's recomputed kinetic allowance, rather than
the Legendre action's value. Its free symbol has one zero at momentum
zero, all other corners at energy at least `2/a`, and leading Dirac
velocity one. Both chiralities and every internal component are retained.

The delta-to-zero comparison is needed only to reproduce the
discontinuous-sign architecture. It is not a prerequisite for the
permitted positive-delta construction. The finite-range Wilson action
`(beta+H_W)/a` can also be chosen directly; it is the uniform
infinite-smoothing endpoint on the bounded Wilson spectrum. Neither
choice inherits exact-sign chiral-projector identities. The physical
influence estimate below controls a declared spacing refinement; a
convergent physical state and interacting records, chiral matter and mass
shells require additional same-action estimates.

### Refinement-uniform influence and an interacting thermal state

The [physical-distance construction](../../model/particle_physical_causality/PROOF.md)
uses the smooth action at fixed delta, or a separately declared finite-range
Wilson action, with `0<=kappa(a)<=kappa0/a`. Onsite electric coefficients
need no common upper bound. The source interaction norm obeys
`W_a<=W_star/a`. For `0<a<=1` and `lambda>1`, complete local records
affected by a remote intervention obey

```
half_trace <= min(1,C a^(-q) exp(-g/a)),
g=mu d-2 W_star T>0,    mu=(lambda-1)/lambda,
C=2(2K+1)(2R+3)^3,     q=p+3.
```

The assigned physical separation is `d<=a dist_graph(X,Y)`; detector
supports lie in a physical cube of half-width R and use at most `K/a^p`
complete intermediate interventions. All outcomes and untouched references
are retained. The exact planner certifies every smaller real spacing,
including support growth. For the finite-range Wilson choice, native bond
norm `30/a` gives `W_star=360 lambda+48 kappa0 lambda^2`. With lambda=2,
kappa0=1/100, T=1/1000, d=4, R=1/10 and two interventions, the complete
record error is at most `2^-41` for every `a<=1/128`. This is a conservative
causal bound in declared length and clock units, not a measured light
speed, a Lorentz transformation or existence of an interacting continuum.
Finite-box, representation and implementation errors add separately on
both compared experiments. Choosing total implementation error at most a
on each side preserves vanishing influence outside the cone.

A state suitable for physical interpretation is a separate requirement.
The [same-H variational witness](../../model/particle_field_state_response/PROOF.md#5-a-concrete-boundary-of-the-spatial-mixed-preparation)
lowers the degree-zero spatial mixed preparation's mean energy by
`750 V/(17a)`: fill both negative-beta spin rows for all fifteen internal
components at every cell, retaining the identical gauge density and
electric/plaquette terms. This gauge-singlet Slater factor has no
intersite coherence. For periodic side at least two, the comparison proves
that the original tracial Fock preparation is not a finite ground state.
The coefficient belongs to the degree-zero Legendre action, not the smooth
action. It neither establishes that the alternative is the interacting
ground nor excludes every continuum interpretation of the mixed family.

The [finite interacting state construction](../../model/particle_field_state_response/PROOF.md)
uses the unchanged charge-six quantum-flux Hamiltonian on all
`16*13=208` states. Its target is `rho_beta=exp(-beta H)/Z` at declared
beta=1/64. A positive rational exponential polynomial gives an actual
faithful `rho_M`, a deterministic canonical source/reference preparation
and a certified half-trace error to that target. There is no postselection
or autonomous cooling assertion. All 169 original two-time flux labels
are retained under an exactly unitary rational approximation to the same
clock. At duration1/64 and four Cayley steps, the evolution allowance is
`1331/6291456`, with preparation error added. The target obeys source
KMS, passivity, entropy/energy identification and thermal response. The
actual polynomial density is stationary and passive, but is generally
not exact Gibbs at beta; its logarithm is not replaced by `beta H`.
Its canonical purification is not stationary under `H tensor I`, and
arbitrary joint operations using that reference are not source-passive.

### Declared spatial geometry and same-state pressure

The [coframe response](../../model/particle_geometric_response/PROOF.md)
declares a constant spatial coframe V, metric `G=V^T V` and volume
`v=det V` on its periodic parent carrier. Its shifts are unitary. Direct
substitution of the primary open shifts instead gives
`H_primary-H_naive=dGamma(beta_D sum_i Pi_i)/(2a)`, where
`Pi_i=I-T_i†T_i`; the physical operator norm is `45 L²/a`.
The [boundary calculation](../../model/particle_primary_gravity_premises/PROOF.md)
therefore requires an explicit boundary-completed coupling before any
same-action transfer to the open primary. Dirac differences couple through
`V^(-T)`, Wilson differences
through `G^(-1)`, and electric and based-plaquette Gram forms through
`G/v`. The magnetic Gram form includes its subtraction
`Q_ij=Re Tr[(I-R(P_i))^dagger(I-R(P_j))]/30-delta_ij`, so at V=I it
retains exactly `Q_ii=-Re Tr R(P_i)/15` and the specified magnetic energy
zero. Identity plaquettes have `Q=-I`; at nonzero kappa their scalar
energy and coframe derivative cannot be discarded. All six strain derivatives
are gauge preserving. Exact noncommuting-link checks retain all sixty
native spin/internal columns at kappa=1/7; this is a coordinate-symbol
control, not a normalizable gauge-field state. The full Fock extension is
algebraic. Its magnetic scalar multiplies the Fock identity, rather than
being second-quantized as part of the matter operator.

At fixed bare regulator length, isotropic `V(s)=exp(s)I` gives
`K_iso=-H_E-H_D-H_B-2H_W`. The extra Wilson term cannot be dropped.
Constant frame rotations have their own covariance identity, without
proving Lorentz boosts.

The [coframe locality estimate](../../model/particle_geometric_response/PROOF.md)
joins this ordinary-Wilson action to the physical-distance construction.
Every cell owns its outgoing quantum links, so the positive electric
quadratic form, including cross terms, remains onsite. The hopping and
magnetic terms give
`W_a(V)<=W_K/a` uniformly for V in a fixed compact subset K of
positive-determinant coframes and `kappa(a)<=kappa0/a`. Bounds on V and
its inverse also control physical-region support growth. The assigned
coframe distance satisfies `d_V<=a b(V) D`, where
`b(V)=max_i ||V e_i||` and D is graph distance. On periodic boxes d_V is
the torus quotient distance, minimized over periodic lifts; an unwrapped
coordinate separation across a periodic seam cannot replace it. This
extends complete-record causal control to the declared geometry family,
without identifying that dictionary as proper geometry or proving a
limiting interacting state.

The same 208-state quantum-flux Hamiltonian and
beta=1/64 target above supply a quantum directional control at the separately
declared kappa=0. Its exact action has no magnetic term. This is an algebraic
reduction with transverse terms removed, not a normalizable reducing
subspace of the full three-dimensional link theory or the kappa=1/7 check.

Complete effects `(I+-K/28)/2` read either its longitudinal or isotropic
strain mean, with all occupations and the conjugate reference retained.
Their values on `rho_M` are exact rational records; the target Gibbs mean
has the separately certified preparation error. For the declared volume,
`p_beta=-<K_iso>_beta/(3v)` and
`dS/ds=-beta^2 Cov_beta(H,K_iso)`. These equilibrium derivatives compare
separately prepared Gibbs states; isolated unitary driving need not follow
that family. Adding `xi log(det V)I` preserves every flat-source record
but changes the pressure by `-xi/v`, exposing why flat records alone
do not determine the geometry coupling or its vacuum stress. Choosing a
coupling is permitted; identifying it with observer-derived proper
geometry requires an additional argument. The ten-probe Einstein consumer
also needs its common-state modular diamond, area and small-ball estimates.

### Actual field position and joint crystal momentum

The [complete translation/readout construction](../../model/particle_joint_translation_readout/PROOF.md)
uses the periodic polynomial field, its finite word cutoff and its
referenced mixed input. Translation moves every matter and link register.
On the actual canonical purification its admitted measurement also moves
the conjugate reference. These joint translations commute with the same
Hamiltonian and preserve the joint Gauss-zero subspace; their Fourier
projectors need not commute with gauge transformations on every ambient
charged input. Their complete sum is identity, with outside-domain and
zero labels retained. The result is whole-system crystal momentum, not
momentum of a selected particle.

These finite global menus are outside the local-record thermodynamic
theorem: the momentum projectors and the full-number normalization
`N^(-1)` do not have fixed local support. Their complete finite laws
do not identify global canonical purifiers across volumes.

The position effects are `F_x=n_x/N` on `N>0`, with `F_vac=1_(N=0)`.
They sample an occupied mode's site while retaining every number sector.
For `V` cells and sixty modes per cell, put `v=2^(-60V)`,
`A=(1-v)/V` and `C=E[sqrt(n_x n_y)/N;N>0]` for distinct sites under
the full initial tracial Fock law. The owning proof gives its exact finite
binomial sum. Position followed by joint momentum has original weights
`Prob(vac,0)=v`, `Prob(vac,k!=0)=0`,
`Prob(x,0)=[A+(V-1)C]/V` and `Prob(x,k!=0)=(A-C)/V`.
The strict bound `A-C>=60*2^(-60V)` proves nonzero momentum backaction
without postselection. Later interacting records use displacement-dependent
correlations; this initial formula is not substituted for them.

The actual coordinate first moment and momentum generator satisfy finite
Robertson uncertainty. Recorded position variance includes the positive
POVM noise and can exceed the first-moment operator variance. The exact
cyclic wrap term and number-dependent Weyl phase are retained; a finite
commutator cannot equal `i hbar I`. These are verified finite operational
observables, without a universal physical `hbar/2` bound or a particle
dispersion identification.

The [regional modular calculation](../../model/particle_regional_geometry/PROOF.md)
uses the faithful Q*=0 Gibbs preparation rather than the Q*=6 ground.
Restriction to one family's six modes has positive-support ranks
`(1,6,2,13)` and 42 retained zero outcomes. Its modular Hamiltonian has
the induced pair term `log(p_dark/p_pair) B_f`; for the first family
the coefficient is `log(13440/41873)`. A gauge-invariant pair/dark probe
detects it. A source-admitted entropy variation supplies that same slope
with a controlled finite remainder. A number-only modular ansatz fails,
and the constructed tracial conditional expectation cannot preserve this
interacting Gibbs state. The scalar-only regional modular flow is trivial
on its retained algebra even though the global dynamics changes scalar
occupation. Thus a global modular clock cannot simply be restricted and
renamed a local geometric clock.
The identity `tau=-beta s` belongs to this faithful Q*=0 Gibbs state;
it does not furnish a modular clock for the Q*=6 ground or the separate
gauge-coordinate field preparation.

## Declared field programs and complete observer outputs

The [finite observer-output package](../../model/particle_observer_outputs/README.md)
calculates complete interacting experiments for specified constructions.
Its quantum-Higgs neutral and charged experiments share one Hamiltonian,
prepared sector ground and clock. Its gauge-coordinate field is a distinct
declared program; their spectra are not pooled into one particle list.

| Construction | Calculated output and error control | Physical domain |
| --- | --- | --- |
| [Quantum-Higgs transition](../../model/particle_observer_outputs/HIGGS_RECORDS.md) | The literal full-CAR vertex has a complete 307-state conserved sector with two bright and 305 dark states. The same vertex prepares its ground. Complete two-time occupation records recover `Delta=sqrt(d^2+16 Tr(Y^dagger Y))` on a stated unaliased window. | A finite gauge-invariant scalar/fermion transition. Spatial hopping and gauge transport are excluded from this local Hamiltonian. |
| [Charged Higgs records](../../model/particle_observer_outputs/CHARGED_HIGGS_RECORDS.md) | The same ground has two addition lines per internal singular channel and one removal energy. All 1,650 addition and 36 removal states are retained, including 1,614 dark addition states. A compensating charged probe keeps its quantum reference and both original event probabilities. | Actual interacting charged response energies and residues. Their physical dispersion, spin and mass dictionary are not inferred from their internal labels. |
| [Gauge-coordinate field](../../model/particle_observer_outputs/FIELD_ROUTE.md) | A gauge-covariant bounded overlap polynomial acts on the full Dirac Fock space of a finite three-dimensional lattice. The same Hamiltonian defines the local conserved current. A representation-word cutoff gives a rational Dyson-tail bound on complete histories, including references. | The coordinate completion, action and apparatus are permitted choices. Both chiralities are retained; a chiral interacting relativistic limit is a separate output claim. |
| [Positive-smoothing field](../../model/particle_smooth_sign_field/PROOF.md) | A convergent interaction expansion controls complete local quantum records at fixed duration, uniformly in volume. Polynomial, volume and finite-word errors total `3/2^20` at duration one. | A declared action with full quantum electric links and vectorlike matter. Zero smoothing is optional; physical spacing refinement and chirality do not follow from this thermodynamic limit. |
| [Refining causal family](../../model/particle_physical_causality/PROOF.md) | The complete local influence bound is uniform in assigned physical distance for `kappa(a)<=kappa0/a`, including polynomially growing apparatus supports and errors on both compared experiments. | A conservative causal premise; no convergent state, Lorentz symmetry or measured light speed is inferred. |
| [Interacting thermal control](../../model/particle_field_state_response/PROOF.md) | The unchanged 208-state charge-six quantum-flux source has deterministic positive-polynomial preparation, a Gibbs target and all 169 complete flux records on one clock. | A finite referenced thermal preparation. It is distinct from the scalar Gibbs state, the spatial mixed family and a physical vacuum. |
| [Constant-coframe response](../../model/particle_geometric_response/PROOF.md) | Six full-source spatial strain derivatives, compact-family causal control, a frame-rotation identity and complete work/pressure mean detectors on the same 208-state thermal control. | Geometry is a declared coupling. The sixty-column algebraic control uses kappa=1/7; the directional thermal reduction uses kappa=0 on a separate carrier. Proper stress and Einstein premises require their own identification. |
| [Joint position/momentum instruments](../../model/particle_joint_translation_readout/PROOF.md) | Complete `n_x/N` position effects plus vacuum, Fourier projectors of all-register/reference translations, original Born joint weights and finite uncertainty on the periodic field. | Gauge-preserving on the joint invariant preparation domain. These whole-field measurements are not single-particle momentum or a physical mass shell. |
| [Executed spinor/flux readout](../../model/particle_observer_outputs/FIELD_ROUTE.md#executed-common-state-spinor-and-quantum-flux-control) | The native `q=6Y=6` singlet on one periodic cell has exact initial flux, current and energy `3,6,11`. Two signed-flux reads retain all 1,369 or 3,721 outcomes at cutoffs 18 or 30. Exact rational Taylor weights enclose the complete finite record law; the cutoff-30 evolution half-L1 error is below `1.155e-35`. | One interacting quantum coordinate on the declared clock. These values are neither GeV masses nor a demonstration of three-dimensional propagation. |

The [readout contract](../../model/particle_observer_outputs/READOUT_CONTRACT.md)
gives full pointer unitaries and complete instruments, including dark,
unused and outside-encoding outcomes. Their exact matrices may be primitive
apparatus services. The alternative coherent-source compiler gives a finite
resource and reference-stable approximation bound; its long tapes are not
claimed as executed. One fixed evolution macro defines the program clock,
with preparation and detector evolution paused or explicitly included.
The finite Born theorem applies to the same effects and state. Consensus
copies their orthogonal records and preserves the calculated weights.

For the field cutoff, `q_m(H_W)` is bounded on every gauge configuration,
including small Wilson gaps. Accuracy relative to the exact sign
requires the stated gap domain or a bound on the state's weight outside it.
The implementation compresses the polynomial of the parent operator,
rather than evaluating a different polynomial after truncation. Its
two-time rotor weights are exact nonnegative rationals without artificial
renormalization. The comparison combines their Taylor remainder bounds
with the separate infinite-flux cutoff allowance.

The remaining common-theory statements concern interacting spatial
propagation, controlled chiral complementary modes, uniform volume and
refinement limits, and the physical detector dictionary. The atlas's
[named-field meanings](../../model/particle_atlas_horizontal/FIELD_NATURE.md)
distinguish matter occupation, gauge polarizations and scalar fluctuations;
they do not identify the finite readout lines as those laboratory species.

<a id="source-tests-for-the-first-five-atlas-inputs"></a>

## Finite tests of declared atlas inputs

These are comparison calculations for specified atlas inputs. Their input
categories are not additional nodes of the construction graph.

The [construction mission](../../model/CONSTRUCTION.md#what-is-being-built)
allows a declared carrier, software, scheduler, memory and observer actions.
OPH does not require proof of its own uniqueness or that these construction
choices are forced. Source-derived output means the complete observer laws
computed from this fixed instance. Conditional theorems constrain its
consequences; matching nature tests those consequences. A setting chosen
using a measured target is labeled as calibration for that comparison.

The [source-selection calculations](../../model/particle_calibration_selection/README.md)
separate a common unit convention from selection of a physical spectrum.
Their exact results reuse the configured native memories, finite gauge
currents, family-band theorem and complete fermion/link construction.

| Input | Derived result | Required physical connection |
| --- | --- | --- |
| Phase | The zero-momentum free sea is phase-flat; nonzero momentum can favor an interior phase. The tested complete finite Higgs sector favors zero phase. | A fixed program with an actual phase/interaction prescription, and computation of its observer output. A dynamical phase register is optional. |
| Attenuation | The free sea and tested Higgs sector favor r=1 under energy minimization. | A declared retained-loss process and its complete output calculation. A fitted r remains a calibration input; those energy tests do not predict it. |
| Capacity response | The exact C-duty coefficient at every attempt horizon has bounded ratios below the lepton targets under a linear mass attachment. | A physical capacity-to-program map and a justified effective observable. Increasing this archive budget alone cannot supply the missing hierarchy. |
| Clock and units | Native record rates and logarithmic generator gaps have unit-free ratios; one shared positive scale cancels. | Common operational clock and state for the charged action; physical energy conversion for GeV values. |
| Families and chirality | Per-event cost selects the smallest positive rank-three band without a dimension window. Complete gauge-invariant unitary dressing supplies canonical fermions. | A declared cooling/kinetic program, its matter-family output, spatial chirality and mirror decoupling in the interacting limit. |

These results concern their named source domains and reference states.
Calculating a configured phase or attenuation does not require a theorem
selecting that input uniquely. Their physical mass consequences require
the interacting detector and species attachment. These input categories include
overlapping program choices and physical identifications; they are not ten
independent scalar parameters. The [particle explainer](OPH_PARTICLE_ATLAS_V3.svg)
distinguishes finite repair, winding examples and detector records from their
physical field interpretation.

The [common executable hypothesis](../OPH_EXECUTABLE_HYPOTHESIS.md) and
[completion DAG](../COMPLETION_PROGRAM.md) place this work inside one target:
the same source must supply observer clocks, defect identity and energy,
an interacting chiral field description and the realized Einstein premises.
In particular, the photon-containing full theory is gapless, and stationary
repair or record entropy is not the geometric entropy variation needed for
gravity. Einstein dynamics and the pure-gauge gap are deferred branches of
the broader target; they are not prerequisites for the active mass comparison.
The distinctions below retain those separate domains.

The [particle mass response](OPH_PARTICLE_MASS_RESPONSE.md) records the
numerical constructions and hierarchy exclusions; the
[construction contract](../../model/CONSTRUCTION.md#completion-target-and-acceptance-criteria)
sets the larger observational target.

## Minimum particle comparison dependencies

### A restricted comparison within the unified graph

The canonical construction and its status are specified only by
[the ten-node graph](../../model/particle_joint_reconstruction/DEPENDENCIES.md)
and [dependencies.json](../../model/particle_joint_reconstruction/dependencies.json).
Nodes 1-3 establish the finite primary program, physical boundary-state
preparation and same-source complete local-record bound. Nodes 4-8 retain
the joint refinement, interacting relativistic reconstruction, species,
spectrum and independent-comparison requirements. A conditional calculation
may supply a premise explicitly, while preserving its calibration and
source ancestry; that calculation does not certify the supplied premise.
Nodes 9 and 10 reuse the common source and observer-preparation conditions
for the pure Yang-Mills and geometric-gravity branches. Each branch needs
its own admitted source instance or transfer; no duplicated preparation
requirements are added.

The narrower table below expands acceptance for a single proposed
observable within these nodes. It is not a replacement construction graph.

A comparison of one particle observable requires a narrower construction
than the full common-theory target. A specified source can supply a
prediction without a theorem selecting it uniquely among all architectures.
Its action, reference, clock and continuation must be explicit, and the
compared observable cannot also set the parameter being tested. The
Einstein construction and the complete Standard Model are not prerequisites
for this restricted comparison.

| Step | Depends on | Required result |
| --- | --- | --- |
| Source, state and clock | Configured source specifications and admitted preparation | One action and reference, an operational time convention, and any volume or cutoff family used in the claim. Reading configured parameters is allowed. |
| Symmetry-resolved complete correlators | Source, state and clock | Creators and detectors with stated charge and spatial transformation laws; original probabilities for every success, failure and zero record; controlled correlation and Gram errors. |
| Spectral sector and dispersion | Complete correlators | A stable pole with nonzero original residue, or the appropriate threshold or resonance observable, together with finite-volume and cutoff errors and the claimed energy-momentum relation. A finite level or positive broad band alone does not supply these conclusions. |
| Common-clock ratio | Two identified spectral observables on that clock | A ratio with propagated uncertainty. The common time conversion cancels; different references or generators cannot be combined by that cancellation. |
| Physical comparison | Spectral observable or ratio and its symmetry assignment | A justified species/observable dictionary and an independent measured quantity with matching conventions and uncertainty. A source energy ratio is not a species assignment. |

The dependency chain is `source/state/clock -> complete correlators ->
spectral sector/dispersion -> common-clock ratio -> physical comparison`.
A single calibrated energy observable can bypass the ratio step. Connected
higher-point functions or scattering enter when the claim concerns an
interaction, width or coupling; a mass ratio does not require a complete
scattering theory. Electrically charged sectors with an unscreened gauge
field can have an infrared threshold rather than an isolated mass pole,
as in [Buchholz's Gauss-law result](https://doi.org/10.1016/0370-2693(86)91110-X).
A stable neutral sector avoids that particular infrared requirement.

The high-field construction resolves this checklist as follows. Its action,
ground and clock are specified, and a finite native preparation has an
explicit error certificate. Its complete local charge instrument is
constructed. A stable analytic lattice band exists for sufficiently small
`J/B`, with exact second-order dispersion and equal conjugate-charge energies.
The equality gives a source energy ratio of one, not a measured mass ratio.
The local massive relativistic refinement is excluded under bounded `J/B`
and bounded rest energy; physical species labels remain unassigned. The
ordered weak-field family avoids this exclusion and is the continuation
target. It requires source spectral concentration or a controlled
shrinking-tube weight bound before physical comparison.

The source constructions supply distinct parts of this chain:

| Construction | Available input to a comparison | Scope of the identification |
| --- | --- | --- |
| [Native reflected-record transfer](../../model/native_record_response/TRANSFER.md) and [short-record reconstruction](../../model/native_record_response/SHORT_RECORD_RECONSTRUCTION.md) | Positive finite spectral weights, exact recovery of selective maps, and Gram/insertion error bounds from complete records. | Transfer energies use the reconstructed record clock. The [native clock comparison](../../model/native_clock_comparison/PROOF.md) prevents silently identifying it with the working carrier generator. |
| [Cubic free field](../../model/continuum_free_limit/README.md) and [ordered-spin field](../../model/spin_canonical_field/README.md) | A controlled free relativistic limit on its specified scaling; a positive-energy canonical spin field with finite propagation on its own reference and clock. | The free construction does not establish a surviving interaction; the ordered-spin result does not establish a sharp relativistic branch. |
| [Bound pair](../../model/bound_pair/SPECTRUM.md) and [moving pair limit](../../model/moving_bound_pair/README.md) | An isolated two-body upper pole with original contact weight; a separate specified scaling gives a moving quadratic continuum sector. | The fixed-spacing pair's energy/curvature ratio disagrees with the proposed relativistic relation. The two-body reference and shifted, reversed edge clock are explicit in the continuum construction. |
| [One-dimensional charge gap](../../model/mass_generation_audit/README.md) and [relativistic charge scaling](../../model/hubbard_relativistic_sector/README.md) | An exact source-frame Hubbard join, interaction-generated charge threshold and controlled massive dispersion on a supplied running family. | The half-filled reference, centered charge clock and one-dimensional trajectory are specified. This does not supply a three-dimensional particle identification or a local-CAR pole. |

### Configured fermion, kinetic and gauge attachments

The [operator dictionary](../../model/standard_model_source_map/FERMION_ATTACHMENT.md)
attaches the source's fifteen internal states to the configured
`C2_spin tensor M15 tensor C3_family` CAR field. Its exact representations are
`(3,2)_(1/6)`, `(bar3,1)_(-2/3)`, `(bar3,1)_(1/3)`,
`(1,2)_(-1/2)` and `(1,1)_1`, with `Q=T3+Y`.
The same permutation intertwines all twelve currents on all occupation
sectors of the ninety-mode Fock space: `[dGamma(T),a_j^dagger]=sum_i T_ij a_i^dagger`.
The auxiliary even masks, odd propagating fields and three declared family
copies have distinct roles; this finite attachment supplies no mass basis.

The [kinetic attachment](../../model/standard_model_source_map/KINETIC_ATTACHMENT.md)
preserves those currents in every Weyl-node frame. The local sine regulator
retains eight valleys. Opposite winding with the same M15 is a mirror field,
not its conjugate hole: `Tr_M15 Y^5=5/9` changes sign under conjugation.
On its free filled-band reference, a complete local add/remove instrument
assigns each full valley cell original weight `1/8`; this is not a shrinking
physical-energy-window residue.
A single-valley projection fails to reduce the full local gauge-current
algebra. The [quantum-link attachment](../../model/standard_model_source_map/GAUGE_ATTACHMENT.md)
proves full-Fock endpoint Gauss identities for the same M15 Weyl hopping,
and the entire linear/cubic representation anomaly tensor vanishes.
Its dressed finite experiment retains twenty ordinary records and one
other-occupation complement. A selected gauge vacuum, mirror removal and
an interacting chiral continuum are separate requirements; nilpotent quantum
links cannot be replaced by identity background links to import a free band.

The [quantum-Higgs join](../../model/particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md)
uses the existing three-state scalar and its full CAR Yukawa vertex, with
all twelve local Gauss identities. Its rank-one family vertex has exact
norm 2, so complete native feedback admits a common interacting gauge/Higgs
Hamiltonian with `L=2V|g|` at V finite cells. No scalar vacuum expectation is
substituted. A quantum-link Wilson probe has a sign-kernel gap and an exact
Ginsparg--Wilson relation on its one-cell Euclidean singlet sector, but no
light zero. Operator-valued projector entries also produce an extra CAR
commutator. These controls separate finite gauge/Higgs interaction from
the unresolved many-body chiral regulator and physical vacuum.

The [conditional Standard Model benchmark](../../model/particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md)
supplies a complete elementary-field comparison when a renormalized chiral
action, physical state and matching scheme are specified. Three copies of
the fifteen-component representation give forty-five Weyl fields. A
canonical Higgs doublet gives the exact quadratic relations
`m_W=gv/2`, `m_Z=v sqrt(g^2+g_prime^2)/2`,
`m_h^2=2 lambda_H v^2` and `m_f,i=sigma_i(vY_f/sqrt(2))`.
The photon and gluon Lagrangian masses vanish. Antiparticle charges are
conjugate; physical mass equality additionally uses the stated CPT
hypotheses. Neutrino mass fields `nu_1,nu_2,nu_3` retain symbolic masses
and a supplied Dirac-singlet or Weinberg interaction, with flavor states
related by their mixing matrix. The eighteen atlas labels include those
flavor channels, rather than three invented flavor masses.

This benchmark does not identify the finite scalar carrier with a Higgs
vacuum or remove the sine regulator's mirrors. A source derivation must
extract temporal and spatial kinetic kernels and self-energy from common-state
correlators, control a gauge-covariant chiral regulator, and match the
resulting charged thresholds, unstable resonances and quark mass schemes.
Defect identification also requires a dynamical charge/current intertwiner;
the compact texture's localized frequencies do not provide a momentum shell.
Measured reference values label the conditional comparison and are not all
inserted into a tree action as independent exact predictions.

The [three-input charged controller](../../model/particle_candidate_atlas/CALIBRATED_TRIPLET.md)
is a separate constructive calibration within the complete attenuator.
Using the supplied electron, muon and PDG2026 tau central masses gives
`r=0.99999669509849`, `theta=0.22222476189459` and common coefficient
`s_GeV=6.7791801525 GeV`, with all sixteen flags retained. Certified rational
controller values realize those three action energies within `10^-36 GeV`;
this is implementation precision, not experimental or theory accuracy.
All three masses are inputs, so there is no held-out prediction. The
one-ratio balanced Koide comparison instead holds tau out. Neither
calibration computes the interacting physical matching or selects its
controller from primitive repair dynamics.

The [feedback limit](../../model/standard_model_source_map/FEEDBACK_LIMIT.md)
constructs the four specified static flavor actions from weak correlated
record pulses at fixed finite volume. It compares the same complete Markov
archive, with original probabilities, to the mean action. Its coherent time
`t=N h` does not include fixed-angle current-gate and repair overhead;
a common implementation clock must supply that scheduling comparison.
The [mass readout](../../model/standard_model_source_map/source_mass_readout.py)
exports rational intervals for `m_f=|v_H| sigma_f(Y)` and their ratios.
These are quadratic action coefficients, with no interacting pole or
physical unit identification. The [same-action Koide test](../../model/standard_model_source_map/KOIDE_ATTACHMENT.md)
excludes `Q=2/3` for each of the four charged-lepton duty policies; the
separate finite GNS balance does not change their source action. A
[coherent alternative](../../model/standard_model_source_map/KOIDE_SOURCE.md)
constructs `Y=C_ex^2/3` with Koide balance. Its
[full-CAR realization](../../model/charged_lepton_response/CAR_KOIDE.md)
lifts the complete analyzer unitary to the enlarged Fock space. Fresh
ancillary vacuum, a specified charged-Higgs kick and a retained quantum
archive give the same physical ninety-mode pairing Hamiltonian on every
occupation sector. Recording vacuum/nonvacuum flags gives complete-output
unhalved trace error at most `108 |v_H|^2 t^2/N`. The ancillary copies carry
the same gauge charges; their flag coordinate is neutral. Pulse clocks,
fresh-mode preparation and analyzer overhead are explicit service inputs.

The [free relativistic theorem](../../model/charged_lepton_response/RELATIVISTIC_MASS.md)
uses this same coefficient, the native Weyl kinetic symbol and a chosen
filled-negative-energy Gaussian reference. It proves
`E_j(p)^2=c^2|p|^2+|v_H|^2 sigma_j(Y)^2` and controlled compact-momentum
lattice and correlation limits. These are charged-lepton family branches
of the configured free theory, not three thresholds of its whole electric
charge sector. All eight regulator valleys and their mirror-current
boundary are retained.

The [phase readout](../../model/charged_lepton_response/PHASE_AND_PREMISES.md)
recovers the unordered action spectrum from one original eight-flag score;
ordinary measured repair maps alone are phase-blind. At phase zero the
action ratios are `1:1:(34+24sqrt(2))`. A
[calibrated comparison](../../model/charged_lepton_response/CALIBRATED_SANITY.md)
uses only the supplied muon/electron ratio to fix a different program phase,
then gives a tau/electron ratio within one quoted standard uncertainty of
the held-out reference. This is a conditional one-input comparison, not
phase selection or exact agreement with every central mass value.

The [mass-response bound](../../model/charged_lepton_response/MASS_RESPONSE_BRIDGE.md)
propagates proved kinetic, mass-correction and separate observable-matching
errors to ratios and Koide's quotient. Exact permitted corrections to the
native coefficient show that Koide need not be preserved. A stable charged
infraparticle can be described by a lower spectral boundary; unstable
charged leptons require an appropriate resonance observable. The selected
interacting gauge vacuum, mirror treatment, physical clock and those
matching errors are not supplied by the free theorem.

### Capacity-dependent charged candidates

The [capacity readout](../../model/charged_lepton_response/P_DEPENDENCE.md)
evaluates the retained MCPR law on four authenticated coordinate branches,
without fitting or repeating the implicit-root calculation. Its phase is
`delta=2/9+zeta(P,u)` and its amplitude is `rho=sqrt(2) exp(-chi(u))`.
The complete C-duty inverse realizes `spec(g Z_C)={m_j/v}` with two
admissible controllers and all disabled records retained. This proves
existence for the declared law; it does not select that law, its rates or
its physical clock from the primitive repair process. The real pixel
coordinate P is distinct from a total cell count or hardware memory cap.
The rounded ratio centers are:

| Conditional source branch | m2/m1 | m3/m1 |
|---|---:|---:|
| Forward coordinates | 206.768300595181 | 3477.365536596023 |
| Measured-alpha endpoint | 206.768282739123 | 3477.365262268937 |
| Supersymmetric proxy, beta=4 | 206.768297922262 | 3477.365495531309 |
| Minimal census, beta=7 | 206.782621310399 | 3477.644439943646 |

Only the measured-alpha endpoint lies inside both quoted individual
CODATA2022 one-standard-uncertainty ratio intervals. It inherits measured
alpha, and the architecture has target-informed provenance. The closure
roots are unique only within their admitted intervals and prescribed
census, scale and response laws. Their positive-root quotient is
`Q=(1+exp(-2chi))/3<2/3`, so copying their phase into the balanced coherent
analyzer would change the action. No interacting correction or joint
statistical significance is supplied by these comparisons.

The [compatibility calculation](../../model/particle_candidate_atlas/COMPATIBILITY.md)
allows a balanced bare coefficient and a P-deformed effective response.
The [complete attenuator](../../model/particle_candidate_atlas/BARE_EFFECTIVE_CONTROLLER.md)
realizes `C_r=aI+r b(zR+conjugate(z)R^2)` and `Y_r=C_r^2/3` with all
sixteen analyzer flags and a full-CAR feedback limit. Supplying
`r=exp(-chi)` implements the stipulated P-law amplitude; it does not derive
chi from passive memory or an interacting self-energy. The unchanged
[native response](../../model/particle_candidate_atlas/NATIVE_RESPONSE.md)
has finite-step noise, and its actual all-lag diffusion does not close on
literal scalar C3 attenuation in the fixed source frame. Common units and
passive field changes preserve normalized masses. A physical correction
test must extract kinetic and mass responses on the same state and clock,
including the phase difference between calibrated and P-dependent programs.

The [phase theorem](../../model/charged_lepton_response/SOURCE_PHASE_SELECTION.md)
shows that reversing EC to CE conjugates the phase and preserves unordered
masses. The exact ordering potential is `I tensor i kappa(R-R^-1)` and
comes with an order-h dissipator. Since it commutes with R, its conjugation
fixes the coherent mass coefficient; it does not choose a branch phase.
The [eighteen-channel candidate map](../../model/charged_lepton_response/CANDIDATE_MAP.md)
retains the fifteen-weight current dictionary and the two charged-action
families, while keeping the eleven protected objects separate. Mass modes
are spectral subspaces of the whole action, not its constituent port
projectors. Physical masses, spins, widths and defect-to-species links are
unassigned pending the same-source interacting observable attachment.
The [particle atlas](OPH_PARTICLE_CANDIDATE_ATLAS.svg)
shows these source results and the separate physical-identification steps.

### Correlation-matrix spectroscopy and a high-field sector

For one positive ground generator, write the connected Euclidean matrix as
`C_ij(t)=sum_n exp(-E_n t) z_in conjugate(z_jn)`. A positive metric
`C(t0)` and independent overlaps with the first N states define the
generalized eigenvalue problem

```
C(t) u_n = lambda_n(t,t0) C(t0) u_n.
```

[Blossier et al., sections 2 and 3](https://arxiv.org/abs/0902.1265)
give effective-energy errors `O(exp[-(E_(N+1)-E_n)t])` for `t0>=t/2`,
under their separated-level and overlap hypotheses. The finite source
transfer and [bounded response hierarchy](../../model/spin_dynamical_memory/SPECTRAL_BRANCH.md)
supply correlation and Gram constructions to which this method applies.
A certified numerical interval also needs the overlap conditioning,
tail constants and measurement errors. Large-time exponential decay can
locate a threshold without proving a pole. Volume and cutoff control are
separate from finite-matrix eigenvalue convergence.

The [high-field source proof](../../model/high_field_particle/PROOF.md)
matches the declared `H_B=J sum(S_x.S_y-1/4)-B sum epsilon_x S_x^1`
to Yarotsky's weak-perturbation theorem. Onsite frame changes give a
product Hamiltonian with simple excitation 1 and interaction norm at most
`3J/B`. For sufficiently small positive `lambda=J/B`, its ground-relative
generator has an isolated nonconstant analytic lattice band. Fixed-field
locality identifies the theorem's boundary convention with the source
ground. The exact [source matrix tests](../../model/high_field_particle/test_algebra.py)
include every virtual three-flip state and subtract the interacting
ground energy. The bulk result is

```
E(k)/B = 1+3lambda+(lambda^2/2)[3-(sum_j cos(k_j))^2]+O(lambda^3),
E_min = B+3J-3J^2/B+O(J^3/B^2),
D^2 E(0) = (3J^2/B) I_3+O(J^3/B^2).
```

The proof's analytic-kernel argument gives a qualitative uniform `C^2`
remainder; neither its constant nor the smallness threshold is numerically
certified. The complete local preparation has transformed Kraus operators
`K_s=|1><0|` and `K_f=|1><1|`. Its original particle-band weight tends
to one as `lambda` tends to zero. Preparation failure and success outside
the band are retained; the infinite-volume band projector is a mathematical
spectral split, not an asserted exact native detector. The
[complete readout](../../model/high_field_particle/readout.py) has an abstract
full-Fock completion preserving fermion grading. The literal native pointer
agrees on the protected spin code; it has a different out-of-code extension.

A [finite residual certificate](../../model/high_field_particle/preparation.py)
avoids the long adiabatic schedule for the eight-site cube at `J=1,B=400`.
The product reference is within half-trace distance `1/197` of its unique
interacting ground. The complete native preparation, pointer and feedback
execute in 167 pulses, retaining 608 first-failure slots. Exact compiler
bounds keep the full output within `1/100` of the exact-ground instrument.
This is a finite preparation result, not an evaluation of the unspecified
infinite-volume band threshold or a particle-pole residue floor.

The unbroken internal `U(1)` gives opposite sublattice charges and equal
folded bands. These even-CAR spin creators do not supply a lepton spinor
or photon helicity assignment. The isolated band satisfies the locality
and spectral hypotheses for lattice wave operators in
[Bachmann, Dybalski and Naaijkens](https://arxiv.org/abs/1412.2970);
nontrivial scattering amplitudes and asymptotic completeness are unproved.

The same source has gap at least `B/8` and propagation speed at most `72J`.
For spacing `a`, a common positive clock factor `s_a`, bounded rest energy
`s_a E_min<=Mbar` and `J/B<=lambda_bar<lambda_*`, the band obeys
`|grad_p[s_a E(k_a+a p)]|<=576 a lambda_bar Mbar`. Every finite pointwise
limit on a connected momentum window is therefore constant. This excludes
a nonzero-speed massive relativistic refinement under that local spatial
and common-clock identification. A fixed lattice cutoff, a different
justified dictionary, or a trajectory with unbounded `J/B` has different
premises. Configured source parameters are admitted inputs; they do not
by themselves assign a physical species.

The specified joint family `p=n^-2, B/J=n^-5` avoids the high-field bound,
but its positive directed energy bands do not prove a sharp shell. Use its
positive scaled clock K and one filter R supported in a fixed interval
`[ell,Lambda]`, with `ell>0`. Put `q=RQ Omega`, `a=RP Omega`,
`h_Q=<q,Kq>`, `chi_Q=<q,K^-1q>`, and analogously for P, with
`c=Im<q,a>`. The [exact Ward defects](../../model/spin_joint_field_limit/WARD_DYNAMICS.md)
are `delta_Q=h_Q-c^2/chi_P` and `delta_P=h_P-c^2/chi_Q`.
Their vanishing tests the two filtered force equations; the filtered cross
moment c must be measured on the same original response. Sharp corrected
oscillators and [two-frequency controls](../../model/spin_dynamical_memory/RESOLVENT_CLOSURE.md)
both satisfy the known static moment and symmetry bounds. The strict
[quantum force slack](../../model/spin_joint_field_limit/QUANTUM_FORCE_SLACK.md)
excludes the naive mean-spin coefficients, without proving a linewidth or
excluding a sharp subbranch. Susceptibility saturation or a shrinking-tube
bound with nonzero original weight is an additional source requirement.

The [weighted inverse-susceptibility test](../../model/spin_joint_field_limit/SPEED_CONCENTRATION.md)
reduces a sufficient concentration test to one dynamical statistic. On the
same positive-field source and a momentum packet away from zero, let
`X=K_B/[J sqrt(2 sum_j(1-cos P_j))]`, `w=||psi||^2`,
`H=<psi,X psi>` and `chi=<psi,X^-1 psi>`. Exact double commutators express H
in static bond/order correlations; equal-time correlations give w. For a
specified speed c, `D_c=H-2cw+c^2 chi>=0`, and original spectral weight
outside `|X-c|<epsilon` is at most `(c+epsilon)D_c/epsilon^2` without an
upper energy cutoff. An actual preparation `psi/Lambda` contributes the
known factor `1/Lambda^2`; no observed success probability is divided out.
The complete joint energy/momentum score and the positive-B gap give finite
chi intervals. The source static inputs and chi require growing-family
enclosures; the test does not establish `D_c/w->0`. A nonzero full-response
defect does not exclude a sharp component with a continuum beside it.
These results leave the canonical premise classifications unchanged.

## The three field constructions and their source boundaries

The sphere occupation field, the relativistic free field and the ordered-spin
field use different observable dictionaries and scaling trajectories. Their
spatial dimensions do not replace one another.

| Construction | Established field result | Identification boundary |
| --- | --- | --- |
| [Native occupation records](../../model/physical_space_time/OCCUPATION_FIELD.md) | A two-dimensional spatial Dirichlet-form limit on the carrier sphere, with its stated observation clock. | The local rank-three carrier frame is not an extra spatial dimension of this particular field. |
| [Centered Weyl/contact source](../../model/continuum_free_limit/README.md) | A free relativistic field in three spatial dimensions and time, retaining all eight valleys, on the proved weak-contact scaling. | The cubic field action, encoding and state are specified; this branch does not supply a surviving interaction or measured massive species. |
| [Ordered protected-spin source](../../model/spin_canonical_field/README.md) | A regular canonical field net in three continuous spatial dimensions and time, positive energy, continuous translations, finite propagation and nonzero gapless spectral weight. | Its microscopic dynamics interacts. Lorentz covariance, a sharp particle branch and nonzero connected continuum scattering are unproved. |

The ordered-spin chain uses one declared two-flavor field extension throughout.
[Protected spin dynamics](../../model/protected_spin_limit/README.md) gives
`H_J=J sum_<xy>(S_x.S_y-1/4)` from the centered source along
`c_E=J E^4, U_E=J E^8`, with complete repair records. This is distinct from the
bare one-bit occupation source. [Finite ordered preparation](../../model/native_ordered_ground/README.md)
constructs the ground of `H_B=H_J-B sum_x epsilon_x S_x^1`; the selected infinite
state takes volume to infinity and then `B` to zero. The canonical field uses
`X=p*x`, `T=J*p*t` and the positive ground generator `K/(J*p)`.
The finite native construction has an explicit precision budget. The selected
infinite state has no proved effective finite-box preparation modulus.

### Interacting matter on one fixed-spacing source

The [local full-CAR construction](../../model/source_family_compatibility/LOCAL_MATTER_READOUT.md)
keeps fixed c,U>0, the centered full ground-state subsequence, odd fermionic
observables and internal spin on the same positive clock. Complete local
observer histories approach the coherent source histories under uniform
global repair attempts. For m attempts, a radius-R collar and k local
interventions, their half-trace error is at most
`min(1,m[6(2R+1)^3/E+min(2,(2k+1)B_R(T))])`, with the source locality
error B_R given in the proof. Remote raw records are traced in this
accessible output and cannot control its later feedback. This comparison
does not identify the complete global archives.

Each local mode on that state has original addition/removal weights 1/2 and a nonzero
connected four-point response at fixed spacing. The
[charged sum rules](../../model/source_family_compatibility/CHARGED_SPECTRAL_BANDS.md)
give each charge an original weight at least `25/100352` in
`35/256<=E_source^2<=2744/5` at c=1/8,U=1 in the cubic bulk, uniformly
in volume. The interval is a positive-energy response bound, not a gap
or particle pole. Complete finite-resolution energy records have the local readout
construction below; an ideal sharp spectral projector is not a local
intervention supplied by this comparison.
This retains interacting odd matter without the spin projection that sends
its charged energies to infinity. It supplies no interacting spatial
continuum or effective ground-subsequence selection modulus.

The [actual composite fermion](../../model/source_family_compatibility/COMPOSITE_CAR.md)
`b=2(n_opposite-1/2)a` extends the response without adding independent
modes. On the specified proper-cubic average of the invariant full ground,
the [magnetic-momentum reduction](../../model/source_family_compatibility/MOMENTUM_RESPONSE.md)
uses the actual commuting two-step magnetic translations. With
`h(k)=c sum_j cos(k_j) Gamma_j`, `S(k)=sum_j cos(k_j)^2`, its first
composite block and force Gram are

```
C(k)=kappa h(k),                 kappa=4 q1,
G(k)=[g0+2g2 sum_j cos(2k_j)] I,
g0=c^2(9/2-12qstar),            g2=c^2 q2,
R*R(k)=[g0-6g2+(4g2-kappa^2 c^2)S(k)] I.
```

Here q1, q2 and qstar are explicit two-site and seven-site spin/eta
correlations of that same ground. Complete local Pauli records measure
them with stated errors. Exact source symmetries remove the other
range-two coefficients and give `0<=G(k)<=18c^2 I`. The scalar Gram does
not fix the frequency-dependent memory or a particle pole. Averaging the
ground over the proper cubic group is explicit; a symmetry-broken pure
state is not identified with this reference.

The [nodal source estimate](../../model/source_family_compatibility/NODAL_MEMORY.md)
strengthens these constraints: at U/c=8 the actual residual Gram is at
least `18c^2/125` on `sum_j cos(k_j)^2<=3/250`. Local source derivatives
also bound the dissipative memory below at each positive imaginary
frequency. Neither the kinetic node nor its two-block compression is
therefore sufficient to identify the full interacting dispersion.

The [every-momentum spectral construction](../../model/source_family_compatibility/SPECTRAL_RESIDUE.md)
uses finite-range CAR moments and uniform exponential tails to define the
actual matrix measure at every magnetic momentum, with the required Bloch
chart gluing. It is weakly continuous; its resolvent and memory are
continuous at each fixed nonreal frequency. This supplies no modulus when
frequency and distance to a node vanish together. At the magnetic node
`k_*=(pi/2,pi/2,pi/2)`, the [exact source little group](../../model/source_family_compatibility/NODAL_CHANNELS.md)
first reduces the response to two scalar signed channels. Site inversion
inherited from the same finite-box and averaging prescription exchanges
them at the same energy; particle-hole symmetry makes the common measure
even. Thus the full sixteen-component nodal measure is `mu I`, with no
additional state average. The eight transported Weyl-valley probes retain
their kinetic winding and share this same interacting response. They do
not define eight independently selected particle poles or mass laws.

More explicitly, `Sigma(i eta)=-i b_eta I`, with `b_eta>=0`,
`alpha=U/2` and `m_eta=eta b_eta`. The original nodal Lorentzian score is
`[(eta^2+m_eta)/(alpha^2+eta^2+m_eta)] I`. Its zero-energy residue is
`m_0/(alpha^2+m_0)`, with `m_0 I=R*1_{0}(L_perp)R`.
The positive total Gram does not determine this zero-energy component.
The same source does prove a fixed-unit lower-energy event: for `c=1,U=8`
and packets supported on `S<=delta0`, where
`delta0=1296/843578984375`, the combined original measure has weight at
least `2592/14340842734375` on `|E|<=sqrt(16-delta0)`.
Zero energy is included. Restoring c multiplies this energy window by c;
the bound does not supply weight in the shrinking unit-hopping window
required by ballistic scaling.

For one fixed selected c=1,U=8 ground, weak momentum continuity also gives
an [original-weight necessity test](../../model/source_family_compatibility/SPECTRAL_RESIDUE.md).
Under c=v/a, normalized packets supported within `aK` of the actual magnetic
node have limsup weight in any bounded physical-energy window at most
`mu({0})=m0/(16+m0)`. Vanishing packet weight outside that cap is allowed.
If m0=0, their finite-energy original CAR response vanishes. The
[source locality and Bernstein bound](../../model/source_family_compatibility/NODAL_LOCALITY.md)
gives the converse for bounded-window weight: a positive nodal atom
retains uniformly positive original weight in some bounded physical-energy
window. A certified atom floor gives an explicit spacing and score-error
budget. This does not supply positive-energy dispersion or a particle.
The [full-CAR spin/eta star](../../model/source_family_compatibility/NODAL_MOMENT_BOUNDS.md)
also gives the all-momentum bound `G(k)<=c^2(31+sqrt(577))/4`. At U/c=8,
the nodal atom and every original packet's ballistic fixed-window limsup
are at most `(31+sqrt(577))/(95+sqrt(577))`. The atom's positivity remains
unproved. These statements fix the
c-normalized reference; they do not transfer to arbitrary a-dependent
ground choices without an additional uniform estimate.

The [default-cube certificate](../../model/source_family_compatibility/CUBE_RESPONSE.md)
has `3/80 I<R*R<61/1000 I` on all sixteen original mode directions of the
full 65,536-state source. Its three response levels span 48 dimensions,
without adding particles. The [complete finite energy table](../../model/source_family_compatibility/CUBE_SPECTRUM.md)
retains both charges and every clock outcome. These spatial preparations
have no open-cube translation-eigenvalue interpretation.

The [local spectral instrument](../../model/source_family_compatibility/LOCAL_SPECTRAL_READOUT.md)
reconstructs the same ground's original positive addition/removal energy
measures using time-evolved CAR operators. Fourier orthogonality makes
all charge/phase outcomes complete on every input. Its radius-R cut has
complete source/reference half-trace error at most `sqrt(2) B_R(T)`, with no factor
from the number of phase labels. A supplied local ground-marginal error
and the declared coherent-service errors are additional. The
[local commutator bounds](../../model/source_family_compatibility/UV_TAIL.md)
give volume-independent exponential energy tails for every normalized
packet, controlling high-energy aliases. At c=1/8,U=1 a tapered clock
with 16384 retained records certifies each charge's original weight above
`75/401408` in `[0,(13/2)pi+1/4]`, on the same full ground. This large
readout is specified analytically; its native tape is not executed.
Bounded Lipschitz spectral scores converge through the stated ordered
local limits, without assuming a sharp pole or preparing an unspecified
ground by declaration.

The [direct Lorentzian readout](../../model/source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md)
uses four complete charge/sign records at each geometric time and two
dummy records. For normalized packets f and g, `r=(1+eta Delta)^-1` and
`M2<=v^2`, its signed score differs from the original
`Q_eta=integral eta^2/(eta^2+E^2) dnu_f` by at most

```
2||f-g|| + Delta v/2 + r^m
+ 2 min(1, B_R(T)+delta_marginal+epsilon_prep+epsilon_service+epsilon_repair),
T=(m-1)Delta.
```

The collar error is independent of the `4m+2` record count, and failures
retain their original probabilities. Algebraic local clock approximations
make the finite score an even algebraic observable. Computable pressure
secants then bound this observable directly over all relevant
translation-invariant ground states, without substituting the ground of
the isolated cut. The strict stopping search terminates when every such
ground has a common score margin; it does not establish uniqueness or an
effective growing-marginal modulus. At default c=1/8,U=1 the local
all-input second-moment identity `M2=35/128` gives the stronger analytic
bound `nu([0,1])>=93/128`. A specified 8194-record readout and pressure
search certify the operational lower bound `>15/32`; that large native
tape and search are not executed. The separate [native dimer execution](../../model/source_family_compatibility/NATIVE_LORENTZIAN_READOUT.md)
retains 8194 proper outcomes and all raw-pointer and first-failure
records. Exact spectral intervals and a complete-channel bound below
`1/100000` certify original weight above 0.595 in `[0,1]` on its full
16-state source. These fixed-source statements do not give a positive
scaling-uniform low-energy probability.

The same cubic reference has [local sixth and eighth nodal moment observables](../../model/source_family_compatibility/NODAL_NEXT_MOMENTS.md),
with 493 and 27965 Hermitian words after its actual state symmetries.
Their selected-ground expectations are undetermined. Complete weighted
word instruments retain both signs; ground-marginal and native-service
errors multiply the observable's coefficient norm and cannot be omitted.
The resulting degree-four Christoffel bound uses those actual moments.
The onsite force identity alone gives no ground-norm equality or charge gap.

The [strong-contact bound](../../model/source_family_compatibility/STRONG_CHARGE_SCALING.md)
excludes every finite-energy original odd response on ballistic c=v/a
with U/c tending to infinity, including fixed positive cell contact
`g=a^3 U`. At U/c=8 the
[complete-response obstruction](../../model/source_family_compatibility/LOCAL_CHARGE_ESCAPE.md)
is stronger than a finite-cube comparison: on the specified cubic
invariant ground average, every normalized spatial packet has limsup
original weight at most `9/17` in any fixed physical-energy window.
The entire unfiltered normalized response therefore has no pointwise
clock limit continuous at zero. This does not exclude a finite-energy
filtered component or an alternative collective field. The concrete
route is the complete filtered readout above, retaining its complement;
it requires a positive original weight in shrinking unit-hopping windows
and controlled state, momentum and interaction limits. Neither the
fixed-source band nor its exponential ultraviolet tail proves that weight.

The [actual energy-dipole identity](../../model/source_family_compatibility/ENERGY_DIPOLE.md)
separates the quadratic source current from its exact composite interaction
and centering contribution. The quadratic current divided by v^2 converges
to geometric momentum on a specified low-momentum window with an O(a^2)
error. The interacting term and the window's preservation on complete
histories require control before this dipole supplies an operational boost.

The signed Hubbard action and its repair use one declared occupation frame.
Changing to the site-dependent Weyl frame must conjugate every probe and
repair operator. Bare repair in the other frame is a different protocol.
The local-history estimate also applies directly to the original Weyl action
with bare original repair: it uses only even seam support and the source
bond norm c. No equivalence of the two repair protocols is assumed there.
The [aligned-hole construction](../../model/source_family_compatibility/ALIGNED_HOLE_DISPERSION.md)
derives an exact pi-flux dispersion and an operational position/momentum
limit in a separate fully aligned reference sector. Its Schrödinger limit
is free and nonrelativistic; its reference energy cannot replace the
centered interacting ground in the charged construction.

### Ordered spin preparation

At strictly positive selection field the
[uniform gap and native ramp](../../model/spin_twist_closure/FINITE_FIELD_GAP.md)
give a stronger preparation interface. An exact Hadamard change of basis
puts the actual source in the Lee--Yang theorem's class, with full gap at
least `B/8`, uniformly in volume and along a nonnegative bond ramp from
the product state. The adiabatic, time-splitting and native pulse errors
have one explicit complete-output budget. The large tape is specified;
the complete small native compilation witness is executed.

The [positive-field marginal theorem](../../model/spin_fixed_field_locality/PROOF.md)
also gives a computable volume modulus. Interpolating only the bonds across
a box boundary preserves that gap. Gaussian approximate spectral flow
and the actual source locality bound control an interior marginal,
independently of the exterior volume. For a support of m sites at distance
d from M crossing bonds, set `g=B/8`, `T=d/(144J)`. Its half trace error is
at most

```
min(1, JM[(mT/2) exp(-d/3)
          +(sqrt(2T/(pi g))+1/g) exp(-gT/2)]).
```

For fixed B, the surface factor grows polynomially while the last exponent
is `-Bd/(2304J)`. Finite open-source marginals therefore converge to one
positive-field state, with a rational box-size planner. Complete local
instruments preserve the bound, including independent reference systems.
This is a state comparison, not an arbitrary-input channel bound.
The constants diverge as B tends to zero; comparison with the selected
zero-field state requires a separate rate.

### A prepared joint source family

The [joint-family construction](../../model/spin_joint_field_limit/PROOF.md)
uses `p=n^-2`, `B=J n^-5`, the positive-field ground and the same H_B clock
at every scale. This is a specified preparation prescription, not a power
law selected by unread repair. It avoids replacing that state or clock
by the ordered zero-field prescription. At each fixed n, the volume,
adiabatic, native clock and complete local Pauli-record errors have a
computable combined bound. Large certified schedules are resource bounds,
not executed simulations.

On a common bounded-history subsequence it gives two regular canonical
pairs, `(Q3,P2)` and `(Q2,P3)`, with
`Q_a(F)=p^2 sum epsilon_x F(px) S_x^a`,
`P2(G)=-p sum G(px) S_x^2`, and `P3(G)=p sum G(px) S_x^3`.
Both symplectic coefficients are the same source order `m integral FG`,
with `m>1/10`; all cross pairs commute. The
[mixed-momentum proof](../../model/spin_joint_field_limit/SECOND_PAIR.md)
first controls auxiliary transverse phases, then the longitudinal phase,
then the mixed Weyl commutator on complete histories. This avoids assuming
the extra canonical pair from a vacuum commutator alone. Even microscopic
translations give common continuous spatial translations. The full internal
U1 acts on both pairs and commutes with the same positive time generator.
The source locality bound gives a finite propagation cone.

The [directed charged-band construction](../../model/spin_joint_field_limit/DIRECTED_CHARGED_BANDS.md)
uses real-even profiles in two opposite caps of relative radius
`0<delta<=1/4` about `+/-k0`, with `rho=|k0|>0`. The linear canonical
test vector has this spatial support; no such support is assumed for the
nonlinear bounded creator. Actual site inversion gives each internal
charge and each cap original probability at least `1/774400`, uniformly
in rho and delta, in the broad energy interval
`[9rho/774400,580800rho]`. All seven creation/cap/complement records remain.
For a successful unnormalized vector psi of weight w, actual spatial
translations obey
`||(U(y)-exp(i sigma k0.y))psi||^2<=min(4,delta^2 rho^2 |y|^2)w`.
A fixed profile has the inherited finite-source approximation, while a
sharp continuum cap requires its resolution and boundary-weight controls.
The full periodic six-site control retains thirteen charge, translation,
energy and complement records; its directed energy variance is positive.
Neither arbitrarily narrow spatial caps nor that finite variance fixes a
continuum energy shell, physical spin or a Standard Model charge.

The positive-field gap supplies quantitative order variance O(1/n) and
bond-coefficient variance at most `(792/n+24/n^6) a^2 V0` for the owning
coefficient/support normalization. Scalar twist coefficients do not imply
closure of the dynamical field algebra. The construction supplies no
unique continuum history state, effective continuum convergence modulus,
Lorentz transformation or sharp particle branch. Equality with the
ordered zero-field history limit is unproved. The exact
[clock comparison](../../model/spin_joint_field_limit/CLOCK_COMPARISON.md)
controls equal-time selection rotations and states the additional
mixed-history estimate required for that equality.

The [actual filtered force](../../model/spin_joint_field_limit/FORCE_READOUT.md)
isolates a quantitative dynamical target. For the positively oriented
source current `Z_p(G)=p^3 sum D_i^pG C^2_(x,i)`,
`i[L_B,p,P2(G)]=-p^-2 Z_p(G)-(1/n)Q3(G)`.
For a contraction energy filter rho supported in [0,M],
`||rho(L_B,p) Z_p(G) Omega_B|| <= p^2(M||P2(G)Omega_B||
+ ||Q3(G)Omega_B||/n)`.
A sufficient source-norm route to the proposed linear force coefficient a is
the stronger cancellation `||rho(L_B,p)[Z_p(G)+p^2 a Q3(Delta G)]Omega_B||=o(p^2)`.
The O(p^2) current bound does not prove it. Complete finite Pauli scores
measure the unnormalized squared residual with all outcomes retained;
the error must also be o(p^4) to resolve this test. This is a proof and
measurement interface, not a measured mass. Here P2/P3 are canonical
field momenta, distinct from the generator of spatial translations.

[Clock comparison and energy continuity](../../model/native_spin_clock_join/README.md)
control complete histories of `H_B` and `H_J` on the same prepared input and
construct finite local energy currents. They do not identify that finite input
with the infinite ordered vacuum. [Relaxed response](../../model/spin_relaxed_response/README.md)
constructs the inverse-energy variational response in that vacuum without
identifying bare twist energy with its full continuum energy.
[Complete shell readout](../../model/spin_shell_readout/README.md) bounds the
original spectral weight near a specified energy-momentum relation, with
clock leakage, aliases, failed preparations and sampling errors retained.
A finite energy band is not a proof of a continuum particle pole.

The [bounded response-closure test](../../model/spin_dynamical_memory/RESOLVENT_CLOSURE.md)
uses the actual positive clock K and source-preparation isometry V.
With `S=z(K+z)^-1`, `A=V*SV`, the positive defect
`D=V*S^2V-A^2=R*R`, `R=(I-VV*)SV`, vanishes exactly when that response
space reduces K. It needs no unbounded energy moment. A dispersion claim
also requires independent identification of A with `z/(w(P)+z)` on a
space reducing the actual spatial translations. Finite Gram/resolvent
records supply a test only within their stated access and error bounds.
Small D alone gives bounded polynomial comparison, not unqualified time
evolution: high-energy leakage needs its own cutoff/tail estimate.
The actual source dimer distinguishes a nonreducing local preparation
from its reducing two-energy span. This is not an actual continuum
vanishing-defect theorem. A two-speed field with two CCR pairs and full
internal U1 shows why those kinematic inputs and scalar twist coefficients
cannot select one field equation.

Two exact interfaces reduce the remaining calculations. The
[source Ward identities](../../model/spin_joint_field_limit/WARD_DYNAMICS.md)
express the optimal two-force defects as `h_Q-c^2/chi_P` and
`h_P-c^2/chi_Q`, using source first and inverse energy moments.
On a fixed positive energy window these give original-weight tube bounds.
Their saturation is an unproved source estimate. The actual
uniform finite-spin Ward relation fixes one susceptibility exactly, but
does not by itself imply a sharp energy distribution. Filtered cross
moments and unfiltered order commutators remain distinct.

For the uniform finite-field Ward pair, write `q=R M3 Omega` and
`M_j=integral E^j dnu_q(E)` on a positive energy window. Its second
force defect is exactly zero, but the first is
`delta_Q=M_1-M_-1^2/M_-3=M_-3 Var_pi(E^2)`, with
`dpi=E^-3 dnu_q/M_-3`. The unfiltered susceptibility fixes only its own
M_-1; filtered moments require their corresponding readout.
The [uniform-mode estimate](../../model/spin_joint_field_limit/UNIFORM_MODE.md)
proves concentration at zero on the scaled clock: for source order m,
`R<=6/(mn)+n^-6`, `w^2<=R` and `Var_pi(E^2)<=R^2/4`.
The fixed-positive-field marginal theorem supplies explicit growing boxes
with m>1/20. This volume-wide mode is distinct from compact field probes
with nonzero momentum. Its proved original success floor `1/(160N)`
depends on the box volume N and supplies no positive continuum residue.
The complete eight-site controls at n=2,3,4 have positive weighted
variance even after charge resolution. The same fixed box has scaled
gap at least `n^2 g0-8/n^3`; it cannot substitute for the growing-box
limit. Momentum-resolved dispersion and original continuum weight are
the additional estimates needed for particle identification.

The [source spectral hierarchy](../../model/spin_dynamical_memory/SPECTRAL_BRANCH.md)
uses complete time-filtered preparations to match `2r-1` bounded resolvent
moments on an r-step response space. For a Lipschitz score f its original
error is at most `w Lip(f)/sqrt(2r-1)`. The hierarchy therefore constructs
controlled continuous spectral scores without assuming that a response
projection is executable or its Ritz values are physical energy levels.
Window boundaries, growing preparation costs and continuum errors remain
explicit. The [threshold route](../../model/spin_shell_readout/THRESHOLD_READOUT.md)
allows a continuous infrared edge if a source tail modulus and the
appropriate physical detector theory are supplied; a pole is not mandatory
for every charged physical sector.

For these spin prescriptions the precise additional requirements are:

1. Select a common preparation/refinement prescription and give an effective
   continuum history error for the growing spectroscopy supports and times.
   The ordered-zero-field route needs its B-to-zero marginal modulus.
   The direct joint H_B route has complete fixed-scale preparation and
   local histories, but no effective convergence rate to its field limit.
   The [interacting repair budget](../../model/source_family_compatibility/DEFECT_BUDGET.md)
   controls complete repair disturbance on a specified full Hubbard
   trajectory from its propagated occupation-defect expectation. It permits
   doublon--hole creation; the selected state and energy/attempt scaling
   remain required inputs.
2. Prove the ground-specific inverse-current tail bound, or another estimate
   that controls the continuum energy and current responses. Finite positive
   relaxed susceptibility does not establish that tail bound.
3. Prove a common dispersion branch, controlled resonance, or infrared
   threshold with suitable asymptotic detectors. The shell route needs a
   uniform shrinking-tube lower bound; the continuous-edge route instead
   needs controlled original spectral tails. Neither source estimate is
   established by a finite readout or a fitted spectral centre.
4. Establish Lorentz covariance for the same multitime field hierarchy and
   determine whether its connected interactions survive the scaling.
5. Construct the corresponding momentum flux and symmetric conserved stress,
   then attach the responsive metric, entropy and area data in the finite
   Einstein composition theorem.
6. Identify massive sectors, charges and detector couplings on this theory;
   fix one unit convention and compare complete predictions with measurements.

The [current sum rule](../../model/spin_current_spectrum/README.md) and
[shell diagnostic](../../model/spin_shell_diagnostic/README.md) give the exact
missing correlation estimates. Neither an additional finite theorem with
those estimates as hypotheses nor a fitted dispersion discharges them.

For the ordered-zero-field prescription, two quantitative interfaces
distinguish the missing estimates from the proved representation.
Prepared-state control requires a computable bound
on `half_trace(rho_(box,B)|_R, omega_selected|_R)` for growing region R and
the observation times used by the experiment. Finite ground preparation
and qualitative maximal-order convergence give no such joint modulus.
For the original broad-window response, the existing shell test sets
`X_p=K/(J s(P))`, where `s(P)^2=2 sum_j(1-cos P_j)`, and
`d_p=<X_p><X_p^-1>-1`. Proving `d_p->0` and one common limiting mean,
independent of packet envelope and direction, would concentrate that
whole response on a cone. Alternatively, a shrinking spectral tube can
retain a fixed positive fraction of the original weight while a continuum
survives outside it. The source current sum rule fixes a total inverse
moment; it does not establish either concentration statement.

The [local spin creation](../../model/defect_localization_readout/PROOF.md)
supplies a finite excitation on the actual selected cube ground. A local Z
leaves the exterior marginal unchanged at preparation and flips a global
internal spin parity preserved by every C/E record and the centered source
Hamiltonian. Complete native preparation and detection give contrast above
`0.67`; a source-clock estimate retains contrast above `0.4` through ten
repairs on its declared short-time trajectory. This is a finite symmetry
sector with controlled spreading, without a topological boundary charge or
a particle dispersion theorem.

The [driven spin response](../../model/spin_stress_response/PROOF.md) gives
`H_a=sum_x a_x e_x`, currents `J_xy=a_x a_y i[e_x,e_y]` and work densities
`dot(a_x)e_x` on the same finite action. Exact continuity and retarded
response use the actual initial state. A separate two-site source instance
prepares its selected `J=1,B=3/8` ground and evolves it under the exchange
clock; its one-site marginal is `(I+(3/5)cos(t)X)/2`, with spectral floor
`1/5`. Three complete 144-pulse experiments implement a nonzero regional
entropy/modular response with every failure retained.

For a faithful reference with floor `mu` and half-trace perturbation
`epsilon<mu`, the general bound
`0<=D(rho||sigma)<=epsilon^2/(mu-epsilon)` proves the finite first-law
differential on these source tangents. The physical exchange-energy response
in the displayed experiment is zero despite nonzero regional entropy
response. The result supplies the entropy/modular equality premise of the
finite first-law theorem. The boost/edge split, a symmetric Lorentz stress,
responsive area and fixed-proper-volume stationarity require additional
maps. The two-site reference is not substituted for the infinite ordered
state of that ordered-spin construction.

### Source symmetry, spectral response and regional geometry

The [source-charge construction](../../model/protected_defect_catalog/PROOF.md)
uses the ordered state and the same `K/(Jp)` clock to attach internal U1
to the transverse configuration history subnet. Adding the second staggered
configuration Q2 costs only the source commutator bound
`||[Q2_p(F),Q3_p(G)]|| <= (p^4/2) sum|F(px)G(px)|=O(p)`.
No second uniform canonical momentum is appended. The existing variance,
energy-form transport and full-history regularity arguments construct
the joint limit. A source symmetry exchanges the opposite-charge fields;
their joint energy--momentum measures are equal and inherit the nonzero
positive-energy annulus weights of Q3. This is a field attachment of a
source charge, with no Lorentz-spin or mass-shell assertion.

The [dynamical-memory construction](../../model/spin_dynamical_memory/PROOF.md)
retains the complete spectral complement of a source-created vector:
`G(z)=1/(z+epsilon-Sigma(z))`, with positive memory
`Sigma(z)=<r,(D+z)^-1 r>`. On both finite coherent charge channels this
memory is strictly positive. Bounded odd preparation and a resolvent filter
have a uniform positive original probability in the actual field scaling.
Their filtered norms, first two energy moments and resolvents converge on
the common bounded-history subsequence, without assuming energy-moment
tightness for the unfiltered field. The continuum memory is defined by
these source limits; a sharp shell or vanishing force residual is not inferred.

The [regional modular construction](../../model/spin_modular_geometry/PROOF.md)
uses the same selected eight-cube ground rather than the separate dimer.
Its two-site marginal has spectral floor greater than 1/41. The certified
modular generator and a source-current perturbation give a nonzero entropy
response, executed with all native failure records. An informationally
complete 32-record regional detector constructs a positive Fisher metric
on density-response tangents. This is not a spacetime metric. The actual
regional modular generator lies outside the weighted isotropic-bond plus
local-field ansatz. A fixed graph-cut area has zero response while this
allowed pulse changes entropy, so that fixed area cannot supply stationarity
for the displayed variations. Responsive geometry and the physical boost
identification are separate requirements.

### Bounded charged preparations and separated pairs

The [bounded charge construction](../../model/spin_charge_dispersion/PROOF.md)
uses the existing infrared proof with the trial vector `(I+tM)G`.
Its exact 256-mesh integral bound and fixed bias `t=10` certify the same
selected magnetization `m>1/10`. This is a lower bound, with no empirical
magnetization input. An [entangled-dimer trial](../../model/spin_charge_dispersion/DIMER_BOUND.md)
strengthens the same order bound to `m>11/100` and the microscopic current
coefficient to `1/120`; the readout keeps the conservative value `m0=1/10`.
These are source inequalities, not a dispersion speed or fitted stiffness.
For an annular unit test, set `x=sqrt(r)Q3(F)` and
`y=sqrt(r)Q2(F)`. The angular charge projection of `f(x)=x/(1+x^2)` gives
bounded opposite-charge creators `B_+,B_-` of norm at most1/2. Their
positive square-root complement makes one complete three-record instrument.
The actual canonical commutator, source moments and charge conjugation prove
that each joint creation-and-band event has original probability at least
`k=1/1757600`, uniformly as `r` tends to zero. Its momentum lies in
`r/2<=|P|<=r` and energy in `[3kr,r/k]` on the same positive clock.
This broad band proves positive excitation weight; it is not a sharp
dispersion law. Complete finite compilations and ordered history limits
retain all other energies, momenta and preparation records.

The [separated-pair construction](../../model/source_charge_scattering/CONTINUUM_PAIRS.md)
uses those same bounded configuration fields. For the positive effect
`A=B_+^dagger B_+`, finite spatial averaging and Jensen's inequality give
`omega(A_R^n)>=omega(A)^n>=k^n`. Excluding nearby packet centres yields
pairwise separated rational tuples with original weight at least `k^n/4`
for every fixed n and every choice of charge signs. A computable cube
contains a suitable tuple; no effective ground-state evaluation or good-point
selection algorithm is supplied. All `3^n` preparation records are retained.
At n=2 this gives charge+2 and charge-2 probability at least `k^2/4`. The energy
form bounds the unconditioned pair energy by `27r/4`; an upper-energy cutoff
`54r/k^2` retains at least `k^2/8` original weight. Zero-energy components are
not excluded by this pair bound. For arbitrary fixed n the energy bound is
`12r(1-2^-n)^2`. Every integer internal-charge sector is therefore realized
and has energy infimum zero, allowing a zero-energy kernel. These labels
do not count species or prove massless particles. Separation concerns the
centres of Schwartz profiles, not disjoint compact operator supports.
No clustering, independent-particle factorization or scattering state is assumed.

On the selected finite cube, the same
[two-charge response](../../model/source_charge_scattering/PROOF.md)
spans each entire 28-dimensional charge+-2 sector. Complete diagonal
projectors and a certified nonzero Hamiltonian walk instantiate the finite
matrix-unit argument. A complete 29-effect detector observes every vector
in that sector; together with the four preparation records it retains116
coarse outcomes. Exact connected work and retarded response verify source
interaction and propagation. Local complete histories have a fixed-spacing
infinite-volume limit on the selected state, with volume taken to infinity
before the selection field vanishes. An asymptotic scattering limit needs
an actual comparison dynamics and an integrable long-time error estimate.

### Local charged inputs and positive multipacket energies

The [configuration theorem](../../model/spin_multipacket_comparison/PROOF.md)
uses the existing P2 Weyl shifts and the configuration U1 to prove that
every finite independent Q3/Q2 vacuum configuration law is absolutely
continuous. It does not assume a second canonical momentum. Nonzero
polynomials in distinct bounded charged configuration creators therefore
have strictly positive vacuum norm. This gives nondegenerate finite Gram
matrices without Gaussianity, clustering or a tensor-product particle space.

The [local creator construction](../../model/spin_multipacket_comparison/LOCAL_CREATORS.md)
cuts off an annular profile in physical space and proves its response
directly from the surviving canonical commutator. The operators have compact
spatial supports, exact charges +-1 and norms at most1/2. With
`d=3m0/1000`, `m0=1/10`, `a=3d^2/2` and `b=2/d^2`, each original joint
event in `r/2<=|P|<=r`, `ar<=L<=br` has probability at least
`9/200000000`. No small unknown probability is divided out.

A smooth positive-energy filter of the same clock gives operators C_+ and
C_- with energy transfer in `[ar/2,2br]`. Analytic Fourier bounds supply a
scale-independent normalization M, retaining probability
`p=9/(200000000 M^2)` in the inner band. Both are raising operators; they
are exchanged by charge conjugation and are not each other's adjoints.
Their positive third effect retains a complete three-record instrument.
The inherited cone gives explicit almost-local commutator and derivative
tails for these operators on every input and reference.

Finite Jensen averaging plus those commutator bounds gives separated
equal-sign n-fold preparations of original weight at least `p^n/4`,
with all `3^n` records. The selected n-charge branch has its entire energy
in `[nar/2,2nbr]`; this interval is not imposed on the complementary records.
Thus every nonzero integer internal-charge sector has positive spectrum
approaching zero. This statement excludes a kernel-only explanation of
that particular construction; it neither counts elementary species nor
assigns rest masses. The filtered proof applies to equal-sign strings;
the unfiltered commuting construction above has the separate all-sign scope.

The [comparison estimate](../../model/spin_multipacket_comparison/POSITIVE_MULTIPACKETS.md)
isolates the scattering dependency. For a candidate dispersion w(k), the
one-packet derivative norm squared is exactly
`integral |E-w(k)|^2 |fhat(k)|^2 dmu_C(E,k)`. It is independent of time.
A fixed broad band alone cannot make this error integrable on an infinite
time tail. Actual dispersion control and comparison-space normalization
are required in addition to the proved inter-packet locality bounds.
The [source dispersion diagnostic](../../model/spin_multipacket_comparison/COMPARISON_DEFECT.md)
identifies the least scalar comparison error as the conditional energy
variance at fixed momentum. Narrowing a window onto a graph of zero
source weight loses its original success probability. On the selected
cube the complete four-character cube-symmetry readout has strictly positive
optimal defects for both local charge preparations. This concerns the
whole chosen prepared vector; it does not exclude a spectral subbranch.

An exact branch-selection criterion uses nested momentum cells K_l and
energy bins I_j: `A_j,l=sum_K max_I mu(I x K)`. Refining momentum first
and energy second gives the integral of the largest conditional energy
atom at each momentum. A positive limit therefore selects a genuine
reducing graph projection, even beside a continuum; repeated selection
resolves all conditional atoms. The finite statistic is stable under
complete probability error. No positive limiting value, smooth velocity,
uniform branch residue or effective source/refinement error follows from
this criterion. Its sharp continuum projectors need not be local or have
a proved native implementation.

### The energy-flux and boost test

The [source stress construction](../../model/spin_boost_stress/PROOF.md)
defines `B_i=sum_x r_i(x)e_x` and `P_i=i[H,B_i]` from the actual local
energy densities. Exact local momentum balances retain their force term,
and coarse energy blocks retain their internal dipole correction. This
`P_i` is an energy-flux candidate; equality with the existing spatial
translation generator is an additional theorem. On the selected finite
cube the candidate boost fails vacuum invariance and its Lorentz-algebra
defect `D_i=i[B_i,P_i]+L` has expectation less than `-0.13115434068`.
A complete native pointer detects the vacuum change. Finite rapidity
errors are controlled by integrals of that same defect, so taking many
small steps does not remove it. A relativistic continuum construction
requires the relevant scaled defect to vanish on actual low-energy field
histories, together with the momentum and coarse-dipole identifications.

The [weighted source energy](../../model/spin_energy_polarization/PROOF.md)
has an exactly cutoff-independent change on compact preparation supports
and an affine spatial translation law on the same clock. Its Q and P
Weyl-kick expectations converge to the spatial first moments of
`6 kappa F^2` and `kappa |grad G|^2/2`, respectively, with explicit mesh
errors and Schwartz-tail bounds. These are microscopic expectation limits;
energy-moment tightness is required to identify them with continuum boost
operators. A finite energy first moment need not survive a weak limit.

For the ordered-zero-field prescription, the
[actual quadratic twist](../../model/spin_twist_closure/PROOF.md)
has operator-norm remainder O(p^2). Its quadratic bond coefficients become
central, time-invariant variables in an enlarged source-history
representation. Their sum is scalar by the ground-energy inequality;
individual scalarity and closure of the original canonical algebra do not
follow. Exact source forces and a positive oscillator comparison retain
this distinction.
The [selected-ground energy bound](../../model/spin_twist_closure/STATE_CONTROL.md)
gives `||(H_Lambda/|Lambda|-e0) Omega|| <= 9J/L` on a cube of side L,
with complete regional energy-tail bounds. It does not control every
regional observable: opposite ordered source phases have identical energy
statistics but one-site half trace distance greater than `11/50`.
The actual P-field drift contains a different current average amplified
by `p^-2`. Qualitative mean-ergodic cancellation of that average supplies
no force closure without the corresponding quantitative correlation bound.

For the candidate energy polarization B_i, the exact source identity is
`[B_i,Q_p(F)]-[L_p,Q_p(X_i F)]=i R_i,p(F)`. The
[current remainder theorem](../../model/spin_twist_closure/POLARIZATION_LIMIT.md)
proves `R_i,p(F) W_p Omega -> 0` for every fixed finite bounded canonical
time-history word. Its derivative test has integral zero, removing the
entire mean-ergodic current component without assuming ergodicity.
The source phase-transport bound extends this cancellation through time.
An actual product-source control excludes operator-norm vanishing on all
states. The theorem concerns the difference of the two commutators; it
does not establish their separate unbounded limits, the P-field action,
vacuum boost invariance or a Lorentz-covariant particle shell.

### Existing theorem connections

The following results remove repeated construction work while preserving
their actual source hypotheses.

| Reused result | What it supplies to this chain | Exact connection required |
| --- | --- | --- |
| [Canonical field](../../model/spin_canonical_field/README.md) and [energy](../../model/spin_canonical_energy/README.md) | One selected-state history representation, strong field regularity, positive energy, locality and nonzero annulus weight | Ground-specific spectral/force estimates determine dispersion; the representation need not be rebuilt |
| [Selected-state colimit](../../Lean/QFT/ColimitSelectedStateGNS.lean) | A representation from an exactly compatible observer-state family | The family and its ground-state attachment are premises; no effective selected-spin marginal convergence rate follows |
| [Finite source correlations](../../Lean/QFT/SourceCorrelationCapstone.lean) and [net adequacy](../../Lean/QFT/StructuralNetAdequacySurface.lean) | Finite joint-state and structural net implications | Their stated hypotheses do not supply the ordered-spin continuum vacuum, time-slice attachment or sharp spectrum |
| [Finite operator generation](../../Lean/QFT/SourceOperatorGeneration.lean) and [native compiler](../../model/spin_operational_position/PROOF.md) | Matrix-unit generation from diagonal projectors and a connected transition walk; complete finite observer instruments | Supply the actual sector projectors, nonzero source edges and admitted coherent controls |
| [Celestial null rays](../../Lean/Geometry/CelestialNullCone.lean) and [native transition test](../../model/native_lorentz_transition/PROOF.md) | Sphere labels for future-null rays; twelve exact port tests determine an invertible conformal Lorentz map | Actual field/clock transitions must realize the map and its normalization |
| [One-boost cone theorem](../../paper/tex_fragments/OPERATIONAL_CAUSAL_SELECTION.tex) | Finite irreducible rotations plus one operational boost direction determine the invariant causal cone | Boosts must transform admissible histories and readouts; a spherical coordinate action alone does not do this |
| [Minimal Einstein composition](../../Lean/Geometry/SourceOrderEinsteinComposition.lean) | Nine source-direction balances suffice for the finite Einstein-form tensor conclusion | Supply the same physical stress, responsive geometry, balances and normalization; no second tensor derivation is needed |

The finite source action `PSL(2,F5)`, the conformal group `PSL(2,C)` of the
celestial sphere, and the internal spin U1 are different groups. Their
representations are reused only through an explicit map. The sphere is
equivalent to the null-ray space, not to the Lorentz group itself.

The spectral hypotheses remain essential in related analytic results.
[Koma's Goldstone theorem](https://arxiv.org/html/2105.04970v2)
bounds filtered packet energies between positive multiples of momentum.
[Michoel--Verbeure, Theorem 9](https://arxiv.org/html/math-ph/0001033)
obtains closed oscillator equations under its single-peak assumption (22),
with additional fluctuation-correlation hypotheses. Neither statement
supplies the required sharp branch on the joint B(p) family. Unscaled
concentration at zero is insufficient: `(delta_p+delta_(2p))/2` converges
to delta_0, while its E/p distribution retains two peaks. The
[two-speed response control](../../model/spin_dynamical_memory/RESOLVENT_CLOSURE.md)
also preserves both canonical pairs and internal U1. These are implication
boundaries, not a proof that the actual source lacks a sharp branch.

### Position, translations and uncertainty

Uniform native repair has mean load transition `I-L/60`. Its smallest
nonzero Laplacian band, `5-sqrt(5)`, has rank three. The
[record-recovered Gram](../../model/relative_observer_geometry/README.md)
`G=4 P_slow` supplies unit frame directions `u_a` up to orthogonal isometry.
For a fixed-total load `n`, define the dimensionless carrier readout

```
x(n)=sum_a n_a u_a,
x(n_after)-x(n_before)=j (u_b-u_a).
```

The second identity is an actual record readout: the seam and committed current
supply each displacement without hidden source loads. Absolute position needs
an initial reference. This three-component readout is a local load-weighted coordinate,
not an injective position label for every many-body state. The calibrated
structural Gram survives the stated exposure correction; an arbitrary biased
scheduler need not retain it as the dynamical slow band. The
[real frame quotient](../../Lean/Screen/PrimitivePortFrameQuotient.lean)
completes its signed control module to Euclidean three-space. It does not
identify a bounded fixed-total trajectory with that full completion or supply
inter-carrier frame gluing.

The [same-readout operator proof](../../model/carrier_position_translation/PROOF.md)
strengthens the recorded displacement identity. For every complete fine branch
`K_h`, including arbitrary source/reference inputs,

```
[X_i,K_h]=d_i(h) K_h,
exp(i q.X) K_h exp(-i q.X)=exp(i q.d(h)) K_h.
```

The exact certificate uses the unnormalized algebraic frame; dividing every
coordinate and displacement by its common circumradius gives the unit frame
above. Dividing the load coordinate by the fixed total gives the centroid.

Each branch is a filtered partial shift. Its acceptance depends on the load;
it is not a unitary translation. For a finite-dimensional position observable,
any unitary obeys `||U*X_iU-X_i-a_i I||>=|a_i|`. The uniform repair average
contracts the slow position by `1-(5-sqrt(5))/60`, rather than translating it.
These exact statements do not prevent a controlled limit on expanding states,
but exclude the literal full-register translation identity.

The protected spin source has a stronger restriction. Its subspace `P` has
one fermion per site, so `N_x P=P`. Every total-load coordinate
`X_f=sum_x f_x N_x` restricts to `c_f P`, where `c_f=sum_x f_x`.
For the source low-band isometry `W` with `||W-P||<=2 tau/U`,

```
||W*X_f W-c_f I|| <= 4 (sum_x |f_x|) (tau/U)^2.
```

This follows by expanding the centered observable, which annihilates `P`
on both sides. Fixed local load coordinates become scalar on the declared
strong-coupling trajectory. They cannot locate an internal spin excitation.
The source's carrier frame and an excitation's spatial detector have distinct
roles; recovering a rank-three frame does not identify their observables.

In the ordered-spin continuum, translations of test functions define the
actual commuting momentum generators. The
[field-generated position proof](../../model/spin_covariant_localization/CONTINUUM_POSITION.md)
constructs a conjugate coordinate on the closed subspace
`H_Q=closure{Q(F)Omega}` of this same representation. The field uncertainty
relation and the existing upper fluctuation bounds give

```
m^2 integral |Fhat(k)|^2/|k| <= ||Q(F)Omega||^2
                            <= 3 integral |Fhat(k)|^2/|k|.
```

With Fourier measure `d^3k/(2pi)^3`, set `T(Fhat/sqrt(|k|))=Q(F)Omega`.
The source Gram operator satisfies `m^2 I<=T*T<=3I` and commutes with
spatial translations. Its polar factor `V=T(T*T)^(-1/2)` is a unitary onto
`H_Q` intertwining those same translations with multiplication by
`exp(-i k.a)`. Therefore

```
X_i=V i partial_i V*,  P_j=V hbar k_j V*,
U(a)*X_iU(a)=X_i+a_i I,
[X_i,P_j]=i hbar delta_ij I  on V Schwartz(R^3),
Delta X_i Delta P_i>=hbar/2.
```

Here `hbar` fixes an action convention. The proof calculates no laboratory
value for it. A complete extension to the entire field space retains the
no-position outcome `I-Pi_Q`; the invariant vacuum gives that outcome with
certainty. A normalized covariant position measure on the full space without
this outcome is impossible, since the vacuum would induce a translation-
invariant probability measure on all of Euclidean space.

This proves a source-defined position/momentum pair, with explicit operator
domains. Its Gram whitening is generally nonlocal. The
[finite source compiler](../../model/spin_operational_position/PROOF.md)
implements computable finite instruments using the declared coherent load,
pointer and refinable six-level controls. The eight-site coordinate detector
executes a 49,228-pulse decoded core with an untouched reference. The full
49,324-pulse encoded detector has an analytic half-diamond error below
`1.53e-10`; execution of that entire encoded tape is not claimed. Every
pointer label and export failure is retained. These services are hypotheses
of the executable observer and do not follow from bare occupation repair.

The [finite-mode position proof](../../model/spin_position_limit/PROOF.md)
uses the same source Gram and bounded sine probes. Finite positive Gram
compressions converge strongly to the polar map, and their complete
Lueders instruments converge pointwise on fixed states and references.
Computable spatial bins and certified finite matrices admit source
compilation. The limit orders finite-source and field refinement, probe
linearization and mode exhaustion; it gives no effective box-size modulus,
infinite-dimensional uniform channel bound or second-moment tail control.
The coordinate is distinct from the total-load centroid.

The [dynamical extension](../../model/spin_dynamical_localization/PROOF.md)
uses `H_resp=closure span{exp(-itL)Vf}` under the same positive generator.
Its history Gram gives a translation-intertwining isometry into
`L2(R3;ell2)` and a compressed covariant position measure. Both this response
space and its dark complement reduce time evolution. For normalized states in `H_resp`, the position-outcome and momentum
standard deviations obey `Delta X Delta P>=hbar/2` when both
moments exist. Time invariance of `H_Q`, moment persistence and a sharp
particle branch are not inferred. The full response-space measure has no
proved operational continuum implementation. The finite cube's invariant
odd-parity enclosure has eight rank-16 coordinate effects and an even dark
complement; it is distinct from the executed rank-eight detector. Its finite
symmetry labels are reflections, not continuous momenta.

The [finite Robertson theorem](../../Lean/EventAlgebra/Robertson.lean)
also applies directly to the actual finite observables. Finite carrier
uncertainty uses its state-dependent commutator, not a nonzero canonical
matrix commutator. Knowledge of the simulator specification does not make
incompatible instruments jointly sharp.

## Observable charge, energy and gravitational response

[Complete local charge histories](../../model/native_local_charge_history/PROOF.md)
attach every decoded occupation snapshot to the actual incident currents.
The fine-branch operator continuity identity retains arbitrary references,
boundary flux and erasures. One complete-instrument error bounds the joint
history; it does not accumulate once per snapshot. This supplies occupation
transport for CP-05, without assigning an electromagnetic charge or particle.

[Positive transfer-energy windows](../../model/native_energy_identification/PROOF.md)
give a complete spectral partition from actual lag records. The finite fair
total-two example has seventeen bins, including overflow, and less than
`4e-11` preparation error in the whole histogram. These are signed statistical
reconstructions of positive probabilities. The sufficient one-percent
sampling budget exceeds `3.2e25` repetitions. The distinct
[controlled energy instrument](../../model/native_energy_measurement/PROOF.md)
avoids this signed-moment amplification by composing sixty-six existing
destructive phase filters. It retains every energy label, processor failure
and final residual, with complete half-diamond error below `1/12000` under
its primitive certificates. All 66 words are individually admitted; their
3,853,872 pulses are a sequential service composition, not a monolithic
repair program. The generic expanded host retains its resource stop.

The [fair regional entropy attachment](../../model/native_gravity_identification/PROOF.md)
uses the actual canonical source, every local number sector, and the finite
entropy theorem. It gives a full-state entropy maximum and a regional first
law. Its modular generator is central on the number-preserving local algebra;
native repair changes observables on that same algebra. Sector reweighting
also changes edge entropy while a fixed graph cut does not. Thus the direct
attempt-clock/modular-clock and fixed-cut/area identifications fail on this
specified dictionary. The [record-filtered reference](../../model/native_noncentral_reference/PROOF.md)
provides a constructive alternative in 71 native attempts. Its conditional
acceptance is at least `1/12` on a recognized nonconstant sector, and its
accepted spectral floor is at least `1/8192`. The actual nonconstant weight
can be zero; rejects and the coherent constant corner remain explicit.
Its natural regional modular actions agree on overlaps. They generate
local charge phases, fix occupations, and supply no geometric boosts.

The total-two reference feeds the controlled energy instrument through the
same native configuration chart. Exact rational products give the complete
law with rejection `57765/65536`, accepted weight `7771/65536`, all 66
energy outcomes and final residual. The accepted reference observation
traces the sixty mixing-record values; the actual archive retains them.
This is a sound complete finite readout path, without a particle assignment.

The existing ideal phase-experiment reference `rho=P/56` gives a separate
exact state/operator-clock join: `-log rho=H/60+log(56)I`. Its owning
preparation retains the finite selector timeout and numerical-density
admission scope; no outward-rounded preparation channel certificate is
imported. This state differs from the diagonal tree-weighted reference by
half trace distance at least `1/336`, witnessed by the actual norm-one
occupation exchange. The two preparations' regional modular and energy-clock
properties cannot be combined as properties of one common state.

For a physical continuation, fix one observer-visible regional algebra,
translation and refinement dictionary, and prove its complete experimental
law comparison with a vanishing error. Its clock must be the clock in the
propagator. Within the occupation dictionary, its two-dimensional continuum cannot
become a three-dimensional spatial generator by an invertible same-generator
renaming. The separate three-dimensional field constructions above retain
their own stated source and state attachments. The spherical null-direction theorem constrains the Lorentz form
once the spatial readout and clock are supplied; it does not supply those
maps or the action. A different field dictionary must be constructed and
tested, with source interventions and controller dependencies retained.

The noncentral reference's state-preserving regional expectations are
completely positive, generally nonmultiplicative. They cannot fill the
homomorphic restriction field in `FiniteCausalObserverNet`: for an exterior
projection with conditional weight `0<p<1`, `E(P^2)=pI` but `E(P)^2=p^2I`.
Use a local algebra net with its actual inclusion maps, or prove an applicable
weaker restriction theorem. Compatible modular phases alone do not produce
causal translations, an interacting chiral limit, physical rest masses or
responsive Einstein geometry.

These results retain the physical requirements below. Reading the simulator
specifications fixes the source law and operational menu; it does not remove
the need to prove which state, clock and observable maps describe one common
physical field.

## Common-source clocks, energy and reference tests

The [local energy construction](../../model/native_local_energy_clock/PROOF.md)
uses the actual binary federation, with `V=240*4^ell` ports, `M=630*4^ell`
seams, `G=sum_e(I-S_e)` and `P=I-G/(2M)`. Its triangle cover proves the
capacity-independent floor `P>=10/21`. For `rho=P/Tr(P)`, the logarithmic
clock `H=-2M log P` approaches `aG` on fixed local observables in the
rho-weighted norm. Fixed complete classical histories obey a separate
total-variation bound of order `O(M^(-1/2))`; `a` tends to
`1/[1-p(1-p)]` at density `p`.
This comparison does not cover arbitrary input states or growing archives.
The antisymmetric currents `i[I-S_e,I-S_f]` conserve `G` exactly.
Native unread repair instead exchanges energy: on the local total-two
reference its mean `G` increase is `37/840`. A total conserved energy
requires an apparatus-energy observable.

The [regional construction](../../model/native_regional_modular_clock/PROOF.md)
takes reductions of this same `rho`. On local number-preserving observables,
its scaled regional modular action approaches
`[G_A+(1-2p)N_boundary]/[1-p(1-p)]`. A fixed nested-cut control has mismatch
`4/3`; a [collar estimate](../../model/native_regional_modular_clock/COLLAR.md)
suppresses the effect of the boundary on interior observables. The literal
348-port window has error below `1.068e-6` at time `0.01`, density `1/2`.
Joint growing-window control requires `v^2/M -> 0` and
`v*abs(q/V-p) -> 0`, compatible port routing and separately bounded
preparation, service and readout errors. This is graph propagation, with
no identification of a Lorentz cone or geometric modular boost.

The [fair ground experiment](../../model/native_fair_ground_readout/PROOF.md)
uses a different, pure reference. Its finite preparation error is below
`1/4000000`; an occupation phase creates exchange energy `50/33` and a
complete measured mean-score contrast above `3/50`. All preparation and
detector failures remain in that unnormalized comparison. Ground-local
clock words have error `O(M^(-1))`, but scalable preparation is unproved.
The mixed reference has negative commutator frequencies and cannot be
the nonstatic positive-energy vacuum. The pure fixed-number ground has
nonseparating local algebras if its nonzero bounded exchange observables
survive the limit. These are state-specific constraints, not a prohibition
on every possible field dictionary.

Exact [coherent occupation controls](../../model/native_coherent_occupation_control/PROOF.md)
combine the declared source generators before exponentiation. They preserve
number at literal dyadic angles and give a bounded addressed service without
changing the frozen interpreter. The [complete correlation experiment](../../model/native_source_chsh/PROOF.md)
has CHSH score above `12/5` with pre-setting aborts and every native record
retained. Coherent access, initial rails, settings and controller remain
explicit. These services preserve cell-charge block structure under their
stated grammar; encoded entanglement alone does not prepare the growing
cross-cell coherent reference. Neither the finite correlation witness nor
the clock estimates selects a chiral relativistic action.

The direct continuation must choose one prepared state and one scaling
of the local algebra, retain its complete controller/readout law, and prove
a nontrivial positive-energy relativistic limit. Only then can exchange
currents be compared to a conserved symmetric stress tensor and regional
entropy to a metric response. Particle sectors must be identified in that
same limit, rather than by matching finite eigenvalues to measured masses.

The [dilute-ground theorem](../../model/native_dilute_ground_limit/PROOF.md)
tests that limit directly. For a bounded number-preserving observable on
`s` ports, its ground variance is at most `4*norm(A)^2*s*q/V`.
Fixed complete local instrument histories differ from their scalar
empty-sector histories by at most `2*sum_j sqrt(s_j*q/V)` in complete
half trace distance. Thus the fixed-`q=2`, fixed-support, bounded-observable
route has vanishing connected signal at every rescaled time. A positive
density, growing observation regions, renormalized fields or an admitted
charged algebra changes the hypotheses and needs a separate construction.
The exact density branch inherits graph-Laplacian frequencies; a literal
cycle control has low-mode energy ratio tending to four rather than two.
The selected federation's diffusive spatial conclusion retains its
spatial-scaling premises. None of these statements excludes every
possible encoding or a hydrodynamic sector within a relativistic theory.
In particular, retaining spectral weight at least `w` with operator norm
at most `K` requires `s*q/V >= w/(4*K^2)`. A candidate field dictionary
must meet this visibility condition before its pole locations can support
a particle assignment.

The [centered ground construction](../../model/source_family_compatibility/GROUND_MARGINALS.md)
gives a finite algebraic preparation procedure and computable pressure-secant
bounds for local marginals. The symmetric default reference has all four
one-site eigenvalues above `1/90000`. Six native tagging pulses on the
prepared full cube retain 36 pointer pairs, 84 additional first-failure
records and the inherited preparation failures. Its one-site half-trace
comparison to every stated symmetric thermodynamic ground is below `13/50`.
This coarse one-site bound is insufficient for small spectral windows;
multisite marginal convergence with the required accuracy is unproved.

The [fixed-source scaling theorem](../../model/source_family_compatibility/MASS_SCALE_SELECTION.md)
separates bounded-window survival from a massive spectrum. For the original
CAR readout on one fixed normalized reference and coupling ratio, with
the stated clock and observable dictionary, a full ballistic joint limit is dilation covariant. It assigns zero original
weight to every nonzero relativistic mass shell. Convergence along the
whole geometric refinement gives discrete covariance: any positive-mass shell forces all
its discrete mass dilates. Neither statement supplies that limit or excludes
every resonance, changing observable dictionary, or controlled finite-cutoff
reconstruction. A proposed massive identification must derive the actual
capacity dependence of its source, state, clock and readout from the
implementation, rather than choosing them to match a target mass.

<a id="source-to-field-requirements-and-status"></a>

## Native repair and ordered-spin comparison domains

This inventory preserves theorems for the native repair and ordered-spin
programs. It is evidence with explicit source boundaries, not another
canonical dependency list. The selected Wilson program and conditional
chiral extension use the [unified graph](../../model/particle_joint_reconstruction/DEPENDENCIES.md).

The [finite action identification](../../model/source_family_compatibility/SOURCE_ACTION_IDENTIFICATION.md)
recovers the configured positive `U/c` from complete two-time collision
records on a known dimer preparation and clock. It supplies no law selecting
the extension, its capacity trajectory or the interacting ground. Source
specification access permits reading these inputs; a physical derivation
must connect that specification to the common reconstructed observations.
The [volume-uniform cubic inverse](../../model/source_family_compatibility/BULK_ACTION_IDENTIFICATION.md)
derives the local collision coefficient `c^2(U^2+25c^2)` from exact signed
paths on every interior edge. Its finite empty-collar estimate allows
arbitrary exterior matter for that statistic. It consumes a known local
calibration preparation, with quantitative marginal and clock errors;
it does not substitute that preparation for the interacting ground.

The decisive spectral estimates differ between the two prepared source
families. For the default centered CAR family, write `K=c K(1,8)` and
take genuine momentum packets near `k_*+a p`, with `c=v/a`, where `k_*`
is a zero of the actual magnetic symbol `h(k)=sum_j cos(k_j) Gamma_j`.
The original Weyl-chart node labels require the corresponding transport
of the probes and Bloch convention. A positive
uniform lower bound on the original Lorentzian score
`q_a=-eta_a Im<f_a,G_00(i eta_a)f_a>`, `eta_a` proportional to a,
would establish retained finite physical-energy weight: for fixed M,
`nu_a(|e|<=M eta_a)>=q_a-1/(1+M^2)`. The
[finite source filter](../../model/source_family_compatibility/CUBE_MEMORY.md)
supplies its complete record meaning. Its two fixed-cube frequencies do
not prove this growing-volume, shrinking-energy estimate. The actual
low-frequency memory, marginal and detector errors must be controlled on
that same family before extracting a momentum-resolved pole or threshold.

For the joint positive-field collective family, the
[source Ward defects](../../model/spin_joint_field_limit/WARD_DYNAMICS.md)
must bound leakage outside a shrinking energy--momentum tube, using the
actual filtered cross moment and both inverse susceptibilities. The
[bounded preparation inequality](../../model/spin_joint_field_limit/ARCTAN_CHARGED_BANDS.md)
then transfers this leakage to an original charge-probability floor without
a growing preparation normalization. The [two-site source inequality](../../model/spin_joint_field_limit/QUANTUM_FORCE_SLACK.md)
proves `liminf(kappa_B-m_B^2)>=61/14400` on this same joint trajectory.
Exact double commutators therefore exclude closing the unfiltered forces
by replacing every longitudinal spin by its staggered mean. A corrected
infrared coefficient or sharp subbranch is not excluded. This bound does
not discharge the filtered susceptibility/slack requirement. A common linear low-momentum shell
would identify a massless branch only; additional massive sectors require
their own surviving weights and dispersions. Neither route transfers its
state or clock to the other family. Relativistic transformations, surviving
interactions and the empirical quantum-number dictionary remain separate
requirements for the full spectrum.

The component [defect-identification criteria](../../model/particle_identification/ROADMAP.md#defect-to-particle-identification-requirements)
expand the creation, protection and observable-sector tests within unified
nodes 6-8. The interfaces below retain the detailed comparison constructions.
The eight free Weyl valleys are a model control, without an assignment to
eight measured particle species.

The table below concerns the native repair and ordered-spin source-to-field
routes. The gauge-coordinate causal, thermal and coframe results above have
their own source and state domains; transferring them into these routes
requires an explicit comparison map. The table concerns physical
identification rather than whether a finite test suite passes. The
[claim registry](../../claims/claim_registry.yaml),
[premise register](../../tracking/premise_register.json) and
[physical identification register](../../claims/physical_identification_registry.json)
retain the exact hypotheses. The 84 premise entries include alternative
theorem branches; they are not 84 jointly required laws. The 21 physical
identifications are undischarged. No final measured particle-mass ratio is
derived by the constructions below.

| Comparison domain | Established result on its stated source | Boundary of a transfer to the selected family |
| --- | --- | --- |
| **Native executable families** | **Static restrictions, exact corner cell topology and finite routed local-history comparisons established; physical capacity family unproved.** The native program declares its repair, memory, scheduler and observer services. [Canonical restriction](../../model/native_capacity_restriction/PROOF.md) gives every local number sector, exact nesting and finite-size bounds for the fair federation. The [weighted local limit](../../model/weighted_local_limit/PROOF.md) controls full local number mixtures for two declared densities of the weighted extension. [Dynamic product restriction](../../model/native_capacity_restriction/PRODUCT.md) preserves complete biased repair/hold records under unused-factor extension. The [local-history limit](../../model/native_local_history_limit/PROOF.md) controls complete fair records with time bins under matching rooted-window and state hypotheses. [Rooted refinement](../../model/native_rooted_refinement/PROOF.md) supplies a full-sequence corner cell topology and certified finite routed windows, with an explicit arithmetic-continuation boundary. | Select one family indexed by physical capacity, including its conserved occupation, preparation, controls and clock. The fair federation is not a proved capacity continuation of the weighted H66 observer; the weighted connected extension is not admitted by that interpreter. Port count, FIFO slots and cosmological capacity are distinct inputs. |
| **Native observable access** | **All-sector local preparation, complete occupation readout, exact dyadic coherent controls and a complete CHSH witness established; universal field-observable access unproved.** The [versioned control service](../../model/native_coherent_occupation_control/PROOF.md) and [complete witness](../../model/native_source_chsh/PROOF.md) retain coherent-access, state, setting and controller premises. The [fair local construction](#local-preparation-complete-readout-and-a-shared-attempt-clock) supplies controlled preparation and a fresh detector on one attempt counter. Declared global occupation resolves infeasible constant alternatives; general-density readout keeps bounded explicit erasures. The [full-number E-tree detector](../../model/native_sector_readout/PROOF.md) reads any declared global sector with complete back-action and a retained empty/full ambiguity when the local number is unknown. [Complete E records](../../model/native_record_response/PROOF.md) give rank-66 population observability. The [coherent clock protocol](../../model/native_clock_observables/PROOF.md) adds complete Pauli instruments and finite Robertson inputs using existing declared services. The [tree detector](../../model/native_effect_coverage/PROOF.md) supplies configuration readout, basis preparation and the inherited tomographic menu. [Rank-one synthesis](../../model/native_effect_synthesis/PROOF.md) constructs bounded rational targets, with complete error below `10^-6` for words passing the exact primitive gates. [Coexistent weak measurements](../../model/native_coexistent_effects/PROOF.md) compose up to eight source filters into one complete noncommuting menu with explicit state changes and retained failures. [General finite spin instruments](../../model/spin_operational_position/PROOF.md) are constructively compiled with arbitrary precision under the declared refinable coherent service; [finite-mode position limits](../../model/spin_position_limit/PROOF.md) retain their ordered convergence boundary. Direct compilation and lazy histories retain their distinct host execution limits. | Extend the proved local mixed-reference and occupation instruments to the required interacting field state and complete field observables. Establish the operational equivalences and effect coverage required by the quantum theorems. The supplied local occupation-wire coupling implements a bounded phase experiment; neither it nor the coexistent family implements arbitrary effects, every coexistent pair or an arbitrary-precision coherent refinement. |
| **Native spatial and clock comparisons** | **Native graph, shared attempt clock and state-specific local Hamiltonian comparisons established; physical spatial limit unproved.** The [local energy theorem](../../model/native_local_energy_clock/PROOF.md) gives a fixed-history comparison for the faithful transfer reference, and [regional collars](../../model/native_regional_modular_clock/COLLAR.md) control boundary-flow errors. [Clock composition](../../model/native_fair_local_clock/PROOF.md) counts fixed local operations and uniform repairs in one program, and bounds the reservation cost in its separate binned comparison. [Spatial scaling](../../model/native_spatial_scaling/PROOF.md) attaches the source routing and affine-coordinate bounds to the actual counts, with distinct global-attempt and incident-clock estimates. The [exchange clock](../../model/native_exchange_clock/PROOF.md) controls the endpoint comparison to the native local action at counted time `T/(2M)`. The complete local-history comparison uses per-seam attempt time `T/M`, aligned time bins and a stated record trace. [Authenticated dataflow](../../model/native_causal_dataflow/README.md) retains the full single-FIFO chain. [Independent repairs](../../model/native_record_causality/PROOF.md) have equal complete channels after a specified reversible record permutation, for fixed completed schedules. [Lorentz transition tests](../../model/native_lorentz_transition/README.md) constrain supplied transitions. | Derive spatial readouts, translations, causal cones and clock normalization on one refinement family. The reconstructed transfer generator and working modular generator are inequivalent. Their proposed same-reference identification has a [complete-output separation bound](../../model/native_clock_comparison/PROOF.md). A physical quotient must preserve observable controller signals; adaptive selectors and timeouts are outside the independent-repair equivalence. The celestial null-ray parametrization alone supplies neither a three-dimensional field net nor this limit. The [centered-source dipole identity](../../model/source_family_compatibility/ENERGY_DIPOLE.md) isolates the interaction-current contribution to a possible geometric boost. |
| **Native reference preparations** | **All-sector local mixed preparation and finite fair and weighted transfer-ground experiments established; common interacting vacuum unproved.** The [fair ground construction](../../model/native_fair_ground_readout/PROOF.md) certifies its preparation, local energy contrast and complete readout; growing ground preparation and a physical vacuum algebra are separate premises. [Fair local preparation](../../model/native_fair_local_preparation/PROOF.md) gives half-diamond error `2^-b` after `11+60b` attempts, with actual sector weights, arbitrary references and a stated archive marginal. It does not select the global canonical input or a pure ground. [Observable warmup](../../model/native_correlation_preparation/PROOF.md) controls fixed-`q` reflected scores from arbitrary quantum input, without claiming full-state convergence. The [repair-only reset bound](../../model/native_sector_readout/PROOF.md#a-bounded-repair-only-reset-obstruction) excludes exact finite deterministic reset of all inputs on a nontrivial sector to one basis target; approximate or separately controlled preparation needs its own theorem. [Transfer-ground reconstruction](../../model/native_vacuum_response/PROOF.md) recovers positive-energy correlations as finite combinations of complete lag experiments; its phase-intervention witness excludes a direct identification with the working tracial state. The [controlled transfer experiment](../../model/native_transfer_control/PROOF.md) prepares and evolves the weighted finite transfer ground using declared coherent services, with complete failures and a certified response. It does not provide the growing fair-source vacuum. The [protected native reference](../../model/protected_reference/PROOF.md) classifies stationary states, bounds sharp-charge preparation and joins actual repair to its modular reference. Its [gauge-singlet corner](../../model/protected_reference/GAUGE_JOIN.md) has exact replacement preparation and a Lean finite-thermodynamic adapter. | Prepare the state used by all field sectors on the same source and clock, with positive-energy and refinement control. The centered finite full-ground projector has an algebraic preparation procedure; actual pressure secants give local-marginal enclosures, with uniform one-site faithfulness at U/c=8, but no unique growing marginal or arbitrary-accuracy modulus is inferred ([source proof](../../model/source_family_compatibility/GROUND_MARGINALS.md)). A native stationary density, a chosen Higgs ground and a physical vacuum cannot be identified without a state-preserving map. |
| **Source response comparisons** | **Complete native spatial matrices, contact response, finite excitation-energy detection and local-history error bounds established; common relativistic response unproved.** The [full centered-cube readout](../../model/source_family_compatibility/CUBE_SPECTRUM.md) encloses complete finite-clock probabilities and true source-energy windows for eight prepared profiles, both charges and both spin components. Its open-box profiles are not translation eigenmomenta. The [composite residual](../../model/source_family_compatibility/CUBE_RESPONSE.md) is strictly positive in every original mode direction. The [local energy currents](../../model/native_local_energy_clock/PROOF.md) conserve the exchange generator and distinguish controlled evolution from native repair heating. The [same-occupation pair matrix](../../model/native_pair_response/PROOF.md) joins complete reflected cross probes, local-window errors and a collision-retaining exclusion comparison. The [native contact response](../../model/native_contact_response/PROOF.md) derives a visible primitive pole below its centered independent-pair comparison, together with a nonvanishing connected-cell leakage bound. [Connected contact dressing](../../model/native_connected_contact/PROOF.md) bounds its negative-parameter local resolvent and proves collapse of the global comparison threshold. The [reflected four-E probe](../../model/native_correlation_sectors/PROOF.md) isolates a degree-two correlation block of the fair connected transfer, with normalized response independent of occupation on the same graph. The [native record resolvent](../../model/native_record_response/PROOF.md) has 65 observable population-relaxation poles with actual back-action. [Correct record reflection](../../model/native_record_response/TRANSFER.md) derives a positive transfer Hamiltonian distinct from the native modular clock. [Local response](../../model/native_local_response/PROOF.md) reduces selected complete histories to their backward dependency cone and controls finite-window errors in global and internal repair counts. | Derive the common spatial response `G_N(z,p)` with controlled state, clock, record and numerical errors. The [all-volume nodal memory bound](../../model/source_family_compatibility/NODAL_MEMORY.md) forbids dropping the omitted composite dynamics near default kinetic nodes. The local comparison uses a stated observation map and fixed selectors; it does not remove observable controller signals. Attempt relaxation, internal repair counts and bounded Abel response require separate connections to relativistic energy and momentum. |
| **Comparison field limits and actions** | **A free relativistic branch and a separate three-dimensional ordered-spin canonical continuum are established on their stated sources; an interacting chiral relativistic theory is unproved.** The [uniform charged scaling bound](../../model/source_family_compatibility/STRONG_CHARGE_SCALING.md) excludes finite-energy original odd response at ballistic U/c tending to infinity, including a finite nonzero bare cell-contact coefficient. It does not exclude a critical intermediate trajectory or source-even collective sector. The [ordered-spin field](../../model/spin_canonical_field/README.md) has continuous translations, positive energy and a finite cone; the [same-source requirements](#the-three-field-constructions-and-their-source-boundaries) identify its additional correlation and state estimates. The [native density construction](../../model/physical_space_time/README.md) has spectral dimension two in its stated dictionary; it does not become a three-dimensional relativistic field by interpreting the transfer logarithm as a Hamiltonian. Bounded quadratic profiles have a conditional Wick limit within that diffusive Gaussian field, with moment control; this does not establish an additional interacting continuum. The [controlled continuum](../../model/continuum_free_limit/README.md) is a free eight-valley Weyl vacuum net on its stated scaling. [Literal gauge transport](../../model/gauge_defect_transport/PROOF.md) is complete on its declared extension. | Prove positive energy, locality, relativistic covariance and joint-record convergence for the same source. Derive its action and demonstrate that the required interactions survive. Gauge dynamics, chiral matter, anomaly conditions and cutoff control must hold on that construction; Lie-algebra classification does not select them. |
| **Comparison spectral sectors and ratios** | **Model sectors identified; measured massive spectrum unproved.** The finite native contact pole has a nonzero complete-record residue, but its unchanged isolated-cell vector fails the connected-capacity persistence test. The global comparison threshold is at most `550/V`; controlled off-spectrum dressing supplies no positive-energy persistence theorem. The [exact kinetic-node map](../../model/source_family_compatibility/KINETIC_NODE_MAP.md) attaches all eight windings to massless particle/hole sectors on the proved free family. The [joint prepared collective branch](../../model/spin_joint_field_limit/ARCTAN_CHARGED_BANDS.md) has bounded opposite-charge preparations with uniform original weight at least 1/387200 in broad energy--momentum bands; a sharp shell is not inferred. The [ordered-spin shell readout](../../model/spin_shell_readout/README.md) certifies finite spectral weight with full detector errors; its nonvanishing shrinking-tube condition is an unproved asymptotic premise. [Massive Dirac shells](../../model/source_mass_shell/README.md) use supplied onsite mass coefficients. Protected memory charges and finite response lines are distinct objects. | Identify stable poles, resonances, thresholds or bound sectors using creation, propagation, quantum numbers and complete detectors. Derive their energies from the selected action/state. Ratios cancel a common unit choice; they do not remove unspecified couplings, preparations or sector assignments. |
| **Source-specific empirical comparisons** | **Not established for the source-derived particle spectrum.** Conditional calibrations and numerical mass comparisons have separate provenance. | Freeze the observable dictionary, calibration inputs and error model, then compare predicted joint records, masses, mixing, widths and scattering where applicable. No target-informed choice may be reported as an independent prediction of that target. |

The component [construction criteria](../../model/particle_identification/ROADMAP.md)
state complete-channel and defect formation-to-readout tests. This table
catalogues their source-specific evidence. Calculations for different
actions, states or clocks require explicit comparison maps before supporting
the selected family's nodes in the unified graph.

### Connected readout, correlation sectors and spatial scaling

The [native sector detector](../../model/native_sector_readout/PROOF.md)
uses one E repair on each edge of an ordered spanning tree. Reading the
source's conserved global total `q` fixes the otherwise undetermined root
occupation. All `V-1` records then determine the initial and terminal
configuration. Each nonzero fine map is
`sqrt(w_h)*|y_h><x_h|`, with its arbitrary-reference output retained. Coarse
effects labelled by the decoded initial configuration are the exact basis
projectors. At the evaluated source levels 0,1,2 the detector has
239,959,3839 commands. A lazy raw-word query and the stock expanded-tree
`RESOURCE_STOP` are distinct execution results.

A local restriction generally does not inherit knowledge of its local
number. On the full binary mixture the same tree distinguishes every
nonconstant configuration, while the all-empty and all-full inputs share
one rank-two outcome, preserving their mutual coherence. The protocol does
not supply their phase measurement or exact repeatable restoration. Any
bounded repair-only controller has a positive no-motion path: for at most
`T` repairs its probability is at least `(1-r+r*p_min)^T`. This excludes an
exact deterministic reset of every input to a chosen different basis state.
The exclusion concerns a nontrivial sector, `0<q<V`; the empty and full
sectors are single states. It permits approximate, heralded or separately
controlled preparation.

A final tree can decode the terminal state of a body of repairs. Reversing
every retained signed current then reconstructs its earlier configurations,
without supplying an initial census. The caller must authenticate the body
program, selector probabilities and complete FIFO/archive. The final readout
changes the terminal state; it is part of the same observed experiment.

The [correlation-sector theorem](../../model/native_correlation_sectors/PROOF.md)
decomposes the actual fair population transfer into orthogonal primitive
orders `r=0,...,min(q,V-q)`, of dimensions
`binom(V,r)-binom(V,r-1)`. Adding unoccupied sites is an algebraic intertwiner,
not an admitted creation command. For disjoint native seams `(a,b)` and
`(c,d)`, the field `phi=(n_a-n_b)*(n_c-n_d)` has no constant or one-body
component. Its variance and lag covariance satisfy

```math
v_q=4\binom{V-4}{q-2}/\binom Vq,\qquad
C_q(t)/v_q=h^*P_2^t h/4,\quad 2\le q\le V-2.
```

The four-E experiment reads this statistic with complete back-action:
the first pair uses outgoing contrasts `d-2j`, the final pair incoming `d`.
Every pair retains all 625 formal raw slots, including zero maps. The
normalized curve is independent of global occupation on a fixed graph,
not independent of graph size or clock choice. Different primitive orders
can share eigenvalues. An initial half-trace error `epsilon` contributes
at most `2*epsilon/v_q` after normalization; fixed small `q` makes `v_q`
decay as `V^-2`. No practical sampling bound or measured mass follows.

The [spatial attachment](../../model/native_spatial_scaling/PROOF.md) treats
the pinned midpoint subdivision and successful injective port routing. It
extends the local halo incidence argument through that combinatorial
family. Its affine cell-centroid readout is an explicit observable,
not a selected physical metric. At level `ell`, the sum of squared crossing
coordinate differences is
`20-4*sqrt(5)+(-40+16*sqrt(5))/(3*2^ell)`, below 12. In the stationary
single-unit sector this gives mean squared displacement at most
`T/(12600*16^ell)` in `T` global attempts. A fixed set of `s` incident seams
and its actual `k`th selection instead have a separate stopped-clock bound
`k^2/(10*s*4^ell)`. Substituting the mean stopping time into the deterministic
bound is not the proof. The literal geometry checks reach level 3; native
expanded evaluation reaches level 2, and the all-level combinatorial result
is conditional on the source routing returning its required postconditions.
A separate specified-sender experiment has a positive first-incident-tick
arrival word of probability `(2M)^(-D)` along a shortest path of length `D`.
Incident ticks therefore do not define a one-hop support cone. These results
constrain the declared coordinate and clock dictionaries, without selecting
a physical metric, light cone or bulk space dimension.

The same-source continuation must supply one prepared reference, complete
field observables and a clock/refinement map that preserve their joint
records. A proposed three-dimensional field must explain how it is related
to the [two-dimensional native density limit](../../model/physical_space_time/README.md),
or define and justify a different source observable dictionary. Uniqueness
of that native density limit retains its separate conductance-convergence
obligation. Derive a controlled action on those observables before assigning the primitive
response blocks to particles. The tree detector, reflected probe and
spatial estimate have compatible native endpoint laws; their distinct
state, clock, number and observation contracts must be checked in any
composed experiment. In particular, the displayed displacement bound uses
`q=1`, where the degree-two contrast vanishes; its spatial bound cannot
simply be assigned to the nonzero degree-two probe. The common-occupation
response below has its own instrument and estimates.

### Prepared pair responses on one source

The [common native response proof](../../model/particle_identification/NATIVE_FIELD_RESPONSE.md)
joins preparation, spatial cross-correlations and a local exchange response
on the fair connected source. Fix `2<=q<=V-2`, independent uniform seam
selection, fair endpoint coins, repair probability one and one native
attempt clock. The comparison reference is `tau_q=I/binom(V,q)`. This
theorem does not identify the fair branch with the weighted H66 observer.

For two disjoint seams in each setting, the
[spatial pair matrix](../../model/native_pair_response/PROOF.md) is

```math
C_{ij}(T)=c_q h_i^*P_2^T h_j,\quad
c_q=\binom{V-4}{q-2}/\binom Vq,\quad v_q=4c_q.
```

The working source retains its chosen `q`; `P_2` is an exact coefficient
representation. Past outgoing contrasts and future incoming contrasts
include the four E reads' back-action and every raw outcome. Different
settings may overlap, so the Gram matrix `h_i^*h_j/4` and explicit smearing
weights are retained. It is not an independently prepared two-unit source.

[Native warmup](../../model/native_correlation_preparation/PROOF.md) prepares
this signed expectation from any initial `q`-sector density, including an
untouched correlated reference. With
`N=binom(V,min(4,q,V-q))-1` and a certified gap lower bound `gamma_0`, `n`
unconditioned warmup attempts give normalized reference-valued error at most
`sqrt(N/v_q)*(1-gamma_0)^n`. This uses a degree-at-most-four scored effect,
not convergence of the full density matrix. All warmup records remain;
the comparison ignores their values without conditioning. The exact graph
bound `gamma_0=2/(M*V*D)` consumes a proved diameter bound `D`. Integer
half-contraction blocks supply sufficient counts; budgets beyond the
existing single-program ceiling are not reported as completed runs.

The [exchange-clock theorem](../../model/native_exchange_clock/PROOF.md)
uses the same source operator `G_q=sum_e(I-S_e)` and
`P_q=I-G_q/(2M)`. At `t=T/(2M)`, every occupation sector satisfies
`||P_q^T-exp(-t*G_q)||<=3/(5T)` for `T>=1`, with zero error at zero.
Thus the preparation allowance adds to this endpoint-response allowance.
For primitive degree at most two and maximum graph degree six, the
rescaled logarithm has the stronger generator-domain statement
`||-2M*log(P_q)-G_q||<=144/(M-12)`. That generator estimate does not hold
uniformly on the full sector: at half occupation its error is at least
`M/16`. Neither exponential comparison changes the actual attempt-count
distribution or implements a unitary clock.

The pair comparison `R=P_2-A_adj/M` retains collisions as absorbing failures.
Its exact positive-entry telescoping formula controls the difference from
exclusion. On a four-port path, the normalized one-attempt difference is
`1/12`. For initially separated ports at distance `d`, its collision
probability is at most
`min(1,T/M,binom(T,d)*(Delta/M)^d)`. These are auxiliary pair-kernel errors,
not a replacement of the source's quantum instrument. For local probe
records, include every observed seam endpoint and retained terminal vertex
in the window roots; the native cone bound contributes `2*p_escape/v_q`
to a normalized score. A host-stop mass contributes `p_stop/v_q` when
stopped outputs are assigned score zero. Success conditioning is excluded.

General quadratic snapshots use a separate observation map: retain the
whole pure-E body tape, append one complete native E tree, decode its
body-terminal mask using the public `q`, and reverse all body currents.
This determines all earlier occupation profiles without an initial census.
The detector's disturbance follows the body. The local-probe cone bound
does not apply to this global tape and detector.

For bounded mean-zero profiles `f,g`, exact removal of constant and
one-body terms gives the primitive quadratic
`W=Y_f*Y_g-<f,g>/V-a*sum_i f_i*g_i*(n_i-q/V)/(V*chi)`, with
`a=(V-2q)/(V-2)` and `chi=q*(V-q)/(V*(V-1))`. Its last term has
mean square at most `a^2*||f||inf^2*||g||inf^2/(V*chi)`. On the owning density
subsequence, `q=V/3` and attempts `floor(1260*16^ell*tau)`, exact canonical
subGaussian bounds supply uniform eighth moments and covariance convergence
to Wick products of the diffusive Gaussian field. The quadratic variables
are generally non-Gaussian; their transfer is the symmetric tensor square
with generator `A tensor I + I tensor A`. This bounded-smearing limit does
not establish an interacting relativistic sector or the same limit for
fixed microscopic seam contrasts. No physical identification row is discharged.

### Complete local histories and bounded quantum menus

The [local-history theorem](../../model/native_local_history_limit/PROOF.md)
uses the fair binary endpoint law, independent uniform seam choices and a
fixed C/E mixture. Its observation retains the ordered complete raw records
on chosen seams, a common time-bin label, the final local source and an
arbitrary untouched reference. Other raw values and exact global tick
addresses are traced; the original archive remains readable. With
`T=M*t`, each bin boundary must satisfy `M*t_j` integer. This is the
per-seam attempt scale; exchange-action time is `t/2`.

For roots `S`, a radius-`r` window with `v` ports and `K` internal seams,
maximum degree `Delta`, and initial joint window/reference error `epsilon`,
the complete local output differs from its rate-one marked-Poisson window
comparison by at most

```math
\min\left(1,\epsilon+
\min\left(1,|S|\binom{T}{r+1}(\Delta/M)^{r+1}\right)
+\min(1,tK^2/M)\right).
```

All endpoint coins, law tags, impossible slots, currents and disturbance
are part of the compared instrument. The Poisson law is a comparison,
not an admitted native scheduler. Exact timestamps or global record lengths
would distinguish it. For the canonical initial state and an independent
reference, `epsilon<=v(v-1)/(2V)+v*abs(q/V-p)`, capped at one. Correlated
references require a joint bound. Stable rooted windows and compatible
states give a conditional local-process limit; bounded-degree compactness
provides a subsequence, without selecting a unique rooted refinement.
The literal cell-zero radius-one checks cover levels 0,1,2. They do not
prove all-radius convergence of the routed source constructor. Finite
Poisson enumeration retains an overflow allowance, and native host stops
are separate outputs.

The [coexistent-measurement construction](../../model/native_coexistent_effects/PROOF.md)
uses the distinct weighted H66 carrier and its declared coherent services.
For each Gaussian-rational `v_j` and rational `s_j` in `[0,1]`, two existing
phase words give success `L_j=I+(s_j-1)P_j` and abort effect
`(1-s_j^2)P_j`, where `P_j=v_j*v_j^dagger/(v_j^dagger*v_j)`.
With `A_j=L_j...L_1`, the effects
`A_(j-1)^dagger*(1-s_j^2)*P_j*A_(j-1)` and residual `A_m^dagger*A_m`
sum exactly to the identity. The effects can be noncommuting. Earlier
measurement disturbance is included; a separate rank-one test cannot be
substituted for a stage of this common instrument.

There are at most eight stages, two individually admitted 58392-pulse
nine-bit words per stage, and 934272 pulses along the surviving path.
Inputs have at most 64-bit rational components. Passing every emitted
primitive's exact stored-array check gives complete half-diamond error
below `1/50000` and gain below `10001/10000`. Each first processor failure,
one observed padding flag per word with 446 internal Kraus rows, and the
final quantum residual are retained. The generic host's six-bit/256-pulse
stop is unchanged; direct tape and lazy instrument queries have their own
domains. The construction supplies a bounded coexistent family, without
universal Born-premise coverage or a growing quantum-service embedding.
It cannot be applied to the local-history theorem's fluctuating-number
window merely because both constructions have source records.

### Contact poles and the transfer-ground dictionary

The [contact calculation](../../model/native_contact_response/PROOF.md)
uses the actual fair twelve-port source at `q=2`, with
`tau=I66/66` and `P=I-G/60`. On its primitive pair space, the centered
independent-pair comparison is `H0=H+2Q*A_adj*Q`, with
`H0>=2*(5-sqrt(5))*I`. The coefficient two follows from source exchange;
it is not fitted. A complete reflected probe on seams `(0,1),(2,4)` sees
all twelve distinct primitive poles. The lowest root of
`x^3-30*x^2+280*x-776` is approximately `5.01827276928`, with normalized
residue approximately `0.0194379829330`, certified between `1/100` and
`1/50`. Its separation below that finite comparison threshold exceeds
`2/5`. The raw four-E variance is `2/33`; 480 unconditioned warmup attempts
give normalized score error at most `33/65536`.

Every normalized vector in the fifteen-dimensional cubic subspace has
connected-federation action leakage outside its original cell greater
than `1/10`, for every valid three-port attachment. The proof covers all
220 attachment triples and checks the actual levels 0,1,2. Thus that
unchanged cell vector is not a capacity-stable eigenvector. A dressed
candidate must supply its residual, spectral window and nonzero accessible
residue: `||(G_l-E_l)f_l||<=eta_l` implies spectral weight outside
`[E_l-delta_l,E_l+delta_l]` at most `eta_l^2/delta_l^2`. It does not supply
those inputs. At growing `V` and fixed `q=2`, raw variance is
`8/[V*(V-1)]`; a normalized line can have vanishing raw weight and a
sufficient independent-shot budget growing as `V^4`.

The [transfer-ground response](../../model/native_vacuum_response/PROOF.md)
constructs `H_M=-2M*log(P_q)>=0` and its unique constant ground vector
`Omega` on the transfer Hilbert space. Complete four-E lag scores recover
`W_ij(t)=<F_i*Omega,exp(-it*H_M)F_j*Omega>/v_q`, a positive-type matrix
correlation. On bounded primitive degree `r`, choose
`b=floor(M/(2*r*Delta))` for `M>=2*r*Delta`; then `P^b>=I/2` and a
rational-complex polynomial in `P^b` approximates this response. Its
truncation and coefficient amplification are capacity independent at
fixed `r,t`, while the required native lag grows with `M`.

For coefficients `a_j` and `A=sum_j(abs(Re(a_j))+abs(Im(a_j)))`, systematic
normalized lag error `eta` contributes `A*eta`; score-zero stops add
`A*p_stop/v_q`. Independent samples have a separate coefficient-amplified
bound. This is signed postprocessing of several complete experiments,
not one acceptance probability or executed transfer unitary. An ideal phase-kick
witness gives a positive response in the transfer ground while bare
C/E-record laws are invariant under the corresponding ideal diagonal phase
kick, even with an arbitrary reference. The terminal working quantum
state is excluded from that blindness statement. A physical dictionary
must therefore supply a different state/intervention/readout map or a
prepared implementation that passes this witness. Recovered positive-energy
correlations alone do not identify the working tracial state with a vacuum.
An executed phase with complete half-diamond error `delta` adds a bare-score
allowance `2*delta`. The concrete order-12, `t=1/64` continuation has
truncation `1/4096`; its sufficient level-zero program has 116424316
attempts per setting, while levels one and two exceed `2^32`. Its
conservative 99-percent independent-sampling budget exceeds `10^17`
repetitions. These bounds do not assert practical execution or a universal
resource lower bound.

### Local preparation, complete readout and a shared attempt clock

Three owning packages provide the finite construction on the fair source:
[preparation](../../model/native_fair_local_preparation/PROOF.md),
[readout](../../model/native_fair_local_readout/PROOF.md) and
[clock composition](../../model/native_fair_local_clock/PROOF.md).
Their `packet --check` producers bind code, proofs and exact source controls.
The local reference is a channel, rather than an assumed input state:
it replaces each nonconstant twelve-port sector by its uniform density,
keeps its actual input weight and reference block, and preserves the entire
empty/full coherent corner. An eleven-E tree and `b` two-sweep blocks
implement it with half-diamond error at most `2^-b`. Every sector's exact
integer minorization is checked. The comparison retains tree records and
traces mixing-record values; every actual archive entry remains available.
Feedback consulting those marginalized values needs a separate bound.

For readout, incident currents undo boundary motion in each observed tree
contrast. The resulting initial occupation pattern is unique except for
constant alternatives. Reading the declared global total removes an
alternative only when it is infeasible: empty requires `q<=V-|B|`, full
requires `q>=|B|`. Consequently a fresh eleven-E root-cell tree is exact
for global `0<q<12` on every feasible local sector, including the empty
one. Sixty-four preparation blocks followed by this tree use 3862 native
attempts and nodes, with comparison error at most `2^-64` for arbitrary
input and untouched reference. No uniform time or occupation-density limit
is inferred from this capacity-independent attempt count.

At general density a one-command uniform E word acquires every chosen
tree edge with failure probability at most `(|B|-1)*2^-floor(T/M)`.
Coverage failures and unresolved constants have explicit erasure flags;
there is no success normalization. The full quantum/classical output
distance from an occupation-label oracle followed by the same instrument
equals their total probability. For canonical global marginal `q=V/3`,
a certified level-six 348-port window with `64M` read attempts, preceded
by 64 preparation blocks, has error below `10^-16`, 3852 nodes and
165154571 attempts. Unitality preserves that canonical marginal through
preparation even with a correlated reference. The preparation theorem
does not establish the initial canonical premise. The generic expanded
host retains its geometry and branch stops.

Both native programs use the actual attempt counter. The separate binned
local-history comparison charges reserved fixed operations by
`1-exp(-K*L/M)<=K*L/M`, together with the selector, cone and joint-state
errors. Its short-time window bound cannot be substituted for the long
uniform-readout bound. Native exact ticks and full global archives are not
identified with continuous Poisson histories.

The [addressed occupation coupling](../../model/native_fair_local_preparation/COUPLING.md)
intertwines every raw C/E branch with the twelve local occupation bits
and exterior. Under the existing supplied proper-code load/export premise,
a fresh tree-conditioned branch admits a beam/phase/inverse-beam experiment
and complete compiled E detector in 305 pulses. Every local occupation
sector, constant outcome, processor failure and rounded number leakage is
retained. The ideal contrast is the nonconstant branch probability `a`;
the certified stored-map score increment is at least `a*9999999/10000000`.
The rounded work register stays in the apparatus. A quantum wire is not
reconstructed from classical records, and no pulse-to-attempt physical
time conversion is supplied. The sector-uniform reference alone is
invariant under these number-preserving controls.

The physical source-to-field obligations concern the interacting prepared
state, autonomous coherent action and common physical clock, with
capacity/refinement control of complete observables. Finite occupation
readout and the conditional phase experiment do not discharge universal
effect coverage or identify particle species.

### Controlled ground, rooted windows and connected contact

The [controlled transfer experiment](../../model/native_transfer_control/PROOF.md)
uses the weighted 66-state source, its reflected scalar transfer
`S=Pi^(-1/2)*P*Pi^(1/2)` and `H=-log(S)`. The source reference is
`pi`, with weights `(25,25,8,...,8)`; its pure transfer ground has
amplitudes `Omega_x=sqrt(pi_x)`. An eleven-E detector reads the initial
and terminal configuration. Record-controlled routing and the existing
nine-bit amplitude compiler then prepare `Omega`. For arbitrary input
and reference, the ideal prepared source factors from the retained archive
and reference; the archive can remain correlated with that reference.

A fixed rational spectral proposal passes independent exact Gram and
source-residual gates of `10^-9`. A polar basis is only a proof device;
66 rank-one phase words on the actual rational rays approximate
`exp(-it*H)` for declared 64-bit rational `|t|<=1`. With the preparation,
phase intervention and final E-tree readout, the complete protocol uses
68 separately admitted coherent words, at most 3,877,451 native pulses
and 22 E commands. Every emitted stored primitive must pass its interval
gate. The complete half-diamond error is below `1/11000` and gain below
`1001/1000`; every processor failure and one padding-abort flag per word
remain. Numerical outputs have bounded normalization error and are not
asserted to be exact probabilities. At `t=1/10`, a phase on configuration
`(0,1)` produces a certified implemented score increment greater than
`1/2000`, after the full error allowance.

This supplies a controlled preparation/intervention dictionary on that
finite service branch. Bare unread repair remains a contracting mixed-state
process. The transfer clock differs from the working modular clock, and
its controlled parameter is not physical elapsed time. Direct nine-bit
compilation and lazy complete trees retain the generic host's word and
branch stops. No growing coherent-service embedding or common physical
vacuum follows from this experiment.

The [rooted refinement](../../model/native_rooted_refinement/PROOF.md)
belongs to the separate fair federation. For the five-face corner at
vertex zero, insertion of child0 gives identical addressed cell balls
whenever `2^ell>r+2`. This proves an all-radius, full-sequence cell-topology
limit. It does not prove a full-sequence limit for the floating port routing.
A local evaluator reproduces the finite source expressions, including
coordinate cleaning and tie rules, without allocating the remote sphere.
Exact binary32 score intervals certify the routed windows at radii
1,2,4,8 through specified levels up to 14; ambiguous decisions are rejected.
The finite degree/level limits and the integer edge-key domain are explicit.

For declared density `p`, the canonical family has
`q_ell=floor(p*240*4^ell)`. Its restriction to `v` window ports retains all
number sectors, with Bernoulli comparison error at most
`v*(v-1)/(2V)+v*abs(q/V-p)`, capped at one. In the radius-four example,
a 348-port window serves as the dynamical neighbourhood for records on
only the thirty root-internal seams and the final twelve-port root state.
Other window records, remote values and exact global ticks are traced by
this observation. An independent reference is allowed; correlated input
requires its joint restriction bound. At level 14, `p=1/3` and per-seam
attempt time `t=1/630`, the complete local-output comparison is below
`10^-6`. The associated `2^28`-attempt program is admitted, but expanded
execution returns `RESOURCE_STOP_UNEVALUATED_CARRIER`. It is a controlled
kernel comparison, not a completed high-capacity run. A routed continuation
must specify arithmetic, normalization, score ties and addresses at every
radius; neither finite routed tables nor the exact cell topology supplies it.

The [connected contact response](../../model/native_connected_contact/PROOF.md)
uses that fair source at `q=2`, with canonical reference and actual four-E
records. An explicit affine quadratic coefficient projected to global W2
proves that the compressed independent-pair comparison has
`0<lambda_min(H0)<=550/V`. The isolated cell's pole above five is therefore
above the connected comparison bottom, even at level zero. This is a
variational upper bound, not an inference from a vanishing lower bound.
It does not exclude dressed resonances or bound bands at fixed momentum.

Keeping all 66 within-cell pair configurations, write the global exchange
operator in core/exterior blocks. Its exact core is
`A=G_cell+N_boundary`, with `B^*B=N_boundary<=2I` and `0<=D<=24I`.
At negative spectral parameter `-s`, `s>0`, the Schur complement is
`F=A+sI-B^*(D+sI)^(-1)*B`. Its degree-`J` exterior-path truncation satisfies
`0<=F^(-1)-F_J^(-1)` and operator error at most
`2*(24/(24+s))^(J+1)/s^3`, uniformly in capacity. Intermediate core states
are not projected back into cell W2. At `s=24,J=2`, the cell-zero controls
at levels 0–2 give normalized response in
`[0.660081931825,0.660515959604]`, width at most the exact `1/2304`
allowance before decimal rounding. The isolated response exceeds the
connected enclosure by more than `1/120`.

With native discount `q_s=2M/(2M+s)` this response equals
`(1-q_s)*sum_n q_s^n*C(n)/v_V`, where `v_V=8/[V*(V-1)]`.
A cutoff at `N` adds at most `q_s^(N+1)`. State preparation and independent
sampling require their own bounds; host stops remain distinct outputs. The discount
is an analysis parameter. The negative-parameter bound does not invert
`D-E` near the positive cell pole and supplies no dressed-particle
persistence or physical mass assignment.

The next source-level constructions are therefore specific: prove an
all-radius routed arithmetic continuation, extend the controlled quantum
services and prepared reference to the same local number mixtures, and
establish a positive-energy connected sector with nonzero complete-record
residue and propagation control. These finite results leave the eight
physical requirements and all 21 identification obligations undischarged.

### Local states across capacity

The [canonical restriction](../../model/native_capacity_restriction/PROOF.md)
uses the actual fair federation with `V=240*4^ell` ports and a declared
conserved global occupation `m`. Its stationary canonical state
`rho_(V,m)=I_m/binom(V,m)` restricts to any `v` ports as

```math
\rho_W=\sum_{k=0}^{v}
\frac{\binom vk\binom{V-v}{m-k}}{\binom Vm}\rho_{v,k}.
```

All number sectors, including zero-weight sectors, are retained. These
restrictions nest exactly. The half trace distance to the canonical
sharp-`q` coarse state is `1-h(q)`, where `h` is the displayed number law;
for an arbitrary state supported in that sector this is a lower bound.
Thus choosing one fixed number independently at each size is not an exact
state-refinement map for a nontrivial global canonical family.

A fixed local window has a different comparison: its distance to the
Bernoulli product state of density `m/V` is at most `v(v-1)/(2V)`, capped
at one, and its distance to vacuum is at most `vm/V`. The exact local law
can be evaluated using integer arithmetic without constructing the full
carrier. At declared level 200, all thirteen number outcomes of twelve
ports are evaluated at a 408-bit population with density `1/3`; the
Bernoulli bound is `66/V`. Any common complete local instrument preserves
this state-error bound. Boundary-crossing repairs require a larger
dependency region or a separate locality bound; correlated references
require a joint input-state bound.

This supplies a static local calculation and a precise state-compatibility
test. It selects no physical relation between `V` and cosmological capacity,
no occupation trajectory `m(V)`, and no exact finite preparation or dynamic
refinement map. The weighted H66 source behind the native record response
has weights `(25,25,8,...,8)`; the growing federation is fair. Their
different laws and references cannot be combined as one capacity family
without an explicit source comparison.

The [dynamic product restriction](../../model/native_capacity_restriction/PRODUCT.md)
supplies an exact growing embedding of the biased repair/hold program into
the admitted independent-cell products. Appended factors can be correlated
with the active source and reference: tracing them out intertwines every
complete adaptive history, including raw records, FIFO exports and matched
stops. The comparison verifies and relabels the known program identifier;
allocation and specification queries remain distinguishable. Returned host
executions require the common evaluated domain, at most 64 cells. This
result covers repair/hold programs, not the separate compiled coherent
clock services. It creates no interactions or connected physical geometry.
A direct same-record, single-attempt replacement by the fair federation
fails: its coin probability is `1/2`, while the biased seam gives `25/33`,
so the complete-output total variation is at least `17/66`.

Reading simulator specifications is permitted. The
[field-attachment audit](../../model/common_source/repair_field_attachment.py),
[staggered field](../../model/common_source/native_staggered_field.py) and
[native feedback matrix](../../model/common_source/native_repair_mass.py)
distinguish supplied controls from outputs of repair. The feedback matrix's
six-row work-register evolution does not itself provide spatial propagation.
Adding a kinetic or Higgs coupling would define an extension whose predictions
must be labelled accordingly; it would not derive that coupling from the
unextended source.

**The construction supplies a common finite defect, gauge-current and
pixel-dependent mass operator, with a complete record service and a controlled
quartic extension on its localized algebra. It does not establish the physical interacting
chiral Standard Model or its measured spectrum.** Distinguish the action chosen by
the producer, the information in accessible records, the inverse justified
by those records, and the physical interpretation being tested. These are
different steps even when one software implementation performs all of them.

The unchanged binary native transfer also supports a
[count-derived localization and spectrum construction](../../model/native_quantum_defects/PROOF.md).
Its faithful one-port phase-twisted functional has an exactly matched full
exterior and positive dimensionless transfer-energy excess. Complete native
histories preserve the exterior comparison while the contrast decays; no
native phase gate or physical preparation is inferred. The same source has
an exhaustive 66-mode spectrum and compact-support classification. Its
[nonquadratic interaction test](../../model/native_modular_dynamics/PROOF.md)
excludes every one-body quadratic in the fixed native occupation dictionary.
These direct source results use reconstructed history functionals. Their
modular dynamics, regional entropy response and physical boundaries are
spelled out in the [finite field-selection guide](OPH_PARTICLE_FIELD_SELECTION.md).

### Connected weighted references and clock normalization

The [weighted extension](../../model/weighted_source_family/PROOF.md) uses
the source federation graph and the literal endpoint C/E rule, with a declared
positive weight at every vertex. It is outside the native program class's
weighted-federation admission rule. For fixed total occupation `m`, its
stationary diagonal reference is `pi_m(S)=product_(i in S) w_i/e_m(w)`, where
`e_m` is the elementary symmetric polynomial. Restriction to `A` retains
every local number with probability
`e_k(w_A)e_(m-k)(w_outside)/e_m(w)`. The restrictions compose exactly. The
faithful product family `p_i=z*w_i/(1+z*w_i)` also restricts exactly for one
fixed fugacity and inherited weights. Its modular action restricts under
partial trace. The recorded current gap `j*log(w_a/w_b)` has no global
normalization or capacity term; no physical mass is identified by that gap.

Repeating the native weight word recovers the H66 reference conditional on
two loads in one cell. Cross-cell repair changes that local number. Replacing
crossing repairs by a complete hold control gives half-diamond error at most
`1-(1-beta)^t` after `t` global attempts. Each source cell has thirty internal
and three crossing seams. Under the uniform seam selector, the probability of no crossing through the
`k`th internal repair is `(10/11)^k`, independently of federation level.
The bound does not vanish at fixed local repair count merely by enlarging the
exterior. The positive local symmetrized generator is derived from the same
population law; its unitary continuation is not a native service or a physical
Hamiltonian. Physical capacity, occupation, fugacity, coherent control and
field scaling remain separate selection and comparison requirements.

The [weighted local-limit theorem](../../model/weighted_local_limit/PROOF.md)
gives a controlled state comparison for `B` copies of the twelve-port weight
word, at declared mean occupation `r=2` or `r=4` per cell. The canonical
total-number state with `m=rB` is exactly the Bernoulli product conditioned
on that total. Its product fugacity is respectively
`(sqrt(1256)-16)/1000` or `(1+sqrt(6401))/1600`; the weights determine the
individual probabilities `p_i=z*w_i/(1+z*w_i)`. For any `k` cells, the full
local number mixture satisfies

```math
\tfrac12\|\rho_{B,r}|_k-Q_r^{\otimes k}\|_1
\leq\epsilon_{B,k}:=
\min\left(1,\sqrt{\frac{k\log(240B)}{2B}}\right).
```

The proof bounds the conditioning probability below by `1/(240B)` and
distributes its relative-entropy cost using cell exchangeability. It works
for every finite `B`, without a central-limit approximation or global
state enumeration. At declared level 200, the one-cell bound is below
`1.673e-60`. The `r=2,B=1` state is H66; its growing one-cell restriction
fluctuates in local number and is not confined to that sector.

An exactly trace-preserving common local instrument inherits this bound,
including every record and the terminal state, with an independent reference.
That instrument must be defined on every retained local number sector;
the H66-only coherent services do not automatically act on this mixture.
A rounded completely positive implementation of trace gain `G` and target
error `eta` instead gives `G*epsilon+eta`, without capping it at one or
renormalizing its outputs. An arbitrary correlated reference needs a joint
input bound. These are static reference estimates: they choose neither the
physical density nor preparation, weighted native admission or a dynamic
refinement map.

### Local histories and finite-window response

The [native local-response proof](../../model/native_local_response/PROOF.md)
attaches an exact backward-cone reduction to the endpoint C/E instruments.
The observation retains selected complete raw records, their positions,
the final local source and an arbitrary untouched reference. Remote raw
values and the terminal exterior are explicitly traced by this observation;
the source itself retains its FIFO and exports. Fixed selectors and matched
input-independent clock metadata are allowed. Every state-dependent record
event seeds both its endpoints, in addition to the final measured region.
Reversing the repair order then includes every unread seam meeting the
accumulated support. The resulting conservative cone suffices to compute the
complete observed channel from its initial reduced state; minimality is not
claimed.

For a finite window `W`, retain the same global selector tape and replace
crossing repairs by identity in the comparison channel. If the cone stays
inside `W`, the two observed outputs agree exactly. Their half trace distance
is therefore at most the cone escape probability, plus any joint initial
window/reference error. This compares controlled local experiments, not a
smaller universe with a renormalized clock. Actual controller feedback can
signal across disjoint repairs: the retained witness has receiving probability
`1/4` versus zero. Adaptive remote selectors, state-dependent archive access
and unmatched stops must include those dependencies and are outside this
local comparison.

For one whole cell, the source graph has a useful fifteen-port window: twelve
cell ports and three external neighbours in distinct cells. Each neighbour
has five further exterior seams and none to the other two neighbours. These
facts are checked for all 420 cells at levels 0, 1 and 2. At other levels the
following formulas require that geometry as a premise. With
`M=630*4^ell` seams, independent uniform seam choices and repair probability
one, the surviving conservative cone
has four states, indexed by its number `n=0,1,2,3` of infected neighbours:

```math
Q_{nn}=1-\frac{3+4n}{M},\qquad
Q_{n+1,n}=\frac{3-n}{M},\qquad
p_{\rm escape}(T)=1-\mathbf1^TQ^Te_0
\leq\min\left(1,\frac{15\binom T2}{M^2}\right).
```

The exponent is the number of global attempts. The fifteen terms count
ordered two-edge escape witnesses. At level zero and thirty attempts the
exact bound is about `0.01499`. For a general radius-`r` window, root set
`S`, maximum degree `Delta` and seam probability `q_max`, a separate
first-exit walk bound gives
`min(1, |S|*binom(T,r+1)*(Delta*q_max)^(r+1))`. These rate examples use
repair probability one.

The cell formula starts from the whole cell and retains record events with
both endpoints there. Observing crossing-record values enlarges the root
set and requires the corresponding cone bound instead.

Stopping immediately at the `k`th internal cell repair gives a different
clock. The internal selection rate is `30/M`. Reversing its geometric waiting
blocks gives the exact surviving matrix

```math
R=30(D-B_0)^{-1},\qquad D_{nn}=33+4n,\quad (B_0)_{n+1,n}=3-n.
```

The error is at most `1-1^T R^k e_0`, bounded by `min(1,k(k+1)/120)`;
for `k=1` it is `677/50061`. Elapsed global time is retained, with no
conditioning on rare observed times or outcome values. The bound is
independent of `M` and does not vanish by enlarging the exterior at fixed
local repair count. For any retained score of magnitude at most one,
the expectation error is at most twice that bound. Its normalized Abel
average, with analysis parameter `0<=q<1`, has exact bound
`2[1-(1-q)*1^T*(I-q*R)^(-1)*e_0]`, also bounded by
`min(2,q/[30(1-q)^2])`.

The proof covers arbitrary correlated references when the compared initial
cone/reference states agree or have a stated joint error. Combining it with
the weighted static limit requires an independent reference or that joint
premise, and uses the declared weighted extension rather than weighted
stationarity of the fair federation. The all-duration kernel and Abel series
are mathematical response families. Returned native executions retain the
host limits of 64 steps and 65536 branches; the expanded four-command control
stops after three. No unbounded executed protocol, physical spatial dimension,
energy or scattering pole is inferred from these estimates.

### Complete finite readout and the clock comparison

The [configuration detector](../../model/native_effect_coverage/PROOF.md)
reads eleven spanning-tree seams with native E repair. Subtracting the
accumulated recorded displacement from each incoming contrast reconstructs
the initial edge difference. Connectivity and total occupation two determine
the initial configuration. Every nonzero fine history is rank one, and the
coarse effects are the 66 computational projectors. The proof retains all
raw slots, terminal configurations and exports. Streamed branch continuation
and an exact lazy evaluator implement the stated domains; the stock expanded
history evaluator has an explicit resource stop and is not identified with
the completed return of the lazy evaluator. The nine-bit coherent words use
the direct source compiler; the generic `execute_word` host has a six-qubit
ceiling.

Existing coherent routing maps the recorded terminal code to a fixed basis
code, supplying ideal repeatable pure preparation, or back to the decoded
initial code, supplying the ideal basis Lueders instrument. Composing this readout with the
existing real/imaginary pair analyzers gives 4291 settings spanning all 4356
Hermitian coordinates, with source compiler and complete-failure bounds.
The rounded maps have bounded normalization error; their actual positive
record functionals retain full rank without renormalizing the outputs.
This replaces the inherited current-only detector's timeout allowance on its
stated execution domain. The coherent service and quantum instrument semantics
are declared; tomography alone does not supply the operational equivalences
and coexistent-effect coverage of the universal Born theorem.

The [rank-one synthesis protocol](../../model/native_effect_synthesis/PROOF.md)
extends the executable measurement menu using the same nine-bit work register,
source pulses, export and tree detector. For any nonzero Gaussian-rational
66-vector `v` and rational `0<=lambda<=1`, with input numerators and
denominators bounded to 64 bits, it targets
`E=lambda*v*v^dagger/(v^dagger*v)`. A padding amplitude completes the
normalized preparation vector. Its inverse word routes the selected output
to a genuine source code before export. Every native record and first
processor failure is retained; the 446 internal padding Kraus rows contribute
to one observed export-abort flag. The target effect does not specify a
Lueders state update.

Integer interval arithmetic constructs 41-bit dyadic angles and exact Walsh
coefficients. At most 23579 existing pulses give full-work-register target
error at most `11/2^41`; certified stored-primitive error and all failures
give complete half-diamond error below `10^-6`, with bounded trace gain,
when that word passes its exact primitive gates. The generator and pulse
bound cover the stated input domain; the numerical gate is checked for each
emitted word and rejects a failed certification. No universal rounding
theorem for the numerical library is assumed.
The concrete unequal-weight projector onto `(3|a>+4|b>)/5` uses an
eighty-pulse word for the displayed pair. It is outside the positive cone of
one balanced analyzer, while this longer admitted word constructs it.
Direct nine-bit compilation and the exact lazy detector retain the separate
generic word and expanded-history resource stops. The finite precision and
pulse cap supply no arbitrarily accurate limit, general positive-effect
compiler or common disjoint implementation of every coexistent pair.
Universal operational Born premises and physical field access are therefore
not discharged by this bounded rank-one construction.

The [clock comparison](../../model/native_clock_comparison/PROOF.md) uses the
actual accessible Y probe and working reference. Its mean is zero at every
working modular time. Identifying the scalar-record transfer generator with
an evolution on that same working state instead gives a mean greater than
`1/2000` at transfer parameter `1/10`. Exact rational logarithm enclosures and
an odd-response Taylor bound establish this finite-time inequality. Including
native preparation, complete readout and a conditional half-diamond allowance
of `1/100000` per proposed clock service leaves output half-trace separation
greater than `1/5000`. This rejects the specified same-reference clock map,
including reference-preserving unitary basis changes. It does not execute the
transfer clock or exclude a different reference/encoding. Bare repair records
with the final working state traced are blind to interleaved diagonal phases;
reading the simulator specification and using coherent analyzers are allowed.

## The defect-to-particle argument

The [identification route](../../model/particle_identification/ROADMAP.md#how-a-protected-defect-becomes-an-observable-particle)
starts with a source operation creating a locally distinguishable state,
then requires protection, propagation, field coupling and complete detection
on that same state. Charge conservation alone does not establish a stable
particle: the charge can spread or be carried by other excitations. A
spatially localized wave packet can also contain many energies without
representing many species.

The [charged cube experiment](../../model/ground_spectral_source/README.md)
provides an explicit finite formation-to-readout chain for its supplied
action. It prepares the full-Fock ground, applies a local odd creation
operation preserving the initial exterior restriction, and measures the
same-action response with complete records and probability error below
`3/5000`. Its protection holds during the declared repair waits; its
78 spectral gaps are not 78 identified particles. The native ordering-memory
charge and the cube number flag are different constructions. Their proposed
identification requires a map preserving creation, motion and records.

The [protected reference theorem](../../model/protected_reference/PROOF.md)
joins the actual biased ordered-cover repair, conserved memory charge,
stationary state and modular-current covariance. The
[complete gauge transport](../../model/gauge_defect_transport/PROOF.md)
supplies a locally covariant edge operation on a declared matter/link
extension, including blocked motion. Their
[joint process](../../model/protected_reference/GAUGE_JOIN.md) has a common
stationary reference and an exactly preparable neutralized singlet corner,
with its field attachment and corner admission declared. These finite
results test parts of the required connection. Neither supplies a source-selected interacting
relativistic Hamiltonian or a physical particle assignment.

## One source, complete records and a common spacetime

The [interacting extension](../../model/interacting_extension/README.md)
supplies a complete finite lepton/Higgs experiment, conditional native
references, a free spatial control, an interacting abelian flux control and
a local chiral-projector electric calculation. These occupy distinct
carriers and states. The [native mass join](OPH_NATIVE_PARTICLE_MASS_JOIN.md)
specifies their inputs; the [mass-response guide](OPH_PARTICLE_MASS_RESPONSE.md)
gives the signed response and its numerical interpretation. None of these
finite calculations identifies a physical interacting chiral field by itself.

Reading the simulator, scheduler, observer, repair and memory specifications
is allowed. To establish an effective field description, identify a source
preparation or encoding `E`, its complete experimental process `S_M`, a
record/state readout `R`, and the field experiment `T_M`, for each admitted
menu `M`. One sufficient comparison is

```math
\tfrac12\|R\circ S_M\circ E-T_M\|_\diamond\leq\epsilon_M.
```

The output includes all records and the remaining accessible state. A
restricted-state or energy-constrained statement is legitimate when its
domain is explicit and preserved by each intervention. A continuum claim
requires compatible menus covering the claimed local observable algebra
and errors tending to zero in the declared refinement. Determinism, equal
matrix dimensions or a match of isolated spectra does not supply these maps.

The [source trace](../../model/particle_identification/ROADMAP.md#the-central-common-model-test)
locates the gap concretely. `native_program.members.record_clock_member`
constructs a native record-clock program; its interpreter applies the
specified endpoint law and record-dependent selectors. It does not emit
the 438/1596-dimensional Higgs Hamiltonians or prepare their chosen ground.
The Higgs producer separately reads exact family matrices, sets its scalar
and interaction coefficients, and diagonalizes the resulting Hamiltonian.
A common-model proof must connect this action, preparation and pulse clock
to one source program, including every failure and continuation record.

### Native event order and observable influence

The [complete native record response](../../model/native_record_response/PROOF.md)
uses the weighted H66 program with its actual uniform seam selector and
tagged C/E records. Unconditional E-record scores measure incoming
occupation contrasts. Their ordered correlations include the repair that
produces each record, and give the exact response
`R(z)=D(zI-P)^(-1)A`. All 65 nonstationary population modes occur as poles.
They are relaxation factors in attempted repairs, with diagonal record
effects; neither full quantum tomography nor particle masses follow.
The weighted occupation means do not close under a twelve-component
linear law: one attempted repair enlarges their span from dimension 12
to 19. The fair source's closed mean equation is therefore not a substitute.

The [reflected-record reconstruction](../../model/native_record_response/TRANSFER.md)
gives a positive transfer on the 66-dimensional Hilbert space generated by
these stationary scalar records. With `Pi=diag(pi)`, it is
`S=Pi^(-1/2) P Pi^(1/2)`, with `(296/495)I<=S<=I`, and `H=-log(S)>=0`.
Actual raw-label reversal makes every scalar cylinder reflection form
positive. Leaving incoming contrasts unchanged under reflection instead
gives an exact negative control. This supplies a finite Euclidean
representation. Its [controlled finite implementation](../../model/native_transfer_control/PROOF.md)
uses separately declared coherent services; a spatial field and autonomous
physical time evolution are not inferred.
The [short-record realization](../../model/native_record_response/SHORT_RECORD_RECONSTRUCTION.md)
uses the configured weighted source to recover all 721 selective population
maps. Complete single-record functions have rank 42; the constant, eleven
single contrasts and fifty-four two-record products have rank 66. Correctly
reflected moments of at most five attempts determine the Gram and insertion
matrices, with exact Gram floor `55756800/2129251819009`. Their inverse
reconstructs arbitrary finite stationary record histories, including zero
outcomes. The [finite algebra theorem](../../Lean/SourceLaw/RecordRealization.lean)
proves insertion uniqueness and all-word similarity under explicit inverse
premises. The exact source certificate supplies those premises. This closes
finite population observability and realization on the configured source;
it does not identify unspecified hidden architectures from five-event data,
recover quantum coherences or supply the common relativistic particle clock.

The working source's modular generator `K=-log(Pi)` is a different clock:
an exact source matrix element has `[S,K]!=0`. The reconstructed ground
direction `sqrt(pi)` is not a preparation of the native full working state.

The [native dataflow calculation](../../model/native_causal_dataflow/README.md)
replays complete bounded H66/H78 programs from classical basis preparations.
It authenticates positive-branch endpoint values, selector reads, FIFO
memory and archive writes, retaining zero-weight branches in the full law.
Last-writer dependencies always point to an earlier emitted event, giving
an acyclic graph and the hypotheses of the finite
[semantic provenance construction](../../Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean).

The endpoint footprint alone is insufficient for adaptive programs. A
first repair on `(0,1)`, followed by `(2,3)` for positive recorded current
and `(2,10)` otherwise, gives final port-3 occupancy probabilities `1/4`
and `0` for preparations `{0,2}` and `{1,2}`. Both distributions retain
all 625 raw histories. The second seam is disjoint from the first in
every branch; the controller carries the influence. With fixed second
seam `(2,3)`, both probabilities are `1/2`.

The full single-FIFO read/write graph is a chain. It is a conservative
operational dependency relation: writing a shared archive does not by
itself prove physical influence. That exact chain cannot be identified
order-preservingly and order-reflectingly with an open Lorentzian region
containing spacelike pairs. This does not exclude a physical Lorentzian
limit. A [complete-record comparison](../../model/native_record_causality/PROOF.md)
supplies a narrower semantics result: disjoint fixed endpoint commands
commute on all H66/H78 input operators after a reversible permutation of
their complete record slots, archive/FIFO resplitting and known program
identifier relabeling. The compared channels have diamond distance zero,
including arbitrary references. Literal chronological records can differ.
The theorem requires completed fixed schedules without a controller
consuming the reordered records; adaptive feedback, timeouts and
overlapping repairs have explicit counterexamples. It does not erase
storage dependencies or select a physical causal quotient.

### Twelve-direction transitions and the field cone

The [finite transition theorem](../../model/native_lorentz_transition/PROOF.md)
uses the actual twelve spatial vectors `u_p` with squared radius
`R²=phi+2`, where `phi=(1+sqrt(5))/2`. Supply the lifts `k_p=(1,u_p)`,
the metric `G=diag(-R²,1,1,1)` and an invertible real linear transition
`L`. The exact native evaluation matrix has rank nine on symmetric
four-dimensional forms and kernel `span(G)`. Therefore

```math
(Lk_p)^T G(Lk_p)=0\ (\text{all twelve }p)
\quad\Longrightarrow\quad L^TGL=\lambda G,\qquad\lambda>0.
```

Normalizing by `sqrt(lambda)` gives a Lorentz isometry. Future orientation
(`L_00>0`) and properness (`det L>0`) are separate requirements. Conformal
factors multiply under composition. After subtracting
`((L^T G L)_00/G_00)G` from `L^T G L`, the exact decoder bounds the
entrywise residual by `(1-sqrt(5)/5)` times the largest absolute null reading.
These are algebraic consequences for supplied transitions and null lifts;
the native source must produce their physical readout, linear overlap
maps and common scale. The finite `PSL(2,F5)` action is not the continuous
spinor group `PSL(2,C)`.

The [field-cone criterion](../../model/field_cone_compatibility/PROOF.md)
tests whether the same cone belongs to a given static wave equation.
Let its spatial principal matrix be `a>0`, with supplied invertible
readback `dx=J dz` and positive clock conversion `dt=alpha d tau`, `alpha>0`.
Its pulled-back tangent
form is `K=diag(-alpha²,J^T a^(-1)J)`. All twelve rays are null exactly when

```math
K=(\alpha^2/R^2)G,
\qquad a=(R^2/\alpha^2)JJ^T.
```

This closes a conditional cone-comparison calculation. It selects no
source clock, mass term, vacuum or interaction. In the fixed dictionary
`x=r n`, an unwarped radial product with one angular coefficient and one
global clock cannot equal this flat cone at both `r=1` and `r=2`; the
angular tangent metric scales as `1/r²`. The exclusion concerns that
dictionary, not curved spacetime, variable coefficients or other frames.

### Acceptance conditions for a common field limit

The [golden-source order limit](../../Lean/Geometry/GoldenSourceCausalLimit.lean)
has an eventual causal-order equivalence away from the null boundary for
its specified `Fin 3 -> Fin(fib n)` populations, neighbourhood read law
and layer durations. Its volume companion and flat-diamond fraction
provide further conditional geometric data. Neither the single-FIFO
program nor a rank-three response Gram proves those population and read
premises. The native occupation-field limit has spectral dimension two
on the sphere; the three-dimensional radial product and active wave
service have additional declared inputs.

A common source construction must therefore produce the observer-visible
event relation, spatial readback and clock, prove their cone/coverage and
count-volume estimates, and show that the field propagator on those same
histories has the matching principal form. It must also control state
preparation, interactions, gauge/chiral structure and refinement. Equal
cones alone cannot determine a mass spectrum or Einstein dynamics. For
gravity the same state must supply stress, geometric entropy variation,
volume constraint and the responsive metric limit.

The [local-history bound](../../model/particle_identification/CAPACITY.md)
allows a finite buffered calculation at fixed spacing, with error
`min(1, epsilon_0 + sum_i(delta_i/2 + eta_i + a_i + zeta_i))` for
initial joint-state error, locality, splitting, coefficient/clock mismatch
and complete instrument error. It requires these inputs uniformly on the
claimed source family. No identification of global capacity with volume,
spacing, copies or regulator cutoffs is made. This gives a precise test
for when a manageable simulation suffices; it supplies no missing
capacity-to-action or state-selection law.

Reproduce the finite spacetime interfaces from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q \
  model/native_causal_dataflow model/native_lorentz_transition model/field_cone_compatibility
python -m model.native_lorentz_transition.certificate --check
python -m model.field_cone_compatibility --check
```

The dataflow certificate is `model.native_causal_dataflow.dataflow.certificate()`;
its tests check the complete probability law and JSON roundtrip. These
finite checks do not certify the conditional continuum hypotheses.

## The common finite particle operator

The simulator carrier, memory rule and coherent pulse service are explicit
design inputs. Their consequences can be derived without proving that bare
overlap consistency uniquely selects that design. A design selected with
knowledge of measured masses retains that calibration provenance even when
its execution reads no measured mass values. The
[source particle map](../../model/completion/SOURCE_PARTICLE_MAP.md) separates
those inputs from the proved operator identities and physical attachments.

| Finite construction | Established result | Physical boundary |
| --- | --- | --- |
| [Certified winding cores](../../model/certified_winding_core/PROOF.md) | Complete finite spectral enclosures certify 16 two-dimensional kernels among 28 supplied winding matrices; 34 supplied contact cases have controlled full-Fock ground pairs and equal addition/removal gaps | The winding texture and contact parameters are supplied; these are not 34 identified species |
| [Radial current algebra](../../model/radial_current_map/PROOF.md) | The conserved sums of onsite even operators of the original radial carrier give eta-pairing su(2), with only u(1) in the number-preserving part | This carrier does not realize the full native Standard Model current algebra through those local operators |
| [Native flavor and memory](../../model/native_flavor_orbit/PROOF.md) | The native current commutant gives the complete canonical same-current family frames; full-record coefficients, mixing and correlated-channel restrictions are calculated | The fixed common-duration example fails the direct same-scale charged-lepton hierarchy comparison; its static mean is not its finite unread channel |
| [Lattice texture energy](../../model/lattice_texture_selection/PROOF.md) | The selected ordered spin state gives positive shell-degree energy bounds; finite energy excludes nonzero degree on sufficiently distant shells | Bare protected-spin C/E repair does not select the texture, and the absolute crossing-energy bound is not a lifetime theorem |
| [Pixel-dependent native memory](../../model/mcpr_native_memory_bridge/PROOF.md) | Two admissible choices of current-probe durations and a source gain realize the retained charged-lepton mass operator, with rational existence bounds and controlled correlated weak-pulse dynamics | The mass-law architecture has calibration provenance; the source coordinates and operational scale require their physical attachment |
| [Defect, current and mass inclusion](../../model/defect_current_bridge/PROOF.md) | One local enlarged carrier has the same defect modes, native current action and pixel-dependent mass operator; a complete record service converges to its full finite quadratic Hamiltonian | Family multiplicity, carrier and coherent service are chosen; occupation does not create winding, and no interacting relativistic continuum follows |

The positive operator inclusion uses the full radial operator \(h\), its
two-column kernel frame \(W\), and the actual native matter representation
\(M_{15}\otimes\mathbb C^3\). The three family copies are chosen, rather than
deduced from the two-dimensional radial kernel. With \(J=i\sigma_2\) and
\(C_x=(-1)^{x_1+x_2+x_3}J\), the source identities give
\(hW=0\) and \(C\overline W=WJ\). For the symmetric all-left Higgs kernel
\(K\), the local creation pairing is \(\Delta=C\otimes\overline K\).
The Nambu inclusion

```math
E=\operatorname{diag}(W\otimes I_{45},\overline W\otimes I_{45}),\qquad
BE=EB_0,\qquad
B_0=\begin{pmatrix}0&J\otimes\overline K\\
 (J\otimes\overline K)^\dagger&0\end{pmatrix}
```

also intertwines all twelve actual native currents. The source permutation
transports their original port representation, so the identification does
not consist of attaching gauge labels to an unrelated spectrum. The induced
current action is spatially uniform; this construction supplies neither
dynamical gauge links nor a local Gauss constraint. The induced
CAR inclusion intertwines quadratic dynamics. On even states, tensoring
with the normalized trace on the complementary modes supplies a recoverable
full-carrier state map; that complement is not a vacuum.

The memory construction supplies a positive family matrix \(M(P)\) whose
eigenvalues are the retained conditional \(m_i(P)/v(P)\). Its construction
uses the complete C-enabled stationary coefficient \(Z_C\), with
\(M(P)=gZ_C\), rather than conditioning away failed records. Choosing
\(Y_e=\sqrt2 M(P)\) at \(H=(0,v(P)/\sqrt2)\) puts those coefficients in
the actual Higgs kernel on the same localized modes.

This coefficient has an explicit record-driven implementation. Fresh
auxiliary current probes generate the complete memory law independently of
the spatial field. On every record apply the chosen local finite-Fock
Hamiltonian \(H_r=H_0+\mathbf1_C gL(P_j)\), where \(j\) is the updated
memory port and \(L\) is the lepton pairing map. Each record preserves the
same core inclusion, and its stationary mean is \(H_0+L(gZ_C)\). If
\(\|H_r\|\le M\), the reduced channel after \(N\) steps of duration
\(\tau/N\) differs from that mean Hamiltonian's unitary channel by at most

```math
\frac{4M\tau}{N(1-p)}+
\left(\frac8{1-p}+4\right)\frac{M^2\tau^2}{N}
```

in unhalved trace norm, including arbitrary initial classical-memory
correlations and untouched finite references. Stationary product memory
removes the first term. The bound is for a fixed finite volume and grows
with volume; it is not a continuum estimate. The finite unread channel
can be nonunitary, and this coherent service is distinct from the original
site-field repair.

For observations confined to the same 90-mode core, the
[interacting extension](../../model/interacting_defect_response/PROOF.md)
replaces the extensive norm by
`M_core<=45(1+g)+2025lambda` for `lambda>=0`, after adding
`lambda(N_W-45)^2`. This gives a radius-uniform complete-record limit on
that algebra, including untouched references. It does not give uniform
control of arbitrary spatial observations or lattice-spacing refinement.
Writing the unread bound as `delta`, the complete final history/core state
has error at most `sqrt(2 delta)` against the same history distribution with
ideal mean evolution. This includes initially correlated classical memory
with its conditional history laws preserved. Thus final record retention has
an `O(N^-1/2)` bound; it does not inherit the unread `O(N^-1)` rate. No uniform
rare-history or intervening field-measurement guarantee is inferred.
The interaction preserves native currents and the core inclusion but has
nonzero spatial tails. Its exact four-mode mean-Hamiltonian slice displays
connected correlations and interaction-dependent fermion poles. Its
Gibbs response supplies an entropy/energy/work identity on the same algebra,
without a geometric modular-flow or Einstein-stationarity conclusion.

The radius-one witness retains 2,430 annihilator modes and all 4,860 Nambu
coordinates. Its core has 90 annihilator modes and 180 Nambu coordinates,
including every zero direction. Those counts are operator degrees of
freedom, not species. The complete bands are \(\lambda(h)\pm s_j(K)\).
Their additive form does not prove relativistic massive dispersion. The
finite inclusion therefore connects the three specified operators without
identifying their states with nature's particles. A common interacting
vacuum, physical clock and mass matching, dynamical defect formation, and
the responsive Einstein limit require additional derivations.

The texture energy cannot be omitted from that test. Identifying the
prescribed radial hedgehog with the unscreened ordered-spin rotation field
would put degree one on every surrounding shell. The proved positive
per-shell bound then forces divergent total energy. This rules out that
specific identification as an isolated finite-rest-energy defect. It does
not exclude a screened defect or another localized topology with its own
derived field action.

The [compact-pair construction](../../model/compact_defect_pair/PROOF.md)
gives one explicit finite-energy alternative, with opposite local degrees
and a constant exterior. Its four localized levels are bounded and certified
away from zero in the displayed example, including in the controlled
infinite-volume limit. The current/pairing inclusion therefore retains a
nonzero `A=W* h W`. The
[compact interacting extension](../../model/compact_interacting_response/PROOF.md)
retains that operator and `lambda(N_W-90)^2` on all 180 localized modes.
Every native memory branch preserves this algebra; its bound
`M_core<=90r+90(1+g)+8100lambda` is independent of exterior volume.
The complete sixteen-level slice spectrum has competing even, mixed and
blocked ground sectors. Exact source intervals place all three retained
lepton slices in the blocked sector at `lambda=1/128`, with two poles
`e1+lambda +/- m` and residues `1/2`. A slice is not the full interacting
vacuum. The memory theorem concerns a fixed texture during its pulse interval.

[Guarded texture repair](../../model/compact_defect_transport/PROOF.md)
retains all twelve records of a continuing law, preserves both shell
degrees, and admits explicit endpoint transport of the band, currents,
pairing and interaction. The guards and coherent control are design
inputs. [The hopping-work law](../../model/compact_band_response/PROOF.md)
differentiates the same projected interaction as the band moves; its
moving-projector term cannot be omitted. Its work/entropy relation does
not turn the hopping parameter into spacetime geometry.

[Local native Gauss constraints](../../model/native_quantum_gauss/PROOF.md)
have a finite 60-state link realization with the actual five matter
multiplets, a static dual reference and complete endpoint readouts.
Its local Ward identity acts on the full finite Fock space. The compact
band is not dressed into this distinct Gauss sector. Controlled texture
protection, finite current conservation and interacting pole calculations
therefore have explicit operator domains, without a common physical
chiral continuum or laboratory mass assignment.

The [pixel feedback calculation](../../model/pixel_mass_closure/PROOF.md)
also uses distinct declared running/scale designs. Each certified root
feeds the retained mass formula and an admissible native-memory inverse.
This closes a conditional coefficient consistency problem; identification
with the electromagnetic response of these field records requires its
own current-response and physical matching derivation.

## What a particle means in this construction

A particle identification is a statement about a continuing observable
response in a specified state and spacetime/clock regime. Preparations,
currents and detector records must consistently exhibit the proposed
excitation, with its quantum numbers, propagation and interactions. A finite
internal basis vector, memory level, channel decay rate or agreement record
does not have that meaning by itself. There need be no fundamental particle
object behind the observers' statistics.

For a stable massive excitation in a suitable vacuum regime, one can seek an
isolated charged response with nonzero residue, the appropriate dispersion,
spin/statistics and current matrix elements. All parts must use the same
state, dynamics, observer clocks and accessible operations. For a stationary
reference state, a finite closed Hamiltonian gives a discrete Lehmann response;
this alone neither
establishes relativistic particles nor selects how many physical species
survive a volume or continuum limit.

Other physical regimes require other mass observables. Unstable excitations
use a specified resonance pole and analytic continuation, or an explicitly
controlled line-shape description. Confined quark parameters are defined by
renormalization and matching to physical observations, rather than by
asserting isolated observable colored-particle poles. Charged sectors with
massless gauge fields require an infrared treatment; a universal isolated
one-particle-pole requirement would be too strong. These qualifications do
not license arbitrary replacement of a failed pole identification by an
unrelated threshold. State the observable and its domain before comparison.
See the [PDG resonance review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-resonances.pdf),
[particle-mass and scheme conventions](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf),
and [Buchholz's Gauss-law obstruction under its charged-sector hypotheses](https://doi.org/10.1016/0370-2693(86)91110-X).

The intended charged-spectrum target is nine running charged Yukawas and
the four independent CKM parameters in compatible conventions. These are
not nine directly observed isolated stable-particle poles. A matching map
must connect them to the actual detector and interaction observables.

Mass–energy equivalence does not select topology uniquely. A real free scalar
field has a contractible field target and trivial vacuum topology, yet its
one-particle states satisfy relativistic massive dispersion; see the
[free-field construction](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S2.html).
This is a counterexample to inferring a topological sector from a dispersion
relation alone. Whether a microscopic realization of an effective particle
is a defect is an additional question about that source realization.

## The hypothesis of particles as localized defects

The [protected-object catalogue](../../model/protected_defect_catalog/CATALOG.md)
contains explicit finite classifications and source-to-observation links.
Its stated domains distinguish complete charge classifications from partial
lists of sectors and prescribed textures. The [candidate-assignment table](OPH_PARTICLE_FIELD_SELECTION.md#candidate-assignments-from-the-existing-catalogues)
uses these results to rank particle hypotheses by their signatures and
falsifiers. Such hypotheses guide calculation while physical assignments
remain subject to the same-source field and measurement tests.

The particle hypothesis is **a localized, persistent
topological defect in the native overlap/transport structure**. Its continued
maintenance by the same repair process would account for its rest energy.
The finite common operator above supplies fermions bound to a prescribed
texture. It does not prove that creating a fermion creates that texture or
that its eigenvalue is the texture's physical rest energy.

The proposed object is a bounded core whose surrounding records approach the
reference regime, but whose boundary transport carries a nontrivial invariant
that ordinary admitted local repairs cannot erase. The obstruction is to
trivializing the transport relative to that reference, not a permission for
observers to disagree about the same measured event. Moving the core must
transport its charge. Allowed creation, annihilation, fusion and decay must
follow from the same source law, rather than from freezing selected seams
or supplying punctures to protect an assigned label.

For stable massive matter the proposed object is a persistent excitation.
Unstable species could be
metastable or internal defect excitations with derived decay channels.
Massless propagation and gauge interactions must also be tested in the same
construction; their meaning cannot be supplied by selecting unrelated
numerical modes.

| Part of the hypothesis | Concrete test |
| --- | --- |
| Defect identity | Exhibit an actual native state/observable sector and a boundary charge invariant under the admitted repairs; derive its coefficient group and transport rules |
| Persistence and localization | Retain the full law and show that the core persists or has a controlled lifetime, with reference-like response outside it; distinguish active recovery from autonomous protection |
| Maintenance and energy | Derive the defect's excess action/energy relative to the reference using a calibrated local clock, with all kinetic terms from that same process |
| Mass and propagation | In a justified local relativistic regime, recover the same rest energy from propagation and dispersion; the rest-mass relation uses the common signal speed, not a conversion fitted per species |
| Particle identity | Obtain spin/statistics, charges, current couplings and any family structure on those same defects or excitations, before comparing mass ratios |

The tested double-load configuration fails this persistence requirement,
including its quantum collision subspace. The actively corrected matching
code supplies persistent information without the proposed topological charge
or maintenance-energy map. Neither object is identified with a physical
particle. The radial common-operator construction has a different specified
carrier and law; it does not reverse either of these exclusions.

## Choose one source object before comparing masses

The different calculated quantities are not interchangeable particle
definitions. Their roles must stay visible in every spectrum table.

| Calculated object | Legitimate role | Missing identification |
| --- | --- | --- |
| Representation components, charges and Casimirs | Internal types and allowed current/Higgs couplings | A propagating source excitation with those quantum numbers |
| Protected bands, memories and family kernels | Persistent information and candidate interaction coefficients | A physical family sector and a derived mass/kinetic interaction |
| Repair contractions and transfer eigenvalues | Attenuation or phase in their stated event/model parameter | A common energy/clock map and the appropriate measured response |
| Finite CAR addition/removal poles | Excitations of the specified interacting finite law, with residues and charge | Physical spatial, gauge, spin/statistics, continuum and mass matching |
| Supplied Weyl/Nambu symbols, staggered/dimer gaps and Hubbard thresholds | Separate conditional field or interaction benchmarks | A derivation from one selected native source object; their numbers cannot be pooled into a spectrum |
| Running charts and calibrated mass coordinates | Parameterization, matching controls and empirical comparison | Independent selection of source couplings and physical readout |
| Topological textures and their bound modes | Certified finite winding sectors, energy bounds and a common current/mass operator on localized modes | Creation of a particle as creation of a protected dynamical defect, with physical rest energy and transport |

A candidate identification must name **one excitation or sector of one
continuing source process** before its mass comparison. Specify the reference
regime, local preparation, support/localization, conserved charge or defect
invariant if present, admitted repairs, propagation and allowed creation/
annihilation channels. The same object must account for currents, detector
response, kinetic normalization and its mass observable. A separate calculation
becomes evidence for that object only after an explicit map preserves the
relevant state, operations and records. Selecting one convenient number from
each benchmark does not supply that map.

Each finite CAR experiment belongs to a named supplied-action benchmark.
The four-mode interacting ring and the chiral quadratic candidate are distinct
laws. The inverse-response diagnostic tests the ring's response;
the spectral census and record tomography test the same chiral
candidate, with tomography restricted to one exactly decoupled component.
They introduce no additional physical particle candidate or mass ansatz and
do not turn these benchmarks into one source. The common radial operator
above supplies a particular finite join with its own declared carrier;
its physical chiral-field attachment requires further premises. No benchmark
is promoted to a Standard Model particle because its spectrum resembles a
measured ratio.

For a topological-defect route, first exhibit the actual native configuration
or sector and an invariant under the allowed local repair/deformation class.
Show persistence, transport and the rules for pair creation, annihilation and
fusion where relevant. Representation labels or a persistent memory index
alone do not establish topology. The defect's spin/statistics and gauge
response must be tested on that same process; they are not supplied by the
word defect.

## When repair maintenance cost can explain mass

A rest mass can be the excess rest energy of a localized excitation over the
reference state in a regime where energy, clock and spatial dispersion are
defined. To identify this with repair maintenance, the source must derive
how the excitation changes an intrinsic action/energy readout or the same
clocked response. Subtract the reference contribution and account for cutoff,
volume and renormalization dependence where required. Then test propagation
and interactions with that same identification.

Counting repairs, timing host code, taking a relaxation logarithm, or declaring
an arbitrary function of a persistent register does not establish this energy.
Different schedules or implementations can alter resource counts without
changing the claimed physical response. A maintenance-cost law must explain
that correspondence, including calibrated local-clock changes, rather than
fit a conversion independently for each species. Attenuation can govern
lifetime or decoherence instead of inertial mass; the response calculation
must decide that distinction.

No established native topological-defect maintenance law yields the measured
masses. The specified memory and coherent field service realizes the retained
pixel-dependent mass operator without establishing that maintenance-energy
interpretation. That interpretation requires the source object and its
cost/energy correspondence; unrelated gap and kernel scans cannot supply it.

### Existing defect theorems and the direct native decision test

The [consensus paper](../../paper/reality_as_consensus_protocol.tex) does prove
a useful defect result: on a fixed graph with positive mismatch weights and
fixed edge offsets b, nonzero cycle holonomy forces a strictly positive
minimum mismatch. Repairs of the vertex assignments cannot eliminate that
frustration while b is fixed. The theorem does not supply an update of b,
localize a moving defect or identify the residual minimum with excitation
energy above its own reference ground state. Its particle interpretation is
a hypothesis requiring those further maps.

The [finite flux criterion](../../code/a5_closure/flux_defect_criterion_certificate.py)
and [native defect gate](../../code/a5_closure/native_defect_stability_gate.py)
also distinguish supplied link/puncture sectors from native records. The full
two-dimensional carrier complex has vanishing H1 while H2 does not; the seam
graph alone has a different cycle space. This is not a proof that all topology
is absent. A nontrivial bundle needs actual transition data. Native load
contrasts are exact differences and have zero face curl. Declared Z6 links
and punctures allow a separate flux experiment, but those are additional
state/update premises, and unrestricted local link changes can erase its pair.

The [native defect decision experiment](../../model/common_source/native_defect_decision.py)
and [complete exact receipt](../../model/common_source/native_defect_decision_results.json)
test a specific local candidate directly: double occupation of one port in
the native total-load-two C/E process. They retain all raw outcomes and compare
the same state transitions, transport and cost observables. The invariant
test uses the support of all native seams; lifetime and activity numbers use
an explicitly supplied independent uniform-seam/fair-coin probe schedule.
That schedule does not replace the production feedback policy.

| Test under the unchanged native law | Exact finite result | Consequence |
| --- | --- | --- |
| Configuration sectors | Twelve transient double-load classes feed one closed irreducible binary class of 66 states | A scalar configuration label preserved by every admitted branch is constant on this 78-state domain; this is not a classification of quantum or complete-history invariants |
| Local double-load persistence | Under independent uniform choice of the 30 seams, survival is (5/6)^n and mean lifetime is six seam-attempt counts | The candidate splits into two unit loads; it is not a protected defect transported by this law. No physical duration is assigned to an attempt |
| Quantum collision-sector persistence | Both complete C/E instruments satisfy sum_h K_h* P_collision K_h = (5/6) P_collision, with no binary-to-collision feed | The same survival law holds for arbitrary quantum inputs and untouched references. There is no persistent subspace wholly supported on these collisions under the stated uniform unread evolution; other quantum/history sectors are untested |
| Continuing binary transport | Nonzero current probability is 5/33 per attempt in the uniform binary stationary law, while every binary event has zero drop of Phi=one-half sum(load squared) | Repair activity and motion are not equal to potential decrease |
| Proposed edge disagreement cost | The retained binary edge cost takes both values 8 and 10 and has native transitions in both directions | It is not a monotone maintenance-energy law of this process |
| Changed-law control | Removing the five seams incident to one port yields binary closed classes of sizes 11 and 55 | An artificially frozen port creates a sector by changing the allowed law; it does not establish protection under the original law |

These results reject this proposed autonomous configuration defect and the
tested automatic cost-to-energy identifications, including quantum persistence
wholly within its collision subspace. They do not exclude other quantum,
history or federated defects, another repair law, a common energy offset, or
the general particle-defect hypothesis. A replacement candidate must identify
its actual source sector and allowed dynamics before any mass comparison.

The [actively corrected code](../../model/common_source/native_gauge_object.py)
is a different persistent object: coherent recovery preserves logical
information and actual port rotations act on it. The fixed matching
code leaks out of its subspace without recovery, by 31/108 for C and 55/108
for E on the normalized code projector under the uniform-seam test. This is
a one-step preparation-specific escape result, not an all-state lifetime.
Its complete corrected
record law has the form rho -> rho tensor p(records), so a fixed cost of
seams/records alone cannot distinguish its logical states. The
[fixed-code calculation](../../model/common_source/parameter_selection.py)
also gives scalar one-port and native-edge quadratic compressions. This is a
restriction on logical-state mass splitting from those costs, not a proof of
zero common mass. Its fixed matching table and the separately derived
equivariant table in [native code transport](../../model/common_source/native_code_transport.py)
must not be silently identified for a nonscalar observable or an additional interaction.

## The inference chain to retain

| Step | Evidence needed | Insufficient substitute |
| --- | --- | --- |
| Source law | A specified preparation, instrument/channel or coherent action, schedule, complete records and controls | Naming a matrix a field or interpreting a solver counter as elapsed time |
| Accessible response | Executed probes with all outcomes, references and continuation after measurement; declared clock and apparatus | A Hamiltonian supplied only to the inverse, successful-branch averages, or expectation values that discard the record tree |
| Identifiability | An inverse on a declared model class, with phase/bandwidth/rank conditions and indistinguishability controls | A fit to one oscillation frequency or a few momenta |
| Joint dynamics | The same response fixes kinetic normalization, gaps, residues and current/Higgs vertices | A spectrum with a separately fitted Gram matrix, label table or vertex |
| Particle content | Charge/current action, visible and dark modes, multiplicities, spin/statistics and controlled propagation | Sorted eigenvalues or an internal irreducible dimension alone |
| Physical regime | Locality/covariance, state/vacuum, gauge constraints, chirality, interactions, limits and errors | A finite free symbol, a supplied background or removal of unwanted branches |
| Measurement comparison | Fixed dictionary, scheme/scale, threshold/width/infrared treatment and joint uncertainties | Reusing empirical inputs as validation or comparing running coefficients directly with pole masses |

A valid detector inverse need not be unique among all conceivable laws. Its
claim must state the class it identifies, the observations it preserves and
the competing laws that remain indistinguishable. A declared apparatus can
be a legitimate design choice. A proof that the apparatus was selected by
native repair is a separate claim.

## Record reconstruction and supplied-action spectroscopy

The [four-mode response calculation](../../model/common_source/particle_identification/RESPONSE.md)
computes the complete interacting finite Lehmann response. It establishes
which encoded charged fields have residue at the gap measured by the
continuing clock experiment. Its global charge is occupation, the ring is
spinless, and its Hamiltonian/clock are specified. The result does not supply
an electromagnetic, color or relativistic identification.

For fermionic interpolators O_a in a stationary finite system, the retarded
anticommutator response has the spectral form

```math
G_{ab}(z)=\sum_n\frac{R_{n,ab}}{z-\omega_n},\qquad \mathrm{Im}z>0.
```

The complete residues, reference state and signed frequencies are part of
the response. Normal and anomalous sectors require their stated CAR/Nambu
convention. In a canonical complete interpolator frame the zeroth moment is
the anticommutator matrix. A restricted detector menu can leave a mode dark;
zero residue in that menu is not a proof the mode is absent. Conversely, a
Hamiltonian eigenvalue with no accessible residue is not an observed line.
Repeated eigenvalues require spectral projectors or another invariant
subspace calculation, rather than basis-dependent eigenvector labels.

Away from singular or ill-conditioned points, the same response determines

```math
\Gamma(z)=G(z)^{-1},\qquad
\partial_z\Gamma=-\Gamma\,(\partial_zG)\,\Gamma.
```

This identity is a finite-response diagnostic, not a continuum reconstruction
theorem. Under an invertible change of interpolators,

```math
O'=S O,\qquad G'=S G S^\dagger,\qquad
\Gamma'=S^{-\dagger}\Gamma S^{-1}.
```

The kinetic and mass
forms must transform together. Generalized poles remain invariant; ordinary
matrix eigenvalues can change. Incomplete response data cannot determine a
missing kinetic form by declaring a convenient inner product.

The [Gaussian CAR source](../../model/common_source/chiral_candidate/SOURCE.md)
executes the declared quadratic action through local native control words,
including fermionic swaps, full work/reference covariances and complete
occupation records. This is stronger than diagonalizing a supplied symbol.
The coefficients, neutral Higgs background, initialization, control services
and physical field identification remain separate inputs. A covariance
representation is exact for each Gaussian branch; averaging distinct
conditioned branches and reapplying Wick's rule to that mixture is generally
invalid. Coherently controlled phases are not recoverable from covariance
alone.

The [protected momentum experiment](../../model/common_source/chiral_candidate/protected_momentum.py)
recovers oscillation frequency and visibility from its declared finite probe
records. Those marginals do not identify every matrix residue, sign, phase,
current or kinetic derivative of a general charged field. The
seven-momentum agreement also leaves local continuations that disagree on
a second shell. Tests of additional preparations, phases, momenta and Higgs
perturbations must be withheld from the inverse and predicted afterwards.

A response reconstructed from record expectations must count the apparatus
that makes those expectations accessible. Parity-even probes need not be
charge invariant. Charged-field coherences require an explicitly allowed
reference/dressing and a charge-compatible preparation and detector, or a
stated additional phase-reference service. Native selection of that service
cannot be inferred merely from tomography being algebraically possible.

### Response inverses and record tests

The [inverse-response diagnostic](../../model/common_source/particle_identification/kinetic_response.py)
computes G inverse and its frequency derivative from the complete finite
four-mode Lehmann data. It validates actual residue positivity/Hermiticity,
the canonical CAR sum and the retarded domain before inversion, and reports
the actual integrated residue separately from the expected coordinate Gram.
Nonunitary coordinate controls preserve generalized poles while moving raw
matrix eigenvalues. The interacting derivative varies with frequency and is
not absorbed into a fitted static mass. This consumes the declared complete
Hamiltonian response; it is not detector tomography or a physical kinetic
metric.

The [charge-resolved spectral diagnostic](../../model/observer_particles/spectral_response.py)
tracks all signed subspaces and their matrix residues for fixed interpolators.
The supplied 180-mode Weyl/Nambu report retains dark and complementary modes.
Numerically unresolved clusters keep their full frequency lists; their weights
are not presented as one exact pole, and numerical ranks are not exact proofs.

The [record two-point experiment](../../model/common_source/chiral_candidate/record_two_point.py)
and [retained receipt](../../model/common_source/chiral_candidate/record_two_point_results.json)
reconstruct one charged L/e-conjugate pair, native family axis zero.
It uses neutral source pulses and a charged reference fermion in a fixed
total-charge-zero code. Complete detector/clock records reconstruct the complex
two-mode transfer and, under an explicit principal-phase bound and canonical
quadratic closure, its generator and simple-pole residues. Joint interventions
recover the z-axis momentum and neutral-background multiplier derivatives.
Neither a dynamical Higgs/current vertex nor an interacting wavefunction
renormalization is inferred.

The retained 304 contexts have 4,864 terminal outcomes with full unnormalized
continuation states. At 64 source splits, withheld momentum/background settings
and a second duration have generator residual at most about 6.40e-8 per unit
of reconstructed clock coordinate, and complete continued-output distance at
most about 4.36e-8. The two
durations repeat the same elementary source word; the composition defect is
about 2.23e-16. These are numerical replay diagnostics, not finite-shot,
outward or physical error certificates. The second duration alone cannot
eliminate unbounded generator aliases; the phase-bound prior is retained.

Controls exhibit opposite pairing signs with identical empty-pair Rabi
marginals, a two-hop law agreeing on training momenta but differing on a
withheld momentum, and wrong conditional states with unchanged record weights.
The complete response detects these differences. The source coefficient,
canonical field class, preparation/control services, reference and programmed
clock remain premises. This completed coherent-service experiment does not
replace native attempt, abort or export histories by unit success. Its mass
coefficient is supplied and is not selected by repair maintenance.

## Repair, coherent evolution and physical energy

The complete native instrument can be contractive while an attached closed
Hamiltonian evolves unitarily. The no-intertwiner result for the tested
conditioned response excludes a
nonzero direct linear identification of the tested conditioned native
response with the unmeasured closed quadratic evolution. It does not rule
out a detector-coupled or environmental effective description; the
[dephasing control](../../model/common_source/chiral_candidate/KINETIC.md)
exhibits the distinction explicitly.

The common defect/current construction has an exact CAR intertwiner for
its chosen local Hamiltonian and each record-specific coherent field pulse.
Its weak-step averaging theorem gives a controlled limiting unitary after
the complete classical memory is traced out. These statements do not equate
a finite unread native channel with a closed unitary, and hence do not
contradict the conditioned-response exclusion.

A contraction modulus describes attenuation in its declared update or event
parameter. The argument of a transfer eigenvalue describes a phase only
with an unwrapping and clock convention. Neither quantity becomes a rest
mass just by taking a logarithm. The physical energy identification requires
the same measured propagation and clock law. Conditioning must retain its
success probability, failure branches and changed subsequent state.

Addition, removal and neutral-transfer spectra also differ. A transfer can
measure a sum of charged excitation costs. Dividing its frequency by two
identifies a single charged gap only when the relevant equality is actually
established. A chemical-potential intervention can change addition/removal
costs while leaving a neutral-transfer frequency unchanged. That is an
identifiability control, not a change in units.

## Gauge, chirality and multiplicity boundaries

The finite representation decoder supplies internal charge/Casimir labels
under its declared spin/exterior lift. The common operator construction
puts the actual corresponding currents on the same localized radial modes
as the pixel-dependent mass kernel. Physical particles additionally require
their controlled propagation and physical spin/statistics. Equivalent carrier
representations on different factors do not identify physical spatial and
internal rotations. A valid conditional spin/matter tensor product likewise
does not prove that the continuing source realizes its physical action.

The staggered free extension retains two Dirac tastes and both chiralities;
the naive three-axis lattice has eight balanced Weyl nodes. Nambu particle/
hole redundancy is a bookkeeping doubling and must be distinguished from
these additional regulator species. Zero, dark and complementary modes must
remain in the census. An arbitrary projection or restriction of momenta is
not a decoupling argument.

The [regulator routes](../../model/common_source/chiral_candidate/REGULATOR_ROUTES.md)
include Euclidean overlap tests on nonconstant nonabelian backgrounds and a
stated locality domain. They do not establish the interacting chiral
measure, physical continuation, common state and observer correspondence.
The finite abelian quantum-link construction supplies a Gauss sector;
the nonabelian background transport calculations use supplied classical links.
Neither the internal anomaly sums nor a free regulator alone close the
physical nonabelian chiral theory.

## Conditions for a physical identification

These observation-specific acceptance tests expand nodes 3-8 of the
[unified graph](../../model/particle_joint_reconstruction/DEPENDENCIES.md)
for one proposed sector. Their source and measurement domains must be fixed
before any comparison.

1. Fix the physical response domain, allowed apparatus, reference state and
   clock. State separately which operations are native, compiled coherent
   services or supplied preparation/readout controls.
2. On that domain reconstruct an informationally sufficient response from
   complete records. Demonstrate ambiguity controls and predict withheld
   detector protocols, momenta and Higgs perturbations. Keep uncertainty and
   conditioning errors explicit.
3. Extract kinetic, mass and interaction response jointly. Retain dark/extra
   modes and invariant projector subspaces, distinguishing unresolved numerical
   clusters from exact multiplicities; check charge/current covariance
   and coordinate-invariant poles. Component checks on a supplied action are
   useful but keep that input visible.
4. Establish the physical chiral interacting regime, state and gauge constraints
   with controlled limits. Recheck what the detector can access in that regime.
5. Apply one fixed physical dictionary and a compatible matching/uncertainty
   contract to independent observations. A changed rule or matching map is a
   distinct candidate with its exposure to comparison data retained.

A positive finite test advances only its declared part of this chain. A
negative test identifies the indistinguishable alternatives, failed candidate
or exact excluded class. Identifying the measured particle spectrum requires
the field/record and numerical comparison contracts to hold in one continuing
construction.
