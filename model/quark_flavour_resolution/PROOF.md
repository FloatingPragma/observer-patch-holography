# Charge-sensitive quark spectra on the full native source

## Exact equal-coupling domain

The input is the full thirty-mode one-cell Hamiltonian and its five-state
quantum scalar, with all three primitive Yukawa coefficients and the scalar
occupation cost equal to one. The preparation is its actual native-current
singlet ground doublet divided by two. The complete elementary quark-removal
maps and spectral measures are those in
[native_quark_response](../native_quark_response/PROOF.md). Every transitive
input hash of that receipt is checked before using its exact rational
intervals and moments.

The up-conjugate and down-conjugate targets both have dimension 26 and the
same complete thirteen-root annihilator, characteristic polynomial, exact
ground moments and positive spectral weights. Equal moments and the common
annihilator give an exact cyclic isometry

`H_u^j c_u|g> -> H_d^j c_d|g>, j=0,...,12`.

Its Gram identity follows for every pair of powers, since the annihilator
recurrence determines all higher moments. Strict positivity at all thirteen
roots gives cyclic dimension thirteen. This is a spectral symmetry of this
particular source and preparation. Its extension to spatial preparations
uses the microscopic symmetries below.

This exchange cannot commute with every native gauge current. The literal
full fermion/scalar hypercharge current Q obeys

`[Q,H]=0`, `[Q,c_uc]=(2/3)c_uc`, `[Q,c_dc]=(-1/3)c_dc`.

A unitary commuting with Q cannot map one of these operators to the other.
The exact current algebra therefore permits charge-sensitive dynamics
without supplying separate quark Yukawa coefficients. The left-quark
removal charge is `-1/6`. All twelve global current commutators with Q
vanish; the scalar contribution to Q is retained. These statements supply
global covariance, without a dynamical gauge field or local Gauss theorem.

## Microscopic symmetry of the complete spatial law

The signed internal quark exchange and particle-hole map verified in
[quark_gauge_dressing](../quark_gauge_dressing/PROOF.md) sends

`c_quark(x,s,i) -> D_x epsilon_st S_ij c_quark(x,t,j)^dagger`.

Its internal antisymmetric pairs are `(0,12)`, `(1,13)`, `(4,14)`, `(2,3)`,
`(5,6)` and `(7,8)`. It exchanges the complete up/down interactions and
preserves the charged-lepton channel and every scalar coordinate. The
literal compact identity `h C+C conjugate(h)=0`, with `C=D epsilon`,
proves invariance of all spatial hopping and texture terms. It maps
`Y -> Y-B`, `B -> -B`, with fixed signs on the remaining eleven native
currents. Thus the global singlet B=0 domain is preserved. This operation
alone relates up removal to down addition.

The additional antiunitary `A c_i A^-1=i c_i^dagger`, with scalar
conjugation, preserves every local native interaction and reverses all
native currents. A general Hermitian spatial one-body term changes sign
under A, since `Tr h=0`. On the actual compact texture, the unitary
`R=(x -> -x) tensor i sigma_x` satisfies `R^dagger h R=-h` on all 686 rows.
Its determinant-one spin rotation preserves every Yukawa pair, and its
site reflection preserves all 343 scalar stiffness entries and the pair
stagger. Functional calculus therefore preserves Omega as well. Hence
`C_full=R A` is an exact antiunitary symmetry of the full joined law.

Composing C_full with the quark particle-hole operation maps up removal
to down removal, with a reflected site and a unitary spin frame. Its
charge map is `Y -> -Y+B`, `B -> B`. An invariant stationary density on
the global native-singlet B=0 domain gives identical complete up/down
removal spectra at reflection-fixed sites, including both compact cores,
or for matched reflected probes. In particular, the mixed lowest spectral
projector on that domain is invariant. The actual one-cell physical ground
has B=0. The baryon sector of the unrestricted spatial singlet ground is
not determined here, so that particular spatial ground is not silently
substituted.

The target baryon charge of either removal is `1/3`. Consequently the
same exchange sends `Y^2` to `(Y-B)^2`, which differs on the charged
targets. The neutral polarization term transforms to `(-Y+B) X_A`;
changing the memory sign cannot remove its B term. This identifies the
actual symmetry broken by both charge-sensitive source controls.

## Four fixed unit-coupling controls

The same exact bare energy measure is transferred by four explicit laws.
No coefficient is selected using the measured up/down ratio.

For `H+Q^2`, a target of charge q has an energy shift `+q^2` and unchanged
weights. The unchanged neutral source preparation is a ground in its
original Q=0 domain. This positive source-current term makes the selected
up record higher than the down record.

For `H+4N_ref+Q Q_ref`, the actual compensating reference has charge `-q`,
so the interaction shift is `-q^2`. The bare source energy is unchanged,
the declared static reference cost is four, and the cross interaction has
energy `-q^2`. Detaching that interaction restores the equal bare spectra.
The shift is an observer-reference binding energy, without an intrinsic
particle-mass identification. No minimum over all charged/reference sectors
is inferred from the neutral source ground.

The positive completion `H+4N_ref+(Q+Q_ref)^2/2` has source and reference
self energies `q^2/2` each and cross energy `-q^2`. Its total charge-field
energy on the compensated sector is zero. Keeping only its attractive
term would change this law. This square is an algebraic control; a local
Maxwell or quantum-link implementation is outside the calculation.

## Intrinsic neutral polarization

Adjoin a neutral two-state memory to the source itself and declare

`H_pol=H tensor I2+I tensor N_A+Q tensor X_A`,
`N_A=diag(0,1)`, `X_A=[[0,1],[1,0]]`.

Every coefficient is one. This is an explicit source architecture, not a
charge-dependent fitted offset. The preparation is the same Q=0 native
ground times the same memory blank `|0>`. The compensating charged reference
retains its original representation and cost four; it does not couple to
the polarization memory. The source dimension is `10*2^30`. Its Hamiltonian
is bounded on the complete carrier: the native scalar-star bound gives
`||H||<=29`, the actual hypercharge norm is at most `13/2`, and hence
`||H_pol||<=73/2`. No occupation sector is discarded to obtain this bound.

On the Q=0 domain, the added memory has Hamiltonian N_A and the native
ground is unchanged. This is a ground statement on that domain, without
a minimization over all charged domains. On a complete charged target,
conservation of Q gives the exact tensor Hamiltonian

`H_target tensor I2+I tensor M_q`, `M_q=[[0,q],[q,1]]`.

There is therefore an explicit intertwiner from the complete native target
times the memory to the extended source target. All native target energies,
their degeneracies and their zero residues survive. Let
`delta_q=sqrt(1+4q^2)`. Each native atom at energy epsilon and weight w gives
two joint atoms

`epsilon+lambda_-, lambda_-=(1-delta_q)/2`,
`epsilon+lambda_+, lambda_+=(1+delta_q)/2`,

with original weights

`w (1+1/delta_q)/2`, `w (1-1/delta_q)/2`.

The blank-memory resolvent and exact Schur selfenergy are

`G_q(z)=[z-H_target-q^2/(z-H_target-1)]^(-1)`,
`Sigma_q(w)=q^2/(w-1)`.

Thus the split is a source response involving an actual neutral degree of
freedom. Both branches are required by its spectral law. In an individual
memory eigenstate the receipt also gives the expectation of N_A and QX_A;
their sum is lambda. These are subsystem expectations of noncommuting
terms, without simultaneous sharp energy records for those terms.

The exact memory identities are `<0|M_q|0>=0` and
`<0|M_q^2|0>=q^2`. Since every original removal measure has weight one half,
the complete unconditioned first transfer moment is unchanged and the
second moment increases by `q^2/2`. The up-minus-down second-moment change
is exactly `1/6`. The different lower poles do not imply a change of the
whole spectral mean. At zero coupling the upper-memory residue vanishes
and the complete original measure is recovered, with the zero branch
retained explicitly.

The polarization interaction is number-conserving in the fermions. Every
literal native annihilator obeys `[c_i,Q]=Y_i c_i`, so

`{[c_i,Q tensor X_A],c_j}=Y_i {c_i,c_j} tensor X_A=0`.

The executable checks all 900 such components. The full anomalous Higgs
contact tensor therefore remains the same unit-coefficient tensor. The
ground scalar Gram is also unchanged. The extra contribution to its
retarded first moment contains one X_A, whose expectation in the blank
memory is zero. Its normalization remains
`Z_H=2-7n/4-3e/4`, approximately `9.58523463744285`, on this unchanged
neutral ground. Higher scalar moments and poles need not be unchanged.
This is a charge-dependent propagator selfenergy with a universal primitive
Higgs contact, without an inserted family Yukawa matrix.

There are 26 joint energy records for each singlet channel and 50 for the
left-quark channel, on complete target dimensions 52, 52 and 104.
Every zero residue is retained. Each channel also keeps
its original failure probability one half and the zero outside-target
outcome. Failure records do not measure source energy, so no zero physical
failure energy is assigned. The summed successful probability is one half
and is never divided out. All nonzero-weight success energies in these
unit controls are strictly positive by exact rational interval checks.

## A fixed one-scale comparison

The bare dominant gap is approximately `0.7815736835381646`. The positive
Q-squared control gives up/down ratio `1.3734054118`. The compensated
reference interaction gives ratio `0.5028308111`, with bare source and
reference energies explicitly separated. The polarization memory gives
up poles `0.4482403502048313` and `2.114907016871498`, with original
attempt weights approximately `0.3965334900670651` and
`0.0991333725167663`. Its down poles are `0.6806484709608330` and
`1.8824988961154963`, with weights approximately `0.4540433107451119`
and `0.0416235518387195`. These four values refer only to the original
dominant atom; the other complete records are retained. The lower-pole
ratio is `0.6585489710600182`.

The comparison uses the frozen running quark parameters in
`common_source/calibrated_matter_results.json`, at their stated common
renormalization scale. Their ratio is approximately `0.4574566638`. Exact
positive interval intersections test one shared energy conversion under
the specified record-zero/lower-pole dictionary. No control fits both
fixed parameters by one scale. In the reference-attraction case, anchoring
the up value gives conditional down value `0.002457966237 GeV`, compared
with target `0.002701766647 GeV`. This is a failed dictionary comparison,
without fitted primitive coefficients, propagated experimental uncertainty
or a physical rest-mass assertion.

The construction proves that the equal-coupling local degeneracy is not
an unavoidable equality of all native-current responses. It supplies
complete dimensionless charge-sensitive spectra. Physical quark masses
require a corresponding propagation, renormalization and observable
identification on the chosen source.
