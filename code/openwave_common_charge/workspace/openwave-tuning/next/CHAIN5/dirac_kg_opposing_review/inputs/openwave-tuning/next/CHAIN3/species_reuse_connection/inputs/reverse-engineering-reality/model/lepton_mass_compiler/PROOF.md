# A complete lepton analyzer compiled into free Dirac mass records

## Literal source and its entire unitary

The address, character and analyzer data are those of
[lepton_cycle](../lepton_cycle/PROOF.md). Two independent sixty-address
registers give a 3,600-dimensional source space. Their prescribed blank
state is S=u tensor u, where u is uniform. The actual native record partition
has fifty blank, five positive and five negative addresses in each cell.
The three analyzer masks have cardinalities 500, 250 and 250. Let phi be
the sum of their normalized indicator vectors divided by sqrt(3).

Let H_S be the tensor product of the two real Householder maps taking the
first address to u, and let H_phi take the first pair address to phi.
Let R be the product of the two reflections about u. Both are involutions,
and RS=S. On character k, D_k multiplies positive-X addresses by
`exp(i q_k)`, negative-X addresses by `exp(-i q_k)`, and all other addresses
by one, with `q_k=theta+2pi k/3`. The fixed phase is theta=2/9. Define the
complete unitary on all 10,800 address/character coordinates by

```
U=direct_sum_k (H_phi D_k R H_S),
Q=|address0><address0| tensor I3,
W=U^dagger Q U.
```

Every factor acts on every input coordinate. This is the unitary completion
of the actual analyzer kick, since its pullback projector is
`R D_k^dagger |phi><phi| D_k R`. No postselected returned-source map is
substituted for U. The code executes U^dagger U and W on every input column
in bounded batches, without allocating the dense square matrix. Unitarity
also follows exactly from the declared factors, independently of the
floating execution residuals.

The blank-to-address-zero block and its positive compression are

```
K00=diag(lambda_k),
lambda_k=sqrt(5/108)+2sqrt(5/216) cos(theta+2pi k/3),
Y=K00^dagger K00=diag(e_k), e_k=lambda_k^2, Tr Y=5/18.
```

The constants follow from sqrt(500)/(60sqrt(3)) and
sqrt(250)/(60sqrt(3)). The three root-of-unity sums prove the trace exactly.
Rational sine remainder bounds and integer-square-root enclosures certify
that every lambda is positive, and `e_1>1/20000`. The numerical coefficients
in ascending order are

```
e_1=0.0000753756987638,
e_2=0.0155854570501,
e_0=0.262116945029.
```

The phase is an inherited coherent program input, with the same retrospective
architecture and character assignment as the lepton cycle. Its coincidence
with the separate weak-current rate 2/9 is not asserted to select a phase
from passive probabilities or elapsed native-attempt time.

## Exact map to the configured charged-family matrix

The independently configured matrix in
[KOIDE_SOURCE](../standard_model_source_map/KOIDE_SOURCE.md) is
`C_ex=sqrt(10/72) I+sqrt(5/72)(z P+conjugate(z) P^2)`, where P is its
actual cyclic permutation. The witness obtains its constant, cosine and
sine coefficient matrices from the existing `exact_source` evaluations at
`z=1,-1,i,-i`. It checks equality to that declared matrix polynomial before
evaluating it at `z=exp(2i/9)`. The existing executable's rational-phase
restriction is respected; this exact polynomial extension does not define
the configured matrix from the native effect.

Let `omega=(-1+i sqrt(3))/2` and `F_jk=omega^(-jk)/sqrt(3)`. Direct exact
matrix arithmetic checks `F^dagger P F=diag(1,omega,omega^2)` and

```
F^dagger C_ex F=sqrt(3) diag(lambda_k),
F E_native F^dagger=C_ex^dagger C_ex/3.
```

The identity is first verified as a polynomial in independent real cosine
and sine coordinates, then evaluated at theta=2/9 and compared to the
actual complete lepton-cycle analyzer. The existing native geometric
frame O is imported and checked through `R O=O P`; the geometric witness
is O F. Gauge currents are identity on this same three-family factor,
so the frame change preserves their representation.

The actual prepared double-load effect has trace 5/18. The separate native
pair preparation has trace 2/9, and its three squared amplitude parameters
are each smaller by 4/5. The executable retains this exact 5/4 preparation
factor. Negative controls reject replacing Y by Y squared, substituting
the native-pair effect or dropping the 5/4 factor, dividing out the active
probability, and conjugating the Fourier frame with the mode labels fixed.
The identity joins the coefficient and its family frame. The configured
eight-flag apparatus and native 3,600-address apparatus have different
complete source laws; the full-CAR implementation below uses the literal
native analyzer. The receipt separates this exact admission from physical
clock, Higgs normalization and observed-particle identifications.

## Full-CAR interaction, before compression

Use the same ninety physical modes per spatial site,
`C2_spin tensor R15 tensor C3_character`, as the native charged-Higgs
attachment. Give every one-particle mode the full 3,600-address coordinate.
The enlarged carrier has 324,000 fermion modes per site. Physical modes
occupy address zero; the other 323,910 modes are supplied in their vacuum.
Let J be this vacuum embedding on the ENTIRE physical Fock space and
P=JJ^dagger. The physical source has every fermion occupation, rather than
only a one-particle or low-occupation input.

The exact full-CAR theorem in
[CAR_KOIDE](../charged_lepton_response/CAR_KOIDE.md) applies to a unitary
analyzer with this rank-three Q. Its hypotheses are checked on this literal
U and the actual native charged-Higgs support. Use U on the entire native
lepton doublet and conjugate(U) on the conjugate singlet, with identity on
other internal/spin factors. Denote this one-particle unitary by S_U and
its exterior lift by Gamma(S_U). The native doublet and singlet projectors
commute with all twelve currents. The two Higgs components retain their
covariance; fixing a neutral scalar orientation preserves only its unbroken
electric charge, not all twelve generators.

Let D be the actual directed neutral charged-Higgs tensor with unit singular
normalization. Declare the neutral background norm w_H=1. In the usual
parameterization w_H=v/sqrt(2), this gives v=sqrt(2) in source units.
The background is supplied; it is neither an empirical weak-scale input
nor a selected interacting vacuum expectation. In the source's
creation-pair convention, set

```
Delta_Q=epsilon_spin tensor
         [conjugate(D) tensor Q+conjugate(D)^T tensor Q],
V=Gamma(S_U)^dagger H_(Delta_Q) Gamma(S_U).
```

The exact CAR identity transforms a pair kernel as
`Delta -> S_U^dagger Delta conjugate(S_U)`. Its oriented doublet/singlet
block is therefore `U^dagger Q U=W`, and its reverse is W^T. Using U on
both legs would instead produce `U^dagger Q conjugate(U)`, a different full
operator. The implementation tests that difference on the actual source
frame. No identity between a second-quantized sum and a sum of
second-quantized unitaries is invoked.

Vacuum compression kills every monomial containing an ancillary creator or
annihilator. All remaining terms give the full operator identity

```
J^dagger V J=H_Y,
Delta_Y=epsilon_spin tensor
          [conjugate(D) tensor Y+conjugate(D)^T tensor Y^T].
```

This equality holds on all `2^(90V_sites)` physical states. The uncompressed
kick has six disjoint two-mode pairing generators per site, each of norm
one. For V_sites sites, `||V||=6 V_sites`. Its mean has norm
`2 V_sites Tr Y=5 V_sites/9`. These are full-Fock norms.

The physical-vacuum subspace is NOT reducing under V. Its leakage is kept
in the subsequent instrument. For a useful exact finite control, each
character block of W has rank one. The span of the physical flag vector
and W times that vector has dimension two, because 0<e_k<1. In a suitable
ancillary-only coordinate frame its matrix is

```
[[e_k,sqrt(e_k(1-e_k))],
 [sqrt(e_k(1-e_k)),1-e_k]].
```

All other address directions are exact zero-coupling spectators. This gives
a sixteen-state Fock control for one spin-oriented physical pair and its
two ancillary modes. The executable retains every physical occupation,
all four ancillary outputs, all 256 histories of four rounds, and an
entangled physical/reference input. A flag-independent number Hamiltonian
is included; it commutes with the analyzer and fails to commute with the
pairing. Its full reduced dynamics, compression and archive error are
checked independently. This numerical active-span control does not claim
execution of the enormous precision examples below.

The 256 fine labels specify archive coordinates. The production instrument
records only the sixteen vacuum/nonvacuum histories. Inside each such
history the output is the coherent sum of its Kraus operators tensored
with their orthogonal ancillary history kets. The receipt gives both the
fine diagonal probabilities and the sixteen coarse probabilities. In the
complete trace comparison only the all-vacuum block overlaps the ideal
output; every other coarse block is positive, so its trace contributes
its full probability with all internal coherences retained.

## Common kinetic law, archive and clock

At finite periodic volume, put the same free nearest-neighbor Weyl hopping
on every family and flag copy, with hopping c/(2a), c=1. The extended
Hamiltonian H0 preserves the ancillary vacuum and is identity on the
family/flag coordinates, so it commutes with S_U. Its physical restriction
is exactly the kinetic term used with H_Y. No kinetic operator or clock
is exchanged at the compression.

A round applies `exp(-ih(H0_extended+V))`, records vacuum versus nonvacuum
ancillary occupation, and exports the entire ancillary output before
supplying fresh vacant modes. The archive retains every charge, coherence
and failure. Only the invariant coarse vacuum flag is measured; no
nonabelian fine-occupation covariance is assumed. The input can be any
physical Fock state with any untouched reference. In particular, it may
be the filled negative band of the target quadratic law. The supplied
ancillary vacuum need not be the ground of H0, and the complete controller
preparation is not claimed to be a stationary interacting ground.

Let N rounds have h=t/N, and L=6 V_sites. In the interaction picture of
H0, the compressed first-order integrals agree with those of H_Y at every
time. Integral remainder and contraction telescoping give the same
complete-archive theorem as CAR_KOIDE:

```
Pr(any nonvacuum flag) <= L^2 t^2/N,
||actual complete output - ideal target with vacuum archive||_1
    <= 3 L^2 t^2/N =108 V_sites^2 t^2/N.
```

The comparison includes arbitrary retained references and the quantum
complements of every flag. No conditioning on a successful reset occurs.
The bound uses committed classical vacuum/nonvacuum flags. Before this
specified measurement, a wholly coherent flag dilation has the weaker
unhalved trace bound `2sqrt(2) L |t|/sqrt(N)`. Retaining each measured
sector's quantum complement does not remove this distinction.
A certified complete dilation error eta per round adds 2N eta to the
trace bound. These bounds require no factor ||H0|| because their Dyson
remainders are taken in its interaction picture.

The sole model clock is accumulated pulse duration t. Analyzer preparation,
uncomputation, measurement and archive transfer are specified services
with separate overhead. Equality with elapsed native-attempt or laboratory
time is not supplied. At each fixed finite volume and time, N tending to
infinity gives the same target unitary generated by H0+H_Y. At finite N,
its noisy recorded evolution has no asserted exact target energy spectrum.

## Dirac propagation and the same rest energy

The native spin-pairing and Clifford identities of
[relativistic_mass](../charged_lepton_response/RELATIVISTIC_MASS.md) apply
to this derived Y, without changing its normalization. On all 180 Nambu
coordinates,

```
h(p)=-sigma.p tensor I45,
K=D tensor conjugate(Y)+(D tensor conjugate(Y))^T,
Delta=epsilon_spin tensor conjugate(K),
H(p)=[[h(p),Delta],[Delta^dagger,-h(-p)^T]],
H(p)^2=diag(|p|^2 I+Delta Delta^dagger,
            |p|^2 I+Delta^dagger Delta).
```

Thus the three charged Dirac fields per valley have rest energies
`m_j=w_H e_j=e_j`. The coefficient is not squared again. Their positive
Nambu rows include two spin states and conjugate charge; Nambu duplication
adds no species. The complete native symbol also retains thirty-nine zero
internal coefficients of this charged-only free action.

With `p_i -> sin(a p_i)/a`, all eight lattice valleys are retained. For
`|p_i|<=P_box`, the exact symbol error is at most
`delta_a=a^2 P_box^3/2`; unitary error is at most |t|delta_a. For the charged
negative-band reference, the two-point bound is
`delta_a(|t|+1/m_*)`, using the certified `m_*=1/20000`. All eight valley
orientations are executed on the complete native symbol. The resulting
free continuum Dirac equation has its standard hyperbolic propagation.
This is a controlled free-field construction, without a chiral mirror
removal or an interacting gauge continuum theorem.

## Complete mass observation and finite precision

The rest observation uses the reservoir-SWAP and filled-reference grammar
of [mass_instrument](../coherent_defect_wave/mass_instrument/PROOF.md), with
its Hamiltonian reconstructed for this source's mass m. No Wilson rest
parameter or spectrum is copied from that package. For one charged Dirac
block at p=0, choose the unit native charged supports and write their
annihilators as `(l_up,l_down,r_up_conjugate,r_down_conjugate)`. The literal
canonical transformation

```
b=(l_up,l_down,r_down_conjugate^dagger,-r_up_conjugate^dagger)
```

turns `m(l_up^dagger r_down_conjugate^dagger
-l_down^dagger r_up_conjugate^dagger+h.c.)` into `b^dagger(m beta)b`.
Its vacuum is the state with the two conjugate right modes occupied.
The code builds all sixteen Fock columns of this transformation and checks
every annihilator, the pairing Hamiltonian and electric charge exactly.
Native electric charge is `2-N_b`; the record's relative Dirac number is
`N_b-2`, with its negative electric charge recorded separately.
Thus `h=m beta` and `P_-=(I-beta)/2` belong to the same pairing source.
The negative-band Slater
projector is constructed by its actual second minors. Its complete
sixteen-state Fock Hamiltonian relative to the sea has energies
`0,m,2m,3m,4m`. The code checks the complete polynomial energy projectors,
all nine nonzero energy/charge effects, the 32-state fermionic SWAP, and
all 36 setting/reservoir/energy/charge records. Thirty records have zero
weight in the preparation and are retained. A successful single-fermion
addition or removal reads m; a failed probe can create a neutral pair at
energy 2m. This distinguishes the rest mass from a pairing-sector gap.
The reservoir mode has electric charge -1. An idle reference of opposite
charge is prepared with the same initial occupation and retained. Every
nonzero source/reservoir/reference record has zero total electric charge
and even total parity. This compensation leaves the thirty-six original
probabilities unchanged; the controlled probe's work is supplied apparatus.

The exact spectral projectors are a declared ideal detector service of
the limiting target. A separate finite-duration Ramsey construction gives
an explicit error transfer for the compiler. Prepare
`(|ground,reference blank>+|positive-energy particle,dual reference occupied>)/sqrt(2)`.
It has fixed even total parity and zero total electric charge. The reference
is idle and retained. At p=0 the ground and particle eigenvectors of beta
are independent of the positive coefficient m; the preparation and the
readout effects do not encode its value. On that two-dimensional span let X be Pauli X and
`Y=i(|0><1|-|1><0|)`; extend them by zero on its orthogonal complement.
The complete four effects `(I±X)/4,(I±Y)/4` include the recorded choice of
X or Y with probability one half. All complement states are retained.
The exact limiting probabilities are

```
(1+cos(mt))/4, (1-cos(mt))/4,
(1+sin(mt))/4, (1-sin(mt))/4.
```

If the complete compiler output error is epsilon in unhalved trace norm,
each Pauli expectation changes by at most epsilon. The complex contrast
changes by at most sqrt(2) epsilon. For epsilon<=1/4 its phase error is at
most 4 epsilon, so the energy readout `atan2(<Y>,<X>)/t` has additive error
at most `d=4 epsilon/t` on the specified small-angle branch. The four
original probabilities each have error at most epsilon/2. This statement
uses the whole source/reference state and its fixed half-probability
settings, without conditioning on a rare successful fermion event.

A ratio can only be reported after resolving its denominator. With exact
source enclosures `[l_i,u_i]`, the implementation requires
`d<1/20000<=m_e` and returns

```
[(l_i-d)/(u_e+d), (u_i+d)/(l_e-d)].
```

At t=1, choosing N=`108 V_sites^2 *10^12` gives epsilon<=10^-12 and
`d<=4*10^-12` in source-energy units. One-site rest readout consequently
requires 108,000,000,000,000 rounds under this conservative bound. A
side-four periodic example has V_sites=64, 20,736,000 enlarged modes per
round and 442,368,000,000,000,000 rounds. The receipt gives the full exported
ancillary-mode count. These are proved resource formulas, not an executed
simulation at that round count or a practical apparatus claim. The
intervals concern exact computed probability laws. Finite-shot noise,
physical detector error and interacting-model error are separate inputs.
Multiplication of the guarded ratio intervals by the frozen electron
anchor gives corresponding GeV intervals. Their error covers the finite
compiler record law; uncertainty of the empirical anchor is separate.

## One empirical energy unit

Only after the source coefficients are computed, use the frozen electron
energy 0.00051099895069 GeV as one scale anchor. The source energy unit is
approximately 6.77935938333 GeV. The three same-source free rest parameters
then become

```
electron: 0.00051099895069 GeV,
muon:     0.1056594144959 GeV,
tau:      1.7769849708134 GeV.
```

The electron value is a calibration input. Muon and tau data do not enter
the source or its conversion factor. The mode labels, phase 2/9 and
three-character architecture retain their retrospective input ancestry.
Converting the unit background literally gives w_H=6.77935938333 GeV;
no empirical weak-scale matching is supplied. No dynamical Higgs kinetic
residue is derived, so Y=e is the compiled pairing coefficient, without an
identification as an absolute canonically normalized Standard Model
Yukawa coupling.
The code proves a concrete effect-to-free-mass dynamics for this declared
compiler. It does not identify these three free fields with observed
interacting leptons, derive their widths, select a Higgs vacuum, eliminate
the other valleys or calibrate a passive native clock.
