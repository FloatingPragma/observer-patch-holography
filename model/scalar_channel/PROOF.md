# Native scalar capacity and current domains

## Source and scope

The law is `A-native-scalar-overlap-v3`, law 2 of the
[construction manual](../construction_manual/CONSTRUCTION_MANUAL.md): the
scalar channel written on ordinary registers. The manual distinguishes this
scalar law from the separate load-and-memory law. No moving charged-owner
Hamiltonian or interacting seven-species dynamics transfers through this
scalar compilation.

[gauge_capacity.py](gauge_capacity.py) supplies exact finite dimension and
matrix witnesses. Its 49-state owner frame and eight-state fermion frame
are independent small counterexamples, not dense executions of the complete
432-bit compiler. The full compiler uses local Pauli words, stabilizer
constraints and explicit fusion-path columns.

## One parity sector and its complete operator algebra

The source has 150 scalar modes: 30 seams and 120 wedges. Twelve fuel and
twelve spectator reference modes give 174 modes in the fixed instrument
graph. Its 432 interaction edges carry 432 encoded bits. For this connected
graph the cycle rank is

```text
432 - 174 + 1 = 259.
```

There are 260 local penalty cycles, comprising 240 star triangles and 20
cell faces. Their rank is 259. The all-positive cycle code therefore has
dimension `2^(432-259) = 2^173`. The graph number operators satisfy
`product_i B_i = I`, since each edge Z occurs twice. Under
`B_i = 1 - 2 n_i`, this fixes total fermion parity to even.

For n modes the full Fock space has dimension `2^n`; each parity space has
dimension `d = 2^(n-1)`. The abstract even CAR algebra is

```text
CAR_even(n) = M_d(C) direct_sum M_d(C).
```

Its restriction to even total parity is one `M_d(C)` summand. It is complete
on that sector and is faithful there. It is not faithful on the direct sum:
the nonzero central projection `(I - parity)/2` maps to zero. The executable
three-mode witness checks that projection on all eight occupations and its
restriction to the four-dimensional even space. Thus “complete even
algebra” in a compilation claim must name the selected parity sector.
Indeed the even monomials span exactly the matrix units between occupation
words of equal parity, giving the two full matrix blocks directly.

Both fuel inputs are available in this same architecture. Start with 120
occupied source modes and set the selected fuel and its decoupled reference
both to zero or both to one. Total occupation is respectively 120 or 122.
Both are even. Number-conserving source/fuel pulses preserve total parity;
they can change source parity by exchanging a particle with fuel. Keeping
the reference allows these experiments without representing the odd total
sector or identifying the full even algebra with its one-sector quotient.

## Capacity of the owner code

One owner is the seven-state space consisting of vacuum, a color triplet
and a Cartesian weak triplet. Its neutral subrepresentation is
`spin0 direct_sum spin1`, of dimension four. This is the tensor product of
two virtual spin halves. On r owners its collective spin-zero multiplicity
is the Catalan number

```text
C_r = binomial(2r,r)/(r+1).
```

The construction selects the first `2^r` ordered singlet fusion paths. This fits
at the declared degrees: `C_10 = 16796 >= 1024` and
`C_15 = 9694845 >= 32768`. Distinct fusion paths are orthogonal by the
Clebsch decomposition. Their explicit coefficient tensors give an isometry
`V_r`; path zero is the product of owner vacua. Arbitrary code unitaries
have the complete raw extension `V_r U V_r* + I - V_r V_r*`.

This construction stores r scalar logical bits collectively. It does not
turn a seven-state owner into seven independent fermionic modes. Seven
unrestricted modes require 128 states, or seven coherent bits, per site.
Even a fixed parity sector requires 64 states. Neither embeds isometrically
in a seven-state owner or in a seven-state owner with one additional bit
(dimension 14). A seven-state alphabet permits vacuum plus six single
excitations; adding a seventh species gives eight states even before
multiple occupancy is admitted. More hardware or an explicit constrained
sector is necessary for a seven-species matter construction.

## Raw currents and the accepted scalar response

Use the inherited owner representation in the displayed vacuum/color/weak
coordinates. Eight traceless Hermitian color matrices, three Cartesian weak
rotation generators, and the color identity span the twelve-dimensional
compact algebra `u(3) direct_sum so(3)`. Its derived algebra has dimension
eleven and its center dimension one. The exact matrix calculation verifies
these ranks and bracket closure. These are coordinates for the declared
representation; they do not reconstruct a response menu from scalar records.

For each current Q, let `Q_total = sum_j Q^(j)` on a star. The selected
fusion paths obey

```text
Q_total V_r = 0
```

for all twelve currents. Color and the central color current vanish on the
neutral owner subspace. The weak total vanishes on its spin-zero paths.
Consequently the old current action and every commutator response it induces
are zero on the accepted scalar code. The twelve-current representation is
faithful on the raw owner space, which contains charged complement states.
Faithfulness on that larger space cannot be asserted for scalar-code states.

The complete 49-state two-owner witness uses code columns

```text
|vacuum,vacuum>,
(|x,x> + |y,y> + |z,z>)/sqrt(3).
```

All twelve collective generators annihilate these columns. A color state
outside the code has central charge one and is unchanged by the complete
lift of a logical X gate. A weak generator acting on only one owner has
squared norm `2/3` on the second codeword. Neutrality therefore uses the
shared star frame. Independent owner rotations are different operations.

## Number symmetry, relative phase and the intrinsic-response boundary

An encoded scalar number N is nonconstant on the scalar code and commutes
with every old collective current. If N is independently admitted into the
current list, its code response adds one direction to the zero old code
response. On a raw representation retaining the faithful old currents, it
adds a thirteenth independent direction. The resulting specified algebra
has center dimension two, rather than one. Independence follows immediately
by restriction to code: every old current vanishes there while N does not.
The two-owner witness realizes this extension by the projector onto its
second singlet column, whose code compression is `diag(0,1)`.

Whether that direction belongs to an intrinsic gauge response requires an
operational premise. Conservation of total number alone does not supply it.
For the source and apparatus jointly, the free law and the declared hopping
pulses commute with total number. On a fixed-number input a global phase
`exp(i theta N_total)` is scalar and changes no outcome law. A mixture of
number sectors with no coherence between sectors has the same invariance.

A phase of source relative to fuel is operationally different. The exact
three-mode witness starts in `|source,fuel,reference> = |0,1,1>` and uses a
balanced source/fuel splitter, a phase, and the inverse splitter. Every
branch has even total number two and reference occupation one. A zero phase
returns occupied fuel with probability one; a source phase pi returns it
with probability zero. A global total-number phase pi/2 reproduces the zero
phase records. All eight terminal occupation outcomes are included, with
their zero probabilities retained. The source-number phase does not commute
with the converter; it is a phase relative to the apparatus.

These facts distinguish a representation, a global conserved number, and a
controlled relative phase. They neither identify all admissible intrinsic
currents nor prove the twelve-dimensional faithful complete-response premise
A1 of the gauge theorem. That premise cannot be transferred from the old
charged response to a neutral scalar code by adding a phase control. The
compiler, its constrained scalar capacity and its declared instruments have
the finite domains above; interactions and the complete intrinsic response
of a larger matter process require their own construction.

## Routed loop instruments and retained preparation records

The routed initialization has 432 data bits and twelve neutral scratch bits,
one at each native center. The following instrument identity holds for an
arbitrary data state and arbitrary untouched quantum reference, on the
accepted owner codes and with all scratch bits initialized to zero.

Write a Hermitian cycle Pauli as `C = product_b C_b`, with one local Pauli
factor for each supporting star. Assign its overall sign to one factor.
Starting at a declared center, prepare its scratch bit with H. Apply the
controlled local factor at each participating center. Neighbor SWAP gates
move the scratch state along the declared route; each destination scratch
is zero before the transfer, and the source scratch becomes zero after it.
A closed route returns the scratch state to its starting center and leaves
every other scratch bit zero. This compiles the joint controlled operator

```text
|0><0| tensor I + |1><1| tensor C.
```

Apply H to the returned scratch and measure its local Z. For reported sign
`s = +1` or `-1`, the induced Kraus operator on data is exactly

```text
K_s = (I + s C)/2,        sum_s K_s* K_s = I.
```

All routing SWAPs and controlled factors are part of the unitary dilation;
no global Pauli measurement is a primitive. For data/reference density rho,
the unnormalized output is `(K_s tensor I_R) rho (K_s tensor I_R)`. The sign
is retained in a separate classical record. The local scratch reset is an
outcome-controlled X after that record is written. It does not erase the
record or reset the data. The identity therefore transfers the complete
instrument, including coherence and reference correlations, rather than
only expectation values on one preparation.

Choose the spanning-tree fundamental cycles `C_1,...,C_259`. They commute.
Each has its own non-tree chord, which no other fundamental cycle uses.
The chord Z anticommutes with its own cycle and commutes with the others
and with all number operators B_i. For sign string s, let F_s be the
product of these chord Z operators at negative signs. The complete corrected
instrument is

```text
rho -> sum_s |s><s| tensor F_s P_s rho P_s F_s*,
P_s = product_j (I + s_j C_j)/2.
```

Tensor identities on any untouched reference are implicit. Corrections can
be applied immediately after their measurements because they commute with
every other cycle projector. Summing all sign records gives a trace-
preserving channel. The routed implementation executes these same maps
using only local scratch reads and neighbor transfers.

Scratch diagnostics extend this instrument to arbitrary scratch inputs.
Let `Z_b^+` and `Z_b^-` be the two local scratch projections. The diagnostic
outcomes are a first negative sign at center j, with projector
`D_j = (product_{b<j} Z_b^+) Z_j^-`, or all positive signs, with projector
`D_ok = product_b Z_b^+`. These orthogonal projections sum to identity.
Every first-negative branch retains `D_j rho D_j`, its probability and its
classical diagnostic record, and halts without reset. The all-positive
branch has blank scratch and proceeds through the loop instrument above.
This defines a complete instrument for coherent dirty scratch and arbitrary
untouched references. The simulator evaluates exact deterministic scratch
reads without changing the tableau; uncertain reads use the local Z
projection. No global predicate discards a coherent dirty input.

For the specified all-zero data-bit initialization, every nonempty product
of distinct fundamental cycles has an X on at least one unique chord.
Its expectation in the initial Z-basis word is zero. Expanding P_s therefore
gives probability `2^-259` for every sign string, without enumerating or
postselecting those strings. The 174 number signs are initially positive
and remain positive. After feedback all cycle signs are positive too.
The number and cycle stabilizers have joint rank 432, so their common
positive eigenspace is one-dimensional. Every sign branch reaches the same
encoded vacuum density while retaining its distinct preparation record.
Uniform signs and the unique vacuum conclusion use this particular blank
initialization; the Kraus identity itself has the arbitrary-state domain.

The executed finite program routes 259 loop reads through 941 neighbor
transfers and then uses 120 local even pair operations to occupy the 120
source transit modes. Before each loop read it records twelve local scratch
diagnostics, giving 3,108 diagnostic reads. These are certain positive signs
on the stated blank preparation, so the 259 loop signs keep their uniform
law. The sign lists, reset actions and final stabilizer state are retained. An
omitted chord correction produces a wrong syndrome whenever its measured
sign is negative, a failure of the specified prepared-code condition. A
missing or wrong acknowledgment halts the protocol with a failure record
and preserves the actual partial data/scratch state, including its
entanglement. It does not substitute the ideal loop projector, discard that
branch, or reset the incomplete apparatus silently.
The acknowledgment controls supply a classical fault at the first
transport. Conditional on that supplied fault and positive diagnostics,
the halted quantum map is the actual unitary prefix of the route. These
controls do not implement a noisy communication channel or infer a fault
from an unspecified message distribution.

The controller, the preparation of blank scratch, the local gates and the
finite record allocation are supplied parts of this finite process. The
routed instrument theorem does not construct an autonomous controller,
unbounded fresh memory, or a work balance for repeated reset and reuse.

## From routed preparation to finite scalar instruments

Let V be the isometry from the even-total-parity Fock space of 174 modes
to the positive-cycle graph code, followed by the declared star isometries.
The number signs and cycle signs specify a unique occupation basis state.
Consequently the routed output with all 120 transit modes occupied and all
other modes empty equals V applied to that occupation state, up to a global
phase. It is the Slater state with occupied frame
`F_0 = [0_(30 by 120); I_120]`. Preparing fuel and its parity reference
together in 00 or 11 preserves even total parity and changes no source
occupation. The reference is then decoupled from the source and converter.

For each allowed one-body matrix h, the number and hopping identities give
`H_compiled V = V dGamma(h)`. Hence their exponentials intertwine, as does
every finite product of the declared ramp steps, pulses and phase settings.
The ramp executes the occupied-frame update `F -> exp(-i delta h) F` at
each stated midpoint. Its output covariance is `C = F F*`; this calculation
does not replace the prepared state by eigenvectors of the final source.
The final eigensystem only measures its distance from the exact filled sea.
The covariance residuals are numerical checks, not interval error bounds.

The executed helper accepts a pure Slater source with even particle number.
It checks that C is Hermitian, has spectrum in [0,1], is a projector, and
has even integer trace, within its declared numerical tolerance. Both paired
fuel/reference inputs then lie in the represented parity sector. A general
fractional covariance, or an odd-rank source projector with these paired
inputs, is outside this preparation domain. An odd source with an odd
number of encoded apparatus particles can lie in the even-total-parity
space, but that is a different preparation. An untouched reference outside
the encoded carrier does not change the carrier's parity constraint.

Write U for the product of the two source-on one-body converter pulses,
the free interval and the selected fuel phase. Let Gamma(U) denote its
number-preserving Fock-space unitary. With blank pointer bits, the terminal
two-pointer isometry is

```text
psi -> |00> tensor (I - n_f) psi + |11> tensor n_f psi.
```

Its complete instrument has Kraus operators
`K_00 = (I - n_f) Gamma(U)`, `K_11 = n_f Gamma(U)`, and
`K_10 = K_01 = 0`. They satisfy `sum_r K_r* K_r = I` on the whole represented
parity sector. For any input rho in that sector and any untouched external
reference R, the retained branch is
`(K_r tensor I_R) rho (K_r* tensor I_R)`. Conjugation by V gives the same
instrument in the accepted raw owner code. This operator identity supplies
the arbitrary-reference statement; it does not assert that a covariance
matrix specifies every conditional many-body or reference state.

For the executed Slater preparations, the occupied fuel probability is
`(U C_total U*)_ff`, where C_total includes the stated fuel and reference
occupations and every other apparatus mode is empty. The four pointer
probabilities are therefore `[1-p, 0, 0, p]` in the stated order. The replay
stores this complete pointer law. It does not publish a full conditional
many-body density matrix. Both pointers are supplied local neutral bits at
the fuel center; their controlled writes are whole-star primitives. This
terminal comparison does not implement transport of a record to a second
spatially separated observer. Preparation reset, phase settings and pointer
duration remain supplied services with the boundaries stated in the ledger.

## Finite resources and work

The architecture supplies 432 encoded bits, twelve additional neutral
routing scratch bits, and two additional neutral pointer bits at the active
fuel center: 446 coherent bits. Raw owner dimensions are those of the star
isometries, rather than 446 independent physical qubits. The finite
classical archive retains 259 loop signs, 3,108 diagnostic signs, the gate
and transport records and the terminal two-bit word for each instrument
setting. These are retained record counts, not a minimum-memory theorem.
The twelve scratch bits return to zero by recorded outcome-controlled
operations. The archive is not reset or reused in this construction.
Initialization of the raw owners and unused registers is a supplied vacuum
preparation. Only the scratch diagnostics and the stated encoded loop
instruments are executed here; a physical owner-reset diagnostic is not
inferred from an energy penalty.

For one-body source h, define the observable H = dGamma(h), with zero bare
fuel and parity-reference energy. A converter pulse switches V on,
evolves for its declared duration with H + V, and switches V off. Conservation
of H + V during that interval implies

```text
Delta <H> = <V>before - <V>after = W_on + W_off.
```

The covariance evaluation `Tr(h C)` gives these expectations on the
declared Slater state. The ledger evaluates both sides independently for
each pulse. Free evolution preserves H. The fuel phase and pointer writes
commute with H because h has source-only support; their source work is zero.
This does not make their controller work zero or assign them a duration.
Every finite history retains these stage values and the total balance
residual. Numerical residuals have explicit gates; they are not certified
interval bounds. The pulse work identity is an externally driven source
balance, not autonomous energy conservation. Preparation-ramp and controller
work are outside this ledger.

## Independent contrast-memory source map

The construction manual also describes a contrast-memory law, implemented as the primary
model in `model/contrast_memory_spherical_v1`. It acts on qutrit loads and seam memories,
not on the scalar fixed-parity space. The independent implementation
`model/contrast_replay.py` checks a complete 99-state seam, retaining all ten
nonblank memory values and the blank value, each with nine load values.

Its one-cell total-two sector has 78 states: 66 pairs of distinct occupied
ports and twelve doubled-load states. These two parts are reducing under
the swap law; no outcome is dropped because its probability is zero on
the selected initial adjacent pair. The complete finite spectral
decomposition includes both parts. Roundoff normalization retains every
unnormalized squared amplitude and requires its sum to differ from one
by less than `1e-12` before normalization.

The explicit port map to the primary carrier is
`(0,1,5,7,11,9,10,8,4,6,2,3)`. It maps all thirty seams. If T is its induced
78-state permutation, the checked source identity is

```text
H_import = T* (H_primary + 12 kappa I) T.
```

The reference shift accounts for the primary sector's energy convention;
it changes only a global evolution phase. The full 99-state seam operator
also agrees before restriction to blank memory. Both identities are checked
at kappa 0 and 8. Each operator is affine in kappa, so the two values prove
coefficient equality for the specified finite matrices. The replay's
floating residuals audit that exact construction identity.

The raw-memory chart also tracks seam orientation. If the port map reverses
an ordered seam, it swaps the two load coordinates and maps every nonblank
record `(c,d)` to `(1-c,-d)`; blank memory stays blank. Denote this 99-state
permutation by S. The write obeys `S_record W = W S_load`, since exchanging
the two output loads exchanges the coin label and exchanging the input
loads reverses the contrast. The record generator commutes with the record
part of S. These identities imply `S h_kappa S* = h_kappa` on the full seam
space, including dirty memory. Applying this chart to each reversed seam
makes the port relabelling a consistent full-memory map. For example,
`(4,6)` maps to the oppositely ordered `(11,10)`. Omitting its memory-label
map gives operator error 9 at kappa 8; the complete chart has only numerical
roundoff residual. The executable retains both coordinate permutations and
the orientation of all thirty seams.

For complete-cell federations with external port partial matching and an
initial adjacent total-two state, the primary binding theorem applies with
`kappa > 3`. It bounds loss of pair support by `40/(2 kappa - 6)^2` for all
times. Its block lower bound separately requires `kappa >= 0`; a negative
coupling admits the uniform doubled-load counterexample of energy
`24 kappa < 12 kappa`. The finite time grid is a numerical screen of this
contrast law. Neither the source map nor that bound identifies the scalar
law with it or transfers the scalar instruments and clocks.
