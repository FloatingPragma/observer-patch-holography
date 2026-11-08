# Complete finite instruments from admitted source controls

## 1. Source services and the observable being implemented

Use the two-flavor spin source, selected finite ground and control interface
of `native_spin_observer` and `native_ordered_ground`. The source's exact
twelve six-level port generators synthesize single-register X, Y, Z and
scalar phases, controlled phase, Hadamard and the existing three-pulse CNOT.
Coherent addressed load/export, blank pointer registers and classical
pointer commitment are the same declared apparatus services. No additional
projective spin measurement is assumed.

One logical spin is encoded as `0 -> down`, `1 -> up`, using five source
pulses. Its inverse decodes onto the up bit with the down bit blank. A
complete ideal spin operation can therefore be executed as decode, a word
on the logical bits and observer ancillas, then re-encode. This costs ten
wrapper pulses per spin. Native errors and their proper-code complement
outputs are retained, rather than projected back into the spin code.
Repairs occur between completed operations, as in the inherited control
protocol. A native pulse count is not a calibration of laboratory time or
the duration of an autonomous spin evolution during measurement.

The literal occupation repair alone is the identity on the singly occupied
spin sector and has no spin-dependent records. The construction uses the
explicit coherent service above; it does not contradict that obstruction.

## 2. A terminating complex rank-one compiler

Let `v` be a nonzero complex rational vector of length `D=2^n`. Clear a
common denominator and write its coordinates as pairs of integers
`a_j+i b_j`. The squared coordinate weights `a_j^2+b_j^2` are exact.
The existing rational square-root/arctangent algorithm computes each
positive amplitude-tree angle to error at most `2^-b`.
Zero-weight nodes have angle zero. All branch decisions use integers.

The phase of a nonzero coordinate is computed from the same algorithm
applied to `a_j^2,b_j^2`, followed by its exact signed quadrant. A rational
enclosure of pi handles the negative-real quadrants. Zero coordinates
receive phase zero. The phase error is at most `2*2^-b`; no undecidable
real-zero test is used. A Walsh transform expresses the phase diagonal as
commuting Z products. Gray parity circuits implement those factors with
the original CNOT, retaining the final CNOT in each cycle.

The amplitude tree uses at most `4D` pulses and its phase diagonal at most
`4D`. Their ideal unitary `A_v` satisfies `A_v|0>=v/||v||`.
The emitted rational-angle approximation obeys

    ||A_tilde-A_v|| <= (n+2) 2^-b.                           (1)

The tree estimate sums the maximum error of each level's mutually
exclusive conditional rotations. The diagonal estimate is its maximum
coordinate phase error. This is an operator estimate on all inputs.

The marked phase `exp(i theta |0><0|)` has the exact expansion

    |0><0| = 2^-n sum_(S subset {1,...,n}) Z_S.               (2)

Its commuting parity word uses at most `4D` pulses. Applying the words in
the order `A_tilde inverse`, marked phase, `A_tilde` implements the rank-one
phase on `v`, to operator error at most twice (1), with at most `20D`
pulses. The special phase pi is the Householder reflection. All global
phases are retained because a later word may be controlled.

`generic.py` implements this compiler for every complete finite complex
integer vector supplied to it. It uses rational angles and terminating
integer interval algorithms; it does not turn a fixed binary64 command
interface into an arbitrary-precision one.

## 3. Every computable finite instrument has an explicit source word

Consider a specified complete finite instrument on `q` spins, with
computable Kraus matrices `K_r` satisfying `sum_r K_r* K_r=I`. Multiple
Kraus operators per classical outcome are assigned separate retained
environment labels. Here m counts all fine Kraus/environment branches,
including different Kraus operators with the same coarse classical label.
For example, a one-outcome depolarizing qubit channel has four such branches.
Pad that fine alphabet to `M=2^ceil(log2 m)` and
include its zero maps. The complete isometry is

    W psi = sum_r K_r psi tensor |r>.                       (3)

Its action on an arbitrary untouched reference is `W tensor I_R`.
Extend (3) from the blank-record subspace to a unitary on `q+log2 M`
bits. Finite Gram orthogonalization constructs such a completion.
For computable entries, a nonzero pivot is selected by searching for a
strict positive rational lower bound; the known rank guarantees a pivot.
An equality oracle for arbitrary computable reals is unnecessary.

A unitary of dimension `D=2^(q+ceil(log2 m))` has a product decomposition
into at most D complex rank-one phases and a final diagonal phase. One
stable elimination step maps a unit column x to `s e_j`, where `s=+1` or
`-1` is chosen with `Re(s x_j)<=1/2`. Overlapping rational tests choose
such a sign without deciding equality at a threshold. With
`v=x-s e_j`, the unitary

    H=I-v v*/(1-s x_j)                                      (4)

sends x to `s e_j`, fixes earlier coordinate axes, and has only one
nontrivial eigenvalue `-conjugate(d)/d`, where `d=1-s x_j`.
The denominator's real part is at
least one half. Thus (4) is a rank-one phase with computable parameters.
Approximate its direction by complex rational coordinates and its phase by
rationals, with independently certified errors. For the eigenvalue lambda,
enumerate rational theta until a refined rational sine/cosine enclosure
certifies `|exp(i theta)-lambda|<epsilon_phase`. Density and computability
ensure termination without a global argument branch at a phase cut. The
phase-error budget is this chord/operator distance. The factor compiler in
section 2 then applies. A positive finite POVM is included by taking the
computable positive square roots `K_r=sqrt(E_r)`; this fixes its finite
Lüders instrument rather than silently choosing a postmeasurement state.

The explicit pulse bound, including spin decoding and re-encoding, is

    N <= 20D^2+4D+10q.                                     (5)

The bound is exponential and is not an efficiency claim. The existing
integrated-port command algorithm gives each ideal pulse complete-isometry
error at most `29*2^-p`. Telescoping on the full source, ancillas and every
first-export-failure register gives at most `29N*2^-p`. A complete-isometry
operator error epsilon bounds half the diamond distance of the committed
instruments by epsilon, including arbitrary reference inputs.

For a requested half-diamond error epsilon, `bounds.py` assigns epsilon/8
to normalized rank-one-vector rounding, epsilon/8 to factor phases,
epsilon/4 to the state-preparation trees and epsilon/2 to source ports.
The factor direction error is at most `epsilon/(32D)`, since changing
a normalized direction by delta changes its rank-one phase by at most
`4 delta`. Each phase error is at most `epsilon/[8(D+1)]`. Choose b so
`2D(n+2)2^-b<=epsilon/4` and p so `29N2^-p<=epsilon/2`.
These are terminating integer choices.

This theorem provides a finite program at each finite register size and
accuracy under the inherited coherent-service contract. The expanded
host evaluator's register and resource ceilings are separate. It supplies
no efficient infinite-ground preparation, no automatic access from a
bare-repair observer and no physical detector locality theorem.
A coarse instrument is recovered by grouping its fine environment labels
in the stated observation map. The complete implementation retains that
environment; its size is never inferred merely from the number of coarse
outcomes.

## 4. The selected cube's coordinate and its complete backaction

Use the unchanged eight-site source at `J=1,B=1/4` and its certified
256-component ground witness. In the Marshall frame, project the local Z
probe onto its eight graph-reflection characters. The integer columns
`c_k` are exactly orthogonal and nonzero. Set

    psi_k=c_k/||c_k||,
    eta_x=8^-1/2 sum_k (-1)^(k.x) psi_k,
    Pi_x=|eta_x><eta_x|,    Pi_dark=I-sum_x Pi_x.             (6)

The companion `spin_covariant_localization` bounds the difference between
this proposal-defined detector and the true source-ground detector using
the full ground certificate. These labels are the finite reflection group;
they are not small continuum momenta or identified particle positions.
If each normalized character vector has error at most d, the eight-column
Frobenius error is at most `sqrt(8)d`. For the Lüders instrument used here,
the rank-one projection rows and dark projection give a dilation error
at most `sqrt(8)*sqrt(8)d=8d<9d`. Thus the retained numerical allowance
`9d` also bounds this instrument, independently of the companion packet's
destructive coordinate-readout convention.

Eight retained integer Householder normals define an exact rational
orthogonal matrix O. The verifier applies O to every full integer `c_k`
and checks each norm exactly. Integer square-root enclosures bound

    delta^2 >= sum_k ||O psi_k-e_k||^2.                     (7)

The proposal generator is not trusted. Changing a normal requires the
complete exact residual check again. The resulting delta is below
`10^-10`. The complete ideal detector built from O differs from (6)
in Stinespring operator norm by at most `3 delta`: the rank-one rows
contribute at most `4 delta^2` in square sum and the complementary
projection at most `4 delta^2`.

Apply O, then the three-bit Walsh transform on the low bits. Compute a
validity flag when the high five decoded field bits are zero. Copy the
three low bits to the coordinate register only when that flag is one.
This labels anchors 0 through 7 as records 8 through 15. All other 248
states receive the single dark record 0. Records 1 through 7 and the
pointer's reserved nonbinary slots are retained as zero ideal outcomes.
Uncompute the Walsh transform and O on the field before committing the
pointer. Therefore each ideal field branch is exactly its projection,
including the **whole dark projection**. The dark subspace is not secretly
measured in a finer basis. This proves the probabilities and postmeasurement
state for every coherent input and reference.

Decode and re-encode the physical spin code around this word, and apply
the actual Marshall phase on both sides. The resulting instrument acts
on the same physical spin variables used by `native_ordered_ground`.
No positive Marshall-frame vector is substituted for its physically
prepared counterpart. There are 49,228 decoded-core pulses and 49,324
physical pulses including these wrappers.

`circuits.py` independently checks each signed amplitude-tree angle,
every fixed gate angle and every literal binary64 integrated-port
coefficient. It adds those errors to `3 delta` and the true-ground detector
error. Their sum is below `10^-6` in half diamond norm, without accepting
only a successful source branch. Consequently the exact finite group's
position covariance has complete-instrument defect at most twice this
error when compared using the ideal group action. An implementation of
the group action must add its own word error.

## 5. Executed domain and continuum boundary

The literal decoded core is executed on all 4096 quantum output rows,
with two untouched reference columns coherently pairing a probed-ground
vector and a dark vector. It retains all sixteen binary pointer labels
and every first-export failure. The full twenty-bit encoded wrapper has
the all-input analytic bound above; this packet does not mark that larger
literal execution as performed. Its code/Marshall attachment is proved
by the exact encoding identities and charged in the error bound.

The finite detector therefore has a concrete source implementation under
the existing service interface. It can be applied to the actual selected
finite state or to any other source/reference input. Preparing an input
approximately adds its complete-state error; the instrument bound itself
does not depend on a state preparation.

Section 3 also supplies the finite observer operations needed by a proposed
continuum approximation. Convergence to a particular continuum position
measure must use the same source states, spatial translations and clock.
This finite implementation theorem does not prove that convergence or
second-moment tail bounds. The selected cube's cyclic coordinate code
leaks under its Hamiltonian, as its existing energy-variance certificate
shows. Its specified input is a local-spin response; attachment to a
protected topological charge is unproved. No particle stability, dispersion shell, rest mass, mass ratio or
experimental particle label follows from the measurement compilation.
