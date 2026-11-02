# Native preparation and charged energy readout of the centered cube ground

This package prepares the actual centered `N=8` cube ground to a proved
error using the existing six-level native source, then executes the
complete four-bin charged energy instrument. The reference is the unique
hidden-spin-zero, eta-zero ground of `K=H-N/2`, at `c=1/8,U=1`. It is not
the original `N=4` spin-two reference. The complete field carrier has
65536 rows; every preparation and clock outcome, including first-export
failures, remains represented. A four-bin experiment does not separately
resolve the 228 visible spectral lines and does not identify physical
particle species or relativistic masses.

## A source-derived target, without an amplitude oracle

The input is the fixed integer witness of
[the centered spectral certificate](../centered_particle_spectrum/README.md).
Its exact symmetry basis and dyadic ground column expand into 4484 nonzero
signed integer amplitudes on the full occupation carrier. There are 4900
balanced configurations, including the exact zeros. `target.py` rebuilds
this expansion, reduces a common integer divisor, and fixes a harmless
global sign. The normalized target is `a/sqrt(sum a^2)`. This normalization
specifies the input program; no output of a native experiment is normalized
or postselected.

The exact residual/Gram certificate is rerun on the complete 32-dimensional
ground symmetry block. Its isolated lowest eigenline, together with the
proved global ground symmetry and uniqueness, bounds the normalized
target's distance from the exact centered ground by approximately
`1.245e-13`. It also supplies an enclosing interval for the ground energy.
The integer amplitudes and this residual bound are part of the receipt.
They are computed from the readable source and a checked witness, not
obtained by querying unknown quantum amplitudes.

## A shorter explicit native amplitude tree

A binary tree on 16 work bits splits each prepared higher-bit prefix into
the next bit's two child weights. At all but the last level use the square
roots of the two exact integer subtree weights; at the last level use the
signed leaf amplitudes. The ideal angle is `theta=-atan2(y,x)`, so
`exp(i theta Y)|0>=(cos theta,-sin theta)` has the desired components.
An empty subtree uses angle zero. All those zero requests are retained.

At a level with k higher control bits, put `n=2^k` and let
`g(j)=j xor (j>>1)` be the cyclic binary reflected Gray code. Compute

```
beta_g = (1/n) sum_a (-1)^(a.g) theta_a.
```

The chronological circuit is `Y(beta_g(j))` on the target, then a CNOT
from the one control bit changed between `g(j)` and `g(j+1 mod n)`.
There is no CNOT for k=0. In control sector a, the accumulated X toggles
conjugate each Y rotation by `(-1)^(a.g(j))`; the complete cycle restores
X to identity. Thus its exact action is `Y(theta_a)`. Keeping the final
cycle CNOT is essential on arbitrary control inputs. A mutation test
removes it and detects the changed full operator.

Each CNOT is the existing phase-exact three-pulse `H,J,H` source word.
Each Y request is one native pulse. Consequently the full preparation has

```
65535 Y pulses + 65534 CNOTs = 262137 native pulses.
```

No 16-qubit unitary is inserted as a primitive. `tree.py` uses a private
192-bit arithmetic context only to propose raw angles. An integer fast
Walsh transform of their exact binary64 representations proposes the
literal beta angles, avoiding a quadratic transform. Independent
certification checks every raw conditional two-vector with exact
rational trigonometric enclosures and inverts the rounded beta transform
exactly. For each level it adds the maximum raw two-vector error and the
maximum reconstructed-angle discrepancy. These are maxima over orthogonal
control sectors, not an erroneous sum over all possible prefixes.
Telescoping the 16 levels gives the signed-state error.

The tree prepares the hidden-frame state. The 184-pulse inverse onsite
frame word is then executed on the full local Fock spaces. Its exact
Gaussian-integer certificate checks vacuum, single-occupation and double-
occupation phases. This yields the state in the original Weyl-coordinate
source carrier used by the subsequent Hamiltonian and local probe.

## Full native generator and outcome accounting

Every pulse is generated from the same twelve response ports on the full
six-level processor. Logical one- and two-bit gates use the existing
2- and 4-dimensional proper codes and coherent load/export service.
All complementary export Kraus outcomes are retained, including zero
outcomes. A first failure stops that history with its event and complement
label; success continues without renormalization.

For efficient certification of many literal angles, let `b_p` be the
actual binary64 base controls at angle one. For each literal angle theta,
read the actual controls `c_p(theta)` and compute the exact rational
residual `delta_p=c_p(theta)-theta*b_p`. Then

```
||sum_p c_p(theta) D_p-theta T||
 <= |theta| ||sum_p b_p D_p-T|| + sum_p |delta_p| ||D_p||.
```

The base error and each port-generator norm are bounded from the exact
`Q(sqrt(5))` source. This identity uses the actual rounded coefficients;
it assumes no floating-point relative-error model. Both generators are
anti-Hermitian, so Duhamel bounds their unitary difference without an
exponential norm penalty. The full export isometry preserves that bound.
Fixed CNOT/frame angle errors, target-state error, and all source control
errors are included separately.

Execution partitions the literal pulse tape into consecutive blocks on
at most three work bits. Each block is contracted from its actual full
native matrices and every first-failure Kraus matrix. A compiled C++
tensor contraction applies that complete local instrument over every
65536-row state, without number, spin, symmetry or spectral projection.
The block is then discarded to save host storage; the immutable program
and event offsets reconstruct every failure vector from its unnormalized
prefix. Failure probabilities are stored for every pulse and every slot;
nonphysical padding slots are identified separately. This is actual
native pulse execution, not replacement by the ideal amplitude tree.
Numerical residuals are diagnostics; rational generator bounds establish
the control claim.

## Centered charge and four-bin clock

At the origin, the phase-exact native charge instrument first copies the
input down-mode occupation to a supplied blank tag, then applies
`gamma=Z_0 X_1`. Its two logical outputs are exactly `c_1^dagger G` and
`c_1 G`, with the CAR sign retained. Reading the tag produces two classical
columns, each of ideal norm squared 1/2. All six raw tag records and any
prior export failures are retained. The tag is reused as a blank clock
bit only after that classical readout; a second clock bit is supplied
blank. This is not a coherent merger of the two charge histories.

The two clock bits are prepared uniformly. Conditional field durations
are one and two, using the existing complete full-Fock Strang circuits
with 16 and 32 steps. The actual source Hamiltonian remains H. On each
controlled duration t, append the explicit number phase

```
exp(+i t N P_clock/2),
exp(+i t n_j P_clock/2)
 = exp(i t [II-ZI-IZ+ZZ]/8)       for each field mode j,
```

and the one-bit ground phase `exp(+i t E_G P_clock)`. Their product with
the controlled H evolution is the desired controlled
`exp(-i t(K-E_G))`. Every exact split term conserves N, so moving this
number phase after the split introduces no extra splitting error.
The actual rounding of these extra source words is included in the bound.
The ground-energy phase uses the certified interval, not an assumed exact
floating eigenvalue.

The number correction cannot be absorbed into one scalar reference for
both charges. Omitting it while retaining the K ground phase shifts the
addition/removal frequencies by 9/2 and 7/2. If the reference scalar is also
changed to the H eigenvalue `E_G+4`, the shifts are instead +1/2 and -1/2.
The countercontrol retains this distinction.

A positive Fourier transform then reads the complete M=4 Fejer energy
POVM. For each charge tag the ideal joint law is

```
p_tag(k)=sum_{i=1}^{228} w_i |sum_{j=0}^3 exp(ij(2*pi*k/4-E_i))|^2/16.
```

The finite spectral certificate encloses every positive weight w_i and
excitation E_i; they are recomputed before interval propagation. Both
charge laws agree by the centered source symmetry. The energies lie
between approximately 0.36168 and 3.88644 in declared source-clock units,
within one 2*pi phase period. These four broad outcomes are a finite-
resolution readout of that measure, not 228 individually resolved lines.
Removal has positive physical excitation energy; its negative sign is
only a signed response-frequency label.

The full clock has 607071 native pulses, and a complete successful history
including preparation and the charge instrument has 869405. All 216 raw
charge/two-clock pointer combinations are listed, including the 208
proper-code zeros, and every earlier export failure remains in the
absorbing stop record. Full execution retains 262144 clock/field rows and
two classical charge columns. The source error plus the proved full-Fock
Strang error is less than 0.03 in complete-state norm and hence in total
variation after the complete pointer measurement. It bounds any coarse
event, with no division by a success probability. Finite-shot sample
requirements are stated separately; no random laboratory shots are
claimed to have been collected.

## Scope and reproduction

The supplied apparatus includes blank registers, a coherent tensor memory,
source control precision, the pulse schedule, classical tag readout and
feedforward, and final pointers. This experiment has no interleaved bare
field repair during preparation or clock evolution. It therefore does not
prove that the centered ground is autonomously selected or protected by
that repair law. Its resource counts and complete outcomes make this
finite controlled experiment concrete under the admitted apparatus.

Run `python -W error -m pytest -q model/centered_native_preparation` for the
strict checks. `python -W error -m model.centered_native_preparation --write`
and `--check` execute a fresh complete preparation/readout in a cold worker,
with a lossless full-row archive and a bound current source manifest.
The new checks and inherited receipt custody are reported separately;
older source receipts are neither modified nor silently recertified.
