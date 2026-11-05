# Native coherent detection of the cubic charged response

Let \(K=H-N/2\) be the specified open eight-site Weyl/contact generator at
\(c=1/8,U=1\), expressed in its exact hidden-spin site frame. Let \(g\)
be its unique half-filled ground. The probe is
\[
C_0=\sum_{y\in\{1,2,4\}}c^\dagger_{y\uparrow}c_{y\downarrow}c_{0\downarrow}.
\]
The source preparation is the retained literal Gray-tree preparation of
`centered_native_preparation`, stopped before its inverse site-frame word.
Its full 65,536-amplitude `hidden_success` column and every tree abort are
authenticated. The inherited 262,137 pulses are not reexecuted by this
package. The tree certificate and target-to-ground error are recomputed.
No autonomous ground-selection or protection against ongoing bare repair
is assumed.

## Exact coherent instrument

In occupation-qubit order the adjacent spin bilinears cancel their
Jordan-Wigner strings and
\[
C_0=\frac18\sum_{y\in\{1,2,4\}}
 Z_0(X_1+iY_1)(X_{2y}-iY_{2y})(X_{2y+1}+iY_{2y+1}).
\]
There are 24 signed unitary Pauli words. Five work bits prepare their
uniform coherent superposition, with zero amplitude on the eight labels
whose neighbor field is 3. Three bits use Hadamards. The two neighbor bits
prepare \((|00\rangle+|01\rangle+|10\rangle)/\sqrt3\).
The variable rotation uses the literal hexadecimal angle in `circuits.py`;
its sine and cosine are enclosed independently against \(1/\sqrt3\) and
\(\sqrt{2/3}\). The complete preparation unitary, including its action on
every other input, is inverted after selection.
For two real Y rotations, their difference satisfies
\((U-V)^\dagger(U-V)=4\sin^2((\theta-\phi)/2)I\). Thus the enclosed
first-column error equals the full operator error and also bounds the
inverse rotation on entangled inputs. The factor two in the certificate
accounts for both occurrences.

For each label, `selection_template` implements its signed Pauli word.
Controlled phases use
\[
\exp(i\theta\prod_j |b_j\rangle\langle b_j|)
=\prod_{S}\exp\left(i\theta 2^{-k}(-1)^{\sum_{j\in S}b_j}
                            \prod_{j\in S}Z_j\right).
\]
The scalar term is retained. Hadamard conjugation implements conditional
flips. Each resulting string is compiled into the declared one-bit and
two-bit native pulses. An independent integer occupation calculation
checks all 65,536 entries of the sum of the 24 signed permutations against
\(8C_0\). Thus the zero selector outcome has exactly the ideal Kraus
operator \(C_0/3\). Unitarity supplies the complete complement without a
success normalization or an imposed two-outcome Lüders completion.

The 1,008-pulse word has 56 bounded-support contraction blocks. Exact
template consumption accounts for every pulse. Rational trigonometric
enclosures bound literal angle errors, and the source's quadratic-field
port enclosures bound its full six-level generators. Telescoping the
complete local instruments, with orthogonal first-failure records, gives
the rational operator error in `control_certificate`. It is below
\(1.6\times10^{-13}\), for every field input and untouched quantum
reference. Adding the authenticated preparation and target errors gives
the stated complete-output norm bound and total-variation bound for any
subsequent readout. Numerical replay admissions are reported separately
and do not define this theorem's constants.

## Complete field, reference and record domain

The probe acts on original field modes \((0,1,2,3,4,5,8,9)\). Their eight
occupation bits form local rows. The other eight field bits remain an
untouched quantum reference, with 256 coherent columns. This is an
occupation-coordinate permutation; all fermionic signs occur in the
probe. It is not a signed tensor-factor identification of odd CAR
operators on disjoint regions.

Linearity permits partitioning the reference basis into groups of at
most 32 columns. The complete 2,097,152 successful joint amplitudes are
reassembled, so off-diagonal reference coherences are retained. Summing
column norms is only the partial trace used to assign record
probabilities. Every first-export failure has its original source Kraus
map and weight. `failure_reference_amplitudes` reconstructs any such
record's full coherent reference columns; rectangular padding is excluded
from the outcome set.

The final native readout has six raw outputs per bit. The archive retains
all \(6^5=7776\) raw selector records in lexicographic order. Exactly
7,744 are zero because successful exports lie in the two-dimensional code
space. All 32 valid binary records are included, alongside the separate
native export failures and inherited preparation aborts.

## Spectral attachment and boundary

`charged_composite_search` independently certifies the complete 138-line
spectral measure of this same \(C_0g\). Its total weight divided by nine
encloses the ideal zero-selector probability. The native execution gives
approximately 0.08261623248916, within the complete preparation/control
error of that certified interval. All other detector outcomes are
retained. This apparatus detects the total cubic response; it does not
resolve its 138 energies. No energy clock is executed here.

The hidden spin label is a symmetry of this finite source. The result
does not identify observed particle species, relativistic spin or rest
masses, a volume-limit pole, an interacting continuum field theory or
Einstein dynamics. The declared finite simulator and observer controls
are readable source inputs, not inferred laws of nature.
