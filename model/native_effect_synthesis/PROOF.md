# Bounded native rank-one effect synthesis

The input is the existing weighted H66 carrier, embedded by its declared
nine-bit code into all 512 working addresses. Coherent requests use only the
X, Y, Z, I, H and J pulses provided by
`observer/source_scalar/compiler.py`, with its existing proper-code load and
complete export services. The final detector is the eleven-repair E-record
tree in `native_effect_coverage`. Coherent gluing, processor blanks and Kraus
trace semantics remain declared. This is a finite observer composition of
those services; the bare repair/hold grammar does not issue the coherent word.

## An unequal-weight projector from existing controls

Fix two original labels a<b. The existing reversible X/CNOT routing takes
their work codes to addresses 2 and 1. One CNOT takes these to 2 and 3. On the
three low bits let

\[
 A=Y_0\frac{I-Z_1}{2}\frac{I+Z_2}{2}.
\]

The four Pauli terms commute. The existing `pauli_rotation` words therefore
implement the four factors of exp(i theta A). On the ordered pair 2,3 its
matrix is [[cos theta,sin theta],[-sin theta,cos theta]], and it is the
identity on the other six low-bit addresses. Undoing the CNOT and routing
gives this rotation on each routed spectator pair. The selected source pair
is reducing; other source and padding inputs are retained, not projected away.

Choose theta with cos theta=3/5 and sin theta=4/5. Selecting final native
label a then has effect

\[
 P_{ab}=\frac1{25}
 \begin{pmatrix}9&12\\12&16\end{pmatrix}_{ab}.                 \tag{1}
\]

The emitted theta is a certified binary approximation described below.
The pair (0,65) uses 80 pulses, with a complete error bound below 2/10^9.
Equation (1) lies outside the single balanced-analyzer positive menu;
that obstruction did not exclude other sequences of admitted source gates.
This construction supplies such a sequence. It does not change the native
repair rule or derive a measurement update from an effect alone.

## A rational rank-one target on the entire source

Let v be a nonzero 66-component Gaussian-rational vector and let
0<=lambda<=1 be rational. Each supplied rational numerator and denominator
has at most 64 bits. The target is

\[
 E=\lambda\frac{|v\rangle\langle v|}{\langle v,v\rangle}.
                                                                  \tag{2}
\]

Choose one existing padding address p. In the same 512-dimensional register
define the unit vector

\[
 |w\rangle=\sqrt{\lambda}\frac{C|v\rangle}{\|v\|}
             +\sqrt{1-\lambda}|p\rangle,                         \tag{3}
\]

where C is the original source embedding. No extra qubit or additional
ancilla service is used. All 512 squared magnitudes in (3) are rational.
The construction emits a binary amplitude tree followed by a diagonal phase
word. Its ideal target W sends the literal work address 0 to w. The selected
measurement word is W inverse followed by existing X routing from address 0
to one genuine source code c(a). This final routing is essential: address 0
can be padding, and a padding-abort coarse record cannot isolate it.

After the entire coherent word, perform the complete source-address export.
Only its source-code success sector enters the native tree detector. The
446 padding Kraus rows contribute to the owning export's single observed
abort flag; this does not add address-resolved padding measurements.
No chart projection occurs between coherent pulses. The selected native outcome has ideal effect

\[
 C^*W|0\rangle\langle0|W^*C=E.                                 \tag{4}
\]

All other source outcomes, processor failures and the padding-abort flag remain
present. In particular (4) neither conditions on source-code success nor
specifies a Lueders update. The actual inverse word consists of reversed
negative native requests, executed independently; it is not an adjoint of a
returned forward floating array.

## Certified binary angles without inverse-trigonometric assumptions

Set Q=2^41. A split with exact left and right squared norms L,R needs an
angle theta in [0,pi/2] satisfying

\[
 \cos(2\theta)=(L-R)/(L+R).
\]

For L=R=0 choose zero. For other splits, `split_angle` bisects integer grid
points theta=n/Q and compares a directed enclosure of cos(2theta) with the
exact rational right-hand side. If adjacent grid points bracket the answer,
their separation is 1/Q. If a cosine enclosure contains the target, its
width delta satisfies 4 delta Q^2<=1. For a,b in [0,pi],

\[
 |\cos a-\cos b|\ge2(|a-b|/\pi)^2.
\]

Indeed use the cosine-difference identity, symmetry about pi/2 and
sin x>=2x/pi on [0,pi/2]. Hence the unresolved split differs from its target
by at most 2 sqrt(delta)<=1/Q. Zero endpoints are handled explicitly.
The initial upper endpoint is the ceiling of pi Q/2, so its possible
overshoot is less than 1/Q and never enters an interior comparison.

Sine and cosine enclosures use directed 120-bit integer interval arithmetic,
the Taylor terms through degrees 97 and 96, and exact integer checks that the
remainders are less than one interval unit for |angle|<=4. Every division
rounds outward. The pi enclosure follows
pi=16 arctan(1/5)-4 arctan(1/239), with 90 alternating-series terms and their
next-term remainders. The identity follows from the tangent addition formula
and the indicated first-quadrant range. No `atan2`, library square root or
floating trigonometric value decides the synthesis certificate.

For a Gaussian rational x+iy, the same split applied to x^2,y^2 determines
the first-quadrant phase; sign and quadrant corrections use the certified
pi grid point. Its angular error is at most 2/Q. Zero coordinates have phase 0.

At each tree level the desired signed angles have the Walsh expansion

\[
 \widehat\theta_m=2^{-c}\sum_{u=0}^{2^c-1}
             (-1)^{u\cdot m}\theta_u .                          \tag{5}
\]

Each term is an existing tensor-Pauli request. The commuting terms exactly
recover every conditional split. There are at most 256 split prefixes at a
level. The full diagonal phase has 512 Walsh coefficients. With 41-bit grid
angles and |phase|<4, every signed-sum numerator has absolute value below 2^52.
Every coefficient is consequently exactly binary64 representable. The code
checks exact Fraction equality after conversion and verifies every inverse
Walsh relation. It retains zero coefficients and their actual service requests.

Each of nine conditional-rotation layers differs from its exact target by at
most 1/Q in operator norm; the diagonal phase differs by at most 2/Q. Telescoping
the unitaries gives

\[
 \|W_{\rm grid}-W\|\le11/Q.                                    \tag{6}
\]

This is a full 512-dimensional operator estimate, not only an estimate on
the prepared column. It therefore applies to the inverse measurement and
arbitrary untouched references. The short real-pair construction uses 1/Q.

## Existing finite resource domains

The amplitude tree has 12291 native pulses. The full nine-bit diagonal word
has 11270 pulses, including its scalar phase. Mapping address 0 to a source
code uses at most 18 more pulses. Thus every rank-one word has at most 23579
pulses, below the unchanged 65536-pulse source-word admission ceiling. Its
description also satisfies the existing payload limit. All registers stay
at nine bits. The native detector has eleven repair attempts, within the
unchanged 64-step host bound.

These are distinct execution interfaces. Generic `execute_word` returns
`RESOURCE_STOP` before executing the nine-bit word because its host accepts
at most six qubits and 256 pulses. The existing direct source compiler executes
the declared nine-bit pulse tape on source/reference amplitude columns. The
generic expanded tree evaluator also retains its three-step resource
stop. The complete native tree instrument is supplied by its exact lazy branch
factorization. Neither resource stop is relabeled as a successful full run.
The bounded angle adapter emits existing public source pulses; it does not
claim that the libm-based state compiler emits identical angles.

## Full primitive and complete-output bounds

For every distinct emitted pulse, including negative and zero requests, the
certificate compares the original stored 6-by-6 array with the ideal full
source pulse at its exact binary64 angle. Every returned real and imaginary
component is converted to an exact fraction. Directed trigonometric intervals
bound the squared Frobenius error, which must not exceed (10^-11)^2.
For H, a directed integer-square-root enclosure also bounds 1/sqrt(2).
No polar projection, normalization or deletion of complement blocks occurs.

The ideal synthesis and resource bounds hold throughout the stated rational
input domain. The numerical implementation bound below is conditional on
`certify_word` passing these exact checks for every emitted primitive. A
failed primitive rejects that word's certificate. This does not assume a
uniform rounding theorem for SciPy or assert that every possible input passes
the stored-array gates.

The Pauli basis and parity constructions also contain fixed `math.pi` angles.
Their binary payloads are separately compared with the rational pi enclosure;
the error is at most 10^-15 before division by 2 or 4. The relevant ideal native
generators have norm 1. Thus g=2*10^-11 bounds each stored primitive against
the required ideal gate, including the fixed-angle error. This explicitly
accounts for the CNOT and basis changes needed by (5). It bounds stored arrays;
physical hardware error is outside this certificate.

For N<=23579, put x=Ng<1. The good-code product has gain at most 1/(1-x),
and its difference from the rational-target unitary is at most

\[
 e=\frac{x}{1-x}+\frac{11}{Q}.                                 \tag{7}
\]

Use 1/Q for the short pair. Each pulse's complete processor complement has
operator norm at most g. Keeping the first failure, with its addressed code
reset and every spectator/reference retained, gives total failure mass at
most F_N=N g^2/(1-x)^2. Subsequent failure branches stop the word. Every
complement row, including an exactly zero one, remains a public address.

After the good-code product, the complete address export and native tree
instrument are exact completely positive trace-preserving maps under their
declared semantics. Their composition contracts the trace norm. For every
source/reference input the complete output comparison therefore has

\[
 \delta_N\le\frac{e(2+e)+F_N}{2},\qquad
 G_N\le\frac1{(1-x)^2}+F_N .                                   \tag{8}
\]

The first quantity is a half-diamond bound and the second a gain bound.
The general construction has delta_N<1/10^6. Rounded maps need not be exactly
trace preserving, and returned positive record blocks need not be normalized.
No survivor normalization is used. An accept-score expectation can change by
at most 2 delta_N. Binary acceptance/rejection is only a coarse-graining of
the complete record channel. Complete first-failure queries preserve the same
arbitrary reference; the numerical selected-row and output-array diagnostics
are independent checks, not extra zero-error statements about floating
state propagation.

## What remains for the operational Born theorem

This gives a bounded-precision source route for all targets (2) in the stated
finite rational-description domain, and an explicit answer to the single-menu
unequal-weight-projector obstruction. It neither compiles an arbitrary
positive effect nor constructs a general Naimark isometry. A whole eigenbasis
can require a different word budget; no such budget is silently enlarged.
Individual rank-one tests also do not implement a common disjoint-event
refinement for every coexistent pair.

`OperationalBornBridge.lean` requires implementation-independent event weights
and universal coexistent refinement. Its prose identifies the further need
for a dense-effect interface when exact uncountable realization is unavailable.
The fixed binary64 alphabet, fixed pulse limit and the nonzero allowances
here do not supply an arbitrarily accurate limit. The source consumes
Kraus trace probability semantics. This result adds an operational synthesis
construction; it does not derive that semantics from repair counting,
discharge universal Born premises or identify relativistic fields or particles.
