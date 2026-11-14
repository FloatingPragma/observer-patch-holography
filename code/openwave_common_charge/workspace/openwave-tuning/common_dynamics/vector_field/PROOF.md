# Autonomous fields on the retained carrier current

The declared source is a local fermion-mode walk with one shared modular
electric register on every actual carrier seam. It extends the proposed
`native-local12-fermion-modes-v1-proposed` construction. The field alphabet
is Z16. The source retains all 2048 matter modes, including 1152 physical
modes, 288 four-bit field registers and three clean arithmetic work bits.
Fermionic statistics, coherent code transfer and the field couplings are
architecture choices. The underlying integer-port balancing law is not
asserted to reproduce this many-mode evolution.

## One evolution and one Gauss constraint

Let B have +1 at the lower-index endpoint and −1 at the other endpoint.
The matter occupation at a port is the sum over its twelve physical event
modes. Onsite mode reflections preserve that sum. A transport pair swaps its
two occupations and increments the corresponding electric register by
`nu_a-nu_b`. Thus n+B F=G modulo16 is invariant under the full retained
matter-current step. The modulus distinguishes every local occupation from
0 through12. It does not by itself distinguish unrestricted regional totals.

Let A be the local face matrix supplied by `cell_complex.py`. Its 208
columns have BA=0. Define the field shift U_f by F→F+A_f modulo16 and the
centered electric value e(F) in −8,...,7. The field factors are

\[
 U_B=\prod_f\exp\{i\beta_f(U_f+U_f^*)/2\},\qquad
 U_E=\exp\{-i\alpha\sum_e e(F_e)^2/2\}.
\]

The total step is V=U_E U_B Gamma(W)_F. All face shifts commute. Every
factor preserves the same modular Gauss constraint and has a complete
unitary action on its declared registers. The coefficients alpha and beta_f
are recorded inputs. The electric observable and each plaquette generator
are kick generators; their sum is not asserted to be a conserved logarithm
of V. The branch-free observable I−(V+V*)/2 is conserved by V. Its conservation
does not identify it as the time generator or a physical rest energy.

The finite register obeys the Weyl clock/shift algebra. Its integer-valued
representative has the boundary correction
[E,U]=U−16|boundary_out><boundary_in|. Integer Gauss and the unbounded U(1)
commutator require a separate cutoff transfer. A normalized uniform state
on an infinite integer flux fibre is not used.

## Source words and complete records

An individual face column has a pivot of coefficient ±1. Modular shears
replace every other supported field coordinate by
F_e−A_e A_p F_p. They turn translation by A_f into one cyclic pivot shift.
A native four-bit Fourier transform, a diagonal cosine phase, and the
inverse Fourier transform implement its exponential. Reversing the shears
gives the full plaquette kick on all field states. The field circuit uses
the actual M6 scalar compiler, the native Toffoli/CNOT routines and complete
coherent code services. The arithmetic test executes all256 clean and2048
dirty inputs of a two-register modular adder. Each local Fourier, electric
and cosine word is executed on all16 columns.

`native_field.complete_word` emits the complete matter, current and field
word. Onsite Gaussian reflections and adjacent mode swaps retain every
fermionic spectator sign. An optional fresh native phase tick follows the
completed word; it uses wire3203 and angle−h/16, giving measured relative
phase h/8. This marks the shared declared time t=n tau0 h. It is not an
identification of elapsed processor time with pulse count. Full controlled-V
spectral readout must control this entire emitted word, including every
current addition and every field phase.

`composition.py` consumes the exact retained source response and literal
binary64 controls. Integer outward trigonometric enclosures certify the
Gaussian preparation angles; phase-exact Pauli templates certify every
spectator string, source primitive and field word. The matter/current layer
contains1172864 native pulses. The one-face field member adds36562 pulses.
Its complete V has source Stinespring norm error at most1.990772e−10; the
homogeneous208-face word with the same coefficients has bound3.466855e−10.
These bounds retain arbitrary quantum references and all export complements.
They assume the declared exact response and coherent code services, rather
than uncontrolled laboratory noise.

The initialization program resets3203 logical source/work carriers and
emits the actual one-particle preparation at port0 with uniform event
amplitudes and F=0. A separate one-hole preparation fills every physical
mode except the initial hole and applies the corresponding real Gaussian
word. Its reference is G=12 ones−delta0, so its relative field current has
the opposite sign. These are separate preparations of the same source;
they are not a simultaneous particle-antiparticle pair or a ground-state
identification.

The full terminal instrument reads all2048 occupation bits and1152 field
bits with actual native pointers. Every pointer ready-reset timeout,
source/export failure and decoder complement remains an outcome. The
minimum reset success probability over all admitted native initial columns
bounds every conditional reset, so the union bound does not assume that
unknown initial carriers are independent. The complete three-step record
law differs from the ideal law by total variation at most3.296486e−5.
`PREPARATION_PROGRAM.json` gives the initialization and readout words;
`COMPOSITION.json` records every error and reset factor. Factorized operator
identities establish the full circuit action. A dense2^3203 state array is
not enumerated.

## Executed interaction member

`autonomous.py` executes alpha=0.17 and beta_0=0.41 with all other face
couplings zero, together with the separate zero-coupling interventions.
It retains the complete one-particle state on all physical positions,
event labels and every reachable Z16 field configuration. No small amplitude
is discarded. The plaquette orbit update includes every one of its16
outcomes. The source is the one-particle restriction of the same full
mode/field program; the one-hole control uses its charge-conjugate
preparation. All reported branches satisfy Gauss and the complete norm.

After three steps the coupled and zero-field programs differ in position
probabilities by total variation0.0004197274. The electric expectation
differs by1.00395795. The same native current decoder reports a nonzero
edge-current response. Both position and electric responses exceed the
complete source/readout error bounds of the two programs: the electric
observable lies in[0,9216], including all modular labels. The records are
joint source/event/field measurements, not a mean field fed to a separate
particle Hamiltonian. The electric step changes later matter interference,
and matter transport changes the same retained field.

The reported current response is an operational transport observable. A
force dictionary additionally needs the common length/time and inertial
mass attachment. Pair forces need one shared two-excitation state and
reference. The one-hole comparison supplies neither premise. The numerical
member has one active face; its results are not relabeled as a homogeneous
three-dimensional run.

## Transverse harmonic response

The homogeneous field's Bloch face symbol is reconstructed from the actual
size-three routed graph, with36 edge and26 face coordinates per cell.
Its zero-momentum incidence has rank11 and its face matrix has rank22.
The three-dimensional joint kernel of B(0) and A(0)* is the harmonic edge
space. Projection of the three coarse edge chains into it gives vectors
w_a with exact Gram matrix (2/3)I. A rational Schur calculation yields the
effective curvature symbol

\[
 (2/3)(|k|^2 I_3-k k^T)
\]

in an orthonormal harmonic frame. There are two nonzero transverse branches;
the third direction is longitudinal. This is derived from the local face
attachment and its actual incidence, rather than inserted as a continuum
execution matrix. Finite Bloch controls reproduce the coefficient and its
lattice corrections.

Write U_e=exp(i theta_e), fixing the sign of the unwrapped canonical angle.
The harmonic field kick with order U_E U_B has
E_next=E−beta A A* theta and theta_next=theta+alpha E_next. On a curvature
mode lambda, cos(omega)=1−alpha beta lambda/2. Its leading transverse speed
per cell and field step is sqrt(2 alpha beta/3). Matching another matter
speed therefore constrains the declared product of couplings and number
of field substeps; it is not a prediction of a fine-structure constant.

The exact infinite-lattice convolution row bound controls the curvature
spectrum on every periodic refinement, before any small-volume folding.
The optional Dirac-matter adapter uses eight field substeps per matter macro,
alpha=1/16 and beta=3/8. Its leading speed is1 and its harmonic ultraviolet
stability bound is158(3/128)=237/64<4. This is a parameter and clock adapter
for a separate matter amendment, not a source equivalence to the base walk.

On a retained state with all admitted plaquette angles bounded by r_f and
leakage probability p, replacing each cosine kick by its quadratic version
has vector error at most

\[
 \sum_f |\beta_f|r_f^4/24+2\sqrt p.
\]

The bound follows from the cosine Taylor remainder and unitarity on the
complement. Multistep error is the sum evaluated on the actual intermediate
states. Finite electric-register wrap errors are additional. A controlled
quantum Maxwell continuum requires a preparation with certified angle and
flux tails, their dynamical bounds, and clock/native-precision refinement.
The harmonic symbol and one-kick bound do not by themselves establish that
full interacting limit. Particle species, observed masses and nuclear
forces are outside the claims of this finite construction.
