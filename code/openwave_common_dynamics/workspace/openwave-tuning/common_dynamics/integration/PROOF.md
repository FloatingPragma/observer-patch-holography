# One massive source with field, position and clock records

The candidate is `native-local12-mass-fermion-plaquette-v3-proposed` at
`h=1/8`, mass parameter zero or one, electric coefficient `.17`, and magnetic
coefficient `.41` on face zero. The other 207 magnetic coefficients vanish.
The prepared state has one particle at vertex zero, uniform twelve-event
amplitudes and zero field. The observation interval contains at most three
complete macros. These settings are architectural inputs.

The frozen source macro is

`V_m = U_E U_B Gamma(S)_F Gamma(R exp(-i h m R))`.

Its native implementation prepends the declared on-site Gaussian mass word
to the complete matter/current/field word. The native clock is a separate
prepared plus state with a tick `exp(-i h Z/16)` after that entire macro.
Thus the ideal joint macro is `V_m tensor exp(-i h Z/16)`. Its measured clock
phase after three macros is `3/64`, and its declared time is `3 tau0/8`.
The clock is not a self-starting particle rest oscillation or a measurement
of elapsed processor time.

The stated time and source comparison concern completed macros. Native
compiler prefixes and global controller/readout access have their own
resource and locality domains. A bound on the completed spatial update
does not prove Lorentz locality of every intermediate encoded pulse or
give a spacetime interpretation to its service duration.

## Allocation and complete program

The source uses wires 0 through 3202. The seven position flags use 3203
through 3209. The clock uses 3210. The pointer's two-bit auxiliary code uses
3211 and 3212, with a separate native ready-carrier reset before every read.
There is no spectral probe in this experiment. The public allocator rejects
overlap among these regions, a different modulus or source, two-particle
input, and mass/clock settings outside the frozen finite members.

`joint_program.complete_program` emits the complete chronological program:
native configuration resets, source preparation, clock preparation, the
mass/field/tick macros, position compute/read/uncompute, every field-bit
read, and clock analysis/read. Every native read includes its reset, six
source pulses and the actual C or E commit on the ready carrier. All
source pulses use the producer's complete coherent loading/export service,
including its complement outcomes. A reset timeout stops the program and
retains the fine tape and final carrier. A source/export failure retains
its original record and the last fully completed macro count. No guessed
position or clock result is substituted after failure.

The public wrapper `verified.complete_program` checks every frozen producer
hash before generation and again on completion, then emits the explicit
terminal decoder declaration. `verified.error_contract` applies the same
guards to the bound. A one-byte source mutation is executed through that
digest admission function and rejected. Producer files are not modified.
The exact central contract and upstream criterion snapshot read by the
audit are retained under `consumed/`, with an explicit original-path map
in the manifest. They are historical context inputs, not imported dynamics.
Producer and imported-source hashes retain their exact executable meaning.

The seven position bits are address parities of mode occupations. On the
one-particle sector they implement the projectors onto position, retaining
all within-position event and field coherences. Compute/read/uncompute
returns the flags to zero on the ideal branch. On a general Fock state
the same circuit measures seven parities. It has no two-particle position
interpretation. All 128 address outcomes are retained; the 32 nonphysical
addresses are decoder complements. Field readout is an additional stated
dephasing intervention on all 1152 field bits. Logical event/padding
degrees not read by this apparatus are retained in the final quantum state.

The pointer's complete local matrix-unit identity supplies the two
Lueders maps on any target bit with an arbitrary untouched reference.
Wire relabelling changes neither its norm bound nor its reset cost.
Therefore positioning its auxiliary code above the flag and clock wires
repairs the alias without introducing another measurement premise.

## Complete-error composition

Let `epsilon_mass` be the frozen complete source/preparation/position/field
bound from `massive_coupled.error_contract`. It includes 3203 blank native
resets, seven blank flag resets and 1159 ready pointer resets. Let `p0`
and `pr` be the minimum success probabilities of the two actual native reset
controllers. Let `delta_clock` be the certified native preparation, tick,
analysis and pointer bound from `phase_clock.exact_word_bound` at the same
`h` and macro count. The joint complete-instrument bound is

`epsilon_joint <= epsilon_mass + (1-p0) + (1-pr) + delta_clock`.

The complete experiment therefore has 3211 blank native resets and 1160
ready native resets. The rank-one conclusive fine reset Kraus operators
give a diagonal success effect with least eigenvalue equal to the minimum
of its diagonal. This operator statement applies to arbitrary coherent
and correlated native inputs. Conditional first-failure bounds compose by
a union bound, without multiplying unconditional marginal success rates.
The original private labels/tapes are retained; they are marginalized only
at the explicitly stated comparison with an ideal source preparation.

The coherent-word errors are operator/Stinespring bounds with all native
export complements. Their tensor extension and the complete measurement
contract preserve arbitrary quantum references. Adding them by the
triangle inequality gives the displayed complete terminal-record bound.
There is no conditional-success normalization and no inference of a
physical noise model from a floating matrix residual.

The optional spectral experiment would require controlling every source
factor, including the field-current additions, mass phases and field
kicks. Controlling only the bare walk on a plus-probe state violates the
same Gauss sector with probability one quarter after one macro. That
rejecting control is separate from this clock attachment, whose clock
does not control the source.

## Independent executions

The audit reconstructs the mass kick through the exact rank-one spectral
projector of the twelve-event reflection. An independent complete sparse
implementation executes the native graph transport, current increments,
one-face magnetic orbit and centered electric phase. It compares every
amplitude and every archived event probability for both mass settings at
steps zero through three. No small-amplitude threshold is used.

The separate current review reconstructs all 288 signed moving-slot
probabilities after the same mass/reflection word and verifies ordinary
mean continuity `Delta<n>+B<j>=0` at each step. It does not replace the
integer current by differences of centered modular field means at a wrap.
The corresponding event/current read is an additional stated instrument;
the terminal position/field/clock readout does not silently measure the
unobserved event index. Omitting the current gives a rejecting control on
the admitted source preparation.

The clock's native one-bit word and full C/E pointer instrument are
executed separately and joined to the same source law. The product is an
ideal source/clock law for the declared product clock preparation. The
actual full native program differs by the complete bound above. A
correlated clock preparation requires its own joint law and is not
represented as this calibration. Finite-shot confidence intervals are a
separate statistical statement.

The audit records an ordered native block program, its digest, exact
repetition counts and maximum wire. It reuses the frozen complete field
macro, the emitted on-site mass word with every wire map, and literal
clock/readout words. Repeating a certified block is exact ordered
concatenation; no intermediate state reset or commuting-factor exchange is
introduced. `complete_program` exposes the primitive expansion, while the
certificate avoids hashing millions of identical pulse records again.
Factorized local operator identities establish the circuit on the full
source. No dense Hilbert space with more than three thousand qubits is
enumerated, and no completed flat-stream hash is claimed.

## Field cutoff and source-family boundaries

For the centered cyclic field of modulus `m`, initial support cap `E0`,
total per-edge matter-motion budget `J` and total absolute magnetic
coefficient `Lambda`, the root cutoff theorem uses
`d=m/2-E0-J>0` and the norm bound

`min(2, 2 sum_(k>=d) Lambda^k/k!)`.

An independent rational partial sum and later geometric remainder bracket
the exponential tail inside the producer's bound. Each low-order term with
`k<d` stays strictly below the positive cyclic wall and agrees with the
integer-field process at every prefix. Electric phases are diagonal
unitaries on both domains, so their unbounded integer generators do not
invalidate the finite ordered-product proof. The constant factor two
bounds the two omitted tails. The same subspace operator bound extends
to references; initial outside-support probability adds at most twice its
square root. The claim concerns complete terminal instruments agreeing
under the centered embedding.

A state at the register wall falsifies uniform convergence. A one-face
coefficient sum cannot bound all 208 active faces. The comparison target
is an integer compact rotor; there is no quantum Maxwell vacuum or
spatial continuum theorem in this cutoff result. Native costs and errors
for modulus 16 do not transfer to a larger register by changing a label.

The Dirac proposal is a different matter source. The audit independently
checks its Clifford relations, free Lie-product bound and actual signed
route currents. It does not substitute its spin result into the scalar
mass source. Its full preparation/native-error contract, physical rotation
apparatus, interacting canonical-field limit and lattice species control
are separate claims. A conditional Pauli `g=2` identity is not a native
magnetic-moment measurement.

The accompanying 31-criterion map records the supporting finite domains
and the required physical connections. No official platform status is
awarded by these integration checks.

## Shared normalization

The units audit consumes the frozen unit-flux basis and verifies it against
the actual cell symbol: its cycle Gram matrix is `(2/3)I`, its three cut
crossings are the identity, and the dual angle basis makes the symplectic
pair canonical. Conditional on the stated canonical homogeneous field
limit, `N` substeps give electric coefficient `A_e=N(2/3)alpha`, magnetic
coefficient `A_b=N beta`, speed squared `A_e A_b`, and dimensionless charge
coupling squared `sqrt(A_e/A_b)`. The latter is gauge coupling, not magnetic
Lande g.

The Dirac speed setting has `(A_e,A_b)=(1/3,3)`, hence speed one and charge
coupling squared one third. A different pair `(1/6,6)` has the same speed
and half that coupling. Unit conversion cannot turn these two dimensionless
couplings into the same prediction. The real static Coulomb coefficient is
conditional on an admitted static/canonical field response, not an executed
two-charge force from the one-face source.

The scalar speed setting has conservative ultraviolet product bound
`79/16>4`. That bound alone gives no full ultraviolet stability certificate;
it does not prove instability. The eight-substep Dirac setting has
`237/64<4`. Definite zero electric field, a finite cycle-character state,
and a quantum Maxwell vacuum have different preparation domains.
