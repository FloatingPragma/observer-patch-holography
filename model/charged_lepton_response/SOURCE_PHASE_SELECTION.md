# Source phase and the C/E distinction

The specified coherent program has phase `z=1`, or `theta=0`. Its exact
positive action spectrum has ratios `1:1:(34+24 sqrt(2))`. A common mass
factor cannot turn this repeated pair into the three unequal charged-lepton
masses. The program is fully specified and its phase can be read without an
inverse experiment. The finite C/E geometry and measured repair records do
not supply a theorem replacing this configured phase by a nonzero value.

[`phase_selection.py`](phase_selection.py) verifies the actual half-turns,
their action on the coherent family operator, the original memory duties,
and the operator class preserved by ordered feedback. The assertions
concern these specified sources and interfaces. They are not a universal
impossibility theorem for phase selection by another source construction.
The [phase readout](PHASE_AND_PREMISES.md) supplies the complete experiment
and interval inverse for a given coherent phase.

## Actual half-turns and orientation

The geometry is read from
[`native_family_record_memory.py`](../common_source/native_family_record_memory.py).
Its C half-turn is the unique native icosahedral rotation exchanging ports
2 and 3. Its E half-turn exchanges ports 2 and 10. Let C and E also denote
their real orthogonal restrictions to the actual three-dimensional family
frame. They obey

    C^2=E^2=I, R=EC, R^3=I, tr R=0,
    CE=R^-1, CRC=ERE=R^-1.                         (1)

The executable checks these identities both as the literal 12-port
permutations and on the family frame. The generated subgroup has six
elements. The full native rotation group has 60 elements. The coherent
source uses an orthogonal change of family frame taking R to the regular
three-cycle P; equation (1) is independent of that change of coordinates.

For a=√(10/72), b=√(5/72) and a supplied z=exp(i theta), set

    C_theta=a I+b(z R+conjugate(z) R^-1),
    Y_theta=C_theta^2/3.                            (2)

Either literal half-turn sends C_theta to C_-theta and Y_theta to Y_-theta.
Their unordered eigenvalues are equal. Reversing the order EC to CE therefore
cannot choose between two different mass ratios. The positive-root chamber
is `|theta|<=pi/12` modulo `2pi/3`; generic interior phases have three distinct
positive action eigenvalues. In this chamber the signed-root identity is
the physical positive-root Koide identity for the configured action.

The three active Fourier outcomes replace z by `z exp(-2pi i m/3)`. Their
unordered spectra are also equal: this operation permutes the three roots.
Each outcome retains its own input-dependent effect and original probability.
Conditioning on an active label does not select a continuous phase and does
not justify replacing all three maps by the flag-zero map.

## Covariance is stronger than a change of orientation

Every Hermitian matrix commuting with R has the form

    H=a0 I+br(R+R^-1)+i bi(R-R^-1),                 (3)

with real a0, br and bi. The C3 commutant has complex dimension three; its
Hermitian part has real dimension three. Equation (1) shows that invariance
under either half-turn forces bi=0. The invariant Hermitian subspace has
real dimension two. With `P_axis=(I+R+R^2)/3`, it is the span of I and
P_axis, so its two-dimensional orthogonal subspace has a repeated eigenvalue.

The same restriction applies directly to the squared action in (2). Its
coefficient of R has imaginary part

    Im beta=(2b/3) sin(theta) [a-b cos(theta)].      (4)

Since a=√2 b>b, half-turn invariance forces sin(theta)=0. Within the positive
root chamber and the fixed R convention, this is theta=0. The phase-pi root
has a negative eigenvalue and does not instantiate positive-root balance.
Thus imposing full C/E covariance cannot select a nonzero phase while
retaining the three distinct balanced eigenvalues.

Full A5 covariance is stronger. One native rotation moves P_axis, as checked
on the actual carrier. An invariant matrix in the two-dimensional space
above must then be scalar. This does not imply that the actual reset-memory
coefficient is scalar. That source includes a distinguished reset port and
does not have the additional covariance premise. Its always-duty coefficient
is real symmetric and does not commute with the actual R.

## Static duty and orientation are different operations

The original source has C event `(coin=0,raw=1)` and E event
`(coin=0,raw in {4,5})`. All 128 raw labels, including zeros, and all 12 old
memory inputs are retained. The base probabilities are p=1/12 and q=1/6.
The destination applies C on a C event, then E on a joint event; other
records and core destinations reset to port 2. Reconstructing the stationary
projector coefficients from these 1,536 labelled rows gives:

| Operation | Exact consequence | What it can discriminate |
|---|---|---|
| Always duty | Positive rank 3, original trace 1 | A specified static action |
| C duty | Positive rank 3, original trace p=1/12 | Distinct from E duty |
| E duty | Positive rank 2, original trace q=1/6 | A zero action singular value |
| CE duty | Positive rank 2, original trace pq=1/72 | A zero action singular value |
| Reverse the word EC to CE in (2) | theta maps to -theta; unordered masses equal | Orientation-sensitive observations, not unordered masses |
| Change the active Fourier label | theta changes by a multiple of 2pi/3 | Root labels, not unordered masses |
| Exchange all C/E event laws and reset wiring | A different complete source law | Requires its own source calculation |

The rank statement for the four positive-rate duties also follows from
[`native_flavor_orbit/algebra.py`](../native_flavor_orbit/algebra.py),
`memory_certificate()`, for every 0<p,q<1. This includes the declared gated
source. The numerical p and q in the table belong to the ungated base
control; they are not substituted for the gated source rates. Rank two
cannot supply three nonzero action masses without an additional contribution.
Rank three alone does not imply the observed mass ratios or Koide balance.
The four current-gated static coefficients fail the balance test in
[`KOIDE_ATTACHMENT.md`](../standard_model_source_map/KOIDE_ATTACHMENT.md).
Neither static duty is identified with the coherent source (2).

## Ordered memory has a potential, with a different operator placement

The actual feedback pulses use real rank-one projectors P_i and generators
`X tensor P_i`. Every ordered unitary word obeys

    [U,X tensor I]=0,
    (Z tensor I) conjugate(U) (Z tensor I)=U.        (5)

This follows for each exponential and is preserved by ordered multiplication.
It does not assume commuting pulses. The executable verifies all twelve
native projectors and an explicitly noncommuting two-pulse word. For a
differentiable unitary evolution through I with the same identities, its
Hermitian generator has the form

    H=I tensor iA + X tensor S,
    A real antisymmetric, S real symmetric.         (6)

Indeed, commuting with X leaves only I and X components; differentiating
(5) makes the I component purely imaginary and the X component real.
Hermiticity gives the stated transpose symmetries. Thus an imaginary family
potential in the I block cannot be renamed as the imaginary part of the
X-block action. At generic phase, X tensor Y_theta violates (6) in this
fixed native frame by (4). Transporting the complete frame transports this
constraint. Equation (6) excludes neither arbitrary alternative dictionaries
nor an isospectral real representation, and it does not assert that the
finite unread memory channel is unitary.

The [all-lag source theorem](../particle_memory_response/PROOF.md) supplies
an actual ordering correction. For the base always duty it gives

    H_ord=I2 tensor A,
    A=i (603/6666724)
      [[0,3+sqrt(5),0],[-3-sqrt(5),0,-2],[0,2,0]],
    G1=-i[H_ord,.]-(1/2) sum D_ij[H_i,[H_j,.]], D>=0. (7)

The same actual R in (1) gives the exact orientation identity

    A=i kappa(R-R^-1), kappa=603(1+sqrt(5))/6666724,
    [A,R]=0.                                        (8)

Consequently, in the source's orthogonal cyclic frame `R O=O P`,
`O^T A O=i kappa(P-P^T)`. The ordering correction supplies an actual
oriented C3 generator, with its coefficient determined by the native law.
Because C_theta and Y_theta are polynomials in R, they commute with A.
Thus conjugation by exp(-it A) leaves both unchanged and cannot turn C_0
into C_theta at a nonzero phase. This positive identification is distinct
from the static always-duty coefficient, which does not commute with R. Equation (8) does not specify
a duration, a branch-phase controller, or a transfer from the I block to
the X block. It therefore does not select a scalar Yukawa phase.

This is a positive result from the entire correlation sum. It is not a
vanishing-commutator argument. The source theorem compares the unread field
on a stationary product preparation, with arbitrary reference, to
`exp(t(G0+hG1))` with a rational O(h^2) bound at fixed t=Nh. Its second-order
comparison is not a comparison of every classical archive to that semigroup.
The complete leading weak limit has the separate archive bound in
[`FEEDBACK_LIMIT.md`](../standard_model_source_map/FEEDBACK_LIMIT.md).

At fixed source gain and t, h H_ord and the h D dissipator both vanish as
h tends to zero. Retaining a second-order term by changing gain changes the
leading term and the diffusion at the same time. Removing the diffusion,
subtracting a leading generator, or converting the resulting potential to
a branch phase requires additional specified controls and their error bound.
Taking a logarithm of the full finite memory channel supplies no shortcut:
the source channel is noninvertible; its recurrent classical block also has
the negative eigenvalue -p(1-q). The slow-block expansion used in (7) has its
own stated preparation and reduction and is not that full logarithm.

## Geometric phase and finite theorem premises

The actual coherent-line geometry in
[`native_mass_holonomy.py`](../common_source/native_mass_holonomy.py) gives
the spin-half Bargmann phase pi/10 around an oriented icosahedral face.
The spin-one coherent line doubles it to pi/5. Integer closed seam walks
have the corresponding integer phase lattice. A separately specified
three-equal-link rule assigns one third of the loop phase to a link;
the positive balanced unordered classes are 0, pi/30 and pi/15. The exact
geometric phase is not automatically the branch phase in (2). Neither the
choice of coherent line and loop nor the equal-link attachment is supplied
by the measured C/E memory law.

The Lean theorem
[`SeamU1HolonomyClassification`](../../Lean/Screen/SeamU1HolonomyClassification.lean)
classifies the actual seam connections modulo gauge by 19 chord holonomies.
Its holonomy map is surjective. Its unique tree-gauge representative is
unique for a supplied connection, not a selection of a connection from
repair data. [`HolonomyInterference`](../../Lean/Screen/HolonomyInterference.lean)
relates a supplied path character to the relative loop phase. It does not
choose that character or the connection.

[`KoideCirculant`](../../Lean/ObserverPatchHolography/KoideCirculant.lean)
proves the balanced quotient independently of phase, under its positivity
and modulus hypotheses. The finite tracial-GNS balance and its operator
attachment are stated in [the premise crosswalk](PHASE_AND_PREMISES.md).
A phase-independent theorem cannot by itself select theta. The conditional
charged-response rigidity construction in
[`derive_charged_source_law_rigidity_conditional.py`](../../code/particles/leptons/derive_charged_source_law_rigidity_conditional.py)
requires additional ranks, trace state, path table, clock and ancestry
premises. Its baseline offset `zeta=theta-2/9` stipulates that phase origin.
The code does not derive those physical premises from the C/E source.

## The retained pixel-dependent MCPR construction

The minimal complete public response (MCPR) law gives a definite conditional
phase from declared pixel capacity P and unified response coordinate u:

    chi=(u^2/512)/(1-u/77),
    zeta=[(P-phi)^2/(21*pi)]/[1-u/27-u^2/135],
    delta(P,u)=2/9+zeta, rho=sqrt(2) exp(-chi),
    r_j=1+rho cos(delta+2pi j/3), m_j proportional r_j^2. (9)

These are the fixed points of the retained affine response equations with
register dimensions 512,77,21,27,5. Their implementation is
[`derive_charged_mcpr_completion_conditional.py`](../../code/particles/leptons/derive_charged_mcpr_completion_conditional.py).
The exact denominator conditions `1-u/77>0` and
`1-u/27-u^2/135>0` give delta>=2/9, with equality precisely at P=phi.
For u>0, chi>0 and rho<sqrt(2). When all roots are positive, their quotient
is `(1+exp(-2chi))/3<2/3`. Thus this retained response includes a deformation
of both phase and balance. The balanced angle below 2/9 obtained by the
separate one-input calibration is not delta in (9). Copying only delta
into (2) is a specified alternative controller; it changes the MCPR
amplitude and is not equality of the two source actions.

There is a finite dynamical realization of (9), rather than only a table
of ratios. The exact inverse construction in
[`mcpr_native_memory_bridge/PROOF.md`](../mcpr_native_memory_bridge/PROOF.md)
selects two admissible branches of complete C-duty probe rates p(P,u),
q(P,u) and a gain g for which the unconditional coefficient satisfies
`spec(g Z_C)={m_j/v}`. All 128 repair records and thirteen-outcome left and
right current probes are retained; the grouped instrument has 512 records
per old memory input. Its correlated weak-pulse theorem realizes the same
coefficient, up to a gauge-trivial family unitary, with a finite error
bound. These chosen rates and gain are part of the configured MCPR source;
they do not equal the fixed pi/5-duty program whose exclusions are above.

The [pixel closure](../pixel_mass_closure/PROOF.md) also certifies two
conditional fixed points, each unique in its stated interval, using this
same charged law. One uses a supersymmetric running proxy with scale
coefficient four; the other uses the minimal census with a specified color
scale coefficient seven. The minimal census with coefficient four is
excluded on the stated domain. None of these statements selects the scale
ansatz, thresholds, response kernel or the MCPR baseline from primitive
repair alone. The forward source-audit decimal coordinates and the
measured-endpoint coordinate branch have distinct ancestry; the latter
inherits measured-alpha input. The retained forward source artifact quotes
`delta=0.2222247634275244202993587776065...` and
`chi=0.00000330489736391485586617533549...`. These are cited conditional
outputs, not values read or rederived by this certificate. Its owning
[receipt](../../code/particles/runs/leptons/charged_mcpr_completion_conditional.json)
and producer identify their source coordinates and provenance.

The MCPR architecture specifies its phase offset, eight registers, path
monomials, signs, amplitude and determinant normalization. Its source files
explicitly record selection after exposure to the charged spectrum, while
its execution reads no charged-mass reference. Readable specifications
therefore determine a conditional source result and a valid finite memory
realization. They do not convert that historically informed architecture
into a prospective source-only derivation of the offset. The distinction
is between an available configured construction and an unproved choice of
that construction, not between existence and nonexistence of any phase law.

A measured event frequency, such as 2/9 in a declared register, is a
probability. Converting it into radians requires an explicit phase-per-event
service and clock; a frequency alone is not a unitary phase. A fixed known
coherent controller can supply such a service and permits the exact action
calculation and readout. The unprovided attachment is a source-derived
controller whose branch phase follows from the primitive dynamics with its
complete outcomes and clock. The geometric and ordered-memory theorems do
not provide this attachment; the MCPR model supplies it conditionally on
its stated architecture. The calibrated chosen-phase comparison is a
separate diagnostic with an empirical input. It is not consumed by this
certificate or promoted to a source-only prediction.
