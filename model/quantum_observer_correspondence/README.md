# Coherent quantum histories on overlapping native clocks

This package extends the single-carrier observer correspondence to full
quantum output states, entangled references and a declared coherent
interaction. The native and effective descriptions agree exactly, including
their local clock marks, feedback, unresolved waits and terminal measurement
disturbance. The result is finite and conditional. It does not identify an
interacting relativistic field theory, physical duration or measured spacetime.

The nearest [diagonal correspondence](../observer_correspondence/README.md)
uses the same 66-state source, thirty-seam schedule and overlapping regions,
but discloses definite initial occupations. **This experiment changes that
preparation/read contract.** It omits the initial disclosure, which would
dephase the coherent probe, and selects feedback from a measured current
instead of an unmeasured occupation. Its 1,442 records are not the separate 5,768-word experiment with an arbitrary coherent input substituted.
Fields here are accessible quantum observables and instruments. They are not
simultaneously assigned classical occupations inferred from an unread state.

The [native operator construction](../native_quantum_history/README.md)
provides the existing signed-difference reduction on all 4,356 operator units.
We reuse it to check the unread source, then split the actual Kraus maps into
remote and current-marked incident parts. The
[joint-clock process](../joint_clock_process/README.md) has complete reference outputs on a separate four-state support with a different schedule. The attachment here is the full single-carrier quantum domain and its exact
finite-cap representation, with independent joint clock/response replay.

## The declared experiment

The source is the binary total-two subspace of twelve occupation ports,
of dimension 66. Native C repairs use an independent uniform choice among
the thirty source seams and a fair placement coin. The coin and seam are
hidden CP alternatives; different native current branches are never added
coherently. All original C branch semantics remain intact.

Observer A accesses `(0,1)` and B accesses `(0,5)`. Their union is `(0,1,5)`.
A local incident clock ticks when the selected seam meets its region,
including zero-current attempts. A mark contains both clock increments and
the three signed local occupation increments. Remote repairs continue while
both clocks wait. The joint record collector and the attempt-cap service are
declared; physical message delivery and simultaneity are not inferred.

1. Prepare `(|ports 0,1> + |ports 0,5>)/sqrt(2)`, with no initial occupation
   read. A separate input entangles these two basis states with a two-state
   untouched reference.
2. Request one union incident, with at most two native attempts. A timeout
   ends native evolution but remains in the final output.
3. Following completion, perform one native C repair on `(0,1)` if the
   recorded change at port 0 is positive; otherwise use `(0,5)`. Retain its
   current and count its incident tick for both observers.
4. Apply the declared ideal union intervention
   `U_gamma = R_(1,5) (-1)^(gamma n0 n1)`, for fixed `gamma=0` or `1`.
   In the pair basis `(|10>,|01>)`, R has columns `(3/5,4/5)` and
   `(-4/5,3/5)` and fixes `|00>,|11>`. The occupation interaction precedes
   the rotation.
5. Request a second union incident with the same cap, retaining any timeout.
   Measure port-1 occupation with both projectors on every final branch,
   including first-interval timeouts. Keep the entire unnormalized
   source/reference density matrix for every record.

The phase and rotation are supplied coherent apparatus, evaluated by exact
rational matrices. They are not inferred native dynamics or contracted from
emitted source pulses. Each successful first interval incurs one phase
operation, including the identity control, and one rotation; every branch
incurs one terminal measurement. They add no native incident ticks. There
are at most four scheduled attempts and one extra native feedback repair.
Expected scheduled attempts are `328/125`; expected feedback repairs are
`16/25`. No apparatus compilation, physical duration or service-error
certificate is supplied.

The phase is a genuine two-occupation operation on this declared sector:
its phase product on occupations `{0,1}` and `{2,5}`, divided by that on
`{0,5}` and `{1,2}`, is `-1`. Every diagonal one-body occupation phase has
ratio `+1`. This rules out that particular reduction to one-body phases;
it does not derive the interaction from repair or establish a field theory.

## Quantum correspondence and the allowed matrix cones

Write A for the remote CP map per native attempt, and B_m for the CP map
of incident mark m. The code constructs their actual C Kraus numerators.
It checks the effects on every matrix unit, including off-diagonal units:

```text
Tr A(E_ij) = (3/5) delta_ij,
Tr sum_m B_m(E_ij) = (2/5) delta_ij.
```

Thus survival through h attempts is `(3/5)^h` for every preparation and
untouched reference. Set p=`2/5`, D=`I-A`, and define

```text
R = p D^-1 = p sum_(n>=0) A^n,       G = D/p,
e = R(rho),                        rho = G(e).
```

The geometric CP series converges in complete trace norm: each A power has
trace effect `(3/5)^n I`. Its sum proves R is CPTP, including reference
extensions. Positivity is not inferred from positive scalar entries of an
inverse matrix. The exact inverse is computed in connected operator blocks
of dimensions 1, 7, 8, 9 and 36. Together they cover all 4,356 matrix units;
the largest off-diagonal signed-difference sector has only ten units before
the remote split. No approximate inverse or spectral truncation is used.

For cap h, the actual and effective incident maps and the timeout map are

```text
K_m,h = B_m sum_(r=0)^(h-1) A^r,
F_m,h = D^-1 B_m (I-A^h) = R K_m,h G,
Q_h = A^h.
```

The completed transfer, without a cap, has effective maps `D^-1 B_m`, which
are CP. **The capped F_m,h are not ordinary positive maps on every positive
effective matrix.** For an effective input `|ports 0,2><ports 0,2|`, the
cap-one mark with B alone ticking and port 5 gaining one has output
population `-1/720` at `|ports 0,5>`. G itself also gives a population
`-1/24` on this positive input. The input lies outside the admissible image.

At each reference dimension k, use the transported matrix cone

```text
C_eff,k = {(R tensor id_k)(rho) : rho >= 0}.
```

This is a proper subcone of the ordinary positive matrices. G is a positive
decoder on this cone; the capped maps preserve it at every reference
dimension. There is no assertion that G is a CP channel on the ambient
matrix algebra, nor that the effective coordinates factor into independent
local quantum systems.

Allowed local intervention branches have number-preserving Kraus operators
supported on the union ports. They commute with each remote Kraus operator
on the full twelve-port tensor space and preserve its total-two sector.
Consequently their superoperators commute with A, R and G. The executable
checks 864 such operator commutators for its actual apparatus, detector and
feedback branches, as well as the apparatus's exact unitarity. Total-two
space itself is not asserted to factor as a union tensor remote space.

Initialize the effective history with R(rho), apply the F maps and the same
admitted interventions, and decode the terminal state with G. Adjacent G/R
factors cancel. Induction gives the identical unnormalized output on every
finite adaptive record word. Extending the identities on all operator units
by complex linearity and tensoring with the identity proves the result for
arbitrary coherent mixed states and finite untouched references. The
retained executions use real rational coefficients; that executable choice
does not restrict the complex-linear map identity to real states.

A cross-boundary rotation on `(1,2)` fails the commutation premise, with
an exact entry `4/375` in `A U - U A` at the superoperator level. Such a
service would need its own transformed intervention dictionary. It cannot
reuse the untransformed local service by declaration. Likewise R(rho) is
essential for nonstationary input: replacing it with rho changes the
decoded timeout state. No boundary reset or selected hidden state is used.

## Complete outputs, entanglement and censoring

All 1,442 output words, including zero-probability words, are retained for
each interaction and reference setting. The interaction-on unreferenced
probe has 243 positive words. The native contraction, transported quantum
representation and independent endpoint/coin/current oracle agree exactly.
The [clock-response audit](CLOCK_RESPONSE.md) extracts a joint observable
from these same outputs and includes matched-marginal and dephasing controls.
Its clock marginal is invariant by construction, so the joint law tests the
correlation between seam records and the detector.
It uses no separately fitted clock or field model.

An additional entangled input is
`(|ports 0,2>|0> + |ports 0,3>|1>)/sqrt(2)`. Its full decoded output agrees
and leaves the unconditional reference marginal `I2/2`. The first-timeout,
detector-zero branch has probability `9/25`. In its partial transpose, the
unit vector `(|ports 0,2>|1> - |ports 0,3>|0>)/sqrt(2)` has expectation
`-17/200`, or `-17/72` after conditional normalization. Separable states
have positive partial transpose because transposition preserves positivity
of each product factor. This exact negative value proves that source/reference
entanglement survives in a retained unresolved wait. The complete output
cannot be replaced by scalar record probabilities or a classical mixture.

For two requested incidents the probability of any timeout is exactly
`1-[1-(3/5)^h]^2`, giving `369/625` at h=2. This remains uniform over all
inputs and references under the declared intervening instruments. On a common
output alphabet that retains timeout flags, the capped and completed
experiments have complete-output trace distance equal to that probability:
the completed instrument adds CP histories omitted by the cap, and the capped
instrument puts their equal total mass into distinct timeout blocks. No
successful subset is renormalized. Larger caps reduce this mathematical
censoring error; they add no physical clock calibration.

This closes a finite coherent/reference correspondence interface. It leaves
spatial refinement, Lorentzian microcausality, interacting continuum dynamics,
physical action normalization, physical source/apparatus identification and
comparison with independent observations as separate obligations. Native repair is the stated instrument; no unitary Hamiltonian evolution is
silently substituted for it. No federation is required for this local result.

```sh
python -m model.quantum_observer_correspondence --check model/quantum_observer_correspondence/results.json
python -m pytest -q -W error model/quantum_observer_correspondence
```
