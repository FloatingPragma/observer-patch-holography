# Icosahedral closure certificate (C19)

Exact representation-theoretic certificate for the icosahedral closure theorem:
the twelve-port screen module, the compact-Lie trichotomy, the rank-five
noncentral corollary, and the oriented face-phase multiplicities.

Referenced as **[C19]** in `extra/compact_proof_of_oph.tex` and as the
closure certificate in `paper/observers_are_all_you_need.tex` and
`paper/recovering_relativity_..._compact.tex`.

## Run

```bash
python3 a5_screen_sm_closure.py          # module identity, trichotomy, face phases, kinetic relation
python3 exterior_sm_completion.py        # exterior matter, anomalies, weak load, deck control
python3 a5_selection_certificate.py      # sharpness data for Cohn-Kumar universal optimality
python3 echosahedral_selector_certificate.py all  # issue #565 source selector, receipt, controls
python3 port_current_inner_certificate.py all     # issue #566 physical port-current algebra, receipt, controls
python3 super_tannakian_matter_lift_certificate.py all  # issue #314 super-Tannakian matter lift, receipt, controls
python3 axis_center_descent_certificate.py all    # issue #567 conditional kernel/lattice receipt and physical no-go
python3 -m unittest discover -s tests -v    # issue #565/#566/#314/#567 regression and adversarial suites
python3 a5_compact_lie_classifier.py     # compact-Lie enumeration
python3 a5_harmonic_decomposition.py     # angular multiplet sequence
python3 bh_log_correction.py             # conditional horizon log-coefficient decision tree
python3 independent_trichotomy_check.py  # independent re-derivation, trusts nothing above
python3 survival_boundary_certificates.py # exact Q0 fixtures and physical-boundary controls
python3 test_audit.py                    # regression suite
```

Requires Python 3.11+ and SymPy. The suite exits 0.

## What is certified (exact; no floating-point fit, no measured number)

| Object | Statement |
|---|---|
| **Echosahedral counting and incidence selector (#565, #625)** | On the declared federation-of-twelve-port-echosahedra branch, twelve primitive central atoms of trace `1/12`, an additional integer atom-counting grammar on the total-12 fiber, and the normalized central-readback Hilbert-Schmidt cost give the unique all-one split with exact quadratic gap `2`. Oriented edge/face incidence and the refinement lineage independently give the unique graph-distance-three antipode, `Aut+ = A5` by a faithful conjugation action on five Klein-four subgroups, six axes, and the exact rank-three Gram frame `G^2=4G`. A covariant half-count candidate exists on the reduced finite-atomic carrier. Its complete A1 operational/refinement data, A2 naturality, and A3 optimizer lift are open, so #625 does not establish full-schema independence. Issue #628 owns the operational source construction. The receipt checks all 60 proper and 60 improper frame determinants, refinement cocycles, arbitrary relabelling equivariance, and the typed negative controls. |
| **Source-bound port-current algebra (#566, #599)** | Given the certified #565 carrier and the hash-pinned semantic response artifact measured from the carrier dynamics (sector structure `1+3+3'+5` with Galois-paired triplets, oriented port-to-vertex frame, potential and rotation response channels with measured signs, defect-port persistence maps), the certificate recomputes every exact artifact claim in `Q(sqrt5)` and derives the map `K : P12 -> u(C^3 (+) C^3)` with twelve-dimensional skew-adjoint commutator-closed image `u(3) + so(3)`, one-dimensional center, derived dimension 11, positive-definite invariant pullback, A5 covariance, an inner A5 action, and naturality along the declared tower and the physical persistence maps. The equivariant intertwiner space is exactly four-dimensional; the construction-model string is control-lane only and the physical source gate is computed from the binding. Nineteen negative controls fail closed. Theorem and proofs: `code/a5_closure/port_current_algebra_certificate.md`. |
| **Source-bound super-Tannakian matter lift (#314)** | Given the hash-pinned source-bound #566 packet: the sixty proper implementers lift exactly to SU(2) with a 120-element lift group carrying a **unique involution** (non-split binary icosahedral double cover — PORT-SPIN-LIFT); the current algebra acts faithfully (rank 12) on the matter carrier `V = C (+) W` and on the selected module; the derived equivariant projector `P_even - P_vac` on the 32-state auxiliary CAR/Fock space has exact rank 15 and realizes `Lambda^2 V + Lambda^4 V` with charge spectrum `{1/6:6, -2/3:3, 1:1, 1/3:3, -1/2:2}`; chirality, conjugation through the invariant top line, all listed realized anomaly traces, even Witten parity (4 doublets), and three exact one-dimensional Yukawa invariant lines are machine-checked; the common action kernel is computed on the genuine simply connected cover `R x SU(3) x SU(2)` and **emitted** as data without forming the quotient. The trace-balanced charge pair `(-1/3, 1/2)` is **derived** by BLOCK-DETERMINANT-BALANCE (anomaly freedom of the realized package forces the balance line, primitive `q = 6Y` integrality and the measured artifact orientation fix the pair), the scalar and channel list are **derived** by the selection scan (admissible scalar charges exactly the conjugate pair `+-3` equal to the weak-block charge), the category typing is **forced** on the realized module by the failing Vec/sVec controls, and the artifact persistence maps are intertwined on the carrier and the Fock realization. The measured spin statistics artifact supplies the transport double cover from source: the exact quaternion lifts of the sixty measured deck rotations, the binary-icosahedral order profile, the centre `{+1,-1}`, the eight-way section-sign exhaustion over all five Klein four-subgroups, and the unique spin structure on the oriented support. The exhaustive 1024-subset anomaly scan (machine-checked in `Lean/Screen/ExteriorSelection.lean`) selects the unordered conjugate rank-15 pair as the unique nonempty chiral anomaly-free selection, with the fermionic-parity grading as an output, and the implementation enumeration forces the Spin/odd-Weyl typing: the gauge centre is excluded by the Lean fermion-parity no-go and the measured `-1` is the unique remaining central implementation. Candidate-class nonemptiness is witnessed without promoting uniqueness; scalar existence and multiplicity are owned by #609. Twenty typed negative controls fail closed. Theorem and proofs: `code/a5_closure/super_tannakian_matter_lift_certificate.md`. |
| **Source-bound axis-centre descent (#567)** | Consuming the pinned #314 receipt and the measured global-form artifact, exhaustive enumeration over the 36 central candidates gives the diagonal common kernel generated by `(omega_3 I_3, -I_2, e^{i pi/3})`; the character lattice has basis `(1,0,-2), (0,1,-3), (0,0,6)` with Smith invariants `(1,1,6)` and exact dual basis `(1,0,0)`, `(0,1,0)`, `(1/3,1/2,1/6)`. Local tensors alone descend through quotients by `1`, `Z2`, `Z3`, and `Z6` (character-residue counts `36,18,12,6`) and cannot select — that stays on the record as the negative control. The selection is carried by measured source data: the order-120 deck action on the incidence-nerve federation, the measured six-axis class group of exactly the kernel order, and the measured flux-sector menu `{0..5}` realized by two-puncture seam witnesses on the base and refined supports, with exact single-puncture impossibility. Realized matter transports single-valuedly through every measured sector (42 monodromy checks) while a fractional singlet obstructs every nonzero sector; the unique global form whose sector menu equals the measured menu is the `Z6` quotient. The Dirac-pairing commutant theorem selects the electric (Wilson) polarization as the unique maximal mutually-local lattice containing the realized lines (12 lattices enumerated, 1 admissible). `h^3=(1,-I_2,-1)` acts trivially and is not fermion parity; screen-scope Spin attachment comes from the #314 source-derived typing with no mixed quotient. Four-dimensional instanton normalization, theta periodicity, monopole dynamics, and laboratory flux measurement stay open named gates. Lean: `Z6Descent.lean`, `ExteriorSelection.lean`; executable receipt schema v4. |
| Vertex module | `chi_P12 = (12,0,0,2,2)`, so `P12 = 1 + 3 + 3' + 5`, multiplicity-free |
| Adjacency spectrum | `det(xI-A) = (x-5)(x+1)^5(x^2-5)^3` → canonical ranks `1,3,3,5` |
| SM adjoint restriction | `ad su(3) = End_0(3') = 3' + 5`; with `su(2) → 3`, `u(1) → 1`, the total is `1 + 3 + 3' + 5` |
| Icosahedral selection | 3 distinct inner products `{-1, ±1/sqrt5}`, spherical 5-design ⇒ **sharp** (`m=3`, strength `2m-1=5`). By Cohn–Kumar (JAMS 20, 2007) it uniquely minimizes every strictly completely monotonic pair cost of squared distance, up to `O(3)`. |
| D-optimal selector | Maximizing `det(F1) det(F2)` at fixed vector/quadrupole trace gives `F1=2 I3`, `F2=(4/5) I5`; the six centered projectors form a regular simplex, hence the unique real ETF(3,6), the icosahedral axes. This is an independent optimality cross-check. The certified echosahedral-federation branch uses the declared counting/cost and incidence packet instead. |
| Compact-Lie trichotomy | Exactly three algebras: `u(1)^12`, `su(2)^2 + u(1)^6`, `su(3) + su(2) + u(1)` |
| **Inner-action closure** | If the `A5` action is **inner**: `dim Z(g) <= 1` (inner autos fix the center pointwise); `Z(g)=0` forces `su(2)^4` whose fixed-space dimension is a multiple of 3, contradicting `dim g^{A5}=1`; hence `dim Z(g)=1`, semisimple dim 11, and `11=8+3` uniquely ⇒ `su(3)+su(2)+u(1)`. **Needs no `W5-NONCENTRAL` receipt.** |
| Angular multiplets | `l=2: 5` (irreducible); `l=3: 3'+4`; `l=4: 4+5`; `l=5: 3+3'+5`; `l=6: 1+3+4+5`. First nonconstant invariant at `l=6`. Frozen forward target FZ-02 with kill bands and timestamped custody; Lean receipt `A5AngularMultiplets.lean`; deterministic receipt via `a5_harmonic_decomposition.py --receipt`. |
| Horizon log coefficient | **Conditional decision tree, not a universal law.** Unconstrained `q`-state cells give `S=K log q` exactly, i.e. `c=0`. Exact 12-port balance: `c=11/2`. 24-slot balance: `c=23/2`. Nonabelian `SU(3)xSU(2)` singlet saddle: `c=11/2`. Full 12-dim gauge singlet incl. `U(1)`: `c=6`. The horizon measure must be derived before `c` is frozen. |
| Rank-five corollary | `W_5` noncentral ⇒ `su(3) + su(2) + u(1)` uniquely |
| Face-phase multiplicities | `m_w = (dim V - chi_V(3A))/3 = (0,1,1,1,2)` on `(1,3,3',4,5)`; minimal nontrivial extension has dimension 3 |
| Missing-four/Higgs no-go | **Rejected.** `End_{A5}(4) = R`, so no `A5`-invariant complex structure exists for a commuting hypercharge `U(1)`. The identity `1+3+3+4+5=16` is dimensional only. |
| Exterior matter witness | On a declared trace-balanced carrier `V=C+W`, the selected non-vacuum even package `Lambda^2 V + Lambda^4 V` is one 15-state chiral Standard Model generation. Its hypercharges, three one-Higgs invariant lines, perturbative anomalies, and weak-doublet parity check exactly. |
| Weak multiplicity | The exterior package contains three colored weak-doublet copies and one lepton doublet, so the per-generation multiplicity is exactly `4`. Its identification with a physical port load requires `PORT-LOAD-TRACE`. |
| Free-deck control | Any declared free `Z6` action on 24 slots has four orbits and a four-dimensional invariant-function space. A physical result requires the deck action and `PORT-WEAK-INTERTWINER`. |

## Why the trichotomy holds (proof sketch)

A compact Lie algebra is reductive, `g = z + [g,g]`, and both summands are
characteristic, hence `A5`-stable. The identity component of the center is a
torus whose exponential kernel is an integral lattice preserved by every group
automorphism, so the `A5`-action on `z` is defined over `Q`. Because `3` and `3'`
are Galois conjugate over `Q(sqrt5)`, a rational submodule contains them with
equal multiplicity. Enumerating rational centers against compact semisimple
dimensions leaves exactly three cases; the other five die because no compact
semisimple algebra has dimension 1, 5, or 7, or because an `A1^4` / `A1^2`
adjoint cannot supply `5` (any homomorphism `A5 → S_4` is trivial, `A5` being
simple of order 60).

## Claim boundary

The exact chain is a finite `A5` module, a compact coefficient bracket,
trace-balanced block integration under coefficient `1/2`, and a six-axis
lattice quotient `Z6`. Central record projectors commute; the full-rank
compact skew-adjoint construction is verified by
`port_current_inner_certificate.py` with the response representation, four
signed coefficients, oriented frame, and physical refinement maps
determined by the semantic response artifact measured from the carrier
dynamics (#599 closure) and recomputed exactly before use. The matter lift
is verified on top of that packet by
`super_tannakian_matter_lift_certificate.py` with a passing physical
source gate: the exact non-split SU(2) double cover cross-checked against
the measured spin statistics artifact, the unordered rank-15 pair selected
by the exhaustive 1024-subset anomaly scan with the parity grading as an
output, the typing forced by the measured centre and section obstruction,
the derived BLOCK-DETERMINANT-BALANCE charge pair, the derived scalar and
channel compatibility, the emitted action kernel, and the candidate-class nonemptiness
witness. Scalar existence and economy stay deferred rows owned by #609 and
never enter the passing gate.
The #567 descent certificate
(`axis_center_descent_certificate.py`) computes the common kernel on every
realized tensor and the character/cocharacter arithmetic of its maximal
effective quotient, and keeps the four-way local non-identifiability on
the record as the negative control. Its physical gate passes at finite
source-model scope on the measured global-form artifact: the deck action,
the six-axis class group, the flux-sector menu with realized-matter
transport consistency, the unique realized-compatible polarization, and
screen-scope Spin attachment. Four-dimensional instanton normalization,
theta periodicity, and laboratory flux measurement stay open named gates.

`survival_boundary_certificates.py` records five exact-small controls at this
boundary: completion non-identifiability from the current source reduct, the
ten-Majorana rank-15 projector, the distinction between a composite response
cubic and a fundamental 1PI cubic, complement-complete refinement with a
hidden-zero-mode counterexample, and a finite exhaustive settlement fixture.
The matching sorry-free logical no-go is
`../../Lean/Screen/PhysicalA5ForcingNoGo.lean`.
These are diagnostic Q0 certificates. They do not emit a
physical screen current, a Spin/locality lift, three chiral families, a scalar
field, or a refinement-stable quantum theory.

The family result is therefore split explicitly. The finite screen algebra
supplies a canonical rank-three candidate band, and the declared completion
criterion selects the least admissible multiplicity inside the window
`3 <= N_g <= 5`. A physical `N_g = 3` statement additionally requires the
rank-45 attachment map and the family descent, exchange, residue, and
refinement receipts listed below.

The exterior calculation supplies a second conditional matter route. It does
not turn coefficient directions into physical currents or select the global
form. It is also not the full even Clifford module: that module contains
`Lambda^0 V` as well as `Lambda^2 V + Lambda^4 V`. Removing the singlet and
excluding other anomaly-free light sectors requires the explicit source-completeness
premise or a source-derived observer-visible discriminator. Selecting `H=W`
as the physical scalar is a separate gate. If the face-phase `A5` action remains
exact on the family multiplicity space, it restricts Yukawa tensors to
`A5`-invariant pairings; general family matrices require a source-derived
breaking, hiding, or forgetting mechanism.

### Exact receipts on the declared echosahedral realization

- **UD12** is exact within the declared realization checked by
  `echosahedral_selector_certificate.py`: the integer total-12 fiber and the
  normalized central readback norm `H(q)=sum q_i^2` are inputs, and
  `H=12+sum(q_i-1)^2` gives the unique all-one minimizer with exact next floor
  `14`.
- **RP-A5** is exact without a downstream representation or measured datum:
  the source-oriented incidence produces the unique distance-three antipode,
  the positive automorphism group is explicitly `A5`, and the distance Gram
  matrix satisfies `G^2=4G`, giving the regular six-axis frame. The theorem,
  data model, equivalence, refinement proof, and countermodels are in
  `code/a5_closure/echosahedral_selector_certificate.md`.

These results use the declared simulator assumption that every local carrier
lineage is a quotient-visible twelve-port echosahedral packet. UD12 also uses
the separately declared integer counting grammar and normalized readback cost.
A covariant half-count candidate exists on the reduced finite-atomic carrier.
Its complete A1 operational/refinement data, A2 naturality, and A3 optimizer
lift are open. Issue #625 therefore carries no full-schema independence
closure, and issue #628 owns the operational source. The results do not prove
that arbitrary OPH carriers must have this type.

- The **conditional PORT-CURRENT-INNER algebraic construction** is verified by
  `port_current_inner_certificate.py`: given the charged-double-triplet model
  and four signed coefficients, the map is full-rank, compact, closed,
  A5-equivariant, inner, and natural along the declared algebraic tower maps.
  These representation data are branch premises, not physical measurements.

- The **source-bound super-Tannakian matter lift** (issue #314) is verified
  by `super_tannakian_matter_lift_certificate.py`: given the pinned
  source-bound #566 packet, the exact non-split SU(2) double cover of the
  sixty proper implementers (PORT-SPIN-LIFT), the faithful current action
  on the matter carrier and the fifteen-state module, the derived
  equivariant projector on the auxiliary CAR/Fock space, realized anomaly
  and Witten checks, exact chirality and conjugation, three Yukawa
  invariant lines, naturality along the declared tower and the artifact
  persistence maps, the emitted action kernel, and the candidate-class nonemptiness
  witness. The charge pair is derived by BLOCK-DETERMINANT-BALANCE, the
  scalar and channels by the selection scan, and the category typing is
  forced by the failing Vec/sVec controls; the kernel emission contract
  and the candidate-class declaration stay typed declarations. Theorem:
  `code/a5_closure/super_tannakian_matter_lift_certificate.md`.

### Conditional receipts and open physical gates

- **PORT-CURRENT-INNER** is source-bound (#599 closure): the semantic
  response artifact measured from the carrier dynamics determines the
  representation, coefficients, frame, and physical refinement maps, and
  the certificate recomputes every exact artifact claim before use.
- **PORT-SPIN-LIFT** rides on the source-bound packet: the non-split
  binary icosahedral cover is forced (unique involution; split lifts are
  impossible on the realized module).
- **BLOCK-DETERMINANT-BALANCE** is derived inside the matter packet:
  anomaly freedom of the realized package forces the balance line,
  primitive integrality and the measured orientation select
  `(-1/3, 1/2)`.
- **AXIS-CENTER-DESCENT** computes the common kernel and the
  character/dual-cocharacter lattices of its maximal effective quotient,
  and attaches them to measured source data at finite scope: the measured
  deck action, the six-axis class group, the flux-sector menu with
  realized transport consistency, and the unique realized-compatible
  polarization. The four quotient choices stay indistinguishable on local
  tensors alone; 4d instanton/theta and laboratory receipts stay open.

### Open receipts

- *(weaker group-level branch only)* **A5-COMMON-ACTION** + **W5-NONCENTRAL**:
  one group-level action shared by ports and gauge reconstruction, plus one
  source-derived repair composition with nonzero projected `W_5` commutator.
  The source-bound route above is stronger.
- **A5-FAMILY-ATTACHMENT** (family corollary only): prove the chiral family
  fiber's local face-corner phase is the restriction of a global `3` or `3'`
  action, independent of port labels, worker IDs, chart choices, and refinement
  presentation.
- **A5-FAMILY-DESCENT**: derive how the selection symmetry is broken, hidden,
  or forgotten before general family Yukawa matrices are admitted.
- **PHYSICAL-GLOBAL-FORM-CONTINUUM (#567 remainder)**: the deck/loop class,
  sector menu, polarization, and screen-scope Spin attachment are delivered
  at finite source scope; the remaining receipts are four-dimensional
  instanton sectors with topological action normalization (before any theta
  periodicity statement), monopole dynamics, and laboratory flux
  measurement (#569 lane).
- **EXTERIOR-PACKAGE-SELECTION**: the exhaustive 1024-subset anomaly scan
  proves the unordered conjugate pair is the unique nonempty chiral
  anomaly-free selection at exterior-module scope
  (`Lean/Screen/ExteriorSelection.lean`). The MGFC-grade exclusion of light
  sectors beyond the exterior module — direct sums, inert doublets,
  vectorlike matter, neutral singlets — is owned by #609.
- **PORT-WEAK-INTERTWINER** and **PORT-LOAD-TRACE**: identify a physical
  four-dimensional screen invariant with the four weak-doublet copies and
  prove that its normalized additive load is the physical `4P` readout.

## Conditional coupling relation

If the physical quadratic gauge kinetic operator on the port module is an
adjacency polynomial `K = f(A)`, color coherence (the color triplet and quintet
sharing one kinetic coefficient) forces, at degree two,

```text
k1 = 3*k2 - 2*k3
```

independent of which inequivalent triplet carries color (both assignments
verified). At degree one it collapses to bare coupling unification. A physical
discriminator requires a source-derived polynomial degree, kinetic
normalization, carrier scale, and complete threshold/RG map. The relation is
not a forward test.

## Formal status

`Lean/Screen/Compact12.lean` formalizes the abstract
`u(3) + so(3)` matrix commutator, Lie laws, dimension, and a matrix
noncentrality witness. `Z6Exact.lean` formalizes the six-axis lattice quotient,
`S2DesignSignature.lean` formalizes the `11/25` arithmetic, and
`UnitSplit12.lean` proves only that twelve positive integer weights summing to
twelve are all one. `A5OPH.lean` (merged 2026-07-23 from the external
cross-audit lane, rebuilt under the project toolchain) formalizes the finite
core of the compact-Lie trichotomy — triviality of every `A5`-action on at
most four objects, the unique partitions `11 = 3+8` and `12 = 3+3+3+3` over
the compact-simple dimension list, the excluded semisimple dimensions
`1, 2, 4, 5, 7`, the characteristic-centre step — plus the paper's
noncentral-quintet witness, the gluing-class quotient with both invariance
clauses, the `A5 ⊄ SU(2)` unique-involution obstruction, and the absence of
`ℤ/6` in `A5`. `A5CharacterField.lean` (#605, closed) proves the
Galois-stability half of the `Q(sqrt 5)` rationality lemma via the doubled
character table over `ℤ[√5]`, with the torus/cocharacter step a declared
hypothesis. `A5SixAxes.lean` (#604, closed) lists the sixty elements of the
six-axis `PSL(2, F5)` action, kernel-checks distinctness, closure,
2-transitivity, and the sharp stabilizer-coset fiber counts, proves the
five-dimensional summand irreducible over `ℚ` by stabilizer-coset averaging,
and closes the dimension-six branch unconditionally: `1 ⊕ 5` has no
three-dimensional invariant subspace and no `3 + 3` invariant splitting.
`A5PortModule.lean` (#604, closed) carries the matching S5 centre receipt on
the twelve-port module: the sixty port rotations (row-for-row identical to
`A5PortAction.perms`) fix exactly the constant line, so a centrally trivial
submodule — the centre under an inner action — has dimension at most
one. `A5Commutant.lean` (#568, closed) proves the four-orbit structure of
ordered port pairs and the exact four-dimensional orbital commutant of the
port action, with entry invariance equivalent to commutation;
`TraceBalancedKernel.lean` (#568, closed) checks the six-element cyclic
kernel of the trace-balanced cover on central parameters with bijective
`U(1)` coordinate; `TrichotomyCases.lean` (#568, closed) assembles the
centre-versus-semisimple enumeration with every hypothesis in the theorem
signature; `UnitSplit12.lean` carries removability witnesses for each
unit-splitting premise. The independent cross-check
`independent_lane_check.py` (this directory) rebuilds the rotation group
from golden-ratio coordinates alone and reproduces the orbit, commutant,
fiber-count, spectrum, kernel, and enumeration data (22 checks, all
passing 2026-07-23). The issue-#565 artifact supplies an exact executable finite
check of the declared integer domain, strict readback cost, antipode, `A5`
action, frame, refinement, relabelling, and countermodels. The issue-#625
artifact records the reduced-carrier half-count candidate and the exact local
A3 expansion:
the Kullback-Leibler second-order Taylor coefficient is `6 I`, while the
Hessian, equivalently the Fisher matrix, is `12 I`. Neither infinitesimal
object supplies the exact discrete physical cost. A full operational,
refinement-natural three-axiom lift of the half-count candidate is open. The issue-#566
artifact supplies the exact executable current lift, closure, innerness, and
moduli proofs on the declared response branch. The issue-#314 artifact
supplies the exact executable spin lift, matter transport, CAR/Fock
selection, chirality, anomaly, kernel-emission, and declared-tower descent
proofs on top of that packet. A Lean port of those
complete finite packets remains available work; the Python receipts are not
being relabelled as Lean theorems. `Phi`, `Theta`, the trace-balanced group,
and physical descent retain their separate support boundaries (the
compact-Lie trichotomy's finite steps are now Lean-checked in `A5OPH.lean`,
with the former paper steps #604/#605 discharged in `A5SixAxes.lean`,
`A5PortModule.lean`, and `A5CharacterField.lean`, leaving the Lie-theoretic
classical inputs — compact-simple classification, reductive decomposition,
exponential surjectivity, torus/cocharacter step — as the declared
remainder); only the conditional current algebra and the
conditional matter lift are certified here, while source binding of the
inherited upstream response premises is tracked in #599.

## Novelty boundary

Finite-group flavor models use triplet representations, and Standard-Model
global quotients and line sectors are established subjects. This certificate
makes no priority claim. Its scoped result is the full gauge-adjoint module on
the icosahedral vertices plus the conditional compact-Lie classification.
