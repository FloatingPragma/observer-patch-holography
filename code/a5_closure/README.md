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
python3 claim_boundary_certificates.py # exact Q0 fixtures and physical-boundary controls
python3 test_audit.py                    # regression suite
```

Requires Python 3.11+ and SymPy. The suite exits 0.

## What is certified (exact; no floating-point fit, no measured number)

| Object | Statement |
|---|---|
| **Echosahedral counting and incidence selector (#565, #625, #628)** | On the declared federation-of-twelve-port-echosahedra branch, twelve primitive central atoms of trace `1/12`, an additional integer atom-counting grammar on the total-12 fiber, and the normalized central-readback Hilbert-Schmidt cost give the unique all-one split with exact quadratic gap `2`. The complete #625 diagram classification fixes additive counting up to one positive unit scale, so the half-unit readback is a units rescaling rather than a half-event implementation. The #628 patch machine generates integer event counts, applies register-level conservative repairs, and derives its move-count settling cost from those dynamics. Every move lowers the load-square Lyapunov function by the exact amount `2(d-1)`; the seam quadratic and the A3 Hessian remain separate comparison objects. Oriented edge/face incidence and the refinement lineage independently give the unique graph-distance-three antipode, `Aut+ = A5` by a faithful conjugation action on five Klein-four subgroups, six axes, and the exact rank-three Gram frame `G^2=4G`. |
| **A1/A2 Lie-type theorem and conditional current witness (#566, #599)** | The twelve-port module is `1+3+3'+5`. A faithful compact commutator-closed response tangent with endogenous proper-`A5` holonomy has centre dimension at most one; its one fixed line excludes the centreless `su(2)^4` branch, so compact classification forces `u(1)+su(2)+su(3)`. The source artifact independently derives `R=-J`, its four sector signs, the oriented frame, and carrier persistence. It does not select a matrix current or reconstruct a bracket. The declared charged-double-triplet fixture realizes the forced Lie type exactly as `u(3)+so(3)`. Its physical source gate is false until ordered current tomography and same-current overlap holonomy are produced. Twenty-four negative controls fail closed. Lean bridge: `Lean/Screen/A2HolonomyBridge.lean`; executable theorem and fixture: `port_current_inner_certificate.py`. |
| **Conditional super-Tannakian matter lift (#314)** | Given the declared charged-double-triplet current fixture, the finite packet checks the non-split binary-icosahedral lift, the faithful action on `V=C (+) W`, the rank-15 exterior projector `Lambda^2 V + Lambda^4 V`, the Standard Model hypercharge multiset, chirality, anomaly traces, Witten parity, three invariant Yukawa lines, and the common central kernel. The exhaustive 1024-subset scan leaves the unordered conjugate rank-15 pair as the unique nonempty chiral anomaly-free selection inside the declared exterior menu. The finite Spin artifact fixes odd-Weyl typing inside this construction. It does not source-select the current representation, the physical matter action, one conjugate representative, a scalar sector, or a continuum Spin attachment. The conditional algebraic gate passes and the physical matter-source gate fails closed. Twenty typed controls pass. |
| **Conditional axis-centre descent (#567)** | On the declared matter table, exhaustive central enumeration gives the diagonal common kernel generated by `(omega_3 I_3, -I_2, e^{i pi/3})` and the maximal effective image `(SU(3) x SU(2) x U(1))/Z6`. Local tensors also descend through the cover and the `Z2` and `Z3` quotients, so they do not select a global form. The six-axis calculation has Smith residue six only after diagonal and zero-sum coefficient relations are declared. The exact generator intertwiner and line-polarization calculations are conditional on that relation system. A complete source character category, a same-source loop-to-kernel theorem, source-selected matter, and laboratory attachment remain open. The physical global-form gate fails closed. Lean: `Z6Descent.lean`, `ExteriorSelection.lean`; executable receipt schema v4. |
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
compact Lie type follows from the complete faithful A1 response and the A2
same-response internal-holonomy premise. `port_current_inner_certificate.py`
checks one exact charged-double-triplet matrix witness with four signed
coefficients and an oriented frame. The semantic response artifact derives
the response involution, sector signs, frame, and persistence maps. It does
not reconstruct the matrix generators, their bracket, or same-current
overlap holonomy, so the physical source gate remains false.

`super_tannakian_matter_lift_certificate.py` verifies the exact non-split
SU(2) double cover, the unordered rank-15 pair selected by the exhaustive
1024-subset anomaly scan, its parity grading, the typing fixed by the
measured centre and section obstruction, the BLOCK-DETERMINANT-BALANCE
charge pair, scalar and channel compatibility, and the emitted action
kernel inside declared current, exterior, and scalar fixtures. Source
selection of the matrix current and physical matter action, scalar
existence and economy, extra-light-sector completeness, continuum Spin
attachment, and laboratory matter identification remain open.
The #567 descent certificate
(`axis_center_descent_certificate.py`) computes the common kernel on every
declared tensor and the character/cocharacter arithmetic of its maximal
effective quotient. The six-axis class group, flux-sector menu, line
polarization, and refinement checks are exact inside a declared diagonal
and zero-sum relation system. The carrier does not select those relations
or identify its order-six loop class with the complete current and matter
kernel. The physical global-form gate therefore remains false, alongside
the open four-dimensional and laboratory attachments.

`claim_boundary_certificates.py` records five exact-small controls at this
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
The complete #625 diagram classification fixes additive counting up to one
positive unit scale. Its half-unit member is the same counting mechanism in a
different unit convention, and the inequivalent `6 I + A` comparison is
excluded from the operational cost cone. The #628 mechanism supplies an
operational event-counting and repair-cost realization while retaining the
positive readback scale as a unit convention. The results do not prove that
arbitrary OPH carriers must have this type.

- The **conditional PORT-CURRENT-INNER algebraic construction** is verified by
  `port_current_inner_certificate.py`: given the charged-double-triplet model
  and four signed coefficients, the map is full-rank, compact, closed,
  A5-equivariant, inner, and natural along the declared algebraic tower maps.
  These representation data are branch premises, not physical measurements.

- The **conditional super-Tannakian matter lift** (issue #314) is verified
  by `super_tannakian_matter_lift_certificate.py`: given the pinned
  conditional #566 packet, the exact non-split SU(2) double cover of the
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

- **PORT-CURRENT-INNER** verifies the declared charged-double-triplet
  fixture. The semantic response artifact determines `R=-J`, its relative
  signs, the oriented frame, and carrier persistence. Ordered source
  tomography, the bracket, and same-current overlap holonomy are open.
- **PORT-SPIN-LIFT** rides on the conditional current packet. The non-split
  binary icosahedral cover and its unique involution are exact on the finite
  Spin fixture. Physical matter and continuum Spin attachment are open.
- **BLOCK-DETERMINANT-BALANCE** is derived inside the matter packet:
  anomaly freedom of the realized package forces the balance line,
  primitive integrality and the measured orientation select
  `(-1/3, 1/2)`.
- **AXIS-CENTER-DESCENT** computes the common kernel and central descent
  congruence of the maximal effective image on the declared matter table.
  The full character lattice also carries the nonabelian highest weights.
  The six-axis group and polarization use a declared coefficient-relation
  system. Source completeness and the loop-to-kernel identity are open, as
  are four-dimensional instanton/theta and laboratory receipts.

### Open receipts

- *(weaker group-level branch only)* **A5-COMMON-ACTION** + **W5-NONCENTRAL**:
  one group-level action shared by ports and gauge reconstruction, plus one
  source-derived repair composition with nonzero projected `W_5` commutator.
  The A1/A2 inner-action theorem above is stronger.
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
artifact classifies complete additive readbacks as positive multiples of atom
counting and records the exact local A3 expansion: the Kullback-Leibler
second-order Taylor coefficient is `6 I`, while the Hessian, equivalently the
Fisher matrix, is `12 I`. Neither infinitesimal object supplies the exact
discrete physical cost. The issue-#628 artifact supplies the finite
record-repair mechanism and its derived settling cost. The issue-#566
artifact supplies the exact executable current lift, closure, innerness, and
moduli proofs for the declared charged-double-triplet fixture. The issue-#314 artifact
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
remainder); only the conditional current algebra and the conditional matter
lift are certified here. Source reconstruction of the current and matter
actions remains open.

## Novelty boundary

Finite-group flavor models use triplet representations, and Standard-Model
global quotients and line sectors are established subjects. This certificate
makes no priority claim. Its scoped result is the full gauge-adjoint module on
the icosahedral vertices plus the conditional compact-Lie classification.
