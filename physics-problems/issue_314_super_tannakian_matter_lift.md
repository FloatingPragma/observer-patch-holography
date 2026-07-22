# Issue #314: conditional super-Tannakian matter lift on declared matter contracts

## Result

This artifact proves the exact algebraic **SUPER-TANNAKIAN-MATTER-LIFT**
construction — including the exact algebraic **PORT-SPIN-LIFT** target —
conditional on the pinned conditional #566 current packet and the declared
matter-lift contracts. It does not close either receipt physically:

> **Theorem (conditional super-Tannakian matter lift).** On a certified
> twelve-port echosahedral carrier lineage, given the hash-pinned
> conditional #566 current packet (a declared charged-double-triplet
> response representation with four signed coefficients), the declared
> trace-balanced exterior matter contract, the fermionic statistics
> contract, and the Spin/odd-Weyl category typing, the source packet
> derives: an exact non-split SU(2) double cover of the sixty proper
> implementers (binary icosahedral, unique involution — the algebraic
> PORT-SPIN-LIFT target); a faithful skew-adjoint Lie-algebra
> homomorphism of the twelve-dimensional current algebra onto the matter
> carrier `V = C (+) W`; a 32-state auxiliary CAR/Fock space with
> derived fermionic parity, super tensor structure, and conjugation
> through the exactly invariant top line; a derived equivariant
> projector of exact rank fifteen selecting the exterior package
> `M1 = Lambda^2 V (+) Lambda^4 V` from that space; exact chirality,
> vanishing realized perturbative anomaly traces, even Witten parity, and
> exactly one invariant line per declared Yukawa channel; and the common
> action kernel on the simply connected cover `R x SU(3) x SU(2)` —
> infinite cyclic with generator `(zeta_6, omega, -1)`, whose sixth power
> is the unit deck translation (one full central turn, not the identity
> on the cover), of residual order exactly six modulo the pure deck
> translations — emitted as data for the downstream global-form descent
> without forming either the central compactification or the quotient.
> The construction
> descends naturally along the declared algebraic tower maps, and the
> realized packet witnesses the declared MAR class nonempty without
> promoting uniqueness.

The algebraic properties are machine-checked. The conditional current
algebra is strictly upstream through the hash-pinned #566 packet, whose
physical source gate is recorded open; the matter-lift contracts — the
trace-balanced block charge pair, the one-scalar choice, the Yukawa
channel list, the statistics and category typing, the kernel emission
contract, and the MAR class declaration — enter as explicit typed branch
premises, not as measurements or physically source-bound data. The
receipt therefore records a passing conditional algebraic gate and an
open physical source-realization gate (`issue_closure_condition`).

No family attachment, scalar potential, pole mass, measured coupling, or
global-form choice is an input. Every proof decision is exact arithmetic
in `Q(sqrt5)`; no floating point occurs.

The exact executable receipt is produced by:

```bash
python3 code/a5_closure/super_tannakian_matter_lift_certificate.py all
python3 -m unittest discover -s code/a5_closure/tests -v
```

---

## 1. Source packet

The source packet is

\[
\mathsf M=\bigl(\mathsf K_{566},\ (y_C,y_W),\ \mathsf S,\ \mathsf C,\ \mathsf E,\ \mathsf R\bigr),
\]

with the following fields.

### M1. Pinned upstream current packet

`K_566` is the issue-#566 source packet, referenced by path and SHA-256
pin together with its stored receipt hash. The verifier requires the
stored receipt to be a #566 port-current receipt certifying exactly the
pinned manifest with a passed conditional algebraic gate and a recorded
open physical source gate (this packet inherits that conditionality),
then **rebuilds** the twelve current generators, the sixty proper
rotations, and all 66 structure constants from the pinned source (nothing
is copied from the receipt). Dependencies are #565 and #566 only.

### M2. Trace-balanced exterior matter contract

`(y_C, y_W)` is an exact rational charge pair on the color and weak
blocks with `3 y_C + 2 y_W = 0` (checked arithmetically from the pair,
`TRACE_BALANCE` otherwise; no redundant declared balance flag is
accepted); the reference values are `(-1/3, 1/2)`. The contract further
declares the one-scalar choice (the weak block itself) and the invariant
Yukawa channel list `(Q S u_c)`, `(Q Sbar d_c)`, `(L Sbar e_c)`. An empty
channel list is rejected (`GAUSS_DATA_EMPTY`). This contract is a typed
branch premise of the conditional exterior lane, not a measurement: the
balance itself is **not** derived here (BLOCK-DETERMINANT-BALANCE stays
open).

### M3. Fermionic statistics contract

`S` types the matter statistics `fermionic_odd` and keeps matter distinct
from the bosonic coefficient records of the screen lane. The typing is
checked against a derived fact: the matter-building operators of the
auxiliary CAR algebra anticommute exactly, so a bosonic-even typing fails
closed (`STATISTICS_TYPING`).

### M4. Spin/odd-Weyl category contract

`C` types the category `spin_odd_weyl_super` with `double_cover: true`,
the selection rule `parity_even_minus_derived_invariants`, and the
realization `operator_projector`. A production manifest accepts only
these values; weaker same-reduct typings are admitted only through the
control path (`allow_control_contracts`) and fail there against derived
facts: `vec` cannot carry the derived grading and cocycle
(`VEC_TYPING`), a declared split lift contradicts the derived
unique-involution witness (`SPIN_LIFT_SPLIT`), and
`representation_arithmetic` is rejected because representation arithmetic
alone is not realization (`REALIZATION_NOT_OPERATOR`).

### M5. Kernel emission contract

`E` requires the common action kernel to be emitted as data
(`emit: true`) and forbids assuming the final global quotient
(`assume_global_quotient: false`, `KERNEL_EMISSION_CONTRACT` otherwise).

### M6. MAR class declaration

`R` declares the class `one_generation_one_scalar_chiral_anomaly_free`
with `promote_uniqueness: false`. Promoting uniqueness inside this packet
is rejected (`MAR_UNIQUENESS_PROMOTION`); the packet only discharges the
nonemptiness precondition.

### M7. Source firewall

The #565 firewall plus a matter extension reject family-attachment,
scalar-potential, pole-mass, measured-coupling, and global-form-choice
tokens (`FORBIDDEN_DEPENDENCY`).

---

## 2. Derived objects

### D1. The upstream current algebra

From the pinned #566 packet the verifier re-derives the oriented frame,
the generators `K_p` on `C^3_E (+) C^3_W`, the sixty rotations `R_g`, and
the exact structure constants. The kernel block transforms in the
Galois-conjugate rotation family `sigma(R_g)`, exactly as in the #566
covariance.

### D2. PORT-SPIN-LIFT

For each proper implementer the verifier solves exactly, in `Q(sqrt5)`,
for a special-unitary `U_g` with

\[
U_g\,(\sigma\!\cdot\!v)\,U_g^{\dagger}=\sigma\!\cdot\!\bigl(\sigma(R_g)\,v\bigr),
\]

using the exact half-angle data: `u^2 = (1 + cos theta)/2` and
`c^2 = (1 - u^2)/|n|^2` both have exact square roots in `Q(sqrt5)` for
every icosahedral rotation (orders 1, 2, 3, 5). The 120-element set
`{±U_g}` is verified closed under multiplication with **exactly one
involution** (`-1`), the binary icosahedral order profile
`{1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24}`, order-four lifts of all
involutions, and irrational order-five spinor characters.

### D3. Matter carrier and transport

\[
V=C\oplus W=\mathbb C^3\oplus\mathbb C^2 ,
\qquad
T(K_p)=\Bigl(A_p\ \text{on}\ C\Bigr)\oplus
\Bigl(-\tfrac i2\,(w_p\!\cdot\!\sigma)+i\,\tfrac{\tau_p}{3}\,\tfrac{y_W}{y_C}\,I_2\ \text{on}\ W\Bigr),
\]

where `A_p` is the source even block (kept exactly, including its trace
part `i tau_p`), `w_p` is the source kernel-block axis, and the weak
central term is the declared trace-balanced redistribution. `T` is
verified skew-adjoint, a Lie-algebra homomorphism on all 66 brackets, and
faithful (rank 12 on the 5-dimensional carrier).

### D4. Auxiliary CAR/Fock space

`F = Lambda^* V` (32 states) with exact creation operators `a_i^dagger`.
The verifier checks all 50 CAR relations, vacuum cyclicity (rank 32), the
parity grading `(-1)^N` (currents parity-even, creations parity-odd), and
the second-quantization derivation identity
`[dGamma(X), a^dagger(v)] = a^dagger(Xv)` on all 60 generator-mode pairs;
with vacuum annihilation and cyclicity this pins `dGamma` as the unique
Lie-homomorphic extension. The wedge product with Koszul signs is the
super tensor structure; the CAR-derivation identity is its exact Leibniz
rule.

### D5. Derived equivariant projector

The joint kernel of the twelve second-quantized currents is computed
exactly: complex dimension `2`, spanned by the Fock vacuum and the top
line `Lambda^5 V`. The top line is exactly invariant **because of the
declared trace balance** (`tr T(K_p)|_V = 0`). The matter projector is

\[
\Pi=P_{\mathrm{even}}-P_{\mathrm{vac}} ,
\]

parity projector minus the derived invariant line in the even sector:
diagonal in the subset basis, idempotent, self-adjoint, commuting with all
twelve currents and with parity, of exact rank **15**.

---

## 3. Main theorem

### Theorem 314.1: conditional super-Tannakian matter lift

Let `M` be a source packet satisfying M1-M7. Then:

#### (a) The algebraic PORT-SPIN-LIFT target, non-split

The sixty proper implementers lift to SU(2) with lift group of order 120
and a unique involution. A split extension `A5 x Z2` would carry 31
involutions; hence the cover is genuinely non-split (binary icosahedral),
and no two-dimensional realization of the kernel-block action exists
without it. Order-five spinor characters are irrational, so the spinor
sector admits no signed register-relabeling realization.

#### (b) Faithful current action on the matter tensors

`T` is a faithful skew-adjoint Lie-algebra homomorphism on `V` (rank 12),
and the second-quantized action remains faithful on the selected
fifteen-state module (rank 12). The 720 conjugation transports
`pi_V(g) T(K_p) pi_V(g)^* = T(K_{g(p)})` hold exactly.

#### (c) Derived super structure

Fermionic parity is `(-1)^N` from the CAR grading; the tensor structure is
the wedge with Koszul braiding (exact CAR anticommutation and Leibniz
identity); conjugation is the wedge pairing into the top line, which is
nondegenerate and exactly invariant; chirality is derived, not declared
(part e). Nothing is imposed by hand: each structure is an exact
consequence of the CAR algebra over the derived carrier.

#### (d) The fifteen-state module from the equivariant projector

`Pi` selects `M1 = Lambda^2 V (+) Lambda^4 V` (15 states) from the full
32-state auxiliary CAR/Fock space. The realized charge spectrum is exactly
`{1/6: 6, -2/3: 3, 1: 1, 1/3: 3, -1/2: 2}` with derived integrality
normalization `6`; the five blocks `Q, u_c, e_c, d_c, L` have dimensions
`6, 3, 1, 3, 2` and scalar commutants (multiplicity-free), and no derived
invariant line survives inside the module.

#### (e) Chirality

Any intertwiner `S` with the dual module maps each charge eigenspace into
the negated eigenvalue; the exact spectra are disjoint, hence
`Hom(M1, M1*) = 0`. The conjugate module is realized concretely in the
opposite parity sector through the top-line pairing.

#### (f) Realized anomalies and Witten parity

On the realized operators, the traces `tr(Y)`, `tr(Y^3)`,
`tr(Y T_a T_b)` (all 64 su(3) pairs), `tr(Y S_i S_j)` (all 9 su(2)
pairs), and the 120 symmetrized `su(3)` d-symbol traces vanish exactly.
The realized weak-doublet count is `4` — even Witten parity (the finite
mod-2 surrogate the exterior packet lists).

#### (g) Nonzero invariant Gauss sector

Each declared Yukawa channel carries an invariant subspace of exact
complex dimension `1` (three lines in total); the forbidden control
channel `(Q S d_c)` carries dimension `0`.

#### (h) Emitted common action kernel

On the simply connected cover `R x SU(3) x SU(2)` — with the central
factor the non-compact `R`, not `U(1)` — the set of `(central turn,
su(3) center, su(2) center)` triples acting trivially on every realized
weight — matter states, carrier modes, and the declared scalar — is
computed by exact congruence enumeration with a derived denominator
bound and periodicity argument. The kernel is **infinite cyclic**,
generated by `g = (zeta_6, omega, -1)`; its sixth power is
`g^6 = (1, 1, 1)`, the unit deck translation (one full central turn),
which acts trivially on every integral weight but is **not** the
identity on the cover. The pure deck translations `<g^6>` form the
kernel's intersection with the central `R` factor, and the residual
modulo them has order exactly **6**. Closure, torsion-freeness over the
translations, and the deck relation are verified element by element.
Charge, triality, and duality are additive over tensor factors, so
triviality extends to every realized matter tensor. The kernel —
generator, deck relation, and residual — is **emitted as data**; neither
the compactification of the central `R` nor the global quotient is
formed (AXIS-CENTER-DESCENT consumes the data downstream, and the
emitted generator and relation determine the kernel image in every
candidate quotient).

#### (i) Declared-tower descent

Every declared algebraic tower map lifts through the frame rotation and
the spin lift to an exterior-power implementer that commutes with `Pi`
and intertwines all twelve currents on both the carrier and the Fock
realization. Physical refinement intertwining is not source-bound here.

#### (j) MAR nonemptiness without uniqueness

The realized packet — fifteen multiplicity-free states, one declared
scalar, exact chirality, vanishing realized anomaly traces — witnesses the
declared MAR class nonempty. Uniqueness is not promoted; the contract
rejects promotion inside this packet.

#### (k) Conditional algebraic gate

Vec, split-sVec, opposite-Weyl, bosonic-statistics, truncated-selection,
full-even-module, empty-Gauss, assumed-quotient, kernel-killing,
representation-arithmetic-only, charge-dead, unbalanced,
uniqueness-promoting, and firewall countermodels fail closed with typed
codes (Section 6).

---

## 4. Proof

### Lemma 314.2: exact spinor lifts exist over Q(sqrt5)

For a proper rotation of order `1, 2, 3, 5` with exact cosine
`1, -1, -1/2, (-1 ± sqrt5)/4`, the half-angle squares are
`1, 0, 1/4, (3 ± sqrt5)/8`, and `(3 ± sqrt5)/8 = ((1 ± sqrt5)/4)^2` is a
square in `Q(sqrt5)`. The scaled sine coefficient satisfies
`c^2 |n|^2 = 1 - u^2` with `|n|^2` the exact squared axis norm; because
the true icosian lift has entries in `Q(sqrt5)`, `c` itself lies in
`Q(sqrt5)` and the exact square-root solver finds it. The adjoint
transport check `U (sigma.v) U^dagger = sigma.(Rv)` on all three Pauli
directions pins `U` up to the global sign, and either sign choice is
verified explicitly. ∎

### Lemma 314.3: the lift is a non-split double cover

The 120-element set `{±U_g}` is closed under multiplication (all 14400
products verified in the set), so it is a group `G~` with a 2-to-1
projection onto the sixty rotations. If the extension split, a subgroup
isomorphic to `A5` would embed in `G~` and carry the fifteen involutions
of `A5`; the verifier counts exactly **one** involution in `G~` (namely
`-1`), so no such subgroup exists and the extension is non-split. The
order profile matches the binary icosahedral group exactly, and every lift
of an order-two rotation has order four. ∎

### Lemma 314.4: the matter transport is a faithful homomorphism

The even block is kept exactly, so the su(3)-plus-central part of every
bracket is inherited from the source algebra. The kernel-block map
`hat(w) -> -(i/2)(w · sigma)` satisfies
`[-(i/2)(u·sigma), -(i/2)(v·sigma)] = -(i/2)((u × v)·sigma)`, matching
`[hat(u), hat(v)] = hat(u × v)`; the central terms commute with
everything, and the trace of a commutator vanishes, so central
coefficients match on both sides. All 66 brackets are checked exactly
against the source structure constants. Faithfulness: the traceless even
parts span su(3) (dimension 8, exact pivots), the lifted kernel parts span
su(2) (dimension 3), and the central lane acts nontrivially since
`(y_C, y_W) != 0`; the flattened images have exact rank 12. ∎

### Lemma 314.5: dGamma is the Lie-homomorphic extension

`dGamma(X)` annihilates the vacuum and satisfies the CAR-derivation
identity `[dGamma(X), a^dagger(v)] = a^dagger(Xv)` (60 exact checks). Two
operators that annihilate the vacuum and have equal commutators with every
creation operator agree on the vacuum-cyclic space (rank-32 cyclicity is
verified), and `[dGamma(X), dGamma(Y)]` and `dGamma([X, Y])` both have
these properties by the Jacobi identity, so `dGamma` preserves brackets. ∎

### Lemma 314.6: the derived invariant sector and the projector

The joint kernel of the twelve `dGamma(T(K_p))` has exact complex
dimension 2. The vacuum is invariant since `dGamma` kills it; the top line
is invariant because `dGamma(X)` acts on `Lambda^5 V` by `tr(X|_V)` and
the declared trace balance makes every transported generator traceless on
`V`. The exact nullspace computation shows there is nothing else. `Pi` is
then the parity projector minus the derived invariant line in the even
sector; both ingredients are canonical (the CAR grading and the derived
invariant sector plus the Fock inner product), so `Pi` is derived rather
than declared,
and its equivariance, self-adjointness, idempotence, and rank 15 are
checked exactly. ∎

### Lemma 314.7: realized package identification

The charge operator is diagonal on the subset basis with additive weights,
so the exact spectrum on `im Pi` is
`{1/6: 6, -2/3: 3, 1: 1, 1/3: 3, -1/2: 2}`; together with triality and
duality this splits the module into the five exterior blocks with the
listed dimensions. Each block has scalar commutant (exact complex
commutant dimension 1 from the joint commutant nullspace), so each block
is irreducible and the module is multiplicity-free. ∎

### Lemma 314.8: chirality from disjoint spectra

An intertwiner `S : M1 -> M1*` satisfies `S dGamma(Y) = -dGamma(Y)^T S`,
so `S` maps the `Y = lambda` eigenspace into the `-lambda` eigenspace of
the dual spectrum. The exact spectra `{1/6, -2/3, 1, 1/3, -1/2}` and their
negation are disjoint, hence `S = 0`. Because the representation is
skew-adjoint (unitary at group level), the conjugate module coincides with
the dual, so this is full chirality. ∎

### Lemma 314.9: conjugation through the invariant top line

The wedge pairing `B(u, v) = coefficient of e_top in u ∧ v` between the
even selection and its parity complement is a signed permutation pairing
(nondegenerate, rank 15). Invariance
`B(dGamma(X)u, v) + B(u, dGamma(X)v) = 0` holds exactly for all twelve
generators because the top line is exactly invariant (Lemma 314.6). Hence
the parity-complement selection realizes the conjugate module inside the
same Fock space, using `Lambda^k V ≅ (Lambda^{5-k} V)* ⊗ Lambda^5 V` with
trivialized top line. ∎

### Lemma 314.10: realized anomaly identities

The module is a direct sum of irreducibles, so each realized trace is the
sum of the per-block contributions of the exterior package; the
machine checks compute the traces directly on the 15-dimensional realized
operators and find exact zeros for `tr(Y)`, `tr(Y^3)`, all
`tr(Y T_a T_b)`, all `tr(Y S_i S_j)`, and all symmetrized
`tr(T_a {T_b, T_c})`. The doublet count comes from the exact su(2)
Casimir diagonal: eight states with nonzero diagonal, hence four doublets,
which is even. ∎

### Lemma 314.11: the emitted kernel is infinite cyclic with residual Z6

Any kernel element `(r, a, b)` of `R x Z3 x Z2` must satisfy
`r q + a t/3 + b d/2 in Z` for every realized weight `(q, t, d)` in the
derived integral normalization (`N = 6`). Integer combinations of the
conditions give `r in (1/(6g))Z` for `g` the gcd of the nonzero integral
charges, and membership depends on the numerator `k = 6g·r` only through
`k mod 6g` (adding `6g` to `k` shifts each phase by the integer `q`), so
the exact enumeration over one fundamental window of residues is
exhaustive for the full non-compact `R` coordinate.

The unit deck translation `(1, 1, 1)` — one full central turn with
trivial center components — always lies in the kernel, since `r = 1`
times an integral charge is an integer; on the cover it is **not** the
identity, so the kernel is infinite. The enumeration finds exactly six
residues, one per numerator residue class, with no torsion residue
`(0, a, b) != (0, 0, 0)`; hence the kernel projects injectively to
`(1/6)Z` and is infinite cyclic, generated by `g = (1/6, omega, -1)`
with `g^6 = (1, 1, 1)` (verified exactly, without reducing the `R`
coordinate modulo full turns) and residual order six modulo `<g^6>`.
Membership checks, e.g. on `Q` (`q = 1, t = 1, d = 1`):
`1/6 + 1/3 + 1/2 = 1 in Z`; on `u_c` (`q = -4, t = 2, d = 0`):
`-2/3 + 2/3 = 0`; on the carrier modes `(q, t, d) = (-2, 1, 0)` and
`(3, 0, 1)`: `-1/3 + 1/3 = 0` and `1/2 + 1/2 = 1`. Charge, triality, and
duality are additive over tensor factors, so triviality on the verified
list extends to every realized matter tensor. Closure,
torsion-freeness, and the deck relation are verified element by
element; neither the central compactification `R -> U(1)` nor the
global quotient is ever formed. ∎

### Completion of Theorem 314.1

Parts (a)-(j) are Lemmas 314.2-314.11 combined with the exact executable
checks; part (k) is Section 6. The current algebra is strictly upstream
through the hash pins of M1; the matter-lift contracts M2-M6 enter as
typed branch premises and are enforced fail-closed; the firewall M7
rejects downstream target data. The receipt keeps the conditional
algebraic and physical source-binding gates separate. ∎

---

## 5. Exact certificate values

| Object | Exact value |
|---|---:|
| spin-lift witnesses / lift group order | `60 / 120` |
| involutions in the lift group | `1` (non-split) |
| lift order profile | `{1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24}` |
| order-five spinor characters irrational | `24 / 24` |
| transport bracket checks / covariance transports | `66 / 720` |
| faithful rank on carrier / on matter | `12 / 12` |
| CAR relation checks / vacuum cyclic rank | `50 / 32` |
| second-quantization derivation checks | `60` |
| derived invariant sector | complex dimension `2` (vacuum, top line) |
| projector rank | `15` |
| realized charge spectrum | `{1/6: 6, -2/3: 3, 1: 1, 1/3: 3, -1/2: 2}` |
| integrality normalization | `6` |
| block dimensions (`Q, u_c, e_c, d_c, L`) | `6, 3, 1, 3, 2`, all commutants `1` |
| chirality | spectra disjoint, `Hom(M1, M1*) = 0` |
| realized anomaly traces | all exactly `0` (incl. 120 d-symbol checks) |
| Witten parity | `4` weak doublets, even |
| Yukawa invariant lines | `1, 1, 1`; forbidden control `0` |
| emitted kernel (on the cover `R x SU(3) x SU(2)`) | infinite cyclic, generator `(zeta_6, omega, -1)`, `g^6 =` unit deck translation, residual order `6` |
| refinement maps intertwined | all declared (`3` in the reference tower) |
| MAR class | nonempty witnessed, uniqueness not promoted |
| conditional algebraic gate | `passed: true` |
| physical source-realization gate | `passed: false`; issue closure remains open |

The receipt also records a sixteen-step `derivation_chain`,
`factor_origins` for every numeric constant, `branch_scope`,
`acceptance_criteria_status` (nine rows, with the source-derivation row
honestly `false`), a `dependency_acyclicity_note`, the
`verifier_command`, and every negative-control outcome.

---

## 6. Conditional algebraic gate and negative controls

### N1-N2. Vec and split-sVec same-reduct typings

`vec` typing is rejected against the derived nontrivial grading and
non-split cocycle (`VEC_TYPING`); `svec` with a declared split lift is
rejected against the derived unique-involution witness
(`SPIN_LIFT_SPLIT`). Both fail against derived facts, not by fiat.

### N3. Opposite-Weyl same-reduct selection

The odd-parity selection realizes the conjugate module; every declared
one-scalar Yukawa channel then has charge sum away from zero, so the
Gauss sector is empty (`YUKAWA_CHANNEL_EMPTY`). Together with N1-N2 this
discharges the issue's requirement that the Vec/sVec and opposite-Weyl
same-reduct controls do not pass alongside the reference.

### N4. Bosonic matter statistics

The derived matter-building operators anticommute; bosonic-even typing
fails closed (`STATISTICS_TYPING`).

### N5-N6. Truncated and over-full selections

`lambda2_only` yields ten states with three weak doublets — odd Witten
parity (`WITTEN_PARITY`). `even_including_vacuum` is the full even
Clifford module, which keeps the trivial vacuum line inside matter
(`TRIVIAL_LINE_IN_MATTER`); this is exactly the exterior lane's
"non-vacuum even package, not the full even Clifford module" boundary.

### N7. Empty Gauss data

An empty declared channel list is rejected (`GAUSS_DATA_EMPTY`).

### N8-N9. Kernel contract

Declaring the global quotient assumed violates the emission contract
(`KERNEL_EMISSION_CONTRACT`). An extra integral-charge singlet scalar
collapses the computed kernel to the pure full-turn deck translations;
a trivial residual kernel cannot satisfy the packet (`KERNEL_TRIVIAL`).

### N10. Representation arithmetic only

A packet declaring `representation_arithmetic` realization is rejected
(`REALIZATION_NOT_OPERATOR`): the fifteen-state module must be selected by
the operator projector receipt on the auxiliary CAR/Fock space.

### N11. Charge-dead package

`(y_C, y_W) = (0, 0)` is trace-balanced but annihilates the central lane;
the current action drops to rank 11 and faithfulness fails
(`CURRENT_ACTION_NOT_FAITHFUL`).

### N12. Unbalanced trace charges

`3 y_C + 2 y_W != 0` is rejected at the contract (`TRACE_BALANCE`); the
top line would not be invariant and conjugation would fail downstream.

### N13. MAR uniqueness promotion

`promote_uniqueness: true` is rejected (`MAR_UNIQUENESS_PROMOTION`).

### N14-N15. Firewall

Family-attachment and scalar-potential/pole-mass injections are rejected
(`FORBIDDEN_DEPENDENCY`).

---

## 7. Formalization-ready finite model

### JSON source schema

```text
schema
current_manifest_path, current_manifest_sha256
current_receipt_path, current_receipt_sha256
exterior_matter_contract.{block_trace_charges,one_scalar,yukawa_channels,extra_scalars}
statistics_contract.{matter_statistics,auxiliary_car_modes,distinct_from_bosonic_records}
category_contract.{typing,spin_lift,selection_rule,realization}
kernel_emission_contract.{emit,assume_global_quotient}
mar_class.{declared,promote_uniqueness}
```

No spin lift, Fock operator, projector, anomaly value, kernel element, or
target representation appears in the input.

### Deterministic derivation order

```text
validate contracts and the matter firewall
hash-pin and gate-check the #566 packet; rebuild generators, rotations, structure constants
solve the sixty exact SU(2) lifts; verify closure, unique involution, order profile
build the matter transport T; check 66 brackets, 720 transports, rank 12
build the CAR/Fock space; check 50 CAR relations, cyclicity, parity, 60 derivation identities
compute the derived invariant sector (dimension 2); build Pi = P_even - P_vac; rank 15
identify the five blocks; check commutants, spectrum, integrality normalization
check realized anomaly traces, Witten parity, chirality, conjugation pairing
compute the three Yukawa invariant lines and the forbidden-channel zero
enumerate and emit the common action kernel on the cover (infinite cyclic; residual order six); do not compactify or quotient
check refinement descent on all declared tower maps
emit receipt and negative controls
```

### Lean-facing signatures

```lean
structure MatterLiftPacket where
  current : PortCurrentPacket
  charges : ℚ × ℚ           -- (y_C, y_W) with 3 y_C + 2 y_W = 0
  channels : List YukawaChannel

def spinLift (P : MatterLiftPacket) : A5 → SU2 := ...
theorem spinLift_nonsplit ... :
  ∀ x ∈ liftGroup P, orderOf x = 2 → x = -1 := ...

def T (P : MatterLiftPacket) : currentAlgebra →ₗ⁅ℝ⁆ u (C³ ⊕ C²) := ...
theorem T_faithful ... : Function.Injective (T P) := ...

def Pi (P : MatterLiftPacket) : FockSpace →ₗ[ℂ] FockSpace := ...
theorem Pi_equivariant_rank ... :
  (∀ p, Commute (Pi P) (dGamma (T P (K p)))) ∧ rank (Pi P) = 15 := ...

theorem kernel_emitted ... :
  actionKernel P ≃ ℤ ∧ generator = (ζ₆, ω, -1)
  ∧ generator ^ 6 = unitDeckTranslation
  ∧ actionKernel P ⧸ deckTranslations ≃ ZMod 6 := ...
```

The Python receipt supplies the exact lift matrices, transport images,
projector, spectra, and kernel table needed to instantiate these
signatures without search inside Lean. The receipt is not being relabeled
as a Lean theorem.

---

## 8. Acceptance matrix

| Issue acceptance item | Discharge |
|---|---|
| fermionic parity, spin lift, chirality, conjugation, tensor product source-derived | open physically: each structure is derived exactly on the branch — parity from the CAR grading (D4), the spin lift solved from the pinned rotations with the non-split witness (D2, Lemmas 314.2-314.3), chirality from disjoint exact spectra (Lemma 314.8), conjugation from the wedge pairing into the invariant top line (Lemma 314.9), the tensor structure from the Koszul wedge (D4) — while the upstream response representation and the matter contract remain declared branch premises, not source-bound data |
| physical current algebra acts faithfully at Lie-algebra level on the matter tensors | the conditional #566 current algebra (physical source binding open there) acts with rank 12 on the carrier and rank 12 on the selected fifteen-state module; 66 exact bracket checks; 720 exact conjugation transports (Lemma 314.4) |
| conditional exterior package realized on the cover; listed perturbative anomalies and Witten parity checked | the package is realized as operators on the cover data (SU(3)-valued even implementers, the SU(2) spin lift, the R-central lane); all listed anomaly traces vanish exactly on the realized operators and the Witten parity is even (Lemma 314.10) |
| common action kernel emitted rather than assumed as a Z6 quotient | the kernel is computed by exact congruence enumeration on the genuine cover `R x SU(3) x SU(2)`: infinite cyclic with generator `(zeta_6, omega, -1)`, sixth power the unit deck translation (not the identity on the cover), residual order six modulo the pure deck translations; verified trivial on every realized weight and emitted as data; `assume_global_quotient` is contractually false, and neither the central compactification nor the quotient is ever formed (Lemma 314.11) |
| declared MAR class proved nonempty before uniqueness is promoted | the realized packet is the witness; `promote_uniqueness: true` is rejected fail-closed (M6) |
| family attachment and scalar-potential or pole-mass claims outside | firewall M7 plus controls N14-N15; the claim boundary lists them as not closed |
| Spin/odd-Weyl category nonempty and source-produced; Vec/sVec and opposite-Weyl same-reduct controls do not both pass | the realized category carries the nonzero fifteen-state odd matter object with faithful action, produced from the pinned conditional packet (source production is conditional on its open source binding); N1-N3 fail against derived facts |
| nontrivial faithful current action and nonzero physical invariant sector; empty Gauss data or zero common kernel cannot satisfy | rank-12 faithfulness plus three exact invariant Yukawa lines; `GAUSS_DATA_EMPTY` and `KERNEL_TRIVIAL` controls N7, N9 |
| 15-state module selected from the full auxiliary CAR/Fock space by a source-derived equivariant projector; representation arithmetic is not realization | `Pi = P_even - P_vac` is derived, not declared, on the 32-state space: equivariant, rank 15 (Lemma 314.6); `REALIZATION_NOT_OPERATOR` control N10 |

---

## 9. Claim boundary

This theorem proves the conditional exact algebraic matter lift for the
declared matter contracts over the pinned conditional current packet:
given those premises, the super-Tannakian matter category exists, is
realized by exact operators, passes the conditional algebraic gate, and
its common action kernel is emitted for the downstream global-form
descent. It does not close **SUPER-TANNAKIAN-MATTER-LIFT** or
**PORT-SPIN-LIFT** as physical source-bound receipts;
`issue_closure_condition.met_locally` is therefore `false`.

The matter-lift contracts — the trace-balanced block charge pair, the
one-scalar choice, the Yukawa channel list, the statistics and category
typing, the kernel emission contract, and the MAR class declaration —
enter as typed branch premises, not as measurements. The upstream #566
packet is itself conditional on a declared response representation whose
physical source binding is open, and this packet inherits that
conditionality.

It does not derive the trace balance from source dynamics
(BLOCK-DETERMINANT-BALANCE stays open), does not source-bind the
response representation or the physical refinement maps, does not choose
the global form (AXIS-CENTER-DESCENT consumes the emitted kernel
downstream), does not promote MAR uniqueness, does not attach families
or count them, does not exclude other anomaly-free light sectors (the
MGFC-grade no-extra-sector boundary), and makes no scalar-potential,
pole-mass, measured-coupling, continuum spin-statistics, or
quantum-field-theory claim. The realized package is a statement about
the constructed matter category on this branch, not an identification
with physical particle content; that identification requires the open
receipts above.
