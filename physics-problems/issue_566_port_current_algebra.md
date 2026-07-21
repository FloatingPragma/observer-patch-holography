# Issue #566: the physical port-current algebra on the echosahedral response branch

## Result

This artifact closes the **PORT-CURRENT-INNER** receipt on the declared
echosahedral response branch:

> **Theorem (physical port-current algebra).** On a certified twelve-port
> echosahedral carrier lineage with declared reversible response data, there
> is an injective source-derived port-to-generator map
> `K_r : P_12,r -> u(H_r)` with twelve-dimensional image on a faithful
> charged response space `H_r`. Its image is a compact-type skew-adjoint
> commutator-closed Lie algebra with one-dimensional central `u(1)` kernel
> and adjoint rank eleven. The map intertwines the derived icosahedral
> action, the induced `A5` action on the image lies in `Int(g)`, and the
> construction is natural along the declared refinement tower. The
> equivariant lifts form exactly a four-parameter family: four independent
> `A5`-equivariant response band scales and one common odd-response sign are
> the open source data.

No Standard Model current, particle assignment, product-adjoint count,
measured coupling, or gauge target is an input. Every proof decision is
exact arithmetic in `Q(sqrt5)`; no floating point occurs.

The exact executable receipt is produced by:

```bash
python3 code/a5_closure/port_current_inner_certificate.py all
python3 -m unittest discover -s code/a5_closure/tests -v
```

---

## 1. Source packet

Fix a refinement stage `r` on one certified carrier lineage. The source
packet is

\[
\mathsf R_r=\bigl(\mathsf E_r,\ \Lambda,\ \varepsilon,\ \mathsf T\bigr),
\]

with the following fields.

### R1. Certified carrier packet

`E_r` is the issue-#565 source packet: twelve primitive central port atoms of
trace `1/12`, oriented icosahedral edge/face incidence, integer defect
readback, and orientation-preserving refinement lineage. Its closure packet
(`code/a5_closure/echosahedral_selector_certificate.py`) supplies, as derived
objects: the unit lines, the unique fixed-point-free antipode `iota`, the
faithful proper action `Aut+ = A5` with all sixty port permutations, the
exact rank-three Gram frame `G^2 = 4G`, and refinement/relabeling
naturality. The #566 manifest references this carrier manifest by path and
SHA-256 pin; the pin is enforced fail-closed.

### R2. Reversible response typing

The manifest types two disjoint operation families:

* **Reversible response automorphisms**: the separately measured proper
  incidence automorphism family of the carrier. These are declared
  reversible (the family is closed under composition and inverse; the
  verifier rechecks group closure through the #565 enumeration) and are the
  only declared current source.
* **Irreversible strict-descent repairs**: coordination-defect repair
  records. These are typed irreversible, carry `defines_currents: false`
  row by row, and are excluded from the current construction. A repair
  declared as a current source is rejected fail-closed
  (`REPAIR_RESPONSE_CONFLATION`).

### R3. Response scales and the odd sign

`Lambda = (lambda_1, lambda_5, lambda_3, lambda_3')` are four exact rational
response band scales, one per isotypic band of the port module
`P_12 = 1 + 5 + 3 + 3'`, and `varepsilon in {+1, -1}` is one common
odd-response sign. These are the declared source response data. The theorem
proves (part g) that this is exactly the freedom an equivariant lift has, so
nothing else is smuggled in.

### R4. Source firewall

The manifest schema rejects downstream target fields: Standard Model,
product-adjoint, electroweak, hypercharge, Higgs, particle-mass,
measured-coupling, and gauge-target data. The construction is a function
only of R1-R3.

---

## 2. Derived objects

### D1. Oriented frame realization in exact coordinates

By the #565 theorem the port Gram frame is realized by twelve unit vectors,
uniquely up to `O(3)`, and the oriented face class reduces the equivalence
to `SO(3)`. The verifier realizes it with the twelve **unnormalized**
icosahedron vertices — the cyclic permutations of `(0, ±1, ±phi)` with
`phi = (1+sqrt5)/2` — whose entries lie in `Q(sqrt5)` and whose common
squared norm is `N = 2 + phi = phi*sqrt5`. The inner products are exactly
`±N` (same axis) and `±phi` (distance one/two), reproducing the #565
distance dictionary.

An assignment `psi : ports -> vertices` is admissible when it is a distance
isometry mapping the declared coherently oriented faces to positively
oriented vertex triples (`det > 0`, decided exactly). The verifier counts
exactly `120` distance isometries and exactly `60` orientation-matched
assignments: one proper orbit. Exact band data is verified equal for an
alternative assignment, so the derived invariants are
realization-independent.

Write `v_p` for the assigned signed vertex of port `p` (so
`v_{iota(p)} = -v_p`), pick one representative port per antipodal axis, and
let `u_1, ..., u_6` be the six axis vectors.

### D2. Band split of a port field

A real port field `f in P_12` splits under the antipode into even axis
coordinates `b_i = (f_p + f_{iota p})/2` and odd axis coordinates
`d_i = (f_p - f_{iota p})/2`. With `c = (1/6) sum_i b_i` and
`b^0 = b - c 1`, the four bands are the constant even line (`1`), the
sum-zero even band (`5`), and the odd frame/kernel bands (`3`, `3'`) below.

### D3. Frame map and the Galois kernel intertwiner

Let `U = [u_1 ... u_6]` and define

\[
Ud=\sum_i d_i u_i=\tfrac12\sum_p f_p v_p,
\qquad
\sigma(U)d=\sum_i d_i\,\sigma(u_i),
\]

where `sigma` is the Galois automorphism `sqrt5 -> -sqrt5` of `Q(sqrt5)`
applied entrywise. `U` annihilates the kernel band `W = ker U` and
`sigma(U)` annihilates the frame band (Lemma 566.4), so the two maps
separate the two odd bands without any projector choice.

### D4. Charged response space and derived implementers

\[
H_r=\mathbb C^3_E\oplus\mathbb C^3_W .
\]

For every response automorphism `g` the verifier solves exactly for the
unique matrix `R_g` with `R_g v_p = v_{g(p)}` for all twelve ports and
verifies `R_g^T R_g = I`, `det R_g = 1`, entries in `Q(sqrt5)`. The
implementer on `H_r` is

\[
\Pi(g)=\operatorname{diag}\bigl(R_g,\ \sigma(R_g)\bigr).
\]

### D5. The port-to-generator map

With `hat` the cross-product generator map (`hat(x)y = x × y`),

\[
K(f)=\operatorname{diag}\Bigl(
\lambda_3\,\widehat{Ud}
\;+\;i\bigl(\lambda_1 c\,I_3+\lambda_5\,\Phi_0(b^0)\bigr),\ \ \
\varepsilon\lambda_{3'}\,\widehat{\sigma(U)d}
\Bigr),
\qquad
\Phi_0(b^0)=\sum_i b^0_i\,u_iu_i^{\mathsf T}.
\]

`Phi_0(b^0)` is symmetric traceless; the even part of `K` is
`i`-times-Hermitian, the odd part is real skew, so `K(f)` is skew-adjoint.

---

## 3. Main theorem

### Theorem 566.1: physical port-current algebra

Let `(R_r)` be a source packet satisfying R1-R4 with all four scales nonzero.
Then, at every stage `r`:

#### (a) Source-defined operators, domain, inner product, pairing, refinement

The domain is the space of real port fields on the twelve primitive central
atoms; the operators are `K(f)` above, built only from the derived frame,
the Galois intertwiner, and the declared scales; the inner product is the
standard Hermitian pairing on `H_r`; the response pairing is the
Hilbert-Schmidt pullback `B(f, f') = -Re tr(K(f)K(f'))`; the refinement
maps are the declared carrier tower maps. Nothing else enters.

#### (b) Injectivity and twelve-dimensional image

`K` is real-linear and injective with `dim_R K(P_12) = 12`. The image is

\[
\mathfrak g=K(P_{12})=\mathfrak u(3)\oplus\mathfrak{so}(3)
=\mathfrak u(1)\oplus\mathfrak{su}(3)\oplus\mathfrak{su}(2)
\]

in the block presentation (even block `u(3)`, dimension 9; kernel block
`so(3)`, dimension 3).

#### (c) Compact-type skew-adjoint closure

Every `K(f)` is skew-adjoint; the image is closed under the commutator with
structure constants in `Q(sqrt5)` (all 66 basis brackets solved exactly in
the span). The pullback form `B` is symmetric, `A5`-invariant, and positive
definite, and decomposes band-scalar with the exact Hilbert-Schmidt band
coefficients

\[
\beta_1=\tfrac14,\qquad
\beta_5=3+\sqrt5,\qquad
\beta_3=5+\sqrt5,\qquad
\beta_{3'}=5-\sqrt5
\]

at unit scales; the two odd coefficients are Galois conjugates. An invariant
positive-definite form makes `g` compact type.

#### (d) Central `u(1)` kernel and adjoint rank eleven

The center of `g` is one-dimensional and is exactly the image of the
constant even port line (`K(1) = i lambda_1 I_3 (+) 0`). The derived
algebra has dimension eleven, so the adjoint representation has rank eleven
with the central `u(1)` as its kernel. The adjoint action is **not**
required to have rank twelve, and does not.

#### (e) Faithful charged response space

`g` acts faithfully on `H_r`, and the central generator acts with charge
`i lambda_1` on the even sector and `0` on the kernel sector, so the
response space is charged: the full twelve dimensions are visible on `H_r`
even though the adjoint sees only eleven.

#### (f) Icosahedral intertwiner, innerness, refinement naturality

For all sixty response automorphisms and all twelve basis fields,

\[
K(g\cdot f)=\Pi(g)\,K(f)\,\Pi(g)^{*},
\]

`Pi` is a faithful homomorphism (all `3600` products checked), and every
`Pi(g)` is the exponential of an element of `g` (sixty exact rotation
normal-form witnesses), so the induced `A5` action on `g` lies in
`Int(g)`. Every declared refinement map is intertwined by `K`, and the
tower cocycle holds by the #565 receipt.

#### (g) Response moduli

The space of `A5`-equivariant real-linear maps `P_12 -> g` has dimension
exactly four. Hence the four band scales together with the common
odd-response sign are exactly the open source data; no fifth parameter and
no per-axis freedom exists.

#### (h) Physical-current gate

The abelian-record model and every rank-deficient scale choice fail the
gate with typed error codes (Section 6).

---

## 4. Proof

### Lemma 566.2: exact frame realization

The twelve coordinate vertices have the icosahedral distance dictionary, so
the #565 isometry enumerator produces all `120` distance isometries onto
them; central symmetry of the vertex set splits them `60/60` by the exact
face determinant sign, and the coherent manifest orientation selects the
`60` proper assignments, which form one `A5` orbit by composition. Each
assignment realizes the canonical Gram frame because inner products depend
only on graph distance.

### Lemma 566.3: derived implementers

For a response automorphism `g`, solving `R [v_i v_j v_k] = [v_{gi} v_{gj}
v_{gk}]` on one spanning triple gives a candidate matrix over `Q(sqrt5)`;
the verifier then checks `R v_p = v_{g(p)}` on **all twelve** ports,
orthogonality, and `det R = 1`. Uniqueness holds because the vertices span.
The assignment `g -> R_g` is a faithful homomorphism (checked exactly on
all `60 × 60` products), and `sigma` preserves orthogonality, determinant,
and products, so `g -> sigma(R_g)` is also a proper faithful homomorphism.

### Lemma 566.4: the Galois kernel intertwiner

The odd coordinate action `rho(g)` of a response automorphism is a signed
permutation with rational entries, so `sigma` fixes it. Applying `sigma`
entrywise to the intertwining identity `U rho(g) = R_g U` gives

\[
\sigma(U)\,\rho(g)=\sigma(R_g)\,\sigma(U).
\]

Thus `sigma(U)` intertwines the odd module with the Galois-twisted vector
action. Its kernel is a submodule; the odd module is `3 (+) 3'` with the
frame band of type `3` and, since a nonzero equivariant map between
non-isomorphic real irreducibles is impossible (Schur), `ker sigma(U)` is
the frame band and `sigma(U)` restricts to an isomorphism of the kernel
band onto `R^3`. Machine check: the six-column matrix `(Ud, sigma(U)d)`
over the odd basis has exact rank six, so the two odd maps jointly separate
and jointly surject.

### Lemma 566.5: injectivity, image, and closure

On the even part, `b -> lambda_1 c I + lambda_5 Phi_0(b^0)` is injective
onto scalar-plus-traceless symmetric matrices (`Phi` is an isomorphism onto
`Sym(3)` since the six projector outer products are linearly independent —
their Gram has the exact simple eigenvalues from `(u_i . u_j)^2 in {N^2,
phi^2}`), so the even image is `i Sym(3)`, dimension six. On the odd part,
Lemma 566.4 gives all pairs of real skew matrices, dimension six. Hence
`rank K = 12` and

\[
K(P_{12})=i\,\mathrm{Sym}(3)\oplus\mathfrak{so}(3)_E\ \oplus\
\mathfrak{so}(3)_W=\mathfrak u(3)\oplus\mathfrak{so}(3),
\]

which is closed under commutators as a block-diagonal sum of matrix Lie
algebras. The verifier does not take this identification on faith: it
checks that the kernel block is exactly real, that the block projections of
the twelve generators have exact ranks `(9, 3)` (so the even image is all
of `u(3)` and the kernel image all of `so(3)`), solves all 66 basis
brackets in the span exactly, and records the structure constants. The
derived-algebra block projections have exact ranks `(8, 3)`; commutators
are traceless, so the eight-dimensional even derived block is exactly
`su(3)` and the kernel derived block is `so(3)`, which pins the type
`u(1) (+) su(3) (+) su(2)` by machine check rather than by construction.

### Lemma 566.6: compact type and the band coefficients

`B(f, f') = -Re tr(K(f)K(f'))` is `Ad`-invariant on skew-adjoint matrices
and `A5`-invariant by (f). Exact symmetric elimination produces twelve
positive pivots, so `B` is positive definite; a Lie algebra of skew-adjoint
operators with an invariant positive-definite form is compact type. Since
the four bands are inequivalent irreducibles, `B` is band-scalar (Schur);
the verifier reconstructs `B = sum beta_band P_band` exactly. At unit
scales the frame-band value on an axis basis vector is
`-tr(hat(u)^2) = 2|u|^2 = 2N = 5 + sqrt5`, and the kernel band is its
Galois conjugate `5 - sqrt5`; the unit and quintet coefficients are `1/4`
and `3 + sqrt5`.

### Lemma 566.7: center, derived algebra, adjoint rank

The centralizer system `[K(x), K_j] = 0` reduces through the exact
structure constants to a linear system whose null space has dimension one
and is spanned by the constant even field; `K(1) = i lambda_1 I_3 (+) 0`
is central. The span of all 66 brackets has exact rank eleven
(`= su(3) (+) so(3)`). Hence `dim Z(g) = 1`, `dim [g, g] = 11`, and the
adjoint kernel is the central `u(1)`: adjoint rank eleven.

### Lemma 566.8: covariance

All four band components of `K` are built from the equivariant expressions
`(1/12) sum_p f_p`, `(1/2) sum_p f_p v_p v_p^T`, `(1/2) sum_p f_p v_p`, and
its Galois twist. For example
`(1/2) sum_p (g.f)_p v_p = (1/2) sum_q f_q v_{g(q)} = R_g (1/2) sum_q f_q v_q`,
and `hat(R x) = R hat(x) R^T` for proper `R`. The verifier checks the
identity `K(g.f) = Pi(g) K(f) Pi(g)^*` exactly for all sixty automorphisms
and twelve basis fields rather than trusting the derivation.

### Lemma 566.9: innerness

The odd bands supply every block pair `(hat(x), hat(y))` inside `g`
(Lemma 566.4, scales nonzero). For each `g`, each block of `Pi(g)` is a
proper rotation `R` over `Q(sqrt5)`; the verifier computes its exact axis
`n` (`ker(R - I)` one-dimensional for `R != I`), its exact cosine
`cos theta = (tr R - 1)/2`, and verifies the Rodrigues normal form

\[
R=I+\frac{s}{|n|}\,\widehat{n}+\frac{1-\cos\theta}{|n|^2}\,\widehat{n}^{\,2},
\qquad s^2=1-\cos^2\theta,
\]

entirely inside `Q(sqrt5)` (the proportionality scalar `t = s/|n|` is
solved entrywise and `t^2 |n|^2 = 1 - cos^2 theta` is checked exactly).
Hence `R = exp(theta hat(n/|n|))` with skew generator along the certified
axis, so `Pi(g) = exp(X_g)` with `X_g = diag(theta_E hat(n_E), theta_W
hat(n_W))` an element of `g`. Therefore `Ad(Pi(g))|_g = e^{ad X_g}` lies in
`Int(g)`. For order-five elements the two block cosines are the Galois pair
`(-1 ± sqrt5)/4`, recorded per element.

### Lemma 566.10: Schur block rigidity

The even-block operator module of `u(3)` restricts to `A5` as
`1 + 5 + 3` (character `chi_3(g)^2`), which contains no `3'`. Exact
character arithmetic gives

\[
\frac1{60}\sum_{g}\chi_{3'}(g)\,\chi_3(g)^2=0 ,
\]

so no equivariant map can send the kernel band into the even block: the
block allocation of the four bands is forced, not chosen.

### Lemma 566.11: refinement naturality

Every declared tower map is an orientation-preserving incidence
automorphism (#565 receipt), so Lemma 566.8 applies verbatim:
`K(c.f) = Pi(c) K(f) Pi(c)^*` is checked exactly per declared map, and the
tower cocycle is inherited from the carrier receipt. The lines, bands,
current algebra, and gate verdicts therefore descend to one cofinal
refinement class.

### Lemma 566.12: the response moduli are four-dimensional

`Ad Pi(g) K(e_i) = K(rho(g) e_i)` by covariance, so the target module `g`
is equivariantly isomorphic to the port permutation module `P_12`, and

\[
\dim\operatorname{Hom}_{A_5}(P_{12},\mathfrak g)
=\frac1{60}\sum_g|\operatorname{Fix}(g)|^2
=\frac{240}{60}=4 ,
\]

the exact Burnside rank of the icosahedral port action (distance classes
`0, 1, 2, 3`). All four irreducible bands have real endomorphism field
`R`, so the equivariant lifts are exactly the four band rescalings; a sign
choice on the odd bands is the residual discrete datum, declared as the
common odd-response sign. This proves that R3 is the complete open source
data.

### Completion of Theorem 566.1

Parts (a)-(g) are Lemmas 566.2-566.12 combined with the exact executable
checks; part (h) is Section 6. Every derived object is a deterministic
function of R1-R3 and the firewall R4 rejects downstream target data. ∎

---

## 5. Exact certificate values

| Object | Exact value |
|---|---:|
| distance isometries / orientation-matched | `120 / 60` |
| image real dimension | `12` |
| block dimensions (even, kernel), verified | `9, 3` |
| derived block dimensions (even, kernel), verified | `8, 3` |
| compact Lie type (machine-identified) | `u(3) (+) so(3) = u(1) (+) su(3) (+) su(2)` |
| derived dimension / center dimension | `11 / 1` |
| adjoint rank | `11` |
| central line | constant even port field |
| band coefficients `(1, 5, 3, 3')` | `1/4, 3+sqrt5, 5+sqrt5, 5-sqrt5` |
| covariance checks | `60 × 12` |
| implementer homomorphism pairs | `3600` |
| innerness witnesses | `60` (orders `1, 2, 3, 5`) |
| order-five block cosines | `(-1 - sqrt5)/4` and `(-1 + sqrt5)/4`, Galois-paired |
| kernel band multiplicity in even block | `0` |
| equivariant lift moduli | `4` |
| refinement maps intertwined | all declared (`3` in the reference tower) |

The receipt also records all 66 structure constants over `Q(sqrt5)`, the
twelve positive elimination pivots, the sixty rotation normal forms, a
twelve-step `derivation_chain`, `factor_origins` for every numeric
constant, `branch_scope`, `acceptance_criteria_status`, a
`dependency_acyclicity_note`, the `verifier_command`, and every
negative-control outcome.

---

## 6. Physical-current gate and negative controls

The gate requires: source-defined packet; injective twelve-dimensional
image; skew-adjointness; commutator closure; compact type; center exactly
one-dimensional on the constant line; derived dimension eleven; charged
faithfulness; covariance; implementer homomorphism; innerness witnesses;
refinement naturality; repair/response distinctness; firewall. Each failure
carries a stable typed code.

### N1. Abelian record model

`K(f) = i diag(f)` on `C^12` with permutation implementers matches the
coefficient module (twelve-dimensional, equivariant, skew-adjoint,
commutator-closed) and **fails the gate**: derived dimension `0`, center
dimension `12` (`CENTER_NOT_ONE_DIMENSIONAL`). Its induced `A5` action is
also non-inner, since `Int` of an abelian algebra is trivial while the
action permutes the records nontrivially. This is the exact sense in which
coefficient classification does not produce physical nonabelian currents.

### N2-N3. Rank-deficient models

Kernel-band scale `0` drops the image to dimension `9`; unit-band scale `0`
drops it to `11`. Both are rejected as `IMAGE_RANK_DEFICIENT`. The second
control also witnesses the charge rigidity: an equivariant lift either
charges the response space through the scalar central character or
degenerates.

### N4-N5. Broken equivariance

A per-axis rescaling of the quintet response and a non-common odd-axis
sign each break `K(g.f) = Pi(g)K(f)Pi(g)^*` (`COVARIANCE_BROKEN`),
witnessing that the four scales must be band-uniform and the odd sign
common.

### N6. Symmetric record pairing

Dropping the `i` on the even response produces symmetric non-skew
operators (`SKEW_ADJOINTNESS_BROKEN`).

### N7. Repair conflation

An irreversible strict-descent repair declared as a current source is
rejected (`REPAIR_RESPONSE_CONFLATION`), keeping repairs and reversible
response automorphisms distinct.

### N8. Firewall

A measured-coupling target injected into the source manifest is rejected
(`FORBIDDEN_DEPENDENCY`).

---

## 7. Formalization-ready finite model

### JSON source schema

```text
schema
carrier_manifest_path, carrier_manifest_sha256
response_model
response_band_scales.{unit_band,quintet_band,frame_band,kernel_band}
odd_response_sign
reversible_response_automorphisms.{source,reversible,defines_currents}
strict_descent_repairs.{irreversible,defines_currents,ledger[*]}
```

No group table, rotation matrix, current algebra, structure constant, or
target representation appears in the input.

### Deterministic derivation order

```text
validate source firewall and repair/response typing
load and hash-pin the certified carrier packet
re-derive antipode, A5 action, refinement tower (issue-#565 verifier)
enumerate distance isometries onto the exact vertex model; keep the 60 proper ones
build K on the four bands; check skew-adjointness
derive implementers; check covariance (720) and homomorphism (3600)
check rank 12; solve all 66 brackets; derive center, derived algebra, adjoint rank
check positive-definite invariant pullback and exact band coefficients
check assignment independence of the band coefficients
verify 60 rotation normal forms (innerness), refinement naturality
check Burnside moduli = 4 and Schur block rigidity = 0
emit receipt and negative controls
```

### Lean-facing signatures

```lean
structure PortCurrentPacket where
  carrier : EchosahedralCarrier
  scales : Fin 4 → ℚ
  oddSign : ℤˣ

def K (P : PortCurrentPacket) : (P.carrier.Port → ℝ) →ₗ[ℝ] u3 × so3 := ...

theorem K_injective ... : Function.Injective (K P) := ...
theorem K_image_closed ... : IsLieSubalgebra (range (K P)) := ...
theorem K_center_rank ... :
  finrank ℝ (center (range (K P))) = 1 ∧
  finrank ℝ (derivedSeries (range (K P)) 1) = 11 := ...
theorem K_covariant ... : K P (g • f) = Ad (Π g) (K P f) := ...
theorem A5_action_inner ... : ∀ g, Ad (Π g) ∈ innerAutomorphisms := ...
theorem equivariant_moduli ... :
  finrank ℝ (P12 →ₗ[A5] u3 × so3) = 4 := ...
```

`Lean/ObserverPatchHolography/Screen/Compact12.lean` already formalizes the
abstract `u(3) (+) so(3)` bracket, Lie laws, dimension twelve, and a
noncentrality witness; the Python receipt supplies the exact rotation
tables, structure constants, and normal forms needed to instantiate the
remaining signatures without search inside Lean. The receipt is not being
relabeled as a Lean theorem.

---

## 8. Acceptance matrix

| Issue acceptance item | Discharge |
|---|---|
| operators, domain, inner product, response pairing, refinement maps source-defined | R1-R4, D1-D5, receipt `source_definedness` |
| closure, compactness, rank, faithfulness, icosahedral intertwiner proved | Lemmas 566.5-566.9, exact structure constants, positive pivots, 720 covariance checks |
| abelian-record and rank-deficient models fail the physical-current gate | N1-N3 with typed codes, stored countermodel bundle |
| coefficient classification distinguished from physical Yang-Mills current realization | N1 separating witness + receipt `classification_vs_realization`; the trichotomy admits `u(1)^12`, the gate rejects it |
| no measured coupling, particle assignment, or Standard Model current as input | schema firewall R4, fail-closed control N8 |
| central `u(1)` kernel, adjoint rank eleven (not twelve) | Lemma 566.7, receipt `closure` |
| repairs distinct from reversible response automorphisms | R2, control N7 |

---

## 9. Claim boundary

This theorem closes **PORT-CURRENT-INNER** on the declared echosahedral
response branch: given the certified carrier and the declared reversible
response data, the physical current lift exists, is essentially unique (four
scales and one sign), and passes the full gate.

It does not derive the response scales, the odd sign, or the existence of
the charged response sectors from raw OPH consensus dynamics; those are
declared source response data on this branch, exactly as the carrier
architecture is declared in issue #565. It does not close block determinant
balance, `PORT-SPIN-LIFT`, the physical `Z6` deck/line descent
(`AXIS-CENTER-DESCENT`), exterior package or matter selection, family
attachment, or any continuum Yang-Mills, coupling, mass, or measured-number
statement. The Lie type `u(1) (+) su(3) (+) su(2)` is a statement about the
constructed current algebra on this branch, not an identification with the
physical Standard Model gauge group; that identification requires the open
receipts above.
