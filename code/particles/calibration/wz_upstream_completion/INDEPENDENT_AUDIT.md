# Independent audit — Physical W/Z upstream mathematical completion

Date: 2026-07-20  
Audited archive: `physical_wz_upstream_mathematical_completion_2026-07-20.zip`  
Archive SHA-256: `e6184ec1a535244c0f805c14e275b085a64c7ed078d780828c038d4f21c86a0d`

## Verdict

> **Post-audit integration note.** The isolated v4.2 copy now applies M1--M8:
> complete QCD gauge fixing, the rank-`n-1` Laurent hypothesis, separate loop
> and Planck symbols, deterministic/stochastic output separation, removal of
> minimality, strengthened conditional reconstruction assumptions, finite-EFT
> field-redefinition restrictions, and residual-specific numerical bounds. It
> also corrects the 9-schema/10-instance count, repairs clock-parent DAG edges,
> makes the template checker unconditionally non-promoting, and adds explicit
> QCD-gauge, Laurent-rank, FJ-orientation, and self-attestation regressions.
> Integrated validation is now 8 symbolic checks and 38 tests. C1--C5 remain
> open as **production** work: the v1 interfaces still do not carry or resolve
> proof-bearing artifacts. The audit below records the defects as found in the
> pristine Pro archive.

The package is a useful, scientifically conservative **draft sufficiency
specification**. Its central strict-pole algebra and its main non-entailment
arguments are sound. It also correctly refuses current OPH-native W/Z
promotion.

It is not yet a complete mathematical completion or an executable scientific
receipt system. In particular:

1. the shipped checker validates fixed templates and trusts producer booleans
   and digest strings; it does not verify any proof-bearing artifact;
2. the schemas are checklist envelopes rather than receipts containing the
   action, matching maps, FJ conversion, BRST residuals, contours, null vectors,
   amplitudes, laws, covariances, or clock calculation;
3. the declared full-SM BRST gauge-fixing action omits the QCD gauge-fixing and
   QCD ghost sector;
4. Theorem 15.1 is false without a rank/nullity hypothesis;
5. the clock formula reuses the already-defined loop factor `h` as though it
   were Planck's constant;
6. the final theorem conflates a stochastic pole law with one deterministic
   pole pair and conflates covariance with a certified enclosure; and
7. “minimal” or “smallest” is not proved and is not true in the literal logical
   sense.

Accordingly, this archive closes **no positive OPH-native W/Z scientific
issue**. It defines a strong implementation target after the corrections below.

## 1. Archive integrity and reproducibility

### What passed

- The archive SHA-256 is the value shown above.
- `unzip -t` passed for every member.
- Every member lies below the single declared top-level directory; no absolute
  path, `..` traversal, or symbolic link was present.
- In a fresh extraction, the manifest verified all 35 listed files, byte
  counts, aggregate byte count, and SHA-256 digests: **35/35**.
- The frozen gauge grid contains 45 distinct IDs and 45 distinct parameter
  tuples. Its internal canonical hash verifies.
- The receipt DAG is acyclic and its internal canonical hash verifies.
- `python3 run_all.py` completed with:

  - 7 symbolic checks passing;
  - the template checker returning
    `THEOREM_STACK_DEFINED__SIMULATION_RECEIPTS_OPEN__NO_OPH_NATIVE_POLE_PROMOTION`;
  - 34 pytest tests passing.

### Reproducibility defect

`run_all.py` rewrites `outputs/validation.log`, which is itself listed in the
package manifest. The rewritten log includes absolute interpreter and checkout
paths. Therefore the officially documented validation command makes the
manifest stale immediately after a successful run.

Required correction: either make validation output go to an unmanifested build
directory, canonicalize path-independent output and regenerate the manifest,
or publish separate immutable `SOURCE_MANIFEST` and post-run
`VALIDATION_MANIFEST` objects.

### Schema-count correction

There are **nine distinct JSON Schema documents**, not ten. The checker has ten
validation roles because W and Z instantiate the same physical-pole schema
twice. Documentation and status JSON should say:

```text
distinct_schema_documents = 9
receipt_instances = 10
```

## 2. Results that are mathematically sound

The following results can be retained after notation and scope cleanup.

### 2.1 Structural non-identifiability

Theorem 1 is correct: fixed gauge group, representations, generation count, and
one Higgs doublet leave an open family

\[
p=(g_s,g,g',m^2,\lambda,Y_u,Y_d,Y_e),
\]

and already

\[
m_W^2=\frac{g^2}{4}\frac{-m^2}{\lambda},\qquad
m_Z^2=\frac{g^2+g'^2}{4}\frac{-m^2}{\lambda}
\]

vary inside that family. Structural SM data therefore do not identify a W/Z
pole pair.

Theorems 2 and 2.1 are also correct: an idempotent normalizer does not select a
law, and one density matrix or MaxEnt representation does not uniquely select a
Lorentzian Hamiltonian, clock, correlation hierarchy, or 1PI functional.

### 2.2 Tree electroweak and anomaly arithmetic

For

\[
V=m^2H^\dagger H+\lambda(H^\dagger H)^2,
\qquad m^2<0,\quad \lambda>0,
\]

the positive tree coordinate

\[
v_F=\sqrt{-m^2/\lambda}
\]

gives

\[
w=\frac{g^2v_F^2}{4},\qquad
z=\frac{(g^2+g'^2)v_F^2}{4},\qquad
m_h^2=2\lambda v_F^2.
\]

The local anomaly arithmetic and the even weak-doublet count are correct for
the declared conventional SM census. The theorem should be titled
“perturbative anomalies plus the conventional SU(2) Witten anomaly” unless a
separate global-anomaly/bordism argument is added for the chosen global
quotient.

### 2.3 FJ parameter identity and finite reparametrization

Differentiating

\[
m^2+\lambda v_F^2=0
\]

gives the correct parameter counterterm

\[
\frac{\delta v_F^{\rm par}}{v_F}
=\frac12\left(\frac{\delta m^2}{m^2}
-\frac{\delta\lambda}{\lambda}\right).
\]

For the declared map orientation

\[
p_L=p_F+\hbar\Delta p^{(1)}+O(\hbar^2),
\]

the strict one-loop conversion rule is also correct:

\[
s_{1,F}(p_F)=s_{1,L}(p_F)
+\Delta p^{a(1)}\partial_a s_0(p_F).
\]

This is a conditional coordinate theorem. It is not an OPH-to-FJ conversion
until an actual complete map and its independently recomputed residuals exist.

### 2.4 Strict charged and neutral pole series

With the declared convention

\[
\Gamma_W^T=s-w+h\Pi_{WW}^{(1)}+h^2\Pi_{WW}^{(2)}+O(h^3),
\qquad h=(16\pi^2)^{-1},
\]

the coefficients

\[
s_{W,1}=-\Pi_{WW}^{(1)}(w),
\]

\[
s_{W,2}=-\Pi_{WW}^{(2)}(w)
-s_{W,1}\Pi_{WW}^{(1)\prime}(w)
\]

are correct.

For the neutral matrix, the strict coefficients

\[
s_{Z,1}=-\Pi_{ZZ}^{(1)}(z),
\]

\[
s_{Z,2}=-\Pi_{ZZ}^{(2)}(z)
-s_{Z,1}\Pi_{ZZ}^{(1)\prime}(z)
+\frac{\Pi_{ZA}^{(1)}(z)\Pi_{AZ}^{(1)}(z)}{z}
\]

are also correct. The off-diagonal product is order \(h^2\), so it must not be
inserted in a strict-one-loop root.

The strict square-root expansion is correct, as is the warning that an exact
square root of a loop-truncated pole introduces higher-order kinematic terms.

### 2.5 Conditional Nielsen and Rouché results

Given a correctly renormalized extended Slavnov--Taylor identity and Nielsen
kernels regular at a simple pole,

\[
\partial_\eta\Gamma=\Lambda\Gamma+\Gamma\widetilde\Lambda
\]

implies

\[
\partial_\eta\det\Gamma
=\operatorname{tr}(\Lambda+\widetilde\Lambda)\det\Gamma,
\]

and therefore gauge-parameter independence of a simple determinant zero. This
is standard conditional mathematics; it does not verify a diagram
implementation. The cited Nielsen paper supports that conditional use.

The Rouché theorem statement is correct. A production receipt still needs an
actual contour, analytic reference, ball inequality, and interpolation bound.

### 2.6 Pure-SM gauge beta coefficients

For `gprime = gY`, three generations, and one complex Higgs doublet,

\[
b_{g'}=\frac{41}{6},\qquad b_g=-\frac{19}{6},\qquad b_{g_s}=-7
\]

is correct. Rejecting the MSSM tuple on an otherwise pure-SM interval is also
correct.

## 3. Blocking mathematical corrections

### M1 — QCD is not gauge-fixed

Definition 8 introduces the SU(3) ghost \(c_s^A\), but Definition 9's
gauge-fixing fermion contains only charged electroweak, Z, and photon terms.
There is no \(F_s^A\), \(\bar c_s^A\), \(b_s^A\), or QCD gauge parameter in
\(\Psi\).

The displayed action is BRST invariant, but it is only **partially** gauge
fixed. It cannot define a complete full-SM propagator/rule/counterterm engine.
This omission does not change a bare strict one-loop electroweak W/Z
self-energy, where no internal gluon occurs, but it blocks the document's
stronger claims of a complete SM BRST and renormalization receipt.

Required correction: add, for example,

\[
F_s^A=\partial^\mu G_\mu^A,
\qquad
\Psi_s=\int d^4x\,\bar c_s^A
\left(F_s^A+\frac{\xi_s}{2}b_s^A\right),
\]

include \((\xi_s,\chi_s)\) if QCD gauge independence is in scope, generate the
QCD ghost action from \(s\Psi_s\), and bind it to the action/rule hashes.

The covariant-derivative sign and the derived relation
\(e=gs_W=g'c_W\) must also be frozen explicitly.

### M2 — The mixed-propagator Laurent theorem lacks a necessary hypothesis

Theorem 15.1 assumes null vectors and

\[
\ell^\dagger\Gamma'(s_p)r\ne0
\]

but does not assume a one-dimensional kernel or rank \(n-1\). It is false as
written. A counterexample is

\[
\Gamma(s)=\operatorname{diag}(s,s),\qquad s_p=0.
\]

For \(r=\ell=e_1\), the displayed derivative is nonzero, but
\(\Gamma^{-1}=I/s\), not the claimed rank-one principal part.

Required correction: assume `rank Gamma(s_p) = n-1`, equivalently one-dimensional
left and right kernels and a simple determinant zero, before applying the
rank-one residue formula.

### M3 — The operational-clock equation uses the wrong `h`

Definition 13 fixes

\[
h=(16\pi^2)^{-1}.
\]

Theorem 16.3 later writes

\[
E_\star=\frac{h\nu_{\rm clk}}{\varepsilon_{\rm clk}}
\]

without redefining `h`. As written, this is dimensionally wrong. It must use
Planck's constant:

\[
E_\star=\frac{h_{\rm P}\nu_{\rm clk}}
{\varepsilon_{\rm clk}}
=\frac{\hbar_{\rm P}\omega_{\rm clk}}
{\varepsilon_{\rm clk}}.
\]

The receipt must bind the SI/J/eV conversion convention and constants. The
loop-expansion symbol should be renamed, for example, \(\kappa=(16\pi^2)^{-1}\).

### M4 — Deterministic and stochastic conclusions must be separated

The final theorem allows either a deterministic source or a stochastic law,
then concludes that there is “one source-defined physical pair.” That conclusion
is valid only in the deterministic mode.

The stochastic conclusion must instead be a pushforward law

\[
\mu_{WZ}=F_\#\mu_{\rm src}
\]

on \((s_W,s_Z)\), with its joint correlations retained. A covariance matrix is
only a second-moment summary; it is not a support enclosure and cannot be
“composed” into one without support or tail assumptions. Certified deterministic
remainders must be propagated separately by interval/ball maps.

Required replacement:

- deterministic branch: one certified pole pair plus deterministic error
  enclosure;
- stochastic branch: one certified joint pole law, covariance/pseudocovariance
  where defined, declared confidence/tail object, and separate deterministic
  theory enclosure.

Complex covariance must be typed either on the real vector
\((\Re s_W,\Im s_W,\Re s_Z,\Im s_Z)\) or as covariance plus
pseudocovariance; a bare transpose formula is otherwise ambiguous.

### M5 — “Smallest/minimal” is unsupported

The package proves sufficiency of one large augmented branch, not minimality.
For a strict W/Z calculation, a correctly bound renormalized 1PI/current packet
can be sufficient without exposing a full ultraviolet action packet. Even
inside the SM lane, W/Z one-loop two-point functions can be parameterized by
mass eigenvalues and flavor invariants rather than all weak-basis matrices.
Independent duplicate engines are valuable assurance, not a mathematical
necessity for existence.

Required correction: replace “smallest sufficient” and “minimal” with
“one explicit sufficient augmented branch,” unless an irredundancy theorem and
a precise ordering of admissible information structures are supplied.

### M6 — The finite-source reconstruction theorem is not yet sufficient

Theorem 4A invokes reflection positivity and continuum Schwinger functions but
does not state enough Osterwalder--Schrader data to reconstruct a unique
Lorentzian QFT. It needs, in an appropriate gauge-invariant observable algebra,
the full relevant OS conditions: Euclidean covariance, permutation symmetry,
reflection positivity, regularity/temperedness, clustering or vacuum
uniqueness, compatibility of the continuum limit, and a precise reconstruction
and uniqueness statement. Gauge fields and chiral measures need their own
careful formulation.

The existence and differentiability of a global Legendre transform also require
gauge fixing or a quotient/restricted source construction; “convexity
sufficient” merely restates the desired conclusion.

Required correction: either expand this into a genuine conditional OS theorem,
or explicitly demote it to a research-program hypothesis. Keep the formal
strict-one-loop validation lane separate.

### M7 — Field-redefinition equivalence is too broad at finite EFT order

The equivalence class currently admits any local invertible analytic field
redefinition. At finite EFT order this can generate operators outside the mask,
change Jacobians, and move terms between retained and discarded orders.

Required correction: restrict equivalence to order-by-order maps with the
Jacobian/local anomaly treatment, transformed sources, induced operators, and
remainder bound included in the common term mask.

### M8 — The universal residual-radius formula is not a theorem

The recommended condition

\[
\operatorname{rad}(R_{256})
\le 8\sum_k r_k+2^{-128}\max(1,S)
\]

has no derived universal constant or dimensional normalization. Input radii can
belong to quantities with different dimensions and sensitivities.

Required correction: derive a residual-specific interval/Jacobian or Lipschitz
bound in normalized coordinates. A fixed safety factor may be a policy only
after its domain and scaling are declared; it is not proof-bearing by itself.

## 4. Checker and receipt audit

### C1 — The checker is a fixed-template linter, not an aggregate verifier

`validate_all()` always loads the ten files named in the hard-coded `TEMPLATES`
map. There is no command-line or API input for a production bundle. The main
routine deliberately exits if those shipped templates ever promote.

This is useful for checking that examples remain non-promotable. It cannot
validate a simulator run.

### C2 — Scientific evidence is self-attested

The promotion predicate consumes booleans such as:

```text
promotion_ready
independent_checker
complete
same_pole
rouche_passed
linearized_identity_passed
covariance_psd_checked
target_blacklist_passed
```

It never resolves the named artifact, recomputes its digest, or independently
recalculates the stated result. An arbitrary nonzero 64-hex string is accepted
as a hash. The checker compares only selected action and term-mask strings.

The prose requires all receipts to share action, census, scheme, FJ, term-mask,
sheet, and source-root identities. The code checks only:

- six action-hash fields; and
- six term-mask fields.

It neither verifies the digest preimages nor enforces the full seven-key
cross-receipt subject tuple.

Required architecture:

1. every producer emits immutable artifacts and a canonical manifest;
2. the aggregate verifier receives an explicit bundle directory;
3. it canonicalizes and hashes every artifact itself;
4. booleans are derived outputs, never trusted inputs;
5. every verifier binds the exact common subject tuple;
6. evidence is produced by a separately versioned checker hash;
7. caller-supplied claims are retained only under
   `unverified_evidence_claims`; and
8. promotion is the verifier's output, never a producer field.

### C3 — The 34 tests overstate adversarial coverage

At least 17 tests set a field to the value it already has in the fail-closed
template, or assert a reason that was already present before the mutation.
Examples include target ancestry, chart lane, matching uniqueness, beta
derivation, Yukawa completeness, FJ-map completeness, gamma5 restoration, UV
cancellation, bare counterterm generation, diagram completeness, Nielsen,
precision nesting, Z-current equality, clock gap, and clock naturality.

Because tests do not begin from a fully admissible synthetic bundle and do not
compare `before` with `after`, these assertions do not show that the individual
mutation caused failure.

Required test pattern:

```text
valid independent synthetic fixture -> zero promotion reasons
apply exactly one mutation            -> exactly the named new reason
restore it                            -> zero promotion reasons
tamper one bound artifact             -> digest/semantic failure
```

The synthetic fixture must contain no physics target data and cannot be used as
a scientific receipt; it exists only to exercise verifier logic.

### C4 — The schemas omit the proof-bearing payloads

#### Action

- No canonical action AST is present, only a hash.
- `minItems: 18` does not enforce the actual 19-entry declared SM census.
- Field IDs need not be unique; representations, multiplicities, activity, and
  canonical normalization are not semantically checked.
- Operator exclusions are arbitrary strings and do not cover the required
  theta, neutrino, higher-dimensional, scalar, and vectorlike choices.
- No renormalized numerical parameter vector is present.

#### Matching

- Only three gauge beta coefficients are represented; the Yukawa, Higgs, and
  mass beta system is absent except behind an opaque hash.
- Intervals need not be adjacent, nonoverlapping, ordered as a list, or cover
  the requested path.
- Threshold enclosures need not lie in their interval or be ordered.
- Decoupling maps, Jacobians, contribution masks, and remainder vectors are
  hashes without resolved artifacts.
- The special pure-SM beta check can be bypassed by changing an unconstrained
  `eft_id` string.

#### Yukawa

- `active_fermions` can be any 12 strings.
- Unitarity, SVD identities, CKM construction, matrix scale, and open-channel
  accounting are not recomputed.
- Entries are bare floating-point values, not exact/interval values.

#### FJ equivalence

- There is no \(\Delta p\) map, parameter list, field map, counterterm map, or
  derivative payload.
- Residual balls and equality flags are not recomputed.
- Engine IDs and hashes may be identical while independence booleans say they
  are independent.
- Coefficient units and normalization are absent.

#### Renormalization and BRST

- No bare maps, counterterm AST, UV coefficient table, ST residuals, restoration
  linear system, basis rank, nullspace, or normalization conditions are
  present.
- Diagram completeness can be asserted with `diagram_count = 1`.
- The 128/192/256-bit fields contain no balls; only booleans.
- The frozen gauge-grid digest is not bound into the BRST receipt.
- Engine identifiers can be equal while the independence flags remain false.

#### Physical pole

- `s_pole` is a point-valued float pair, not a complex ball.
- No contour geometry, reference determinant, boundary balls, or interpolation
  proof is present.
- No null-vector components, residual balls, derivative ball, or Laurent
  denominator ball is present.
- The physical channel is an arbitrary nonempty string; there is no amplitude,
  current, external-state, or residue artifact/hash.
- `argument_count = 1`, `same_pole`, and `rouche_passed` are declarations.

#### Source law and covariance

- No deterministic source point values are present.
- No stochastic support, outcomes, weights, continuous density, covariance
  matrix, or correlation object is present.
- Jacobians are unresolved hashes with no dimensions or values.
- Nonstochastic errors are names without bounds.

#### Clock

- No Hamiltonian, eigenstates, transition operator, gap derivation, or
  calibration calculation is present.
- The source-parent list accepts an arbitrary single string.
- Planck/SI/eV conversion is not represented.
- The dependency DAG gives `SOURCE_CLOCK_1` no incoming source/action/law/
  no-target parent edges.

### C5 — Required production schemas

The simulator needs resolved, canonical objects for at least:

1. `renormalized_parameter_packet_v1` — exact/ball values, units, scheme,
   scale, basis, and source root;
2. `canonical_action_ast_v1` — fields, operators, coefficients, conventions,
   exclusions, and canonical digest;
3. `eft_interval_map_v1` — complete beta monomials, thresholds, finite maps,
   Jacobians, and vector remainders;
4. `fj_conversion_map_v1` — every parameter/field/counterterm shift and every
   tree derivative;
5. `identity_residual_bundle_v1` — exact or ball-valued UV, Ward, ST, Nielsen,
   and FJ residuals;
6. `diagram_universe_v1` — canonical universe and produced-record set with a
   set-difference proof;
7. `pole_contour_receipt_v1` — contour, sheet path, determinant balls,
   argument/Rouché proof, root ball, derivative ball, and uniqueness;
8. `physical_current_amplitude_v1` — physical external states/currents,
   amplitude implementation, pole equality, and nonzero residue ball;
9. `source_law_v1` and `source_covariance_v1` — actual law and derived moments,
   not flags; and
10. `operational_clock_calculation_v1` — source Hamiltonian, selector, gap ball,
    calibration, constants, ancestry, and unit conversion.

## 5. Exact scientific status

The strongest defensible status after this package is:

```text
structural_non_entailment_theorems = substantially complete
strict_pole_series_theorems = complete for the declared convention
conditional_nielsen_and_rouche_theorems = complete after stated hypotheses
upstream_sufficiency_specification = draft; corrections required
production_receipt_schemas = not complete
production_aggregate_verifier = absent
external_SM_validation_engine = absent
OPH_to_renormalized_action = absent
OPH_to_FJ_vev = absent
source_matching_certificate = absent
source_law_and_covariance = absent
independent_gauge_BRST_receipt = absent
physical_current_WZ_pole_receipts = absent
operational_clock = absent
OPH_native_dimensionless_pole = false
OPH_native_physical_GeV_pole = false
```

The package's existing negative status is safe. Its positive line
`mathematical_sufficiency_stack_defined = true` should be narrowed to
`draft_sufficiency_stack_defined = true` until M1--M8 and C1--C5 are resolved.

## 6. GitHub issue closure implications

This audit made no GitHub mutation.

### Closure-ready independently of this archive

- **#521 — Source-only Higgs pole pipeline:** closure-ready only on its
  published negative branch, using the separate fixed-\(P\), target-free,
  hash-bound two-completion witness in
  `the retired staging workspace/issue-521-negative-closure/`. This upstream archive's broad
  structural non-entailment theorem is supportive context but does not replace
  that issue-specific witness.

### Correctly closed negative result

- **#590 — Force exactly the Standard Model:** keep closed on its accepted
  delimitation/no-go branch. This package strengthens the distinction between
  structural SM type and a numerical renormalized action. It must not be cited
  as a positive action or common electroweak carrier.

### Must remain open

- **#547 — Common screen/electroweak load carrier:** no physical common carrier
  or `v_chart -> v_F` theorem was produced.
- **#545 — Physical pixel readback/root:** no target-clean source root and no
  D10-to-`SM_MSbar_FJ` matching certificate was produced.
- **#334 — Source-derived physical clock:** only a schema and conditional unit
  formula were supplied; the source clock calculation is absent, and the
  formula itself needs the Planck-constant correction.
- **#575 — Lorentzian event manifold:** the archive lists event/frame/cone
  requirements but instantiates none.
- **#574 — Cyclicity and modular intersections:** `MGNS-1` remains a requirement,
  not a receipt.
- **#503, #573, and #592:** no required campaign/source receipts were produced.
- Any positive Physical W/Z, FJ, BRST, matching, covariance, or physical-pole
  issue remains open.

### Reopen/re-scope recommendations unchanged

- **#32:** reopen; no physical threshold/decoupling/scheme-matching artifact
  was produced.
- **#34:** reopen; no target-clean D11 source producer was produced.
- **#525:** reopen or explicitly supersede by #575; the archive documents the
  H3/event-base defect but does not repair it.
- **#308:** reopen or re-scope to a conditional theorem; the independent mixed-
  GNS premise is still absent.

## 7. Acceptance gate for the simulator

The external-SM validation lane may be built before OPH source closure, but it
may be called complete only when all of the following are true:

1. a complete imported SM parameter packet is hash-bound to a canonical full-SM
   action, including QCD gauge fixing;
2. two actually independent rule/self-energy/FJ implementations produce raw
   diagram and coefficient artifacts;
3. a third verifier recomputes digests, diagram set equality, UV cancellation,
   Ward/ST/Nielsen residuals, and the full FJ coordinate map;
4. strict one-loop order masks are enforced, with no accidental resummation or
   two-loop A-Z product leakage;
5. complex balls contain every scalar integral, sum, root, derivative, contour
   inequality, and engine residual;
6. W and Z each have a declared continuation path, one certified simple root,
   and the rank-\(n-1\) Laurent hypothesis;
7. a gauge-invariant charged-current amplitude and a gauge-invariant
   neutral-current amplitude contain the same respective pole with nonzero
   residue balls;
8. all hashes are recomputed from artifacts and the common seven-field subject
   tuple is identical across every receipt;
9. the mutation suite starts from an admissible target-free synthetic verifier
   fixture and proves one-cause/one-failure behavior;
10. source law, covariance, deterministic errors, and clock remain hard-false in
    the external validation lane; and
11. target comparison is read-only post-processing and cannot influence any
    source or validation artifact.

Only after that external lane passes should an OPH-produced action/matching/
law/clock packet be substituted. Doing so isolates simulator correctness from
the still-open OPH source theorem.

## 8. Reference check

The cited primary literature supports the standard conditional ingredients:

- [Gambino and Grassi, Nielsen identities and complex-pole gauge
  independence](https://arxiv.org/abs/hep-ph/9907254);
- [Dittmaier and Rzehak, tadpole/VEV scheme distinctions](https://arxiv.org/abs/2203.07236);
- [Grassi, Hurth, and Steinhauser, practical algebraic
  renormalization](https://arxiv.org/abs/hep-ph/9907426).

Those references do not certify this software, instantiate the OPH source
packet, prove the archive's claimed minimality, or repair the missing
proof-bearing receipt payloads.
