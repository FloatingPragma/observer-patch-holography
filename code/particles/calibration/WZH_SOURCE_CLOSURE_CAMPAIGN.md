# W/Z/H source-closure campaign

## QFT tier status after the July 20 theorem audit

This campaign implements producers, not theorems already discharged on paper.
The conditional QFT-Q1 local-action implication, QFT-Q2-E determinant-line
criterion, QFT-Q2-H Hamiltonian soundness theorem, formal QFT-Q3 BV/ST
restoration implication, strict fixed-parameter W/Z coefficient algebra, and
separate QFT-Q4 OS/resonance implications are now frozen in
`q1_q4_wz_theorem_oracle/`.

The dependency graph is not linear: QFT-Q2 and QFT-Q3 are parallel descendants
of QFT-Q1. Strict W/Z algebra is a QFT-Q3 leaf and has no promotion edge to
QFT-Q4. This campaign remains open until it supplies the source-selected action
and normalization, one exact QFT-Q2 object, regulator-specific QFT-Q3 matching
and restoration evidence, two independent identity engines, dressed physical
current amplitudes, numerical freeze and covariance, source law and clock, and
the independently required QFT-Q4 observable tower and analytic sheet.

## Current result

The repository has exact conditional composition theorems, not a completed
source-only pole prediction. The D10 evaluator proves the algebra after
QT1--QT5 are assumed. The D11 evaluator proves the declared split formulas
after their coefficients are assumed. Neither evaluator supplies the missing
finite carrier characters, and neither mass chart is a BRST-complete complex
pole calculation.

The expected D10 certificate
`code/particles/runs/calibration/d10_ew_quotient_path_certificate.json` is
intentionally absent. No QT or DS gate may be marked closed until the evidence
below is generated from independently frozen carriers.

## What kind of computation this is

This is exact, proof-producing finite enumeration, not Monte Carlo simulation.
The D10 and D11 carriers must be observer-like self-reading patches: bounded
local state, typed ports and boundaries, readback, persistent records,
feedback or repair moves, and a public evidence bundle. A generic graph chosen
to reproduce the desired coefficients is not an OPH derivation.

## Final theorem boundary

The generic implication stack is complete. Model-extension counterfamilies
prove that the present structural axioms do not emit the D10 coefficients, D11
split, absolute scale, or physical poles. A formal weighted path realization is
not evidence by itself: any finite polynomial can be encoded one monomial per
path. The carrier is non-vacuous only when its transitions, admissibility,
quotient, response map, normalized measure, mismatch signs, exhaustive census,
rigidity/MAR gap, and no-target DAG were fixed independently of the desired
output.

The remaining source-emission packets are exactly:

1. `C_clk`: factorized electromagnetic/electron/nuclear/atomic clock packet;
2. `C_10`: independently weighted and rigid D10 carrier;
3. `C_11`: D11 split-character and rigidity carrier;
4. `P_pole`: BRST-complete W, neutral photon/Z, and Higgs two-point kernels.

Once supplied, certified clock refinement, deterministic RG/matching, Rouché
pole isolation, Nielsen gauge independence, pole-displacement enclosures, and
mass/width error propagation require no stochastic simulation. Runtime target
separation still does not prove historical blindness; a prospective frozen
claim and disclosure are separate requirements.

## Phase 0: freeze the carriers

Freeze the pre-repair D10 and pre-split D11 transition systems before running
an enumerator. Their manifests must hash the slots, source rules, transitions,
quotient generators, protected observables, return/depth convention, and
admissible-deformation class.

The carriers must not contain measured W/Z/H values, the desired D10 incidence
table, `1/3`, a claimed projector trace `1/6`, the fibre formulas, the D11
residual polynomials, or their desired output coefficients. A forbidden-payload
scan is a guard, not a proof that the carrier is physically justified.

## D10 certificate: QT1--QT5

1. Enumerate every legal raw path with exact arithmetic. Distinguish edge
   length from return count, start from singleton paths, and record every
   rejection reason.
2. Quotient presentation redundancy. Keep orbit members, stabilizers, rules,
   and canonical representatives. A second canonicalizer should reproduce the
   orbit root, but agreement is regression evidence rather than completeness.
3. Build the exact response matrix without prefiltering to `q_2,q_n`. Prove
   those vectors independent, prove image equality, and emit left-nullspace
   witnesses for forbidden channels.
4. Derive the primitive bidegree and four-slot measure. Transitivity gives
   equal weights only after invariant support and total probability one are
   separately proved.
5. Derive the charged and neutral path characters. Use exactly one declared
   weighting convention: either a raw group average or groupoid-cardinality
   weights with stabilizers. Do not quotient and then divide by the same orbit
   size again.
6. Construct the actual center action. On a representation where the physical
   Standard Model quotient center acts trivially, the projector is the identity
   and its normalized trace is one. A trace of `1/6` requires a separately
   declared D10 center-label/transport representation satisfying
   `dim(V^Z6)/dim(V)=1/6`. Group order six supplies the averaging coefficient,
   not that trace. The source amplitude and subtraction sign are independent
   derived obligations.
7. Derive the charged and neutral signs from an oriented mismatch first
   variation. Construct the fibre vectors, Gram matrix, and residual pairing
   before applying the generic convex minimizer lemma.
8. Complete a consistent finite deformation census. Either prove no
   output-changing quotient direction survives, or use an independently
   justified target-free selector with a positive winner gap. A new MAR vector
   cannot be introduced only because it selects the desired law.

## D11 certificate: DS1--DS5

The D11 certificate must prove:

1. one strict source branch and a frozen pre-split carrier;
2. exactly two independent physical response directions;
3. explicit characters deriving `rho_HT`, `R_T`, `R_H`, `pi_y`, `pi_lambda`,
   and `-16/9` rather than storing them as rules;
4. oriented descent, canonical normalization, and positivity;
5. rigidity or an independently justified positive-gap selector over the
   complete admissible deformation class.

The existing D11 JSON is a conditional formula receipt and does not satisfy
DS1--DS5.

## Trusted checking architecture

Use

```text
frozen carrier -> untrusted generator -> exact certificate -> small checker
```

The checker must recompute path legality and exhaustion from the carrier,
orbit equivalence, response sums and rank, the center representation and trace,
mismatch descent, fibre elimination, and the deformation census. JSON Schema
conformance is only structural validation. The final formal proof should bind
the carrier and certificate hashes and use no unreviewed axioms. Independent
generation is a reproducibility check, not a substitute for the trusted
checker.

## Remaining non-enumeration gates

Even valid QT and DS certificates would emit running/chart coordinates, not
physical masses. Full promotion still requires:

- a unique interval-certified source root on one branch;
- a source clock packet. If the source emits a dimensionless clock gap
  `epsilon_clk`, then `E_star = h nu_clk / epsilon_clk`; using a measured clock
  as calibration must remain labeled as such;
- frozen beta functions, scheme, thresholds, matching maps, and interval
  remainder bounds;
- BRST-complete physical inverse two-point blocks for W, Z, and Higgs,
  including mixing fields, Ward/Slavnov--Taylor and Nielsen identities,
  declared sheets and contours, nonzero Laurent couplings in dressed
  gauge-invariant current amplitudes (and the separate Higgs scalar amplitude), widths, and
  uncertainty bounds;
- a prospective same-branch dependency DAG and claim manifest.

The pole convention is

```text
s_B = (M_B - i Gamma_B/2)^2,  M_B > 0,  Gamma_B >= 0.
```

Until all of these pass, the fail-closed status remains
`source_incomplete_no_end_to_end_wzh_pole_prediction`.
