# P Derivation

This directory is a clean paper-math implementation of the OPH `P <-> alpha`
closure experiment.

The goal is to avoid public-facing ambiguity between the paper
equations and the larger `code/particles` calibration stack. The code here is
therefore built directly from the equations stated in:

- `paper/deriving_the_particle_zoo_from_observer_consistency.tex`
- `paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.tex`
- `paper/observers_are_all_you_need.tex`

## What is implemented

For a trial pixel constant `P`, the code reproduces the paper D10 forward map:

1. `M_U(P) = E_P * exp(-2*pi) * P^(1/6)`
2. `E_cell(P) = E_P / sqrt(P)`
3. Solve the one-dimensional D10 pixel-closure equation for `alpha_U(P)`:
   `ellbar_SU(2)(t2) + ellbar_SU(3)(t3) - P/4 = 0`
4. Build the source-locked electroweak anchor
   `a0(P) = alpha_em^-1(m_Z^2; P)`

The closure step inserts that paper-side `P -> alpha` map into the outer
equation:

`P = phi + alpha * sqrt(pi)`

The first-principles chain is:

```text
phi from golden-ratio entropy balance
sqrt(pi) from boundary Gaussian normalization
P -> M_U(P) -> alpha_U(P) -> alpha_i(m_Z;P) -> a0(P)
a0(P) -> A_T(P) by Ward-projected U(1)_Q Thomson transport
P = phi + sqrt(pi) / A_T(P)
alpha(0) = 1 / A_T(P)
```

The solver has no default reference constant.

`emit_p_closure_trunk.py` is the compressed five-equation trunk emitter. It
packages the code path as:

```text
P -> M_U -> alpha_U -> alpha_i(m_Z) -> a0(P) -> alpha_in(P) -> P
```

and writes `runtime/p_closure_trunk_current.json`. That artifact is the
canonical audit surface for the simplified chain, but it is not the certified
particle root. Source-only promotion requires a populated source spectral
measure payload, same-scheme remainder, and an interval-level fixed-point
certificate. The public display row can also use the OPH plus empirical hadron
closure surface defined in `../../docs/HADRON.md`.

## Terms

Canonical naming, used on every surface in this directory: `P` is the
measured-endpoint working value of the pixel constant, located from the
measured fine-structure constant; `P_fwd` is the forward-model closure
point, the pixel value implied by the solver's own contracted output. The gap
between them is the open hadronic loop residual tracked on the P-closure
issue ([#545](https://github.com/FloatingPragma/observer-patch-holography/issues/545)).

`P ratio (pixel size)` means the dimensionless screen-cell area

`P = a_cell / l_P^2`

so it is the outer pixel area measured in Planck units.

`inside minimal observation` means the electromagnetic coupling that acts as the
smallest durable observation scale inside the emitted universe. In the strict
D10 core this is the source anchor `a0(P) = alpha_em^-1(m_Z^2; P)`. In the
full low-energy extension it is the Thomson-limit coupling after internal
running.

Informally, the closure program is:

- the outer description says how far the realized pixel sits above the
  golden-ratio equilibrium, via `P = phi + alpha * sqrt(pi)`
- the inner description says what electromagnetic observation scale that same
  pixel emits through the OPH chain
- the derivation solves for the fixed-point `alpha` where those two
  descriptions agree, and only then computes `P`

## Readout modes

Two alpha readout modes are supported:

- `thomson_structured_running`
  This is the default. It uses the D10 core for
  `a0(P) = alpha_em^-1(m_Z^2;P)`, then computes
  `Delta alpha^-1(P)` from the internal Stage-5 charged-spectrum continuation
  with the exact one-loop fermion transport kernel
  `K_f(Q^2;m_f) = (2 N_c Q_f^2 / pi) * integral_0^1 x(1-x) log(1 + Q^2 x(1-x)/m_f^2) dx`
  together with the quark screening factor `1 - N_c alpha_s(P)/pi`.
  No paper endpoint, frozen ratio, or imported charged bundle is inserted.

- `thomson_structured_running_asymptotic`
  Legacy comparison path. It keeps the older high-energy asymptotic
  approximation
  `Delta alpha_f^-1 ~= N_c Q_f^2 (log(Q^2 / m_f^2) - 5/3) / (3 pi)`.
  This is preserved for auditability only.

- `mz_anchor`
  Uses the source anchor `a0(P) = alpha_em^-1(m_Z^2; P)` directly. This is
  mainly a debugging surface.

The important claim-boundary caveat is:

- `paper_math.py` is zero-insert with respect to paper-side target values:
  it contains no hard-coded reference `P`, no hard-coded Thomson endpoint, and
  no imported compare bundle from `code/particles`.
- `thomson_structured_running` is therefore the cleanest internal closure
  experiment in this directory, and it uses the exact one-loop fermion kernel
  rather than the older asymptotic log expansion.
- It is a continuation beyond the strict theorem-grade D10 core, because
  the final low-energy transport law is being modeled by the internal
  Stage-5 structured-running ansatz rather than by a closed theorem.
- `fixed_point_witness.py` emits a numerical witness, not an interval
  certificate. It samples the declared `alpha -> alpha` map and records local
  finite-difference slopes, but it does not prove a Banach contraction bound or
  interval-wide uniqueness.
- `fixed_point_certificate.py` emits a local numerical contraction certificate
  for the implemented map. It is a stricter machine artifact than the witness;
  the formal interval-arithmetic certificate for the full map is
  `interval_contraction_certificate.py` (see "Interval certificates" below).
- A separate pending hardware note reports an optical-cavity check of the same
  fixed-point geometry; this is treated as corroborating engineering evidence.

## Full derivation claim boundary

`FULL_DERIVATION.md` records the complete derivation contract and the endpoint
audit packet. `THOMSON_TRANSPORT_THEOREMS.md` records the theorem suite and its
source-payload rule. The short version is:

- the D10 source map and the outer fixed-point witness are implemented
- the measured input endpoint is
  `alpha^-1(0) = 137.035999177(21)` with
  `P=1.630968209403959324879279847782648941...`; this is the measured input
  that locates the working `P`, never a solver product
- the solver output is
  `alpha^-1 = 136.994835177412937295289429464436...`
  (default exact one-loop readout, converged precision-100 rerun;
  interval-certified unique fixed point, and by the domain-global certificate
  below the only fixed point of the readout map on the declared physical
  domain `alpha^-1` in `[100, 200]`), defining the
  forward-model closure point `P_fwd`
- the difference between the solver output and the measured input endpoint is
  the open hadronic term: the loop residual tracked on the P-closure issue
  (#545); no target or scorer exists
- the missing term is a source-only Thomson transport contribution of
  `0.0411639995870627047105705355631...` in inverse-alpha units
- at the measured-endpoint working pixel
  `P=1.630968209403959324879279847782648941...`, the endpoint package gives
  `a0(P)=128.307965473286248209961108741756716187...`,
  `Delta_required(P)=8.728033703713751790038891258243283813...`, and
  `Delta_source_residual(P)=0.041465861005223389053448715357314044...`
- the residual is equivalent to a required Ward-projected quark-screening factor
  `S_required=0.895400132647658797805800283181670641...`, or
  `c_Q=0.658025759927155435638230170232360050...` in the parameterization
  `S=1-x+c_Q x^2`, with `x=N_c alpha_3(m_Z;P)/pi`

Run the audit after producing a report:

```bash
python3 alpha_gap_audit.py --report runtime/full_p_alpha_report_current.json
python3 thomson_endpoint_package.py --report runtime/full_p_alpha_report_current.json
python3 screening_invariant_no_go.py
python3 thomson_endpoint_interval_certificate.py
python3 transport_theorem_manifest.py --report runtime/full_p_alpha_report_current.json
python3 measured_endpoint_calibration.py
```

## Interval certificates

`interval_contraction_certificate.py` emits
`runtime/p_interval_contraction_certificate_2026-07-14.json`, the stage-2
contraction certificate of the basin-then-contract protocol for the P
coordinate (P-closure issue #545).

What is certified:

- Existence and uniqueness of the closure-map fixed point on a stated alpha
  interval, by direct Banach with a mean-value (centered) form: the interval
  evaluation proves `sup |g'| <= L < 1` on `I` and
  `g(I) subset g(mid(I)) + g'(I)(I - mid(I)) subset interior(I)`.
- Two readout maps carry certificates: the declared solver mode
  `thomson_structured_running` (fixed point near `alpha^-1 = 136.9948`), and
  the gauge-width (#545) mixed map `thomson_structured_running_plus_gauge_width`
  that adds the unified gauge width `alpha_U(P)` to the inverse-alpha readout
  (certified fixed point `alpha^-1 = 137.035660136946577...`, the gauge-width fixed-point value).
- The full declared chain is evaluated in binary interval arithmetic
  (`mpmath.iv`, outward rounding on every elementary operation; private
  contexts, global mpmath precision untouched). The two implicit solves
  (tree-level `m_Z`, D10 pixel closure for `alpha_U`) are enclosed by verified
  sign-change brackets with sign-definite interval residual derivatives
  (existence, in-bracket uniqueness, and C^1 dependence by the implicit
  function theorem). The map derivative is enclosed by forward-mode dual
  arithmetic over intervals with the implicit nodes handled by the same
  theorem.
- The SU(2)/SU(3) edge-sum truncation tails are bounded by explicit geometric
  majorants and added one-sidedly to the interval enclosures, so the
  certificate covers both the declared cutoffs (`su2_cutoff=120`,
  `su3_cutoff=90`) and the infinite-cutoff edge sums.
- The repository Decimal pad backend (`interval_backend.py`) is not used by
  this artifact; the `backend` field records the mpmath.iv basis, so the
  directed-rounding caveat of that module does not attach here.

What is not certified:

- The declared one-loop RG/matching conventions, the tree-level `m_Z`
  closure, the Stage-5 continuation masses, and the exact one-loop kernel are
  certified as declared numerical structure, not as physical endpoint
  theorems.
- No relation to the measured fine-structure constant is certified. The P
  coordinate carries a stage-2 contraction certificate for the declared
  map; the stage-3 landing verdict is unchanged: the source fixed point
  remains outside the measured basin (P-closure issue #545), pending the hadronic
  transport term.
- Global uniqueness of the inner roots (`alpha_U`, `m_Z`) over the solver's
  full scan windows is not claimed; the inner-root certificates hold on the
  verified brackets, which contain the roots the declared scan-and-bisect
  procedure selects.

Run and test:

```bash
python3 interval_contraction_certificate.py --mp-dps 60 --iv-dps 60 \
    --su2-cutoff 120 --su3-cutoff 90 --half-width 0.000004 --refine-passes 10
python3 -m pytest test_interval_contraction_certificate.py -q
```

### Global uniqueness on the declared domain

`global_uniqueness_certificate.py` emits
`runtime/p_global_uniqueness_certificate_2026-07-14.json`, the global
at-most-one supplement to the stage-2 certificate. The declared physical
domain is the solver scan window of `paper_math.solve_closure`:
`alpha in [0.005, 0.01]` (`alpha_inv in [100, 200]`), containing both
certified fixed points. The domain is covered by an adaptive interval
subdivision (256 pieces per readout map, adaptive bisection with a depth cap
and a deterministic work budget; the production sweep needed no bisection);
on every piece the stage-2 forward-mode interval AD chain certifies
`sup |g'| < 1` (worst piece `L <= 0.3041` source mode, `L <= 0.3039`
gauge-width mode; cutoffs `su2_cutoff=120`, `su3_cutoff=90` with the edge-sum
tail majorants folded in, so the verdicts also cover the infinite-cutoff
sums; dps 40; exceptional set empty). `|g'| < 1` on the convex domain gives
at most one fixed point per declared map by the mean value theorem; combined
with the stage-2 existence certificate (same maps, same cutoffs), each map
has exactly one fixed point on the declared domain. Fallback verdict tiers
(sign-definite residual derivative, root-free residual enclosure, monotone-run
synthesis with certified endpoint residual signs) are implemented and tested;
the production sweep did not need them. The claim boundary of the stage-2
certificate carries over unchanged; the inner-root scan-window caveat above
is unaffected.

```bash
python3 global_uniqueness_certificate.py --mp-dps 40 --iv-dps 40 \
    --su2-cutoff 120 --su3-cutoff 90 --initial-pieces 256
python3 -m pytest test_global_uniqueness_certificate.py -q
```

### Maximal-domain uniqueness extension (GAP-A7)

`global_uniqueness_extension_certificate.py` emits
`runtime/p_global_uniqueness_extension_certificate_2026-07-17.json`, the
maximal-domain supplement that closes proof-spine gap GAP-A7. The declared
maps read alpha only through `P = phi + alpha*sqrt(pi)`, and every point of
the maximal analytic domain carries a window-consistent triple
`(alpha_U, L = ln(mu_U/m_Z), P)` inside the declared solver windows
(`alpha_U in [0.02, 0.08]`, `L in [0, 50]`). An envelope lemma from the
window bounds and the `m_Z`-closure identity certifies
`1/alpha_3(m_Z) >= 4`, positive quark screening, a global `P` range, and a
positive inverse-alpha readout floor on the whole domain. An adaptive
interval sweep of the full pixel window (constraint propagation on the
closure equations; 78 pieces at the production parameters; exceptional set
empty) certifies for every piece one of: the `m_Z` window is infeasible, the
closure-forced alpha value is disjoint from the certified readout enclosure
for both readout maps, or the alpha value lies inside the declared physical
interval. Every fixed point of either map on the maximal analytic domain
therefore lies inside `alpha in [0.005, 0.01]`; the exterior fixed-point set
is empty, and the composition with the stage-2 existence certificate and the
at-most-one certificate gives exactly one fixed point per map on the maximal
domain. The sweep quantifies over every window-consistent inner-root
selection, so the exterior exclusion covers the declared scan-and-bisect
selection. The certificate also records the certified outer alpha bound of
the domain and point-located domain edges under the declared grid semantics
(right edge inside `(0.4125, 0.4137)`, left edge inside
`(-0.5594, -0.5578)`).

```bash
python3 global_uniqueness_extension_certificate.py --mp-dps 40 --iv-dps 40 \
    --su2-cutoff 120 --su3-cutoff 90
python3 -m pytest test_global_uniqueness_extension_certificate.py -q
```

The closed printed-pair identity is enforced by
`test_printed_pair_identity.py`: the runtime trunk and full-report artifacts
are regenerated at solver precision 100 with enough outer iterations that
`alpha_root = (P_fwd - phi)/sqrt(pi)` holds to at least 30 significant digits
(stated test tolerance 1e-30 relative; measured 3.1e-36 report, 1.2e-38
trunk). Regeneration commands:

```bash
python3 emit_p_closure_trunk.py --precision 100 --max-iterations 120
python3 -c "import json; from paper_math import build_report; \
    report = build_report(precision=100, mode='thomson_structured_running', \
    su2_cutoff=40, su3_cutoff=30, scan_points=60, max_iterations=110); \
    open('runtime/full_p_alpha_report_current.json','w').write( \
    json.dumps(report, indent=2, sort_keys=True) + chr(10))"
```

The convergence depth of the internal bisections scales with `--precision`
(`max(32, precision + 8)` steps), so the emitter default precision keeps its
fast profile and the identity-grade artifacts are produced with the explicit
flags above.

## Empirical endpoint insert

`measured_endpoint_calibration.py` emits
`runtime/measured_endpoint_calibration_current.json`. This is the OPH plus
empirical hadron closure surface for tables, plots, and public numeric surfaces
that need a single fine-structure value. The source-only spectral payload stays
separate in the audit ledgers.

The measured input endpoint surface is

```text
alpha^-1(0) = 137.035999177(21)   (measured input; locates the working P)
alpha(0)   = 0.007297352564331425...
P          = 1.63096820940395932487927984778...
```

The JSON keeps audit guard fields for theorem tooling and carries row class
`oph_plus_empirical_hadron_closure`. Consumer-facing prose leads with the
source-only value and the closure status: the forward map contracts to
`alpha^-1 = 136.9948...`, the residual to the measured basin is the open hadronic
term (P-closure issue #545), and the measured endpoint
`137.035999177(21)` is the measured input. The empirical hadron-closure row is a
labeled closure-status display:

```text
source-only forward value 136.9948...; the gap to the measured input endpoint is the open hadronic closure term (P-closure issue #545); the OPH plus empirical hadron closure row states closure status, never the OPH result
```

The source computation is blocked because the required object is a
nonperturbative Ward-projected hadronic spectral measure, not a scalar that can
be recovered by higher precision on the existing formulas. A source computation
requires a production finite-volume gauge-field run, vector-current levels and
residues, current normalization, continuum/chiral/finite-volume/current-matching
budgets, a spectral pushforward, quadrature and tail bounds, and a directed
interval fixed-point certificate. The corpus has a schema and contract
for that payload. The empirical hadron closure path uses a separate
`e+e- -> hadrons` payload class through the policy in `../../docs/HADRON.md`.

`EMPIRICAL_HADRON_SCHEME_BRIDGE.md` records the direct PDG/CERN
`Delta alpha_had^(5)(M_Z)` diagnostic. Raw PDG/CERN hadronic-running values are
dimensionless denominator shifts, not additive OPH inverse-alpha packets. With
the current OPH source anchor and lepton packet, direct insertion of the PDG
row `Delta alpha_had^(5)(M_Z)=0.02761` gives
`alpha^-1(0)=136.382895072695577...`. Hitting the NIST/CODATA comparison value
with that same hadron row would require a same-scheme source-anchor bridge of
`0.6350718999845777629...` inverse-alpha units, or a hadronic denominator shift
of `0.0322443435578872888...`, outside the PDG/CERN range. Therefore raw PDG
hadronic running is an empirical diagnostic until a source-side scheme bridge or
same-scheme Ward-projected endpoint map is supplied.

## Hierarchy certificate handoff

The electroweak-hierarchy proof bundle in `../particles/hierarchy` uses the
same two branch surfaces recorded here. On the measured-endpoint working branch
(canonical name `P`; bundle label `P_C`) it reads

```text
P_C = 1.630968209403959324879279847782648941
alpha_U(P_C) = 0.041124336195630495
v/E_star = 2.0199803239725553e-17
```

On the forward-model closure branch (canonical name `P_fwd`; bundle label
`P_source_audit`) it reads

```text
P_source_audit = 1.63097209585889737696451390350695562847912625483895268486516
alpha_U(P_source_audit) = 0.0411242474418166851408899338896597194
v(P_source_audit)/E_star = 2.0198114078576330591337426182670149e-17
```

The bundle certifies the local `P -> alpha_U -> v/E_star` hierarchy lane by a
declared dependency graph and a Krawczyk inclusion for the `R_U` interval. The
public Thomson endpoint transport gate remains the same source-side hadronic
spectral object described above.

## Usage

From `reverse-engineering-reality/code/P_derivation/`:

```bash
python3 derive_p.py --mode thomson_structured_running --precision 40
```

The `derive_p.py` command is the human-facing fine-structure CLI. It prints a
colored progress display by default, including the scan stage, the fixed-point
bracket, each bisection step, the source anchor, transport packet, empirical
closure packet, and CODATA/NIST comparison. Use:

```bash
python3 derive_p.py --help
python3 derive_p.py --color always
python3 derive_p.py --no-hadron-closure
python3 derive_p.py --json --output runtime/report.json
```

The default display row adds the ledger closure packet
`0.041164012378350542050005414212806737971` in inverse-alpha units, the
constant hard-coded in `derive_p.py`. This is
the OPH plus empirical hadron closure surface, not a source-only hadron theorem.
That packet constant predates the converged precision-100 reruns; the converged
forward gap is `0.041163999587062704...` inverse-alpha units (printed-pair identity,
closed).
Increase `--max-iterations` when increasing `--precision` and treating extra
digits as audit-relevant.

If you only want the shortest possible entrypoint that prints the current
candidate fine-structure value, use:

```bash
python3 minimal_alpha.py --precision 40
```

To run the compact PDG/CERN hadron-input diagnostic:

```bash
python3 fine_structure_fixed_point_demo.py
python3 fine_structure_fixed_point_demo.py --compare-alpha-inv 137.035999177
```

To emit a full JSON report:

```bash
python3 derive_p.py --mode thomson_structured_running --precision 40 --json
```

To save the report:

```bash
python3 derive_p.py --mode thomson_structured_running --precision 40 --output runtime/report.json
```

To emit the fixed-point witness:

```bash
python3 fixed_point_witness.py --mode thomson_structured_running --precision 40 --output runtime/fixed_point_witness.json
```

To emit the compressed P-trunk artifact:

```bash
python3 emit_p_closure_trunk.py --output runtime/p_closure_trunk_current.json
```

The default trunk emitter uses the paper-compression asymptotic structured
running profile so it is fast enough for routine ledger refreshes. To run the
stronger exact one-loop continuation profile, pass:

```bash
python3 emit_p_closure_trunk.py --mode thomson_structured_running --precision 40 --su2-cutoff 120 --su3-cutoff 90 --output runtime/p_closure_trunk_exact.json
```

To include a reference inverse-alpha value in an audit report, pass it
explicitly:

```bash
python3 fixed_point_witness.py --mode thomson_structured_running --precision 40 --compare-alpha-inv <external-alpha-inv> --output runtime/fixed_point_witness.json
```

To emit the local numerical contraction certificate:

```bash
python3 fixed_point_certificate.py --mode thomson_structured_running --precision 40 --output runtime/fixed_point_certificate.json
```

For a fast smoke check of the witness plumbing, use the electroweak-scale anchor
debug path:

```bash
python3 fixed_point_witness.py --mode mz_anchor --precision 10 --su2-cutoff 6 --su3-cutoff 4 --scan-points 8 --max-iterations 3 --sample-points 1
python3 fixed_point_certificate.py --mode mz_anchor --precision 10 --su2-cutoff 6 --su3-cutoff 4 --scan-points 8 --max-iterations 3 --interval-half-width 0.0001 --derivative-step 0.0001 --sample-points 3
```

## Output

The CLI prints:

- the solved `alpha` and `alpha^-1`
- the implied `P`
- the source anchor `a0(P)`
- the D10 point at the closure solution
- the internal structured-running decomposition used to reach the Thomson limit

The JSON report also includes every alpha-space closure step.

The witness JSON additionally records:

- `claim_status = numerical_witness_not_interval_certificate`
- the sampled local slopes of the closure map
- any explicitly supplied reference inverse-alpha value
- the pixel ratio implied by that reference value

The certificate JSON additionally records:

- `claim_status = numerical_local_contraction_certificate` when the sampled
  interval brackets the fixed point and all sampled finite-difference slopes
  have absolute value below one
- the explicit alpha interval and endpoint residuals
- the maximum sampled contraction slope and per-sample contraction margins

## License

This code surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE).
