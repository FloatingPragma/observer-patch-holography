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
particle root. Promotion requires a populated source spectral measure payload,
same-scheme remainder, and an interval-level fixed-point certificate.

## Terms

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
  for the implemented map. It is a stricter machine artifact than the witness,
  but still not a formal interval-arithmetic proof of the full map.
- A separate pending hardware note reports an optical-cavity check of the same
  fixed-point geometry; this is treated as corroborating engineering evidence.

## Full derivation status

`FULL_DERIVATION.md` records the complete derivation contract and the endpoint
audit packet. `THOMSON_TRANSPORT_THEOREMS.md` records the theorem suite and its
source-payload rule. The short version is:

- the D10 source map and the outer fixed-point witness are implemented
- the public fixed-point readout is
  `alpha^-1(0) = 137.035999177(21)` and
  `P=1.630968209403959324879279847782648941...`
- the default exact one-loop readout gives
  `alpha^-1 = 136.994835164621649457949994585787...`
- the missing term is a source-only Thomson transport contribution of
  `0.0411640123783505420500054142128...` in inverse-alpha units
- at the public endpoint pixel
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

## Measured endpoint calibration insert

`measured_endpoint_calibration.py` emits
`runtime/measured_endpoint_calibration_current.json`. This is the display
surface for tables, plots, and public numeric surfaces that need a single
fine-structure value, with the source spectral payload recorded separately in
the audit ledgers.

The endpoint surface is

```text
alpha^-1(0) = 137.035999177(21)
alpha(0)   = 0.007297352564331425...
P          = 1.63096820940395932487927984778...
```

The JSON keeps audit guard fields for theorem tooling. Consumer-facing prose
should lead with the fixed-point readout and keep source-audit status in
ledger sections:

```text
fine-structure endpoint uses OPH fixed-point readout; finite source audit recorded in code/P_derivation ledgers
```

The source computation is blocked because the required object is a
nonperturbative Ward-projected hadronic spectral measure, not a scalar that can
be recovered by higher precision on the existing formulas. A source computation
requires a production finite-volume gauge-field run, vector-current levels and
residues, current normalization, continuum/chiral/finite-volume/current-matching
budgets, a spectral pushforward, quadrature and tail bounds, and a directed
interval fixed-point certificate. The current corpus has a schema and contract
for that payload, but not the production payload itself.

## Usage

From `reverse-engineering-reality/code/P_derivation/`:

```bash
python3 derive_p.py --mode thomson_structured_running --precision 40
```

If you only want the shortest possible entrypoint that prints the current
candidate fine-structure value, use:

```bash
python3 minimal_alpha.py --precision 40
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
