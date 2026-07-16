# Ward-Projected Hadronic Transport Payload (Generator G1): Status

Protocol correction: `V1_PROTOCOL_ERRATUM_2026-07-16.md` supersedes the old
emit-and-compare procedure. Historical V1 is exploratory and cannot be
promoted. Its target file and runtime artifact remain immutable. The repaired
`run_bracket.py` is source-only; `score_bracket.py` validates the sealed source
artifact and corrective target in a separate process. The current result is
always `NOT_EVALUABLE`: v3 is inactive, the source bracket is an uncertified
singleton diagnostic, and `Delta_EW` is open.

## 1. Declared contract implemented

Per `code/P_derivation/THOMSON_TRANSPORT_THEOREMS.md` (Theorems 3, 5, 6) and
`code/P_derivation/SOURCE_SPECTRAL_THEOREM.md`:

- Emitted object: the source transport packet in inverse-alpha units,
  `Delta_source(P) = Delta_lep(P) + Delta_had(P) + Delta_EW(P)`, with the
  hadronic part as the once-subtracted dispersion moment of a Ward-projected
  `U(1)_Q` spectral density against the kernel `mZ^2/(3*pi*s*(s+mZ^2))`.
- Scheme: same subtraction as `a0(P) = alpha_em^-1(m_Z^2;P)`, family
  `d10_running_tree`, R-ratio normalization (massless parton density
  `N_c Q^2`).
- Screening coordinates: `x = N_c alpha_3(m_Z;P)/pi`; the declared residual
  form is `S = 1 - x + c_Q x^2`; the harness reports
  `S_hadronic = Delta_had / Delta_quark_naive`, conditional
  `S_QEW_effective = (Delta_had + Delta_EW) / Delta_quark_naive`, and
  `c_Q = (S_QEW_effective - (1 - x)) / x^2` for every payload. The two S
  coordinates are not identified merely because the open branch currently
  assigns `Delta_EW = 0`.
- Evaluation point: the internal fixed-point root of the implemented chain,
  `P* = 1.63097209585889737696451390350695563...`, D10 point rebuilt with
  `paper_math.PaperMathContext` (precision 40, su2 cutoff 60, su3 cutoff 45).
  At this point `alpha_3(m_Z;P*) = 0.118335861957773`,
  `a0(P*) = 128.308268057988`, `x = 0.113002424253719`.
- `Delta_EW`: declared zero branch, unproven (Theorem 4 open). It enters
  every payload as an explicit `declared_zero_branch_unproven` field.

No CODATA/NIST alpha, measured hadronic cross section, PDG hadron datum, or
empirical endpoint interval enters the source computation. Section 5 records
the superseded historical comparison and supplies no evidential weight.

## 2. Internal inputs used

Quark and lepton masses: `mass_source = internal_stage5_continuation`
(`paper_math.diagonal_quark_masses`, `charged_lepton_masses` at the rebuilt
D10 vev). Ratios `m/mZ_run` at `P*`:

| species | m/mZ_run    | display GeV (lane clock candidate) |
|---------|-------------|------------------------------------|
| u, d    | 4.08061e-05 | 0.003737                           |
| s       | 1.46902e-03 | 0.134545                           |
| c       | 8.81411e-03 | 0.807270                           |
| b       | 5.28847e-02 | 4.843618                           |
| t       | 1.90385     | 174.370 (excluded from the sum)    |
| e       | 5.57418e-06 | 0.000511                           |
| mu      | 1.15257e-03 | 0.105562                           |
| tau     | 1.93840e-02 | 1.775353                           |

alpha_s: four-loop MS-bar running of the chain's own `alpha_3(m_Z;P*)` with
three-loop decoupling, thresholds at the internal Stage-5 `m_c`, `m_b`,
`m_t` (engine ported from the source-only lane
`code/particles/qcd/derive_lambda_qcd_source_transmutation.py`).

Lambda_QCD scale: dimensionless ratio `Lambda3/mZ` from the source-only
transmutation lane artifact
`code/particles/runs/qcd/lambda_qcd_source_transmutation.json`
(central 3.6556e-03, lane threshold-sweep interval
[3.4884e-03, 3.8187e-03]). Only the ratio is consumed, so the unclosed
clock candidate cancels.

## 3. Computed variants and the development bracket

Declared grid (58 runs, `run_bracket.py`, stated in advance, not tuned):

- `parton_free`: free-quark one-loop dispersion, Stage-5 masses. Reproduces
  the chain's naive quark sum exactly (`S_hadronic = 1`).
- `pqcd`: parton density times the massless MS-bar R-ratio series
  `1 + a + c2 a^2 + c3 a^3` (`a = alpha_s(s)/pi`) above the declared IR
  cutoff `s0 = (k Lambda3)^2`; grid: Lambda3 {lane_lo, lane_central,
  lane_hi} x k {2, 4, 8} x below-cutoff {free parton, zero support} x
  truncation order {1, 2, 3}.
- `constituent`: dressed masses `M_q = sqrt(m_q^2 + (kappa Lambda3)^2)`,
  kappa = 1 declared, Lambda3 over the lane interval.

Chain reference at `P*` (implemented screen `1 - x`):
`Delta_lep = 4.309397144817`, `Delta_quark_naive = 4.934816164436`,
`Delta_quark_screened_impl = 4.377169974609`,
`Delta_impl_total = 8.686567119425`.

Key historical rows (Delta values in inverse-alpha units; `Delta_source total`
means `Delta_lep + Delta_had + Delta_EW`):

| variant                        | S_hadronic | Delta_had | Delta_source total |
|--------------------------------|----------|-----------|--------------|
| parton_free                    | 1.000000 | 4.934816  | 9.244213     |
| pqcd o3 k=4 free (central)     | 1.038029 | 5.122481  | 9.431878     |
| pqcd o3 k=4 zero (central)     | 0.655241 | 3.233492  | 7.542889     |
| constituent kappa=1 (central)  | 0.659106 | 3.252567  | 7.561964     |
| min: pqcd o1 k=8 zero (hi)     | 0.557779 | 2.752536  | 7.061933     |
| max: pqcd o2 k=2 free (lo)     | 1.054347 | 5.203008  | 9.512406     |

Development bracket over the full grid:

- `Delta_had`: [2.752536, 5.203008]
- `Delta_source total`: [7.061933, 9.512406], width 2.450472
- `S_hadronic`: [0.557779, 1.054347], width 0.496568
- implied `c_Q`: [-25.78, +13.11]
- residual `R_Q` vs the implemented screen: [-1.624634, +0.825838]

Historical artifact: `runtime/ward_projected_payload_bracket_current.json`
(deterministic, content-hashed; full-grid wall time 0.36 s plus 3.8 s for
the D10 evaluation-point rebuild). It is retained byte-for-byte and is not an
artifact of the repaired schema. New emissions default to
`runtime/ward_projected_payload_source_bracket_current.json`; protocol runs
should write to a sealed path outside the source checkout.

## 4. Historical v2 precision diagnostic, not a scoring contract

V2's scalar pass rule is algebraically invalid: it compares a source total
with a residual, assigns one target to two different completed maps, and
treats the open `Delta_EW` remainder as zero. The following figures are
preserved only as historical scale diagnostics. They define no valid pass or
falsification threshold.

- Historical v2 tolerance field: 2.1e-8 in `Delta_source` inverse-alpha units.
  Its implied relative scale on `Delta_had` was 4.0e-09.
- Development bracket width: 2.450472 inverse-alpha units. Ratio to
  tolerance: 1.17e+08.

Comparison against the standard literature:

- Data-driven dispersive determinations, which integrate measured
  `e+e- -> hadrons` R-ratio data, reach
  `Delta alpha_had^(5)(mZ^2) = 0.02766(7)` (Davier et al. 2019) and
  `0.02761(11)` (Keshavarzi-Nomura-Teubner 2019). In inverse-alpha units the
  hadronic transport is about 3.8 with uncertainty 0.010 to 0.015, a
  relative precision of 2.5e-3 to 4e-3. That sits a factor of roughly 5e5 to
  7e5 above that historical scale, and the method consumes measured hadronic
  cross sections, which the source-only lane forbids.
- Lattice QCD(+QED) first-principles hadronic vacuum polarization: BMW 2020
  (Borsanyi et al., Nature 593, 51) gives `a_mu^LO-HVP = 707.5(5.5)e-10`
  (0.8% relative); the 2024 BMW/DMZ update and the Mainz/CLS window program
  reach roughly 0.4% to 1%. Direct lattice computations of the running to
  `mZ^2` via the Adler function (Ce et al. 2022) carry about 0.5% to 1% on
  the hadronic transport. That sits a factor of roughly 1e6 to 2e6 above the
  historical scale. Published lattice scale setting consumes one measured
  dimensionful hadronic input, so even this route is only
  hadron-input-free after replacing scale setting with a source-emitted
  scale.
- Scaling: shrinking a 0.5% lattice error to 4e-9 relative requires a 1e6
  error reduction, a naive 1e12 statistics factor, and simultaneous control
  of continuum, finite-volume, QED/isospin, and scale-setting systematics at
  the 1e-9 relative level. Current world programs (Fermilab/CERN g-2 theory
  initiative targets) aim at the 1e-3 relative level.
- Perturbation theory cannot substitute below about 2 GeV: at
  `sqrt(s) = 1 GeV` the internal chain gives `alpha_s/pi ~ 0.15`, the known
  R-ratio terms shift the density at the several-percent level per order,
  and the region below 2 GeV carries an O(1) fraction of the transport
  moment. The spread across the declared IR-cutoff and support branches
  (this bracket) is the quantitative image of that failure.

Statement about the historical v2 scale: no current method reaches a width of
2.1e-8 inverse-alpha units (4e-9 relative on
the hadronic transport) without measured hadronic input. No current method
reaches it with measured hadronic input either; the best existing
determinations of the hadronic transport sit five to six orders of magnitude
above that historical scale. A first-principles determination at that width
would need a lattice QCD+QED program with total error budgets
roughly six orders of magnitude beyond the current state of the art, with
scale setting emitted by the source chain instead of a measured hadron mass.

## 5. Superseded non-blind development comparison

This session read canon documents that quote the target-side diagnostics, so
the following comparison carries no blind weight. Canon scales seen in
`THOMSON_TRANSPORT_THEOREMS.md`: `S_required(P_C) = 0.895400132648`,
`c_Q(P_C) = 0.658025759927`, `Delta_missing(root) = 0.041163999587`.

- The historical `S_required` lies inside the development `S_hadronic` bracket
  [0.557779, 1.054347]. The bracket width in `S` units (0.4966) exceeds the
  tolerance in `S` units (4.3e-9) by the same 1.2e8 factor.
- Containment carries no promotion weight. Theorem 6 of the canon file
  establishes that the current corpus invariants do not determine the
  spectral moment; a bracket that contains the required scale is consistent
  with that no-go and does nothing to close it.

## 6. Repaired emission and scoring protocol

1. Build an isolated emitter environment containing only source code and
   allowed source-side inputs. Exclude `falsification/`, score artifacts,
   target-bearing prose, and the network.
2. Run `run_bracket.py` at its source-derived internal Stage-5 P. It exposes no
   P override and rejects measurement-located or target-supplied P values.
3. Write the sampled singleton diagnostic to a sealed output path. Any grid,
   node-count, or cutoff change creates a new protocol. Do not describe its
   variant extrema as numerical, theory, or interval-certified error bounds.
4. Record the file hash and embedded `content_sha256` before transfer. Make no
   post-emission edit.
5. In a separate environment, run `score_bracket.py` on the sealed artifact
   and `hadronic_closure_target_2026-07-16_v3.json`. The scorer fails closed on
   target activation, hash, schema, coordinate-kind, P-domain, provenance,
   certification, or `Delta_EW` gate failures.
6. Expect `NOT_EVALUABLE`. V3 is permanently inactive. A detached successor
   may permit a real score only after blind or audited clean-room method
   selection, external anchoring before payload work, a preregistered
   three-way decision policy, and an eligible target-blind P-domain object.
   CL-1 and CL-2 must then be solved and certified independently; no scalar
   containment shortcut is valid.

## 7. Files

- `payload_harness.py`: contract, evaluation point, closed-form parton
  moments, deterministic quadrature, `emit_delta_source`.
- `spectral_modules.py`: alpha_s engine, R-ratio series, declared variant
  modules, synthetic test module.
- `run_bracket.py`: target-free declared grid and source-bracket emission.
- `score_bracket.py`: separate fail-closed target validator and future scorer
  gate; it currently issues no closure verdict.
- `test_payload_harness.py`: source harness and bracket tests.
- `test_score_bracket.py`: source-contamination, immutable-hash, schema,
  coordinate, evaluation-point, and historical V1 rejection tests.
- `runtime/ward_projected_payload_bracket_current.json`: immutable historical
  V1 artifact.
- `V1_PROTOCOL_ERRATUM_2026-07-16.md`: corrected protocol and permanent
  non-promotion ruling.
