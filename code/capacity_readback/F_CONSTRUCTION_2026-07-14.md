# Construction Attempt for the Capacity Readback Map F (2026-07-14)

Candidate construction run for generator G2 against
[F_READBACK_SPEC.md](F_READBACK_SPEC.md). Target: a concrete candidate for the
self-consistent screen-configuration count `|Omega^sc_N|` derived from the declared
structure only, the induced readback map `F = Cap_read o Obs o nf`, and certified
fixed points per branch. Anti-fitting discipline: every underdetermined reading is
enumerated and carried as a branch; no branch is selected by its landing; all fixed
points are recorded before any comparison with a reference capacity.

Status header, written after execution: **no branch lands near a reference capacity.
CL-7 stays open. Outcome class (b) with one named gap item (G2-GAP-1).** The
constructed candidates are recorded here as excluded.

## 0. Inputs and blindness cone

The construction consumes:

- `P` only as the certified forward closure enclosure from
  [`code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json`](../P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json),
  `P ∈ [1.6309720958588973769645139031692..., 1.6309720958588973769645139038446...]`
  (display value `P_fwd = 1.630972095858897`). No measured endpoint, no SL-3
  measured-endpoint root `1.6309682094`.
- Integers supplied by declared theorems: 12, 24, 6, 20, 60, and the constant `pi`.
- No measured Lambda, no `3.31e122`, no `N_EW = pi*exp(6*pi/(P*alpha_U))` value, no
  `alpha_U` (Section 4 records the one reading that would import it and excludes it
  under V-08).

## 1. Derivation, numbered steps

### Step 1. Object and units

`F(N) = Cap_read(Obs(nf(U_N)))`, closure equation `N_CRC = F(N_CRC)`, count
representation `Pi(N) = |Omega^sc_N| * e^(-N)`, selector
`N_star = MAR argmax_N [log|Omega^sc_N| - N]`.
Source: `paper/recovering_relativity_..._compact.tex`, Definition
`def:self-closure-density` and Remark `rem:self-closure-counting-target`;
`paper/tex_fragments/OBSERVERS_SYNTHESIS_SECTIONS.tex`, global self-closure section
(finite readback form, contraction condition `0 <= F' <= kappa < 1`). All capacities
in nats, `log dim H_partial,N = N`.

### Step 2. Screen decomposition: cells and the P/4 budget

`def:self-closure-density` supplies the equal-area chart count

```
K_cell = A_screen / a_cell = 4N / P.
```

Each cell carries `P/4` nats of cut density. The same definition carries an explicit
warning: `P/4 ~ 0.408` nats is an edge/cut density and is stated as such; it is not
`log dim H_cell` for an autonomous cell factor, and a fixed elementary carrier would
need extra hypotheses (homogeneous cellulation, isomorphic cell algebras, equivariant
overlaps, refinement preservation, carrier selection).

**Branch point BR-0 (menu size 2).** Reading (i): treat the `K = 4N/P` cells as
independent record carriers of `P/4` nats each for counting purposes. Reading (ii):
refuse the product structure, in which case no configuration count is derivable at
this step and the run ends in outcome class (c). Reading (i) is the only executable
reading; reading (ii) is recorded as the fallback interpretation of the same warning.
The run carries (i) and registers the warning as an inherited hypothesis on every
branch below.

Consistency check for reading (i): the raw product count is
`(e^(P/4))^K = e^(K*P/4) = e^N = dim H_partial,N`. The cell reading reproduces the
declared screen dimension exactly. This is compatibility, and it is the reason
reading (i) is executable at all.

### Step 3. Repair closure: the Z6 reserve. Branch axes BR-1, BR-2

Membership clause 1 of `Omega^sc_N`: closed under repair. The declared repair pricing
is the protected Z6 reserve
(`paper/screen_microphysics_and_observer_synchronization.tex`, Theorem
`thm:z6-reserve-trace` with declared branch input `l_shared = P/4`;
Definition `def:z6-reserve-projector`): reserve mean `tau_q(Z6) = P/24` per
accounting unit, Poisson reserve exponent `e^(-P/24)`. Every number downstream
inherits the `l_shared = P/4` branch input and the MaxEnt-maximal-symmetry named
hypothesis (remark after `thm:icosahedral-screen-sieve`).

**Branch point BR-1, reserve semantics (menu size 3).**

- **BR-1a (poisson):** repair-closed survival factor `e^(-P/24)` per unit. Equivalent
  to subtracting `P/24` nats of readable capacity per unit; the multiplicative and
  subtractive readings collapse to one branch since
  `e^(P/4) * e^(-P/24) = e^(P/4 - P/24)`.
- **BR-1b (presence):** survival factor `1 - P/24` per unit. The presence reading of
  the same reserve trace: `epsilon` is a conditional trace of the Z6 projection, a
  presence probability, and the receipt-consistent survival is `1 - P/24`
  (`extra/chi_nu_collar_survival_presence_correction.md`; Lean module
  `CollarGatePresence.lean`; the corpus carries both readings as of 2026-07-14).
- **BR-1c (none):** reserve slots hold repair records that remain countable states;
  no exclusion from the count.

**Branch point BR-2, reserve attachment (menu size 2).**

- **BR-2a (per cell):** one reserve unit per equal-area cell, `K` units.
- **BR-2b (per shared edge):** the reserve is priced on shared cuts; a locally
  six-valent cellulation (`thm:icosahedral-screen-sieve` background hypothesis) has
  `3K` shared edges, so `3K` units.

Resulting linear coefficient `rho` in `log|Omega^sc_N| = rho*N + ...`:

| rho branch | semantics x attachment | rho |
|---|---|---|
| R1 | poisson x per-cell | `1 - (4/P)(P/24) = 5/6` (exact) |
| R2 | presence x per-cell | `1 + (4/P) ln(1 - P/24) ≈ 0.8273998401` |
| R3 | none (attachment collapses) | `1` |
| R4 | poisson x per-edge | `1 - (12/P)(P/24) = 1/2` (exact) |
| R5 | presence x per-edge | `1 + (12/P) ln(1 - P/24) ≈ 0.4821995204` |

Menu size of the rho axis after collapse: 5.

### Step 4. Self-reading constraint: the 12-port load. Branch axis BR-3

Membership clause 4: the configuration's own horizon record surface reads back
capacity `N`. The declared readout surface is the twelve-port screen sieve
(`thm:icosahedral-screen-sieve`): the invariant screen load is `X = log(N/pi)`,
read locally as `X/12` per port; the oriented 24-slot repair register
`R_24 = P_12 x {+,-}` is bookkeeping of the same surface and creates no independent
carrier (`def:oriented-24-slot-register`).

**Branch point BR-3, count effect of the readback record (menu size 4).**

- **X+ (port multiplicity):** the port orbit is defect structure over and above the
  `e^N` bulk; filling it consistently at load `X` multiplies the count by
  `e^X = N/pi`. Contribution `(c,d) += (+1, -ln pi)`.
- **X0 (delta constraint):** the readback clause is a membership constraint on defect
  data outside the bulk count; no multiplicity, no cost. `(0, 0)`.
- **X- (record cost):** the bulk must dedicate `X` nats to store the readback record;
  count factor `e^(-X) = pi/N`. `(-1, +ln pi)`.
- **X= (write+check cost):** the oriented register requires write and check copies,
  cost `2X`. `(-2, +2 ln pi)`.

A fifth reading exists: equate the per-port read to the inner electromagnetic
observation step of the D10 lane, which forces
`log(N/pi)/12 = pi/(2 P alpha_U)` and a constant readback at
`pi*exp(6 pi/(P alpha_U))`. This is the CL-3 electroweak-bridge expression, listed
by name in the spec blindness cone (Section 3 of the spec). Excluded pre-evaluation
under V-08; recorded as branch CAP-B in Step 7 and never executed.

### Step 5. Observer/checkpoint subfederation. Branch axis BR-4

Membership clause 2: support at least one stable observer/checkpoint subfederation.

**Branch point BR-4 (menu size 3).**

- **K0 (unmarked):** an at-least-one property absorbs into the bulk count for large
  `K`; contribution `(0, 0)`.
- **K1 (one marked chain):** one marked host cell among `K = 4N/P`; factor `K`.
  Contribution `(+1, +ln(4/P))`.
- **K2 (chain + checkpoint):** marked host plus marked checkpoint cell; factor `K^2`
  (leading order). Contribution `(+2, +2 ln(4/P))`.

Readings with one chain per port (k = 12) or per face-corner flag (k = 60) import
structure the membership clause does not declare; recorded, not executed.

### Step 6. Symmetry quotient. Branch axis BR-5

The no-marked-point rule makes the twelve ports indistinguishable; the face-corner
bundle supplies twenty faces `F_20 ~ A5/C3` and sixty face-corner flags forming a
free transitive A5-set (`thm:icosahedral-face-corner-bundle` and its boundary
remark).

**Branch point BR-5 (menu size 3).**

- **S- (A5 quotient):** configurations counted up to the icosahedral action; divide
  by `|A5| = 60`; `d += -ln 60`.
- **S0 (already quotiented):** `nf` already delivers one representative per class;
  no factor.
- **S+ (flag multiplicity):** the sixty flags label distinct record anchorings;
  `d += +ln 60`.

### Step 7. Assembled count and candidate families for Cap_read. Branch axis BR-6

Steps 2 through 6 give, per branch,

```
log|Omega^sc_N| = rho*N + c*log N + d,
c = c_BR3 + k_BR4  in  {-2, ..., +3},
d = d_BR3 + k_BR4*ln(4/P) + d_BR5.
```

No declared structure produced an `N log N` term or any other tail: every declared
multiplicity and cost in Steps 3 through 6 is per-cell (linear in `N`), per-load
(`log N`), or global (constant). This is recorded now since it fixes the reachable
fixed-point scale before any evaluation.

The spec does not supply a formula for `Cap_read`; it supplies the growth form of P3
(a declared normalization times the log of the sector count) and the port readout
surface. **Branch point BR-6, Cap_read family (menu size 4).**

- **CAP-L (log-count readback):** `F(N) = 1 * log|Omega^sc_N| = rho*N + c*log N + d`.
  Normalization 1 is the only non-tunable choice: both sides are nats of `log dim`
  (spec Section 1.1). Sub-branches: the full (rho, c, d) lattice,
  5 x 4 x 3 x 3 = 180 rows.
- **CAP-P (port-inversion readback):** the sector reconstructs `N` by inverting the
  port read: `F(N) = pi * exp(12 * x_read(N))` with `x_read` the surviving per-port
  load. Sub-menu (6): survival multiplicative on the load per port
  (`x = s*X/12`, `s in {e^(-P/24), 1-P/24}`), per oriented slot pair
  (`s in {e^(-P/12), (1-P/24)^2}`), or additive reserve cost per port
  (`x = X/12 - P/24`) or per oriented pair (`x = X/12 - P/12`).
- **CAP-K (cell-count readback):** the sector reconstructs the chart count and reads
  `F(N) = (P/4) * K_readable = s*N`, `s in {e^(-P/24), 1-P/24, 5/6, 1/2}` (menu 4).
- **CAP-B (bridge readback):** the excluded constant reading of Step 4; barred by
  V-08 before evaluation.

### Step 8. P4 stationarity per branch

`l(N) = log|Omega^sc_N| - N = (rho-1)*N + c*log N + d`, `l'(N) = (rho-1) + c/N`.

- `rho < 1, c > 0`: unique interior stationary point `N_star = c/(1-rho)`, `l'' < 0`.
- `rho < 1, c <= 0`: `l' < 0` everywhere; the MAR argmax sits on the left boundary of
  the admissible interval; no input-free interior selector. P4 fails structurally.
- `rho = 1, c > 0`: `l' > 0` everywhere; argmax at the right boundary. P4 fails.
- `rho = 1, c < 0`: `l' < 0` everywhere. P4 fails.

So P4 admits an interior selector only on `rho < 1, c > 0` branches, where
`N_star = c/(1-rho)`. Largest value over the whole lattice: `c = 3`, `rho = 5/6`,
giving `N_star = 18` nats; smallest interior value `c = 1` on R5, `N_star ≈ 1.93`.
This is a pre-evaluation structural bound, recorded before any fixed point is
computed.

## 2. Execution registry (declared before comparison)

DECLARED-BEFORE-COMPARISON MARKER. Everything above this marker was written before
any branch was evaluated. The registry below lists every computed fixed point. The
comparison table appears only after the registry, per the anti-fitting discipline.
Implementations: [`F_candidate_capL.py`](F_candidate_capL.py),
[`F_candidate_capP.py`](F_candidate_capP.py),
[`F_candidate_capK.py`](F_candidate_capK.py); certificates in `runtime/`.

### 2.1 CAP-L lattice (180 rows)

Fixed-point condition `(1-rho)*N = c*log N + d`, stable root certified by the toy
machinery's interval Banach check (self-map enclosure plus derivative enclosure
`F' = rho + c/N`), `mpmath.iv`, 40 decimal digits, P as the certified enclosure.

Summary of the 180 rows (full per-row records in
`runtime/F_candidate_capL_certificates.json`):

180 rows executed. 98 rows carry a certified fixed point (interval Banach
certificate: centered mean-value self-map enclosure, `sup|F'| <= L < 1`, certified
enclosure of width `<= 2e-12`, most `2e-25`). 48 rows have no positive fixed point
(`F(N) < N` on all of `(0, inf)`). 18 rows are rejected with no contraction
(`rho = 1, c >= 0`: `F' >= 1` everywhere). 15 rows carry only an unstable fixed
point (`|F'| >= 1` at the root; fails P2/A4). 1 row is the identity (excluded by
P5). All certified fixed points lie in `[1.4686, 1452.33]` nats.

| rho branch | c | rows | certified | no fixed point | rejected/unstable | certified N_CRC range (nats) |
|---|---|---|---|---|---|---|
| R1 (rho=5/6) | +3 | 3 | 3 | 0 | 0 | 49.6 .. 113.7 |
| R1 (rho=5/6) | +2 | 6 | 5 | 1 | 0 | 24.66 .. 89.23 |
| R1 (rho=5/6) | +1 | 9 | 5 | 4 | 0 | 24.6 .. 67.47 |
| R1 (rho=5/6) | +0 | 9 | 5 | 4 | 0 | 12.25 .. 49.07 |
| R1 (rho=5/6) | -1 | 6 | 4 | 0 | 2 | 2.183 .. 24.49 |
| R1 (rho=5/6) | -2 | 3 | 2 | 0 | 1 | 2.542 .. 10.31 |
| R2 (rho=0.8274) | +3 | 3 | 3 | 0 | 0 | 46.94 .. 109 |
| R2 (rho=0.8274) | +2 | 6 | 5 | 1 | 0 | 23.01 .. 85.69 |
| R2 (rho=0.8274) | +1 | 9 | 5 | 4 | 0 | 23.48 .. 64.93 |
| R2 (rho=0.8274) | +0 | 9 | 5 | 4 | 0 | 11.83 .. 47.38 |
| R2 (rho=0.8274) | -1 | 6 | 4 | 0 | 2 | 2.163 .. 23.82 |
| R2 (rho=0.8274) | -2 | 3 | 2 | 0 | 1 | 2.526 .. 10.14 |
| R3 (rho=1) | +3..0 | 27 | 0 | 8 | 19 | none |
| R3 (rho=1) | -1 | 6 | 4 | 0 | 2 | 3.142 .. 1452 |
| R3 (rho=1) | -2 | 3 | 2 | 0 | 1 | 3.142 .. 24.33 |
| R4 (rho=1/2) | +3 | 3 | 2 | 1 | 0 | 18.95 .. 29.87 |
| R4 (rho=1/2) | +2 | 6 | 4 | 2 | 0 | 7.637 .. 24.59 |
| R4 (rho=1/2) | +1 | 9 | 5 | 4 | 0 | 5.022 .. 20.06 |
| R4 (rho=1/2) | +0 | 9 | 5 | 4 | 0 | 4.084 .. 16.36 |
| R4 (rho=1/2) | -1 | 6 | 4 | 0 | 2 | 1.491 .. 9.964 |
| R4 (rho=1/2) | -2 | 3 | 2 | 0 | 1 | 1.936 .. 5.762 |
| R5 (rho=0.4822) | +3 | 3 | 2 | 1 | 0 | 18 .. 28.59 |
| R5 (rho=0.4822) | +2 | 6 | 4 | 2 | 0 | 7.084 .. 23.58 |
| R5 (rho=0.4822) | +1 | 9 | 5 | 4 | 0 | 4.736 .. 19.3 |
| R5 (rho=0.4822) | +0 | 9 | 5 | 4 | 0 | 3.943 .. 15.79 |
| R5 (rho=0.4822) | -1 | 6 | 4 | 0 | 2 | 1.469 .. 9.678 |
| R5 (rho=0.4822) | -2 | 3 | 2 | 0 | 1 | 1.914 .. 5.644 |

Closed-form landmarks inside the registry: `capL.R3.Xm.K0.S0` has
`F(N) = N - log N + log pi` with certified fixed point exactly `pi`;
`capL.R3.Xq.K1.Sp` tops the registry at `pi^2 * 240/P ≈ 1452.327` with
`L ≈ 0.99934` (the slowest certified contraction of the run).

P4 record: `p4_coherent_rows = 0`. Every certified row with an interior selector
registers a discrepancy `N_CRC != N_star` (both O(1)-O(10^2) nats); every `c <= 0`
row fails interior stationarity structurally. No branch satisfies P4 as coherence;
all would enter the ledger as registered-discrepancy candidates even before the
basin comparison.

### 2.2 CAP-P family (6 rows)

| branch | reading | map | status | N_CRC (nats) |
|---|---|---|---|---|
| capP.s_poisson_port | `s = e^(-P/24) ≈ 0.9343005` | `pi (N/pi)^s` | certified | `pi` (exact) |
| capP.s_presence_port | `s = 1 - P/24 ≈ 0.9320428` | `pi (N/pi)^s` | certified | `pi` (exact) |
| capP.s_poisson_pair | `s = e^(-P/12)` | `pi (N/pi)^s` | certified | `pi` (exact) |
| capP.s_presence_pair | `s = (1 - P/24)^2` | `pi (N/pi)^s` | certified | `pi` (exact) |
| capP.add_slot | `x = X/12 - P/24` | `N e^(-P/2)` | no positive fixed point | none |
| capP.add_port | `x = X/12 - P/12` | `N e^(-P)` | no positive fixed point | none |

Every multiplicative degradation of the per-port load `X/12 = log(N/pi)/12`
collapses the readback fixed point to `N = pi`: the unique solution of
`(N/pi)^s = N/pi` with `s != 1` is `N = pi` regardless of `s`. The port-inversion
family cannot place the fixed point anywhere else; the reserve strength only sets
the contraction rate. P4: registered discrepancy on every certified row
(`N_CRC = pi` differs from every interior `N_star` of the lattice, range
`[1.93, 18]`).

### 2.3 CAP-K family (4 rows)

| branch | readable fraction s | map | status |
|---|---|---|---|
| capK.s_poisson | `e^(-P/24)` | `F(N) = s N` | no positive fixed point |
| capK.s_presence | `1 - P/24` | `F(N) = s N` | no positive fixed point |
| capK.s_nat_share | `5/6` | `F(N) = s N` | no positive fixed point |
| capK.s_edge_share | `1/2` | `F(N) = s N` | no positive fixed point |

The cell-count readback is linear through the origin with `s < 1` on every
reserve reading: a contraction whose unique fixed point is `N = 0`, outside the
admissible interval. The whole family is excluded; the P3 bracketing pair
`F(a) >= a` does not exist for any `a > 0`.

## 3. Comparison table (written after Section 2 was frozen)

Reference capacities, used here for the first time in this document and never inside
any candidate cone: the Lambda-located SL-4 basin `3.31e122` and the CL-3
electroweak-bridge value `3.53e122` (`log` values `282.06` and `282.13` nats of
`log N`; as capacities, `~3.3e122` nats).

Comparison record: `runtime/F_construction_comparison_2026-07-14.json`
([`compare_references.py`](compare_references.py), the only file of the run whose
cone contains the reference values). Landing criterion
`|log10(N_fp / N_ref)| < 1` against either reference.

| family | rows | certified fixed points | N_CRC range (nats) | log10 shortfall to 3.31e122 | landed |
|---|---|---|---|---|---|
| CAP-L | 180 | 98 | 1.4686 .. 1452.33 | 119.36 .. 122.35 | 0 |
| CAP-P | 6 | 4 | pi | 122.02 | 0 |
| CAP-K | 4 | 0 | none | n/a | 0 |
| total | 190 | 102 | 1.4686 .. 1452.33 | >= 119.36 | 0 |

No branch lands. The best certified fixed point of the whole run
(`capL.R3.Xq.K1.Sp`, `N_CRC ≈ 1452.33` nats) sits `10^119.36` below the SL-4
basin and `10^119.39` below the CL-3 bridge value. Every certified branch is
excluded as a candidate for the physical readback map at the reference scale.

## 4. Outcome

**Outcome class (b): no branch lands; CL-7 stays open; the constructed candidates
are recorded as excluded.** The exclusion is informative in three ways.

1. **Scale exclusion.** Every count assembled from the declared 12/24-port,
   Z6-reserve, and P/4-budget combinatorics has the form
   `log|Omega^sc_N| = rho*N + c*log N + d` with `rho` bounded away from 1 by the
   reserve fraction (`1 - rho in {1/6, 0.1726, 1/2, 0.5178}`) and `|c| <= 3`,
   `|d| <= 7.3`. Both displays of the target then sit at
   `N_star = c/(1-rho) <= 18` and `N_CRC <= 1452.33` nats. A reference-scale fixed
   point `~3.3e122` requires either `1 - rho ~ 1e-122` or a structure-supplied
   tail (an `N log N` term with a stationary point near `exp(282)`), and the
   declared structure supplies neither. This holds independently of every branch
   choice in the run.
2. **P4 exclusion.** Zero rows are P4-coherent: the readback fixed point and the
   count-density stationary point disagree on every certified branch. The two
   representations of the target separate as soon as `Cap_read` is anything other
   than the exact log-count at the stationary point, and no declared normalization
   makes them meet on this lattice.
3. **Fixed-point rigidity of the port inversion.** The port-inversion family pins
   `N_CRC = pi` for every degradation strength: the readback scale is set by the
   `pi` in the invariant load `X = log(N/pi)`, and reserve costs cannot move it.
   Any port-based landing at `~1e122` must inject `~282` nats of additive per-load
   structure, which is Step 4's excluded fifth reading (the CL-3 bridge exponent).

**Gap item G2-GAP-1 (new).** The declared structure underdetermines the count at
exactly one load-bearing point: no theorem supplies a per-cell or per-port record
multiplicity of order `exp(6*pi/(P*alpha_U))`, equivalently a mechanism producing
`1 - rho` at the `1e-122` scale or an `N log N` tail with stationary point near
`exp(282)`. The only corpus object of that order is the CL-3 electroweak-bridge
exponent, whose value is barred from the evaluation cone by the spec (Section 3)
and whose derivation from the source side is precisely the CL-3 resolution path
("derive the EW-bridge/capacity connection"). Until a source-side theorem ties the
port load `X = log(N/pi)` to the inner observation step of the D10 lane, every
count built from the declared combinatorics closes at O(1)-O(10^3) nats. G2 (this
map) and the CL-3 connection are one missing theorem, seen from two rows.

**Ledger action: none.** No branch passes A7, so
[CLOSURE_LEDGER.md](../../CLOSURE_LEDGER.md) row CL-7 is untouched; this document
and the runtime certificates carry the exclusion record. Promotion of any future
candidate requires: a count with a structure-derived reference-scale stationary
point, P4 coherence or a registered discrepancy row, a directed-rounding backend
(Arb/MPFI class) for the certificates, and a blind A7 landing.

Wall-clock for the full run (three families, 190 rows, comparison):
about 80 seconds on the development machine; each family script is deterministic
(byte-identical artifacts across reruns, checked by SHA-256).
