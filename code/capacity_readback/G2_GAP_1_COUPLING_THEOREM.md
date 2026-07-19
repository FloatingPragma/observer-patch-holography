# EW Capacity Identification Reduction

Conditional algebraic reduction for coupling the twelve-port screen load
`X = log(N/pi)` to the D10 observation step. It compares an independently
constructed public stable-record fixed point with the electroweak bridge only
after the physical load receipts hold. It does not construct the capacity map.
The direct producer contract is [F_READBACK_SPEC.md](F_READBACK_SPEC.md).

**Status: conditional electroweak identification theorem. The displayed fixed
point belongs to a target-tainted affine candidate and has no CL-7 promotion
weight. Its exact algebra survives as a comparison after an independent direct
public-record `F` exists.**

The assumption-removal split is:

- the exterior theorem discharges the mathematical multiplicity `4`; identifying
  that multiplicity with the physical coefficient `beta_EW` still requires
  `PORT-LOAD-TRACE` and the declared D10 normalization;
- load additivity is discharged by unitary invariance, direct-sum additivity,
  and `L_P(C)=P`;
- `PORT-LOAD-TRACE` must identify that additive load with the physical port
  readout;
- no `24/6=4` screen inference is admitted because a physical free `Z6` action
  on the bookkeeping slots has not been constructed;
- the independent public-record equalizer, stable code, whole-fiber
  scalarization, capacity-extension naturality, greatest fixed point, and
  refinement receipts remain open;
- an averaging carrier affects dynamics and fluctuations. It is not required
  to derive the algebraic fixed-point location once the two physical readouts
  are equal.

## 1. Statement

**Theorem (target-coupled affine candidate; conditional on CP-1, CP-2, CP-3).**
Let `P` be the certified forward pixel-closure fixed point and let
`alpha_U(P)` be the unified gauge width fixed by the D10 pixel-closure equation
`ellbar_SU2(t2) + ellbar_SU3(t3) = P/4` on the realized product-gauge branch. Under
CP-1 through CP-3, the target-coupled candidate equation `N = F_EW(N)`, with
the twelve-port screen load `X = log(N/pi)` read through the D10 observation
step, has a unique fixed point on the admissible
interval, and that fixed point is exactly the electroweak bridge capacity

```
N_CRC = pi * exp( 6*pi / (P * alpha_U(P)) ).
```

Equivalently, the tick-projection bridge balances at the transmutation count,

```
Pi_EW(P, N_CRC) = beta_EW * P,        Pi_EW(P,N) = 24*pi / (alpha_U(P) * log(N/pi)),
beta_EW = N_c + 1 = 4,
```

and the bridge residual vanishes, `B_EW(P, N_CRC) = alpha_U*log(N_CRC/pi) - 6*pi/P = 0`.

At the certified forward enclosure `P_fwd` (display `1.630972095858897`) with the
certified `alpha_U(P_fwd)` enclosure (display `0.041124247441816685`), the
candidate fixed-point enclosure is

```
N_CRC in [ 3.532131543418935831227120163546618376916e122,
           3.532131543418935831227120721623402237107e122 ]      (relative width 1.58e-25).
```

The theorem is parametric in the branch pair `(P, alpha_U(P))`: on the public endpoint
branch (`P_star = 1.6309682094...`, `alpha_U = 0.041124336...`) the same closed form
gives `3.5323546226...e122`, the value recorded by
`R_EW_global_capacity_certificate.json`.

## 2. Declared sources

| object | source |
|---|---|
| direct public-record map, closure `N = log M_0(U_N)`, count kernel | [F_READBACK_SPEC.md](F_READBACK_SPEC.md); `correctable_public_record_capacity.py` |
| twelve-port sieve, load `X = log(N/pi)` read as `X/12` | `paper/screen_microphysics_and_observer_synchronization.tex`, `thm:icosahedral-screen-sieve`; oriented register `def:oriented-24-slot-register` |
| D6 radius identity `r_CRC/ell_star = (N/pi)^(1/2)` | `code/particles/hierarchy/certificates/R_N_global_repair_tick_certificate.json`, `proved_by_certificate` |
| D10 pixel closure, `alpha_U(P)`, `t_tr = 2*pi/(beta_EW*alpha_U)`, `beta_EW = N_c+1 = 4` | `code/particles/runs/calibration/d10_ew_forward_transmutation_certificate.json`; `code/P_derivation/paper_math.py` (pixel residual, line 588); `code/P_derivation/FULL_DERIVATION.md`; `extra/compact_proof_of_oph.tex` |
| repair tick `-log|g_*'| = log(N/pi)/48`, round count `m_rep = 24` | `R_N_global_repair_tick_certificate.json`; `R_m_rep_24_certificate.json` |
| tick projection `Pi_EW = 24*pi/(alpha_U*X)` and bridge algebra | `R_EW_tick_projection_certificate.json` (derivation chain steps 1 to 3, 5) |
| EW share of screen depth `Gamma_EW = (P/12)*log(N/pi)` | `extra/compact_proof_of_oph.tex`, QCD-free hierarchy witness |
| contraction carrier, `lambda = 1/2` free | `R_EW_global_capacity_certificate.json`, derivation chain step 7 |
| seed fixed point `pi` of the pure port inversion | [F_CONSTRUCTION_2026-07-14.md](F_CONSTRUCTION_2026-07-14.md) Section 2.2; [F_candidate_capP.py](F_candidate_capP.py) |
| `beta_EW = 4` relation | conditional D10 relation stated in the compact paper; no registered comparison target |
| certified `P_fwd` and `alpha_U(P_fwd)` enclosures | `code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json` (`interval_diagnostics.alpha_u`) |

## 3. Notation

All capacities in nats, `N = log dim H_partial,N`. Load `X = X(N) = log(N/pi)`.
`t_tr = 2*pi/(beta_EW*alpha_U(P))` is the D10 dimensional-transmutation depth,
`-log|g_*'| = X/48` is the per-tick screen contraction, `ell_cell = P/4` is the
per-cell cut density, `x_EW(P) = 6*pi/(P*alpha_U(P))` is the D10-side load.

## 4. Derivation, numbered steps

Each step carries the label PROVEN (follows from a cited declared theorem by the stated
reasoning) or PREMISE (declared here, carried in Section 5).

**S1 (PROVEN, definitional). Object.** The direct contract first constructs the
stable public record code size `M_pub(q)` on each terminal quotient, returns the
set of values across the complete terminal fiber, and permits a scalar `F` only
after whole-fiber agreement. Membership in that fiber does not require
`M_pub=D`. The closure equation is imposed only after this independent
scalarization. Source: [F_READBACK_SPEC.md](F_READBACK_SPEC.md).

**S2 (CONDITIONAL THEOREM). Twelve-port read.** The invariant screen load `X = log(N/pi)` is read
locally as `X/12` across twelve indistinguishable ports
(`thm:icosahedral-screen-sieve`). The oriented 24-slot register is bookkeeping of the
same surface and creates no independent carrier (`def:oriented-24-slot-register`).
Inherits IH-1a (strict unit splitting), IH-1b (inverse pairing and D-optimal
tomography), and IH-2 (`l_shared = P/4` branch input).

**S3 (PROVEN). Seed and load coordinate.** The D6 radius identity
`r_CRC/ell_star = (N/pi)^(1/2)` is proved by the global repair-tick certificate, so
`X = 2*log(r/ell_star)` and `X = 0` exactly at `N = pi`. The zero-coupling
port-inversion candidate has certified fixed point exactly `pi` on all four degradation
branches of the CAP-P run. The seed-times-inflation form `N = pi*e^X` is the radius
identity itself; composition at the level of the load coordinate is therefore licensed
by declared structure. The stronger claim that `Cap_read` reconstructs capacity through
this inversion is CP-2 (step S9).

**S4 (CONDITIONAL D10 STEP).** The pixel-closure equation
`ellbar_SU2 + ellbar_SU3 = P/4` fixes `alpha_U(P)`; the forward transmutation law then
fixes the transmutation depth `t_tr = 2*pi/(beta_EW*alpha_U(P))`. The exterior
certificate gives four weak-doublet copies per generation. Setting the physical
D10 coefficient `beta_EW` equal to that multiplicity requires `PORT-LOAD-TRACE`
and the declared normalization. This gives `beta_EW=4` on that conditional branch
(`d10_ew_forward_transmutation_certificate.json`, theorem block; compact proof).
Inherits IH-3 (declared D10 branch conventions).

**S5 (PROVEN ON THE PRODUCT-ADJOINT BRANCH). Global repair tick.**
`|g_*'| = (N/pi)^(-1/48)`, so the per-tick
contraction is `-log|g_*'| = X/48`, with `48 = 2*m_rep` and
`m_rep = 2*dim(su(3)+su(2)+u(1)) = 24` derived by the representation-to-spectrum
round-count theorem. This derivation does not use the twelve-port screen, the global
quotient, or the MAR matter realization (`R_m_rep_24_certificate.json`). The equalities
`12 = 8+3+1` and `24 = 2(8+3+1)` therefore do not identify the sieve ports or their
oriented bookkeeping register with gauge generators or the product-adjoint repair
register. Inherits IH-4 (repair-tick declared items).

**S6 (CONDITIONAL ON HIERARCHY-SCREEN-READOUT). Identity I1, the tick-projection
formula.** The bridge premise permits the screen load and product-adjoint tick law to
enter one readout. The tick projection is the ratio of the transmutation depth to the
per-tick contraction:

```
Pi_EW(P,N) := t_tr / (X/48) = (2*pi/(beta_EW*alpha_U)) * (48/X)
            = 4*pi*m_rep/(beta_EW*alpha_U*X) = 24*pi/(alpha_U*X).
```

The coefficient arithmetic is exact: `(2/4)*48 = 24`. This reproduces derivation-chain
step 3 of `R_EW_tick_projection_certificate.json`, discharged there. Origin note for
the factor `pi`: the `pi` in `24*pi` is the one-loop `2*pi` of the transmutation
exponent `t_tr`, carried through the ratio; the numerator decomposes as
`4*pi*m_rep/beta_EW`. No declared theorem assigns a per-slot phase of `pi` to the
24-slot register; the candidate decomposition "24 slots times `pi` per slot" appears
nowhere in the sources and is unused here.

**S7 (PREMISE CP-1). Identity I2, the balance condition.** The coupling asserts

```
Pi_EW(P, N_CRC) = beta_EW * P.
```

Three displays of the same condition, equivalent by exact algebra (S8):

1. tick balance: the projected transmutation depth equals one global repair tick per
   electroweak transmutation channel per local pixel-area unit,
   `Pi_EW = beta_EW * P`;
2. depth balance: the electroweak share of screen depth equals the transmutation
   depth, `Gamma_EW = t_tr`, where
   `Gamma_EW = beta_EW * ell_cell * X/12 = (P/12)*X` (compact proof, QCD-free
   hierarchy witness);
3. per-port coupling: the per-port load equals the inner observation step of the D10
   lane, `X/12 = pi/(2*P*alpha_U)` (the fifth reading of construction Step 4).

The corpus consumes this condition as the declared resonance target; the
tick-projection certificate scope note states that a deeper geometric derivation of the
product target belongs to a separate strengthening theorem. It is carried here as
CP-1. Structural evidence, recorded and distinct from a derivation: in the
preregistered sweep, all Standard Model gauge triples fail structurally (no
pixel-residual root exists on any window, 20 diagnosed chains), and `beta_EW = 3`
chains close more than three orders of magnitude below the weak scale, while
`beta_EW = 4` chains reproduce the weak-scale pair.

**S8 (PROVEN). Equivalence arithmetic.** Exact chain:

```
Pi_EW = 4P  <=>  24*pi/(alpha_U*X) = 4P  <=>  alpha_U*X = 6*pi/P
        <=>  X = 6*pi/(P*alpha_U)  <=>  N = pi*exp(6*pi/(P*alpha_U)).
```

Coefficients exact (`24/4 = 6`); the `Gamma_EW` and per-port displays reduce to the
same exponent (`12*2/beta_EW = 6`, `12/2 = 6`). Numerical witness at 60 digits, at the
certified `P_fwd` and `alpha_U(P_fwd)` (Section 8): all four residual displays vanish
to below `5e-60`.

**S9 (PREMISE CP-2). Readback form.** `Cap_read` on the coupled branch returns the
port-load inversion, `F(N) = pi * exp(X_read(N))`. Licensed components: the nat units
on both sides under the archived candidate conventions and the seed identity (S3).
Premise component: the selection of the inversion family among `Cap_read`
candidates (branch axis BR-6).
Exclusion support: the recorded run certifies the CAP-L and CAP-K families excluded
at reference scale under every declared branch choice, leaving CAP-P as the only
executable family whose coupled extension reaches the bridge scale. Exclusion narrows
the menu; a derivation of the form is the discharge obligation.

**S10 (PREMISE CP-3). Re-emission dynamics.** The re-emitted load is the average

```
X_read(N) = (1 - lambda)*X(N) + lambda*x_EW(P),      lambda in (0,1),
```

with `lambda = 1/2` recorded as the free averaging weight
(`R_EW_global_capacity_certificate.json`, step 7 note). CP-3 is dispensable for the
location of the fixed point: at any fixed point of a CP-2 map, `X_read = X`, and CP-1
pins the balanced load at `X = x_EW` regardless of the dynamics. CP-3 supplies the
constructive contraction tested in the archived candidate checks A3 to A5. The
`lambda = 1` constant reading is the historically excluded CAP-B branch.

**S11 (PROVEN, conditional on CP-1..CP-3). Composition and Banach certificate.** Under
CP-1 through CP-3 the readback map is affine on the load coordinate,
`C(y) = (1-lambda)*y + lambda*x_EW`, with Lipschitz constant exactly
`1 - lambda = 1/2 < 1`. By the Banach fixed-point theorem the fixed point is unique and
equals `x_EW(P)`, independent of `lambda`; hence `N_CRC = pi*exp(x_EW(P))`, which is S8's
bridge capacity. Interval certificates: [F_candidate_coupled.py](F_candidate_coupled.py),
artifact `runtime/F_candidate_coupled_certificates.json`.

**S12 (record). Certified enclosure.** At the certified `P_fwd` enclosure and the
certified `alpha_U(P_fwd)` enclosure (both from the P interval contraction
certificate, `iv_dps = 60`):

```
x_EW  in [281.032552984795966015533961145, 281.032552984795966015533961261]
N_CRC in [3.532131543418935831227120163546618376916e122,
          3.532131543418935831227120721623402237107e122]
```

relative width `1.58e-25`, `mpmath.iv` outward rounding at 60 decimal digits.
Theorem-grade promotion requires a directed-rounding backend (Arb/MPFI class), matching
the spec's promotion gate.

## 5. Premise ledger

Premises declared by this theorem:

| id | assertion | what would discharge it |
|---|---|---|
| CP-1 | Balance condition `Pi_EW(P, N_CRC) = beta_EW*P`: readback self-consistency holds where the projected transmutation depth equals one global repair tick per electroweak channel per pixel-area unit (equivalently `Gamma_EW = t_tr`; equivalently `X/12 = pi/(2*P*alpha_U)`) | a geometric counting theorem for the product target `beta_EW*P` from declared screen structure, the strengthening theorem named in the tick-projection certificate scope note |
| CP-2 | `Cap_read` on the coupled branch is the port-load inversion `F(N) = pi*exp(X_read(N))` (seed and units proven, family selection premised) | a uniqueness theorem for the inversion form among `Cap_read` families. Reduction record: [CP2_INVERSION_FORM_REDUCTION_2026-07-17.md](CP2_INVERSION_FORM_REDUCTION_2026-07-17.md) decomposes CP-2 into load mediation plus read consistency, classifies the surviving freedom exactly (the load gauge `phi`, pinned at the seed and, under P4 plus CP-4, at `x_EW`), and proves a P4-coherence argument alone cannot pin more than the fixed-point location |
| CP-3 | The re-emitted load is the `lambda`-average of the screen-side and D10-side loads, `lambda in (0,1)` (dispensable for the fixed-point location, load-bearing for the constructive contraction) | a derivation of the averaging carrier, for instance from the write/check orientation split; any declared contraction toward balance suffices |

Inherited named hypotheses of cited theorems (declared in their sources rather
than independently here):

| id | hypothesis | source |
|---|---|---|
| IH-1a | feasible integer charge domain and strict unit-splitting cost | `thm:icosahedral-screen-sieve` |
| IH-1b | inverse-port pairing and source-side D-optimal vector/quadrupole selector | `thm:icosahedral-screen-sieve` |
| IH-2 | shared-cut entropy density `l_shared = P/4` as declared branch input | preamble of `thm:z6-reserve-trace` |
| IH-3 | D10 branch conventions: pixel-closure form, one-loop MSSM coefficients, tree-level `m_Z` closure, Stage-5 continuation (certified as declared) | `p_interval_contraction_certificate_2026-07-14.json` claim boundary |
| IH-4 | repair-tick declared items: homogeneous positive-root one-tick normal form; the archived `pi/rho_read^2` model is retained only as an independent geometric estimator | `R_N_global_repair_tick_certificate.json`, `declared_not_derived`; [F_READBACK_SPEC.md](F_READBACK_SPEC.md) |

## 6. Certified execution and spec compliance

Executed candidate: `F(N) = pi*exp((1-lambda)*log(N/pi) + lambda*x_EW(P))`,
`lambda = 1/2`, [F_candidate_coupled.py](F_candidate_coupled.py).

| test | result |
|---|---|
| A1 finite executability | pass; deterministic agreement across reruns |
| A2 refinement stability | executed as inapplicable: no finite-cutoff `F_r` family is constructed for the coupled map; the obligation travels with CP-2/CP-3 |
| A3 self-map enclosure | pass on the load interval `[x_EW - 1, x_EW + 1]`; capacity-coordinate check passes on `[0.9, 1.1] * N_CRC` |
| A4 contraction | pass; derivative enclosure thin at `1/2` on the load coordinate (`L = 0.5827` on the capacity-coordinate interval) |
| A5 fixed-point enclosure | pass; enclosures of S12, `x_EW` contained in the certified box |
| A6 blindness (V-08) | fails by construction, recorded before any comparison: the candidate definition contains the bridge exponent `6*pi/(P*alpha_U)`. It cannot serve as a blind CL-7 landing test |
| A7 landing | informative comparison only, no verdict: fixed point over the SL-4 display `3.31e122` is `1.0671` (`log10` offset `0.0282`). This restates the CL-3 residual (bridge capacity above the Lambda-located capacity by 6.6 percent in the ledger's units); it is the single live test if CP-1 through CP-3 discharge |

Properties: P1 pass on the certificate interval (physical-grammar totality travels with
CP-2/CP-3); P2 pass (`F' > 0`); P3 pass (`L = 1/2`); P4 open obligation, recorded (the
count representation of the coupled membership clause is unconstructed, no coherence
claim is made; the reduction of the obligation to premise CP-4 under the
per-configuration reading of membership clause 4, and its degenerate closure under the
aggregate reading, are recorded in
[P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md](P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md));
P5 pass for `lambda in (0,1)`.

## 7. Consequence for the ledger geometry

CP-1 through CP-3 close the algebra and contraction analysis of the displayed
target-coupled affine candidate. They do not close CL-7 because its evaluation
cone contains the electroweak bridge exponent. CL-7 requires the independent
public-record equalizer and stable code, complete-fiber scalarization,
capacity-extension naturality, refinement, and greatest-fixed-point receipts
in the primary specification.
Once that independent fixed point exists, this record supplies a conditional
comparison: the physical load receipts identify it with the electroweak value,
while a mismatch rejects the bridge identification without changing `F`.

## 8. Verification record

Exact coefficient checks (rational arithmetic): `(2/beta_EW)*2*m_rep = 24`,
`24/beta_EW = 6`, `12*2/beta_EW = 6`, `12/2 = 6`, all at `beta_EW = 4`, `m_rep = 24`.

Numerical witness, `mpmath` at 60 significant digits, at the certified enclosure
midpoints `P_fwd = 1.63097209585889737696451390350695562985390`,
`alpha_U(P_fwd) = 0.041124247441816685140889933889659717292128290516...`:

```
N            = 3.5321315434189358312271204430193741e122
X            = 281.03255298479596601553396120312297
Pi_EW        = 6.5238883834355895078580556140278225
beta_EW * P  = 6.5238883834355895078580556140278225
Pi_EW - 4P                =  -6.2e-61
X/12 - pi/(2*P*alpha_U)   =   0.0
(P/12)*X - t_tr           =   5.0e-60
alpha_U*X - 6*pi/P        =   1.2e-60
```

## 9. Canon actions

This record closes no ledger row by itself. CL-7 stays open until a source-only
direct-capacity receipt satisfies [F_READBACK_SPEC.md](F_READBACK_SPEC.md). The
target-coupled execution remains a comparison and negative-control artifact.
