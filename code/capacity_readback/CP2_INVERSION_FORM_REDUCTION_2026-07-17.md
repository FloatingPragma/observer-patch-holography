# CP-2 inversion-form reduction

> **Archived target-coupled reduction.** The canonical capacity producer is
> `F_hat_r(D)=M_pub(D)` in [F_READBACK_SPEC.md](F_READBACK_SPEC.md). Neither
> port-load inversion nor `pi/rho_op^2` defines `F`. This record is retained as
> a negative-control and electroweak comparison analysis only.

Companion to [G2_GAP_1_COUPLING_THEOREM.md](G2_GAP_1_COUPLING_THEOREM.md) and
[F_READBACK_SPEC.md](F_READBACK_SPEC.md). Target: GAP-A3, the discharge
obligation of premise CP-2 (the port-load inversion form of `Cap_read` on the
coupled branch). Machine checks:
[cp2_p4_premise_reduction.py](cp2_p4_premise_reduction.py), artifact
`runtime/cp2_p4_premise_reduction_check.json`.

**Status: reduction record. CP-2 stays open. The premise decomposes into two
named components (LM, RC); the surviving freedom without RC is classified
exactly (the load gauge); the fixed-point location is forced across the whole
freedom under P4 coherence plus one named count-side premise. No ledger row
moves.**

## 1. Setting

CP-2 asserts that `Cap_read` on the coupled branch returns the port-load
inversion, `F(N) = pi * exp(X_read(N))`, with `X = log(N/pi)` the invariant
screen load. The proven ingredients (coupling theorem S3): the load chart
`X: (0, inf) -> R`, `X(N) = log(N/pi)`, is a strictly increasing bijection
supplied by the D6 radius identity `r_CRC/ell_star = (N/pi)^(1/2)`, with
inverse `X^(-1)(x) = pi * e^x`, and the zero-coupling readback has certified
fixed point exactly `pi`. The premise component named by S9 is the selection
of the inversion family among `Cap_read` candidates (branch axis BR-6).

## 2. Decomposition of CP-2

Two structural conditions, stated separately so that neither hides the other:

- **LM (load mediation).** On the coupled branch, `Cap_read` consumes exactly
  the re-emitted invariant screen load: `F(N) = G(X_read(N))` for some
  strictly increasing `G` on the read range, with no sector datum other than
  the load entering `G`. Support: membership clause 4 reads through the
  horizon record surface whose declared invariant is the load
  (`thm:icosahedral-screen-sieve`); the two constructed non-load families of
  the 2026-07-14 run (CAP-L log-count, CAP-K cell-count) close at or below
  `10^3` nats or have no fixed point. LM is a statement about the
  construction, not about the function `N -> F(N)` alone: since `X_read` is
  injective for `lambda < 1`, any map factors set-theoretically; LM says the
  factorization is the mechanism.
- **RC (read consistency).** The reconstructed capacity carries the load that
  was read: `X(G(x)) = x` on the read range.

**Lemma (gauge parameterization; exact).** Under LM, `phi := X o G` is
strictly increasing and `F(N) = pi * exp(phi(X_read(N)))`. The correspondence
`G <-> phi` is a bijection between LM candidates and strictly increasing load
gauges `phi`. RC holds iff `phi = id`. Proof: `X` is a proven bijection.

**Theorem CP2-A (uniqueness under LM + RC).** Any candidate satisfying LM and
RC equals the port-load inversion `F(N) = pi * exp(X_read(N))`. Proof: a right
inverse of a bijection is the inverse, and `X^(-1)(x) = pi * e^x` by the chart
theorem. Reading rule for this theorem: RC sits one proven-chart step from
the conclusion, so CP2-A relocates the content of CP-2 into RC; it is a
decomposition, not a discharge.

## 3. Sharp independence result

**Theorem CP2-B (the surviving freedom is exactly the load gauge).** In the
class defined by LM, the spec axioms P1, P2/P3 (monotone contraction on a
certificate interval), P5, and the seed theorem (the decoupled member fixes
`pi`, equivalently `phi(0) = 0`), the candidate set is exactly

```
{ F_phi = pi * exp(phi o X_read) :  phi strictly increasing, phi(0) = 0,
  (1 - lambda) * sup phi' <= kappa < 1 on the certificate interval }.
```

The inversion form is not unique in this class. Witnesses, machine-certified
(artifact, part A2/A3): the power gauges `phi_s(x) = s*x` for
`s in {1/2, e^(-P/24), 1 - P/24, 0.99}` (the CAP-P degradation branches of
the recorded run) and the affine-in-capacity gauges
`phi_a(x) = log(a*e^x + 1 - a)` for `a in {0.3, 0.7}` satisfy every axiom of
the class, violate RC with a certified residual, and place the coupled fixed
point off balance by certified separations between `0.71` and `187.4` nats of
load. The shift gauge `phi(x) = x + 1/10` is excluded by the seed
normalization (named principle), independently of RC.

## 4. Conditional forcing of the fixed-point location

**Theorem CP2-C (P4 pins the location, not the gauge).** Assume CP-1 (balance
load `x_EW = 6*pi/(P*alpha_U)`), CP-3 (re-emission
`X_read = (1-lambda)*X + lambda*x_EW`), the class of CP2-B, and

- **CP-4 (stationarity at balance).** The coupled count density
  `l(N) = log|Omega^sc_N| - N` is differentiable on the admissible interval
  and carries the derivative-sign certificate with its unique zero exactly at
  load balance `X(N) = x_EW`.

Then for a candidate `F_phi`, P4 coherence (MAR argmax of `l` equals the
fixed point of `F_phi`) holds iff `phi(x_EW) = x_EW`, and every P4-coherent
candidate has the same unique fixed point `N_CRC = pi * exp(x_EW)`, the
bridge capacity, independent of `phi`. Proof: CP-4 places the argmax at
`X^(-1)(x_EW)`; if `phi(x_EW) = x_EW` then `X = x_EW` solves
`X = phi((1-lambda)*X + lambda*x_EW)`, unique by the Banach bound
`(1-lambda) * sup phi' < 1`; conversely a fixed point at `x_EW` forces
`phi(x_EW) = x_EW` by direct substitution.

Machine witness (artifact, family `balance-pinned sine gauge`): the gauge
`phi(x) = x + (1/10)*(x_EW/pi)*sin(pi*x/x_EW)` differs from the identity by
about `8.9` nats at `x_EW/2`, satisfies every class axiom, pins
`phi(0) = 0` and `phi(x_EW) = x_EW` exactly, and its coupled fixed point is
certified inside the bridge box.

Consequences, stated exactly:

1. Under P4 plus CP-4, the load-gauge freedom is harmless for the fixed-point
   location and for the bridge value `N_CRC = pi*exp(6*pi/(P*alpha_U))`.
2. The freedom stays load-bearing for the off-fixed-point dynamics: the
   contraction rate at the fixed point is `(1-lambda)*phi'(x_EW)`.
3. A P4-coherence forcing argument alone cannot discharge CP-2: P4 constrains
   the gauge at exactly one point. The discharge route named in the proof
   spine closes the value content of CP-2 and leaves the form content (RC) as
   a separate obligation.

## 5. Negative controls

| family | named verdict | certificate |
|---|---|---|
| inversion `phi = id` | RC residual certified below `1e-40` on the probe grid; coupled fixed point in the bridge box | artifact A1, A2 row 1 |
| power `phi_s(x) = s*x`, four reserve readings | live gauge members; RC violated (residual `1-s`); fixed point off balance by `5.6` to `187.4` load nats; excluded by P4 plus CP-4 | artifact A2/A3 |
| affine capacity `F0 = a*N + pi*(1-a)` | live gauge members; RC violated; fixed point near `x_EW + 2*log(a)`, certified off balance | artifact A2/A3 |
| logarithmic / log-count (CAP-L) | outside LM; certified fixed points at or below `1452.33` nats, `p4_coherent_rows = 0` | `runtime/F_candidate_capL_certificates.json` (recorded run) |
| averaging in capacity (`lambda = 1` constant reading) | CAP-B branch, excluded by spec P5 and barred by V-08 | spec Section 3; coupling theorem S10 |
| shift `phi(x) = x + 1/10` | violates the seed normalization | artifact A2 |
| balance-pinned sine gauge | live gauge member with `phi(x_EW) = x_EW`; RC violated; fixed point certified at the bridge box | artifact A2/A3 |

## 6. Derivation attempt record

Routes attempted for a full discharge, with outcomes:

1. **P4 forcing from the recorded decoupled lattice.** Fails: on the
   2026-07-14 lattice every family registers a P4 discrepancy, including
   CAP-P itself (fixed point `pi` against interior stationary points in
   `[1.93, 18]` nats), so decoupled P4 coherence selects nothing.
2. **Clause-4 pointwise pinning.** Membership clause 4 plus the chart theorem
   pin `G` only on the set of member capacities; with the coupled count
   unconstructed this is at most the two points of CP2-C. No declared theorem
   extends the pinning to the read range.
3. **RC from declared structure.** No declared theorem asserts that the
   sector's reconstruction inverts the certified chart. The former IH-4 model
   `pi/rho_read^2` is superseded as a producer and survives only as an
   independent geometric estimator. RC therefore has no role in constructing
   the direct public-record map.

## 7. Discharge obligations

CP-2 discharges when either (a) RC is derived from declared structure (a
theorem that a self-reading sector's reconstructed capacity carries the load
it reads, for instance from the write/check orientation split of the 24-slot
register), together with the standing structural support for LM, or (b) the
coupled count is constructed and CP-4 becomes a theorem, which by CP2-C
closes the value content and leaves RC as a calibration statement with no
effect on `N_CRC`. Until then CP-2 stands, reduced to LM plus RC with the
surviving freedom named.
