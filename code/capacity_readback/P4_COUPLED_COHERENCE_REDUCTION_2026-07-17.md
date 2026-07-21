# P4 count-density coherence for the coupled readback: reduction

> **Archived target-coupled reduction.** Section and property references below
> refer to the superseded candidate specification. The canonical producer is
> the direct public stable-record capacity in
> [F_READBACK_SPEC.md](F_READBACK_SPEC.md).

Companion to [G2_GAP_1_COUPLING_THEOREM.md](G2_GAP_1_COUPLING_THEOREM.md).
Historical target: GAP-A8, the P4
obligation for the coupled readback (the MAR argmax representation and the
fixed-point representation select the same `N`). Machine checks:
[cp2_p4_premise_reduction.py](cp2_p4_premise_reduction.py), artifact
`runtime/cp2_p4_premise_reduction_check.json`.

**Status: reduction record. Under the per-configuration reading of membership
clause 4, P4 coherence is equivalent to one named count-side premise (CP-4);
under the aggregate reading, coherence holds degenerately by support
collapse. GAP-A8 stays open; no ledger row moves.**

## 1. The four objects in one frame

All objects live over the admissible capacity interval `I`, related by the
proven load chart `X(N) = log(N/pi)`.

- **Finite count object.** `A_r(N) = |Omega^sc_{r,N}|`, the terminal normal
  forms at cutoff `r` and capacity `N` satisfying membership clauses 1 to 4
  of `def:self-closure-density` (repair-closed, observer/checkpoint
  supporting, recovered-package carrying, own-surface readback of `N`).
- **Limiting density.** `l(N) = lim_r log A_r(N) - N`, on the same cofinal
  refinement as `F_r -> F` (spec P1), where the limit exists.
- **Selector and its carrier.** The MAR argmax
  `N_star = argmax_{N in I} l(N)`, non-uniqueness resolved by the MAR
  (Minimal Admissible Realization) lexicographic order; the count
  representation consumes a smoothness/interior-maximum assumption the MAR
  order does not license (`def:self-closure-density`), and the accepted proof
  carrier is the derivative-sign certificate `H_N = l'` with a unique
  positive-to-negative sign change (spec Section 1.3).
- **Fixed-point readback.** `N_CRC = F(N_CRC)` with the coupled map of
  CP-2/CP-3, `F(N) = pi * exp((1-lambda)*X(N) + lambda*x_EW)`, unique fixed
  point `pi * exp(x_EW)` by the Banach certificate (coupling theorem S11).

## 2. Reading dichotomy of membership clause 4

Clause 4 admits two readings, and the P4 obligation is a different statement
under each.

**Aggregate reading.** The readback is a function of `N` alone: every
configuration at capacity `N` reads back `F(N)`. Then clause 4 is
all-or-nothing: `Omega^sc_N` is empty unless `F(N) = N`.

**Theorem P4-A (support collapse).** Under the aggregate reading,
`supp(A) = Fix(F)`, so for any readback map with a fixed point the MAR argmax
over the support equals the fixed point, and P4 coherence holds degenerately.
For a map with no fixed point (the recorded CAP-K family) the selector is
empty. Machine witness (artifact B1): the coupled readback residual is
certified nonzero at every off-balance probe (load offsets `0.1` to `30`) and
encloses zero at balance, so the clause-4-exact support is the fixed-point
box. Consequence: under this reading the count representation carries no
information beyond the fixed point; the two representations are one display,
and the one-capacity rule is satisfied vacuously.

**Per-configuration reading.** Readback varies across configurations at fixed
`N`; `A_r(N)` counts the configurations whose own surface reads `N`; `l` is
nondegenerate. This is the reading under which the count representation is an
independent check, and the remainder of this document addresses it.

## 3. Transfer theorem

**Theorem P4-B (coherence transfer).** Assume CP-1 to CP-3 (so `F` has the
unique fixed point `N_CRC = pi*exp(x_EW)` on `I`) and assume `l` is
differentiable on `I` and carries the accepted derivative-sign certificate
(unique positive-to-negative sign change of `H_N = l'` on `I`). Then

```
MAR argmax of l  =  N_CRC      iff      H_N(N_CRC) = 0,
```

that is, P4 coherence for the coupled readback is equivalent to

- **CP-4 (stationarity at balance).** `l'(N) = 0` exactly at load balance
  `X(N) = x_EW`.

Proof. The sign-change carrier gives a unique interior maximum at the zero of
`H_N`; the maximum sits at `N_CRC` iff the zero does. The converse is direct.
With the stronger concavity bound `-M <= l'' <= -m < 0` (spec Section 1.3)
the same equivalence holds, and on a discrete capacity grid the grid argmax
lies in the two grid neighbors of the continuum stationary point.

Equivalent display. For the coupled map the readback residual carries the
balance sign, `sign(F(N) - N) = sign(x_EW - X(N))` on `I` (load coordinate:
`C(y) - y = lambda*(x_EW - y)`), so CP-4 states exactly that the count
gradient matches the readback residual in sign,
`sign(l'(N)) = sign(F(N) - N)`. A theorem of the form
`F(N) - N = c(N) * l'(N)` with controlled `c(N) > 0` implies CP-4 and is the
constructive route to it.

**Theorem P4-C (status of CP-4).** CP-4 is not derivable from declared
structure: the count object for the coupled membership clause is
unconstructed (coupling theorem Section 6, P4 row), and every constructed
count of the 2026-07-14 run violates it, with interior stationary points
`c/(1-rho) <= 18` nats sitting a certified `>= 121` decimal orders below the
bridge capacity (artifact B5, matching `p4_coherent_rows = 0` of the recorded
run). P4 for the coupled readback is therefore exactly as open as CP-4: the
transfer theorem gives no more, the controls below show no less.

## 4. Controls

Machine-certified in the artifact (`mpmath.iv`, 60 decimal digits, outward
rounding; deterministic, SHA-256 checked across reruns):

| control | content | result |
|---|---|---|
| B1 | aggregate reading: readback residual nonzero off balance, zero at balance | support collapse certified |
| B2 | positive: CP-4 carrier `l = -(X - x_EW)^2/2`; derivative-sign certificate; argmax box against fixed-point box | agreement certified below `1e-8` load nats |
| B3 | location clause dropped: stationarity at `x_EW - 5` and at `0.85*x_EW` | divergence certified (`>= 4.99` and `>= 42.1` load nats) |
| B4 | sign-carrier clause dropped: remote maximum 20 nats below balance, height 300 | argmax leaves the fixed point (value gap `>= 100`) while balance stays stationary to below `1e-80` |
| B5 | recorded-lattice counts, interior stationary points `c/(1-rho)` | separation `>= 121.29` decimal orders from the bridge capacity |

B3 and B4 are the required negative controls: each clause of CP-4 is
load-bearing, and removing either one separates the two representations.

## 5. Discharge obligations

GAP-A8 discharges when either (a) the aggregate reading is adopted as the
declared semantics of clause 4, which closes P4 degenerately and simultaneously
demotes the count representation from independent check to display (a
declaration, not a theorem, and it must be recorded as such in the spec), or
(b) the coupled count object is constructed on the per-configuration reading
and CP-4 is proven for it, which by P4-B closes P4 and by CP2-C
(CP2_INVERSION_FORM_REDUCTION_2026-07-17.md) also forces the fixed-point
location across the load-gauge freedom of CP-2. Until one of these lands, P4
stands as an open obligation reduced to CP-4.
