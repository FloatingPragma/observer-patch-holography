# Canonical support wiring beside the record-metric route

This is the paired diagnostic for [#776](https://github.com/FloatingPragma/observer-patch-holography/issues/776).
The [receipt](support_wiring_receipt.json) includes both routes and their
controls. Its [byte-exact mirror](../source_net_causal_poset/support_wiring_receipt.json)
is in the causal-poset package. The [frozen specification](../../code/support_wiring/SPECIFICATION.md)
was committed as `b013a019` before production results. All numbers below are
finite declared-experiment results, not acceptance targets.

## Executed architecture

The public geometry primitives at revision
`7faa47b5cf00b42f6bf6b3e95ed4f7eb64f4239f` reproduce every previously captured
L3 face and W12 port pair. The twelve-port carriers have thirty intra-carrier
seams. Every sweep visits every intra-carrier and glued seam once using exact
pair means. Four sweeps execute at each of L3, L4, L5, with the committed
childwise-constant cell pullback between levels, applied to each port label.
All 307,200 refinement writes explicitly read their parent-port version.

The retained trace has **8,063,280 write events**, comprising 15,360 initial
writes, 307,200 refinement copies and 7,740,720 writes from 3,870,360 mean
actions. Of those actions, 644,760 are glued-seam attempts. Dyadic numerators
are retained exactly in three unsigned 64-bit limbs; the final denominator is
`2^96`. The two destination writes of each mean action are simultaneous and
incomparable. Intra-carrier ancestry is retained, not silently skipped.

The cumulative signed load includes preparation and telescopes to the current
load. Each writer's own twelve loads give its rank-three position; cell centres
are a second annotation. No geometric cone supplies any provenance edge.
The joins use a declared common port labelling. The committed transport is
the commutative cell scaffold, not a derived physical frame transport.

## Wiring census and canonical controls

| Level | Carriers | W12 glued seams | W3 glued seams | W12 unused ports | Antipodal-label matches, W12 / W3 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 3 | 1,280 | 7,650 | 1,920 | 60 | 3,926 / 1,356 |
| 4 | 5,120 | 30,690 | 7,680 | 60 | 17,498 / 5,796 |
| 5 | 20,480 | 122,850 | 30,720 | 60 | 74,836 / 23,804 |

There are sixty cells touching the twelve pentagonal vertices at every level;
each has eleven W12 neighbours and exactly one unused port. Every other cell
uses all twelve ports. W3 uses three ports per cell. The receipt also reports
port-direction dot products after transformation to the declared surface
frames. Literal antipodal labels and geometric direction agreement are
different statistics; neither is imposed as a success condition.

For `n=300`, the declared four-cell response-kernel samples give:

| Level | W12 slow-band share range | Matched W3 range | Isolated |
| --- | ---: | ---: | ---: |
| 3 | 0.900704–0.927947 | 0.640916–0.641966 | approximately 1 |
| 4 | 0.900660–0.928524 | 0.640916–0.641966 | approximately 1 |
| 5 | 0.900648–0.927897 | 0.640916–0.641966 | approximately 1 |

These are sample ranges, not bounds over all carriers. All five probe times
`n=1,5,30,100,300`, all matrices, eigenvalues and per-cell shares are retained.
The historical L6 W3 64-cell sample has n=300 median 0.641256 and range
0.635905–0.642201. The isolated historical n=300 kernel agrees with `4P_slow`
to the archived tolerance. Its complete five-time readout is included too.
Regulator, preparation, schedule and budget differences are explicitly labelled.

Every canonical mean conserves the component total and satisfies the exact
quadratic descent identity. Repeated forward or reversed complete sweeps
converge to the preserved component mean; the finite-dimensional proof and
its limits are in the [code README](../../code/support_wiring/README.md).
The four-sweep endpoints differ by up to 0.386131, 0.362616 and 0.367399 for
W12 at L3, L4 and L5. They have **not** settled to an exactly common endpoint.
The historical floating W3 L6 run was also budgeted, not terminated. Historical
integer nearest-agreement confluence is labelled separately and does not
supply any provenance result here.

## Paired order readouts

Diamonds are inclusive order intervals between the predetermined tips. Both
placements select an anchor before the counts are computed. Every interval
retains both placements and their spatial extents. Replotting the same interval
cannot change its order fraction. All W12 intervals have executed past and
future on a closed S2 support with no spatial cutoff.

| Anchor rule | Top | Phase gap | Events N | Strict pairs C | Ordering fraction | MM dimension |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Writer readback | L3 sweep 2 | 8 | 9 | 32 | 0.888889 | 1.18605 |
| Writer readback | L3 sweep 3 | 16 | 75 | 1,149 | 0.414054 | 2.24830 |
| Writer readback | L4 sweep 2 | 41 | 2,594 | 705,947 | 0.209908 | 3.10501 |
| Writer readback | L5 sweep 2 | 74 | 24,772 | 51,237,452 | 0.166999 | 3.38441 |
| Cell centre | L3 sweep 2 | 8 | 9 | 32 | 0.888889 | 1.18605 |
| Cell centre | L3 sweep 3 | 16 | 75 | 1,149 | 0.414054 | 2.24830 |
| Cell centre | L4 sweep 2 | 41 | 2,802 | 760,376 | 0.193766 | 3.20318 |
| Cell centre | L5 sweep 2 | 74 | 27,369 | 55,882,840 | 0.149213 | 3.52077 |

The receipt includes successive count ratios, longest-chain heights, event-kind
censuses, cell coverage and complete interval membership arrays. These are raw
count-growth readouts, not a fit of a volume exponent or a physical clock.
In particular, the small single-level intervals do not establish a stable
2+1 reading.

The record-metric q=13 route executes **1,176,764 reads** and 10,985 writes,
with its original neighbour predicate, self reads, preparations and four-layer
read law. Its exact counts agree with the previously committed family:

| Lag | Events N | Strict pairs C | Ordering fraction | MM dimension | Spatially interior? |
| --- | ---: | ---: | ---: | ---: | --- |
| 1 | 2 | 1 | 1 | 1 (endpoint convention) | Yes |
| 2 | 181 | 359 | 0.0220381 | 5.75966 | Yes |
| 3 | 360 | 16,966 | 0.262550 | 2.82788 | Yes |
| 4 | 1,529 | 102,990 | 0.0881646 | 4.14933 | **No: clipped** |

The old q13 receipt already flagged lag four as clipped. Calling its 4.15
estimate an interior result would be incorrect. The lag-one endpoint `r=1`
is reported as `d=1`; the historical estimator returned null there. Neither
choice changes any count. The contrast between the routes is measured under
different laws, event granularities and clocks; it does not isolate wiring as
the sole cause of the gap, select the metric read law, or derive M1.

The requested second audit added the q13 family's existing lower-dimensional
scientific controls under an explicit [post-run declaration](../../code/support_wiring/CONTROL_ADDENDUM.md).
They execute the same metric threshold and read law with one or two spatial
axes, retaining all writer versions and integer values (316 and 20,196 reads).
The complete four-lag readouts are in `q13_controls.json` and the paired receipt:

| Spatial axes | Lag | Events N | Strict pairs C | MM dimension | Spatially interior? |
| --- | ---: | ---: | ---: | ---: | --- |
| 1 | 1 | 2 | 1 | 1 | Yes |
| 1 | 2 | 9 | 15 | 2.24012 | Yes |
| 1 | 3 | 16 | 66 | 1.87199 | Yes |
| 1 | 4 | 29 | 202 | 2.00658 | **No: clipped** |
| 2 | 1 | 2 | 1 | 1 | Yes |
| 2 | 2 | 39 | 75 | 3.98565 | Yes |
| 2 | 3 | 76 | 986 | 2.47979 | Yes |
| 2 | 4 | 205 | 4,610 | 3.04457 | **No: clipped** |

Their flat-space reference fractions are 1/2 and 8/35, corresponding to MM
dimensions 2 and 3. These are comparisons, not acceptance thresholds. Exact
counts at all eight intervals agree with the historical source family; the
lag-four agreement with a reference dimension does not remove its clipping.

## Contract and verification

| Issue requirement | Retained evidence and check |
| --- | --- |
| Full L3–L5 canonical architecture | `trace/`, pinned geometry; independent exact replay of every attempted seam and copy |
| Wiring, antipodes, pentagonal deficits, W3 comparison | `wiring` in the paired receipt; full incidence, optimal assignment and connectedness checks |
| Confluence and slow-band share beside historical controls | `controls.json`, `kernels.npz`, `historical_L6`; forward/reverse finite residuals and independent full-field Gram checks |
| One-level and joined-tower provenance, both placements and growth | `provenance.json`, `interval_members.npz`; reverse descendant count independently reproduces all strict-pair counts |
| Same q13 readouts and scientific controls | `q13.json`, `q13_reads.npz`, `q13_controls.json`, both `q13_control_d*_reads.npz`; exact metric predicates, complete readback metric identities, versioned read laws and counts |
| Independent verifier, tests, causal-poset mirror | `code/support_wiring/verify.py`, mutation tests, dedicated CI, byte inventory and mirrored receipt |

Run `python code/support_wiring/verify.py` for the full check. The verifier
imports neither producer. The mutation suite rejects changed laws, rounded or
changed means, stale writers, omitted seams/phases, incorrect joins, altered
readbacks, counts, kernels, confluence residuals and q13 clipping flags, including
corruptions whose custody hashes have been recomputed. Fresh reproduction
commands and the analytic confluence argument are in the code README.

The Exit is the paired receipt and controls above. It is not a target dimension,
an M1 derivation, a source-selection theorem, a physical clock or a continuum
limit. The September 10 partial integer-law log and its draft numerical prose
are not used as evidence.
