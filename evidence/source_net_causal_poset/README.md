# The OPH causal poset

This package carries the expected causal poset of Observer Patch Holography at
finite regulator, produced by the carrier stack, together with every property
of it that is proved, computed, or declared. The receipts are byte-exact
mirrors of the simulator receipts `data/exact/source_net_causal_limit_receipt.json`
and `data/exact/carrier_source_net_receipt.json` (producers
`oph_exact/source_net.py` and `oph_exact/carrier_source_net.py` in
`oph-physics-sim`, each with an independent verifier and a test file). The
carrier event logs at `q = 5` and `q = 8` are mirrored under
`carrier_source_net_logs/`.

## Definition

Fix a Fibonacci number `q = F_n` with successor `p = F_{n+1}`, the golden
ratio `phi`, and `L = 2/sqrt(phi + 2)`.

1. Records. For every triple `b` in `[0, q)^3` the integer record
   `z(b) = (b2 - m1, b2 + m1, b3 - m2, b3 + m2, b1 - m3, b1 + m3)` with
   `m_i = -floor(b_i phi)` lies in the conservative seam-current image `D_6`
   of one twelve-port carrier.
2. Carriers. One exact twelve-port icosahedral carrier per site holds `z(b)`
   as its port loads on the six antipodal axes.
3. Position. Each carrier reads its position from its own loads through the
   rank-three response, `x = 2 P_slow N`. The readback metric equals the
   paper's source-axis contraction `s(b) = L (xi_{b1}, xi_{b2}, xi_{b3})` with
   scale exactly one, an identity in `Q(sqrt 5)` checked on every pair.
4. Reads. In round `j + 1` every carrier reads the record of every carrier
   whose readback lies within `a_q = L/sqrt(q)`, including itself, and writes
   `q_i(j + 1) = 1 + sum` of the values read. `K_q = ceil(sqrt q)` rounds.
5. Order. Precedence is generated from the read-after-write provenance of
   those reads with no declared parent lists. It equals the layered order
   `(j, s) <= (j', t)` iff `d(s, t) <= j' - j`, with `d` the graph distance
   on the neighbour graph.

The population, the read law, the round duration `Delta_q = a_q/c` and one
counted event per carrier and round are supplied, as in the source paper.

A separate [routed-read extension](routed_read_law/README.md) now executes the
complete q=13/q=21 logical menus on declared W12 L4/L5 support. It retains all
intermediate events and costs and verifies a semantic refinement under a
supplied local M1 feedback law. The original receipts below remain the
declared-family records; the new routing does not derive M1 or turn its
auxiliary events into spacetime-volume counts.

## Properties

| Property | Value | Status |
| --- | --- | --- |
| Precedence is a partial order | reflexive, antisymmetric, transitive | proved, `Lean/Geometry/SourceNetLayeredOrder.lean`, `layerPrec_isPartialOrder` |
| Precedence equals the reflexive transitive closure of one-round reads | both inclusions | proved, `layerPrec_iff_reflTransGen` |
| Graph-distance characterization | `(j,s) <= (j',t)` iff `d(s,t) <= j' - j` | proved, `walk_iff_walkDist_le`; concrete bridge `precedes_iff_layerPrec` |
| Width | equals the site count `q^3` | proved, `width_eq_card_sites`; computed at every `q` |
| Height | every event at round `j` has height `j` | proved, `height_eq_layer` |
| Cone sandwich | `displacement <= k(a - 2h)` implies reachability; reachability implies `displacement <= k a`, for `k >= 1`, an `h`-covering population in a convex window, and `a > 2h` | proved, `Lean/Geometry/SourceNetCausalCone.lean`; computed cone checks at every `q` with zero violations |
| Interior-diamond dimension, 3D family | 2.63, 2.93, 4.15, 3.61, 4.10, 4.09 at `q = 5, 8, 13, 21, 34, 55` | computed, `vertical_intervals[-1].myrheim_meyer_dimension` |
| Ordering fraction, 3D family | 0.306, 0.241, 0.088, 0.138, 0.092, 0.093 against the flat value 1/10 | computed; the constant 1/10 is proved, `Lean/Geometry/OrderingFractionFourDimensional.lean`, `orderingFraction_eq` |
| Layer parity | even `K` approaches 1/10 from below, odd `K` from above | computed |
| Moving-tip diamonds (boost) | 0.250, 0.135, 0.121, 0.115, 0.106 at `q = 8, 13, 21, 34, 55`; dimension 2.89, 3.64, 3.77, 3.83, 3.93 | computed, `moving_tip_interval` |
| Control families | 2D family reads 3.04, 2.78, 3.01, 3.02 and 1D family 2.01, 1.94, 2.00, 2.01 at `q = 13, 21, 34, 55` against the proved constants 8/35 and 1/2 | computed; constants proved, `myrheimMeyer_three`, `myrheimMeyer_two` |
| Count law | normalized count against `pi c^3 T^4/24`: 0.193 against 0.198 at `q = 13`, 0.141 against 0.147 at `q = 34`, 0.172 against 0.177 at `q = 55` | computed; the paper's error bounds are recorded and exceed the volume at these `q` |
| Count clock | `(|I|/|J|)^{1/4}` against the model-time ratio: 1.705 against 2 at `q = 13`, 2.088 against 2 at `q = 34`, 1.907 against 2 at `q = 55` | computed; enclosure proved, `Lean/Time/SourceCountClockEnclosure.lean`, `finite_clock_enclosure` |
| Readback metric | equals the source metric with scale exactly one | computed exactly in `Q(sqrt 5)`, `readback_metric` |
| Neighbour relation from readbacks | equals the family's relation edge for edge at `q = 5, 8, 13, 21, 34` | computed, `neighbours.equals_source_net_digest` |
| Provenance order | equals the layered order on the centre diamond at every `q`; derived rank equals the round | computed, `provenance` |
| Intervention | a `+1` at the centre changes exactly the future cone in every round | computed, `intervention` |
| Theory replay | 28 structural fields byte-equal to the source paper's replay at `q = 5, 8, 13` | computed, `rer_cross_check` |
| Operation costs | reads per round `2E + n`; at `q = 34`: 157,922,880 reads, 275,128 writes, 2.40 GB under the declared byte model | computed, `operation_costs` |
| Continuum reconstruction | On past- and future-distinguishing Lorentzian spacetimes, chronological order determines the conformal geometry; an identified metric volume measure fixes the conformal factor | imported continuum theorem; not a finite-poset existence theorem |
| Selection of the population and read law by native repairs | | is work in progress |
| Physical clock and source-selected continuum | not supplied by these receipts; the declared-family analytic limit is in the source paper | distinct from finite numerical diagnostics |
| Isotropy of link directions, interval abundance profile, curvature estimators, S2 support wiring readout | | is work in progress |

## Independent checks

`build_causal_poset.py` rebuilds the poset at a given `q` from the definitions
above with the standard library and numpy only (no simulator code), and
compares 99 quantities per level with the mirrored receipts: records, readbacks,
the metric identity, the neighbour graph and its digest, the strict pair count
of the centre diamond, the ordering fraction, the dimension, the count clock,
the read chains, the intervention, the operation costs, and, when the stored
log is present, the provenance order generated from the log alone.

```bash
python3 evidence/source_net_causal_poset/build_causal_poset.py --q 5 8
python3 evidence/source_net_causal_poset/verify_causal_poset_archive.py
```

Recorded output of the generator at `q = 5, 8`:

```text
q=5 dim=3: 99 comparisons, all agree
q=8 dim=3: 99 comparisons, all agree
CAUSAL_POSET_REBUILT_AND_EQUAL_TO_RECEIPTS
```

`verify_causal_poset_archive.py` checks every manifest digest, strict JSON
decoding, the receipt schemas and pins, the nonclaim flags, and runs the
generator at `q = 5` and `q = 8`.

## Boundary

The population, the read law, the round duration and one counted event per
carrier and round are supplied. Their selection by accepted native repairs,
the identification of the round with a physical clock, and the continuum
limit are outside this package. Finite runs at `q <= 55` do not demonstrate
the asymptotic statements of the source paper; they read the finite values
above.
