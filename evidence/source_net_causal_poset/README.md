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

Two receipt generations are retained. The current generation,
[`source_net_causal_limit_receipt_2026-09-24.json`](source_net_causal_limit_receipt_2026-09-24.json)
and [`carrier_source_net_receipt_2026-09-24.json`](carrier_source_net_receipt_2026-09-24.json),
counts every strict pair exactly at every level through `q = 89` (704,969
sites, ten rounds) in both lanes; the property table below reads from it. The
first generation, `source_net_causal_limit_receipt.json` and
`carrier_source_net_receipt.json` of 2026-09-09 (stratified estimates at
`q = 55`, carriers through `q = 34`), stays under its original names because
the support-wiring, crossing-read, feedback-transport, native-program and
postdiction-ledger packages pin its bytes. Every row shared by the two
generations is identical; the sampled `q = 55` readings of the first,
0.0926 (SE 0.0003) and 0.1057 (SE 0.0015), enclose the exact 0.0925 and
0.1077 of the second. Each generation keeps its own frozen generator,
`build_causal_poset.py` for 2026-09-09 (pinned by the packages above) and
`build_causal_poset_2026-09-24.py` for the current receipts; the verifier
rebuilds both.


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
   `q_i(j + 1) = 1 + sum` of the values read. `K_q = ceil(sqrt q)` rounds. Levels
   `q = 5, 8, 13, 21, 34, 55, 89` (125 to 704,969 carriers).
5. Order. Precedence is generated from the read-after-write provenance of
   those reads with no declared parent lists. It equals the layered order
   `(j, s) <= (j', t)` iff `d(s, t) <= j' - j`, with `d` the graph distance
   on the neighbour graph.

The population, the read law, the round duration `Delta_q = a_q/c` and one
counted event per carrier and round are supplied, as in the source paper.

## Properties

| Property | Value | Status |
| --- | --- | --- |
| Precedence is a partial order | reflexive, antisymmetric, transitive | proved, `Lean/Geometry/SourceNetLayeredOrder.lean`, `layerPrec_isPartialOrder` |
| Precedence equals the reflexive transitive closure of one-round reads | both inclusions | proved, `layerPrec_iff_reflTransGen` |
| Graph-distance characterization | `(j,s) <= (j',t)` iff `d(s,t) <= j' - j` | proved, `walk_iff_walkDist_le`; concrete bridge `precedes_iff_layerPrec` |
| Width | equals the site count `q^3` | proved, `width_eq_card_sites`; computed at every `q` |
| Height | every event at round `j` has height `j` | proved, `height_eq_layer` |
| Cone sandwich | `displacement <= k(a - 2h)` implies reachability; reachability implies `displacement <= k a`, for `k >= 1`, an `h`-covering population in a convex window, and `a > 2h` | proved, `Lean/Geometry/SourceNetCausalCone.lean`; computed cone checks at every `q` with zero violations |
| Central-diamond dimension, 3D family | 2.63, 2.93, 4.15, 3.61, 4.10, 4.09, 4.07 at `q = 5, 8, 13, 21, 34, 55, 89` | computed, `vertical_intervals[-1].myrheim_meyer_dimension`; every one of these outermost diamonds has `continuum_diamond_inside_cube = false` (`central_diamond_inside_cube_3d` in the manifest), so none is an interior-diamond reading; the largest interior diamonds, `vertical_intervals[-2]`, read 4.04, 5.02, 2.83, 4.38, 3.55, 3.84, 3.92 |
| Ordering fraction, 3D family | 0.306, 0.241, 0.088, 0.138, 0.092, 0.093, 0.094 against the flat value 1/10; exact all-pairs counts at every level | computed; the constant 1/10 is proved, `Lean/Geometry/OrderingFractionFourDimensional.lean`, `orderingFraction_eq` |
| Layer parity | even `K` approaches 1/10 from below, odd `K` from above | computed |
| Moving-tip diamonds (boost) | 0.250, 0.135, 0.121, 0.115, 0.108, 0.102 at `q = 8, 13, 21, 34, 55, 89`; dimension 2.89, 3.64, 3.77, 3.83, 3.91, 3.97 | computed, `moving_tip_interval`, exact all-pairs |
| Control families | 2D family reads 3.04, 2.78, 3.01, 3.02, 3.02 and 1D family 2.01, 1.94, 2.00, 2.01, 2.01 at `q = 13, 21, 34, 55, 89` against the proved constants 8/35 and 1/2; their count clocks use the cube and square roots (`clock_exponent`) | computed; constants proved, `myrheimMeyer_three`, `myrheimMeyer_two` |
| Count law | normalized count against `pi c^3 T^4/24`: 0.193 against 0.198 at `q = 13`, 0.141 against 0.147 at `q = 34`, 0.172 against 0.177 at `q = 55`, 0.163 against 0.165 at `q = 89` | computed; the paper's error bounds are recorded and exceed the volume at these `q` |
| Count clock | `(|I|/|J|)^{1/4}` against the model-time ratio: 1.705 against 2 at `q = 13`, 2.088 against 2 at `q = 34`, 1.907 against 2 at `q = 55`, 2.051 against 2 at `q = 89` | computed; enclosure proved, `Lean/Time/SourceCountClockEnclosure.lean`, `finite_clock_enclosure` |
| Readback metric | equals the source metric with scale exactly one | computed exactly in `Q(sqrt 5)`, `readback_metric` |
| Neighbour relation from readbacks | equals the family's relation edge for edge at `q = 5, 8, 13, 21, 34, 55, 89` | computed, `neighbours.equals_source_net_digest` |
| Provenance order | equals the layered order on the centre diamond at every `q`; derived rank equals the round | computed, `provenance` |
| Intervention | a `+1` at the centre changes exactly the future cone in every round | computed, `intervention` |
| Theory replay | 28 structural fields byte-equal to the source paper's replay at `q = 5, 8, 13` | computed, `rer_cross_check` |
| Operation costs | reads per round `2E + n`; at `q = 34`: 157,922,880 reads, 275,128 writes, 2.40 GB under the declared byte model; at `q = 89`: 21,845,550,850 reads, 7,754,659 writes, 424 GB | computed, `operation_costs` |
| Continuum reconstruction | On past- and future-distinguishing Lorentzian spacetimes, chronological order determines the conformal geometry; an identified metric volume measure fixes the conformal factor | imported continuum theorem; not a finite-poset existence theorem |
| Selection of the population and read law by native repairs | | not supplied by these receipts |
| Physical clock and source-selected continuum | not supplied by these receipts; the declared-family analytic limit is in the source paper | distinct from finite numerical diagnostics |
| Isotropy of link directions, interval abundance profile | odd multipoles vanish, `l = 4` power 0.036, 0.023, 0.0060, 0.0016 at `q = 8, 13, 21, 34`; intervals of 64 to 1024 events match a flat sprinkling within 14 percent, links and small intervals do not (layered order) | computed at `q <= 34`, [extension of 2026-09-25](#extensions-of-2026-09-25) |
| Chain counts `C_2`, `C_3`, `C_4` and interval spectrum | central dimensions 4.076, 4.016, 3.971 at `q = 144`; central 25-point CDF distances 0.075810, 0.024725, 0.032342 at `q = 55, 89, 144` | sampled point estimates, no uniform monotone improvement or joint statistical decision; complete [refinement interface](manifold_refinement/README.md) retains all regions and controls |
| Curvature estimators | | not supplied by these receipts |

## Independent checks

`build_causal_poset.py` rebuilds the poset at a given `q` from the definitions
above with the standard library and numpy only (no simulator code), and
compares 99 quantities per level with the mirrored receipts: records, readbacks,
the metric identity, the neighbour graph and its digest, the strict pair count
of the centre diamond, the ordering fraction, the dimension, the count clock,
the read chains, the intervention, the operation costs, and, when the stored
log is present, the provenance order generated from the log alone.

```bash
python3 evidence/source_net_causal_poset/build_causal_poset_2026-09-24.py --q 5 8
python3 evidence/source_net_causal_poset/build_causal_poset.py --q 5 8   # the 2026-09-09 generation
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

## Extensions of 2026-09-24

Two receipt-level extensions sit beside the current generation and are checked by the
package verifier without simulator code:

- `flrw_record_density_readout_2026-09-24.json` (generator
  `flrw_record_density_readout_2026-09-24.py`): the FLRW count measure
  `M_sigma = sum sigma_j^4 n_j` read off the per-layer counts of the current generation's
  vertical diamonds for four supplied scale-factor profiles (constant, de Sitter, radiation,
  matter, normalised to `sigma_0 = 1` and `sigma_K = 2`). It records the sandwich, the
  expanding count clock with its enclosure and the proper-time ratio of a comoving worldline,
  and the redshift reading between congruent reference diamonds at the emission and reception
  epochs. Every identity of `Lean/Geometry/SourceNetConformalRecordDensity.lean` holds at every
  level through `q = 89`. The verifier recomputes every mass, clock, enclosure and redshift
  reading. The profile is a supplied datum; the orders are the flat comoving orders; no
  physical clock, scale or curvature is attached.
- `source_net_causal_limit_q144_sampled_2026-09-24.json` (driver
  `q144_driver_2026-09-24.py`): the record-metric lane at `q = 144` (2,985,984 sites, twelve
  rounds) with stratified sampled pair counting on the three-dimensional family and exact
  counting on the controls. Readings: ordering fraction `0.0955 +- 0.0002` at `K = 12`, with
  the even-`K` ladder `0.0600, 0.0803, 0.0887, 0.0932, 0.0955` at `K = 4, 6, 8, 10, 12`;
  Myrheim-Meyer dimension `4.055`; largest interior diamond `3.946`; moving tip
  `0.1024 +- 0.0004`; count clock `1.959` against the model-time ratio `2`; controls `3.017`
  and `2.004`. These are sampled estimates with standard errors; the exact family receipt and
  its claim boundary (`q <= 89`) are unchanged, and the level has no theory replay row.

## Extensions of 2026-09-25

Two lanes read the order with more than one number; both are receipts of the simulator with
independent verifiers and tests there, and the package verifier checks each of them here
against the exact family receipt (event counts, pair counts, tip sites, edge counts), the
generator pins and every derived field.

- `manifold_observations_receipt_2026-09-25.json` (generator
  `manifold_observations_2026-09-25.py`, the simulator module byte for byte; q = 8, 13, 21,
  34; 722 s): Myrheim-Meyer dimension on interior diamonds of three sizes at the centre and at
  up to sixteen off-centre tips, moving-tip diamonds at four rapidities, interval-abundance
  profiles against fixed-count sprinklings of flat diamonds, the Benincasa-Dowker action,
  link-direction multipoles against a uniform reference (distinct directions counted exactly
  in `Q(phi)`), count clocks and a refinement trend table; every abundance profile is the
  exact between-count matrix. At q = 34: the dimension is the same at every off-centre tip of
  the same size to within 0.04; the four moving-tip diamonds read 0.108 to 0.125 (dimensions
  3.73 to 3.92) against 0.145 for the five-layer vertical diamond; intervals of 64 to 1024
  events are as abundant as in a flat sprinkling within 14 percent while links are 23 percent
  of related pairs against 8.6 percent; the odd
  multipoles of the link directions vanish and the `l = 4` power falls 0.036, 0.023, 0.0060,
  0.0016 over the four levels; the controls read 3.01 and 2.00. The family is a layered order,
  not a sprinkling: links, small intervals and the action read the lattice scale and are
  controls; the bulk statistics are the ones the limit governs.
- `source_net_manifold_sampled_q55_q89_2026-09-25.json` (generator
  `source_net_manifold_sampled_2026-09-25.py`, the simulator module byte for byte; 3,451 s):
  chain counts `C_2`, `C_3`, `C_4` with standard errors from uniformly drawn events and
  partners (breadth-first searches with the exact cone pruning of the causal-limit lane), the
  dimensions they invert to under `C_k = N^(k) chi_k(d)`, and the interval spectrum, the
  distribution of `|I(y,z)|/N` over related pairs against the continuum `(tau/T)^d` laws in
  d = 2 to 6, on the centre diamond, the moving-tip diamond and fourteen off-centre diamonds
  displaced by `L/8`, with the two- and one-dimensional controls. At q = 89 the sampled `C_2`
  is `5.40e10 +- 0.15e10` against the exact `5.5109e10`; the dimensions from `C_2`, `C_3`,
  `C_4` are 4.094, 4.034, 4.017 (moving tip 3.88, 3.96, 3.93; controls 3.007, 2.997, 3.014
  and 1.99, 2.02, 2.04); the spectrum is closest to the flat `3+1` law with maximum CDF
  difference 0.025 (0.076 at q = 55) and mean 0.0103 against `1/105`; the `C_3` dimension
  over the fifteen six-layer diamonds is 4.078 with standard deviation 0.016 and the
  count-volume coefficient 1.0627 with standard deviation 0.0013.
- `source_net_manifold_sampled_q144_2026-09-25.json` (same generator; 13,188 s): the same
  readouts at q = 144 at half the sample sizes. The three-dimensional site graph was built
  chunk-parallel with the module's per-chunk arithmetic; its neighbour digest equals the
  sampled q = 144 family receipt, and the chunked build equals the sequential build at
  q = 21 and 34. The sampled `C_2` of the centre diamond (4,605,145 events) is
  `9.95e11 +- 0.38e11` against the sampled strict pair count `1.0128e12 +- 0.0025e12`; the
  dimensions from `C_2`, `C_3`, `C_4` are 4.076, 4.016, 3.971 (moving tip 3.99, 4.00, 4.03;
  controls 3.04, 3.02, 3.00 and 2.00, 2.01, 2.03); the spectrum is closest to the flat `3+1`
  law with maximum CDF difference 0.032; the `C_3` dimension over the fifteen eight-layer
  diamonds is 4.061 with standard deviation 0.045 and the count-volume coefficient 1.0224
  with standard deviation 0.0027.

The [refinement interface](manifold_refinement/README.md) combines these retained
sampled receipts in one generated per-statistic table. It separates analytic
chain/mean references from sampled CDF/moment references, exposes reported
marginal errors and unavailable covariance, and retains all 96 regions and all
120 outcomes of the separate small clock control. Its verifier is included in
the archive check. The q144 graph driver and original large-run draws remain
unavailable; archive arithmetic verification does not reconstruct them. The
descriptive failure of uniform monotone improvement is not a significance
test or a rejection of continuum convergence.

## Boundary

The population, the read law, the round duration and one counted event per
carrier and round are supplied. Their selection by accepted native repairs,
the identification of the round with a physical clock, and the continuum
limit are outside this package. Finite runs at `q <= 89` do not demonstrate
the asymptotic statements of the source paper; they read the finite values
above.
