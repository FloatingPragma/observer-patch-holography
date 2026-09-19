# Issue 776: declared paired experiment

This specification is written before the production runs. No dimension is an
acceptance threshold. The contract is the Objective, Deliverables and Exit in
https://github.com/FloatingPragma/observer-patch-holography/issues/776.

## Frozen inputs and execution

* Use the public simulator revision
  `7faa47b5cf00b42f6bf6b3e95ed4f7eb64f4239f`. Pin the four transitive geometry
  source files by bytes. Its unpublished `support_wiring.py` is not an input.
* Capture levels 3, 4, 5 (1280, 5120, 20480 triangular cells). W12 means all
  cells sharing a mesh vertex. W3 means cells sharing an edge. Assign ports
  separately for each graph using the public geometric endpoint assignment.
  Reproduce every existing L3 W12 face and port pair before using larger levels.
* The twelve port labels and thirty seams are those of the committed L3
  fixture. Every sweep attempts every intra-carrier and glued seam once.
  Greedily colour the sorted thirty template seams, using the first available
  colour at each endpoint. Execute those colour classes in order on all cells,
  then all glued seams (a matching). Simultaneous disjoint pairs have no causal
  edges between them. Four sweeps per level; no stopping on a desired result.
* Initial cumulative load at port `i` is
  `SHA256("oph776:20260919:" + decimal(i))[0] mod 8`. Each operation replaces
  `(a,b)` by `((a+b)/2,(a+b)/2)`, including identity operations. Use exact dyadic
  arithmetic, not nearest-integer agreement. The cumulative signed repair load
  includes preparation: initial load plus the sum of signed update increments.
  It therefore telescopes to the current port load. It is not an event counter.
* A mean action has two simultaneous write events, one per destination port.
  Both read the same two previous writer versions. Record initial writes,
  every intra-carrier write, every glued write, and every refinement copy.
  The event IDs serialize custody, not additional precedence. Bind each chunk
  and its metadata in a SHA256 chain rooted in the specification and geometry.
* Execute one tower L3 -> L4 -> L5; its L3 prefix is the single-level run.
  At each join use the committed cell refinement's childwise-constant pullback
  on all twelve coordinates. Each child-port COPY event reads the current
  parent-port writer. This is the declared componentwise extension in common
  port labels, not a derived physical frame transport. Retain the area weights
  and check conditional expectation after embedding. No hidden join edges.
* After each simultaneous phase annotate each event with its owner's current
  rank-three position `x = loads @ F`, where `F F^T = 4 P_slow`, and separately
  its unit-sphere cell centre. Positions are annotations; computing a position
  does not add reads to a seam event. One global orthonormal basis fixes the
  declared identification of local rank-three spaces.

## Order and placement protocol

Strict precedence is the transitive closure of authenticated read-from edges,
including intra-carrier and join events. A diamond is the inclusive order
interval `[a,b]`, not an externally imposed Euclidean light cone. Count `N` and
strict comparable pairs `C`; report `r=2C/(N(N-1))` and invert
`Gamma(d+1) Gamma(d/2)/(2 Gamma(3d/2))`. Undefined cases remain null.

Use two predetermined anchor-selection rules, applied to the snapshot after
the first L3 sweep: (a) the carrier whose readback is nearest the mean carrier
readback, (b) the carrier whose centre is nearest the north pole. Ties use cell
ID. This permits the placements to select different intervals without making
geometry into causality. Retain both placements for every selected interval.
The same fixed interval must give exactly the same order statistics in both
placements; that is a control, not a missing measurement.

For each anchor, the bottom is its port-0 writer after L3 sweep 1. Single-level
tops are its port-0 writers after sweeps 2 and 3. Tower tops are the first
ordered descendant's port-0 writers after sweep 2 at L4 and L5. These tips have
executed past and future. The S2 support has no spatial boundary; report the
interval's cell coverage and spatial extent under each placement. No tangent
patch, spatial cutoff, or excision of pentagonal cells is used. Report all
interval sizes against the actual phase gap and longest-chain height, with
successive count ratios. Do not fit a continuum volume law to four points.

For the record-metric q=13 family, reconstruct the committed golden-coordinate
neighbour predicate and execute its `value(t+1)=1+sum(neighbour values(t))` law
for four layers, including self reads. Retain its complete versioned read
graph and values. The central inclusive intervals for lags 1..4 supply the same
`N,C,r,d` and growth readouts. Report the existing exact geometric interior
test; clipped intervals cannot be used as interior results. State the different
event granularity, population and scheduling conventions beside comparisons.

## Controls and confluence

For all three levels compare W12 with matched W3 and isolated carriers.
Finite canonical runs use the same preparation and four sweeps, forward and
reversed phase orders. Report conservation, quadratic descent, distance to the
component mean, and the disagreement of finite schedule endpoints. Distinguish
finite schedule dependence from the common fixed point: connected components
have a unique prescribed-mean fixed point, and repeated fair cyclic schedules
converge to it. Arbitrary unfair schedules and finite exact termination are not
claimed. Keep historical integer-law results separately labelled.

Compute normalized response kernels using the archived definition
`T=I-L/D`, `D=2|seams|/carriers`, `K_n=12 Q [T^(2n)]_cc Q / trace(...)`,
at n=1,5,30,100,300. These are expectation-operator probes, not the trajectory
of the particular phase schedule. At each level and wiring use four cells:
first and last defect cells, first and last regular cells, sorted and deduped.
Report each sample, not population-wide extrema. Compare with the isolated
analytic kernel and the explicitly L6, 64-cell historical W3 sample.

## Exit and verification

Produce a receipt containing wiring census, exact canonical trajectory custody,
confluence diagnostics, response kernels, both provenance placements, q13
readouts, controls, source pins and claim boundaries. Mirror the receipt into
`evidence/source_net_causal_poset/`. An independent verifier must reconstruct
current writers, check exact means and copies, recompute geometry incidence,
all measured interval counts, readbacks, confluence summaries and kernels.
Mutation tests must reject stale writers, changed means, omitted operations,
bad joins, altered counts and broken custody. Record executed commands and
results. A partial trace, a numeric placeholder or a desired dimension cannot
close the issue. No M1 derivation, source selection, physical clock or continuum
limit is claimed.
