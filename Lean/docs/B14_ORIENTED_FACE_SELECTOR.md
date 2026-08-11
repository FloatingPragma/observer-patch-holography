# B14 oriented-face bracket discriminator

## Exact incidence result

`Screen/OrientedFaceBracketSelector.lean` uses the twenty oriented faces in
the pinned Reynolds port frame. `face_bracket_eq_sixty_r13` kernel-checks the
complete tensor-coordinate identity

\[
B_{\mathrm{face}}=60R_{13}.
\]

The resulting alternating bracket is not a Lie bracket:
`jacobi_failure_witness` gives one exact Jacobiator coefficient equal to
`-1`. The independently replayed Python certificate checks the full 2640
independent-coordinate census: 240 nonzero coefficients, split into 120 at
`+1` and 120 at `-1`.

## Conditional three-norm discriminator

On the same 792 upper-triangular structure-constant coordinates, the exact
certificate compares the face bracket with the complete classified compact
`P`, `F`, and `G` families under three supplied coordinate-edit rules:

| edit | `G` | `F` | `P` |
|---|---:|---:|---:|
| total absolute (`L1`) | `30(sqrt(5)-1)` | `60` | `60` |
| Hilbert--Schmidt squared (`L2^2`) | `(615-123sqrt(5))/22` | `(615+123sqrt(5))/22` | `45` |
| worst coordinate (`Linfinity`) | `(5-sqrt(5))/10` | `sqrt(5)/5` | `1/2` |

The `F` value in the `L1` row is an infimum over its compact stratum, not an
attained minimum. Exact primal-dual Python replay proves the 792-coordinate
optimizations and their feasibility conditions. Lean checks the serialized
radical values and strict family order; `three_norm_unique_nearest_G` proves
that `G` wins every one of the three encoded comparisons.

## Carrier-induced invariant-metric phase diagram

`code/b14_jacobi/invariant_metric_phase.py` removes the basis dependence of
the quadratic rule.  The pinned carrier splits multiplicity-free into sectors
`1 + 3 + 3' + 5`; the commutant of the port action has exact dimension four
and is spanned by the four symmetric spectral projectors, so the complete
family of invariant carrier inner products is the positive sector-scale cone
`alpha P_1 + beta P_3 + gamma P_3' + delta P_5`.  Every such metric induces a
Hilbert--Schmidt metric on the bracket space that is diagonal across the
fourteen channels (all cross pairings vanish identically), and the squared
distances from the face bracket to the classified compact families are the
exact three-term Laurent forms

```
d_P^2 = (15-3s)/beta + (15+3s)/gamma + ((15-3s)/2) beta/delta^2 + ((15+3s)/2) gamma/delta^2
d_F^2 = ((15+3s)/2) gamma/delta^2 + (15+3s)/gamma + (60+12s)/(11 beta)
d_G^2 = ((15-3s)/2) beta/delta^2 + (15-3s)/beta + (60-12s)/(11 gamma)
```

with `s = sqrt(5)`; `alpha` never appears.  The certified phase diagram:

- `P` is strictly excluded for every invariant carrier metric (both gap
  tables are coefficient-positive);
- every sector-balanced metric (`gamma = beta`) selects `G` strictly, with
  gap `3 s beta/delta^2 + (90 s/11)/beta`;
- every metric with `beta/delta` in `[1/50, 6]` selects `G` for all `gamma`
  and `delta`;
- the `F` region is the exact side `d_F^2 < d_G^2` of the displayed tie
  surface; the witness `(beta,gamma,delta) = (8,1,1)` proves that region is
  nonempty, but there is no single threshold in one sector ratio over the
  full three-scale cone;
- conjugating `sqrt(5)` and swapping `beta <-> gamma` maps `d_G` to `d_F`
  exactly: the `F`/`G` asymmetry is the Galois image of the sector swap.
  Global reversal of every face sends the face bracket to its negative and
  leaves all three distances unchanged because the compared compact families
  are centrally symmetric, so face orientation does not select this branch;
- a channel-diagonal invariant metric on the bracket space that is not
  carrier-induced reverses the balanced-point selection, so
  carrier-inducedness is load-bearing (carrier-induced metrics force the
  normalized `t_pp_to_p` and `t_pf_to_f` weights to share one monomial).

`verify_invariant_metric_phase.py` independently replays the commutant
system, the projector membership, every sample point through the raw carrier
metric, and every sign fact.  `Screen/OrientedFaceInvariantMetric.lean`
proves the phase consequences as quantified real theorems over the closed
forms (`dG2_lt_dP2`, `dF2_lt_dP2`, `balanced_unique_nearest_G`,
`box_unique_nearest_G`, `F_wins_at_witness`) and ties the reference point to
the pinned selector distances.

## Open boundary

The face incidence is pinned input.  For the quadratic rule the metric choice
is now characterized rather than declared coordinate data: the comparison is
robust over the complete carrier-induced invariant class, with balance in the
two three-dimensional sectors as the stated sufficient condition and an exact
three-scale tie surface beyond it.  The nearest-point repair rule, the
restriction to carrier-induced metrics, and the `L1`/`Linfinity` coordinate
rules remain declared discriminator choices, and none of them is
source-derived.  Neither a minimum-distance Jacobi repair law nor a physical
bracket producer is selected by the source.  The comparison ranges only over
the compact locus classified under the three named compact-Lie inputs; it is
not a classification of the full Jacobi variety or its noncompact
components.  Ordered source tomography, source bracket selection, and
same-current holonomy remain open under issues #705 and #697.
