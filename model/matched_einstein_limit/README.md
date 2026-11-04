# From regional entropy balance to Einstein dynamics

Ten explicit observers reconstruct every component of a symmetric spacetime
tensor. Their exact rational decoder amplifies the largest measured error by
at most `17/8`, and that constant is optimal for this frame. Applied to one
common family of regional entropy, area and matter-stress readouts, this
proves a conditional Einstein limit with quantified small-ball errors.
The [proof](PROOF.md) states the physical inputs; the
[Lean theorem](../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean)
checks the tensor algebra, error bounds, scaling and conservation consequence.

The required regional readouts must come from the same admitted state and
refinement, in the same proper-time normalization. The area is that of actual
small geodesic balls, the stress is the matter source in the action, and the
entropy is regional generalized entropy on its specified variation space.
If their errors are smaller than the fourth power of the proper radius,
and their common area coefficient tends to a positive spacetime constant,
the reconstructed fields obey `G_ab+Lambda g_ab=8 pi G_N T_ab`, with
`G_N=1/(4 eta)` in units `hbar=c=1`. Independent stress conservation and the
Bianchi identity make Lambda constant on each connected component.

This proves the implication without separately assuming null-tensor balance.
It supplies no executable source satisfying the regional hypotheses. The
[finite reciprocal process](../reciprocal_area_dynamics/README.md) has a
derived global thermal area coefficient and responsive geometry; its two
geometry labels and thermal entropy do not instantiate this regional limit.
The particle and spin packages likewise supply no missing regional input.

The nonlinear conclusion requires actual finite-ball deficits for the
admitted states. A first variation about one reference state gives only the
linearized equation. The finite rational profiles in the receipt are
algebraic checks, without source or continuum identification.

```sh
python -m pytest -q -W error model/matched_einstein_limit
python -m model.matched_einstein_limit --write
python -m model.matched_einstein_limit --check
python -m model.matched_einstein_limit --lean-check
```

The Python receipt checks rational arithmetic and binds the Lean sources.
The final command requires the repository's Lean toolchain and executes the
kernel audit; the Lean build checks the same module through its umbrella import.
