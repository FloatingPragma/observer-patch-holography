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
candidate executable model [A-contrast-memory-spherical-v1](../../model/README.md) has no regional entropy, area
or stress readouts, so it does not instantiate this regional limit.
The particle and spin packages likewise supply no missing regional input.

The nonlinear conclusion requires actual finite-ball deficits for the
admitted states. A first variation about one reference state gives only the
linearized equation. The finite rational profiles in the receipt are
algebraic checks, without source or continuum identification.

```sh
python -m pytest -q -W error code/matched_einstein_limit
PYTHONPATH=code python -m matched_einstein_limit --write
PYTHONPATH=code python -m matched_einstein_limit --check
PYTHONPATH=code python -m matched_einstein_limit --lean-check
```

The Python receipt checks rational arithmetic and binds the Lean sources.
The final command requires the repository's Lean toolchain and executes the
kernel audit; the Lean build checks the same module through its umbrella import.

## Reconstructing an Einstein–scalar field on a controlled source

Use three orthogonal radar arms and a local scalar meter on a homogeneous diagonal Bianchi I spacetime. The source evolves the coupled Einstein and Klein–Gordon equations with a quadratic potential. All probes use the same comoving proper clock. Infer the mass squared from the measured scalar equation, the coupling from averaged axial null contractions and the vacuum term from the time-time component at one appointment on one launch. Freeze all three, then predict the equations at three other appointments on two independently initialized launches. Compare radar arm spacings and differencing intervals of 1/20, 1/40 and 1/80. Withhold the evolution state and source coefficients from the decoder. Accept the reference reconstruction if at the finest spacing its largest absolute scalar-equation and Einstein-component errors are both below 0.005 source units and both error envelopes fall at each refinement. Turn off the field's contribution to geometry while retaining the same dynamical scalar and radar readers, and require a held Einstein-component error above 0.05 with the frozen coupling. This is a controlled test of what the observations can identify on homogeneous solutions of a declared action. The field value, proper clock, three orthogonal comoving neighbours and reflection of light are explicit instrument assumptions.

At the finest spacing the observer recovered $m^2=0.699938$, $8\pi G=1.600919$ and $\Lambda=0.119853$, in source units, against the separately recorded inputs 0.7, 1.6 and 0.12. Across six held records the largest Einstein-component error was 0.0002785 and the largest scalar-equation error was 0.0000511. The potential fixes the vacuum energy convention by setting $V(0)=0$. An arbitrary additive potential constant is degenerate with the vacuum term. The control fixes the scalar, its velocity and two initial expansion rates, adjusts the third to the vacuum constraint, and keeps the scalar evolving and gave an Einstein-component error of 0.5884. A secondary forward check initialized the reconstructed equations at the first held appointment on each new launch and forecast the three arm lengths and scalar reading at each of the two later appointments. The largest absolute error in those normalized readings was 0.0000442. The test uses homogeneous diagonal metrics, comoving observers, ideal returned signals and a scalar meter with a declared calibration. It identifies coefficients within that family. Native Interval Echo still requires the corresponding field, stress and geometry instruments on one state. [Source and decoder](observer_field.py) · [Integer clock and scalar records](observer_field.json)

Replay from the repository root with `python code/matched_einstein_limit/observer_field.py`.

## Does measured quantum stress predict independent optical records?

With six quantized scalar modes on the semiclassical box, rotating the occupied modes changes the pressure at fixed initial geometry and energy. One calibration predicts held diagonal curvature to 0.00000559 in box units. A forecast initialized from field quadratures and radar at one appointment predicts later radar ranges to 0.00000999. In the separately quantized Wilson–Maxwell box a prepared photon carries measurable directional stress and leaves the weak optical response invariant. Its geometry is prescribed. [Quantum scalar and radar](quantum_scalar_radar/README.md) · [Photon stress and independent Maxwell probes](quantum_stress_control/README.md)
