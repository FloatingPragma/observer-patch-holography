# From regional entropy balance to Einstein dynamics

Ten explicit observers reconstruct every component of a symmetric spacetime
tensor. Their exact rational decoder amplifies the largest measured error by
at most `17/8`, and that constant is optimal for this frame. Applied to one
common family of regional entropy, area and matter-stress readouts, this
proves a conditional Einstein limit with quantified small-ball errors.
The [derivation from patch records to gravity](PROOF.md#from-patch-records-to-gravity) explains the physical inputs, the exact homogeneous null-balance result, and the domain of the galactic profile. The [matched small-ball proof](PROOF.md#the-common-source-interface) states the limiting estimates; the
[Lean theorem](../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean)
checks the tensor algebra, error bounds, scaling and conservation consequence.

The [galactic-outskirts example](PROOF.md#why-the-outskirts-of-a-galaxy-could-show-an-excess) explains what a local record retains and how a compact exterior makes an additional source visible in orbits.

The required regional readouts must come from the same admitted state and
refinement, in the same proper-time normalization. The area is that of actual
small geodesic balls, the stress is the matter source in the action, and the
entropy is regional generalized entropy on its specified variation space.
If their errors are smaller than the fourth power of the proper radius,
and their common area coefficient tends to a positive spacetime constant,
the reconstructed fields obey `G_ab+Lambda g_ab=8 pi G_N T_ab`, with
`G_N=1/(4 eta)` in units `hbar=c=1`. Independent stress conservation and the
Bianchi identity make Lambda constant on each connected component.

On 80 Interval Echo cells, an independently calibrated clock survey has a realizable full-head entropy maximum at fixed survey volume. Adding ordinary exchange energy or weakening correlations between the two heads increases entropy and reduces the surveyed area at the same coefficient. The correlation change preserves each head's records and ordinary energy. The state is deliberately prepared. A fifth station and a second clock menu expose the limits of its spatial interpretation. A source-selected family of proper balls carrying the conserved physical stress would connect this result to gravity. [Finite proof](source_implication/full_algebra/README.md) · [Signed measurements](source_implication/signed_transfers.json) · [Fifth-station test](source_implication/full_geometry/README.md)

On the same connected board, 31 return-time probabilities determine the complete enrolled two-head state and independently predict both material-force readings. The twelve mean ranges determine the surveyed tetrahedron. A physical change at the volume-matched equilibrium of the reused radar preserves all twelve means, both marked forces, both head marginals and the ordinary energy while increasing entropy. The timestamps therefore let observers measure quantities that the surveyed tetrahedron misses. [Clock readout and finite counts](source_implication/full_algebra/reuse/CLOCK_HISTORY.md) · [Fixed-geometry counterexample](source_implication/directional_stress/connection/EQUAL_FORCE_ENTROPY.md)

Independent vacancy and excess measurements recover the same surveyed area. The modular charge also extends to a conserved native observable during flight. A local current derivative agrees exactly with the signed area response. A separate reciprocal early-arrival metric on the very same state gives the opposite area sign for unpairing at fixed volume. These tests distinguish the scalar equilibrium identity from the unproved physical identification of proper null curvature and stress. [Same-source proof and controls](PROOF.md#direct-null-balance-tests-on-the-same-equilibrium)

For nested proper balls, an exact finite shell identity shows when positive entropy increases force smaller volume at the same proper radius. The active four-station survey admits independent changes of volume at fixed entropy and energy. On its calibrated stationary preparation, a five-by-five optical screen detects the two-head correlation, and a sixth station excludes the homogeneous curved-chart interpretation of the survey ranges. [Shell theorem and count comparison](fixed_volume/README.md#when-a-smaller-area-implies-a-smaller-volume) · [Same-source experimental record](PROOF.md#direct-null-balance-tests-on-the-same-equilibrium)

The matched small-ball implication obtains the tensor balance without a separate null-balance assumption. The conditional homogeneous Hamiltonian and the full-head experiment use different geometric inputs.

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

A separate read-only candidate instantiates the three finite small-ball estimates on homogeneous straight sheets. Existing native clock publications choose the source. Finite rods, regional memory entropies and the conserved sheet action supply independent inputs at the same volume. The full native modular-flow and gravitating-matter identifications remain open. Replay the specific family with `PYTHONPATH=code:. python -m matched_einstein_limit --published-records`. [Source, bounds and countertests](../../model/candidate_laws/interval_echo_partition_repair/published_strings/README.md)

## Which inputs does Interval Echo supply?

The quantifier in the Einstein result concerns one family of physical regions, states, clocks and matter readouts. The table gives the strongest conclusions of the listed experiments for each input.

| Input | Finite result and physical scope |
| --- | --- |
| G1. A common spacetime metric | Both mobile charges attain the same bridge velocity. Its dual distance is exactly twice the cell-hop count. Changing either stored head preparation leaves these distances fixed. Proper time, responsive spatial volume and a common causal metric across refinements require independent readings. [Velocity and its source bound](source_implication/common_front/velocity/README.md) |
| G2. Geometric modular flow and the matter stress it measures | On a stationary cell, spoke and rim energy readings reproduce eight angular components in a protected sector. The coefficients depend on the radial position of the seams. Their identification with the geometrically normalized modular generator on the shrinking physical regions has no established transfer. [Angular energy identity](source_implication/modular_orbit/angular_energy/README.md) |
| G3. Independent conservation of that stress | On a twelve-cell cut containing both source heads, native energy and a directed shift charge have exact local balance equations. The shift has a marked-contact force with a nonzero stationary trace moment. On all 119 unmarked bridges, the instantaneous energy current equals 40 times load current on the cell-uniform orbital sector. This sector leaves both source variations invisible at fixed phase and leaks under free flight. A symmetric conserved tensor on a common proper frame requires identifying additional contributions or controlling an effective limit. [Flux and force proof](source_implication/directional_stress/translation/README.md) |
| G4. Stationarity at fixed proper volume | A four-station survey has a unique global entropy-plus-area maximum and a contracting ensemble procedure to prepare it. Repeated use on the same heads selects a different state with a certified entropy-increasing change invisible to all twelve directed mean ranges. [Connected radar test](source_implication/full_algebra/reuse/invisible.json) The product maximum persists for arbitrarily many surveyed copies. A different faithful preparation at the same volume has a signed slope above 7.47 nats per unit odd fraction. [Global theorem and counterpreparation](source_implication/full_algebra/PROOF.md) On a different, protected clock and two-word source, all-radius stationarity along the correlation variation forces the normalized linear collar signal to vanish. Even exact first-order stationarity on all balls in one static rest foliation leaves directional curvature undetermined. [Routed-ball reference theorem](source_implication/modular_orbit/relative_entropy/README.md) · [Rest-ball counterexample](rate_geometry/README.md#static-null-curvature-from-a-directional-rim-and-clock-array) |
| G5. Direction coverage and small-ball errors | A positive pole weighting on the twelve-port angular instrument restores all six even components. On a known static geometry, measured directional rim areas with errors smaller than radius to the fourth power and clock ratios with errors smaller than radius squared determine the even contractions. The native port readings still require a proper-area calibration. [Directional construction and counts](rate_geometry/static_directions.json) |
| G6. One positive coupling for ordinary and collar matter | The four direction-summed and direction-differenced tests on the 80-cell source exclude even a fitted common rescaling of the proposed sharp-image and bridge-source response. Independent repetitions can distinguish the source preparations. The physical curvature response at the ordinary coupling is an additional measurement. [Signed comparisons](source_implication/jacobi_screen/source_response.json) · [Repeated source readings](source_implication/common_front/velocity/copies/results.json) |

## Reconstructing an Einstein–scalar field on a controlled source

Use three orthogonal radar arms and a local scalar meter on a homogeneous diagonal Bianchi I spacetime. The source evolves the coupled Einstein and Klein–Gordon equations with a quadratic potential. All probes use the same comoving proper clock. Infer the mass squared from the measured scalar equation, the coupling from averaged axial null contractions and the vacuum term from the time-time component at one appointment on one launch. Freeze all three, then predict the equations at three other appointments on two independently initialized launches. Compare radar arm spacings and differencing intervals of 1/20, 1/40 and 1/80. Withhold the evolution state and source coefficients from the decoder. Accept the reference reconstruction if at the finest spacing its largest absolute scalar-equation and Einstein-component errors are both below 0.005 source units and both error envelopes fall at each refinement. Turn off the field's contribution to geometry while retaining the same dynamical scalar and radar readers, and require a held Einstein-component error above 0.05 with the frozen coupling. This is a controlled test of what the observations can identify on homogeneous solutions of a declared action. The field value, proper clock, three orthogonal comoving neighbours and reflection of light are explicit instrument assumptions.

At the finest spacing the observer recovered $m^2=0.699938$, $8\pi G=1.600919$ and $\Lambda=0.119853$, in source units, against the separately recorded inputs 0.7, 1.6 and 0.12. Across six held records the largest Einstein-component error was 0.0002785 and the largest scalar-equation error was 0.0000511. The potential fixes the vacuum energy convention by setting $V(0)=0$. An arbitrary additive potential constant is degenerate with the vacuum term. The control fixes the scalar, its velocity and two initial expansion rates, adjusts the third to the vacuum constraint, and keeps the scalar evolving and gave an Einstein-component error of 0.5884. A secondary forward check initialized the reconstructed equations at the first held appointment on each new launch and forecast the three arm lengths and scalar reading at each of the two later appointments. The largest absolute error in those normalized readings was 0.0000442. The test uses homogeneous diagonal metrics, comoving observers, ideal returned signals and a scalar meter with a declared calibration. It identifies coefficients within that family. Native Interval Echo still requires the corresponding field, stress and geometry instruments on one state. [Source and decoder](observer_field.py) · [Integer clock and scalar records](observer_field.json)

Replay from the repository root with `python code/matched_einstein_limit/observer_field.py`.

## Does measured quantum stress predict independent optical records?

With six quantized scalar modes on the semiclassical box, rotating the occupied modes changes the pressure at fixed initial geometry and energy. One calibration predicts held diagonal curvature to 0.00000559 in box units. A forecast initialized from field quadratures and radar at one appointment predicts later radar ranges to 0.00000999. In the separately quantized Wilson–Maxwell box a prepared photon carries measurable directional stress and leaves the weak optical response invariant. Its geometry is prescribed. [Quantum scalar and radar](quantum_scalar_radar/README.md) · [Photon stress and independent Maxwell probes](quantum_stress_control/README.md)
