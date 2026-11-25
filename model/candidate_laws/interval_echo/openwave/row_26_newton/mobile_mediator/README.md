# One mobile mediates an orientation-sensitive interaction between separated neutral sources

On the 20-cell Interval Echo federation, a freely moving excess or vacancy can visit the neighborhoods of two separated, clamped excess–vacancy dipoles. Those joint visits change a connected Ramsey spectrum. In a sufficiently dilute equilibrium preparation the sources have an attractive free-energy ordering across four noncontact cell separations. Interchanging the excess and vacancy of the second dipole reverses that ordering. This sign reversal is independent of the charge of the mediator. [Preparation and readout](PROTOCOL.md) · [All 40 source panels](results.json)

The table uses an excess mediator and the two gateway orientations. The thermal entries are the connected one-mobile heat trace at inverse native temperature `β=1/8`. The first-order interaction potential per unit mobile activity has the opposite sign and equals eight times the magnitude. The lower and upper trace endpoints agree to the digits shown.

| Cell distance | Mobile hops between affected neighborhoods | First mixed moment order | Thermal trace, parallel gateways | Thermal trace, second dipole reversed |
|---:|---:|---:|---:|---:|
| 2 | 2 | 6 | +1.03487 × 10⁻⁹ | −2.47557 × 10⁻⁹ |
| 3 | 5 | 12 | +9.92460 × 10⁻²² | −1.56269 × 10⁻²¹ |
| 4 | 8 | 18 | +2.92557 × 10⁻³⁵ | −4.06017 × 10⁻³⁵ |
| 5 | 10 | 22 | +3.60509 × 10⁻⁴⁵ | −4.78143 × 10⁻⁴⁵ |

At separation two the ordinary return meter resolves mediation without a temperature or activity assumption. At wait eight the connected real Ramsey amplitude is −0.00940. The inverted pointer records give −0.00942 with a simultaneous 99% radius of 0.00119 over all 640 binary banks. At separations three through five the same real pointer at wait eight leaves zero within its counting interval. The tables include the unresolved readings. The ground-energy sum also depends on separation. With an excess mediator and parallel gateways it is +16.05491, +16.05497, +1.13213 and +16.05580 at the four noncontact distances. Even a completely immobile load gives +17, +17, +1.5 and +17 through its choice of lowest orbital on each arm. A ground-state subtraction alone therefore cannot establish transport. [Complete counts and spectral intervals](results.json)

## Native generator and return measurement

The clamped footprints stay the same in the four source settings 00, 10, 01 and 11 at each separation. Relative to its own stationary zero-mobile background the source-conditioned one-mobile Hamiltonian is

\[
 H_{ij}=H_0+iV_A+jV_B,\qquad H_0=\operatorname{diag}
 \left[(\kappa+\tfrac12)d+\kappa(6-d)((1+q)^2-1)\right]-\tfrac12 A .
\]

Here `q=±1` is the mobile's excess relative to unit load, `d` is its original port degree and `A` is adjacency on the 236 available mobile positions. The local diagonal on a seam is `e(l,r)=κ(l−r)²+1[l≠r]/2`. Activating a pin of load `s` adds `e(1+q,s)−e(1,s)−e(1+q,1)+e(1,1)` at each unpinned neighbor, with no change elsewhere. Applying the full native seam and degree terms to every allowed occupation independently checks these columns. The off-diagonal exchange between unpinned neighbors is −1/2. Blank memories are invariant. The two neutral sources have individual sharp constituents in this test. The earlier pair-recoil preparation allows coherent internal motion throughout each source cell, with an empty intervening mobile sector. [Native contrast and seam action](../../../../contrast_memory_spherical_v1/law.py) · [Earlier confined pairs](../../../../../experiments/material_inertia/pair_recoil/FORCE_PROTOCOL.md)

Let `Δf=f₁₁−f₁₀−f₀₁+f₀₀`. A uniformly sampled initial port, with a coherent evolved and reference arm and a final quadrature pointer, reads the real or imaginary part of `tr exp(−itHᵢⱼ)/236`. The probability for the plus outcome is `(1+component)/2`. The stationary background phase is referenced separately on each arm. This protocol supplies the conditional Ramsey control and an exchange-time standard. It measures a response of the mobile return spectrum to the two sources. A local elapsed proper time or null signal would require its own material clock and propagation readout.

The connected spectrum directly detects transport. The two diagonal perturbations have disjoint supports in every enrolled geometry. If the off-diagonal hopping is removed, each mobile position feels at most one source and the mixed trace vanishes identically for every wait and temperature. With hopping restored, expand the trace as a sum over closed walks. A mixed term must include both potentials. If their supports are `dₘ` hops apart, its earliest order is `n=2dₘ+2`, with one insertion at each end and a shortest walk each way. If `Nₐᵦ` is the exact count of shortest walks, cyclicity gives

\[
 \Delta\operatorname{tr} H^n
 =\frac{n}{2^n}\sum_{a,b:\,d(a,b)=d_m}(2V_A(a))(2V_B(b))N_{ab}^2.
\]

Every smaller moment cancels. The exact integer sums at distances two through five are 3844, 15376, 61504 and 61504 for parallel gateways and an excess mobile. Reversing the second dipole gives −4216, −16864, −67456 and −67456. The trace identity survives subtraction of any source-independent energy zero. At cell distance one the clamped gateways block the direct route and the first mixed order is 20. This nearby geometry is included in the data. [Walk-count algorithm](measure.py) · [Independent native matrix-power check](test_measure.py)

## Static potential on the specified ensemble

Take an equilibrium ensemble containing zero or one mobile of the chosen charge. The dimensionless activity `z` and inverse temperature `β` are common to every source arm. After removing each static background, `Zᵢⱼ=tr exp(−βHᵢⱼ)`, and the connected mobile contribution is

\[
 \Omega_{\rm conn}(z)=-\beta^{-1}\Delta\log(1+z Z_{ij}).
\]

Thus `Ω′conn(0)=−ΔZ/β`. On the κ=8 board all four Hamiltonians have nonnegative spectrum, as the row-sum lower bounds verify, hence `Zᵢⱼ≤236`. Using `0≤x−log(1+x)≤x²/2` gives the finite remainder

\[
 \left|\Omega_{\rm conn}(z)/z+\Delta Z/\beta\right|
     \le 2z\,236^2/\beta .
\]

The receipt takes `z=10⁻⁵⁵` so that even the smallest signed gap in the table survives this bound. At that single positive activity the parallel-gateway potential is lower at each shorter noncontact separation; the reversed-gateway potential is higher at each shorter separation. This establishes a finite dilute equilibrium comparison conditional on the specified ensemble and clamp placements. The count experiment uses one prepared mobile and does not sample the dilute ensemble. Other orientations and placements of the dipoles on the nonuniform board require independent readings.

An equal-energy control uses two adjacent degree-five ports behind each gateway. The first port is the least-numbered regular neighbor of the gateway. The second is its least-numbered regular neighbor in the same cell that is not adjacent to the gateway. Exchanging the two constituent loads then keeps the isolated static energy of the target exactly equal to the source's 116.5 exchange units. Directly summing the eight unit-contrast seams, one contrast-two seam and the two reference units proves this number in both orientations. With these pins the exact first mixed orders on the four noncontact shells are 10, 16, 22 and 26. At `β=1/8` and `z=10⁻⁷²`, reversing the target again reverses the ordering of the finite potential at each adjacent pair of enrolled separations, for both mobile charges. Thus a difference between the two *static* source energies is unnecessary for the sign change. [Equal-energy source and certified readings](matched_energy.json) · [Construction](matched_energy.py)

The thermal interval computation propagates a triangular matrix of four blocks. With `Mᵢⱼ=β(128I−Hᵢⱼ)`, start the baseline block at the identity and the two single differences and mixed difference at zero. The mixed derivative block has homogeneous generator `M₁₁` and sources `−βV_B` times the first single difference and `−βV_A` times the second. The mixed trace is computed without subtracting four large nearly equal traces. A second recurrence with absolute entries bounds the total unsigned path weight. The maximum row sum of that generator is at most 24 and each row has at most nine entries. The first omitted order after 180 terms and its geometric tail bound the truncation error. [Outward numerical enclosure](measure.py)

For binary64 round-to-nearest arithmetic, a path through 180 sparse multiplies/divides, the 236-entry diagonal trace and the final sum incurs a relative absolute-majorant factor smaller than 2⁻⁴⁰. This follows from `γₖ=ku/(1−ku)` with unit roundoff `u=2⁻⁵³` and `k≤20·180+236+180`. The implementation allows 2⁻³⁹ times the computed majorant to cover both recurrences. Even allowing a smallest-normal absolute rounding error at every elementary operation, the total propagated underflow contribution is bounded by `236·180²·40·exp(25)·2⁻¹⁰²² < 10⁻²⁴⁰`. The remaining common exponential factor is enclosed by a rational Taylor bound. Static energies use a full real eigensystem. The measured Frobenius residual for `HQ−QD`, with a dot-product rounding allowance and an independent lower bound on the least singular value of `Q`, bounds the lowest eigenvalue error.

## What survives averaging over unpolarized sources?

Independently reverse each equal-energy dipole with probability one half. The mobile samples either excess or vacancy with probability one half. Each clamped source still has energy 116.5 and zero net load contrast. In units of the doubled Hamiltonian, the two connected single-pin contrasts seen by an excess are 31 and −34. A vacancy interchanges them. A mobile position adjacent to c constituents of the first dipole therefore sees mean contrast −3c/2. The average for the second dipole is independent.

On any connected allowed mobile graph with separated source neighborhoods, write cₐ and cᵦ for the numbers of adjacent constituents at a pair of vertices, and Nₐᵦ for the number of shortest mobile paths joining them. At the earliest mixed order, the orientation-averaged signed path weight is exactly

$$
  \frac94\sum_{a,b:\,d(a,b)=d_m} c_a c_b N_{ab}^2 > 0 .
$$

The first nonzero order is even. The mixed heat trace is thus positive and its leading dilute potential negative at sufficiently small positive inverse temperature on each such finite graph. This follows from the native contrast and the propagating return walks. The random orientation preparation supplies the uniform weights. The integer coefficients in the four equal-energy geometries are 9, 36, 144 and 720.

At the specified finite β=1/8 the eight-arm mean heat traces at cell distances 2, 3, 4 and 5 are approximately 1.31884 × 10⁻¹⁸, 4.40746 × 10⁻³², 1.76715 × 10⁻⁴⁶ and 1.07323 × 10⁻⁵⁶. The connected potential increases towards zero with separation for the unpolarized source. This ordering holds both for averaging the conditional connected free energies and for equilibrating the fair degenerate internal states first. In the second case the logarithm acts on the averaged heat trace. Both prescriptions have the same leading trace and the same absolute finite-activity error bound above, with z=10⁻⁷². The signals quoted here are population predictions on the specified dilute ensemble. [Eight-arm enclosures](unoriented.json) · [Native calculation](unoriented.py)

The additive zero of potential can be removed without fitting it. If U(r)=u−C/r in these station units, the finite difference Cᵣ=r(r+1)[U(r+1)−U(r)] must have the same value at each r. Using U per unit activity gives disjoint intervals for C₂, C₃ and C₄, respectively 6.33041–6.33044 × 10⁻¹⁷, 4.23115–4.23118 × 10⁻³⁰ and 2.82742–2.82745 × 10⁻⁴⁴. Random orientations recover attraction on these preparations. They leave the scale and radial dependence required for Newton gravity unresolved. A graph crossing count and clamped source energy would each require independent calibration as proper distance and inertial mass. [Independent path and potential tests](test_unoriented.py)

## Gravity scope

For κ=0 with the same alphabet, swaps, pinning and board, both orientations have positive mixed thermal traces. The reversal effect therefore tests the native contrast law. On κ=8, two equally neutral dipoles of the same constituent inventory can reverse the sign of the measured mobile-mediated contribution. Their independently summed static source energies above the filled host are +133 and +101 native units. The pins occupy ports of different degree, so the two orientations have different energies, both positive. A calibration of free inertial masses would require released sources. The thermal interaction decays across the enrolled noncontact separations much faster than a Newton potential. A universal attractive mass interaction, a measured inverse-square strength, and a gravitoelectromagnetic limit would require additional source and field results. The measured Ramsey signal also leaves the common clock and light geometry required by OW-27 open. The reference code and active model were inspected at RER `3cc56709a8aa29222b3287de33c06eab8408687a`; the official criteria were read at OpenWave `0df73dc6be3bde9f0b8ced6044effc79f0d92bd6`. [Official comparison criteria](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#simplest-test-per-criterion) · [Existing clock comparison](../../row_27_metric/free_source_clock/README.md)

Run from the RER repository root with the model dependencies.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.measure --write
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.matched_energy --write
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.unoriented --write
python -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_26_newton/mobile_mediator
```
