# A two-port phase hand crosses the equilibrium source

A vacancy or excess carries a readable launch phase between two disjoint pairs of material ports on the 80-cell equilibrium. At a four-unit flight time the two receiver quadratures recover a reserved $\pi/4$-radian preparation in **[0.7771, 0.7946]** radians for a vacancy and **[0.7497, 0.8242]** for an excess. As the flight duration changes, the received records also distinguish this phase transfer from the ticking of the isolated material hand. Both the all-launch and the arrival-conditioned records test that distinction. [Every simulated terminal count and interval](results.json) · [Native propagation and clock reader](measure.py)

## Preparing, transporting and reading the hand

Hold the same two stationary head memories in the calibrated equilibrium on marks $(0,1)$ and $(38,39)$. The emitter is the blank seam $(5,9)$ and the receiver the disjoint blank seam $(40,45)$. A fixed-charge vacancy or excess occupies one of two rails. On either isolated dimer the exact two-orbital Hamiltonian is $8I+(I-X)/2$. Its symmetric and antisymmetric eigenstates have energies 8 and 9. Starting at the first endpoint and evolving for $\phi$ gives the input density $(I+Z\cos\phi+Y\sin\phi)/2$. The relative phase of the two *energy* modes is the clock reading. Both standards use the same native unit gap. [99-state contact and rest calibration](../PROOF.md) · [Material rails](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/README.md)

Release the dimer and open the entire existing board for $t=2,3,4$ exchange units on separate preparations. At the appointed time isolate the receiver pair. Read $Z$ by terminal occupation, or use a quarter-cycle exchange and the appropriate endpoint sign for $Y$. Nothing measures the individual receiver port before this rotation. The outside outcome has score zero in either basis. The collection window, pulses, phase orientation, source preparation and times are supplied. No phase or source read from the flying load is re-prepared at the receiving dimer. The two original head memories are traced over. The source carries no additional mobile constituent. Earlier endpoint-sign measurements at a single occupied port lose this two-position coherence; the existing five-pixel optical study uses a different input ensemble and tests beam-image response. [Sharp arrival fiber](../two_way_phase/README.md) · [Existing coherent aperture](../../jacobi_screen/coherent_ports/README.md)

Denote the two real unconditional terminal means by $y_\phi(t)$. Opposite preparations at $0,\pi$ and $\pi/2,3\pi/2$ give the two columns of a matrix $M(t)$ by half-differences. Dephasing the input energy modes gives $I/2$ for either source label and independently measures the offset $b(t)$. Quantum linearity predicts
\[
 y_\phi(t)=b(t)+M(t)(\cos\phi,\sin\phi)^\mathsf T .
\]
Every $M(t)$ has nonzero determinant on both probe charges. Inverting it at each preassigned dwell recovers the held $\pi/4$ input. The measured intervals for the two separately acquired erased labels contain zero difference on every quadrature. At zero flight duration the disjoint receiver has no hits in any setting. The isolated rest dimer, measured on the same schedule, gives unit determinant and rotates at its unit energy gap. These measurements use independent ensembles, including every miss. None of the fits uses the held launch angle.

## Transported phase and elapsed time are different measurements

If a single ticking hand has unknown phase advance $\theta(t)$, any fixed nonsingular linear calibrations $A,B$ of its two ends and any positive gain $a_t$ predict $M(t)=a_t B R(\theta(t))A$, where $R$ is a plane rotation. Its orientation, the sign of $\det M(t)$, is fixed even for a varying clock rate. The all-launch data give the following simultaneous 99% ranges.

| Native charge | $\det M(2)$ | $\det M(3)$ |
| --- | ---: | ---: |
| Vacancy | $[-1.130,-1.106]\times10^{-12}$ | $[2.742,2.769]\times10^{-12}$ |
| Excess | $[-1.789,-1.331]\times10^{-17}$ | $[1.036,1.064]\times10^{-15}$ |

Integer propagation and the source enclosure verify both sign changes independently of counts. Continuity therefore gives a time strictly between two and three at which this two-quadrature transfer loses at least one input component. The fixed-gap material clock continues ticking in the rest control.

A detector can also divide each signed mean by its *own* observed probability of arrival, allowing the fraction received to vary with the launch phase, time, and charge. Form $M_c(t)$ from these conditional means by the same opposite-input differences. Its measured determinants stay positive at all three times. For a transferred rotation with fixed calibration, the comparison $L=M_c(t+1)M_c(t)^{-1}$ would satisfy
\[
 \frac{(\operatorname{tr}L)^2}{\det L}=4\cos^2(\theta(t+1)-\theta(t))\le4 .
\]
This follows by conjugating $R(\theta(t+1)-\theta(t))$ by $B$; an arbitrary scalar gain cancels. At $2\to3$ the vacancy counts give **[8.58, 10.15]**. For the excess the two successive equal-duration comparisons give **[3.49, 6.61]** and **[1.11, 1.38]**, excluding any constant phase rate at their joint precision. Independently enclosing each native conditional mean places the first comparison in **[9.2758, 9.4121]** for a vacancy and **[4.6654, 4.7131]** for an excess, strictly above four for both species. Source-dependent diffraction, filtering and mixing of position modes can change these quadratures. The experiment establishes delivery of phase data on the retained source. Registering a moving observer's proper time would require a transported standard whose elapsed phase and trajectory are independently calibrated.

## Error bounds and replay

Each input, dwell, charge and receiving basis uses $2^{48}$ simulated starts under the active model. The mixed incidence at four units is about $2.36\times10^{-5}$ (vacancy) and $1.42\times10^{-7}$ (excess). A joint local port read distinguishes positive, negative and outside outcomes. All three counts are stored. Dividing by the number sent and conditioning on actual arrivals give separate analyses. The isolated control uses its own records. Stable command durations and ideal native shutters are assumed.

The integer half-tick propagator encloses each retained amplitude block within $10^{-19}$ for both left reversals. The reference state's trace-normalized Hermitian centre exceeds $I/50$ blockwise by exact rational principal minors. Its Hilbert--Schmidt error is below $2\times10^{-7}$. Thus each positive event probability differs from its central value by at most $10^{-5}$ relatively. Each of the selected inter-dimer matrix elements has modulus below 0.05, as verified before contracting. Each output array entry sums 64 products of two such elements and one density entry. The remaining two-port products have at most 16 terms. The code adds $5\times10^{-13}$ absolute allowance per event for the block contractions, numerical conversion and finite sampling grid. With $x=\log(2\cdot1024/0.01)/N$, inversion of the Bernoulli Bernstein bound gives a half-width $\sqrt{2\hat p x}+4x$ for an observed event fraction $\hat p$. A union bound covers all stored readouts simultaneously. Interval arithmetic propagates these errors through the two inverses. The tests also compare every nonzero transfer block against a separate direct exponential. [Integer kernel](../../full_algebra/native_precision.py) · [Cholesky-free principal-minor check and independent flight](test_measure.py)

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -B -W error -m matched_einstein_limit.source_implication.material_clocks.transported_modes.measure
python -m pytest -q code/matched_einstein_limit/source_implication/material_clocks/transported_modes
~~~
