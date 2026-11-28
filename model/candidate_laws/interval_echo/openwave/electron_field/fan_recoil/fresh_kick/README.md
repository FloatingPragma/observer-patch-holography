# A material kick on the captured face-0 vacancy

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The same coherent 1045 arrival used in the [fresh reference witness](../../fan_path/local_launch/WITNESS.md) passes the charge-and-face gate on **1.584370%** of original launches. After a local bridge split, equal and opposite material kicks give opposite flow into the terminal fan. Halving one bridge rate leaves the balanced packet's formation mean unchanged and halves its initial recoil. Independent contact calorimetry and a full-flight charge counter both resolve the change. [Fixed protocol](PROTOCOL.md) · [Complete simulated records](results.json)

| Simulated detector statistic | Native bridge | Half-rate bridge |
| --- | ---: | ---: |
| Formation mean in neutral-ring clock units | [9.96222, 10.04464] | [9.96226, 10.04465] |
| Initial current, half-difference for opposite π/2 kicks, in clock units | [3.99263, 4.00274] | [1.99578, 2.00144] |
| Free terminal-owner half-difference divided by native dwell 1/128 | [3.88350, 4.11905] | [1.87777, 2.11334] |
| Local quadratic mass from the π/2 work and lattice phase, in declared units | [0.24480, 0.25576] | [0.48220, 0.52063] |

All measurement ranges are outward rounded from the simultaneous 99% count intervals. Each setting has $2^{30}$ original launches. The simulation records gate failures separately and forms conditional intervals from the actually accepted trial count. The signed motion survives when the face read and the old upstream records are left unobserved after the gate. The phase supplied on a sharply localized arrival before the splitter produces a common phase. A phase on the internal face words is common to the two spatial branches after the splitter. Both controls and the zero-kick panel have zero-compatible initial current. The incoming bridge head itself labels the branch in this sector, so phasing that head has the same action here as a material-port phase.

## Continuity and the matched formation read

The original face gate preserves the two receiving words without resolving them. The arrival program reconstructs all of its Pauli moments on the untouched reference, including the correct vacancy chart, before any recoil operation. The physical bridge 1010–1045 is a collection of exact rate-four swap blocks indexed by the same upstream and fan words. A native quarter exchange followed by the declared correcting port phase puts the two bridge amplitudes in a real equal superposition. The full isolated bridge operator annihilates that balanced packet.

Let $P$ count charge in the receiving owner after release. Only the last bridge contributes to $i[H,P]$. On its paired configurations a quarter exchange followed by sharp charge detection gives a signed gateway-minus-entry mean $Y=2\,\mathrm{Im}\sum_e L_e^*R_e$. Its independent bare transfer rate $r$ sets $\dot P=rY$. A fixed extra relative phase of π/2 gives the real quadrature. These two measured quadratures determine the signed impulse phase. A separate cohort waits under the entire free generator and counts the receiving owner directly. The continuity estimate uses the operator bound $\|[H,[H,P]]\|\le 2\|H\|r$ to bound the finite time approximation.

Independent copies record the odd parity on each physical material link and each written ring by comparing a full isolated exchange against idle. The resulting energy contributions are $2r$ times the odd material weight on a rate-$r$ link and one times the odd ring weight. Both occupied ports have zero degree potential. Summing these counts and subtracting the matched filled-field counts gives the reported formation means. The reference undoes the balanced split and removes the vacancy via the declared local compensated converter. The compartment retains the same fan head words, and the independent parked excess has zero cost in this convention. The calculation checks every included local ring effect against its neutral sector. The extra preparation feeder lies outside the free-flight degree calculation.

The rate change preserves the formation energy because it removes half of a bridge term whose expectation is exactly zero on the balanced packet. It changes the observed current and the independently read kinetic increments. At a held kick of phase $\theta$ the contact-work relation follows from the balanced amplitudes and predicts $r(1-\cos\theta)$. Its predictions from the two quadrature reads and the separate bare transfer meter agree with the calorimeter intervals. The observed probability guards include the archived Taylor bounds and the conditional denominator allowance. The count radius $\sqrt{2\hat p h}+4h$, with $h=\log(2\cdot8192/0.01)/n$, follows by solving the two Bernstein inequalities using variance at most the true Bernoulli probability.

## What the recoil comparison fixes

In the declared one-link coordinate, treating the measured phase as canonical momentum gives two incompatible speeds for a single anchored square-root law. The native π/4 and 3π/4 kicks require disjoint intervals for squared speed, **[37.2158, 43.3094]** and **[32.3353, 33.6313]**. In the lattice coordinate $p_{\rm lat}=2\sin(\theta/2)$, the same work is quadratic and the independently calibrated rate predicts the local quadratic coefficient. This coordinate dependence is expected on a finite bridge. A band inertia inferred from a uniform-field wave number requires separate control of spatial translation and material geometry. [Primary lattice-field comparison](https://journals.aps.org/pr/abstract/10.1103/PhysRev.117.432)

The balanced packet is a supplied, localized zero-current preparation with an energy variance of about 45.30 in native units squared. Its independently measured formation mean provides a controlled reference for this kick. The finite bridge result leaves physical momentum, signal speed, a stationary electron rest energy and a mass-linked intrinsic clock to separate measurements.

From the repository root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_recoil.fresh_kick.measure
~~~
