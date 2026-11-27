# How much of the static interaction travels through the material?

Two observers can measure a separation-dependent interaction between charged material probes on Interval Echo. A mobile vacancy can visit both probes even when their local contacts share no port. At three links, observers in 20- and 80-owner federations resolve signed static energies of **−0.0037862 ± 0.0000300** and **−0.0037806 ± 0.0001206** per unit vacancy fugacity at inverse temperature 1/8. The isolated attractive core has an independently measured rms radius of 1.004 links. A single calibrated local interaction and the independently measured propagation of uncharged material predict every near-core thermal comparison. [Ground-state measurements](records.json.gz) · [Dilute-gas measurements](virial.json.gz) · [Held charge-strength comparisons](confirmation.json.gz)

## Why the choice of reservoir matters

The source and probe are guarded loads, with a vacancy denoting material sign −1 and a double denoting +1 relative to an ordinary load. All unguarded contacts and degree-reference terms on the sphere act. A separate ordinary-filled specimen supplies the reference energy at each pair of sign settings. A clock reads the cost of adding one mobile vacancy. Reversing the two signs gives the signed combination of the four readings with coefficients +¼ for equal signs and −¼ for opposite signs. The known direct interaction of adjacent charged guards enters the total energy separately.

In the fixed-inventory experiment a supplied selector prepares the lowest vacancy addition line, mixing levels within 2⁻³⁰ energy units. Local ternary reads locate the mobile vacancy. With only the attractive source mounted, the 99% rms-radius intervals are **[1.00394, 1.00429]** and **[1.00398, 1.00480]** graph links. The two-link signed static energies on the 20- and 80-owner boards are −0.46783 and −0.45297. If the probe changes sign while the neutral-probe vacancy profile is held frozen, independently measured probabilities give −3.99515 and −3.99683. Allowing the material to relax changes the energy. The variational principle bounds every measured energy difference between the work evaluated on the final profile and on the initial profile.

At 32 links on the 80-owner board the signed energy reaches +3.91440. A guard blocking *vacancy exchange* between the two halves preserves +3.91441. Each side has its own addition energy, and a single shared vacancy selects the cheaper side. When the charged neighborhoods are disjoint across the cut, additivity of the two independently stocked side energies gives zero signed interaction. A full native shutter also removes the cross-boundary diagonal energies. At this separation the single-vacancy and independently stocked measurements with the full shutter both agree with zero. For two isolated wells whose addition energies increase under the same sign change, four times the signed work of one shared vacancy equals the overlap length of the two wells’ energy intervals. The result follows by taking the cheaper energy at each sign corner. Resource competition can produce an apparent distant interaction without a vacancy ever crossing the boundary. [Checked finite identities](StaticWork.lean) [Frozen experiment](PROTOCOL.md) · [Acquisition](measure.py)

A separate observer reads the local interaction from the relaxed address frequencies. On a connected stock the unique lowest eigenvector is positive. Its amplitudes follow by taking square roots of the measured probabilities. Sum its eigenvalue equation over the probe's adjacent ports, using the neutral-probe generator at the same source setting. For a resolved line this yields an interval for the local contact energy change as the probe is charged. With an excess at the source and a vacancy at the probe, the 80-owner intervals include −17 exchange units at every sampled distance, from [−17.093,−16.908] at one link to [−17.090,−16.914] at 32. Across the native and zero-contrast controls 73 of 96 inversions have a bounded interval. The other cases have unresolved neighbor amplitudes or admit more than one line. The comparison to the ground-state frozen-profile forecast uses the half-unit spectrum of the registered sharp contact clock. [Count decoder](readback.py)

## What propagates at a fixed chemical potential?

Let z denote the fugacity of a vacancy referred to the ordinary-host addition-energy zero, common to every source sign. At inverse temperature β the grand potential starts with the ordinary energy followed by −z Tr exp(−βH)/β. Here H is the complete native one-vacancy addition generator. This first derivative is exact for the native gas with arbitrary higher vacancy sectors. Finite-density predictions would require the higher coefficients. A fresh uniform mixture of all allowed vacancy addresses interrogates the native generator in two Ramsey quadratures at 192 real-clock durations. A bounded Fourier approximation to the Boltzmann weight turns the complete characteristic-trace records into a thermal energy coefficient. The observer needs the clock counts, normalized inventory size and a spectral enclosure. [Thermal protocol](VIRIAL.md) · [Spectrometer and error bounds](virial.py)

| Distance in links | 20 owners, β = 1/8 | 80 owners, β = 1/8 | Independent neutral-propagation forecast |
| ---: | ---: | ---: | ---: |
| 1 | −3.64209 ± 0.00010 | −3.64211 ± 0.00037 | −3.64210 |
| 2 | −2.05635 ± 0.00010 | −2.05633 ± 0.00037 | −2.05634 |
| 3 | −0.0037862 ± 0.0000300 | −0.0037806 ± 0.0001206 | −0.0037820 |

Each entry is an exchange energy per unit fugacity. The interaction without mobile vacancies contributes to the constant term in z. At three links the neighborhoods of the probes are disjoint. Parking the hopping makes the predicted signed coefficient exactly zero there. The native response at three links therefore detects transmission beyond the common-neighbor geometry. The changed-contrast material and every distance through the farthest port are in the full record. All 48 ordinary-window thermal readings agree with their calibrated predictions; its native three-link intervals alone are too wide to resolve the signal.

Independent single-seam Ramsey measurements give both local signed contact gains, the vacant-port diagonal cost, the hopping amplitude and the private degree energy. The native signed gains are −17.0000000010 and +15.4999999978, each within 0.0000000332 exchange units at the stated simultaneous coverage. These values and the surveyed blank adjacency predict the complete thermal comparison, including repeated visits of the vacancy to either source. The held test uses fresh counts and an energy interval fixed by the row sums of each generator. The same gains predict 12 of 12 readings across two temperatures, three separations and both boards. Multiplying the second probe's two gains by either 0.9 or 1.1 is excluded in eleven settings. Both alternatives are unresolved on the 80-owner board at three links and β = 1/16. All four source settings share calibration uncertainties, and the forecast allowance uses the sum of their absolute bounds. [Held experiment](confirmation.py)

There is an all-distance sign test. In the one-vacancy word basis, write the neutral generator as D−T with D diagonal and the exchange entries of T nonnegative. In the expansion of exp(−βH), a closed walk has a positive neutral weight. Its signed trace difference multiplies that weight by

$$
(e^{-c_+u}-e^{-c_-u})(e^{-c_+v}-e^{-c_-v}),
$$

where u and v are its dwell times on the two charge neighborhoods and c₊, c₋ are the local signed gains. Each factor has the same sign for positive dwell because c₊ exceeds c₋. Walks with zero dwell at either neighborhood cancel. Every connected stock admits a closed walk visiting both, with positive weight on an open set of dwell times. Thus the signed dilute energy coefficient is strictly negative at each finite β on the connected host, and vanishes exactly across disconnected neighborhoods. The finite numerical comparison resolves only the ranges shown by its bounds.

These records distinguish a transmitted response of dilute charged material, its intrinsic local interaction strength and a distant energy from sharing a single vacancy. A ten-percent vertex change at the core would spoil the held forecasts. Identifying a running gauge coupling requires a physical gauge-field propagator and independently identified gauge charges. The full local running-coupling grade is Partial. [OpenWave criterion and M5, M10 decisions](https://github.com/openwave-labs/openwave/blob/main/MODELS.md#simplest-test-per-criterion)

The selector, material guards, phase reference, controlled holds and number ensemble are supplied. Counts are seeded synthetic draws from calculated native-law probabilities. The two-register thermal traces require 2⁴² starts per quadrature, and the high-resolution test uses 2⁵⁰. Sampling and finite-clock allowances cover the stated panels at simultaneous 99% confidence. They quantify measurement resolution for these model calculations.

From the repository root, use one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_23_coupling.relaxed_charge
.venv/bin/python -B -W error -m "$p.measure"
.venv/bin/python -B -W error -m "$p.virial"
.venv/bin/python -B -W error -m "$p.confirmation"
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling/relaxed_charge
```
