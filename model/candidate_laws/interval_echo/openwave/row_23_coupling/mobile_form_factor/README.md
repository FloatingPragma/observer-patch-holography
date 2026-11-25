# Momentum transfer from a mobile material core

An observer can resolve the shrinking size of a traveling two-load core by scattering an address probe. On the tested 80-cell guide, a single measured constituent response and an independent record of the separation between the loads predict all 213 nonforward momentum bins, at rest and in motion. A probe that visits two successive ports before recombination can scatter twice. Its measured recoil resolves this additional contribution. [Protocol](PROTOCOL.md) · [All outcomes](records.json)

## Core and probe

Use the previously enrolled pair of doubled loads on a 214-port closed guide. The two native contact matchings, winding preparation and finite energy filter are those of the [moving-pair experiment](../../row_21_boost/README.md). The pair evolves for 16 periods before the guide is held for its terminal measurement. The two launch orientations have equal prior weight. Each setting enrolls the complete ensemble, including failures at the original source filter. Of the three prior guides only the 80-cell guide has an off-route vacancy beside every interrogated port. The receipt retains the missing addresses on the other two. A ring wavelength is expressed in local guide-port spacings.

The observer supplies a uniformly coherent pointer over the guide addresses. For pointer branch x a controlled native seam exposure of π/32 joins the double or vacancy at x to its vacant off-route neighbour. Sharp reads of all off-route occupancies resolve whether the double spilled. The common all-vacant outcome preserves coherence among the probe paths. A second study repeats the operation at x+1 before recombining the pointer. An early spill ends that launch. A Fourier read gives momentum transfer q=2πj/214 in the surviving runs. The contact strength is the registered value κ=8. The guide schedule, hold, coherent address multiplexing, spill measurement and Fourier transform are supplied apparatus. The simulation uses their declared conditional operations. This is a response of material doubles to an addressed contact. [Native seam gate](../../row_30/experiment.py) · [Source and scattering calculation](measure.py)

## Which part follows from the core size?

An independent bank measures the shorter port separation d. The distance of each equal load from their midpoint is d/2. The rms radii from these [endpoint records](../../row_21_boost/results.json) are **6.588 spacings at rest and 2.376 at winding 13**. The original source filter admits about 5.10% and 2.60% respectively. No source shape is fitted to the scattering counts.

Let a be the no-spill amplitude of the occupied contact. The vacuum amplitude is one. For j≠0 the constant background cancels in the Fourier read. A single sharp double then puts probability g/L² into each bin, with g=|a−1|² and L=214. Summing the independent single-double bank gives g=3.99277. For two loads, tracing over the separately retained material words gives

$$
P_1(j)=\frac{g}{L^2}\,[2+2\langle\cos(qd)\rangle].
$$

Every measured nonforward bin agrees with this forecast within simultaneous 99% sampling allowances. The apparent amplitude relative to a pointlike double-strength source is the square root of the ratio of the measured intensity to 4g/L². Its first resolved 10% fall moves from bin 3 at rest to bin 7 for the narrower moving core. Propagating the size and count intervals places the crossings at **qR between 0.385 and 0.584**, and **between 0.392 and 0.517**, respectively (endpoints rounded for display). Both brackets overlap the earlier core-window edge of 0.5. This test locates the roll-off to a measured interval. The single-load coefficient accounts for the entire one-contact spectrum when the measured core shape is included.

## Repeated scattering in one observation

For material occupations n at x and m at x+1, both zero or one, the no-spill amplitude for a pointer branch is

$$
a^{n+m}=1+(a-1)(n+m)+(a-1)^2 nm .
$$

The additive forecast omits the final term. At the alternating read j=L/2 the two equal visits cancel in the linear term for **every** material word. Only adjacent loads can give the alternating pointer outcome. Its probability per admitted source is g² times the measured adjacency fraction divided by L². At rest the detector counted **1,987,186** alternating outcomes, and in motion **517,322**, each on 2⁴⁰ original launches. The corresponding intervals per admitted source are [3.5355,3.5585] × 10⁻⁵ and [1.7983,1.8213] × 10⁻⁵. The independent predictions are 3.5480 × 10⁻⁵ and 1.8083 × 10⁻⁵.

Over the whole spectrum let A=2+2cos(q). The additive forecast in units of L²P is A g [2+2〈cos(qd)〉]. The exact correction in these units is Pr(d=1)[g²+A g·2Re(a−1)]. The single-load spill probability s determines 2Re(a−1)=−g−Ls. This additional algebraic reduction uses two independent measured quantities in place of the known source gate phase. The multiple-visit forecast agrees with **all 213 bins** in each native preparation. Adding the two one-contact amplitudes alone misses 213 bins at rest and 211 in motion. On the full zero-contrast law the inferred g falls to 0.002409, both alternating counts are zero and neither preparation resolves a failure of additivity at this exposure. The vacant source returns to the forward bin with probability one.

The raw banks are multinomial simulation of the complete terminal law. Each contains every Fourier count, both distinct spill counts and the source failures. Binomial intervals condition only on initial admission and include spills in the denominator. The independent size bank gives Hoeffding bounds for cosine scores in [−1,1] and squared separations in [1,(L/2)²]. A union bound over at most 4,000 tails sets a total 1% sampling allowance. The program retains both exact discrete outcomes and the comparisons. The one-contact shape factorization is an instrument identity. The locations of the roll-off, the incidence of adjacent loads and the resolved size of the double-visit correction depend on the prepared material.

OpenWave asks whether coupling changes with scale with an onset at the core radius. Its M5 precedent credits a bounded classical form factor; its M10 row records an earlier Wilson-dressing comparison. Here a moving core produces momentum-dependent scattering from actual registered contacts. Independent size and single-load controls explain that dependence using one coefficient. A running gauge charge would need an independent interaction that remains scale-dependent after these corrections. The local OW-23 result stays partial. [Current criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L127) · [M10 decision](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L459)

From the RER root, run

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_23_coupling.mobile_form_factor.measure
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling/mobile_form_factor
```
