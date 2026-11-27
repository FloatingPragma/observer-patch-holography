# A pole in the neutral material's two-point response

A local flip on the polar ring has a positive spectral line at each of its independently timed one-island energies. On the ten- and twelve-dimer rings at native contrast, these isolated lines together carry **at least 99.0%** of the unnormalized local flip's unit spectral weight, with simultaneous 99% confidence. Their measured strengths and earlier frequencies forecast 28 signed responses at other owners within 0.00023 in probability. Free transit and local readings use the same commissioned material as the neutral mass and recoil experiment. The Hilbert space and Born instrument belong to the candidate. [Prior physical circuits, emission and charge](../README.md) · [Frozen two-point test](PROTOCOL.md)

| Dimers and contrast | Lowest certified N times a momentum residue | Largest error on opposite-kick response | Largest error assigning all modes the rest frequency |
| --- | ---: | ---: | ---: |
| 10, native | 0.9868 | 0.00023 | 0.55763 |
| 12, native | 0.9857 | 0.00017 | 0.37361 |
| 10, half contrast | 0.9588 | 0.00090 | 0.71957 |
| 12, half contrast | 0.9586 | 0.00122 | 0.24071 |

The residue is the probability mass of one isolated energy in one momentum sector of the ground two-point spectral measure. The normalization in the table makes a sharp local flip distributed equally over N momenta have weight one in each column of a momentum scan. The observer checks a joint spatial and temporal output after a single flip in coherent superposition over address and time. Every preparation refusal, spatial output and off-frequency exit is counted. On fresh copies the emitter gives the origin one of two opposite local quarter turns. An observer two or half a ring of dimers away reads its own dimer's exchange parity. Subtracting the two probabilities isolates the imaginary part of the vacuum correlation. Switching off the dimer exchanges retains the static interactions and eliminates the remote response, to within 0.00024 in all 56 controls. The exact distant commutator then vanishes. [Complete count banks](records.json) · [Derivation, apertures and assumptions](PROOF.md)

The same pole counts bound what an observer can see after arbitrarily long waits on a fixed ring. The connected local return amplitude exceeds 0.31 in magnitude at arbitrarily late times on the ten-dimer ring and 0.28 on the twelve-dimer ring, at the coverage of the pole measurements. These are forecasts from the finite spectrum. Establishing a scattering limit requires control as the material grows. [Long-time return bound](PROOF.md#what-persists-at-arbitrarily-late-times)

## When do two separate flips look like two independent arrivals?

Fresh emitters flip two dimers simultaneously. The clock reads ordinary energy relative to the initial all-lower word. At opposite ends, 98.49% and 98.55% of admitted native launches lie within one unit of twice the previously timed rest frequency on ten and twelve dimers. For adjacent launches, 98.18% and 97.85% instead enter the window about 48 exchange units, the adjacent-word energy. At half contrast the far double-energy fractions exceed 95.45% and the adjacent lower-energy fractions exceed 93.28%. The maximal false acceptance from energies more than two units outside the respective center is 0.00644. All rates retain separate cold-refusal counts. [Prospective pair test](PAIR_PROTOCOL.md) · [Energy and spatial histograms](pairs.json)

At times 0 through 32 on both native rings, two independent single-flip tapes forecast the opposite-end pair's unordered address histogram within 0.021 in total variation. At 64 the discrepancies grow to 0.502 and 0.272. The adjacent preparations already miss the tolerance of 0.05 at time 8 on both lengths. The two-site and half-ring launches also appear separately in the tapes. At half contrast the first failing appointments for separations one, two and half the ring are 8, 2 and 32; at native contrast they are 8, 8 and 64 on both ring lengths. A deviation from this distinguishable independent-arrival forecast can include exchange interference or interactions. These records give a finite, separated double-energy channel and limit how long independent arrivals describe it. Infinite-volume pole persistence, a relativistic spectral measure and asymptotic scattering in PR-97 require further evidence. [Independent-tape decoder](pairs.py) · [Finite residue lemma](PoleWeight.lean)

From the repository root, set the numerical thread count to one and run the two sibling modules and their focused tests.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.experiments.material_inertia.pair_spectrum.neutral_guide.polar_ring.two_point
python -m "$p.measure"
python -m "$p.pairs"
python -m pytest -q model/experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/two_point
~~~
