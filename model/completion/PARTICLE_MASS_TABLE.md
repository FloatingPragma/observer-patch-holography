# OPH mass and energy readout

The table separates proved continuum results, supplied quadratic designs and finite excitation gaps. **No measured particle species or calibrated physical mass has been derived.**

| Sector | Result | Status |
| --- | --- | --- |
| [Canonical free continuum](../ballistic_free_boundary/README.md) | Eight Weyl valleys, E=v·norm(p), rest energy 0 | Proved when U=o(c^(2/7)); no species assignment |
| [Supplied quadratic masses](../native_mass_channels/README.md) | E_rest=abs(m0−norm(m)) and abs(m0+norm(m)) | Four positive and four negative levels per generic branch; coefficients are inputs |
| [Actual ordered-spin continuum](../spin_canonical_energy/README.md) | Nonzero positive energy bands at every small momentum | Sharp rest mass and dispersion remain undetermined |

For a relativistic branch with limiting speed v, mass would be E_rest/v². Finite lattice gaps and a zero spectral threshold alone do not establish that relation. No conversion to eV or MeV is supplied.

## Exact coefficient designs

These are explicitly chosen quadratic Hamiltonians h+B, not predictions or fits. Their formula is not the spectrum of the U=1 interacting centered cube below.

| (m0,m1,m2,m3) | Rest energies | Larger/smaller |
| --- | --- | --- |
| (0, 0, 0, 1) | 1, 1 | 1 |
| (2, 0, 0, 1) | 1, 3 | 3 |
| (9/4, 5/4, 0, 0) | 1, 7/2 | 7/2 |

Choosing m0=(A+B)/2 and norm(m)=(A−B)/2 produces any chosen pair A≥B≥0. In particular, m0=(q+1)/2 and m1=(q−1)/2 produce energies 1 and q for any q≥1. The free coefficients therefore do not predict a mass ratio. This single uniform family has at most two distinct positive rest energies, so a direct assignment cannot reproduce three distinct charged-lepton masses. This statement does not exclude other OPH sectors. The ordered spin expectation has not been identified with these coefficients by a self-consistency proof.

## Centered finite charged spectrum

The first six of 228 visible lines per local spin/charge branch follow. Reference: unique centered N8 ground, K=H−N/2, c=1/8, U=1, open 2×2×2 cube. Intervals use the dimensionless source clock; all nonzero residues remain in the [full 228-line table](../centered_particle_spectrum/lines.csv).

| Line | Excitation energy interval | Relative to first |
| --- | --- | --- |
| 0 | [0.36168271180597, 0.36168271180599] | [1.000000000000, 1.000000000000] |
| 1 | [0.40229304405600, 0.40229304405602] | [1.112281651636, 1.112281651637] |
| 2 | [0.41517426333336, 0.41517426333338] | [1.147896346110, 1.147896346111] |
| 3 | [0.43819649680927, 0.43819649680929] | [1.211549467269, 1.211549467270] |
| 4 | [0.46219220757884, 0.46219220757887] | [1.277894111308, 1.277894111309] |
| 5 | [0.47442977795619, 0.47442977795621] | [1.311729210354, 1.311729210355] |

Treating these finite gaps as charged-lepton rest energies under one positive scale is excluded by the retained ratio certificate. This excludes that dictionary, not all continuum or composite sectors.

## Charged composite spectrum

The same centered ground has 138 positive lines in a cubic probe with hidden spin 3/2 and eta 1/2. This response is orthogonal to every linear CAR response above. Every weight is retained in the [full composite table](../charged_composite_search/lines.csv). The block does not exhaust all charged composites, and hidden spin has no Lorentz-spin assignment.

| Line | Excitation energy interval | Relative to first |
| --- | --- | --- |
| 0 | [0.40246851577010, 0.40246851577013] | [1.000000000000, 1.000000000000] |
| 1 | [0.42208975632848, 0.42208975632850] | [1.048752237229, 1.048752237230] |
| 2 | [0.43899074292156, 0.43899074292158] | [1.090745550820, 1.090745550821] |
| 3 | [0.47993821750041, 0.47993821750043] | [1.192486365255, 1.192486365256] |
| 4 | [0.48578552247412, 0.48578552247414] | [1.207014967480, 1.207014967481] |
| 5 | [0.49403565036905, 0.49403565036907] | [1.227513783093, 1.227513783094] |

Independent source calculations of C* C and C* [K,C] enclose total weight in [0.74354609240642, 0.74354609241091] and conditional mean excitation energy in [0.52537938729092, 0.52537938735399]. These source-clock energies have no physical mass or SI calibration.

## Finite ordered-spin reference

Here J=1 and B/J=1/4 on the open eight-site cube. The [complete certified 256-state readout](ordered_spin_energy_readout.json) has 68 separated enclosure groups; grouping does not assert exact degeneracy. This is a positive-field finite reference, not the infinite zero-field ground.

| Positive group | Excitation energy interval | Relative to first |
| --- | --- | --- |
| 1 | [1.092284677676, 1.092284677698] | [1.000000000000, 1.000000000019] |
| 2 | [1.442323754758, 1.442323754779] | [1.320465062092, 1.320465062137] |
| 3 | [2.166719725499, 2.166719725521] | [1.983658445221, 1.983658445278] |
| 4 | [2.543094000272, 2.543094000294] | [2.328233703355, 2.328233703419] |
| 5 | [2.665699560833, 2.665699560855] | [2.440480595638, 2.440480595704] |
| 6 | [2.729285256683, 2.729285256705] | [2.498694078944, 2.498694079012] |

The [filtered spin-packet readout](../spin_packet_readout/README.md) additionally gives certified conditional excitation energies and success probabilities for this same finite reference, including all rejection outcomes. Its continuum packet construction has an explicit finite resource bound; it does not determine a relativistic rest mass.

## Reproduce

```sh
python -W error tools/build_particle_mass_table.py --check
python -W error -m pytest -q tools/test_particle_mass_table.py
```

The [machine-readable table](particle_mass_table.json) retains all 228 linear charged lines and 138 composite lines, the exact coefficient examples and source hashes. Its builder checks reviewed proof outputs and performs exact rational arithmetic; it does not rerun every inherited proof producer.
