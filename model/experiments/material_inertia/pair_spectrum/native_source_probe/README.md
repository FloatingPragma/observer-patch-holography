# A native source reads regional neutrality of the moving pair

A material source built from the native ternary loads now measures coarse signed charge on the same 24-line bound sector after one million exchange units. It replaces the installed phase coupling in the earlier cell-field experiment with explicit native contact, reference and routing words. Across all 60 regions and three sharp launches, the largest **simultaneous 99% count-derived upper bound is 0.1031%** for an uncompensated regional excess or vacancy. The bound applies after the original spectral filter and a recorded meter-site eligibility test. The literal OW-19 criterion, “a stable neutral massive state,” remains locally partial and officially untested: electromagnetic identification and a source-independent physical inertial mass have not been established. [Count-only inference and all intervals](count_readback/README.md) · [Complete source records](results.json)

| Launch owner | Initial band-filter admission | Smallest meter eligibility among the 60 regions, after filtering | Largest native regional charge upper bound | Matched-marginal control at launch cell, lower bound |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 13.2693% | 60.5379% | 0.1031% | 46.02% |
| 5 | 13.2725% | 61.2870% | 0.1005% | 27.23% |
| 10 | 13.2754% | 99.4861% | 0.1003% | 31.63% |

The 60 regions are every graph-owner ball of radius zero, one or two. Every region is measured on fresh starts. No distant-arrival or same-cell flag is used. Original filter refusals, all four nonzero eligibility records and all nine final source words remain in each 14-bin count bank. Each setting starts 67,108,864 launches. Eligibility can exclude substantial population; the table does not assign the same charge bound to that excluded population. Its charged-region mass is separately reported in the source model. The native event-probability error allowance is at most about 4.12 × 10⁻⁶, separately from the sampling intervals.

At the three original launch cells, an additional count-only bound assigns every meter refusal to a possibly charged event. It bounds the **entire spectrally admitted band** below 0.0955% regional charge, with simultaneous 99% coverage across those three cells. This stronger statement concerns those placements; it does not remove the eligibility condition for all 60 regions. [Whole-band inference](count_readback/README.md)

## How the material source measures charge

A coherent neutral source has two arms, loads 20 and 02. Its cross rail successively contacts the targets in a region. The full native refocusing word implements the diagonal seam cost, including the extreme 0/2 block. Relative to a unit-filled target, that cost is

~~~text
V(n,0) − 17/2 = (31/4)q² + (65/4)q
V(n,2) − 17/2 = (31/4)q² − (65/4)q,      q = n − 1.
~~~

At dwell 8π/31, the common q² phase is exactly 2π on every target word. It therefore vanishes even on a coherent data superposition. The remaining source-controlled phases add the region's signed load. A second setting doubles the dwell. Native source preparation, its positive-wait inverse and a native reference wait π/64 supply the X and Y analyzers. The full 27-state pulse replay verifies the action; an abstract Ramsey gate is not substituted for these operations. [Circuit](circuit.py) · [Derivation and observation domain](PROOF.md)

For each region, source rails lie outside it on a real contact path with degree-five reference-service ports. Native endpoint transpositions route all targets to that path and restore every intervening unknown load. An independent permutation check covers all 3,600 target routes. A separate 243-state native propagation checks concrete compensated positive and negative standards, a host, and both internal neutral-pair orientations, including actual Y-reference pulses and two-target interaction. Its largest probability discrepancy from the local source calibration is 1.78 × 10⁻¹⁴. [Physical calibration and routing](physical.py)

The incoming pair may occupy the meter sites. A native signed-contrast echo, followed by its inverse, checks those sites before source insertion. Within the one-excess/one-vacancy domain, zero contrast uniquely means that both sites are unit-filled. A supplied fresh neutral 20 source replaces those known loads, preserving total load. Its preparation, the controller, protected reference timing, final occupation read and programmed extended spatial profile are apparatus. This is a source-induced charge measurement; an autonomous long-range field is not derived. [Protocol](PROTOCOL.md)

## Controls and interpretation

A full-strength hard-core comparator preserves the two one-port marginals while removing their bound correlation. Its launch-cell charged fractions are tens of percent under the same probe, as are those after four additional exchange-only units from the identical arrived native vector. Thus global balance and a zero mean cannot explain the observed small regional response. The independent count inversion reads only X counts and six standard count banks, jointly covering 1,086 banks. It retains random eligibility denominators and calibration uncertainty.

Omitting the extreme-block refocusing gives about 7.78% source leakage on sharp nonzero target charges at the first dwell, and 26.27% at doubled dwell. Host leakage remains numerically zero. This probe control retains the full ternary contact evolution and demonstrates sensitivity to a concrete omitted operation. It does not establish that every allowed probe can be compiled autonomously.

The existing flight and recoil results concern this same stable moving band. This experiment adds a microscopic material-source realization of its small coarse signed-load response. Physical electrical charge, a continuum length standard and physical inertial mass still need independent calibration. It supplies no full OpenWave closure or filing clearance.

## Replay

From the RER root, with its existing virtual environment:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.native_source_probe.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -m pytest -q model/experiments/material_inertia/pair_spectrum/native_source_probe/test_probe.py
.venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.native_source_probe.count_readback.measure
~~~

The headline replay recomputes the 57,360-configuration source sector, all regional laws, native pulse calibration, routing and seeded counts, then compares the read-only receipt. Use `--output /tmp/new-native-source.json` to acquire a new exclusive output instead. The five tests check the exact cost decomposition and stroboscope, actual native Y/inverse, all routes and compensated standards, complete counts with full-strength controls, and rejection under optimized Python. The independent count-only replay and tests are documented in its own [README](count_readback/README.md).

A three-site meter gives every spectrally admitted launch a direct signed read. The operator bound on separated owners covers the entire 24-mode spectral cluster throughout free evolution. [Whole-band regional charge](regional_census/README.md)
