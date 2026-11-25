# A massive collective wave in commissioned Interval Echo material

An observer can prepare a change in dimer orientation, read it from the loads on a bond and measure its dispersion while the connecting memories remain coherent and unread. Using one held material seed on three board sizes, independent local coupling readings predicted all 19 enrolled frequencies within **0.000585** dimer-clock units. Opposite packet phases reversed both measured wall-density centroids. The maximum discrepancy from the group forecast is **0.03427** for normalized centroids and **0.00382** for the reported signed, unnormalized density moment. The distance is measured along the commissioned strand. [Held counts](holdout.json) · [Packet records](packets.json) · [Detector audit](../submission/row_30_cautions/README.md)

## What travels

A native blank dimer holds one unit of load. Its orientation says which port holds that unit. A second set of collars, written by the registered repair, couples orientations across real seams. An unsatisfied bond is visible wherever adjacent orientations agree in the strand convention. This physical boundary between opposite alternating patterns can travel while each dimer retains its one unit. On each odd-length strand the observers prepared exactly one such boundary, with a spatial phase that they read back on fresh copies. The finite packet has a small admixture of higher odd quasiparticle numbers. [Exact physical reduction and tests](experiment.py) · [Packet implementation](packets.py)

| Initial material, cells | Held strand, dimers | Maximum frequency error | Square-root / parabolic RMS shape error |
| ---: | ---: | ---: | ---: |
| 20 | 20 | 0.000526 | 0.0228 |
| 80 | 26 | 0.000584 | 0.0318 |
| 320 | 29 | 0.000567 | 0.0480 |

The predictors used the same rest and speed coefficients from independent, local bench meters. The absolute column uses every measured mode. The shape column subtracts each strand's first frequency on both sides to remove the common rest offset. Its ratios remained below 0.077 at the upper 95th-percentile endpoints of a joint, shot-noise-only calibration bootstrap. With ordinary Bloch wave number in the continuum formula, the largest corresponding held errors were 0.01157, 0.01667 and 0.01126. The exact finite relation uses the nearest-neighbor spatial difference. [Frozen test and tolerances](HOLDOUT.md) · [Bench meter](coupling_meter.json) · [Uncertainty calculation](uncertainty.json)

For the ordinary wave number $k$ measured in radians per dimer, eleven of the held lines lie in $0<k\leq 0.50$. Their maximum absolute errors against $\sqrt{m^2+c^2 k^2}$ are 0.001959, 0.002262 and 0.001608 on the three board sizes. This cutoff was applied to the archived counts after the internal validation. Its test on fresh inputs is specified in the proposed registration. The finite-lattice relation gives a bound on the continuum approximation for the same branch,
\[
0\leq \sqrt{m^2+c^2 k^2}-\sqrt{m^2+4c^2\sin^2(k/2)}
\leq \frac{c^2 k^4}{24m}.
\]
Indeed $0\leq k^2-2(1-\cos k)\leq k^4/12$, and the two positive frequencies sum to at least $2m$. In units of distance given by dimer count divided by the independently calibrated speed, the long-wave relation reads $\omega^2=m^2+k_{\rm scaled}^2$. This is a statement about the reconstructed single-wall branch; the directly observed pair lines are sums of its frequencies. [Ordinary-wave-number check using archived observer readings](continuum.py)

The original nine boards supplied 55 modes on strands of 13–49 dimers spanning 3–16 cells. On the independent seed all three held boards passed the locked tests, including predictions of pair frequencies left out of the single-mode reconstruction. Terminal load records in cells separated from both source choices changed by 18–37 percentage points across the original nine preparations. With the contrast coupling removed the largest corresponding change was 0.19 percentage points. The native pair peaks had amplitude at least 0.4994 across the nine boards. The largest changed-law amplitude at their frequencies was 0.00055. [Original records](results.json) · [Same-setting controls](holdout.json)

At flight contrast 0.5, a third arm using the same prepared stock and source settings still propagated to remote cells. Its measured pair peaks moved below one clock unit, while the native high pair peak remained near 30. This separates free transport from evidence for the particular native mass. [Propagating changed-law control](propagating_control.json)

A separate phase-kick experiment measures a fastest packet speed of 0.9903 on the 49-dimer strand. The maximum group speed of the band is one and its low-wave-number coefficient corresponds to a speed of four. On a three-arm junction, a local flip registers walls together on two outgoing arms. The straight-strand single-wall forecast passes the straight-arm controls and misses the eight exclusive turn readings. [Packet velocities and energy response](packet_speed/README.md) · [Junction experiment](fork/README.md)

The positive frequency is reconstructed from observable pair sums. The even wall parity algebra forbids creating a single odd quasiparticle from an even vacuum with these detectors. An odd open strand admits a physical single-wall preparation. On the 13-, 29- and 49-dimer strands its independently certified one-quasiparticle weights exceeded 99.73%, 99.34% and 98.85%. A finite propagating one-dimensional matter wave with an exact lattice gradient and a low-wave-number Klein–Gordon limit is the supported assignment. This experiment leaves open an unconditioned spherical vacuum branch, equality of this strand ruler with the independently ranged ambient distance, and a relativistic three-dimensional field or particle identification. [Finite-sector verification](test_collective.py) · [Independent raw-count decoder](decode_example.py)

Two of the twelve surveyed material boards contain a multicell ring. The original ring calculation preserves a global flip sector when it generates counts. A local measurement of every factor of the first Pauli product erases that sector. The archived values 15.000003 and 14.999997 describe the sector-preserving calculation. On those two exposed rings, a new two-ensemble local measurement recovers 15.000082 and 14.999985. The detector audit includes both the comparison against the original count banks and a replay of the new acquisition. [Original ring calculation](rest/README.md) · [Instrument comparison and revised readout](../submission/row_30_cautions/README.md#the-ring-source)

Run from the repository root. The default commands compare new results with the saved receipts. An optional output path writes a separate receipt.

```sh
P=model.candidate_laws.interval_echo.openwave.row_30_collective
python -m "$P.formation"
python -m "$P.coupling_meter"
python -m "$P.experiment"
python -m "$P.holdout"
python -m "$P.uncertainty"
python -m "$P.propagating_control"
python -m "$P.packets"
python -m "$P.decode_example"
pytest -q model/candidate_laws/interval_echo/openwave/row_30_collective
```

The [compressed sample](example_counts.npz) preserves the positive count in every outcome setting for the first three pair lines and their spatial scans on the first board. The complementary count follows from the recorded exposure. The decoder uses NumPy, the count bank and the registered read instructions alone. Each other raw count array has its own SHA-256 digest in its receipt and can be regenerated from its seed. Use `--sample-counts /path/to/example.npz` with the experiment module to archive the representative bank anew.

A signed spatial impulse reconstructed from even two-time products tests the full lattice wave operator and distinguishes matched hopping. Two separately resolved modes test its low-momentum continuum approximation. [Impulse, uncertainty and records](impulse/README.md)
