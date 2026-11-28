# A massive collective wave in commissioned Interval Echo material

Observers can read a moving change in dimer orientation from the loads on each bond while the connecting memories remain coherent and unread. On three held materials, spatial phase measured on even clock holds predicts the profiles on odd holds. The source launch wavelengths fail that same comparison on every board. Two independent local coupling meters forecast all **22** enrolled constituent frequencies on the lattice ruler to within **0.000575** dimer-clock units. The ordinary continuum formula predicts the **14** lines with $0<k\leq0.50$ to within **0.003165**. Distance is counted along the commissioned strand. [Held validation records](../submission/row_30_confirmation/validation.json) · [Comparison protocol](confirmation/TASK.md)

A home observer can read all four possible values of one connector before launch. A second cell distinguishes the two sources by a local wall reading within each recorded connector value. Resolving only the coin halves the source difference at every flight time across that connector. [Connector records and propagation](committed_heads/README.md)

## What travels

A native blank dimer holds one unit of load. Its orientation says which port holds that unit. A second set of collars, written by the registered repair, couples orientations across real seams. An unsatisfied bond is visible wherever adjacent orientations agree in the strand convention. This physical boundary between opposite alternating patterns can travel while each dimer retains its one unit. On each odd-length strand the observers prepared exactly one such boundary, with a spatial phase that they read back on fresh copies. The finite packet has a small admixture of higher odd quasiparticle numbers. [Exact physical reduction and tests](experiment.py) · [Packet implementation](packets.py)

| Board cells | Strand dimers | Frequencies | Maximum lattice error | Low-wave-number lines | Maximum low-wave-number continuum error | Held-profile residual / protocol ceiling |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 23 | 6 | 0.000574677 | 3 | 0.001410236 | 422.56 / 554.07 |
| 80 | 36 | 8 | 0.000551539 | 6 | 0.003164416 | 860.42 / 1137.74 |
| 320 | 30 | 8 | 0.000561597 | 5 | 0.002917395 | 786.58 / 963.37 |

The local meters give a bulk rest frequency $m=15.00051$ and a squared-gradient coefficient $c^2=16.00156$. The lattice prediction uses the nearest-neighbor spatial difference. After subtracting each strand's first frequency to remove the common rest offset, its RMS residual is 4.34%, 8.90% and 8.35% of the matched lattice tangent's residual on the three boards. The absolute-frequency tolerance alone admits the matched tangent. The spatial residual uses all held receiver addresses and scans, with ceilings set by the comparison protocol. [Independent reconstruction of clocks and residuals](../submission/row_30_confirmation/headline.py)

For ordinary wave number $k$ in radians per dimer, the finite-lattice relation bounds the continuum approximation for the same branch,
\[
0\leq \sqrt{m^2+c^2 k^2}-\sqrt{m^2+4c^2\sin^2(k/2)}
\leq \frac{c^2 k^4}{24m}.
\]
Indeed $0\leq k^2-2(1-\cos k)\leq k^4/12$, and the two positive frequencies sum to at least $2m$. In units of distance given by dimer count divided by the independently calibrated speed, the long-wave relation reads $\omega^2=m^2+k_{\rm scaled}^2$. The single-wall frequencies are reconstructed from the observed pair lines. The larger-wave-number continuum discrepancy reaches 0.01664 clock units. [Held wavelength and frequency records](../submission/row_30_confirmation/validation.json)

Across nine independent propagation fixtures, changing the source changes the terminal load records by 18–37 percentage points at cells separated from both source choices. With the contrast coupling removed the largest corresponding change is 0.19 percentage points. The native pair peaks have amplitude at least 0.4994 across the same nine boards. The largest changed-law amplitude at their frequencies is 0.00055. [Propagation records](results.json) · [Matched controls](holdout.json)

At flight contrast 0.5, the same prepared stock and source settings produce remote transport. The measured pair peaks fall below one clock unit, while the native high pair peak is near 30. This isolates dependence of the measured rest scale on the flight coupling. [Propagating changed-law control](propagating_control.json)

A separate phase-kick experiment measures a fastest packet speed of 0.9903 on the 49-dimer strand. Opposite packet phases reverse the signed remote wall-density centroids. The respective largest errors against the group forecast are 0.03427 after normalization and 0.00382 in the unnormalized signed moment. On a three-arm junction, a local flip registers walls together on two outgoing arms. The straight-strand single-wall forecast passes the straight-arm controls and misses the eight exclusive turn readings. [Packet velocities and energy response](packet_speed/README.md) · [Packet centroids](packets.json) · [Junction experiment](fork/README.md)

The even wall parity algebra forbids creating a single odd quasiparticle from an even vacuum with these detectors. An odd open strand admits a physical single-wall preparation. On the 13-, 29- and 49-dimer strands its independently certified one-quasiparticle weights exceed 99.73%, 99.34% and 98.85%. For the tested transverse-field Ising strand, $\omega(k)=2\sqrt{J^2+h^2-2Jh\cos k}$ has maximum group speed $2\min(J,h)$ and squared-gradient coefficient $4Jh$. At $J=8$, $h=1/2$, the maximum group speed is one and the low-wave-number coefficient gives $c=4$, in the same dimer units. Extending the continuum comparison to a three-dimensional field requires a vacuum preparation and calibration of the strand ruler against the ambient distance. [Finite-sector verification](test_collective.py) · [Independent raw-count decoder](decode_example.py)

On two surveyed multicell rings, local measurements on two ensembles recover rest frequencies of 15.000082 and 14.999985. Measuring the individual factors of the source Pauli product erases the global flip sector. The sector correlation is reconstructed by adding measurements of the source product and its product with the signed flip. [Instrument and ring count replay](../submission/row_30_cautions/README.md#the-ring-source)

On a six-dimer preparation both uniform materials pass the spatial read. All 30 mixed first-write histories fail even when each word is supplied to the reader. [Complete conditional census](sharp_growth/README.md)

An independently chosen electric route permits a different free material. Opening every interior contact without dimer guards preserves a 45-mode wall band on the full 44-dimer route. The remote phase contrast falls within \(2.8\times10^{-6}\) of zero when the inter-dimer contacts are cut. On the three measured prefixes every held interior line rejects all three tested endpoint-calibrated relations, including the exact lattice square root. The degree-reference field remains part of that native comparison. The long-wavelength dispersion of this wall requires its own measured gradient and clock. [Released material, direct readout and bound](../row_30_native_wall/README.md)

The held validation replay is specified in the [bound comparison protocol](../submission/row_30_confirmation/README.md#reproduce-from-the-frozen-checkout). The commands below replay the complementary formation, transport, and detector studies from the repository root. An optional output path writes a separate receipt.

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
