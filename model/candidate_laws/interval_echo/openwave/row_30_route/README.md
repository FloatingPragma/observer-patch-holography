# Changing the memory route fills the interior gap

On the 320-cell federation, the repeated four-port path has **zero** energies between −5 and −3. The shortest internal routes through the same 194 gateway pairs have **62**. The new modes collectively put **48.28%** of their probability in the middle half of the path, averaged equally over modes. Both preparations still transport a sharp local load to the next echosahedron. The gap measured for the periodic material therefore depends on how the collars connect the gateways. [Protocol](PROTOCOL.md) · [Experiment](experiment.py) · [Full spectra, weights and per-port counts](results.json)

## Spectrum and location

Both experiments keep the same law, hardware, cell order, entrance and exit ports, native degree references and free service on all seams. Each active path collar couples a move to a written coin reversal. Empty writes on other seams close the one-load path sector. A reachable configuration contains the prefix of reversed collars up to the load. Its location on the path specifies a unique physical port. The periodic preparation has four ports in every cell. Breadth-first search for the shortest internal route gives 23 two-port, 152 three-port and 19 four-port cells, including all endpoints. Every edge in each resulting path is a native seam. Energy zero follows the same reduced one-load convention as the enrolled channel, removing the constant contribution of the prepared active collars. [Native column and path checks](../row_30_autonomous/channel.py)

| Spectrum of the complete finite path | Repeated four-port route | Shortest internal route |
| --- | ---: | ---: |
| Native ports | 776 | 578 |
| Energy at the bottom of the spectrum | −13.656816 | −15.796128 |
| Mode immediately below −4 | −5.656854 | −4.002869 |
| Mode at or immediately above −4 | −2.342922 | −3.791681 |
| Modes in the open interval (−5,−3) | 0 | 62 |
| Lowest and highest energy among those modes | — | −4.996547, −3.079896 |
| Sum of their middle-half probabilities | 0 | 29.932762 |
| Average middle-half probability per mode | — | 48.2786% |

The middle half comprises the central 388 or 289 consecutive port centers, respectively. Of the 62 shortest-route modes, 26 put more than half their probability there. The largest probability for one of those modes at any one port is 41.26%. The largest *sum over the 62 distinct modes* at one port is 0.44215, at port 2437. No single mode puts more than 19.38% in the first and last cells combined. The nearest mode above −4 has energy −3.791681 and middle-half weight above 99.9999%. These diagnostics distinguish filling the interval in the interior from adding boundary states. Individual modes may still be localized inside the path.

## Free arrival at the next cell

We used the same entrance port 2506 in the middle cell 208 for both sharp launches. Preparing the corresponding prefix of flipped coins puts the launch in the enrolled one-load sector. The receiving observer covers the next cell 211. The table shows its full terminal occupation, together with simulated detections out of 262,144 trials per setting. All other addresses give a quiet read at that receiver.

| Native exposure | Periodic route probability | Observed arrivals | Shortest route probability | Observed arrivals |
| ---: | ---: | ---: | ---: | ---: |
| 0.5 | 0.05140412 | 13,298 | 0.14622296 | 38,426 |
| 1.0 | 0.08104288 | 21,284 | 0.05980564 | 15,618 |

The two free solvers apply the native tridiagonal generator by spectral decomposition and independent matrix exponential action. They agree on the evolved vectors to better than 10⁻¹⁰. The all-mode spectrum and spatial weights are generator diagnostics. The adjacent receiver counts concern an unrestricted sharp pulse, which contains many energies. They do not establish transport by any specific in-gap mode. This sensitivity study supplies a preparation control for the measured periodic-channel dispersion. [Four-port measurement](../row_30_autonomous/PROTOCOL.md)

From the RER repository root run the following with the model's NumPy and SciPy environment.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_route.experiment
```

The replay compares new results with the complete archived receipt. Use `--output NEW.json` to retain a separate copy.

An independent free control starts with every collar blank and proves that this subspace is invariant for arbitrary loads. Its single-unit angular frequencies are affine in the squared graph wave number on all three registered spheres. [Blank-collar certificate and its scope](BLANK.md)

A uniform local commissioning service also yields exact coordination-dependent bands on six boards. A spatial return measurement on the largest one-write domain rules out a single rest frequency for its constant wave on every board. [Locally written medium and tagged transport](../row_30_formed/README.md)

### Does a remote observer actually receive the in-gap response?

An eigenmode in a spectral interval might reside almost entirely near its source. We tested a sharp launch at a gateway with a frequency filter centred on the original gap and a load detector in a separate cell. The temporal filter uses controlled native free holds, with the trimer as its supplied frequency reference. The registered interaction acts throughout each flight. The filter can reject a launch. We retain all rejections and every terminal cell. For a source in cell 208 and a detector four route cells away in cell 206, the shortest route gave **2,220 joint filter and detector successes in 4,194,304 launches** at the gap centre. The repeated four-port route gave zero successes anywhere on the same number of launches. At the separate passband setting the detector in cell 206 clicked 417 times on the four-port route. [Instrument and exploratory choices](PROTOCOL.md) · [All four-port and shortest-path outcomes](receiver.json)

For a frequency mismatch δ the probability of passing all 8,192 balanced timing reads is cos(δ/16) raised to the power 16,384. The native spectral enclosure is [−16,8] in these units. At the gap centre of −4 the maximum probability of passing for any excitation confined outside (−5,−3) is therefore at most 1.25 × 10⁻¹⁴. This bound also controls the squared norm of any out-of-window contribution before a terminal location read. The centre-frequency response at a distinct cell on the shortest route is physically resolvable. The filter itself permits motion during its coherent holds. These data measure finite-distance transfer in the tested guides. Identifying an extended band at arbitrary distances or a particular group velocity would require a different measurement. [Receiver calculation](receiver.py) · [Independent filter checks](test_receiver.py)

Replay with `python -B -m model.candidate_laws.interval_echo.openwave.row_30_route.receiver`.
