# A material clock sends its own frequency to another station

On a stationary source in the 20-owner Interval Echo simulation, two clocks made of unit loads report a **near-to-far frequency ratio of 32.01586** and a **far-to-near ratio of 0.0312348**. Each station receives the other station's clock modulation by detecting the load that actually moved in that clock. Double loads using the same bonds, collar, guide, receiving positions and appointments report a ratio of one in both directions. A common gravitational lapse for these two calibrated materials fails. [Prospective settings](PROTOCOL.md) · [All finite probabilities, seeded counts and readings](results.json.gz)

| Station pair | Unit receiving ratio | Simultaneous 99% interval | Double receiving ratio | Simultaneous 99% interval |
| --- | ---: | ---: | ---: | ---: |
| Near emitter, far receiver | 32.01586 | [31.7960, 32.2362] | 1.000007 | [0.99005, 1.01006] |
| Far emitter, near receiver | 0.0312348 | [0.0310597, 0.0314112] | 1.000003 | [0.99005, 1.01006] |

Each entry divides the remotely recovered modulation frequency by the receiving station's locally read clock frequency. Each histogram contains **1,024 appointments with 16,777,216 simulated starts apiece**. Of the 17,179,869,184 starts in the near-to-far unit arm, 1,501,746 recorded arrival at the destination. The weak signal includes the low success probability of the fast near hand. In the reverse direction 1,531,829,017 starts recorded arrival. Closing exchange only during transport eliminates every remote arrival in all sixteen parked links, while the sending hand still moves. Reversing the stationary source's coin at fixed separated energy returns both materials to unit transfer. A zero-direction sham and zero contrast also give unit transfer. [Causal prediction and stationary-transfer identity](PROOF.md)

## What the stations do

Two guarded source contacts join a parked unit at the pole to belt positions 2 and 8. For the ordinary source, split a sharp double between the pole and a source endpoint by the native partial exchange, set its relative phase with the pole's private reference and write the standing pair with one native echo. The outgoing material word is sharp $(1,1)$ and the coin-even signed memory remains on the seam. Routing swaps exchange these sharp singles with shuttered stock while forming the second write. The complete route list retains the two protected core seams. Applying a supplied relative coin phase during formation gives the equal-energy coin-odd control. The source is stationary through hand preparation and flight. Its separated contacts each have energy 32 in both coin preparations. The source preparation, shutters, sharp echo read, local storage, independent clock replicas and appointment calendar are supplied. [Gate and preparation checks](native.py)

The near hand runs on belt bond 2–3 and the far hand on 4–5. Each runs from a sharp load on its first site. A native echo on the bond reads its signed difference and the inverse restores the word. On a bright read that same load starts at the second site and propagates on all ten links. A dark read makes a silent trial. A fair unsignaled choice between two flight durations broadens the remote signal. At reception an observer can recover the final position from nine signed reads along a path of blank perimeter seams. The source memories are on the two pole contacts. The complete receiver alphabet contains the silent outcome and ten distinct port outcomes.

The observer fits a constant and two quadratures, with free frequency, on the first 512 appointments. The second 512 test the frozen fits. Independent sharp-flight calibrations and the measured sending-clock response predict the whole held modulation. Across all sixteen open links, the largest discrepancy against the prediction using independent local clock counts at the two emission epochs was **0.379 of its simultaneous error allowance**. Using only the prefix fit of the sending hand and the impulse calibration, the largest held absolute error on the weak unit near-to-far link was $9.11\times10^{-6}$ in arrival probability. A reference dimer reads one externally supplied appointment step as 0.062317 exchange units with interval [0.060724, 0.063871]; the nominal value is 1/16. The frequency ratios cancel this common scale. All stochastic data are PCG64 multinomial samples from simulated Born probabilities. [Count-only reader](readout.py) · [Independent raw-contact propagation tests](test_link.py)

## OW-27

The known local near-clock detuning predicts a reciprocal rate change only if the signal reaches the receiver with measurable modulation. A stationary two-delay guide preserves any surviving tone and can suppress its amplitude. The simulation resolves that tone, tests a forecast from other records, and rejects a universal interpretation using the second tag. The full OW-27 criterion remains **Open**. These records carry a material intensity modulation. A gravitational redshift or Shapiro delay still requires independently characterized optical propagation and exterior proper clocks on the same gravitational source. The semiclassical scalar comparison supplies proper clocks and null rays on an imposed metric sector. [Original belt and source census](../README.md) · [Matched scalar and geometric radar](../../../../../../code/matched_einstein_limit/quantum_scalar_radar/README.md)

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_27_metric.frequency_link
nice -n 15 python -B -O -W error -m "$p.measure"
nice -n 15 python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_27_metric/frequency_link
```
