# What can a material probe outside the moving quark candidate measure?

On the registered Interval Echo seam, an observer in an adjacent cell can detect a freely launched period-three wall with a stationary probe. The signal changes when the preparation changes, including when three loads rearrange into a compact cluster. Two independent probe signs and a unit calibration give different charge estimates. This experiment establishes native material coupling across a real cell boundary. OW‑12's fractional *electric* charge remains unassigned.

## What the observers measured

The source is the [ten-port belt wall](../quark_winding/belt_walls/README.md). Three excess loads on a host with record 1 make a period-three word with one gap of four ports. Four loads make a word with one gap of one port. Both produce a one-third-turn change of density-wave phase. For each wall we launched equal amplitudes on two sharp positions with relative phase +i or −i, let the belt evolve unaided for one exchange unit, and read all ten wall positions plus an unclassified bin. Reversing the launch took the probability at position 3 from 0.7201 to 0.0678 for the three-load source, and from 0.7179 to 0.0600 for the four-load source. The probes were prepared independently as record 0 or 2 in a medium of record 1.

**Stationary exterior probe.** The source cell's belt faces gateway port 3; its bridge leads to port 48 in the adjacent cell. The observer holds a probe at 48 in one interferometer arm and puts it in matched storage in the other. The five seams between the belt and port 3 remain open during readout. No test material crosses into the source. The table gives the phase after two exchange units for strength 8 and the registered record-zero degree correction. All entries are exact finite predictions rounded to five places. Separate detector counts and joint belt/gateway records accompany each one. [Registered exterior measurements](exterior_native.json)

| Source | Negative probe, radians | Positive probe, radians | Gateway occupied in the reference arm |
| --- | ---: | ---: | ---: |
| Empty host | 0 | 0 | 0% |
| Single excess unit | +0.09274 | +0.01103 | 0.111% |
| Three-load wall, after free flight | −0.02279 | +0.01086 | 0.160% |
| Three-load compact control | −0.01569 | +0.28032 | 0.592% |
| Four-load wall, after free flight | +0.10409 | +0.02084 | 0.255% |
| Four-load compact control | −0.24344 | −0.00136 | 26.929% |

The negative-probe sine quadrature on the two moving walls had count differences −4,530 and +17,530 out of 262,144 shots. The empty control gave −98. Blocking all five source-to-gateway exchanges while keeping their native diagonal energies removes the exact phase for either wall. A unit prepared on the gateway still gives the predicted signed phase in the same held instrument. The signal therefore has a material route. Gateway occupation at the terminal reading does not specify how long material visited the gateway during the experiment. An experimental idle-slot reference centered at record 1 also produces an exterior signal with different amplitudes. [Experimental reference service](exterior.json)

**Common and differential polar probes.** One test qutrit at each pole covers all ten belt sites equally. This common instrument returns exactly three or four unit responses for the corresponding source, even during free flight. A separately prepared native triangle with one excess gave one excess in every joint record. Subtract its mean density, 1/3 per port, from the whole belt, and the two charges become −1/3 and +2/3. Exactly the same subtraction gives exactly the same results for the compact controls. In a second instrument, a single qutrit interferes between the two poles. Its sine quadrature on the three-load wall reverses when we reverse the flight: ±0.26789 for negative probes, ∓0.29139 for positive probes at exposure 1/32. The four-load wall moves while that quadrature stays zero. This detector resolves the occupation of the two interlaced five-port neighborhoods. [Joint polar readouts](results.json)

**Mobile exterior probe.** A third probe starts at port 48 and can move onto port 3 across their actual native seam. After two exchange units the probability to reach port 3 is 70.81% with the belt empty. On the freely launched three-load wall it falls to 0.150% or 0.141% for negative or positive probes, respectively. On the four-load wall it falls to 0.039% or 0.043%. With the bridge guarded the probe stays at 48. This independently measures scattering at the pole. The stationary experiment above tests for a response without probe transit. [Endpoint and simultaneous wall counts](results.json)

## Why these readings differ

For a blank native seam with material records a and b, its diagonal energy is κ(a−b)² plus 1/2 when the two records differ. The exchange matrix element between distinct swapped records is −1/2. After subtraction of the uniform probe self-energy, a signed probe has one extra energy g for each neighboring excess unit. Direct projection of the full 99-state seam gives g = +15.5 for a record-0 probe and g = −17 for a record-2 probe at κ = 8. These are calibrated material contact energies. Changing to exchange alone gives −0.5 and −1, respectively.

On the two poles together the extra generator is exactly g times total excess N. It commutes with every belt exchange. A common Ramsey reading therefore has visibility exp(−i g N t) for every state in a fixed-N sector, including every nonwinding word. The differential apparatus measures the two complementary local occupations instead. With an open gateway the stationary instrument has conditional generator H + g P, where P projects onto actual excess material on gateway 3. The measured overlap between evolution under H and H + g P includes coherent intermediate transport and probe backaction. Neither phase division by the unit standard nor subtraction of the chosen crystalline background yields a probe-independent fractional strength. The companion [flux experiment](../row_12_flux/PROTOCOL.md) probes transport by an imposed U(1) flux. These measurements instead test how the source influences independent material probes.

## Reproduction and scope

From the checkout root, run the following commands (Python with NumPy, SciPy, NetworkX and pytest).

```sh
python -m model.candidate_laws.interval_echo.openwave.row_12_charge.experiment
python -m model.candidate_laws.interval_echo.openwave.row_12_charge.exterior
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_12_charge
```

The default commands recompute all data and compare the archives. Use `--write` on either module to regenerate its own archive. The matrices act on all configurations for the enrolled number of loads, without a wall-band truncation. The 28 belt preparations and 30 stationary exterior preparations per reference service include opposite launch phases, held one-port translations, empty, unit, compact, prepositioned-unit and exchange-only controls as applicable. Every interferometer setting records 262,144 seeded synthetic shots. Common polar counts have shape [quadrature][bright/dark]; differential and stationary counts additionally resolve the terminal source bin. Bin 10 contains all other belt strings. Exterior bin 11 records an occupied gateway. Mobile counts use [near/far][source bin]; settings run at times 1, 2, with a held and open bridge at each time. [Frozen procedures](PROTOCOL.md) · [Instrument implementation](experiment.py) · [Exterior chamber](exterior.py) · [Independent checks](test_experiment.py)

The geometric belt, crystalline preparation, coherent launches and storage arms, ideal regime of repeated sector measurements, inter-cell shutters and readout timings are apparatus. Held-belt motion at fixed load number is identical up to an energy offset for either idle-slot reference. Once source matter enters the gateway, the two reference services differ and their outcomes are reported separately. The experiment includes one native cell, one adjacent probe and specified finite exposures. It establishes neither an electric field traveling through an independently empty region nor a quark flavour or mass. The independently measured one-third slip and material polarization remain properties of these particular observer records. [OpenWave's OW‑12 test](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L114)
