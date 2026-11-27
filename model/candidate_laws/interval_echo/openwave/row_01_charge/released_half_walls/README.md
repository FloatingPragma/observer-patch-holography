# One ordinary load releases two mobile polarization walls

Fill one vacancy of a sharp alternating ring and let every ring contact evolve. Two stations can then detect separated changes of alternating order, each carrying **half the response per distance of an ordinary load** to a calibrated potential gradient. The two initial defects were adjacent. No pin pulls them apart. A second experiment on the emitted source distinguishes the half-unit polarization from the integer number of loads in a sharply bounded region. [Frozen design](PROTOCOL.md) · [Receiver records](records.json) · [Complete terminal word banks](words.npz)

| Source at exposure 3 | All-start probability of both remote walls | Same source with flight contrast off |
|---|---:|---:|
| 14 ports, inject or remove one | 0.426487 | 0.015890 |
| 16 ports, inject or remove one | 0.423120 | 0.009202 |
| 14 ports, inject or remove a double | 0.433357 | 0.015890 |
| 16 ports, inject or remove a double | 0.431376 | 0.009202 |

The event requires exactly two same-sign phase slips on the entire ring, one in each of two disjoint three-bond gates. Both original wall bonds must have returned to their alternating readings. Independent native trit meters at two stations read every occupation. On the 20-cell board, 1,789,126 of 4,194,304 ordinary positive launches passed the two gates at time three, compared with 66,814 on the changed law. The corresponding positive counts were 1,788,341 and 1,788,169 on the two larger deposits. All four registered geometries, both signs, seven holds and all failed-gate words are retained. The three fourteen-port circuits have different physical addresses and the same reduced generator. The two-wall probability at time three, before applying any location gate, is 0.98333 on fourteen ports and 0.97981 on sixteen. With the contrast off it is 0.09158 and 0.04997. [Actual ring contacts, generators and evolution](ring.py) · [Native heads and individual rail codes](readers.py)

To identify the charge carried per displacement, a second preparation places both walls apart. Two neighbouring words differ by one native load exchange. A moving wall advances two bond spacings while the other stays fixed. A separate single-load preparation on those contacts advances just one site spacing. A linear site potential couples to the actual ternary load through two independently calibrated contact references. The observed relative phase divided by the measured displacement gives the following ranges over two held placements and four geometries. The calibration used sharp loads zero, one and two alone. [Calibration and gradient interferometers](readers.py)

| Moving interface | Gradient response per bond spacing in ordinary load units |
|---|---:|
| Equal occupied neighbours on the 0/1 host | +0.49980 to +0.50039 |
| Equal vacant neighbours on the 0/1 host | −0.50037 to −0.49975 |
| Equal occupied neighbours on the 0/2 host | +0.99969 to +1.00086 |
| Equal vacant neighbours on the 0/2 host | −1.00049 to −0.99957 |

A simultaneous 99% Hoeffding bound over the 1,296 calibrated and held fringe quadratures gives a conservative counting-only radius of 0.0067 for the wall-response ratios. Each bright fraction has radius 0.00122 before conversion to phase and division by the common ordinary calibration. The ordinary comparison reads one or two units per site spacing as prepared. Reversing the applied potential reverses the fringe phase. With the potential parked, every inferred wall response lies within 0.00064 of zero. The gradient reads phase work at a controlled impulse. The separate initial release measures propagation under the native ring generator.

## Which electric interpretation survives?

One voltage probe weights every material site in its window equally. Another sums adjacent-bond polarizations by giving each end site half weight and each interior site full weight. At a full turn on the accepted single-load pairs, the equally weighted probes both return bright. The tapered probes both return dark. In separate level-zero positive banks, 1,788,830 and 1,790,658 launches respectively produced those joint accepted outcomes, with zero accepted outcomes in the other three ports. Every bank contains 4,194,304 launches. On the double-stock source the full turn returns bright for both profiles. The samples are seeded multinomial word and joint-probe draws, or binomial quadrature draws, from the computed native probabilities. The angular and count resolutions describe the simulated sampler. [Both probes, acceptance and all other outcomes](records.json)

The two-valued order slip is visible directly in the alternating pattern. The sign and magnitude of its transported polarization depend also on the local occupation and the background material. Both walls from adding one ordinary load carry positive half-unit polarization even though their changes of bulk phase have opposite orientation. The static high-energy medium, parked exterior exchange, probe profiles, clock and impulse electrodes are supplied. There is no independent dynamical electromagnetic field on this source.

**OW-01 is Partial on this evidence.** The result establishes mobile, quantized material polarization and a concrete division of the ordinary port unit in that transport measurement. A sharp region still returns an integer port count. An integer topological electric charge and the exclusion of physical half-charge daughters require an additional identification. The voltage response alone establishes no general electric-force premise. [Mathematical scope and uniform retention bound](PROOF.md) · [OpenWave criterion and M10 decision](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#simplest-test-per-criterion)

Two independently prepared mobile probes resolve the whole local contact detuning of an accepted wall. The same detector counts an ordinary load increment on each of two consecutive bonds. Fixing its total on known stocks assigns half a unit to each isolated wall for every additive local detector in the tested class. Leaving the ring mobile during readout also changes the probability of observing the two walls. [Calibrated material probes](material_probe/README.md) · [Finite calibration theorem](charge_identification/PROOF.md) · [Independent reconstruction and alternatives](probe_audit/README.md)

## Replay

From the repository root, use one numerical thread.

~~~sh
pkg=model.candidate_laws.interval_echo.openwave.row_01_charge.released_half_walls
folder=model/candidate_laws/interval_echo/openwave/row_01_charge/released_half_walls
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 15 .venv/bin/python -B -m "$pkg.measure"
nice -n 15 .venv/bin/python -B -m pytest -q "$folder"
(cd Lean && nice -n 15 lake env lean -j1 -DwarningAsError=true "../$folder/WallCharge.lean")
~~~
