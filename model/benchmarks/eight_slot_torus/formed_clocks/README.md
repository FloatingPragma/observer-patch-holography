# What time can a formed pair tell?

An ordinary head-on collision can start a stationary material clock on the torus
law. Two incoming trains take the native update and meet at a capture port.
The capture rotates their incident amplitude into a matched pair and leaves a
marker on the local wire. A supplied instrument then directs the pair's own
hops. After loading, the free contact update advances the relative phase of
its two conserved rails. A station reads this phase from four local outcomes.
Its elapsed-time reader receives counts from two appointments, with the
birthplace, birth time and appointment labels removed.

### From an incoming pair to an observer record

| Question | Result on the declared apparatus |
| --- | --- |
| How quickly can a known birthplace load a persistent clock? | **26 native services**, plus one for the ordinary approach. Every successful natal write transfers, including an arbitrary rail state entangled with a reference. A shell-transport bound rules out every deterministic on-site controller through service 25. The first write itself succeeds with probability one half for a sharp head-on input. |
| Must the apparatus know where or when a pair formed? | One identical table at every site handles arbitrary natal superpositions in 35 services. The pair's local age controls the table. It uses both routing parities to transmit all amplitudes isometrically. The best success probability a single coherent port can guarantee at every birthplace for the designated positive-parity plaquette is one half. |
| Does a working source keep the mode populated? | On the 81-site carrier, all 81 wires received writes from one ordinary launch ensemble. By tick 193, 46.27% of original launches had a mark, 45.57% had reached the completed program, and one root detector received 11.97% in its four late ports. All write powers and unmarked outcomes enter the replay. |
| What do the counts measure? | In mature material both routing parities satisfy $U_r^2=e^{2i\mu_r}I$, where $\mu_r=\varphi_r/2$. Two services reproduce the full local spatial state and advance the signed hand by $-\pi/2048$. Eight hidden even lapses from 2 to 7,926 services each had exactly one solution in the booked range. Six of eight single-station panels also recovered every lapse on their own. The two weakest stations retained 20 and 80–81 candidates. |
| Does the result survive motion? | Equipped rail-blind controllers carried a loop through six mature native hops and back. The six count estimates of its advance relative to the parked hand ranged from 0.00114686 to 0.00115238 radians. All six intervals contain the forecast of 0.00115049 radians. The apparatus uses a direction-selective filter, local controls and stationary destination counters. |
| What can a moving detector read with the plates off? | On a spreading ordinary birth the moving aperture resolves both signed phase advances. The low branch meets its complex-amplitude forecast at all seven visits. The high branch misses it at every visit, with errors up to 32.93%. |
| What happens when the native age returns to zero? | The loaded two-parity material state revives after four cycles of 8,192 services. A returned, occupied marker cannot release this dark state. A repeating loading table produces 12.665% capture-active weight on its first restart. A supplied five-state program digit preserves four passive cycles before it reenables the plates. |
| Can one receiver tell which cycle passed? | At a calibrated natal corner the total admission separates four laps, and the rail counts separate 4,032 mature even ages on each lap. All eight blinded acquisitions selected the correct epoch modulo 32,768 services. The saturated opening leaves four lap aliases. Rail dephasing leaves 4,032 age aliases. [Acquisition and finite-shot bound](clock_epoch.py) |

The known-place loader consumes coherent local operations on the matched
slots. The origin-independent device also reads the bound age as a local
control. The age-blind terminal uses just one enrollment bit for an early
window, four rail ports, and explicit misses and refusals. The time comparison
uses an independent free optical channel at the observing site. Its private
pulse returns every two native services. A reversible supplied counter counts
those returns. The elapsed decoder uses the source rail splitting, the
matching parity of the endpoints and an 8,191-service upper bound. Learning
the rate from the counted returns also admits the alias displaced by
$\pi$ radians per service. The two rates give identical even-appointment
counts. [Circuit, bounds and transport](CLOCK_LOADING.md) ·
[Blind comparison and controls](clock_comparison.py)

The enrollment condition has a measurable consequence. With births
continuing and a terminal that reads every arrival regardless of its age,
all 81 stations deviate from the mature two-service relation on both
registered intervals. Some apparent increments have the opposite sign.
The mature part obeys the relation to within $3.5\times10^{-17}$ in the
complex quadrature per original start. The discrepancy is accounted for by
young pairs and intervening births. [Complete unfiltered test](clock_unfiltered.py)

### Separating the apparatus from the material

With the local plates off, an ordinary capture has zero overlap with the
protected loop. With the plates on, changing one split leaves about 6.84% of
the natal weight outside the target. A single rail-blind loading calibration
on the simpler 28-service circuit gives worst-input fidelity 0.999954. A local
two-state memory carrying one excitation can prepare and read the clock
through charge-preserving exchanges. Tracing over the memory removes the
hand while keeping the same local click power.
[Source operations and controls](clock_loader.py) · [Original counts](loaded_clock_counts.npz)

Ordinary births also give two distinct clocks without the loading program.
A fixed two-site aperture on a long axial strip reads the dispersing bright
wave. Its signed 257-service advance was 0.197177 radians against a forecast
of 0.197117, at an admission of only $8.16\times10^{-6}$ per original start
at the last appointment. A single-rail point detector on rings of five and
six sites reads 0.390631 and 0.385518 radians per service. All 24 early fits
predicted the held probability records within $3.0\times10^{-5}$ per start.
Three broad five-site settings missed the registered frequency-precision
threshold. These lines depend on ring size. [Dispersive aperture](clock_modes.py) ·
[Single-rail receiver](formed_clock_readout.py)

The result concerns finite observer records on prepared carriers. The
protected modes have flat spatial eigenphases; assigning their rail beat to
a physical particle mass would also require a mobile material dispersion
and a source-matched physical standard. No measured masses, couplings or
clock standards set the controls. Multi-particle interactions, freely moving
proper time and comparison with physical particle data have separate source
requirements. OpenWave's physical row decisions belong to its validators.

### Replay

The first command runs short literal checks of both the lower bound and the
full two-constituent circuit. The other commands regenerate the frozen
sources, verify the original acquisitions and print the decisions, including
unfavorable controls. Run at the repository root with
[the model dependencies](../../model/requirements.txt).

```sh
python -B -m pytest -q model/benchmarks/eight_slot_torus/formed_clocks/test_clock_loader_optimal_literal.py
python -B -m simulator.formed_clocks.clock_comparison
python -B -m simulator.formed_clocks.clock_budget
python -B -m simulator.formed_clocks.clock_ballistic
python -B -m simulator.formed_clocks.clock_cycle
python -B -m simulator.formed_clocks.clock_epoch
python -B -m simulator.formed_clocks.clock_stream
python -B -m simulator.formed_clocks.loaded_clock
python -B -m simulator.formed_clocks.clock_transport --check
python -B -m pytest -q -W error model/benchmarks/eight_slot_torus/formed_clocks
```

The comparison booked 216 words of $2^{36}$ original launches each.
Acquisitions are simulated multinomial counts from a fixed finite source.
The complete original-start likelihoods, count words and local optical
timestamps are available offline. [Comparison acquisition](clock_comparison.npz) ·
[Stream acquisition](clock_stream.npz)

A separate acquisition used a single prebooked receiver in each of two
fixtures and one sixteenth as many launches. Each receiver independently
decoded all eight lapses. With free evolution alone, neither receiver
gave the correct unique value at any of the eight visits. With rail
coherence removed, all 4,096 even candidates remained possible at every
visit. The active, passive and dephased controls together cost
927,712,935,936 simulated launches. [Resource calculation and acquisition](clock_budget.py)

Bernhard Müller's selected-law runs separately establish a heralded
eight-arm preparation, a birth-accessible bright sector and a moving
two-pole clock on an ordinary capture-dark input. Those preparations and
their clocks have distinct source maps. The torus replays use the native
reflected routing, contact phases, capture and opening specified in the
[simulator](../SPEC.md). The records of those runs are the OpenWave
internal-clock and protected-pair reports.
