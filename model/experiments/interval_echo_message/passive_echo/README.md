# What can an observer learn about a message from a passive echo?

On the native 20-cell Interval Echo screen, a receiver can predict its probability of reading a distant bit from separate round trips on the same material signal. On the held-out sender in cell 15, the receiver and sender jointly certified **57 and zero** returns at six positive clock ticks for its two encoding ports, in 1,073,741,824 simulated launches per port. Those counts predict **50.00371–50.01474%** correct answers in a three-tick one-way experiment. A separate bit bank gave **50.01152%**. The marked echoes thus resolve an advantage on this source. Across all three sender stations and all four deadlines, the native single-launch receiver falls below the specified 52% threshold. [Protocol and disclosed pilot](PROTOCOL.md) · [Complete simulated records](records.json)

## Preparing, ranging and decoding

The receiver keeps a double at port 2 and uses a freely propagating single as its signal. In one-way trials the sender chooses a bit after setup with a native exchange inside its cell; the destination is one of a fixed pair of source ports. In range trials the single starts on the receiver's gateway, port 1. All 240 ports and all native contacts remain active during flight. The two loads interact and exchange on one board. The supplied apparatus checks the single local question “Does port 2 hold two units?” after each unit of free exposure. A yes advances the clock, a no expires the trial. Four registered message deadlines require one through four successive positive readings. Reading this one port is a more selective clock than the two-port clock of the routed-distance study. Exposure per service is supplied, and a recorded positive tick and a unit of physical time have different calibration requirements. [Signal and native Hamiltonian checks](../contact_deadline/README.md) · [Comparison with routed distance](../ranged_deadlines/README.md)

At the one-way deadline, the receiver reads whether the single is at its gateway. It guesses zero on a gateway click and one on every other outcome, including expiration. The stored terminal bank also distinguishes the other eleven ports, quiet, and which interrogation expired. To range a selected bit port, a remote observer reads whether the single is there after the same number of appointments. This local mark is the only intervention on the signal. Free flight continues without a return pulse. The logs of a marked remote visit, a subsequent home click and a live clock certify the round trip when the observers compare records. The one-way decoder receives no remote log. The return calculation retains every clock expiration and every unmarked outcome. On the far three-tick pair, 171,733 and 177,982 home clicks followed a remote miss. Counting these as round trips would obscure the marked counts of 57 and zero. The release ports for the first two stations lie in the same sender cell. [Acquisition](measure.py)

## Why a two-way reading forecasts a one-way probability

In the joint port basis the native free Hamiltonian is real and symmetric. The local tick projector P is real and diagonal. One time step conditioned on a tick acts by V = P exp(−iH) P on a release with its double at port 2. V is complex symmetric. Its nth power has the same transpose symmetry, even though probability leaks out whenever the clock fails. Write a for the complete word with the double at 2 and the single at 1, and b for the word with the same double and the single at the marked sender port. Then

$$
p_n(a\to b)=|(V^n)_{ba}|^2=|(V^n)_{ab}|^2=p_n(b\to a).
$$

On a clock-surviving branch a positive remote mark leaves exactly the word b. The two-way marked return probability is therefore the square of the one-way gateway probability. If the two bit ports have return rates r₀ and r₁, the predicted balanced accuracy is (1 + √r₀ − √r₁)/2. This equality uses the same stationary generator and the same successful clock projector on both legs. The simulator propagates the marked and unmarked branches separately; the forecast reads only independent round-trip counts. For zero simulated returns the uncertainty calculation supplies a positive upper limit. [Count-only forecast](check.py)

| Bit-zero source | Tick | Simulated echo counts for bit zero / one | Accuracy inferred from echoes (99% interval) | Accuracy on separate one-way simulated launches |
| ---: | ---: | ---: | ---: | ---: |
| 82 | 1 | 684,958 / 3,493 | 51.16658–51.17873% | 51.17296% |
| 72 | 2 | 1,350 / 58 | 50.03821–50.05044% | 50.04339% |
| 182 | 3 | 57 / 0 | 50.00371–50.01474% | 50.01152% |

Every row uses 2³⁰ launches per source and service, including all expirations. The direct and inferred intervals overlap at every enrolled setting for both laws. The simultaneous 99% intervals use exact binomial endpoints with a union bound over the three services, three stations, two release labels, four deadlines and two flight strengths. Halving the exchange on all seams at fixed contrast makes the near gateway reader reach 52.47672% at tick two and 52.96970% at tick four. All no-encoding control intervals contain chance. This comparison shows how the occupied clock hand, signal transport and interrogation interact; weakening a hopping rate need not reduce the chance of a usable clocked decision.

The primary forecasting rule was fixed before the counts. A separate analysis asks whether another classical reading of the complete single-launch local transcript could meet 52% under the native law. There are at most seventeen different local records. For each possible deterministic decision table, the balanced empirical accuracy is the mean of 2 × 2³⁰ independent correctness indicators conditional on the input labels. Hoeffding's inequality and a union bound over every table and all 24 settings give a joint two-sided allowance of 0.00687 percentage points. Summing the larger of the two bit counts in each bin gives the empirical optimum. Its largest native value was 51.8488%, at port 82 and one tick. The simultaneous upper limit is below **51.856%** even allowing a separate decoder at every station and deadline. A randomized classical decoder is a mixture of deterministic ones.

These records test reciprocal transfer and complete-message performance in one finite clocked service. They provide no calibrated clock distance, unique propagation speed, boost law, or bound for other detectors or coherent measurements. OW-21 and OBS-03 retain their partial local assessments. [Observer evidence](../../../evaluation/EVIDENCE.md#obs-03)

From the RER root set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and VECLIB_MAXIMUM_THREADS to 1 and run

```sh
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.passive_echo.measure
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.passive_echo.check
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.passive_echo.verify_flight
```

The final command reads the raw board incidence and independently evolves a 240 × 240 ordered-pair array with a degree-80 Chebyshev expansion. It compares 624 held-out terminal bins to simultaneous binomial prediction intervals on both laws, including unmarked return addresses. This independently checks the forward propagator and recorded counts. Both computations use binary64 arithmetic, the same material assumptions and the same local projector. [Independent propagator](verify_flight.py)
