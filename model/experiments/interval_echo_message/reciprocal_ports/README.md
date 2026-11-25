# Predicting a one-way message from passive echoes

An observer can calibrate the original one-load terminal reader using successful returns to the sender. The calibration needs the port at which the load turned around. On the 20-cell board at exposure four, direct one-way detection in the receiving cell gave **42.52%** for bit zero. Address-labelled echoes predicted **42.54%**. Taking the square root of a single coherent whole-cell echo gave **13.66%**. Discarding the port labels after a sharp read gave **17.91%**. These are three different experiments on the same flight. [Enrollment and exposed settings](PROTOCOL.md) · [Counts](records.json)

## What the sender and receiver record

Prepare the original sharp releases at the two bit ports. At the midpoint the receiver either finds a particular occupied port or records quiet. On every arrival the unit flies freely for the same duration again, and the sender asks whether it returned to the exact launch port. To calibrate the reader, use only the twelve counts with both an addressed midpoint arrival and a successful return, each divided by all starts. The two legs use equal supplied intervals of native exposure. The acquisition uses no clock register, and the emitter and receiver can compare their records after each run. [Native two-leg acquisition](measure.py)

For the one-load blank sector on the full board, the Hamiltonian is real symmetric in the port basis. Its unitary \(U\) therefore satisfies \(U^T=U\), including the common contrast phase. For source \(s\), the one-way address probability is \(p_r=|U_{rs}|^2\). A sharp read of \(r\) leaves that same local port state at the turnaround. The two-part record has probability \(e_r=|U_{rs}|^2|U_{sr}|^2=p_r^2\). This equality requires matched passive legs. It can fail on a directed three-port cycle that carries \(s\) to \(r\) with certainty and sends \(r\) to the third port.

The one-way forecast is \(\sqrt{e_r}\) for each port and \(1-\sum_r\sqrt{e_r}\) for quiet. In a sharp measurement that suppresses the address labels, the aggregate return probability is \(\sum_r p_r^2\). A single yes-or-no arrival read preserves the local superposition instead. Its return probability is \(\left|\sum_r U_{rs}^{\,2}\right|^2\), so phases as well as relative port weights enter that echo.

## Comparison with direct reads

Each of the 40 panels has three separately sampled banks of \(2^{22}\) starts. They cover both release ports, the 20- and 320-cell boards, five appointments, and native and half-strength exchange during both legs. An independent incidence-based Taylor evolution generates the fresh one-way challenge distributions. Across 520 challenge bins, simultaneous 99% exact binomial intervals all overlap the address-resolved echo forecasts. The aggregate square-root forecast falls below the lower limit for local arrival in 22 panels for a coherent whole-cell read and 19 panels for a sharp read with labels suppressed. The check retains all panels, including zero exposure. [Count-only forecast and simultaneous limits](check.py)

On 20 cells at exposure four, a decoder trained solely on the echo coincidences recovered 5,829,314 of 8,388,608 balanced challenge bits. At half exchange it recovered 5,214,625 using its own echo calibration. The symmetry predicts reciprocity for both laws. The timing and arrival distributions depend on the exchange strength. The result calibrates terminal probabilities and supplies no distance, speed, or free observer clock.

Run these modules from the repository root to replay the banks and their forecasts.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m model.experiments.interval_echo_message.reciprocal_ports.measure
.venv/bin/python -B -W error -m model.experiments.interval_echo_message.reciprocal_ports.check
~~~
