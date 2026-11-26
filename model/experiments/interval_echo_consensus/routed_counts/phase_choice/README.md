# Choosing a phase reading from ordered repair counts

In two driven histories, a sender counted the strict repairs that a cell network recorded. At owner zero a separate apparatus read the sender's decision from ten material rails and used it to choose a phase pulse. Re-reading the same reports backwards changed the resulting local address frequencies. Dropping transitions across replenishment also reversed the choice in both histories. The implemented pairing supplies a phase convention for the declared apparatus. [Frozen counting rule](PROTOCOL.md) · [Live handoff protocol](TRANSFER_PROTOCOL.md)

## What the observers count

Each history contains 96 preparations of nine vacancies and nine doubles in 18 ports on the 20-owner board. Sixty random pair offers follow each preparation. A blank off-tree collar writes for the two opposite extreme loads; every quiet offer counts too. The two runs wrote 680 and 682 times in 5,760 offers apiece. With routing exchange disabled both produced only quiet reports. The raw logs retain the input inventory, all pair indices in lexicographic order of indices into the enrolled port list, and both outcomes for every offer. A window contains 20 offers. Its bucket is the number of writes, capped at seven. [Source, counts and hardware checks](measure.py) · [Complete reports](results.json)

For the ordered window counts C, choose the lexicographically first triple maximizing the magnitude of C(a,b)C(b,c)C(c,a) − C(a,c)C(c,b)C(b,a). A zero magnitude refuses. Every cyclic rotation preserves the product difference; reversal negates it. Thus it suffices to test the 56 strictly increasing triples, including ties between different triples. The tests independently enumerate all 512 ordered triples. Transposing a table preserves every magnitude and hence its first maximizer. Reversal changes the selected sign on every strict table. The symmetric pooled table has zero gap. [General Lean proof](../../../../../Lean/Thermodynamics/RepairCurrentOrientation.lean) · [Independent recount](test_choice.py)

| Ordered history | Selected triple | Signed product difference | With replenishment transitions omitted |
| --- | --- | ---: | ---: |
| First | (0, 2, 5) | −8,728 | +388 on (0, 1, 2) |
| Second | (0, 1, 5) | −2,902 | +544 on (0, 1, 2) |
| Historical B12 count table | (3, 4, 5) | +1,239,691,068 | Different source and bucket definition |

Backward reading transposes the observed counts. This is an operation on the public history. Replenishing the driven source is a physical input. The historical selector was designed after the B12 run; its reproduced sign serves as an archival cross-check. The other two histories and their window rule were reserved before acquisition. [Historical payload](../../../../../code/thermodynamics/repair_current_orientation/repair_current_payload.v3.json)

## What the receiver measures

A controller serially compares the public repair reports, sets the five-dimer code at owner 13 and dispatches the registered transfer. Owner zero interrogates the delivered code against zero and double load standards. An encoded null and a malformed transmission both stop the phase read. The receipt distinguishes the two through the actual twenty-bit readback. With an intact link two native strict repairs and endpoint exchanges prepare the phase stock from the delivered loads. Their pulse program uses only the decoded word and the known source stock. Its two spent collars stay outside the later contact schedule. The phase trial continues on those 240 ports. Compiling the earlier repair history and initializing the transmitting register are supplied controller services. [Live circuit and stock conversion](transfer.py) · [Complete transmitted, received and terminal records](transfer.json)

Under the declared pairing a positive gap selects 3π/2 of addressed contact and a negative gap selects π/2. On the one-excess code the primitive contact angle θ, starting from port zero, prepares coordinates proportional to (cos(θ/2), i sin(θ/2)). Adding the selected dwell before the common pointer read gives address-zero probability (1+s sin θ)/2, where s is the received sign. Writing Uₛ for the analyzer, the selected effect Uₛ†|0⟩⟨0|Uₛ is exactly the projector onto (|0⟩+si|1⟩)/√2 in the calibrated code. This identifies the physical yes effect with the corresponding conjugate projector of the finite theorem. Changing both sign conventions consistently exchanges the names. At θ=0 and π the two analyzers have identical probabilities. The other 16-outcome pointer codes and malformed label are retained.

For the live challenge at θ=π/6 the positive selector gave 6,092–6,187 address-zero reads in 8,192 simulated launches per bank. The negative selector gave 2,021–2,111. At θ=π/3 the ranges were 7,605–7,676 and 546–567. Negative angles exchanged the bright and dark forecasts. All 32 live banks were inside the predeclared simultaneous 99% allowance of 0.02366. For the remaining 64 settings the receiver issued a no-read on every launch, including every bridge-cut control. The 16 separate primary calibration/held banks use a separate joint 99% allowance. Statistical coverage assumes independent simulation draws under the calibrated instrument. Entropy production, an absolute imaginary direction and additivity on every quantum effect are separate questions.

From the RER root, export OPENBLAS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1 and OMP_NUM_THREADS=1, then replay both modules.

~~~sh
base=model.experiments.interval_echo_consensus.routed_counts.phase_choice
nice -n 15 .venv/bin/python -B -O -W error -m "$base.measure"
nice -n 15 .venv/bin/python -B -O -W error -m "$base.transfer"
~~~
