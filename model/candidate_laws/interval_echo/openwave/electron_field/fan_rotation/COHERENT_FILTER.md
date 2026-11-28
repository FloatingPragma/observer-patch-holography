# Finite coherent threshold measurement on the 320-cell guide

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

The observers can turn the two-level energy admission condition on the asymmetric three- and four-fan routes into a finite controlled-clock measurement. Writing just a binary decision and then checking that the clocks reset retains the two field amplitudes. The proof bounds the change in the joint arrival record and in the separated-field records after this finite admission. It holds for both charge signs. [Executable](coherent_filter.py) · [Receipt](coherent_filter.json) · [Path and detector](../fan_path/PROTOCOL.md)

## Fixed clock and decision

Use the native rate-4 intercell terms and the physical 320-cell routes from the path test. Use the same source and the same receiving hub and two outer-ring modes. At the receiver, the complete comparison also tests that every other fan is in its vacuum. Independently park the material bridges to calibrate the separated charge and field reference.

The exact count at 4.25 is two, and the exact count at 4.53125 is also two. The parked threshold lies between 4.5 and 4.53125. Thus the projection below the parked threshold is the same as the projection below 4.5. The generator is positive semidefinite and its integer row-sum bound keeps all phases away from the wrapped boundary at the chosen clock exposure.

Use 41 coherent phase-estimation registers, each with 18 bits and an exposure step of $1/64$. The unitary in each clock is the native flight $e^{-iH/64}$ under a supplied control. Convert each signed Fourier bin to an energy with the opposite sign, and vote yes if it represents an energy between −1 and 4.5. The reversible comparator uses integers −651 through 2933 inclusive in this convention. Write a single yes bit when at least 21 registers vote yes. Reverse the votes, counter and all phase-estimation circuits, and measure just the yes bit and the binary clock-reset outcome. No fine energy is recorded. Proceed to the free wait and physical receiver only when the decision is yes and the clocks reset.

On an energy $E$, write $f(E)$ for the probability that the majority would vote yes if read before reversal. The yes-and-reset amplitude is $f(E)$ because it is the overlap of the clock with its projection onto the yes outcomes, evolved back to the initial clock. The no-and-reset amplitude is $1-f(E)$. These two records have probabilities $f(E)^2$ and $(1-f(E))^2$. The third record combines both decision bits when clock reset fails and has probability $2f(E)(1-f(E))$. The three effects sum to one. This preserves coherence between retained energies on the accepted branch. The small replay check explicitly applies the clock unitary and its inverse, including the reset-refusal contribution.

## Uniform error bound

An incorrect decision bin for one register is at least $1/32$ native energy units away in the circular phase coordinate. This is more than 20 of the $N=2^{18}$ Fourier intervals. Let $d$ be the shortest phase distance as a fraction of a turn. The geometric sum for a clock bin has probability at most $1/(4N^2d^2)$ since $|\sin(\pi d)|\ge 2d$ for $0\le d\le 1/2$. There are at most two grid bins with $Nd$ in each unit interval. The full incorrect-label probability is therefore at most $\frac12\sum_{k\ge20} k^{-2}\le1/38$. The votes are independent conditional on an energy eigencomponent. The probability of at least 21 bad votes is at most $2^{41}/38^{21} < 1.469\times10^{-21}$.

Thus the accepted operation $f(H)$ differs from the exact projector $P$ in operator norm by at most $\eta=2^{41}/38^{21}$. The proof includes the entire spectrum. It carries onto arbitrary superpositions and states with an untouched reference. The independent inertia, dyadic residual and source/receiver overlap checks replay on the native fixture.

## Arrival and separated-field records

Write $a$ for the verified minimum singular value of admission from the two-input source to $P$, and $b$ for the minimum singular value from $P$ to the receiver with the other fans in vacuum. For every wait, the intervening native unitary on $P$ preserves lengths. Ideal arrival therefore has amplitude at least $ab$ for every input. With the finite clock, the all-launch probability of admission followed by the complete joint arrival is at least $(ab-\eta)^2$. It is strictly positive and the successful two-input map has full rank at every wait. The probability of the admission alone is at least $(a-\eta)^2$.

The conditioned vector after actual admission differs from the ideal conditioned vector in norm by at most $\Delta=2\eta/a$. This follows by normalizing both vectors and bounding the change in their norms. The conditioned joint-receiver probability is at least $(b-\Delta)^2$. For a detector of field excitation at least $d$ links from the charge, let $Q_d$ denote its ideal probability bound. The same finite-admitted ensemble gives separated-field probability at most $(\sqrt{Q_d}+\Delta)^2$, capped at one. These bounds concern records from genuinely accepted runs and do not normalize by an assumed perfect filter.

| Observers' bound for either sign | Three fans | Four fans |
| --- | ---: | ---: |
| Minimum yes-and-reset probability | $2.27710\times10^{-9}$ | $8.94813\times10^{-9}$ |
| Minimum complete arrival probability per launch | $3.03362\times10^{-16}$ | $1.08546\times10^{-19}$ |
| Minimum complete arrival probability given admission | $1.33222\times10^{-7}$ | $1.21306\times10^{-11}$ |
| Maximum separated-field probability given admission, by distance | 5.42471%, 1.44662% | 2.66073%, 0.421923%, 0.074051% |
| Minimum magnitude of the receiver-only negative partial-transpose eigenvalue, per launch | $3.71985\times10^{-14}$ | $1.68766\times10^{-16}$ |

For the receiver-only read, keep a normalized maximally entangled input and trace the other fans after the local detector. The two unnormalized outputs before tracing are vectors of norm at most one and differ by at most $\eta$. Their outer products differ in trace norm by at most $2\eta$. The partial trace contracts this distance, and partial transposition on the receiver preserves Frobenius norm. The minimum eigenvalue of the receiver-only Choi partial transpose therefore moves by at most $2\eta < 2.937\times10^{-21}$. The executable recomputes the ideal moment witness and subtracts this error to obtain the two positive magnitude bounds shown above. Thus the actual finite protocol has a negative partial-transpose eigenvalue at every free wait with the other fans unobserved.

The protocol supplies controlled access to the entire fixture generator and its inverse. It uses 738 phase qubits, 41 vote bits, a six-bit reversible counter and two binary record bits. There are 1,476 controlled power calls in both directions, with maximum single-block exposure 2,048 and sum of absolute exposures $10747863/32 = 335870.71875$ native time units. The receipt counts the Fourier and classical reversible gates and gives an explicit literal comparator construction bound. These costs are logical controls, with multi-controlled gates counted before decomposition and controlled evolutions treated as calibrated primitives. The proof assumes exact primitives. An operator-norm implementation error on the full clock operation would add to $\eta$ in each bound. The finite clock and successful read establish a simulated nonzero two-input transfer in this supplied apparatus. Physical implementation and a dissociation measurement beyond the finite parked reference require their own records.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
    python -B -W error model/candidate_laws/interval_echo/openwave/electron_field/fan_rotation/coherent_filter.py
~~~
