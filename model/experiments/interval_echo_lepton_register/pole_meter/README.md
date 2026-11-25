# Which formation energy do the three charged preparations share?

All three Fourier outputs of the marked-edge vacancy detector have a component at **the same formation line** on the freely evolving twenty-cell board. A distant compact vacancy provides the phase reference. For 4,194,304 simulated starts per setting, the signed bright-minus-dark frequencies were 0.121319, 0.231638 and 0.313892. A separate calculation from the cell geometry gives pole weights 0.121906, 0.231008 and 0.313753. At the specified analyzer phase these three preparations yield the familiar lepton ratios in detector *intensities*. Released from those same sharp outputs they have broad, overlapping energy distributions. [Complete counts and controls](results.json) · [Native source and reader](../README.md)

| Starting state | Common-line weight from geometry | Bright minus dark, per start | Formation-energy centroid | Spectral standard deviation |
| --- | ---: | ---: | ---: | ---: |
| Character 0 (brightest intensity) | 0.121906 | 0.121319 | 36.164 | 5.912 |
| Character 1 (weakest intensity) | 0.231008 | 0.231638 | 37.723 | 6.189 |
| Character 2 (middle intensity) | 0.313753 | 0.313892 | 37.363 | 6.135 |
| Independently prepared home compact mode | 1 | 1 | 35 | 0 |

The final two columns are spectral moments of the registered generator, in the storage convention of the calibrated reference. A separate short-wait material-clock read gives the transfer-energy slopes 1.136, 2.705 and 2.348 relative to the compact line. Finite-wait biases are bounded by 0.023, 0.030 and 0.028 units respectively before sampling. The local timer observes a dimer after sixteen elementary steps and recovers a step of 0.015635 against the scheduled 1/64. Tripling the exchange coefficient during free evolution moves the compact mode to 41 units and gives transfer slopes −1.856, 2.863 and 1.742. [Pulse and energy tests](test_pole.py)

The same instrument on the neutral receiver's charged magnetic doublet resolves two positive common-line weights. An additional native herald transfers an arbitrary doublet input into the free compact qubit with an input-independent success probability of 0.00131910. The inherited spatial turns leak from this energy sector. [Joint magnetic and formation read](magnetic/README.md)

## An energy read with one vacancy on both arms

Use a balanced addressed seam contact to split the same vacancy between one home rail and one remote rail, with reversible swaps to bring the two ports together. A home analyzer run backwards prepares the chosen character. The remote observer prepares a compact vacancy by rotating the sharp rail into a protected cell vector on a matching chart. Every federation seam then acts freely for the same exposure on the two arms. Reverse both preparations, compare the two sharp rails through the same routed contact in X or Y, and retain a third outcome for finding the vacancy anywhere else. The test independently propagates the complete 240-port state and replays the full-word echo restriction and routed contacts. A separate phase setting counts lost coherence and a shifted spectral window tests the frequency discrimination.

The two arms have identical total charge at every stage. The common line is 35 exchange units with the declared filled host and degree-six storage reference. Refitting the stored sample-to-reservoir counts against their independent material timer gives 34.99430 ± 0.00474 exchange units for that transfer (one standard error). Adding a common chemical-potential term to the entire vacancy sector shifts every line together and cancels in these clicks. The model does not determine a conversion to keV from this read. [Storage experiment](../../../candidate_laws/interval_echo/openwave/row_04_compact/README.md) · [Energy reference theorem](../../../candidate_laws/interval_echo/openwave/row_03_clock/neutral_reference/README.md)

Let $s$ denote the released character, $r$ the distant eigenstate with energy $E_r$, and $a(t)=e^{iE_rt}\langle s|e^{-iHt}|s\rangle$. Free evolution preserves $r$ and $\langle r,s\rangle=0$. After a relative phase $\phi$, the bright and dark probabilities are $(1+|a|^2\pm2\operatorname{Re}(e^{i\phi}a))/4$. The third probability is $(1-|a|^2)/2$. Their difference gives the real or imaginary return amplitude against the same material reference. Detuning the window by nine exchange units or dephasing the arm before flight gives a signed mean within 0.00069 of zero for every input. The compact positive control yields a bright exit at every wait on the resonant setting. [Protocol](PROTOCOL.md) · [Executable](measure.py)

For the common pole choose $t=|i-j|/4$ with $i,j$ independent uniform integers from 0 through 511. Averaging the real quadrature gives a positive spectral window with value one at the reference. At offset $\delta$ its response is $|\sum_{j=0}^{511}e^{ij\delta/4}/512|^2$. Each trial samples the complete delay-and-outcome law with the selected delay unread in the final tally. The simulated banks retain bright, dark and every other outcome as a single complementary bin. For the twelve randomized-wait means together the 99% Hoeffding half-width is 0.001927. Finite-window contamination is between zero and 0.006105 in each resonant mean. A common-line weight is therefore bounded below by its measured mean minus those two errors, and above by its measured mean plus the sampling half-width. No fine-frequency fit enters this interval.

## Why the same line occurs at every analyzer phase

Let $P$ be the block-diagonal sum over cells of the exact confining projectors. The charge Hamiltonian $H=L/2+8\operatorname{diag}(2d-6)$ and these projections obey $(H-35)P=0$. The program verifies the identity and $P^2=P$ in integers. It certifies $(H-35)^2+P-I/25>0$ on the full board using a dyadic Cholesky factor and an exact, strictly diagonally dominant positive remainder. On the orthogonal complement, every energy therefore differs from 35 by more than $1/5$. This identifies the entire 40-dimensional spectral projector. The Gershgorin interval $[32,54]$ places every other response outside the main peak and its aliased copies. The geometric-series bound gives the quoted leakage.

At the two endpoints and two incident tips of the marked edge, the tip block of $P$ is $\bigl(\begin{smallmatrix}2&-1\\-1&2\end{smallmatrix}\bigr)/6$ and the endpoint rows vanish. The three Fourier modes therefore have exact weights
$$
w_k=\frac{2-\cos(2\theta+4\pi k/3)}9.
$$
Their total is $2/3$ and each lies between $1/9$ and $1/3$. Even an ideal, perfectly resolved energy measurement from a *single* preparation can correctly guess the character with equal-prior probability at most $1-(\sum_k w_k-\max_k w_k)/3$. This is at most $8/9$ at every phase and 0.882363 at the nominal setting. A numerical read grouping eigenvalues within $10^{-8}$ gives 0.574483 on the nominal sources. This finer diagnostic depends on the grouping resolution. Repeated preparations can distinguish different weight distributions on shared energies. [Lean inequalities for the pole and discrimination](../../../../Lean/SourceLaw/ChargedRegisterPole.lean)

Run from the repository root with single-threaded BLAS.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_lepton_register.pole_meter.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_lepton_register/pole_meter
```
