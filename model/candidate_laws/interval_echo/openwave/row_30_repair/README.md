# A home detector tests the zero-momentum matter wave after local repair

An observer can reject a sharp rest frequency without choosing a zero of energy. Prepare the uniform one-load wave on a connected region, let it evolve freely, undo its preparation, and ask whether the load comes back to the launch port. A single temporal frequency gives a certain return, regardless of its value. On 27 settled sharp-repair preparations the home detector missed **3.63–5.60%** of launches at an appointment of 1/16 exchange unit and **36.48–59.11%** at 1/4. Exchange alone, on identical written material with the same source, clock and detector, gave zero misses. Each joint clock-and-detector setting represents 4,194,304 simulated launches. [Protocol](PROTOCOL.md) · [Counts, intervals and control settings](results.json)

## Where the spread in frequency comes from

Start each port at 0 or 2. The local strict repair turns an adjacent 02 or 20 into 11 and writes its signed contrast on the collar. Both endpoints are then ineligible for a further strict repair. Hence the number $r_v$ of signed writes incident to a given port is either 0 or 1. We use the existing local priority service, with a 1/4 service probability. After settlement, apparatus parks the residual loads and brings in one probe. It preserves the deposited sign and coin registers. This preparation uses neither a wavelength nor an ordered channel. Twenty-four histories were previously recorded for the radar experiment. A further seed was reserved for this test on each of the three board sizes. [Original repair histories](../row_28_acceleration/README.md) · [Earlier empty-collar preparation](../row_30_formed/README.md)

The 99-state native operator on a written magnitude-two collar acts on its coin and sign alone throughout the zero- and one-probe sector. The matter cannot cross it. The blank seams transmit. Each supplies the usual exchange and contrast costs, including the fixed blank degree references. Factoring out the freely evolving memory gives the one-matter operator
$$
H=\tfrac12 L+\kappa(6 I-R), \qquad R_{vv}=r_v ,
$$
where $L$ is the Laplacian of transmitting contacts. The zero-momentum source used here is constant over one transmitting component. It can be synthesized from the home port by native splitters and phases. On each transmitting bond its two endpoint amplitudes coincide. The source, reversal, staging and isolated clock are apparatus inputs. [Local factorization and experiment](experiment.py) · [Source synthesis and spatial bond read](../row_30_formed/experiment.py)

For any eigenstate of this spatial observable, $(H-\langle H\rangle)\psi=-\kappa(R-\langle R\rangle)\psi$. Thus its squared frequency width is exactly $\sigma_H^2=\kappa^2\operatorname{Var}_\psi(r)$. For the uniform wave an observer can forecast the width from counts of touched and untouched sites in the component. **Every one of the 940 components in the 27 histories contained both kinds of port.** The forecast widths squared on the 27 tested largest components ranged from 9.5822 to 14.7654 in exchange units squared.

At a short time $t$ the return loss lies between $\sigma_H^2 t^2(1-W^2t^2/12)$ and $\sigma_H^2 t^2$, where $W=14$ bounds the spectral width for the registered coupling $\kappa=8$. This follows by applying the cosine remainder bound to two independent energy draws. The exchange Laplacian divided by two has spectrum in $[0,6]$, since every transmitting port has degree at most six. The diagonal shift from the writes spans at most eight. The home observer jointly reads its dimer-clock bit, calibrated to $\sin^2(t/2)$ on the stated arc, and its return bit. All short-time forecast intervals overlap the measured intervals with simultaneous 99% counting coverage across all settings. These are conditional simulation intervals. The range of losses is across distinct preparations. The changed law removes the contrast term while retaining the same transmitting exchange.

| Fresh seed 301121 | Ports in largest component | Once-written ports | Home misses at 1/16 | Home misses at 1/4 |
| --- | ---: | ---: | ---: | ---: |
| 20 cells | 120 | 92 | 182,009 | 1,846,163 |
| 80 cells | 612 | 453 | 195,349 | 1,987,992 |
| 320 cells | 504 | 384 | 184,763 | 1,855,206 |

A pure joint eigenline of $H$ and $L$ must be supported at one value of $r$. For a family in one such class, $E=\kappa(6-r)+\lambda/2$, where the independently specified squared wave number is $k^2=\lambda$. For any fixed reference, subtracting $m^2+c^2\lambda$ from the squared frequency leaves a quadratic with leading coefficient 1/4. Thus an exact massive Klein–Gordon relation cannot hold at three distinct spatial eigenvalues within the class. For finite widths the triangle inequality gives $\kappa\sigma_R\leq\sigma_H+\sigma_L/2$. The candidate might still support a different long-wave envelope with a coarser spatial observable or another excitation. This result supplies a measurement of why the registered one-load test on this independently formed medium misses the *rest* line. The positive prepared waveguide result and its identified range of accuracy stand on their original preparation. [Prepared four-site channel](../row_30_autonomous/README.md) · [Public criterion](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L136)

## Reproduction

Run from the RER repository root. The default command rebuilds the histories and compares the recorded counts. The tests also propagate the probe with a separately assembled dense matrix reconstructed from archived repair signs.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_repair.experiment
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_repair
```
