# A spectrometer's excitation stays on two ports in free flight

An OPH observer can start with one sharply placed load and blank sample memories, prepare either of two previously specified wavelength analyzers, and release it on a connected board. At **every free exposure and every real contrast strength**, the probability of finding the load on port 9 or 11 is at least **0.999997291223**, counting every unsuccessful preparation. Both ports have an exit. Changing the relative sign gives returns of 0.461707 and 0.398056 at exposure 37 on the 20-cell board. [Original spectroscopy task](../row_30_selective/PROTOCOL.md) · [Complete calculated probabilities](results.json) · [Preparation and tests](DESIGN.md)

## Why these memories and these states?

In [the public donor at `44f7680df`](https://github.com/FloatingPragma/observer-patch-holography/tree/44f7680df55b4e68b89c301896df3e10237f046d/model/candidate_laws/interval_echo/openwave/row_30_selective), the question was whether a material excitation had both a resolved frequency and a resolved wavelength. Each physical seam was given an arrow pointing out of its earlier port in the ascending address order. The source started at port 0. A hop reverses the participating arrow. The spectrometer tested 32 wavelength analyzers against every one of 31 distinct energies. Its primary 90% overlap criterion returned a single match: squared wave number 1, in a two-dimensional space at energy 4 above the chosen energy origin when $\kappa=8$. Both analyzer columns, 14 and 15, are used here. No home-detection probabilities select among them. The arrow rule, spectral criterion, pulse archive and published match precede the [M10.4 pair design](https://github.com/FloatingPragma/observer-patch-holography/blob/e62ead94322952b250918d4dc78137d3f82f34d0/model/candidate_laws/interval_echo/openwave/row_04/PROTOCOL.md).

The excitation is one unit **together with its coherent change to the written arrows**. Its two components end at physical ports 9 and 11, with opposite amplitudes. The first analyzer reaches their common parent by the path $0,5$ and the second by $0,1,5$. These different preceding arrow records are part of the state. The complete reachable material on the donor board has 32 configurations in three cells. No bond is shuttered during flight. This object, its matter source and its two-port measurement differ from the two-load, eleven-neutral-write candidate in [OpenWave #749](https://github.com/openwave-labs/openwave/pull/749#pullrequestreview-6060088147).

## From a blank sample to the donor's material

One controllable local fuel converter writes the historical symmetric superposition of the two contrast signs of magnitude one. Its coin points out of the specified endpoint. Reversing the mount reverses endpoints, sign and coin together and preserves the full native seam law. A fuel measurement distinguishes a completed conversion from a blank seam. In the matched apparatus the occupied fuel rail gets the sample's native degree-reference energy. On the exact two-state input the addressed generator is

$$\begin{pmatrix}\kappa&G\\G&\kappa\end{pmatrix},\qquad G=64.$$

Both diagonal entries follow from the local operators, on different registers. The fuel is a preparation accessory. Its gain and pulse area $\pi/4$ come from the [collar reader](../../../contrast_memory_spherical_v1/collar.py). The transition probability is exactly $1/2$, independent of $\kappa$. A refusal leaves both endpoints vacant and the seam blank. A converter on bare, zero-energy fuel instead gives $\sin^2(\pi r/4)/r^2$, $r=\sqrt{1+(\kappa/2G)^2}$. This is at least $1/4$ on $|\kappa|\le128$, including all the magnitudes in the historical coupling sweep. No stronger uniform bound on this bare instrument is assumed. Other sample contacts are isolated during preparation.

The single load is kept on the board throughout. Two ordinary sharp swaps park it out of the way; a third, measured hop on a provisionally reversed arrow returns it to the root and sets that arrow to its required orientation. Every conversion occurs with vacant endpoints. All ten ordered choices of working triangle around port 0 lead to the same installed material at each size. The final free sample contains only the original load, native degree references and the material specified by the spectrometer. [Formation sequence](source.py) · [Full 99-state seam and 198-state fuel checks](local.py)

## A source that does not require the correct coupling as input

The donor's old pulse durations depend on the material frequency. At $\kappa=8$ they prepare the two rays with protected weights above 0.999999999997. At $\kappa=1/2$ the same waits give two-port probabilities of **0.328879 and 0.378952** at exposure 37. At zero the load never leaves its initial port. Supplying the ideal wave after this failed preparation would change the experiment.

Here an addressed native echo, a binary read of blank memory, and the inverse echo measure a projector $C_a=|p+a\rangle\langle p+a|/2$ on one allowed parent–child hop. The other child passes unchanged through its no branch. Reading the coin afterward resolves a sharp transfer, with probability $1/2$ from either end. Bounded such transfers reach the specified parent. In the coherent part of the preparation let $Q=I-|p\rangle\langle p|$ denote the local center-vacancy query against a blank degree reference. In chronological order query $C_a,C_b,I-Q,C_a,I-Q$. For the recorded answers yes, no, no, no, no the unnormalized output from $p$ is

$$ Q(I-C_a)Q(I-C_b)C_a|p\rangle=(|a\rangle-|b\rangle)/4.$$

Accepting no on the first query instead gives the positive-sign ray up to global phase, with the same $1/8$ probability. Each first refusal is recorded. A sharp read of the two terminal head bits places a refused source at $p,a$ or $b$; a bounded return can fail and then terminates the trial. The complete first-refusal probabilities, resolved by location, are in the receipt. Preparation uses no numerical value of the coupling. The observation appointments are supplied native exposures; no conversion to laboratory time is claimed.

With retry cap $N$, $E$ physical seams and $\ell=1$ or $2$ prefix hops, the original-start detection probability for the negative-sign source obeys

$$P_{\rm home}(t)\ \ge\ 1-(7/8)^N-(E+\ell+N)2^{-N}.$$

There are $E$ potential write timeouts, one commissioning hop, $\ell$ prefix hops and at most $N-1$ returns. The separate $(7/8)^N$ term accounts for no successful five-query word. For bare fuel replace $E2^{-N}$ by $E(3/4)^N$ on the stated strength range. The bound may be uninformative at small caps; the branch recurrence computes actual probabilities there. At $N=96$ the dominant loss is about $2.709\times10^{-6}$ from exhaustion. On the largest board the matched lower bound clears 0.94 first at $N=22$ and 0.98 at $N=30$. The conservative bare-fuel lower bound clears them at 43 and 47, respectively. The budget controls yield; an accepted source has the same signed ray for any positive cap.

## Why the load stays, and what escapes

Write $A$ for the adjacency of actually flippable histories, $d$ for their degrees and $r=6-\deg_{\rm physical}$ for the original blank reference at the occupied port. Apart from the common seam energy, the free Hamiltonian on this sector is

$$H_\kappa=(\kappa/2)F,\qquad F=2\,{\rm diag}(r)-{\rm diag}(d)-A.$$

This formula applies at $\kappa=0$ without division. Each target consists of two terminal histories with the same parent and the same reference. Their signed difference $D$ satisfies $FD=D$ and $({\rm diag}(d)-A)D=D$ in integers. Their amplitudes cancel in the parent row; all other rows off the two tips vanish. The port observation contains the whole selected two-dimensional eigenspace. Any normalized preparation of fidelity $f$ within that space, including a mixture or one entangled with a reference, has home probability at least $f$ at every exposure. Thus persistence depends on the donor-selected material and the phase the source delivers, rather than the particular method of compiling that phase.

| Home probability at $\kappa=8$, exposure 37, including failures | Donor column 14 | Donor column 15 |
|---|---:|---:|
| Reconstructed signed analyzer | 0.999997291224 | 0.999997291224 |
| Opposite phase, same acquisition budget | 0.461706674836 | 0.398056302547 |
| Same prepared state, dephased at release | 0.730851983030 | 0.699026796886 |
| Same populations, rate-matched independent swaps | 0.374998776339 | 0.374999202737 |

The positive-phase controls cross the source cell boundary with probabilities 0.073053 and 0.044729 at that appointment. At zero coupling both flight laws are stationary and their agreement carries no discrimination. The two alternate orderings of the donor do not reach the named incoming histories. The donor also reported a staggered definition of spatial phase, with five spectral matches instead of one; it is separate from the primary wavelength measurement.

| Native cells | Written seams | Complete reachable configurations | Worst all-exposure lower bound, bare fuel, magnitude at most 128 |
|---:|---:|---:|---:|
| 20 | 630 | 32 | 0.999997290585 |
| 80 | 2,520 | 70 | 0.999997288669 |
| 320 | 10,080 | 70 | 0.999997281006 |

Only the first board supplied the donor's original spectral selection. The other two test the same fixed preparations on different physical carriers. All physical port probabilities are retained using a sparse address list with exact zeros on unreachable ports. Source refusal categories are disjoint from terminal address reports. Summing the tip bins gives the observer's two-port measurement; a simulator path is never an input to that reader.

In a separate source control, the same five interrogations prepare both signs at **every** pair of accessible children with a degree-five common parent. On 20 cells there are 14 such forks. All twelve outside the donor's two fall below 94% at strength 8 and exposure 37, despite the same exact preparation yield. The lowest signed return is 0.099957. On each larger board there are 36 forks; six exact signed eigenrays include **four not enrolled by the donor**. For example a signed wave on ports 45 and 47 at the same energy gives zero reading on 9 and 11. The persistence bound concerns the two transferred rays, not the full refined energy band. [All-fork test](forks.py) · [Source branches at five budgets](budgets.json)

## A second route using the original pulse programs

A separate [clock experiment](CLOCK_TEST.md) tests whether the new five-query program is essential to localization. An observer times a local marked contact and reads the first resolved transfer on its increasing half-cycle. Its counts, with no detector or target-state information, set the duration unit for **all 32 original spectroscopy words**. Across sixteen signed strengths from magnitude ½ to 128, all 512 prepared waves pass the donor's 98% fidelity criterion. The lowest simultaneous fidelity lower estimate is 0.99993747 under the declared binomial calibration model and numerical allowance. Every enrolled two-port return at exposure 37 exceeds 0.99999962 conditional on material preparation. The zero-strength calibration correctly gives no clock. The principal five-query result includes that case without extrapolating a clock through zero. These are counts from simulated calibration trials and probabilities for the terminal flight. The [receipt](clock.json) also records the numerical norm drift separately.

## Reproduction

From the repository root with the `model/requirements.txt` environment:

```sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.experiment
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.budgets
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.transfer_clock
python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_04_spectral
```

The first command compares the pinned donor, reconstructed material, complete reachable sectors, source failures and loss controls with the archived probabilities. `--record PATH` creates a first record at a fresh path. The [choice ledger](CANDIDATE.md) distinguishes conditions fixed by the old spectrometer from the preparation apparatus and disclosures for a candidate decision.
