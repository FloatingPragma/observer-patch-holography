# Reading a bit at a neighbouring cell

At supplied exposure 4, one terminal reading recovered 45,559 of 65,536 late-chosen bits on a 20-cell Interval Echo screen and 47,252 of 65,536 on a matched 320-cell station. Both totals include every quiet receiver. At exposure 0.25 the sampled accuracy intervals include chance. A separate stationary clock reached both requested tick counts. Its hand occupies the receiving contact, and the candidate exchange moves that hand out of its prepared register sector. Neither deadline produced a message answer. [Frozen protocol and source hashes](protocol.json) · [Complete measurements](receipt.json)

## Preparation and receiver

The apparatus supplies a sharp unit load on a common sender port, blank seam memories, the reference vacancies and a prior appointment for release and readout. At setup the receiver is empty. After setup a balanced bit selects one of two addressed sender seams. The same π pulse of the native reduced seam term takes the common load to the selected port. Each encoding uses one contact inside the same sender cell. The supply and gate window are apparatus inputs. The simulation calculates the two pulse columns once, then samples independent, freshly prepared releases. No detector watches between release and the terminal gate. [Sender and detector](channel.py)

| Cells | Sender cell | Common port | Bit 0 / bit 1 ports | Receiving cell and contact |
| ---: | ---: | ---: | --- | --- |
| 20 | 6 | 78 | 82 / 74 | cell 0, seam 82–1 |
| 320 | 3 | 40 | 38 / 39 | cell 0, seam 38–1 |

The station map fixes the first receiving external seam, the facing sender port and its lowest-numbered internal degree-five neighbour. These are matched one-link installations. Their physical separation across refinements is uncalibrated. The 320-cell calculation needs two complex columns of length 3,840 (120 KiB together) and an operator with at most 24,000 stored entries. Both boards retain all their links throughout free flight.

For bit $b$, let $u_b(s)=\exp(-isH)u_b(0)$ after the encoding pulse. The cell records one of its twelve ports with probability $|u_{b,j}|^2$, or quiet with probability $1-\sum_{j=0}^{11}|u_{b,j}|^2$. The complete columns give $T=\tfrac12\sum_y|p_0(y)-p_1(y)|$ and optimal success $(1+T)/2$ for this classical instrument with equal priors. This last number is a model-informed benchmark. The receiver learns a separate table for each appointed window from 32,768 labelled preparations per bit. It guesses the more frequent calibration label at the reported port or quiet, breaking ties toward zero. Its only inputs are a plain outcome and the frozen public table. The simulator draws outcomes, and the scorer joins decisions to truth. The receiver has no file, RNG, trial index or amplitude input. [Decoder](decoder.py) · [Acquisition](acquisition.py)

## All terminal windows

Each row has **32,768 evaluation trials per bit**, 65,536 in total, separate from calibration. Detection probabilities are unconditional and ordered by bit. Accuracy and its simultaneous 95% bounds are percentages. The final column counts quiet trials answering correctly and incorrectly. Total errors equal 65,536 minus the correct count. Every receiver guesses, so coverage is one and abstentions are zero. The receipt gives both 13-outcome laws, integer counts by bit and outcome, and calibration counts for every window.

| Cells | Exposure | Detection, 0 / 1 | TV | Optimal % | Correct | Sample accuracy % [bounds] | Quiet correct / wrong |
| ---: | ---: | --- | ---: | ---: | ---: | --- | --- |
| 20 | 0 | 0 / 0 | 0 | 50.000 | 32768 | 50.000 [49.265, 50.735] | 32768 / 32768 |
| 20 | 0.25 | 0.015050 / 0.000059 | 0.014991 | 50.750 | 33232 | 50.708 [49.973, 51.443] | 32765 / 32301 |
| 20 | 1 | 0.146347 / 0.008510 | 0.137837 | 56.892 | 37393 | 57.057 [56.322, 57.793] | 32463 / 27838 |
| 20 | 2 | 0.270953 / 0.028102 | 0.242851 | 62.143 | 40816 | 62.280 [61.545, 63.016] | 31858 / 23810 |
| 20 | 4 | 0.425134 / 0.035519 | 0.389615 | 69.481 | 45559 | 69.518 [68.782, 70.253] | 31652 / 18861 |
| 320 | 0 | 0 / 0 | 0 | 50.000 | 32768 | 50.000 [49.265, 50.735] | 32768 / 32768 |
| 320 | 0.25 | 0.015050 / 0.000059 | 0.014991 | 50.750 | 33230 | 50.705 [49.969, 51.440] | 32766 / 32304 |
| 320 | 1 | 0.147071 / 0.008521 | 0.138551 | 56.928 | 37404 | 57.074 [56.338, 57.809] | 32509 / 27873 |
| 320 | 2 | 0.285686 / 0.028772 | 0.256914 | 62.846 | 41222 | 62.900 [62.164, 63.635] | 31870 / 23416 |
| 320 | 4 | 0.472289 / 0.031855 | 0.440434 | 72.022 | 47252 | 72.101 [71.365, 72.836] | 31705 / 17221 |

**Stage A: CHANNEL_RESOLVED** at both stations. The fixed test required the lower bound at 4 to exceed 52% and the upper bound at 0 to fall below 52%. At 0.25 the upper bounds on the tested decoder's advantage over chance are 1.443 and 1.440 percentage points. At zero the disjoint support gives quiet with probability one analytically. The grid gives no first-information time.

The 20-cell exact distributions on this grid were examined in development, including the weak 0.25 window. The geometry, grid, sample size and decisions were frozen locally at `dc573d9d6` before the evaluation sample. The receipt reproduces the first evaluation counts. Enrolment, label assignment, calibration, detector and label-permutation streams are specified in the protocol. The calibration labels are public only in the training cohort. For fixed calibration and balanced labels the correct indicators are independent bounded draws conditional on the assignment. Across 30 comparisons, the ordinary Hoeffding half-width is $\sqrt{\log(60/0.05)/(2N)}=0.007355$. For fixed guesses and permuted balanced truth, accuracy differs from 1/2 by the displacement of the fraction of positive guesses in a sample of $N/2$ without replacement from its population fraction. The resulting half-width is 0.010401 (rounded). A true 53.5% accuracy misses the 52% lower-bound threshold with probability below 0.00048 under the stated sampling model.

## What changes the reading

- **No encoding.** Both labels select the same sender pulse. Identical laws give zero total variation and exactly chance balanced success for any fixed table. Fresh sham transcripts scored 49.744–50.237%. Independently permuting truth on the native transcripts gave 49.789–50.295%. Every simultaneous control interval contains 50%.
- **Contrast removal.** An occupied port of degree $d$ contributes $d\kappa$ from incident contrasts and $(6-d)\kappa$ from references. Thus $H=L/2+6\kappa I$ in the recovered one-load sector. At κ = 8 the scalar is 48. The literal generator difference has zero residual after subtracting that scalar, and removing it changes the terminal outcomes by at most $2.23\times10^{-16}$. These data resolve exchange transport. Attribution to contrast binding needs a sector with interacting loads and a matched complete detector. The elapsed-word filter is idle here.
- **Moving alternative.** Uniformly halving the exchange strength preserves incidence, locality and symmetry and gives $H'=L/4+48I$. At fixed exposure 2 its law equals the native exposure-1 law. The largest changes in any outcome at exposure 2 are 0.124606 and 0.138615, above the fixed 0.02 sensitivity gate. Rescaling the exchange rescales this audit coordinate.
- **Planted signal and integrity.** At zero exposure, two distinguishable loads placed directly at receiver ports 0 and 1 in separate preparations give TV and planted-reader success one. Tests reject label-dependent quiet removal, a missing complement, an injected label, a permuted detector address and replacement of the native distribution by the slower law. Only the physical comparisons enter the channel interpretation. [Focused tests](test_message.py)

The producer uses sparse exponential actions, records the raw state mass minus one and checks the quiet branch against the mass on all exterior ports. It rejects negative probabilities and never rescales columns. The largest absolute native mass residual is $1.12\times10^{-15}$ or smaller. An independent implementation builds oriented incidence from the pinned edge file, forms $B^TB/2$, drops the common phase and uses an order-18 Taylor recurrence with steps of at most 1/8. In exact arithmetic each step has operator-norm truncation error at most $e^{0.75}0.75^{19}/19!$ since $\|L/2\|\le6$; over 32 steps the accumulated bound is below $2.4\times10^{-18}$. Both calculations use floating point. The largest observed probability difference is $1.56\times10^{-15}$ or smaller. The $10^{-10}$ per-outcome numerical tolerance allows for roundoff in both implementations and implies at most $6.5\times10^{-10}$ uncertainty in the classical optimum if each entry is within tolerance. It is an engineering allowance, with no certified bound on propagation tails. [Independent calculation and count check](headline.py)

## What the stationary clock measured

**Stage B: JOIN_INCOMPLETE.** The 20-cell receiver uses the stationary fair double on rails 2 and 3, the five-position hand in cell 0, its nine-bit counter in cells 4, 1 and 3, and the published digit on rails 54 and 55 read against vacancy 58. Each clock probe records zero at its supplied reset appointment. A joint run would release the message after this appointment. Each independent preparation requests a detector reading on reaching its named count, capped at eight local opportunities. The clock publishes zeros as well as ones. With a one-tick deadline it read `0,1`; with a two-tick deadline it read `1,0,1`. Both requests reached the readout interface. Both returned unavailable, with zero commitments and zero joint message trials. These are isolated material records. No conversion from the supplied exposure was used. [Material attempt](clock_join.py) · [Clock publication](../interval_echo_clock/PUBLICATION.md)

At release a double occupies hand port 1. The sharp sender at port 82 has occupation 1 or 0 according to its bit. On this particular seam the implemented local generator has matrix element $\langle n,2|h|2,n\rangle=-1/2$ for both $n=0,1$. No other seam has this exchanged configuration as its output from that input. The intact free law therefore leaves the fixed ready-hand sector at first order. Protecting even this one contact changes the pilot's one-load receiver law by 0.270953 in a single outcome at exposure 2. An additional isolated-gate check applies the first π contact of a hand advance to the free signal at exposure 2. This exchanges readings at ports 7 and 10, changing a port probability by 0.007111 for bit 0 and 0.001462 for bit 1. These are diagnostics of the proposed composition.

A joint source map requires the signal and doubles on the same load registers, the free and addressed gate windows on that state, and a detector that distinguishes the signal in the occupied clock patch. A controller would have to order publication, counter read and terminal detection, carry all unavailable and overflow branches, and commit the estimate locally. Protecting or moving the hand changes the tested apparatus and calls for recomputing the channel. The optical ranging apparatus uses a reciprocal guided double and a triggered reflection. Its calibration has no identified transfer to this free signal. [Optical itinerary](../interval_echo_optics/README.md)

## Relation to existing experiments

At OPH `6462c2494a6c4f328f730004362af4962cb89954`, spatial tomography pools multiple local ensembles to infer source mixtures, and its forecast uses a local survey. The contact-message experiment distinguishes a capture or release decision from blocks of scalar returns with supplied slots. The torus OBS-03 experiment includes a successful early addressed relay. Here each guess uses one terminal outcome at one appointed station. The clock comparison, relay menu, independently ranged deadlines and transfer to a physical speed remain unresolved. Local OBS-03 keeps its open grade. [Spatial experiment](../interval_echo_space/README.md) · [Contact message](../contact/README.md) · [Observer criterion](../../evaluation/EVIDENCE.md#obs-03)

The frozen acquisition uses OPH `6462c2494`. A full replay at `0ebb9ba9d6280fec258a7d65e2966c13352dfa8d` reproduces the deposited receipt byte for byte. Three pinned files differ. The candidate module adds particle-identification metadata, the specification describes filled-background defect experiments, and the evaluation records particle and dispersion evidence. The one-load generator, carrier incidence, local reader and stationary hand have their original definitions. The verifier accepts either complete source manifest. [Source checks](channel.py) · [Defect comparison](../interval_echo_encounter/README.md)

The upstream dispersion experiment resolves a band on a prepared 194-cell record channel using clock collars and a local matter standard. Its detector reads an ensemble phase at supplied windows. This leaves the single-trial deadline on the ready hand at contact 82–1 unresolved. [Prepared channel](../../candidate_laws/interval_echo/openwave/row_30_autonomous/README.md) · [Clock readout](../../candidate_laws/interval_echo/openwave/row_30_arrows/CLOCK.md)

OpenWave adopted Interval Echo at `4fe10f9a` and awarded partial credit on row 4 for persistence of a prepared pair. Local OBS-03 concerns a separate receiver-clock measurement. [Maintainer decision](https://github.com/openwave-labs/openwave/pull/745) · [M10 specification](https://github.com/openwave-labs/openwave/blob/753232cac4976ab2ca208eb99d73fe7ac88ef333/openwave/xperiments/m10_oph/research/m10_theory_canonical.md)

## Offline replay

From the repository root, using Python 3.14.4, NumPy 2.4.4 and SciPy 1.17.1 in the recorded run:

```sh
export OPENBLAS_NUM_THREADS=1
python -B -O -W error -m model.experiments.interval_echo_message.headline
python -B -O -W error -m model.experiments.interval_echo_message
python -B -m pytest -q model/experiments/interval_echo_message
```

The headline rebuilds all ten native distributions and re-scores the counts without importing the producer. The full replay also regenerates calibration, control counts and the isolated material records. Verification reads the deposited reference. For first-record creation, pass `--record PATH` to the full replay with an unused destination. With one BLAS thread the headline took 1.48 s and full replay 2.74 s, both with exit status zero. The two implementations share the incidence data and protocol and were written in this study. External review, other stations, the 80-cell carrier, multi-load flight with clock hardware, unrestricted relays and a continuum limit were outside the measured domain.
