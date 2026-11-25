# What fixes the eleven written seams?

On the M10.4 board, eleven writes are the **minimum necessary to protect the filed two-load interference state** at its address, among all blank/neutral assignments. Eight suffice for its one-load orbital. This identifies what each write does. The [OpenWave review](https://github.com/openwave-labs/openwave/pull/749#pullrequestreview-6060088147) kept row 4 at partial because the authors selected the pattern with persistence in view. Neither of two material prescriptions older than that selection reproduces the result. The uniform commissioning service used for row 30 can produce the pattern, but at the filed address it does so about once per 13.3 million boards.

## Which contacts does the state need?

Let $a=157$, $b=156$, $c=166$ and $d=163$. The two tips $a,c$ have exactly two common neighbors, $b,d$. A blank seam exchanges one load with an empty endpoint. A neutral write switches that term off and stays at its address under free evolution. The private degree terms retain their original values when a contact is written. On a binary word $S$, each transmitting seam with one occupied end adds $|S\rangle-|S\mathbin{\triangle}e\rangle$ to twice the exchange column and one to the diagonal contrast. The degree reference adds $\sum_{v\in S}(6-\deg v)$ to the latter. These columns include every destination outside the source and the home cell.

For the one-load ray $|a\rangle-|c\rangle$, the eight tip contacts that have no partner must be written. At each shared neighbor, both incident tip contacts must either transmit or both block. This gives equal diagonal terms too. For the two-load ray $|ab\rangle-|bc\rangle$, moving the center load along any of its three other contacts gives distinct destinations whose amplitudes cannot cancel. Those three contacts must also be written. Matching the two preparation routes and the two exit routes gives four possible pairs of switch settings. The opposite sign adds the two amplitudes at an open common destination.

| Common eigenray of exchange and contrast | Passing masks out of 32,768 | Fewest writes | Passing with both preparation and both exit contacts open |
|---|---:|---:|---:|
| One load, minus | 32 | 8 | 8 |
| One load, plus | 8 | 12 | 0 |
| Two loads, minus | 4 | 11 | 1 |
| Two loads, plus | 1 | 15 | 0 |

The four two-load solutions use 11, 13, 13 and 15 writes. The unique minimum is the filed eleven. This classification imposes **no routing-task clauses, eigenvalue or persistence threshold** in the operator test. It asks for this specific ray, whose sign and support were selected for the result. Other persistent states and other preparations are outside the classification. [Integer action and all mask counts](masks.py) · [Independent exit-by-exit check](test_source.py)

The center by itself is mobile. At zero contrast and exposure 37 a lone center load returns to port 156 with probability **0.121087**, while the antisymmetric tip load has probability one on its two ports. Two independent copies, read with the same unordered geometric seam criterion, give **0.233285** home; the ideal two-load dark ray gives one. At exposure 150 the independent-copy probability is 0.229043. The mechanism claims persistence of joint records relative to the written material. It requires no attractive contrast energy. The independent copies are a constituent control, distinct from the acquired two-load state and its finite failures. [Constituent reads](flight.py) · [Home-port and pair-address records](results.json)

## How often does a different experiment produce it?

The [spectroscopy commissioning rule](../../row_30_formed/PROTOCOL.md) takes two fresh fair readings per empty contact and writes if both are positive. The write probability is $p=1/4$. This is an addressed service on the native law. Free evolution of an all-blank board alone has zero probability of making even the first neutral mark. The spectroscopy rule predates these test flights and postdates the M10.4 design.

The probability at the original address of eleven writes and four openings, **allowing arbitrary material everywhere else**, is $p^{11}(1-p)^4=81/4^{15}$. Requiring the other 615 contacts also to be blank gives $1.09714\times10^{-84}$ for the entire filed board. On the local event, the same pulses, read, retry and dark columns apply regardless of exterior writes. The joint probability of that event and a home read at any exposure is at least $7.50427\times10^{-8}$, including acquisition exhaustion. Other materials can also produce home reads.

| Cells | Eligible local layouts | Bounds on probability of at least one matching layout anywhere |
|---:|---:|---:|
| 20 | 104 | 0.00000760025 to 0.00000784547 |
| 80 | 384 | 0.0000280621 to 0.0000289679 |
| 320 | 1,504 | 0.000109905 to 0.000113458 |

The upper bound sums the single-event probabilities. The lower subtracts each pairwise intersection once, treating shared contacts with conflicting requirements as impossible. These are bounds on local motifs, not on a chosen launch's home detection or on whole-board connectivity.

Even if the write count is correlated, a source invariant under **every permutation of the physical seams** is a mixture of uniform $k$-write placements. Among $N$ contacts its probability for this specified local mask is bounded by

$$\max_k\frac{\binom{k}{11}\binom{N-k}{4}}{\binom{N}{15}\binom{15}{11}}.$$

On 630 seams the sharp maximum is **0.000168818** (rounded), attained at 462 writes. This exact counting consequence covers any exchangeable write budget. Sources that use the local geometry or the launch history can break permutation symmetry. The bound does not apply to them.

## What the observer reads without selecting the material

There is an exact obstruction on the unconditioned spectroscopy source at zero contrast. With probability $9/16$ both preparation contacts are blank and the eight-attempt exhaustion is $(364/729)^8$. With probability $3/16$ only the contact from the center to tip $c$ transmits. Each occupied read then succeeds with probability $\cos^2(\pi/3)=1/4$; three successive reads accept only $1/64$ per attempt. In the remaining quarter the contact to $c$ is written and the sharp initial pair never leaves during acquisition. Thus **at every free exposure**,

$$P(\text{home})\le P(\text{accepted})=
1-\tfrac9{16}(364/729)^8-\tfrac3{16}(63/64)^8
=0.832521775011\ldots.$$

This bound uses the donor's coin law and the filed source budget. It requires neither flight propagation nor a sampling assumption about the 128 archived banks. [Exact source branches and check against every bank](admission.py)

In the [fixed flight test](PROTOCOL.md), all 128 successive material banks were retained. At zero contrast and exposure 7 the mean home probability is **0.230586506** (range 0.013125 to 0.566747 rounded); same-material random swaps give **0.166572611**. The one-sided 99% Hoeffding upper limit is 0.364710 under independent material draws, including a $10^{-7}$ numerical allowance. These are computed detector probabilities per material, not counts of simulated particle detections. Chebyshev propagation uses an analytic truncation tolerance of $10^{-17}$; a separate exponential checks the first bank and each historical transfer. The allowance is separate from that tail and from sampling variation. The rational admission ceiling is a separate result.

The following two neutral masks were fixed for purposes other than persistence **before the pair layout**. Only their masks are transferred. Load preparations, parked station residents and instruments of the donor tasks are not identified with the M10.4 pair. All entries use strength ½ and exposure 37. The maximum is the largest eigenvalue of the home-read effect on the two source words, allowing perfect preparation at arbitrary relative phase.

| Material in its original location | Acquired pair | Filed input supplied | Best possible input on the two words | Random swaps on actual input |
|---|---:|---:|---:|---:|
| [First external scattering contact](https://github.com/FloatingPragma/observer-patch-holography/blob/44f7680df55b4e68b89c301896df3e10237f046d/model/candidate_laws/interval_echo/openwave/row_30_memory/PROTOCOL.md) | 0.440925 | 0.440925 | 0.442701 | 0.009313 |
| [Tidal station's 160 sham guards](https://github.com/FloatingPragma/observer-patch-holography/blob/d37d200a794b49446433ac5c7306f241fed18295/model/experiments/interval_echo_einstein_station/native_tidal_response/PROTOCOL.md) | 0.120211 | 0.159186 | 0.171024 | 0.065917 |
| Blank control | 0.440929 | 0.440929 | 0.442706 | 0.009313 |
| Filed eleven | 0.995945 | 0.995945 | 1 | 0.123196 |

The tidal rule selects the center address for a matched-ray experiment. It closes both preparation contacts there and leaves both required external contacts open. The two archived spectroscopy materials have their own [matched flights](../commissioned_medium/README.md). None of these tests establishes a source independent of persistence that both selects the filed material and passes the pair measurement. A candidate for that role needs its original selection record and its actual complete output, including cases that do not prepare the pair. A geometry-dependent service can be tested this way without requiring it to reproduce the filed mask.

**A two-contact test for another source.** Let $u$ be its frequency of both pulse contacts open, and $v$ its frequency of only the center-to-$c$ pulse contact open, before placing the loads. The other 628 contacts can be arbitrarily correlated. A necessary condition for 94% inclusive home detection at zero contrast with this acquisition is
$u(364/729)^8+v(63/64)^8\leq0.06$.
Feed all four counts of an independently fixed material panel to `seam_origin.admission --panel BOTH_OPEN ONLY_Z ONLY_X NEITHER`. Failure excludes the **panel** at every exposure; passing only justifies doing the flights. Finite panel frequencies need a separate sampling argument for a source-population claim. The blank control illustrates the distinction: it passes the source test and loses the pair in flight.

## Replay

From the repository root, with Python, NumPy and SciPy:

```sh
python -B -m model.candidate_laws.interval_echo.openwave.row_04.seam_origin
python -B -m model.candidate_laws.interval_echo.openwave.row_04.seam_origin.admission
python -B -m model.candidate_laws.interval_echo.openwave.row_04.seam_origin --full
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_04/seam_origin
```

The first command verifies the exact census on all three boards. The second derives the admission ceiling. Full replay checks the immutable donor source blobs, regenerates all 128 two-bit histories, propagates all 28,680 two-load configurations per material and compares each terminal bin with the archive. The pulse and terminal instruments reuse the [native material transfer](../commissioned_medium/experiment.py). The raw-board test constructs its own pair graph, pulse and two propagations for the two pre-layout donors and the first ensemble bank. Verification reads the receipt. `--record PATH` creates a separate first record and refuses an existing path.
