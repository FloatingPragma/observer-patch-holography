# How many screen records can one or two mobile loads reach?

Six connected twelve-port cells retain six possible owner parity words for one mobile unit load and sixteen for two. These counts are attained at generic positive times by the unswitched blank-seam Hamiltonian. At the registered time eight, two outside stations resolved five words in 65,536 simulated launches for each inventory. The finite probability of a record matters when inferring accessible capacity from an observed histogram. [Frozen protocol](PROTOCOL.md) · [Full two-station receipt](records.json.gz)

## Capacity at fixed binary inventory

Let a connected closed board have $C$ cells, twelve binary ports per cell and exactly $K$ unit loads. A reported one is an odd cell population. If the word has $r$ ones, each odd cell requires at least one load and at least one vacancy, and the sum of its parities is $K\bmod 2$. Conversely, put one load in each odd cell. The remaining $K-r$ loads can be distributed in pairs. Each odd cell has five available pair places and each even cell six, giving $6C-r$ places. Thus every word satisfying
\[
  r\equiv K\pmod{2},\qquad r\le\min(C,K,12C-K)
\]
has a material configuration. The number of possible words is
\[
  N(C,K)=\sum_{\substack{0\le r\le\min(C,K,12C-K)\\ r\equiv K\,({\rm mod}\,2)}}\binom{C}{r}.
\]
For $K=1$, this gives $C$. For $K=2$, it gives $1+\binom C2$. For a whole closed 20-cell board these numbers are 20 and 191. At 80 cells they are 80 and 3,161, and at 320 cells they are 320 and 51,041. These whole-board numbers are exact consequences of the finite source law. The measured flight here encloses six cells. Even at the largest binary inventory support for a fixed $K$, total parity restricts the accessible words to a single parity class. Sampling or selecting a physical screen with independently variable inventory is an additional input to PR-124.

## Reachability under simultaneous free evolution

On a blank seam with occupancies 01 or 10, the off-diagonal exchange element is $-1/2$ and the diagonal element is $\kappa+1/2=8.5$. Blank 00 and 11 are stationary on that seam. Summing **all** open seam terms on the region therefore gives the binary exclusion graph Hamiltonian. Its diagonal is $8.5$ times the number of unlike neighbouring pairs and each permitted particle-vacancy exchange has matrix element $-1/2$. No sequential choice of transfers is made during flight. [Native columns and sparse Hamiltonian](measure.py)

The configuration graph at fixed $K$ is connected. Choose a spanning tree and fix a leaf to its desired load. When its present load differs, find the nearest particle (or vacancy) of the required kind elsewhere in the remaining tree. Interior vertices of that path have the opposite kind, so it can move to the leaf by successive adjacent exchanges. The leaf can then be retired. The remaining occupied count equals the target count on the smaller tree, which proves the claim by induction. Empty and full inventories each have a single configuration.

For sharp source $a$ and target $b$, let $d$ be their graph distance and $w>0$ the number of shortest walks. A walk of length $d$ from $a$ to $b$ cannot include a diagonal step. Its Hamiltonian weight is always $$-1/2$^d$. Every shorter off-diagonal contribution vanishes. The first nonzero term in the amplitude is consequently
\[
  \langle b|e^{-itH}|a\rangle=\frac{w}{d!}\left(\frac{it}{2}\right)^d+O(t^{d+1}).
\]
All amplitudes are analytic, and there are finitely many targets. For a fixed finite board and input there is a common sufficiently small positive interval on which all amplitudes are nonzero; at other times only a discrete exceptional set can remove configurations. This gives the capacity formula for the ordinary free evolution. The probabilities of its allowed words depend on the initial arrangement and the elapsed time. The receipt stores the exact shortest-walk count and distance separately for every sharp configuration. The six-cell graphs have 72 and 2,556 configurations, 185 and 12,950 undirected exchanges, and maximum distances eight and seventeen from the enrolled stocks. [Independent combinatorial and graph checks](test_inventory.py)

## What the external stations observed

Owners 18 and 19 on the twenty-cell board each compare the complete postflight sharp load word against their own references. Their routes are opened only after the coherent flight and reversed after the reading. The second station sees the sharp word projected by the first. For each available parity class an additional full two-head route execution checks the compact readout representation at all 72 internal ports and the vacant exterior sentinel. These simulated launches use the registered PCG64 stream; the mathematical support proof uses the generator itself. [Station routes and ternary calibration](../apparatus.py)

| Inventory | Native flight time | Parity words recorded | Computed probability of the allowed words absent from that sample |
| ---: | ---: | ---: | ---: |
| 1 | 0 | 1 | 0 |
| 1 | ¼ | 3 | $7.91\times10^{-9}$ |
| 1 | 2 | 4 | $1.11\times10^{-6}$ |
| 1 | 8 | 5 | $1.05\times10^{-5}$ |
| 2 | 0 | 1 | 0 |
| 2 | ¼ | 2 | $5.43\times10^{-6}$ |
| 2 | 2 | 4 | $1.05\times10^{-5}$ |
| 2 | 8 | 5 | $1.11\times10^{-5}$ |

The diagonal-only control retained the same native diagonal, cells, alphabet, input and times. Every launch returned its starting parity. The largest absolute frequency residual across the 1,024 registered probability bins was 0.000362 or less. The simultaneous 99% Hoeffding allowance was 0.009660. The missing masses in the table come from finite numerical evolution; the archive reports values below $10^{-12}$ as a threshold. The deposited full-configuration counts permit independent reconstruction of every parity bin.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.screen_count.free_inventory.measure
nice -n 15 .venv/bin/python -B -m pytest -q model/experiments/interval_echo_consensus/screen_count/free_inventory
~~~
