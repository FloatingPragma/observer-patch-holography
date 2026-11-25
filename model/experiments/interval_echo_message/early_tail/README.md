# How early can the neighbouring port distinguish the two releases?

On both existing passive message stations, the receiving gateway has different click probabilities for the two bits at every free exposure between zero and 1/36, excluding zero. An observer who guesses bit zero on a gateway click and bit one on every other outcome has success probability strictly above one half throughout that interval. At zero both inputs give the quiet outcome with certainty. The infimum of exposures with nonzero discrimination is therefore zero for these two stations. The appointments use the original free-exposure parameter. The size of the advantage and the time recorded by a working clock are separate measurements. [Preparation and complete receiver](../README.md)

Let $a$ and $b$ be the sender ports for bit zero and bit one, and $r$ the receiving gateway. In the one-load basis of the selected blank seam, the contact is $(8+1/2)I-X/2$. Its native π pulse swaps the launch onto the selected sender port. Subtracting a scalar from the one-load Hamiltonian gives $K=L/2-3I$, where $L$ is the unweighted graph Laplacian. Since every vertex has degree at most six, $0\le L/2\le 6I$ and $\|K\|\le3$. The registered boards have $K_{ra}=-1/2$ and $K_{rb}=0$. For either launch the zero-order amplitude at $r$ vanishes. The integral remainder in $e^{-itK}=I-itK+R$ obeys $\|R\|\le9t^2/2$ by unitarity. Thus for $0<t\le1/36$ the near amplitude has magnitude at least $t/2-9t^2/2$, and the far amplitude has magnitude at most $9t^2/2$. Squaring and subtracting gives
$$
P(r\mid 0)-P(r\mid 1)\ \ge\ t^2(1/4-9t/2)\ \ge\ t^2/8.
$$

At exposure 0.01 the click probability difference is 0.000024998 on both screens, giving this fixed decoder a success probability of 0.500012499 for equal priors. These are calculated probabilities. The bound concerns a single adjacent receiver in a finite continuous-time unitary channel. A detection threshold, an observer-clock deadline, and a finite-capacity relay can each give a different operational timing question.

From the repository root, run `.venv/bin/python -B -W error -m model.experiments.interval_echo_message.early_tail.check`. It checks the native sharp source and matrix entries on 20 and 320 cells, compares the complete receiver probabilities to an independent reconstruction from the incidence files, and tests the bound at three exposures. [Checker](check.py)

## How many launches buy a 52% chance of recovering one bit?

At free exposure 0.01, repeat the same unknown balanced bit using independent preparations. Give the receiver all twelve ports of one cell at the end of each launch. The table bounds the number of launches needed to reach the original protocol's 52% mean decoding accuracy. The lower limit applies even to a joint quantum measurement across all received preparations. The upper limit uses a detector that stores a single classical answer to the question “Did any launch put the load in my cell?” Distances here count port links from the nearer source to the cell. A receiver clock would have to be ranged separately.

| Board | Port distance | Receiver cell | Necessary launches | Sufficient launches |
| --- | ---: | ---: | ---: | ---: |
| 20 cells | 1 | 0 | 1,601 | 2,001 |
| 20 cells | 2 | 15 | 256,017,068 | 320,021,335 |
| 320 cells | 3 | 2 | 23,041,339,242,324 | 28,801,674,052,905 |

The bit ports and launch preparation come from the two registered stations. For each board and port distance one through six the checker takes the first cell with unequal source distances, and reports an empty selection when none qualifies. It calculates the probabilities for both sources and both flight laws. These large repetition counts follow from the probabilities. No such ensemble was drawn.

Write $p_b$ for the probability of any load arriving in that cell given bit $b$. After tracing out all other cells, the receiver holds a vacuum word of weight $1-p_b$ and a one-load block of weight $p_b$. For $N$ independent repetitions both receiver states contain the identical all-vacuum word with a common weight of at least $(1-\max_b p_b)^N$. Removing that common term bounds the trace distance by $1-(1-\max_b p_b)^N\le N\max_b p_b$. With equal priors the optimal success probability is one half plus half the trace distance. Thus a rigorous upper bound $u$ on both arrival probabilities implies that reaching 52% requires at least $\lceil 0.04/u\rceil$ launches.

For the converse, guess the nearer bit if any load arrives, and the farther bit otherwise. If $\ell$ bounds the nearer arrival probability below and $v$ bounds the farther one above, success is at least

$$\frac12\left(1+\frac{N\ell}{1+N\ell}-Nv\right).$$

This follows from $(1-\ell)^N\le(1+N\ell)^{-1}$ and the union bound. Each sufficient entry uses $N=\lceil 0.05/\ell\rceil$ and clears 52%. Halving the exchange rate during flight while holding the two released port states and reader fixed leaves even the optimal receiver below 52% at these same counts. Identical released states would give exactly 50%.

The checker verifies the candidate's centered generator against the integer port Laplacian, then evaluates its degree-22 Taylor polynomial with integer numerators and a common denominator. The integral remainder bounds the vector error by $\epsilon=(3|t|)^{23}/23!$. Projection onto any receiver cell changes the squared probability from the polynomial value by at most $2\epsilon+\epsilon^2$. The calculation assumes a terminal read, a fresh copy for each launch and a passive flight. Interrogations during flight, relays and interacting clock hardware generate different records. [Exact arithmetic and complete board output](launch_budget.py)

Run `.venv/bin/python -B -W error -m model.experiments.interval_echo_message.early_tail.launch_budget` from the repository root.
