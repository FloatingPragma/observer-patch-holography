# Cooling stored prefix records and free doubled loads

A longer observation produces a growing native prefix alphabet with one uniform cooling bound in attempt-cost units. Applying the same ground-projection criterion to the physical flight generator has a different outcome on growing connected degree-six federations. Each calculation uses the finite spectrum or complete record probabilities of its specified experiment. [Executable comparison](measure.py) · [Results](results.json) · [Independent bounds checks](test_cooling.py)

## A cofinal family of observed word lengths

An observer prepares a sharp double on either rail of the native fair source in cell zero. It applies the registered half-contact and sharp read on each attempt and stops at the first zero or after $m\ge2$ attempts. The half-contact has probability $1/2$ for either output from both sharp starting rails. The read leaves another sharp rail preparation. The conditional probability of a word of length $k$ is therefore $2^{-k}$ on a continuing source. The observer retains the complete prefix and its cost $E=k\log2$. Histories refine by truncation to the shorter cap. These maps compose, preserve the unique one-bit ground record $0$ and cannot raise its cost or that of any other record. The identity on the finite depths $m\ge2$ is cofinal and order-reflecting in the same unbounded observation-depth order. Each step increases the number of public labels from $m+1$ to $m+2$. Depth two gives exactly $0,10,11$ with the finite record-cost model's calibration. [Registered fair source](../../../candidate_laws/interval_echo/clock.py) · [Three-state finite theorem](../../../../Lean/Thermodynamics/NativeRecordFourLaws.lean)

Write $r=2^{-\beta}$. There is one word at every length $1$ through $m-1$ and two at length $m$. After dividing the partition function by the ground weight, the excited-to-ground ratio is $R_m=\sum_{j=1}^{m-2}r^j+2r^{m-1}$. For $\beta\ge1$, $R_{m+1}-R_m=(2r-1)r^{m-1}\le0$. Hence the exact worst-case off-ground mass is $R_2/(1+R_2)=2r/(1+2r)$. With $j=k-1$, the thermal entropy is $\log(1+R_m)+\beta\log(2)\mathbb E[j]$. Each excited $j$ has multiplicity at most two, so
$$
S_m\le\log(1+2r)+2\beta\log(2)\frac{r}{(1-r)^2}.
$$
Both bounds tend to zero uniformly in depth. At $\beta=1$ the fair source law pushes forward exactly under truncation. For other multipliers the energy-compatible truncation maps need not push the Gibbs law to the Gibbs law. For example, the relative weight of the two children of the old terminal all-ones word multiplies by $2r$ when the cap grows.

At positive integer $\beta=b$, agreement of $b$ independently read prefixes implements the conditional Gibbs distribution. Each proposal has success probability $Z_{m,b}=\sum_h 2^{-b|h|}\ge2^{-b}$. A capped controller returns its old record on a failed proposal and retains a timeout. After $M$ proposals the residual old-record weight is $(1-Z_{m,b})^M$. Unbounded retry terminates almost surely. The expected native cost per accepted sample is $2b(1-2^{-m})/Z_{m,b}\le2b\,2^b$, uniformly in depth. This follows by conditioning the retry cost on the first proposal, using the per-prefix expected attempt count $\sum_{k=1}^m2^{-(k-1)}$. The cofinal sequence of realized integer multipliers suffices to exhibit the cooling limit. The analytic envelope also applies to positive real multipliers $\beta\ge1$. The observer supplies the stopping rule, memory and replica protocol. The measured cost is the number of recorded half-contacts times $\log2$.

The code executes every underlying tape for four depth/replica settings from both sharp starting rails. At depth three, three replicas agree on $0,10,110,111$ on respectively 64, 8, 1 and 1 of 512 tapes. All 438 disagreements remain in the outcome table. [Exact outcome counts](results.json)

## The native flight ground mode

With one double and all seams blank, the physical generator on the port label is $24\kappa I+L/2$, including degree-reference contacts. Here $L$ is the connected board's unweighted Laplacian. A scalar energy shift cancels from all Gibbs probabilities. The unique ground mode is constant. Every shifted eigenvalue is in $[0,6]$ because the maximum degree is six. On $N$ ports this gives, for every fixed $\beta\ge0$,
$$
p_0\le[1+(N-1)e^{-6\beta}]^{-1}.
$$
In fact there are exactly three external ports per twelve-port cell, giving $\operatorname{Tr}(L/2)=21N/8$. Jensen's inequality on the $N-1$ excited levels gives the stronger denominator $1+(N-1)\exp[-21\beta N/(8(N-1))]$. Thus along any unbounded connected federation family with this degree structure, the supremum of off-ground mass at every finite $\beta$ is one. Each finite connected member cools to its ground mode as $\beta\to\infty$. Since the largest eigenvalue is $p_0$, the total entropy also satisfies $S\ge-\log p_0\ge\log[1+(N-1)e^{-6\beta}]$. The ground state has zero entropy. Ground-projection mass and this total entropy concern the whole one-double sector. Energy density, individual local marginals and physically rescaled generators have separate criteria. [Source reduction](../FLIGHT_CONFLUENCE.md)

The deposited 20-, 80- and 320-cell board generators pass an exact sparse-matrix identity and connectivity check. Direct diagonalization at $\beta=1$ gives ground probabilities $0.026310727686$ and $0.006578142421$ on the first two boards. The unbounded claim is a graph theorem under the stated growth conditions. These finite checks do not identify a macroscopic physical regulator. [Flight spectrum evaluations](results.json)

Run the replay and regression from the repository root.

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m model.experiments.interval_echo_consensus.refinement_cooling.measure
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/refinement_cooling
