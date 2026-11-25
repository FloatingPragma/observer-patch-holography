# Counting native marker histories

One round of two fair Interval Echo reads produces three marker positions with probabilities $(1/2,1/4,1/4)$. Two of the four equally likely source tapes leave the marker at its starting port. Counting those tapes gives a reference for longer histories and for readers who see only part of a history. Conditioning on an observed first departure realizes the unique relative-entropy minimum at every tested depth. Replacing the departure constraint by a mean-position constraint gives an irrational optimum for the same apparatus. A bounded fair-read program can approximate that second law, and an adaptive program samples it exactly with three reads on average.

## Apparatus and counting measure

Cell zero contains a sharp double at marker port 0 and another at source port 2. Other loads start at zero. The source receives a half contact on ports 2 and 3 and resolves port 3 against vacancy 4 by echo, read and undo. Both sharp input orders give a fair outcome from the registered 99-state seam Hamiltonian. In each round, the first bit conditionally exchanges marker ports 0 and 1, and a fresh bit conditionally exchanges ports 1 and 5. Three vacancy comparisons at the appointed time publish the one-hot marker record. The labels $0,1,2$ denote ports $0,1,5$. The apparatus supplies the sharp preparation, isolated contacts, appointment schedule, classical feedback and Born read. All its operations lie within one owner on each board. [Native registers](../../../candidate_laws/interval_echo/clock.py) · [Experiment](measure.py)

For a deterministic record map $f$ consuming at most $N$ native fair reads, extend every stopped branch with its unused bits. A record $h$ has probability $R(h)=|f^{-1}(h)|/2^N$. A further unused bit doubles the numerator and denominator. Coarsening by any $g$ gives $g_*R=(g\circ f)_*U$ for the uniform complete-tape law $U$. PR-05's uniform step is a bit choice here. The marker law counts all bit choices that yield the same observation. In label order its one-round transition matrix has rows $(1/2,1/4,1/4)$, $(1/2,1/4,1/4)$ and $(0,1/2,1/2)$. [Finite pushforward, fibre counts and padding in Lean](../../../../Lean/InformationProjection/NativePathReference.lean)

## Projection and refinement

Write $A$ for the first marker departure. Retain the full record of every run, including the probability $1/2$ of rejection. The admitted law is $P=R(\cdot\mid A)$. The observer selects a prefix length by stopping at the first fresh fair one, with the longest length taking the residual probability. At depth $L$ this gives weights $w_j=2^{-j}$ for $j<L$ and $w_L=2^{1-L}$. Each prefix contains the first marker read, so $A$ is visible to every reader.

For any joint competitor $Q(A)=1$, the ratio $P_j/R_j$ on the admitted prefix support equals $1/R(A)=2$. Directly substituting in relative entropy gives

```math
\sum_j w_j D(Q_j\Vert R_j)
  = \log 2 + \sum_j w_j D(Q_j\Vert P_j).
```

The full-prefix term has positive weight, so equality at the minimum requires $Q=P$. The argument also applies to the complete source-bit prefixes when the reader retains those records. Deleting the final round pushes the fine reference and optimizer to their coarse laws. Every feasible coarse law $q$ has a feasible fine lift $q(h)R(h,x)/R(h)$ on reachable prefixes. Two adjacent terminal reader weights add to the coarse terminal weight. Refinement here extends the observation time on the same local program. Deploying it on larger boards checks the local incidence in each deployment. The state-selection dependency is conditional maximum randomness relative to the declared reference. [Premise registry](../../../../claims/axiom_registry.yaml)

| Board cells | Rounds | Source tapes | Reachable marker histories | Rank of separate time marginals | Rank with full marker history |
|---:|---:|---:|---:|---:|---:|
| 20 | 1 | 4 | 3 | 3 | 3 |
| 80 | 2 | 16 | 8 | 5 | 8 |
| 320 | 3 | 64 | 21 | 7 | 21 |

At depth three the marker history leaves 43 directions of the complete source distribution invisible. A feasible controller that keeps the marker stationary after the first round incurs 0.779791 nats of excess marker cost. A second control uses a $\pi/3$ source contact. A positive read triggers a native swap to restore the input orientation, giving independent probabilities $3/4,1/4$ on every draw. Its conditional marker law incurs 0.279125 nats of excess at depth three, versus 0.327030 nats for its complete source record. [All branch records and costs](results.json)

## Exact sampling at a fixed mean

A separate preparation sends the marker to the sum of two measured fair bits. Its mean label is one. For the original three-label reference $R=(1/2,1/4,1/4)$, the law minimizing $D(Q\Vert R)$ at that same mean is

```math
P=(a,1-2a,a),\qquad a=(4-\sqrt 2)/7.
```

The ratio $P_i/R_i$ is proportional to $(\sqrt 2)^i$. Thus every feasible $Q$ averages $\log(P_i/R_i)$ to the same value as $P$, which gives $D(Q\Vert R)-D(P\Vert R)=D(Q\Vert P)$. The two-bit preparation has excess 0.129446 nats. A bounded fair source with $k>0$ admitted tapes can produce only masses $m/k$. The integer $7m^2-8mk+2k^2$ is nonzero because both polynomial roots are irrational. On $0\le m/k\le 1$, factoring about $a$ and using $|7(m/k+a)-8|\le8$ gives $|m/k-a|\ge 1/(8k^2)$. The obstruction includes finite postselection.

An adaptive reader uses the same physical fair source and a classical integer controller. The bits specify a nested dyadic interval for a uniform number $U$. Output zero for $U<a$, two for $U>1-a$, and one between the thresholds. The controller reflects its interval into $[0,1/2]$ and tests the signs of $7x^2-8x+2$ at its two endpoints using integer arithmetic. The polynomial is strictly decreasing there. A cell that straddles the root needs another bit. For every $n\ge1$ exactly two of the $2^n$ depth-$n$ cells straddle a threshold. The algorithm terminates almost surely and $\mathbb E[T]=1+\sum_{n\ge1}2^{1-n}=3$. This exact law assumes an indefinitely available fair source and a controller that can extend its integer prefix without a fixed limit.

The finite replay stops at 12 reads, labels both unresolved branches, and publishes the middle marker on those branches. Its law is $(1513,1070,1513)/4096$, its mean is exactly one, and its total variation distance from the optimum is bounded by the unresolved mass $1/2048$. The expected number of reads is $3071/1024$. The receipt recomputes every stopping prefix through the native fair register and every terminal marker through the physical port comparisons. [Stopping tree and unresolved mass](results.json)

## Coherent inputs

Sharp permutations and diagonal data reads, even when controlled by the independent fair register, yield diagonal effects on the marker basis. A complete diagonal read on $d$ sharp labels leaves a $d(d-1)$-dimensional real off-diagonal kernel on trace-zero Hermitian perturbations. On the actual two-rail seam, $(|20\rangle\pm i|02\rangle)/\sqrt2$ have identical sharp probabilities. A native half contact on the **data** followed by the same sharp read produces opposite certain outcomes. Removing the contrast term preserves this two-rail test. The path experiment therefore certifies classical record separation. General coherent state selection requires corresponding phase-sensitive instruments and their joint constraints. [Native generator comparison](measure.py)

From the repository root run `python -m model.experiments.interval_echo_consensus.native_path_reference.measure` to check the stored receipt, and `python -m pytest -q model/experiments/interval_echo_consensus/native_path_reference` for the five focused tests. From `Lean/`, run `lake env lean -j1 InformationProjection/NativePathReference.lean` to check the pushforward and sampling theorems.
