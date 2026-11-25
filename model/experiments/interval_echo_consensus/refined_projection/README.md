# Information projection for two marked faces

A prepared Interval Echo apparatus puts two distinguishable marks in descendants of one registered face. Its output law is the unique relative-entropy minimizer given the measured probabilities that the marks occupy the same ancestor at each available resolution. Forgetting the last child address sends every fine optimizer to the coarse optimizer exactly. Both directions of feasible-set transport hold. The proof works at every finite depth of the four-child construction. The complete native executions cover the 20 → 80 → 320 cell boards.

## Physical records and inputs

The preparation fixes base face 0. Each parent has three corner children and a central child. The parent of cell $a$ is $\lfloor a/4\rfloor$. The receipt verifies the existing cell-dual certificate against both pairs of deposited boards, including connectivity of each four-child region and the actual crossings of every parent interface. Iterating this subdivision keeps every descendant region connected. Connectedness suffices for the marker transport.

The source has a prepared sharp double, a material comparison bit and two marker emitters. Its five local half-contact readings per decision give one fair selector and two independent two-bit child addresses. Each conditional bit probability is $1/2$ on either preceding sharp state. The source takes all five readings even when a branch discards the second address. While the two addresses share a parent, selector one sends both marks to the first reported child. On selector zero the two reported choices stand. After separation both choices stand on either selector reading. The controller writes the updated agreement bit. No target agreement probability enters its configuration.

Native contact swaps put the marks on rails 0 and 6 of their designated cells. Each local detector compares its rail to a vacant port 10 by echo read and undo. A joint appointment collects the two complete sets of flags. At an intermediate depth the marks visit the first leaf of each selected ancestor on the physical board. The following decision consumes the observed ancestor addresses and the retained comparison bit. Undoing the transport restores every other rail. All reported probabilities refer to these flag records. The source controller, mark routing, common root and grouping by ancestor are design inputs. [Sharp-register implementation](../../../candidate_laws/interval_echo/clock.py) · [Geometric parent certificate](../../mobile_memory/refinement.py)

The independent reference uses the same five-bit source, marked rails and observations with both child addresses taken fresh throughout. Its ordered-pair table $U_L$ is uniform with mass $16^{-L}$. It pushes forward to $U_{L-1}$. This provides a counting reference for the declared address records. The reference on general repair histories is a separate premise.

## Feasible families and their optimizer

Let $E_j$ mean that the two marked faces have the same ancestor at depth $j$, with $E_0$ certain. Once $E_j$ fails, its descendants cannot agree. Conditional on $E_{j-1}$, the source copies with probability $1/2$, and the two fresh uniform child choices agree with probability $1/4$. Hence $s=\Pr(E_j\mid E_{j-1})=5/8$ and $\Pr(E_j)=s^j$. This is the observable supplied by the apparatus.

At depth $L$ take **every** probability table on ordered leaf pairs satisfying $\Pr(E_j)=s^j$ for $1\le j\le L$. Competitors can have arbitrary individual face frequencies, histories and within-stratum correlations. For $0\le k<L$, the stratum of exactly $k$ common child digits has mass $s^k(1-s)$ and $4^k\cdot12\cdot16^{L-k-1}$ outcomes. The all-equal stratum has mass $s^L$ and $4^L$ outcomes. Within each stratum the common-parent controller treats all leaf pairs equally. Call its resulting table $P_L$. The log ratio $\log(P_L/U_L)$ is constant within each stratum. For any feasible $Q$, summing this log ratio under $Q$ therefore gives the same value as summing it under $P_L$. It follows exactly that

```math
D(Q\Vert U_L)=D(P_L\Vert U_L)+D(Q\Vert P_L).
```

Strict positivity of $P_L$ makes equality with the minimum equivalent to $Q=P_L$. Each coarser table is a deterministic grouping of the same marks, so readers agree on all shared records. A jointly separating menu includes the full pair and reports its projections to each resolution $1,\ldots,L$. For the positive reader weights, stop at the first fresh fair one, with resolution $L$ receiving the residual probability. Thus $w_{L,j}=2^{-j}$ for $j<L$ and $w_{L,L}=2^{1-L}$. The two finest weights add to the preceding terminal weight. The native table minimizes every term of the weighted sum of marginal relative entropies. The full-resolution term with positive weight makes the minimizer unique.

## Refinement of arbitrary feasible competitors

Let $C$ discard the final address digit of each mark. Every finer feasible table pushes to a coarser feasible table because all older events are preserved. For a coarse feasible $q$, lift it by attaching the native conditional child table $K$ to each ordered parent pair. For equal parents this table gives $5/32$ to each equal child choice and $1/32$ to each unequal choice. For separated parents each of the 16 choices has mass $1/16$. The lift $qK$ has $C(qK)=q$ and $\Pr_{qK}(E_{L+1})=s\Pr_q(E_L)=s^{L+1}$. This proves surjectivity on the full feasible family, including asymmetric competitors. Taking $q=P_L$ gives $C P_{L+1}=P_L$.

The cost of refinement has an exact expression as well. Put $I=s\log(5/2)+(1-s)\log(1/2)>0$. For every fine feasible $Q$ and $q=CQ$,

```math
D(Q\Vert U_{L+1}) = D(q\Vert U_L) + s^L I + D(Q\Vert qK).
```

Indeed, the correction $\log(16K)$ is $\log(5/2)$ for equal children of equal parents, $\log(1/2)$ for unequal children of equal parents, and zero for separated parents. Their masses are fixed by $E_L,E_{L+1}$. This gives the identity by expanding the three divergences on the support of $Q$. The lift uniquely minimizes the fine cost at fixed $q$. The coarse addresses together with the single finer comparison bit are a sufficient statistic for the fine optimizer versus its reference, since the two conditional laws agree given that statistic. Coarse addresses alone lose the positive, explicitly measured cost $s^L I$.

A changed controller permits copying after the marks have separated. It preserves both uniform one-marker marginals and *every* registered ancestor-agreement probability. Its excess objective at full resolution equals

```math
I\sum_{j=0}^{L-1}(1-s^j).
```

At step $j+1$, previously separated parents have mass $1-s^j$. Their correlated child table costs $I$ more than the independent one under the native optimizer. All other conditional tables coincide. The excess is strictly positive from depth two onwards. Thus agreement statistics alone fail to certify that an arbitrary preparation attains the projection.

## Executable evidence

| Sphere cells | Complete five-draw histories per menu | Distinct physical address pairs | Measured ancestor agreement | Changed-wiring excess, nats |
| ---: | ---: | ---: | ---: | ---: |
| 80 | 32 | 16 | $5/8$ | 0 |
| 320 | 1,024 | 256 | $25/64$ | 0.117281818017 |

The three menus use the same candidate contact and routing. The branch enumeration exactly weights all sharp-instrument outcomes, including unused draws, and checks the native half-contact Born probabilities against the two-input operator reduction. The level-two panel executes 10,240 half-contacts per menu. The regression also checks the next abstract address-table pushforward independently of any higher-board execution. The theorem covers this prepared marked-face observer family. Unrestricted covers, unconstrained free Interval Echo evolution, Einstein null balance and gravitational stress require their own source maps. [Source and verifier](measure.py) · [Complete receipt](results.json) · [Regression](test_projection.py)

From the repository root, run

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m model.experiments.interval_echo_consensus.refined_projection.measure
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/refined_projection
```
