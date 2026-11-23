# When finite repair makes regional entropy stationary

Two overlapping observers share a fair bit and each makes a private error on one quarter of preparations. Independent private repair coins select a unique maximum of the three-bit regional entropy, with zero first variation in **all four feasible joint directions**. A different preparation gives the two owners identical record tables and preserves the complete instantaneous cut, yet its regional entropy changes at **−log 5 nats per unit change in the joint error probability**. An independent seam pointer reads this change in the material correlations. [Measurement design](PROTOCOL.md) · [Finite counts](results.json) · [Coarse regions](coarse.json)

Each owner spans cell zero and one of its neighbours. The overlap contains the sharp digit $X$, and the private cells retain $Y$ and $Z$. The two owner distributions in $(X,Y)$ and $(X,Z)$ have entries $(3,1,1,3)/8$ in binary order. These are the unique accessible minima of relative entropy to the registered uniform reference with the supplied fair-source and quarter-mismatch constraints. The independent preparation reproduces the original projection's joint outcome probabilities after summing over the two extra fair decisions used by this apparatus. The additional sharp scratch rails are returned to fixed values before the terminal survey. The material region thus has eight orthogonal terminal words. Source settings and discarded feedback history lie outside the terminal record algebra. [Registered projection](../../../candidate_laws/interval_echo/projection.py)

Write $E_L=X\mathbin{\mathrm{xor}}Y$ and $E_R=X\mathbin{\mathrm{xor}}Z$. The four constraints on the eight joint probabilities specify normalization, $P(X=1)=1/2$, and $P(E_L=1)=P(E_R=1)=1/4$. At the product of the two private conditional distributions, the joint modular generator is
\[
-\log p_*=\log(32/9)+(\log3)(E_L+E_R).
\]
Its expectation is fixed on the entire feasible family. For any feasible table $p$, $S(p_*)-S(p)=D(p\Vert p_*)$, which gives uniqueness and stationarity. On this classical record algebra the shared-bit sector split is $S(X,Y,Z)=\log2+S(Y,Z\mid X)$. The variation therefore lies in the conditional regional entropy. This split and the material entropy refer to the specified accessible algebra.

The comparison uses two sharp bits $A,B$ produced by fair native contacts. One controller gives both owners the same error $AB$. The other gives the left owner $AB$ and the right owner $A(1-B)$. Varying their mixture physically changes $\theta=P(E_L=E_R=1)$. The four error probabilities are $(1/2+\theta,1/4-\theta,1/4-\theta,\theta)$, independent of $X$. On their faithful interior the finite first law gives
\[
{dS\over d\theta}=\log{(1/4-\theta)^2\over(1/2+\theta)\theta}.
\]
At $\theta=1/16$ this derivative vanishes. At $\theta=1/8$ it equals $-\log5$. Both owners' complete record tables stay at their constrained entropy minima on the entire variation. Each of the three cells contains six units of material load. The full instantaneous distribution of both endpoints on all five outgoing contacts is fixed. The echo pointer is undone after each sharp read. Every function of that terminal cut, including a proposed area in any normalization, has zero variation.

| 320-cell survey | Derivative at the reference | Exact central difference | Central difference from fresh counts |
| --- | ---: | ---: | ---: |
| Regional entropy, independent repair | 0 | +0.00843 | +0.00242 ± 0.01485 |
| Regional entropy, correlated repair | −1.60944 | −1.61195 | −1.58417 ± 0.01722 |
| Seam mismatch, correlated repair | −2 | −2 | −1.99233 ± 0.00954 |

The two-sided step is $1/64$ and each setting uses 4,194,304 simulated preparations. Errors are one propagated standard error. The joint-word and seam surveys use independent ensembles. The seam instrument brings the two private rails to adjacent positions by native swaps, reads the echo pointer and reverses the routes. The probability of a sharp mismatch is $1/2-2\theta$. Equal sharp rails contribute zero expected native contact energy, while unequal double rails contribute 32.5 units. The corresponding slope of the addressed comparator energy is −65 native units per unit $\theta$.

## Connected coarse regions

On the 20-, 80- and 320-cell boards, connected regions of 18, 66 and 279 cells contain 6, 22 and 93 disjoint three-cell repair preparations. Each group uses private output rails only in its interior. For every group the code tests all eight retained words against its routes, per-cell load and outer cut. Relative to the supplied independent preparations across groups, the entropy of the complete regional word adds. Its derivative at independence vanishes. At $\theta=1/8$ it is respectively −9.65663, −35.40763 and −149.67773 nats per unit $\theta$, or $-(\log5)/3$ per enclosed cell. The mean seam disagreement retains slope −2. All joins between groups are parked during the guarded repair service. [Connected regions and terminal-word checks](coarse.py)

This gives a finite stationarity mechanism on the materially prepared product reference and a counterexample to transferring entropy stationarity from the two-owner record algebra alone. Proper spacetime volume, its independently responding boundary area and a calibrated vacuum on shrinking geodesic balls are additional geometric inputs to the fixed-volume Einstein premise. [Einstein small-ball identity](../../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

Replay from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_einstein_station.repair_stationarity.measure --output /tmp/repair-stationarity.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_einstein_station.repair_stationarity.coarse --output /tmp/repair-coarse.json
.venv/bin/python -m pytest -q model/experiments/interval_echo_einstein_station/repair_stationarity
~~~
