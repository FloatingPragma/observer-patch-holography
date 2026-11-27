# Three radar deadlines on the same retained source

At the exact four-station entropy equilibrium, the 12 original mean ranges and both one-head marginals leave six independent directions of two-head pairing undetermined. Varying these pairings controls **four earlier readings independently**: the area and volume reconstructed at deadlines two and four. The original radar closes at time six. During the previously established signed collar intervention the source can keep both shorter areas and both shorter volumes exactly constant, while its entropy increases and its late area contracts at the reference. These are nested time cutoffs at the same stations. Interpreting them as proper radii would need independent geometry. [Interval rank certificate](certificate.json) · [Original signed interventions and retained controllers](../README.md)

## What a second observer reads

The monitor checks every half exposure. It retains every missed launch or return. A recipient with deadline $c$ assigns one-way score $\min(n,2c)/4$ to a completed return at tick $n$, and $c/2$ to a timeout. The early recipient can read the prefix without intervening in the continuing cycle or changing the terminal feedback record. The observer applies the same clock-gain correction and tetrahedral reconstruction separately to each deadline. The three areas are areas of three measured tetrahedra built from the same four physical receivers. [Integer propagation and early scores](kernels.py)

Two archived endpoint banks already have all the required event times. They compare the commissioned positive and negative twirl steps, $10^{-5}$ on each side. Scoring the same trials at each deadline gives the following secants, in squared and cubed radar units per unit change in twirl probability.

| Deadline | Area, simultaneous 99% interval | Volume, simultaneous 99% interval |
|---:|---:|---:|
| 2 | +0.001609 to +0.001699 | +0.000171 to +0.000183 |
| 4 | +0.015381 to +0.015902 | +0.004037 to +0.004143 |
| 6 | −0.019033 to −0.017758 | −0.000183 to +0.000191 |

These count-only intervals cover the registered digital sampling laws on the two prepared endpoints. All 72 clock means enter one Hoeffding union bound, including their reuse across deadlines. The sample contains $2^{66}$ synthetic launches per route, endpoint and common observation menu. No extra trials were generated. The rank theorem below uses separately bounded native probabilities at the exact reference. [Counts reader](records.py) · [Frozen source banks](../observer_banks/records.json)

## What the native clock effects determine

Perturb the three-by-three Pauli pairing coefficients in each of two reversal blocks. These 18 real coordinates preserve the parity weights and both per-block local coin marginals. All are available locally at the faithful reference. Let $R_c$ be the linear part of the twelve directed means on this slice. The certified inverse of a twelve-column minor of $R_6$ fixes an exact six-dimensional kernel. For each of its six free coordinates take the unit step in that coordinate and the unique compensating change in the pivot coordinates. These physical tangent columns form $Z$ and $R_6 Z=0$ exactly. The interval calculation constructs this kernel for the actual native effects, including their uncertainty. [Source algebra and preparability](../../../PROOF.md)

For each of the two earlier deadlines, compose $R_c Z$ with the exact first derivatives of its area and volume. The combined four-by-six matrix $D$ has a four-column minor with a right-inverse residual norm below **0.00374** after row scaling. This proves differential rank four at the actual selected equilibrium. In particular these readings do not descend to four functions of the original twelve ranges alone.

For the error calculation, the first-hit recurrence for the shorter deadlines uses subsets of the original success terms and smaller absolute score weights. The cap-six operator error $7.82\times10^{-16}$ bounds each native deficit. Hermitian projection and $2^{-56}$ rounding give total operator error below $9\times10^{-16}$. The reference's Hilbert–Schmidt radius bounds each range using its nonconstant deficit. All geometric gradients are enclosed on those ranges. For the projection step, write $P$ for the twelve pivot columns of $R_6$. Combining the earlier row with $P^{-1}$ before bounding perturbations retains cancellations between the surveys. The reader calculates the center with rational arithmetic and verifies the final minor by outward intervals. [Projection proof and bounds](certificate.py) · [Frozen native deficits](deficits.npz)

## A held earlier area at fixed volume

Let $X(s)$ be the existing signed collar-erasure path, compensated by the right-head rotation so that $V_6$ is fixed. Seek $Y(s)=X(s)+Zz(s)$, with $z(0)=0$. The implicit-function theorem on the invertible four-column minor supplies four components of $z(s)$ keeping $(A_2,V_2,A_4,V_4)$ equal to their reference values for all sufficiently small $s\ge0$. The twelve cap-six ranges agree **pointwise** for $Y(s)$ and $X(s)$ because they are affine functions of the state. Thus these paths also have identical $A_6$ and $V_6$ throughout. The added term preserves the original local marginals, parity weights and parked energy. Faithfulness is open, and the earlier strict-majorization proof supplies a separately programmed retaining controller for each nearby target. The path uses calibrated choices of source populations and feedback.

At the reference, stationarity of $S+\eta A_6$ on the $V_6$ level gives $dS(Zz')=0$. Hence the entropy derivative of $Y$ is the independently certified collar slope between 0.23107940 and 0.23113300 nats per unit, and its late area derivative is $-S'/\eta<0$. Both earlier area derivatives vanish. Conversely, the same rank result permits arbitrary sufficiently small smooth changes of the four earlier readings with the original twelve means held fixed. The result establishes local existence. It gives no lower bound on an experimentally accessible path length. The twirl interval table concerns a separate commissioned path and tests the practical effect of assigning the different deadlines to real model records.

Run from the repository root with the shared scientific dependencies.

```sh
export PYTHONPATH=code:. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=matched_einstein_limit.source_implication.full_algebra.reuse.retained_equilibrium.nested_ranges
nice -n 15 python -m "$p.kernels"
nice -n 15 python -m "$p.certificate"
nice -n 15 python -m "$p.records"
nice -n 15 python -m pytest -q --pyargs "$p"
```
