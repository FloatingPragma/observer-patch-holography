# Five product terms measure the fixed-volume entropy derivative

On the four-station survey, an observer can read the signed first variation
of head entropy at the calibrated maximum using five two-head Pauli products,
together with individual-head and parity readings. Every reader using four
such terms misses a variation of the same physical preparation. The signed
observable measures an entropy response in nats. Identifying a regional
stress that gives this response requires a separate physical measurement.

Use the 80-cell, six-unit survey with all timeouts included. Its ordinary
exchange calibration sets the positive gain `η`. The prepared state
`X*` uniquely maximizes `S + ηA` at the measured survey volume. The
source supplies a neighborhood of this state in the 31-dimensional
trace-one algebra of two head coins with two reversal parities. Thus on
every source variation with `dV = 0`, its entropy derivative is
`dS = −η dA`. Here both area and volume use the four stations’ projected
radar distances. The result concerns this choice of preparation, clock
and volume. [Source, calibration, reachability and stationarity](../PROOF.md)

## Which observables count

The parity-resolved records have two blocks, each containing two qubits.
A product term in block $b$ is a signed multiple of
$(u\cdot\sigma)\otimes(v\cdot\sigma)/2$ times the projector onto that parity.
Axes and weights may differ for every term and each parity. Each
individual-head reading and the parity contrast have arbitrary coefficients.
The count refers to terms in this factorization. The same means can be
formed from the nine laboratory Pauli-pair readings in each sector on
separately prepared copies, with the event weighted by its recorded parity.
Equivalently the complete timestamp tomography of the same source determines
each block expectation. [Native analyzers](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/README.md) · [Complete timestamp tomography](../reuse/README.md)

Let $U_{ij}=\sigma_i\otimes\sigma_j/2$, with $i,j\in\{x,y,z\}$.
These nine matrices are orthonormal in Hilbert–Schmidt norm and have zero
partial traces. In each parity set
$$
 A_b[i,j]=dA(U_{ij}^{b}),\qquad V_b[i,j]=dV(U_{ij}^{b}).
$$
A product term contributes a rank-one matrix to these coordinates.
Equality of two linear functionals on the fixed-volume tangent space
means they differ by a scalar multiple of $dV$. The multiplier is
common to both parities because the experiment constrains a single survey
volume. Consequently the number of product terms needed for exact
fixed-volume equality is
$$
 \min_{t\in\mathbb R}\bigl(\operatorname{rank}(A_e-t V_e)
                        +\operatorname{rank}(A_o-t V_o)\bigr).
$$
For the lower bound one may restrict to the 18 joint directions, where
all one-head and parity readings vanish. For the upper bound, take a
singular-value decomposition in each block and supply the one-head
coefficients and parity contrast of $-\eta(dA-t\,dV)$. The common trace
coefficient drops out. The recorded even-sector $ZZ$ volume derivative
exceeds $0.03$. The implicit function theorem and the source’s strictly
positive preparation weights therefore integrate every fixed-volume
tangent in this 18-dimensional subspace to an exactly volume-preserving
curve with both signs, keeping each block’s marginals and weight fixed.

## What the native range records determine

For every real $t$, each pencil $A_b-t V_b$ has second singular value
greater than $4\times10^{-8}$, and at least one has smallest singular value
greater than that same bound. Thus a four-term allocation either leaves one
block with at most one term, or gives both blocks at most two. For each
real price the distance from such a reader to the pair of pencils, in the
joint Frobenius norm, exceeds the stated bound. Minimizing this distance
over the price projects onto the fixed-volume tangents. Multiplication by
the independently calibrated gain $\eta>23.1042361$ gives a discrepancy
**greater than $9.24\times10^{-7}$ nats per unit Hilbert–Schmidt tangent**
on at least one signed preparation change for every four-term reader.

The even determinant has opposite signs at $t=2.8$ and $2.85$.
It vanishes between them, where the two ranks are exactly $(2,3)$.
The odd determinant changes sign between 0 and 1, giving the alternative
allocation $(3,2)$. **Exactly five terms suffice.** For direct use, the
signed reader freezes the first price at $2.822544426614$ and records
the five weights and axes, the six one-head coefficients in each block
($XI,YI,ZI,IX,IY,IZ$, each divided by 2), and the coefficient of
$(I_4,-I_4)/\sqrt8$. Its certified worst-case rounding, propagation,
and state-localization error on the fixed-volume tangent is less than
$3.106\times10^{-7}$ nats per unit Hilbert–Schmidt norm. The exact
five-term statement uses a root of the exact pencil. [Reader coefficients and interval certificate](results.json)

The certificate propagates the directed range-effect errors through the
actual tetrahedral area and volume at the source enclosed by the upstream
$2\times10^{-7}$ state bound. It encloses the range gradients by
integrating their interval Hessians from the rational center. For a
$3\times3$ pencil matrix $M$, every two-by-two minor bounds
$\sigma_1\sigma_2$ from below, and
$\sigma_3\ge |\det M|/\|\operatorname{adj} M\|_F$.
The two compact charts $t=4s$ and $t=4/s$, $|s|\le1$, cover all
prices. In the reciprocal chart the certificate uses $sA_b-4V_b$. Division by $|s|\le1$ preserves the lower bound. Directed interval
arithmetic verifies the bounds on 37 and 4 subintervals. Assigning the
same range panels to both parities introduces a common singular price
and the verifier rejects the control. [Certificate program](measure.py) · [Independent range replay and control](test_measure.py)

The rank result supplies an observable target for the complete linear
response of this finite source. Physical null balance further requires an independently measured material
flux and selection of a common source on shrinking proper balls.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/correlation_rank
```
