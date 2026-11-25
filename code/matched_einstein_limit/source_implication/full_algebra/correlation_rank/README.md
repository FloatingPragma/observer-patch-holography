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

## A material measurement on the same preparation

At the exact constrained maximum there is a source curve on which the
survey volume and all twelve registered regional translation means are
constant, while the entropy increases. Both one-head marginals in each
parity and the parity weights stay fixed. Per unit Hilbert–Schmidt tangent
length the entropy slope lies in **[0.00159319, 0.00159487] nats** and the
area slope in **[−0.000069030, −0.000068956]** survey area units.

The regional experiment launches a vacancy or excess through the material
chart and measures the translation phase, the outgoing and incoming cut
contributions, and the force due to the two stationary marks after four
ticks. Their exact momentum balance supplies the summed boundary flux and
the total momentum rate as well. The material operators come from the
regional Hamiltonian and its partial translation. The radar reconstruction
uses separate copies of the same two-head preparation and a six-unit clock.
The translation operators have no entropy or area coefficients.
[Regional experiment and balance law](../../directional_stress/translation/README.md)

For each charge, compress the four independent material operators onto
the joint Pauli directions in each parity. Stack the eight rows with the
measured volume derivative to form a matrix $B$ with 18 columns. Their
rank is nine on the exact source. With the published rational rescalings
of the rows, an approximate right inverse $K$ gives
$\|I-BK\|_F<0.0000525$. Choose a central tangent $d$ in the direction of
decreasing area after projecting out those nine rows. For the true
derivatives,
$$
 D=d-K(BK)^{-1}Bd,\qquad BD=0,\qquad \|D-d\|_2<0.0000716 .
$$
Subtracting a fixed linear combination of the constraint rows from the
area derivative preserves its value on $D$. Interval evaluation of the
result on $D/\|D\|_2$ gives the signed bounds above with the ordinary gain.
All 18 basis directions have zero marginals, and the calibrated physical
preparation contains a neighborhood in that slice. The submersion theorem
therefore integrates $D/\|D\|_2$ to a local curve of positive preparations
with exactly equal volume and exactly equal material means on both arms.
The constrained local level surface has dimension nine.
[Native effects and exact-kernel construction](translation.py) · [Enclosures and preparation direction](translation.json)

Any scalar computed solely from those twelve material means is constant
on this curve. It cannot give the measured entropy slope, even with freely
chosen coefficients or a nonlinear response function. Spatially resolved
material records can contain information lost in these regional summaries.
The physical null-balance premise further requires an independently read
stress and a common family of shrinking proper balls.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.correlation_rank.translation --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/correlation_rank
```
