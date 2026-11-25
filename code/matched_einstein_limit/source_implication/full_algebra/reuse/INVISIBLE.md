# Equal directed ranges and different entropy at the reused-radar equilibrium

At the volume-matched stationary state of the uniformly repeated 80-cell radar survey, a physically preparable perturbation increases the stored-head entropy while **all twelve directed mean ranges, both head marginals and the parity weight remain exactly fixed**. The ranges here are the registered six-unit means, each pooling the two enrolled probe charges with equal weights and counting timeouts. Thus their tetrahedral area and volume stay fixed as well. A parity-blind joint Pauli $Y\otimes Y$ read distinguishes states along this family. [Stationary source](../reuse_equilibrium.py) · [Certificate](invisible.py) · [Values and tangent coordinates](invisible.json)

## Finding a direction invisible to every range

In each parity block allow the nine orthonormal Hermitian products $\sigma_a\otimes\sigma_b/2$, with $a,b\in\{X,Y,Z\}$. Every product has zero trace and zero partial traces. On the two unnormalized parity blocks their 18 independent coefficients form a Euclidean coordinate space, with norm equal to the Hilbert–Schmidt norm of the state displacement. This choice preserves even the separate one-head coin marginals in each parity block. Let $R_{ij}=\operatorname{Tr}(E_i D_j)$ for the twelve complete directed mean-range effects. A vector in $\ker R$ preserves every range on a straight line, at finite displacement.

Each certified native effect is within $\epsilon$ in operator norm of a rational matrix. Since $\|D_j\|_1=2$, the rational range matrix $\bar R$ satisfies $\|R-\bar R\|_{\rm op}<30\epsilon$. The verifier uses a floating singular-value decomposition only to propose a rational right inverse $K$ and a signed entropy direction. It evaluates
$$
d=w-K\bar Rw,\qquad
b=\|I-\bar RK\|_F+30\epsilon\|K\|_F<10^{-6}
$$
with rational products and outward norm bounds. Consequently $RK$ is invertible, $R$ has rank 12, and its kernel has dimension six. The exact physical direction
$$
 D=d-K(RK)^{-1}Rd
$$
satisfies $RD=0$ and $\|D-d\|_2<4\times10^{-7}$. The bound follows by applying the inverse norm estimate $1/(1-b)$ to $\|Rd\|\le\|\bar R d\|+30\epsilon\|d\|$. The receipt gives the center $d$ in the basis $XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ$, with the even block first. The physical effects and this formula determine the correction.

## Resolving the entropy and a different joint coin read

The calculation uses the integer propagation of the 79 terminal records. If $\delta$ bounds the approximation error of each terminal Kraus operator, the direct sum of these operators differs from the exact Stinespring isometry by at most $\sqrt{79}\,\delta$. Tracing out the record and probe then costs at most $2\sqrt{79}\,\delta+79\delta^2$ in Hermitian induced trace norm. Applying this sharper estimate to the same dyadic channels gives a uniform channel error below $1.91\times10^{-7}$. The certified Hilbert–Schmidt contraction of 0.92 gives a conditional fixed-state error below $2.53\times10^{-6}$ in each parity. [Instrument, projection and Taylor error](README.md) · [Conditional contraction](../reuse_equilibrium.py)

All conditional spectra along the segment from the rational center to the exact fixed state are bounded below by 0.18. On this segment the resolvent formula for the matrix logarithm bounds its Hilbert–Schmidt change by the state error divided by 0.18. Bessel's inequality bounds the error in the nine orthonormal correlation coordinates by the same value. Combining both blocks costs a factor $\sqrt{2}$. The log parity weights in the unnormalized blocks contribute scalar identities and vanish on these traceless directions. An independent rational logarithm bound at the centers supplies the remaining numerical error. Pairing this score with $d$, including its certified correction to $D$, proves a strictly positive entropy derivative exceeding $10^{-5}$ nats per unit displacement.

The exact volume-matched parity weight lies between 0.134 and 0.1346. At every weight in this bracket the identity coefficient exceeds the sum of all absolute nonidentity Pauli coefficients by more than 0.20 in both native preparation tables (retained and dephased). The eigenstates of each product Pauli are mixtures of product Pauli eigenvectors in its two eigenspaces. These inequalities therefore put both tables in the interior of the native product-lottery cone. A displacement $tD$ changes their coefficient sums by at most $12|t|\|D\|_2$, so every $|t|\le0.01$ is preparable. Along the shorter segment from $t=0$ to $10^{-7}$ the physical eight-dimensional state has minimum eigenvalue above 0.024. The magnitude of the entropy Hessian there is at most $\|D\|_2^2/0.024$. The positive first derivative therefore gives a **finite** entropy increase greater than $7\times10^{-13}$ nats at $t=10^{-7}$. This interval certifies a finite-model difference. Resolving such a small entropy gap by counts requires its own statistical analysis.

For a parity-blind joint $Y\otimes Y$ measurement the mean changes linearly, with slope $2(D_{\rm even,YY}+D_{\rm odd,YY})>1.53$. The plus probability changes by more than $7.65\times10^{-8}$ at the entropy-certified displacement. Comparing the two preparable endpoints $t=\pm0.01$ gives a more readily distinguished plus-rate difference exceeding 0.0153 while preserving every listed mean range. The larger endpoint comparison establishes readout separation; the finite entropy-increase bound uses the shorter one-sided segment.

## Six additional scalar means recover the missing coordinates

On the eighteen-coordinate slice, measure the parity-blind products **XZ, YX, YY, YZ, ZY and ZZ** on independent copies. This means measuring the indicated Pauli axis on the left and right coin and multiplying the two recorded signs. The parity weight and the two one-head marginals *in each parity block* are supplied in specifying this slice. Fixing only parity-pooled marginals permits additional unknowns.

The enrolled parked coefficient interrogation resolves the left reversal parity without changing the coin states within either outcome. Admit one material excess to each of the two isolated marked contacts. On the left the native 99-state Hamiltonian restricts to $8P$ or $9P$ according to the parity outcome, where $P=I-|B\rangle\langle B|$ and $|B\rangle=(|01\rangle+|10\rangle)/\sqrt2$ in the rail-and-coin basis. An exposure of $\pi/8$ on the even outcome and $\pi/9$ on the odd outcome gives $I-2P=-(Z\otimes Z)\operatorname{SWAP}$ on both. The right head is even and uses $\pi/8$. Apply the calibrated blank-dimer $X,Y,Z$ analyzers to the outgoing rails, average all parity outcomes and retain both product signs. This supplies the same product effect on both parity blocks with unit detection probability in the ideal controlled apparatus. It is a terminal destructive read, with no postselection. [Parked parity interrogation](../../action_response/README.md) · [Native contact reduction](../../material_clocks/PROOF.md) · [Rail analyzers](../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/measure.py)

For product index $j$ its additional row is exactly $2(e_j,e_j)$ on the eighteen orthonormal correlation coordinates. Append the six rows to $\bar R$ and call the augmented matrix $\bar A$. Their contributions are exact, so $\|A-\bar A\|\le30\epsilon$. An independently rounded rational inverse proposal $L$ satisfies
$$
\|I-\bar A L\|_F+30\epsilon\|L\|_F<2.369\times10^{-6}.
$$
Consequently the physical $18\times18$ observation matrix is invertible, with smallest singular value above $9.8966\times10^{-9}$. The original six-dimensional null space proves that fewer than six additional independent scalar means cannot suffice on this slice. This minimality statement counts scalar means. A multi-outcome instrument can report several independent means at once.

Use any known reference with the same parity weight and conditional one-head marginals, such as the product of the two marginals within each block. Subtract its predicted measurements and solve the augmented linear system. Add the resulting correlation coordinates to the reference to recover the full state on this slice, hence its entropy. The demonstrated Euclidean measurement-error to Hilbert–Schmidt state-error amplification is at most $1.011\times10^8$. The large factor matters if attempting a finite-shot state reconstruction. Autonomous entropy stationarity requires a separate dynamical result.

The state at the center is stationary under the equal-lottery reuse channel. Its neighbors are native feasible preparations. The result excludes stationarity against all such variations of entropy plus any differentiable function of the twelve pooled mean ranges at this source. Eighteen native charge-pooled return-time bins strictly before the cutoff suffice to recover all eighteen coordinates. The native clock retains information that twelve means discard. Other physical stress observables have their own effects. [Clock-history certification](CLOCK_HISTORY.md) [Regression](test_invisible.py)

From the repository root, run

    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.full_algebra.reuse.invisible
    OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra/reuse/test_invisible.py
