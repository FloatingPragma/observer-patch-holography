# How much of a formed collar does a regional observer read?

Two armed seams acquire correlations from a freely travelling excess or vacancy. At unit filling the parked heads have total energy 16 in every joint state. The observer can combine heads from the same preparation or from independent preparations. Both readings give exactly the same single-head density matrices. The correlation changes the modular first variation and can change the return of a second material probe. [Native formation](../mobile_collar_source/README.md) · [Protocol](PROTOCOL.md) · [Full outcome laws](results.json)

## An exact finite comparison

Write \(\sigma\) for the formed two-head density matrix, \(\pi=\sigma_A\otimes\sigma_B\) and \(D=\sigma-\pi\). With pairing probability \(f\) the reference is \(\rho_f=\pi+fD\). Vary the source toward stronger pairing at fixed marginals. Per unit transfer the observed tangent is \(fD\). This specifies a direction in density matrices, with the opposite sign for a transfer toward weaker pairing. At the formed endpoint the inward variation is implementable simply by mixing in independent heads. The logarithm is evaluated on the fixed unperturbed reference. All four prepared sources have full rank.

At \(\pi\), the modular operator is a sum of two one-head terms. Thus its expectation has zero derivative on **every** fixed-marginal perturbation, including coherent ones. The identity \(S(\pi)-S(\rho)=D_{\rm rel}(\rho\Vert\pi)\) makes the product the unique entropy maximum at these marginals. At the paired reference the signed first variation is
\[
 E(f)=-f\,\operatorname{Tr}D(\log(\pi+fD)-\log\pi).
\]
Let \(R_t=(\pi+tI)^{-1}\), \(J=\int_0^\infty\operatorname{Tr}(D R_t D R_t)\,dt\), and \(b=\lambda_{\max}(\pi^{-1/2}D\pi^{-1/2})\). For a fixed source and every \(0<f<1\), the finite bounds are
\[
 \frac{J}{1+bf}\ \leq\ \frac{-E(f)}{f^2}\ \leq\ \frac{J}{1-f}.
\]
Indeed, \((1-x)\pi\leq\pi+xD\leq(1+bx)\pi\). After adding \(tI\) and inverting, the resolvent at \(x\) lies between \(R_t/(1+bx)\) and \(R_t/(1-x)\). The quadratic form \(\operatorname{Tr}(D R D R)\) has positive matrix \(R^T\otimes R\) on vectorized \(D\). The corresponding squared bounds apply to the curvature of the logarithm at \(x\). Integration from 0 to \(f\) gives the displayed bounds. For \(D\ne0\), \(J>0\). Two independent uniform draws of indices from a bank of \(N\) copies give \(f=1/N\) exactly. Keeping the same index gives \(f=1\) at every bank size.

| Formation board | Forming load | \(J\), nats | Strengthening slope at \(f=1\), nats |
|---|---|---:|---:|
| 20 cells | vacancy | 0.00500 | −0.001924 |
| 20 cells | excess | 0.01046 | −0.004173 |
| 80 cells | vacancy | 0.11146 | −0.200825 |
| 80 cells | excess | 0.13209 | −0.293848 |

## What a clocked ball changes

A uniform launch hits a ball with probability \(f\) equal to its measured fraction of the host ports. The lab supplies paired heads on a hit and independent heads on a miss. If the regional algebra includes the hit bit, its density is \(f\sigma\oplus(1-f)\pi\), with first variation \(f E(1)\). If the regional algebra contains just the heads, the first variation is \(E(f)\). Both prescriptions keep every trial. Tracing out the bit changes the information accessible to the observer.

For the 327,680-cell survey, doubling the inner clocked range changes the enrollment from 1,158 to 4,630 ports. The flagged slope has doubling exponent 1.9994. The four unflagged exponents are 3.9773–3.9954. The smallest 80-cell survey gives only 2.61–2.99 for the unflagged slopes. These calculations combine separate source copies and a native range census. A fourth-power source on a proper ball would require the corresponding physical observation algebra, \(f/\ell^2\to a\), and a compatible family with controlled \(J\) and \(b\). At fixed source this gives \(E(f)/\ell^4\to-a^2J\) exactly. The sign reverses on weakening the pairing. Attractiveness requires an independent curvature response to the physically selected variation. [Clocked-ball census](../proper_ball_ruler/coarse.json) · [Dark-sector normalization](../../../../cosmology/oph_dark_matter_paper.pdf)

## Energy and round trips on the same collar

A fresh material probe starts at the inner endpoint of one marked seam. A separate Ramsey comparison can read its complex return amplitude against a phase reference. The reciprocal test interrogates the other inner endpoint after time \(t\), then checks for return at \(2t\). If \(B\) and \(C\) are the outgoing and incoming four-by-four propagator blocks on the heads, the complete probability is \(\operatorname{Tr}[(CB)^\dagger CB\,\rho]\). The receipt includes the miss on each leg. The phase reference, projective relay and appointments are supplied apparatus. These are calculated population probabilities.

For any initial source-independent probe, the first two energy moments are identical between the paired and independent arms. Write the native operator as a head-independent term plus one contact per head. The contacts have disjoint spatial supports, so both cross products vanish in the square. Every surviving head term depends on at most one marginal. This proof holds at every board size with these separated contacts. A single positive rescaling of the whole Ramsey time would multiply the nonzero variance by the square of that rate, even allowing an arbitrary energy offset. Equal variance fixes that rescaling to one, and equal mean fixes the offset to zero. The later complex return amplitudes differ.

On the 80-cell board with an excess forming load and excess probe, the complete return probabilities at two units per leg are 0.038839 for paired heads and 0.053209 for independent heads. Their initial energy variance is 5.26308 in both arms. This isolates a delayed material sensitivity at equal ordinary head energy. An inertial mass, a distant static force, a proper clock and an optical null path require separate identifications.

Replay both source panels and the receipt tests from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.interval_echo_einstein_station.collar_pairing.measure
python -m pytest -q model/experiments/interval_echo_einstein_station/collar_pairing
~~~
