# Predicting heat flow from equilibrium records

Observers can predict the response of a guarded Interval Echo conductor from its heat-current fluctuations with both reservoirs at the same temperature. With four energy registers between the reservoirs, the prediction is **1/204** per appointment per unit inverse-temperature bias. The equilibrium measurement gives **0.0049010** after the finite-window correction. Reversing a small temperature difference in independent experiments gives **0.0049051** before the finite-bias correction. A disconnected conductor has measurable short-time fluctuations and zero steady response. [Contact protocol](PROTOCOL.md) · [Complete contact counts](records.json) · [Diffusion protocol](RING_PROTOCOL.md)

## Reading heat

Use one guarded dimer, with native energies 32 and 33, at each address. A fair mixture of sharp starting rails, the energy decoder, a timed BC leak and an inverse decoder give the three original-proposal probabilities
$(1/2,w/2,(1-w)/2)$ for lower, upper and refusal. The survival factor is $w=\cos^2(t/2)$. Its temperature is read from the upper-to-lower odds by a separate thermometer, $\beta=-\log w$ in unit-gap conventions. Five pulse settings have independent calibrations. The three-rail computation certifies these probabilities before sampling. Preparing the initial chain uses first acceptance at each site with every refusal counted. During transport each boundary offer uses exactly one fresh proposal, with refusal counted as an idle appointment. A contact swaps the two data rails by native exchanges on blank routes. The energy lost by the left source is its initial sign minus its outgoing sign, measured in unit gaps. Every completed swap conserves the sum of the two endpoint energies. A randomly selected address controls which bond can exchange; observations use the signs and the appointment tags. These choices describe a laboratory apparatus on the registered law. [Native path and swap proof](../thermal_contact/README.md) · [Pulse, route and census code](measure.py)

## Finite transport identity

Let there be $n$ interior registers and let $c_i$ be the probability per appointment of a completed exchange on bond $i$ (including the end contacts). Each interior rate is one uniform tag probability $q$. The end rates are $q(1+w_L)/2$ and $q(1+w_R)/2$. If $m_i$ is the mean upper population at site $i$, the mean current on a bond is exactly $c_i(m_i-m_{i+1})$, oriented to the right, with reservoir means $p_L=w_L/(1+w_L)$ and $p_R=w_R/(1+w_R)$. This identity follows from a complete swap and holds also for correlated interior sites. At stationarity all the currents coincide, so the signed influx is

\[
 J=\frac{p_L-p_R}{R},\qquad R=\sum_{i=0}^n c_i^{-1}.
\]

At matched sources put $p=w/(1+w)$, $v=p(1-p)$ and let $M=I-H$, where $H$ is the path Laplacian on the interior sites with the two endpoint rates to the sources. The product Bernoulli state is stationary, since every accepted exchange swaps independent identically distributed signs. Conditional means about this state evolve by $M$. For the left source's signed energy loss $j$ the equal-time fluctuation is $C_0=2c_0v$. Immediately after that appointment its covariance with the new interior signs is $c_0v$ at site one and zero elsewhere. The conditional mean next influx is minus $c_0$ times the endpoint deviation. Hence

\[
 C_k=\mathbb E(j_0j_k)=-c_0^2v(M^{k-1})_{11},\qquad k\ge1 .
\]

The substochastic walk with matrix $M$ reaches a reservoir from every site with positive probability in a bounded number of steps. Thus its powers decay geometrically, including when only the left contact is active. To solve the steady Dirichlet problem with left and right populations one and zero, take $m=c_0 H^{-1} e_1$. The left flux is both $c_0(1-m_1)$ and $1/R$. With one contact cut this flux is zero. Finally $dp/d\beta=-v$. Derivatives of the rates in $1/R$ multiply the vanishing equilibrium population difference, giving the exact identity

\[
 \tfrac12 C_0+\sum_{k\ge1}C_k
   =v(c_0-c_0^2(H^{-1})_{11})
   =\left.\frac{dJ}{d(\beta_R-\beta_L)}\right|_{0}.
\]

After lag $K$ the discarded signed sum is $-c_0^2 v(M^K H^{-1})_{11}$. It matters here. For four sites the 64-lag sum is 0.00557597 and the omitted term is −0.00067493. Cutting the right contact gives a 64-lag sum of 0.00462065 with an omitted term of −0.00461971. The four-site first-lag correlation was −0.00195322; treating successive transfers as independent would give about 0.02083 for the intact conductance. This edge-current identity uses a signed transfer on an appointment. The finite theorem for the resolvent pairing of state observables uses a different observable. [Finite equilibrium theorem](../../../../../../Lean/Thermodynamics/GreenKubo.lean) · [Thermal contact and the equilibrium heat-current relation](https://arxiv.org/abs/1302.4538)

## Measured response and uncertainty

At baseline the reservoir odds are 1/2. The exact coefficients for chains of one, two and four registers are 1/48, 1/66 and 1/204 in gap units per appointment per unit thermal bias. A forward and reverse pair at independently calibrated odds near 9/16 and 4/9 gave the raw slopes **0.0208552, 0.0151651 and 0.0049051**. The second pair used 5/8 and 2/5. The right-cut four-site response at this larger spacing was **0.0000010**. The calculation at each finite spacing uses the actual admission probability of each source. The test bounds the warmup error after 1,024 appointments by $c_0\sqrt n\,\rho(M)^{1024}$. It bounds thermometer error separately from heat-current error and compares the two-sided response to the equilibrium prediction including the calculated finite-bias shift.

Each transport bank uses exact multinomial propagation of counts of indistinguishable words. Its $2^{37}$ independent simulated chains occupy at most 16 interior words per initial current class. These are compressed draws from the operator-certified process. The receipt stores every comparison histogram and totals for all offer types, including idle and refused events. With fewer than 1,024 scalar comparisons, define $h(N)=\sqrt{\log(2\cdot1024/0.01)/(2N)}$. Hoeffding and a union bound give simultaneous 99% coverage. The 64-lag sum has a conservative sampling allowance of $(128+1/2)h(N)$, below 0.00087; summing the marginal bounds uses no independence between successive lags. The receipt includes total offers across the 1,024 warmup appointments and final read, original-preparation refusals, and the final counts by word and offer. [Independent checks](test_transport.py)

## Diffusion without a temperature gradient

A separate four-site ring contains exactly one upper sign in a sea of lower signs. It uses the same three interior exchange routes and the same tag probability per offered bond as the conductor, including a closing route in place of the boundary contacts. Reading consecutive upper-sign positions gives an unwrapped signed displacement, including wraparound. Conditional on any current position, the next displacement is +1 or −1 with probability $q=1/8$ each, or zero otherwise. The increments are independent because their conditional distribution is constant for every prior history. The diffusion coefficient in squared bond spacings per appointment is therefore $D=\langle X_t^2\rangle/(2t)=q$, and equals the half-weighted zero-lag current fluctuation plus the sum of the positive-lag correlations. Eight-step counts on 16,777,216 fresh simulated preparations gave **0.125037**. The parked ring gave zero at every step. In the four-register conductor the three interior resistances are then independently predicted by this $D$, and the two end resistances by the thermometer admissions. Their combination with the measured equilibrium susceptibility predicts the driven conductance within its joint error. The displacement is unwrapped along the supplied cyclic order; the finite ring's occupied-address distribution has bounded variance. [Ring source](ring.py) · [All displacement and current bins](ring.json)

The results exhibit diffusion and linear thermal response on a scheduled chain of guarded energy registers. A bulk conductivity of freely interacting matter requires its own current, distance and time calibration; identifying a macroscopic material coefficient in SI requires a physical unit map.

From the repository root, with one numerical thread:

```sh
pkg=model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.heat_transport
.venv/bin/python -B -m "$pkg.measure"
.venv/bin/python -B -m "$pkg.ring"
.venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/refinement_cooling/dimer_band/heat_transport
```
