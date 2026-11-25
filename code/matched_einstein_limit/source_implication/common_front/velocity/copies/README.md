# Can repeated observers see the source in a common geometry?

Independent native repetitions detect both ordinary energy and the loss of pairing at this station. They also measure a shorter mean flight distance on both held-volume paths. The two material probe species infer incompatible changes from that distance after ordinary calibration. The exact velocity distance itself is unchanged for every number of independent copies. These results concern the stationary 80-cell source, time-four flight and the registered geometry measurements below. [Native source](../../../full_algebra/PROOF.md) · [Signal-speed geometry](../README.md)

## A material detector for lost pairing

Use the installed marked contacts for a dwell $\pi/8$. Load an excess on the first port of each material readout rail. Read the left output in $Z$ and the right in a uniformly chosen axis $j=X,Y,Z$ using blank-seam exchange and phase operations. Multiply the two signs. Repeat on a second independently prepared source with the **same** chosen axis, and record one click if the two products agree. Every run has a binary outcome. The apparatus uses the right contact's even reversal sector in every preparation. The left marked contact may occupy either parity. [Native readout gates and position analyzers](../../../../../../model/experiments/interval_echo_einstein_station/collar_pairing/phase_geometry/measure.py)

On the mobile contact the exchange frequencies on the active and spectator components differ by 8 in the even sector and by 9 in the odd sector, in native time. At the registered dwell the even gate transfers the coin to the rail up to a known local sign convention. The odd left rail sees $L_o=\mathrm{diag}(1,-\cos(\pi/8))$ on the head. This follows by diagonalizing the middle two-port exchange block; its probability of changing the material output port is $\sin^2(9\pi/16)$ on input coin one and zero on input coin zero. The left even analyzer sees $L_e=Z$. The right analyzers read $\sigma_j$ with the calibrated sign. Each run uses these same physical gates regardless of its stored parity. [Contact matrices and changed-contrast control](test_copies.py)

Let $\chi_b$ be the subnormalized stationary head state in sector $b$ and $\chi_R$ its right marginal. Write $g_j=\sum_b\operatorname{Tr}\chi_b(L_b\otimes\sigma_j)$, $l=\sum_b\operatorname{Tr}\chi_b(L_b\otimes I)$, $v_j=\operatorname{Tr}\chi_R\sigma_j$ and $h=lv$. Erasing pairing with probability $s$ replaces each $g_j$ by $(1-s)g_j+sh_j$. Any compensating right-unitary rotation applies the same orthogonal rotation to the resulting vector. Conditional on $j$, the two independently prepared products agree with probability $(1+g_j^2)/2$ at the reference. Therefore the changed coincidence probability is exactly

$$
p_s=\tfrac12+\tfrac16\|(1-s)g+sh\|^2 .
$$

On the full stationary enclosure, the unconditional native meter gives $p_0\in[0.53487160,0.53487186]$ and $p'_0\in[-0.069789961,-0.069789459]$. The change at the archived $s=2\times10^{-7}$ endpoint lies in $[-1.395799071,-1.395789050]\times10^{-8}$, including any installed right rotation. Switching off squared contrast during transfer while keeping the same frozen detector erases this distinction. An observer can infer a pairing change by repeated material position reads. [Enclosed probabilities](results.json)

We also exercised a parity-resolved implementation with complete synthetic counts. A parked left exchange parity interrogation with a relative integrated exchange pulse $\pi$ resolves $I-Q,Q$ without disturbing each sector internally. On odd outcomes the apparatus returns an independent fair sign. On even outcomes it uses the two rails just described. Thus the same calculation uses $L_e=Z,L_o=0$; the new slope is $[-0.044620156,-0.044619954]$. For the ordinary preparation a single parity read suffices. Its success probability increases from $0.202437$ toward $1/2$ with slope $[0.297562763,0.297563565]$; unpairing preserves it exactly. [Parked coefficient interrogation](../../../action_response/README.md) · [Binary records](results.json)

## How many independent sources?

For a register with known direction of change and physical probability gap at least $d$, draw $n$ blinded Bernoulli readings and $n$ contemporary reference readings. Declare a change when the signed difference in proportions exceeds $d/2$. Hoeffding's inequality bounds each of the false-alarm and missed-change probabilities by $\exp(-nd^2/4)$. The receipt takes $n=\lceil19/d^2\rceil$ with a rational, downward-rounded $d$, and supplies seeded reference, unchanged and changed binary counts. A given comparison uses only two cohorts. Each parity observation consumes one preparation; each collar coincidence consumes two. The counts are synthetic acquisition data from the model. The error probabilities and source enclosures follow from the stated independent-preparation experiment. [Sampling implementation](measure.py)

| Source lottery $s$ | Ordinary observations per cohort | Parity-resolved collar coincidences per cohort |
|---|---:|---:|
| 1 | 215 | 38,217 |
| $1/8$ | 13,734 | 694,970 |
| $2\times10^{-7}$ | 5,364,589,020,743,195 | 238,580,832,913,902,450 |

There is also a lower bound allowing an arbitrary collective quantum measurement of all $M$ tested preparations and arbitrary unchanged ancillary inputs. For a faithful reference $\chi$, changed head $\tau$, and $\Delta=\tau-\chi$, put $Q=\operatorname{Tr}(\Delta\chi^{-1}\Delta)$. For each positive measurement effect $E$, Hilbert–Schmidt Cauchy–Schwarz gives $|\operatorname{Tr}E\Delta|^2\le\operatorname{Tr}(E\chi)\operatorname{Tr}(E\Delta\chi^{-1}\Delta)$. Sum over outcomes. The measured classical $\chi^2$ divergence is at most $Q$, and $1+Q=\operatorname{Tr}(\tau\chi^{-1}\tau)$ multiplies under independent products. A binary decision with both errors bounded by $\epsilon<1/2$ has classical divergence at least $(1-2\epsilon)^2/[\epsilon(1-\epsilon)]$. Consequently

$$
M\ge\frac{\log[1+(1-2\epsilon)^2/(\epsilon(1-\epsilon))]}{\log(1+Q)} .
$$

For the raw lotteries $\Delta=sD$; the receipt propagates the station enclosure through $Q/s^2$. For the archived installed rotation $\theta$ it independently uses $\|\Delta\|_{\rm HS}\le s\|D\|_{\rm HS}+|\theta|$ and $\chi\ge0.019999 I$. With $\epsilon=0.01$ and $s=2\times10^{-7}$, even these conservative bounds require over **2.25 million ordinary tested preparations** and **26.17 million collar tested preparations** at the installed angles. For a smooth fixed-volume continuation the score bound $-B\chi\preceq\dot\chi\preceq B\chi$ implies $\operatorname{Tr}(\dot\chi\chi^{-1}\dot\chi)\le B^2$. Thus the necessary number is $\Omega(s^{-2})$. If one proposes a head perturbation of order $\ell^4$ as a way to encode a shrinking-ball coefficient, this entails $\Omega(\ell^{-8})$ preparations on that proposal. A macroscopic ball does not by itself prescribe that head family or scaling. [Source-score bounds](../certificate.json) · [Finite adjusted interventions](../../../signed_transfers.json)

## What geometry do these repetitions supply?

Let $H_M=\sum_i H_i$ and register the intensive observable $\bar F=M^{-1}\sum_iF_i$ for $M$ independent copies with the same cell positions. The terms in $i[H_M,\bar F]=M^{-1}\sum_i i[H_i,F_i]$ commute. A product of bridge states at a maximizing edge attains the same extreme eigenvalue on all copies. Therefore $\|i[H_M,\bar F]\|=\|i[H,F]\|$ for every $M$, even when arbitrary entangled initial states are allowed in taking the norm. The dual cell distance remains $2d_{\rm cell}$ and any fixed area prescription using these distances alone is insensitive to the source. This argument uses independent copies, a common position assignment and intensive readouts. Interacting distributed heads and relational coordinates derived from other records are distinct constructions. [Single-copy commutator proof](../README.md)

A different experiment measures **mean travel** by pooling sharp terminal cell clicks. Release either charge from the same sharp port 0, evolve freely for four exchange-time units, then register terminal cell $C$. Draw a shell $k$ uniformly from 1 through 4. Score one precisely when $d_{\rm cell}(0,C)\ge k$. Eight times the success probability is exactly $\mathbb E[\min(2d_{\rm cell}(0,C),8)]$ in the velocity-dual distance units. Independent repetitions converge to that value. The fixed-point propagator encloses all 3,840 terminal amplitudes in each parity sector for four launch columns. The signed variations use the original station and the compensators of the original radar volume. The ratio enclosure cancels the common rotation between the two interventions. [Certified radial flight](radial.py) · [Independent shell and log-slope replay](../results.json)

| Carrier | Mean capped distance | Ordinary distance slope | Collar distance slope | Collar/ordinary |
|---|---:|---:|---:|---:|
| Vacancy | [0.04057541, 0.04057556] | [−0.113033, −0.111310] | [−0.045984, −0.045290] | [0.40521410, 0.40892915] |
| Excess | [0.05926173, 0.05926192] | [−0.040398, −0.039779] | [−0.022045, −0.021826] | [0.54357806, 0.55130120] |

Suppose the two meters respond linearly with separate, fixed gains to a common scalar. Calibrate each gain on the ordinary transfer. Their held collar/ordinary ratios must then coincide. The certified intervals are disjoint. For $a=0.40892915$ at the nearest vacancy endpoint and $b=0.54357806$ at the nearest excess endpoint, any fitted common fraction incurs at least $(b-a)/(b+a)>14.1362\%$ relative error on one charge. This follows by minimizing the maximum of the two relative deviations; the optimum lies between the endpoints. Both mean flights really shorten, and independent repetitions resolve this. A universal proper radius or a common curvature coefficient has yet to follow from these readings. This finite test covers the specified time, sharp preparation, terminal meter, fixed survey volume and static gains. [Outward interval receipt](radial.json)

From the repository root,

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.velocity.copies.measure --write
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.source_implication.common_front.velocity.copies.radial --write
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/common_front/velocity
~~~
