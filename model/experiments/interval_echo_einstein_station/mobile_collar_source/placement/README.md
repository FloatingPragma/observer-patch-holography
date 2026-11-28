# Do two freely formed source heads select the attractive sign?

The sign on the original home route is predictable from how the transmitter was armed. The 20-cell route begins with opposite sharp coin values. The 80-cell route begins with equal values. On **every directed bridge of all three boards**, at native exposures 1, 2 and 4, both load signs produce positive covariance between the two events “this head flipped”. Consequently, reading both original coin values on the same physical load standard gives positive covariance for equal initial coins and negative covariance for opposite initial coins. This accounts for the reversal between the two old home placements. It also predicts the sign on the other source placements before reading their final coins. [Original transmitter](../PROTOCOL.md) · [Prospective placement test and known outcomes](PROTOCOL.md) · [Subsequent census and source comparison](FOLLOWUP.md)

For original binary head readings $X,Y$ define $C=P(11)-P(X=1)P(Y=1)$. The older sharp Pauli covariance is $4C$. On the separately calibrated sharp contact, erasing the correlation while preserving the individual load distributions raises the ordinary energy by $17C$ at full admission. These computations follow the mobile vacancy or excess freely on every seam of its 240-, 960-, or 3,840-port carrier. Two physically addressed, initially independent memories and a localized defect are the supplied source stock. Their initial values direct both armed hops along the selected route. The native law preserves blank memory during free propagation. [Original formation operator](../measure.py) · [Energy readout and opposite physical orientations](../../../interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/README.md)

## Equal weights on source opportunities

Each owner has three outgoing intercell bridges. Give every directed bridge one vote and retain the inherited lowest-address internal neighbour at each end. The paired erasure comparison holds the placement fixed in both arms, so its mean energy change uses the **mean within-placement covariance**. Pooling all raw records before making the independent control introduces covariance between different placements' one-head means. On 80 cells, for example, the pooled vacancy covariance is +0.021324, while the local paired mean is only +0.000283. Neither board size nor the preparation weights have a measured physical-radius interpretation. [Full directed-bridge records and separate pooled contrasts](extended.json)

| Owners | Positive / negative placements at exposure four, either charge | Mean $C$, vacancy | Mean $C$, excess |
| ---: | ---: | ---: | ---: |
| 20 | 24 / 36 | −0.00163102 | −0.00246809 |
| 80 | 144 / 96 | +0.00028301 | +0.00008013 |
| 320 | 656 / 304 | +0.00697250 | +0.02711138 |

The 80-cell mean is particularly sensitive to source placement. Suppose each orientation class remains internally uniform but a fraction $f$ of launches starts with equal coins. Writing $A$ for their positive group mean and $-B$ for the other group mean, the sign changes at $f=B/(A+B)$. At four the thresholds on 80 cells are 58.56% for vacancy and 59.69% for excess. The inherited carrier supplies 60% under equal directed-bridge weights. The equivalent thresholds are 47.01% and 47.94% on the 20-cell carrier, where the actual share is 40%, and 43.01% and 25.16% on 320 cells, where it is 68.33%. These break-even values permit testing a separately determined source weighting. [Class populations and covariances](extended.json)

Opening the choice of arm admits five internal neighbours at each end. All 1,500 placements on 20 cells were run with equal weight. At exposure four the same 536 placements give positive covariance for either source charge; 964 give negative covariance. The means are −0.00346809 and −0.00520566. Independently sampling 240 of the 6,000 placements on 80 cells gives 64 positive and 176 negative for either charge. Sampling 240 of 24,000 on 320 gives 52 positive and 188 negative. The four larger-board sample means range from −0.006077 to −0.009716. Their conservative simultaneous 99% finite-population radii are 0.0711 and 0.0722 on the two boards, so the sign of either complete larger-board arm mean at four remains open. The flip-parity prediction correctly classifies every sampled arm at all three clocks and both source signs. [Held arm samples, every raw outcome probability and initial weights](extended.json) · [Initial 320-cell sample retained before its follow-up census](results.json)

## A short-time sign law for all 31,500 arm placements

Let $s_1,s_2$ be the initial sharp coins and $F_1,F_2$ their flip indicators. Since $X=s_1+(1-2s_1)F_1$ and similarly for $Y$,

$$ C=(1-2s_1)(1-2s_2)\operatorname{Cov}(F_1,F_2). $$

Write the four route ports as $u,a,b,v$. On every placement, among the two remote-arm endpoints only $b$ is reachable from $u$ in two hops, and only via $a$. The program verifies this on every arm pair on the three actual carriers. Thus the unique three-hop walk that flips the distant head has native amplitudes $-4,-1/2,-4$ and flips the near head as well. The first two Hamiltonian powers have zero projection onto a flipped distant coin. The third has precisely one such entry, $-8$ at $v$ with both coins flipped. For either defect sign the squared leading amplitude is $16t^6/9$. Real Hamiltonian coefficients and real sharp inputs make each outcome probability even in $t$. The product of the single flip probabilities starts at order eight. Hence $\operatorname{Cov}(F_1,F_2)=16t^6/9+O(t^8)$. [Native columns and two-hop census](measure.py) · [Independent Hamiltonian-power test](test_placement.py)

There is a uniform finite-time bound. Subtract 40 times the identity from the position and head generator. Its maximum absolute row sum is at most 24 for a vacancy and 29 for an excess, including private degree energies. Use the common norm bound 32. At $0<t\leq 2^{-20}$ the unitary Taylor remainder after three powers is at most $(32t)^4/24\leq t^3/24$. For both the remote-flip and joint-flip projections the cubic amplitude has norm $4t^3/3$, while the near-flip amplitude is at most $32t$. Consequently $\operatorname{Cov}(F_1,F_2)/t^6$ is between $(4/3-1/24)^2-(32\,2^{-20})^2(4/3+1/24)^2$ and $(4/3+1/24)^2$. This is strictly positive. The time interval sets the domain of this sign theorem. [Exact rational enclosures](extended.py)

| Owners | All-arm placements, equal / opposite initial coins | Certified interval for the all-arm mean $C/t^6$, rounded outward |
| ---: | ---: | ---: |
| 20 | 536 / 964 | [−0.619, −0.396] |
| 80 | 1,616 / 4,384 | [−0.933, −0.709] |
| 320 | 5,904 / 18,096 | [−1.016, −0.792] |

Thus a freely moving defect creates same-direction flip records on every transmitting arm at short times. Under equal weights on physically available arms, all three finite boards give a **negative** early mean energy change on the common original readers. Constraining the arms to the lowest-address choice reverses the sign of that early mean on 80 and 320 cells. A distribution over placements is part of the physical source claim. [Complete counts and signed rational means](extended.json)

## Can a physical preparation select a positive average?

Write **both heads to the same fixed sharp native coin** before the defect is launched. Two independent preparations of 0 give the product state 00. Two independent preparations of 1 give 11. Each initial covariance is zero. Use the same physical load reader on both heads, vary the placement without consulting its transmission orientation, and include every final outcome. Both preparations produce a positive average covariance for either sign of the defect at all three clocks, on every population in the table. This is a conditional source construction. The freely evolving law accepts the two armed heads as initial material. [Protocol for common physical arming](COMMON_HEADS.md) · [Full extension on the 20- and 320-cell carriers](common_heads.json)

| Population at exposure four | Vacancy 00 | Vacancy 11 | Excess 00 | Excess 11 |
| --- | ---: | ---: | ---: | ---: |
| 20 cells, all 60 inherited directed bridges | +0.00246550 | +0.00246732 | +0.00323934 | +0.00323243 |
| 80 cells, all 240 inherited directed bridges | +0.00243917 | +0.00244192 | +0.00314982 | +0.00314988 |
| 320 cells, all 960 inherited directed bridges | +0.00573087 | +0.00499264 | +0.01605713 | +0.01605825 |
| 20 cells, all 1,500 bridge and arm placements | +0.00217950 | +0.00218411 | +0.00293275 | +0.00293738 |

The mean has a definite sign even though some individual locations reverse it. On the 20-cell all-arm vacancy census, 00 gives 602 resolved positive, 58 resolved negative and 840 unresolved placements; 11 gives 864, 34 and 602. The largest negative 00 result is −0.0000297. With opposite fixed heads, 01 and 10, both mean covariances at four are negative for each of these same populations and both defect signs. A matched preparation therefore has a measurable effect apart from the reader convention. Under the condition of equal source opportunities and equal delivery and erasure factors, multiplying these means by 17 gives the separately calibrated contact-energy gain per opportunity. Location-dependent admission or exposure would weight the source differently. [Every sharp-input cohort](extended.json) · [Held additions to the complete censuses](common_heads.json)

Finally, two initially **independent fair coins** also have zero covariance. Mix all four evolved joint distributions at each placement before measuring their connected part. On the 20-cell inherited census the means are −0.00115776 and −0.00083786 at exposure four. On the 80-cell census they are +0.00188565 and +0.00618342; on 320 they are +0.00040426 and +0.00282578. Across all 20-cell arm choices vacancy and excess give opposite signs, −0.00006677 and +0.00076736. The two equal fixed standards give a repeatable way to construct a positive ensemble mean; a population of sources with unmeasured polarization, placement, age and delivery weights has no established sign from these results. Identifying what physically arms the heads in astronomical matter, where erasure exchanges energy and how the associated stress is measured by orbital and lensing probes remains open. [Complete source populations and control weights](common_heads.py)

## Controls, uncertainty and reproduction

The complete one-defect, two-armed-head sector has $4N$ states for $N$ ports. The native blank edges, both ternary load signs and both possible armed-edge orderings come from the 99-state local generator. With contrast removed all coin-sensitive blocks vanish. Uniform unit material leaves the two written heads at energy 16 times the identity. Replacing a formed joint state by its within-placement product gives exactly zero covariance. Each table entry uses all final positions and all four coin records. Slow propagation at smaller order could give an apparent sign on nearly insensitive alternative preparations, so the counts explicitly retain an unresolved category. [Full local-operator projection checks](../test_mobile_collar_source.py) · [Native propagation, controls and error calculation](measure.py)

The current source is based on repository revision `9ade2801a7d0b5ebe16a9fc01923b93926246833` with the existing worktree option to retain the scalar in the source Hamiltonian. Subtracting the same scalar on every state changes only the global phase. We used an order-256 Chebyshev expansion on the centered matrix divided by 32. Three sets of integer-argument Bessel coefficients come from rational partial sums with bounded remainders. The ignored Chebyshev tail is below $10^{-20}$. For IEEE binary64 arithmetic, at most seven matrix entries per row and absolute matrix norm at most one bound each local recurrence error by $2\times10^{-14}$. The Chebyshev polynomial of the second kind at step $j$ has norm at most $j+1$; accumulating the recurrence, coefficient, summation and probability-read errors gives an absolute radius below $5.70\times10^{-8}$ on each binary covariance. We compute the simulation outcome laws directly. The only sampling intervals here concern placements drawn without replacement. Relative to the certified 80-cell home route the independently applied exponential and the sharp read reproduce the earlier result. [Original integer certificate](../../../interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/formation_certificate.json) · [Independent full-state evolution and census checks](test_placement.py)

From the repository root, the first three commands reproduce saved results through each acquisition stage.

```sh
python -m model.experiments.interval_echo_einstein_station.mobile_collar_source.placement.measure
python -m model.experiments.interval_echo_einstein_station.mobile_collar_source.placement.extended
python -m model.experiments.interval_echo_einstein_station.mobile_collar_source.placement.common_heads
python -m pytest -q model/experiments/interval_echo_einstein_station/mobile_collar_source
```
