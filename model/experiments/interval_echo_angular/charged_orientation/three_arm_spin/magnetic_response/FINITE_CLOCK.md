# Current-meter bound for the original admitted records

The same original clock acceptance and circular preparations used by the [bond-current meter](PROTOCOL.md) give an upper limit on every native destination-probability difference at every free-evolution time. The limit applies after both source and clock acceptance. On the full original-start denominator, multiply it by the admission probability and keep both refusal bins.

| Guide links per arm | Blank-contact difference, upper bound | Written-contact difference, upper bound |
| ---: | ---: | ---: |
| 19 | $2.152\times10^{-6}$ | $2.017\times10^{-7}$ |
| 20 | $2.339\times10^{-6}$ | $2.193\times10^{-7}$ |
| 21 | $2.552\times10^{-6}$ | $2.392\times10^{-7}$ |

The calculation uses the admitted material chart, the original source and clock, free evolution under the native generator, the calibrated isolated contact pulses, and a parked donor. The bounds are for the expected finite record frequencies. A sampled comparison also needs its own uncertainty interval. The bounds concern vacancy transit currents on this chart. An implication about the exterior magnetic field requires a separate source-to-field theorem with an exterior observer read.

## From the clock to a uniform record bound

Let $H$ generate the native free evolution, $c$ be the [published clock center](../recoil_search/THRESHOLD.md), and $F_N(x)=N^{-1}\sum_{k=0}^{N-1}e^{-ik\tau x}$ with $N=2^{36}$ and $\tau=1/32$. The position and six recorded head bits index the admitted chart. A row-sum bound checks $\tau|E-c|<3<\pi$ throughout its spectrum. The geometric sum and $\sin y\ge2y/\pi$ for $0\le y\le \pi/2$ give

$$|xF_N(x)|\le\frac{\pi}{N\tau}.$$

The existing dyadic enclosure bounds the selected line's detuning and source overlap. If $\delta$ is the detuning ceiling and $b$ the overlap floor on the normalized source row, centering the sum at its midpoint gives an accepted amplitude floor

$$B=b\left(1-\frac{(\tau(N-1)\delta/2)^2}{2}\right)>0.$$

Every centered summand has real part at least the parenthesized quantity. The two orthogonal source rows have equal clock throughput by the checked arm symmetry, so this is also a floor for either normalized circular input. Thus for either normalized admitted circular packet $\psi$,

$$\|(H-c)\psi\|\le \varepsilon=\frac{\pi}{N\tau B}.$$

This norm is invariant under free evolution. For an orthogonal projector $P$, its mean derivative is twice the signed imaginary part of $\langle\psi,(P-I/2)(H-c)\psi\rangle$. The operator norm of $P-I/2$ is $1/2$, giving $|d\langle P\rangle/dt|\le\varepsilon$.

For each written bond, the projector onto the corresponding head bit being one changes only on that bond. Every admitted forward move flips that bit in the same direction. Its derivative is therefore the oriented current through that bond up to sign, so the current magnitude is at most $\varepsilon$.

The circular input is an eigenvector of the three-arm cycle. The calculation checks the exact rotation and reflection identities on the original integer source columns, and both the clock and the flight commute with that cycle. The three corresponding oriented blank currents are therefore equal at every free time. Cut the guide after ordinal $j$ on each arm. The projection onto the hub and the first $j$ sites of every arm has exactly three blank crossings and at most one written crossing per arm. Applying the derivative bound to that projection gives $3|J_{\rm blank}|\le \varepsilon+3\varepsilon$. Every blank bond belongs to one such cut. Divide the blank and written current bounds by their respective transfer magnitudes $1/2$ and $4$ to obtain the limits for the two destination probabilities.

The replay verifies the original integer preparation identities, exact permutation symmetries and microscopic bond assignments, every cut, and the quarter-pulse effects on the native local terms. Once the inherited overlap and energy enclosures are applied with an outward numerical pad, the scalar inequalities use rational arithmetic, $\pi<22/7$, and outward conversion for reporting. No wait-time grid is used in the bound.

From the repository root, run

~~~sh
python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.magnetic_response.finite_clock
~~~
