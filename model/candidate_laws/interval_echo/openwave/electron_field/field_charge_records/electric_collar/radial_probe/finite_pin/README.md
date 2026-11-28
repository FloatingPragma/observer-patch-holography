# A connected spectral shift with finite charge pins

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

The four-fan arrival at hub 1045 supports a connected near–far energy difference when both loads can make virtual excursions from finite traps. At native ring rate, increasing both pin penalties from 64 to 512 reduces the difference from about $-0.003862$ to $-0.00006903$ in ring-clock energy units. The exact finite theorem bounds it quadratically in inverse pin strength on every enrolled head word. The unfiltered energy mean of the original preparation has zero connected difference, and the finite shift is visible in the lowest reachable line of a supplied ideal spectrometer. [Protocol](PROTOCOL.md) · [Full receipt](results.json) · [Previous clamped and dynamical measurements](../README.md)

## Matching the four measurements

Pool the untouched reference outcomes of the native four-fan 30.8125 arrival at face 0. The face pair is 334/365 or its charge complement. The face filter admits about 1.58437% of original launches. The independent material probe sits at 1048 or 1046, two or three intact-board hops from 1045. The two probe ports have the same degree. Calibrate each pin separately by acquiring both phase quadratures of an isolated off-pin load with the actuator on and off. Then run four banks with both units, source only, probe only and neither. For the removed-source bank move the arrived source back to 1010 with the reserved bridge and shutter it. For the removed-probe bank use a separate same-rate feeder to the common exterior port 1055 and shutter it. Neither removal changes the enrolled collar heads. Separate reservoir phases and either unit's isolated energy cancel in $U=E_{11}-E_{10}-E_{01}+E_{00}$. The radial score is $U_{\rm far}-U_{\rm near}$.

The meter reads the lowest reachable energy of each bank. The two input words belong to one connected component in each bank. The connected off-diagonal entries are nonpositive, so Perron positivity and the full-rank arrival density give nonzero overlap with each stationary line selected this way. The receipt retains each calculated line energy, its numerical enclosure and the nominal probability per original launch. Everything outside the selected line contributes to the remainder, including original face refusal. No energy apparatus was built, and the lowest-line weights are calculated from the declared ideal measurement. Using the energy mean of all spectral outcomes conditional on the same initial face cut would instead give exact zero, because every bank's original preparation lies in the same matched pinned block.

For the positive source and probe, near word 45 and far word 47 differ only by the head flip of the admitted native direct probe hop. This preserves the probe load–head sector. This example was selected from the complete comparison of both signs, 64 exterior words and the pin settings. All intervals below enclose numerical error in the ring-clock energy comparison.

| Pin penalties, source and probe | Connecting ring rate | $U_{\rm far}-U_{\rm near}$ |
| --- | ---: | ---: |
| 64, 64 | 1/2 | $[-0.003861951,-0.003861936]$ |
| 128, 128 | 1/2 | $[-0.001045631,-0.001045607]$ |
| 256, 256 | 1/2 | $[-0.000271248,-0.000271204]$ |
| 512, 512 | 1/2 | $[-0.000069068,-0.000068983]$ |
| 64, 64 | 1/4 | $[-0.004135120,-0.004135105]$ |
| 64, 64 | 0, reachable component | $0$ exactly |

At native 64, the nearest both-present line has nominal weight $9.75\times10^{-6}$ per original launch. At native 512, the worst outside-pin probability of the stationary levels over the complete banks is below 0.027%. Swapping the probe sign and interchanging the two equi-degree port labels is an exact symmetry of every term, with the sharp head word unchanged. It reverses the same-word radial score by exchanging its endpoints. With the near address and sharp word 45 held fixed, the two opposite probe signs give a **sum** $U_{++}+U_{+-}$ in $[0.0031055934,0.0031056103]$. A pure charge-product energy would give zero for that sum. Joint charge conjugation and head complementation also hold term by term.

## What the pin sweep establishes

In the subspace with every retained load at its pin, let $A$ be the unpenalized operator. Its coupling to the off-pin block is $K=[K_s\ K_p\ 0]$ in the source-away, probe-away and both-away ordering. For every arm and reached component the executable checks that the four pinned word sets coincide. Source retention adds $12I$ to $A$ and gives $K_sK_s^T=48I$. The probe's $K_pK_p^T$ is the same diagonal with or without the source, with entries 0 or 16. Therefore the lowest level of

$$
A_{\mathrm{eff}}=A-\frac{K_sK_s^T}{\lambda_s}-\frac{K_pK_p^T}{\lambda_p}
$$

has the same source-only offset $12-48/\lambda_s$ at both probe radii and in both probe presence banks. This proves cancellation through first order in inverse pin strength even when $\lambda_s\ne\lambda_p$.

Here is the explicit remainder. Write the off-pin block as $D+\Lambda$, with $D\succeq0$ from the original graph Laplacians and contact energies. Put $L=\min(\lambda_s,\lambda_p)$, let $e$ be the least diagonal of $A$ and take $d\ge\|D\|$. The full lowest energy $E$ lies in $[0,e]$ and $e<L$ on all measured arms. The Schur complement at $E$ is positive semidefinite with a kernel. Thus $E$ is the lowest eigenvalue of $A-K(D+\Lambda-E)^{-1}K^T$. The identity

$$
(D+\Lambda-E)^{-1}-\Lambda^{-1}
=-(D+\Lambda-E)^{-1}(D-E)\Lambda^{-1}
$$

bounds its difference from the lowest level of $A_{\mathrm{eff}}$ by $\|K\|^2(d+e)/[L(L-e)]$. Summing four such bounds controls $U$, and summing the independent radial bounds controls the measured difference. The same off-pin equation gives the bound $\|K\|^2/(L-e)^2$ on outside-pin probability. The bound and exact first-order cancellation give the infinite-pin null. Cutting the connector also gives exact null for every pin strength, since the actually reachable components factor into disjoint internal and external operators with Cartesian-product reachability.

## Weak coupling and visibility

On the same selected pair of radial sectors at pin 64, weakening the connector to $1/65536$ gives $U_{\rm far}-U_{\rm near}\approx-0.00752913$. All positive rates use the same reachable chart. The limiting difference on that chart is $-0.0075366085$ with numerical error below $4.8\times10^{-9}$. At zero connector the line actually reachable from the same initial words has exact connected score zero. The two nearest weak-coupling lines with the probe present flow to closed components inaccessible from the original preparation. Their limiting floors lie at least 0.26129 and 0.25375 below the cut-accessible floors for the probe-only and both-present banks, respectively.

The norm of the removed connector operator is $2g$ at rate $g$. The limiting cut floor on the coupled chart increases by at most $2g$. If the cut-accessible floor exceeds this continued energy by a certified $d_g>0$, projection of the eigenvalue equation bounds the eigenvector's norm on the original component by $2g/d_g$. Multiplying its square by an upper bound on face admission carries the original launch denominator through the line comparison. At $g=1/65536$ each of the two lowest lines occurs with probability below $3\times10^{-10}$ per original launch. The sharp lowest-line selection therefore retains a nonzero energy difference in the weak limit while the probability of reading the relevant lines vanishes.

This is a finite virtual-motion energy of the supplied collar and traps. It supplies neither a stationary force measurement on a freely responding probe nor a long-range law or an independently established elementary charge strength. Rows 1 and 20 remain on **HOLD**.

From the repository root, run

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.electric_collar.radial_probe.finite_pin.measure
~~~

The independent audit reselects the connected head-and-load component by adjacency search and bounds the energy of each arm by exact rational row ratios of a positive trial vector. It checks the native 64 and 512 endpoints, the half-rate and cut controls, and the fixed-word polarity comparison. [Exact-ratio audit](audit.py)

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.electric_collar.radial_probe.finite_pin.audit
~~~
