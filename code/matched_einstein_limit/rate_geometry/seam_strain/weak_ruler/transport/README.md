# A position clock on the traveling probe

A mobile launches from port 0 and is read at port 39, or travels in the reverse direction. At an intermediate appointment a local detector tags whether it occupies the first interior port of the marked chain, 1 or 38 respectively. Both original heads and the bridge strain remain coupled during flight on the four-tick pumped, 80-cell source. A separate bank measures arrival with the intermediate detector absent. The experiment includes both vacancy and excess probes, both intermediate outcomes, and every terminal miss. Times use the native exchange rate.

Can the local mark be interpreted as a passive reading of where the carriers that reach the receiver spent their time? With mark time 4 and terminal time 5.5, the following signed readings are inferred per terminal arrival on the reference source.

| Carrier | 0 to 39, mark 1 | 39 to 0, mark 38 |
| --- | ---: | ---: |
| Vacancy | −0.06047 | −0.06661 |
| Excess | −0.06458 | −0.07755 |

The signs also hold at times 2 and 4.5, on the full ordinary rail change and on the held preparation with independently re-paired heads. All 24 population enclosures lie strictly below zero. All hardware outcome probabilities are nonnegative. [Positive event banks and interval bounds](results.json)

## How the clock is read

Let \(P\) select the marked port, \(F\) the receiver and \(U\) the free propagator, with all internal states retained. Restricted to the launch and reception registers write

\[
D=F U(T),\qquad V=F U(T-s)P U(s),\qquad B=D-V .
\]

For the input state \(\rho\), write \(p(X)=\operatorname{Tr}(X\rho X^\dagger)\). The unmarked arrival probability is \(p(D)\). A sharp mark followed by the receiver gives two exclusive arrival branches \(p(V),p(B)\). The signed joint position clock and its normalized reading are

\[
  W=\operatorname{Re}\operatorname{Tr}(V\rho D^\dagger)
  =\tfrac12[p(D)+p(V)-p(B)],\qquad w=W/p(D).
\]

A variable-strength binary pointer gives the same reading with arbitrarily small disturbance to the total arrival rate. Its two arrival amplitudes are \(aV+bB\) and \(bV+aB\), where \(a^2=(1+\epsilon)/2\), \(b^2=(1-\epsilon)/2\). For the two positive event probabilities \(p_+,p_-\), \(W=[p(D)+(p_+-p_-)/\epsilon]/2\) at every \(0<\epsilon\le1\). The unconditioned terminal disturbance is proportional to \(1-\sqrt{1-\epsilon^2}\). A classical passive tag on pre-existing visits with the same calibrated response would give \(W=\Pr(\text{visit and arrive})\ge0\). Here the coherent through-site and outside-site alternatives interfere. The sharp reader changes the receiver statistics. [Two-time quasiprobability and its measurement](https://doi.org/10.1103/PhysRevA.93.022123)

The discrepancy persists on intervals of time. At fixed terminal appointment, differentiating the Heisenberg site projector gives \(|\partial_s W|\le \sqrt{p(D)}\,\|[H,P]\|\), by Cauchy–Schwarz on the input state. The integer cut Hamiltonian bounds \(\|[H,P]\|^2=\|PH(1-P)\|^2\). Every tested reading remains negative throughout a shutter of half-width \(3\times10^{-5}\) ticks. The individual certificates often allow longer windows. Integrating against any positive affine time weight on this window preserves its sign.

Dephasing the two path alternatives gives a quantitative control. If a physical pointer preserves the diagonal path contributions and scales the real interference by \(g\in[0,1]\), then \(W_g=p(V)+g[W-p(V)]\) and \(p_g(D)-p(D)=2(W_g-W)\). Restoring a nonnegative reading requires a **relative increase in arrivals of at least \(-2w\)**. On the reference source at times 4 and 5.5 this exceeds 12.09% for every charge and direction. An invasive positive clock can be implemented and its changed arrival statistics measured.

## Held source and independent checks

Holding the five-station radar volume fixed by the independently measured right-head rotation gives a separate infinitesimal test. In the forward vacancy experiment at times 4 and 5.5 the ordinary change in \(w\) is 0.00544139 per unit admixture. Calibrating it against the entropy of the retained heads and strain predicts +0.00314014 to +0.00314024 for re-pairing. The measured re-pairing slope lies between −0.000294933 and −0.000294928. All eight ordinary-calibrated residuals exclude zero. The clock measures conditional propagation on this material. Its role in a null-geodesic description would require an independently identified ray, proper clock and screen, together with a physical stress measurement. Coarse ensembles and different crossing instruments can have other responses. [Volume enclosures](../hinge.json) · [Independent source logarithm](../../live_modular.json)

The two appointment pairs were selected from the reference half-tick survey. Each pair was then used on both held preparations. Finite bounds use integer Hamiltonian blocks and trace-norm enclosures of the original pumped source. The simulation also draws \(2^{40}\) trials per instrument and setting. A simultaneous 99% Bernstein interval from the retained counts resolves every sign. These are simulated event banks. Independent propagation evaluates the finite-strength pointer at four strengths. [Implementation](measure.py) · [Replay and tests](test_transport.py)

From the RER root, with the model environment active, run

~~~sh
PYTHONPATH=code:. python -m matched_einstein_limit.rate_geometry.seam_strain.weak_ruler.transport.measure --write
PYTHONPATH=code:. python -m pytest -q -W error code/matched_einstein_limit/rate_geometry/seam_strain/weak_ruler/transport
~~~
