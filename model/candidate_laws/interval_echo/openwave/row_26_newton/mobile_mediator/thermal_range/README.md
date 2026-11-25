# Temperature and range of the neutral-dipole interaction

Two independently randomized neutral dipoles lower the first-order mobile free energy at every finite inverse temperature \(\beta>0\) whenever a native mobile can visit both. The strength decays at least exponentially in their mobile hop separation, uniformly in temperature on the degree-five/six boards. The attractive part obeys the same range bound at finite activity when the fair internal states equilibrate before taking the logarithm. On the spherical midpoint family the entire connected energy is smaller than inverse cell distance asymptotically, including the positive contribution from competition for the mobile. A fixed activity and an arbitrary temperature schedule therefore cannot produce a nonzero Newton tail through this mediator. [Four-arm preparation and interferometric spectrum](../PROTOCOL.md) · [Native one-mobile Hamiltonian](../measure.py)

## Signed response at every temperature

Use the finite blank-memory source with the two clamped dipoles of the equal-energy experiment. Let \(G\) be adjacency on the unpinned sites, \(D_0\) the native diagonal relative to each arm's stationary zero-mobile background, and \(H_{ij}=D_0+iV_A+jV_B-G/2\) for \(i,j\in\{0,1\}\). Here \(V_A,V_B\) are diagonal. Let \(d\) be the shortest hop count between their disjoint neighborhoods, \(s\) an upper bound on the size of either neighborhood, and \(Z_{ij}=\operatorname{tr}\exp(-\beta H_{ij})\). Use \(\Delta Z=Z_{11}-Z_{10}-Z_{01}+Z_{00}\). Bars below average the two independently fair dipole orientations and, if mixed, the mobile charge. The four arms share a mobile graph and clock. The Ramsey return spectrum supplies the conditional frequencies and weights. The heat trace additionally uses a specified thermal weighting of that spectrum. [Native spectroscopy and thermal enclosures](../README.md)

Expand the heat trace in hops around its diagonal part. Each closed path has a nonnegative weight from the hopping amplitude \(1/2\) and its holding times in \(D_0\). Write \(L_A=\int_0^\beta V_A(x_t)\,dt\) and \(L_B\) likewise. The path's connected weight is its baseline weight times \((e^{-L_A}-1)(e^{-L_B}-1)\). A visited neighbor of one pin sees either \(31/2\) or \(-17\) native energy units, with mean \(-3/4\). Both values interchange with the sign of the mobile. If \(k_A(x)\) counts adjacent A pins, a fixed path has mean \(L_A=-(3/4)\int k_A(x_t)\,dt\). Jensen's inequality gives \(\mathbb E e^{-L_A}\ge\exp((3/4)\int k_A(x_t)\,dt)\). Given the path and mobile charge, the two orientation choices are independent. Their averaged connected factor is nonnegative for every path and strictly positive on paths spending positive time in both neighborhoods. A walk to the other source and back has a positive measure of such holding times for every \(\beta>0\). Summing establishes \(\overline{\Delta Z}>0\) and strictly positive single-source increments for nonempty neighborhoods. Suppressing the off-diagonal exchange leaves only constant paths and gives \(\Delta Z=0\) identically on disjoint supports. The positive heat trace gives a negative dilute interaction energy. Comparing energies at distinct separations requires additional information about those separations.

## Bound independent of board volume

The blank one-mobile diagonal at degree \(k=5,6\) and charge \(q=\pm1\) is
\[
D_0=(8+1/2)k+8(6-k)((1+q)^2-1).
\]
It ranges from \(69/2\) to \(133/2\). A vertex adjoining both constituents of a single dipole receives their sum \(-3/2\). Thus \(\lVert V_A\rVert_\infty,\lVert V_B\rVert_\infty\le M=17\). Disjoint supports give \(D_0+uV_A+vV_B\ge b=35/2\) for all \(0\le u,v\le1\). Each source touches at most \(s=12\) available sites. The hopping row sum is at most \(\gamma=3\). These constants follow from local degrees and the native seam. They apply to every finite extension having these local conditions.

Put \(K_t=\exp[-t(D_0+uV_A+vV_B-G/2)]\). Differentiating the trace in both couplings and integrating over the unit square expresses \(\Delta Z\) as the integral of
\[
\beta \int_0^\beta\operatorname{tr}(V_A K_t V_B K_{\beta-t})\,dt .
\]
Entrywise, the positive path expansion bounds \(K_t(x,y)\) by \(e^{-bt}[\exp(tG/2)]_{xy}\). Each kernel between the two supports requires at least \(d\) hops. Anchor the trace at a site in A. For hop counts \(k,l\), summing the product over every intermediate site gives \((G^{k+l})_{xx}\le6^{k+l}\), which also bounds the sum restricted to B. Integration of \(t^k(\beta-t)^l/(k!\,l!)\) gives \(\beta^{k+l+1}/(k+l+1)!\). For \(m=k+l\) there are \(m-2d+1\) admissible splits. Taking absolute values gives the single-orientation bound
\[
 \frac{|\Delta Z|}{\beta}
 \le sM^2\,\beta e^{-b\beta}
 \sum_{m\ge 2d}\frac{(m-2d+1)(\gamma\beta)^m}{(m+1)!}.
\]
The same estimate holds after averaging. At fixed temperature the right-hand side is at most \(sM^2\beta e^{-(b-\gamma)\beta}(\gamma\beta)^{2d}/(2d)!\). Uniformly in temperature, use \(e^{-b\beta}(b\beta)^{m+1}/(m+1)!\le1\) in each summand and sum the geometric series:
\[
 \sup_{\beta>0}\frac{|\Delta Z|}{\beta}
 \le \frac{sM^2}{b(1-\gamma/b)^2}
       \left(\frac{6}{35}\right)^{2d}.
\]
The native energy zero is essential for the temperature-uniform bound. Shifting all mobile energies changes the physical activity associated with the zero-mobile arm. For sources confined to cells a cell-graph distance \(R\) apart, \(d\ge R-2\), since each neighborhood extends by at most one cell and a mobile hop crosses at most one seam. Any nonzero attractive tail proportional to \(1/R\) in cell units violates the bound at sufficiently large \(R\) at fixed activity. Applying that conclusion to proper distance would require its calibration from clocks and signals.

## Exact effect of the zero-or-one preparation

Let \(Z_0=\overline{Z_{00}}\), \(a=\overline{Z_{10}}-Z_0\), \(b_1=\overline{Z_{01}}-Z_0\), and \(c=\overline{\Delta Z}\). All three increments are positive. For the ensemble that equilibrates the internal states before taking the logarithm, the connected grand potential is \(\Omega=-\beta^{-1}\Delta\log(1+z\overline Z)\). The exact ratio multiplying \(-1/\beta\) inside the logarithm is
\[
 \frac{(1+z\overline Z_{11})(1+zZ_0)}
 {(1+z\overline Z_{10})(1+z\overline Z_{01})}
 =1+\frac{z[(1+zZ_0)c-zab_1]}
 {(1+z\overline Z_{10})(1+z\overline Z_{01})}.
\]
Since the denominator is at least \((1+zZ_0)^2\), the attractive part satisfies \(\max(0,-\Omega)\le zc/[\beta(1+zZ_0)]\le zc/\beta\). Thus the exponential range bound extends to every fixed \(z>0\) for this prescription. The term \(-zab_1\) records competition for the single mobile slot.

### The reference at infinite separation

The positive competition energy also decays on the specified ensemble. Invert the partition ratio. Its denominator is \((1+zZ_0)(1+z\overline Z_{11})\ge(1+zZ_0)^2\). Dropping the negative term in its numerator and using \(\log(1+x)\le x\) gives
\[
 \max(0,\Omega)\le \frac{z^2ab_1}{\beta(1+zZ_0)^2}.
\]
For either single-source increment, the trace derivative bound gives \(a,b_1\le \beta sM e^{-(b-\gamma)\beta}\). At fixed \(\beta\), the no-hop contribution \(Z_0\ge n e^{-133\beta/2}\), with \(n\) available sites, yields a positive energy at most a constant times \(n^{-2}\).

For fair mobile charge the bound can be made uniform over all temperatures. Each twelve-port cell has three external gateways and nine degree-five sites before pinning. Deleting the four pins leaves at least \(n/2\) degree-five mobile sites. Half the charge weight is a vacancy. A vacancy waiting at one such site contributes \(e^{-69\beta/2}\) to its heat trace. Consequently \(Z_0\ge (n/4)e^{-69\beta/2}\). Put \(x=(zn/4)e^{-69\beta/2}\); \(1+x\ge x^{2/5}\) for all \(x>0\). Then
\[
 \max(0,\Omega)\le (sM)^2 4^{4/5}z^{6/5} n^{-4/5}\,\beta e^{-7\beta/5}.
\]
The supremum of the last temperature factor is \(5/(7e)\). On the spherical midpoint family the number of available sites is \(n=240m^2-4\) at subdivision frequency \(m\). The greatest cell distance from a chosen corner is at most \(6m\). Hence every cell separation satisfies \(R\le12m\). Combining the two signs gives
\[
 \lim_{R\to\infty} R\,\sup_{\beta>0}|\Omega(R,z,\beta)|=0
\]
for each fixed positive activity on this family. This supplies the infinite-separation reference for the full connected energy. No additive constant permits a nonzero \(1/R\) coefficient. The comparison concerns the same four source arms and excludes the proposed Newton tail even if the temperature varies across the distance comparisons. [Native port incidences](../../../../../../executable_model/refinement.py) · [Midpoint volume and radius proof](../../../../../../experiments/interval_echo_einstein_station/collar_scaling/station_clock/README.md)

### Finite activity crossover

At \(\beta=1/8\), a constant path at any affected site gives \(a,b_1\ge e^{-133/16}(e^{3/32}-1)>3/160000\). The single-coupling derivative is bounded above in magnitude by \(\beta s M\) since \(\|K_\beta\|\le e^{-(b-\gamma)\beta}<1\). Thus \(a,b_1\le51/2\) and \(0<Z_0\le236\) on the measured board. In all four cases these bounds and the archived heat intervals imply \(ab_1>Z_0c\). There is exactly one sign change, at \(z=c/(ab_1-Z_0c)\). The table gives the connected energy signs determined by outward enclosures.

| Cell separation | Mobile hops | \(\Omega\) at \(z=10^{-72}\) | \(\Omega\) at \(z=10^{-35}\) |
|---:|---:|---:|---:|
| 2 | 4 | Negative | Negative |
| 3 | 7 | Negative | Negative |
| 4 | 10 | Negative | Positive |
| 5 | 12 | Negative | Positive |

These are signs relative to the connected zero of energy. A force direction would compare neighboring energies. The second activity was selected after inspecting the four archived heat traces to distinguish the two regimes. A many-mobile bath has a different partition function. Freely responding dipoles, proper distance, inertial calibration and other propagating modes require separate source evidence for the full Newton criterion. [Archived oriented traces](../unoriented.json) · [Endpoint and crossover bounds](results.json) · [Calculation](bounds.py) · [Independent graph and partition-function checks](test_bounds.py)

Run from the repository root with the model dependencies:

~~~~sh
python -m model.candidate_laws.interval_echo.openwave.row_26_newton.mobile_mediator.thermal_range.bounds --write
python -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_26_newton/mobile_mediator/thermal_range
~~~~
