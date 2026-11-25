# A conserved inverse-square profile on a programmed Interval Echo guide

A compact attractive contact has a positive stationary mode on each finite cubic guide. One contact setting, derived from the infinite-lattice binding threshold before measuring any radial profile, gives an inverse-square *occupation density* on increasing interior scales. The unadjusted contrast-eight two-body interaction binds much more tightly. Both statements concern a material excitation. Calling its measured energy an additional gravitational source would require independent evidence against counting that same energy twice.

## Generator, source and measurement

Use the [routed cubic two-load interaction](../../../model/experiments/seam_program/INTERACTING_LIMITS.md). At zero total momentum, in the exchange-even sector, the relative coordinate ranges over the sites of an even periodic cubic grid of side $n$, with the coincident site removed. If $L_*$ is the nearest-neighbor graph Laplacian after removing that site, and $P$ projects onto its six neighbors, the generator in native relative-exposure units, after subtracting a constant, is
$$K=L_*-2k P.$$
The active physical seam has contrast 8. The addressed interaction phase sets the effective $k$. At this stage a pair-separation measurement is relational. Each individual particle in a translation-invariant two-body state has uniform position marginal.

A distinct realization parks a binary source rail at the excluded address. On the other $n^3-1$ logical sites a single mobile load hops with coefficient $a=1/2$. Each degree-$d_x$ site carries the single-load contrast reference $8d_x$. Restore its one missing neighbor reference on the six sites beside the hole. A refocused seam against a vacant helper supplies $(8+1/2)t$ phase per occupation. Duration $8\delta/(8+1/2)$ realizes the missing $8\delta$ in a slice of duration $\delta$. This makes the background generator $48I+aL_*$. Between the parked rail and each neighboring site the verified native occupation word has action $\exp(i17\tau n_{\rm source}n_x)$. Take $\tau=2ak\delta/17$. The logical Hamiltonian is then
$$H=48I+a(L_*-2k n_{\rm source}P).$$
An empty source gives a uniform ground mode on the punctured grid. An occupied source gives the mode below. The state can therefore be localized relative to a physical parked load in a single laboratory. The source's position, the cubic routing, its energy reference and the counterpulse are supplied. The derivation applies to the generated Hamiltonian in the product-formula limit. Six disjoint hopping matchings on an even cube, a diagonal compensation layer and the source phase suffice. At most one mobile load contributes, so each layer's operator norm is bounded independently of side. The first-order splitting error over logical time $T$ is $O(T\delta)$ by the pairwise commutator bound. Elapsed proper time and physical Euclidean ranges of the guide require their own calibration. [Native occupation word](../../../model/experiments/seam_program/range_contact.py) · [Full two-load reduction checked against the relative generator](../../../model/experiments/seam_program/test_kinetic_resonance.py)

The measurement replays the occupied and sham phase words at the threshold setting on the physical 20-cell carrier, giving maximum four-column errors below $4.8\times10^{-12}$. The spatial grids in the radial table are logical-generator calculations. A spectral filter or a prepared eigenstate selects the enrolled state. Time evolution alone need not prepare its ground state.

## Exact stationary state on each finite grid

For the *unpunctured* periodic Laplacian $L$ let $g_m=(L+m^2 I)^{-1}|0\rangle$ and denote the origin and any nearest-neighbor value by $g_0,g_1$. Here $m>0$ is a dimensionless inverse length parameter. Positivity follows from the positive walk expansion of the resolvent. The origin equation is $(6+m^2)g_0-6g_1=1$. At a neighbor of the missing site the deleted contribution to the Laplacian is $g_1-g_0$. Consequently the restriction of $g_m$ to all nonzero sites, normalized in the ordinary squared norm, is an eigenvector of $K$ of eigenvalue $-m^2$ exactly when
$$g_0=(1+2k)g_1.$$
There is one positive root for every $k>0$. For the continuous-time walk with generator $-L$, $g_1/g_0=\mathbb E_e[\exp(-m^2 T_0)]$ by the first-hit decomposition of the resolvent. The finite torus is connected, so this ratio decreases strictly from one to zero. The resulting positive eigenvector is the unique ground state by Perron–Frobenius on the connected punctured graph. Its reality gives zero stationary bond current. The full sparse matrix, independent diagonalization and direct unitary propagation verify the identities on a small grid. The disabled-source control evolves the same state into a different position distribution. [Implementation](measure.py) · [Independent checks](test_mode.py)

## One coupling for increasing domains

Let $G$ be the infinite cubic Green function at zero energy in the Laplacian convention above. Then $G(0)=0.25273100986\ldots$, $G(e)=G(0)-1/6$ and $G(x)\sim1/(4\pi |x|)$. The last coefficient follows from the quadratic minimum of the Fourier symbol $2\sum_j(1-\cos p_j)$; subtract a smoothly cut off $1/|p|^2$ and integrate its remainder by parts. The heat-kernel integral for $G(0)$ is implemented in the cubic program. Fix
$$k_*={1\over 12[G(0)-1/6]}=0.96826781269\ldots$$
on every grid. This selects the infinite-volume binding threshold through local spectroscopy. It uses no velocity, lensing or halo-profile data. [Derivation of the two-body threshold](../../../model/experiments/seam_program/INTERACTING_LIMITS.md#finding-the-threshold-from-the-native-operator)

For fixed $m>0$ the periodic resolvent converges to the infinite one as $n\to\infty$ by ordinary Fourier sums. In infinite volume $g_1/g_0$ at $m=0$ equals $1/(1+2k_*)$ and decreases strictly for $m>0$. The finite roots $m_n$ must therefore tend to zero. At those roots the origin identity gives $g_0=[6+m_n^2-6/(1+2k_*)]^{-1}\to G(0)$. For each fixed lattice vector $x$, subtracting $g_m(x)$ from $g_m(0)$ inserts a factor $1-\cos(p\cdot x)$ in the Fourier sum. Its ratio with $m^2+2\sum_j(1-\cos p_j)$ is bounded uniformly in $m$ for that $x$. Removing a small neighborhood of zero and then restoring it justifies the Riemann limit. Thus $g_{m_n}(x)\to G(x)$ pointwise. Writing $p_n$ for the normalized ground-state probability, for any lattice axis $e$,
$$\lim_{r\to\infty}\lim_{n\to\infty}{p_n(re)\over p_n(2re)}=4.$$
The argument also works on fixed rays of other lattice directions. It proves the iterated limit and permits a sufficiently slowly growing diagonal sequence. The sampled growing-radius sequence below is a separate check. Since $G$ is not square summable, the normalization diverges, so one normalized packet's probability at a fixed site tends to zero. A finite cutoff or a specified increasing total occupation is necessary for a nonzero halo normalization on an increasing domain.

| Guide side | Axial distances | Density ratio at $k_*$ | $n m_n$ |
| ---: | ---: | ---: | ---: |
| 64 | 4 and 8 | 3.96557 | 1.95593 |
| 128 | 5 and 10 | 4.00872 | 1.95078 |
| 256 | 8 and 16 | 3.98296 | 1.94826 |
| 512 | 11 and 22 | 3.99264 | 1.94701 |

At effective $k=8$ the side-64 ratio at 4 and 8 is $3.87\times10^9$. The inverse length parameter tends to a positive value (about 3.396) because the infinite-volume bound state lies below threshold. Positive mass in the resolvent gives exponential spatial decay. This is the unadjusted two-body contact on the logical guide. The parked reduction with that same effective $k$ uses its declared parked-source realization.

## Which energy is conserved?

For any mobile wave put $p_x=|\psi_x|^2$, $\phi=H\psi$ and $e_x=\operatorname{Re}(\bar\psi_x\phi_x)$. The latter is the expectation of the symmetrized operator $(|x\rangle\langle x|H+H|x\rangle\langle x|)/2$. Number and this site energy satisfy exact lattice continuity, with outward currents
$$J_{xy}=-2\operatorname{Im}(\bar\psi_x H_{xy}\psi_y),\qquad
J^E_{xy}=-\operatorname{Im}(\bar\psi_x H_{xy}\phi_y+\bar\phi_x H_{xy}\psi_y).$$
They are antisymmetric by Hermiticity. Differentiating the densities under $i\dot\psi=H\psi$ gives the negative sum of the corresponding outward currents. The energy calculation uses $i\dot\phi=H\phi$ as well. Tests exercise these identities on complex moving packets with nonzero currents. In the stationary mode $e_x=(48-a m_n^2)p_x$. With the parked reference the conserved energy therefore has the same radial profile. A joint reading of the sharp total eigenenergy followed by position realizes these statistics on that mode.

An independent local energy allocation puts half of the measured bond gradient energy $a|\psi_x-\psi_y|^2$ at each endpoint and assigns the on-site potential energy there. It differs from $e_x$ by exactly $-(a/2)(L_*p)_x$. The two regional readings differ by boundary flux. Away from the source in the iterated limit $G(x+e)/G(x)\to1$ at large $|x|$, so the bond-gradient part is $o(p_x)$ relative to the retained rest reference. Changing the reference by a multiple of total load changes the local coefficient by the same multiple. The number equation and the Hamiltonian energy equation do not fix the physical gravitational zero of energy, pressure, momentum flux in a Lorentz frame, or the response of proper geometry. The two-arm optical comparison can also change upon dephasing at fixed population and energy. The [material stress tests](../collar_identifiability/refinement/composition/README.md) distinguish these preparations.

From the repository root:

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code .venv/bin/python -m matched_einstein_limit.coherent_halo.measure
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code .venv/bin/python -m pytest -q code/matched_einstein_limit/coherent_halo
~~~

## A normalized halo on a fixed spatial volume

Keep the same contact \(k_*\) and let the spacing of the programmed guide be \(1/n\) on the unit periodic cube. The normalized probability measures of the parked-source ground states converge in total variation after spreading each site's mass uniformly over its voxel. The limiting density is \(F_\mu^2/\|F_\mu\|_2^2\), where \(F_\mu\) is the periodic Green function of \(-\Delta+\mu^2\) and \(\mu=\lim n m_n\) is positive. No new contact or mass is fitted. Here is a derivation.

Write \(\lambda(\theta)=4\sum_i\sin^2(\theta_i/2)\) and \(G=\int_{[-1/2,1/2]^3}1/\lambda(2\pi t)\,dt\). Use the centered reciprocal cells \(C_j=j+[-1/2,1/2]^3\) on the torus of side \(n\), and put \(a_n(x)=1/[n^2\lambda(2\pi x/n)]\). The massless, zero-mode-deleted sum differs from its integral, in this rescaled normalization, by
\[
 B_n=\sum_{j\ne0}\left[a_n(j)-\int_{C_j}a_n(x)\,dx\right]-\int_{C_0}a_n(x)\,dx.
\]
It has a finite limit \(B\) obtained by replacing \(a_n(x)\) by \(a(x)=(4\pi^2|x|^2)^{-1}\) and summing over \(\mathbb Z^3\). Indeed \(\lambda(2\pi x/n)\) is bounded below by a constant times \(|x|^2/n^2\) on the centered torus. The first and second derivatives of its rescaled denominator \(n^2\lambda\) are \(O(|x|)\) and \(O(1)\). The Hessian of \(a_n\) is therefore \(O(|x|^{-4})\) away from zero, also across the periodic faces using the nearest image. The midpoint linear error integrates to zero. The cell errors away from the origin have a summable \(C/(1+|j|)^4\) majorant. The central singularity is bounded by \(C/|x|^2\), which is integrable.

For each fixed \(v>0\), separating the zero frequency gives
\[
 \lim_n n[g_{v/n}(0)-G]
 = B+v^{-2}+\sum_{j\ne0}\left[
 (4\pi^2|j|^2+v^2)^{-1}-(4\pi^2|j|^2)^{-1}\right].
\]
The differences are dominated by \(C_v/|j|^4\). This continuous function is strictly decreasing from positive to negative infinity. For the negative limit sum the terms with \(1\le |j|\le v/(2\pi)\), whose absolute values are at least \(a(j)/2\). Denote its unique zero by \(\mu\). At the scattering pole the exact binding condition is \(g_m(0)=1/(G^{-1}+m^2)\). Replacing \(G\) by its right side at \(m=v/n\) changes the scaled mismatch by \(O(v^2/n)\). The opposite signs at any two trial values bracketing \(\mu\), continuity and the previously proved unique finite root give \(nm_n\to\mu\).

Before removing the coincident site, the step field with values \(n g_{m_n}(x)\) has discrete Fourier coefficients \(1/[n^2\lambda(2\pi j/n)+(n m_n)^2]\). Their squares have a summable uniform bound \(C/(1+|j|)^4\). On each fixed frequency the voxel interpolation of the exponential converges in \(L^2\) to the continuum exponential. Parseval controls the remaining frequencies uniformly because the step images of the grid characters are orthonormal. The fields converge strongly in \(L^2\) to \(F_\mu\). Dropping the origin costs \(g_{m_n}(0)^2/n\) in squared norm. Normalizing and applying Cauchy–Schwarz to the difference of squares proves the density limit.

Periodizing \(e^{-\mu |x|}/(4\pi |x|)\) gives \(F_\mu\), since this sum solves the periodic massive Green equation and has the Fourier coefficients above. Near the anchor \(F_\mu=1/(4\pi r)+O(1)\). The limiting occupation in a ball of radius \(r\) is \(r/(4\pi\|F_\mu\|_2^2)+O(r^2)\). This fixes the inner slope relative to the total occupation of this one mode and the box size. In an increasing physical box that coefficient changes. The proof describes the interval between a localized source and the screening or box scale. A different physical geometry or a galaxy's baryonic normalization needs independent evidence. [Fixed-volume numerical check](fixed_box.json) · [Discrete Fourier and puncture check](fixed_box.py)

## An independent ball ruler on the same ordinary source

Use the parked guide with the anchor vacant. With degree restoration its one-mobile Hamiltonian is \(48 I+L_*/2\), so the uniform mobile mode has energy 48. An unobserved lottery of strength \(t\) selects this mode or an empty guide with the load parked outside the measured region. The measured local number and symmetrized energy both vary. A second group of rails runs the same outward clock pulses in both preparations. All marks, exits, radii, areas and volumes read by these clock rails are identical under the lottery. This is an explicit separable, programmed clock control. It supplies the chart and the conversion between radius and lattice steps.

In a ball away from the puncture, take the two first exterior distances \(d_\pm\) along one axis and weight occupied positions by \(2\pi d_+d_-/(2\ell)\). The calibrated material switch has gain \(C(0)=0.0645423913551\) nats. Use switch scale \(\Gamma=48/C(0)\) to read the established energy per mobile. The local switch probability is bounded by \(\Gamma\pi\ell(1+h/\ell)^2<1\) on these grids. For the mesh \(h=\ell/N,\ \ell=1/(2048N)\), let \(J_N\) be the sum of the measured weights times voxel volume divided by \(\ell^4\). Exact row sums compute it without allocating a macroscopic state. Both exits have error at most \(h\); Riemann convergence on the unit ball gives \(J_N\to8\pi^2/15\). At \(t=1/2\) the actual terminal entropy slope divided by \(\ell^4\) is \(C(t\Gamma\ell^4 J_N/(1-h^3))\Gamma J_N/(1-h^3)\). It tends to \(48(8\pi^2/15)=252.661873\). The area derivative vanishes at every resolution. At 128 steps per radius the two normalized derivatives are 256.29746 and zero. No finite positive area coefficient cancels them. The signed energy calibration and independent ball definition can therefore both succeed on a chart whose geometry is insensitive to matter. All energies here use the guide's specified exposure unit. [Matched-source control and data](stationarity.json) · [Clock sum and native calibration](stationarity.py)
