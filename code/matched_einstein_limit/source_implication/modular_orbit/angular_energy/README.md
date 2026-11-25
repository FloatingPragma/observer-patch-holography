# An energy representation of twelve regional entropy responses

On one serviced Interval Echo cell, the angular part of a regional entropy variation can be measured with native seam energies. This holds simultaneously on all twelve oriented stars for **every perturbation of the constant mode and first spatial triplet**. A monotone ball kernel on the same seams gives the wrong sign for an equal-energy pair of preparations. These stars have a common reference prepared by a finite service. They supply a finite test of a proposed physical matter input. The Einstein equation also requires an independently measured geometric response on the same proper regions. [Angular source and detectors](../../../../../model/experiments/interval_echo_angular/README.md)

## Source and identity

Use the registered reference on the zero-, one- and two-load words of one cell, prepared by four palindromic sweeps of 30 native seam filters and averaged over all 60 proper rotations. Its acceptance is 0.0126583. Each star contains a pole and its five neighbors. Let \(K_n\) be its marginal negative logarithm, lifted to the parent algebra after subtracting the vacuum reading. For a traceless normalized-state tangent \(D\), \(\operatorname{Tr}(D K_n)\) is the exact entropy derivative of that star at the common reference. The reference is faithful on the full word space. Every sufficiently small tangent supported on the one-load slot is therefore admissible. Its perturbed one-load density can be realized by a lottery of native pure one-load preparations. Number-sector probabilities and the mixtures on the two other sectors can be retained by conditioning the original input lottery on its load count. [Source filters and regional trace](../../../../../model/experiments/interval_echo_einstein_station/cap.py) · [Native one-load preparation](../../../../../model/experiments/interval_echo_space/shell_source.py)

Let \(P\) project onto the constant one-load mode and the three first-band modes. This sector is invariant while the internal seams evolve with common strength and the external gates stay closed. Let \(H_s(n)\) be the sum of the five native seam energies incident on the pole, \(H_r(n)\) the sum on the five edges around its rim and \(H\) the sum on all 30 internal seams. On one load each edge has energies 8 and 9 on its symmetric and antisymmetric endpoint modes and zero on a missed endpoint. The following operator equality holds with the same four real constants at all poles,
\[
 P K_nP=P\big[cI+\gamma H+uH_s(n)+vH_r(n)\big]P.
\]

| Coefficient | Certified interval, rounded outwards |
|---|---:|
| \(u\), spoke weight | \([0.3528196218134,0.3528196218137]\) |
| \(v\), rim weight | \([0.6979178648817,0.6979178648821]\) |
| \(\gamma\), full-cell energy weight | \([-0.0491057198812,-0.0491057198795]\) |
| \(c\), common offset | \([1.99389651448,1.99389651456]\) |

At fixed total cell energy the first two local terms give the entire entropy variation. For an arbitrary one-load perturbation in the protected sector, subtract the average over the twelve pole readings. Each seam occurs in two spoke sets and two rim sets. The centered regional variation is therefore exactly the centered variation of \(uH_s+vH_r\). This supplies all three dipole and five traceless quadrupole coefficients even for perturbations that change total energy. The calibration uses the modular coefficients at one pole. Direct replay on all poles gives a maximum compressed operator error of \(5.83\times10^{-14}\) nats. Repeating the identical energy prediction on a reference without the rotation average gives a discrepancy of \(0.02392\) nats. [Replay and full count data](results.json)

## Proof and interval certificate

Use unit port directions \(\mathbf n_i\) from the first triplet of the cell. The columns of \(B=[\mathbf1/\sqrt{12},\mathbf n/2]\) are orthonormal and \(\mathbf n_i\cdot\mathbf n_j=q=1/\sqrt5\) across a seam. At pole \(n\), the real symmetric compressed operators \(K_n,H_s(n),H_r(n)\) commute with the fivefold stabilizer. The scalar–vector block is consequently \(b n^\top\), and the triplet block is \(c_T I+d\,nn^\top\). If \(a\) is the scalar entry, the nontrivial monopole on normalized differences is \(m=a-c_T-d/3\). Thus \((m,b,d)\) determines such an operator up to an identity. Rotations carry the same coefficients to the other eleven poles.

For an individual seam with endpoint projections \(x,y\) along the pole, write \(A_\parallel=[17(x^2+y^2)/2-xy]/4\) and \(\tau=(17-q)/4\). Its axial coefficients read
\[
 a_e=4/3,\quad m_e=(q-1)/12,\quad b_e=8(x+y)/\sqrt{48},
 \quad d_e=(3A_\parallel-\tau)/2 .
\]
They follow from \(h_e=(17/2)(|i\rangle\langle i|+|j\rangle\langle j|)-(|i\rangle\langle j|+|j\rangle\langle i|)/2\). The five edge types have \((x,y)=(1,q),(q,q),(q,-q),(-q,-q),(-1,-q)\), up to order. Each spoke and each rim consists of five edges of its type. Solve the two dipole and quadrupole equations using these two types. The resulting positive weights are \(u,v\) in the table. The full-cell energy has no angular coefficient by symmetry and has monopole \(30m_e\ne0\). Its coefficient and the identity then fix the remaining two scalar equations.

For the source enclosure, the filter amplitudes in its two one-load seam modes are \(2^{-8/32}\) and \(2^{-9/32}\) per local step. The entries in the port basis are nonnegative, so rounding up or down in every product bounds each matrix entry. Use denominator \(2^{75}\), propagate one sweep \(F\), form \((F^\top F)^8\), sum the rotations and trace out the other six ports. Positive common normalization cancels in the vacuum-centered log. The unnormalized star single-particle matrix comprises a pole, its uniform ring mode and two real Fourier doublets. For pole diagonal \(p\), pole–ring entry \(v_0\), ring diagonal \(z\), ring neighbor \(w\) and second neighbor \(y\), the first block is
\[
 \begin{pmatrix}p&\sqrt5 v_0\\\sqrt5 v_0&z+2w+2y\end{pmatrix}.
\]
The transverse triplets have half their norm squared in the first ring doublet with eigenvalue \(z+(\sqrt5-1)w/2-(\sqrt5+1)y/2\). The constant and axial triplets project into the displayed block as \((1,\sqrt5)/\sqrt{12}\) and \((1,1)/2\). Scalar interval logarithms then give \(m_K\simeq-0.17415148,\ b_K\simeq 6.55201374,\ d_K\simeq-0.23106153\) and certify the coefficient signs. [Integer and interval calculation](certificate.py)

In a changed law that turns off squared contrast, every exchange seam annihilates the constant one-load mode. All radial energy combinations on that law have exactly zero scalar–triplet dipole. The same local preparation filter is implemented in this control. The contrast shift of eight on the always singly occupied filter ancilla contributes a common phase per pulse, and its phase at the \(\pi/2\) read-basis gate is one. Its normalized success operation and hence the prepared reference stay equal. The measured nonzero \(b_K\) distinguishes the two couplings. [Control and independent propagation](measure.py)

## A reference stationary during native flight

The same construction works for a reference stationary under the complete internal cell Hamiltonian on **all 79 words**. After the filters and the rotation lottery, leave the internal seams open. For each pair of distinct energy levels in the same load sector, let the apparatus wait \(\pi/|E_i-E_j|\) on a fair bit and otherwise proceed immediately. The off-diagonal \(ij\) block is multiplied by \((1+\exp[-i(E_i-E_j)t])/2=0\) at that appointment. The waiting channels commute, preserve degeneracies and leave diagonal populations unchanged. The parent has no coherences between load sectors. The complete finite procedure therefore produces \(\bar\rho=\sum_E P_E\rho P_E\) exactly in the ideal timing model. After preparation, uninterrupted internal cell evolution leaves this full reference fixed. The original populations, closed external gates, rotation lottery and timing service are construction inputs.

Isolating the roots of the square-free characteristic polynomials of \(2H\) gives four distinct one-load energies and eighteen distinct two-load energies. Six plus 153 independent fair appointments suffice with a combined maximum dwell below 2281.78 native time units. All durations are explicit computable functions of the native integer Hamiltonian. For an appointment with timing error \(\delta t\), the trace-norm error of the fair waiting channel on a normalized input is at most \(\delta t\,\inf_c\|H-cI\|\), using Duhamel's formula on the unitary branch. Here \(\delta t\) denotes the absolute error. Contractivity allows these errors to be added over the appointments.

There is also an interval proof that the logarithmic signs survive this change of reference. Let \(L\le R\le U\) be the entrywise integer enclosures of the unnormalized real rotation sum and let \(\mathcal E_H\) be its energy dephasing. In each load sector, the symmetric rotation-invariant matrices have disjoint integer indicator matrices as an orthogonal basis. Their coefficients obey integer linear equations for commuting with \(2H\). Rational elimination gives the exact squared Hilbert–Schmidt distance from \(L\) to this commuting subspace. The energy dephasing is the Hilbert–Schmidt projection onto the full commutant. Thus
\[
 \|R-\mathcal E_H R\|_1
 \le\sqrt{79}\,[\|U-L\|_{\rm HS}+
       \mathrm{dist}_{\rm HS}(L,\{X:[X,H]=0,\ X\ {\rm rotation\ invariant}\})].
\]
The partial trace is trace-norm contractive. For either positive regional block, a lower eigenvalue bound \(a>e\) and an operator perturbation of size \(e\) give \(\|\delta\log\|\le e/(a-e)\) by the resolvent integral. Including the vacuum subtraction bounds the operator change of every single-particle star read by \(4.285\times10^{-5}\) nats. This bound uses the original source and exact native commutator equations. It is independent of numerical eigenvalue-gap estimates. [Commuting-matrix certificate and finite waiting protocol](stationary.py)

At the dephased reference the additional rim coefficient lies in \([0.3450764,0.3451201]\), the quadrupole in \([-0.231148,-0.230975]\), and the separating tangent still gives at least \(0.0677743\) nats. A separate simulated angular read gives \(-0.23055\pm0.01016\) nats for axial minus transverse. The replayed source has commutator norm \(8.77\times10^{-21}\) in native units. Changing from a prepared spatially symmetric reference to one that the internal Hamiltonian also preserves therefore leaves the angular test intact. Time stationarity concerns evolution under that Hamiltonian. Fixed-proper-volume entropy stationarity concerns a different family of source variations. [Computed receipts and separate stationary-source counts](results.json)


## Which bulk kernels the records admit

For a positive interior ball weight that decreases towards an outer boundary, the more central spoke carries at least as much weight as a rim seam. The Euclidean midpoint distances from the pole on the slow-mode chord chart have these two squared values, \((1-q)/2\) and \(3(1-q)/2\). For a fixed pole prepare the axial triplet and an equal mixture of the two transverse triplets. Their total energies are identical because the first band is degenerate. This contrast reads \(d_K\) directly. The simulated records below use \(2^{20}\) independent preparations per setting and arm. Standard errors condition on the specified reference.

| Axial minus transverse | Exact-source prediction, rounded | Simulated difference ± one standard error |
|---|---:|---:|
| Regional entropy derivative, nats | −0.231062 | −0.2403 ± 0.0102 |
| Five spoke energies | +7.940983 | +7.9304 ± 0.0099 |
| Five rim energies | −4.345492 | −4.3534 ± 0.0090 |
| Total energy | 0 | Exact by degeneracy |

Any nonnegative inside-star weights \(u'\ge v'\) give a nonnegative contrast \(7.940983\,u'-4.345492\,v'\), of the opposite sign from the regional entropy derivative. This result applies to the two midpoint classes on this chart. A boundary stress, different regional embedding or different choice of reference can alter the bulk comparison.

The finite operator split gives a quantitative boundary correction if the spoke normalization is held fixed. Take a constant interior weight \(u\); the exact additional rim term is \((v-u)H_r\), with \(v-u\in[0.3450982430680,0.3450982430687]\) nats per native energy unit. An inward-increasing candidate with the same spoke weight has a bare rim weight at most \(u\), so fitting the two angular coefficients requires at least this much extra rim weight. The independence of the two angular responses makes this coefficient unique in the two-class basis. Applying this finite, prepared-source correction to a gravitational stress would require the physical source-to-curvature comparison.

There is also a separation for **arbitrary nonnegative weights on all 30 seams**, without any hypothesis about radial symmetry. Choose \(\lambda,\mu\) so the linear functional \(L=m-\lambda b-\mu d\) vanishes on each spoke and rim. The exact axial formulas give \(L<0\) on the other three edge classes, while interval evaluation gives \(0.0678624358315<L(K_n)<0.0678624358338\) nats. For a physical tangent in the basis of the constant mode and triplet set
\[
 D_{00}=1,\quad D_{0i}=-\lambda n_i/2,\quad
 D_{ij}=-\delta_{ij}/3-\mu(3n_i n_j-\delta_{ij})/2 .
\]
It is traceless, and \(\operatorname{Tr}(D O)=L(O)\) for every individual seam and the regional charge. Thus any nonnegative seam kernel misses this response by at least \(0.0678624358315\) nats per unit tangent. A scalar added to the kernel has zero response. The common-energy term accounts for the sign in the measured operator identity.

These coefficients concern a finite regional modular reading and an independently measurable set of local energies. Their normalization was inferred from the reference. The twelve physical directions are the cell's measured directions. A claim about a continuum of physical null rays would require identifying their boost normalization, a proper regional metric, and the curvature response on the same scaling family. [Matched regional theorem](../../../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

Run `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.modular_orbit.angular_energy.measure --write` from the repository root.
