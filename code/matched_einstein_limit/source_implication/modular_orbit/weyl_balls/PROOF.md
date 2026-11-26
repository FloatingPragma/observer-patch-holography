# A common finite Weyl reference for the ball entropy and stress readings

Consider a free, massless, two-component fermion at unit signal speed. Give it a stationary, spatially homogeneous, charge-symmetric excitation with bounded momentum support. At arbitrarily small resolutions one can read its regional entropy and its independently specified hopping energy on a finite menu of shrinking, overlapping balls in finitely many moving frames. For the constructions below,

\[
 \Delta S(B_R,u)=\frac{8\pi^2}{15}R^4\,T_{uu}+o(R^4),
 \qquad D(\rho_{B_R,u}\Vert\sigma_{B_R,u})=O(R^7).
\]

Here \(\sigma\) and \(\rho\) are the actual vacuum and excited **observer marginals** of one machine at a given resolution. The excitation has fixed amplitude as \(R\) shrinks. The proof concerns the matter input for flat tangent diamonds. The area of an observer's ball in a responding geometry is a further observable.

## Signal and measuring instrument

On an antiperiodic cubic torus of side \(L=Na\) with even \(N\), use the number-conserving nearest-neighbour Weyl hopping. Its one-particle matrix is

\[
 H_a(p)=\boldsymbol{\sigma}\cdot d_a(p),\qquad
 d_{a,j}(p)=\sin(ap_j)/a,\qquad p_j=2\pi(m_j+1/2)/L .
\]

Let \(P_\pm(p)=(I\pm\hat d_a\cdot\boldsymbol{\sigma})/2\). The signal reference fills precisely \(P_-\), a single stationary ground state shared across all regions. The small-momentum group speed fixes the clock coefficient of \(H_a\). The other lattice cones lie outside the occupied perturbation window as \(a\) decreases. The theorem tests smooth long-wavelength records.

There is a finite preparation from any input at precisely half filling. First read the commuting hopping-mode occupations. Each outcome has \(N^3\) particles among \(2N^3\) modes. Antiperiodicity on an even grid excludes zero energy, so exactly \(N^3\) modes have negative energy. Choose an unordered pair of distinct modes uniformly and exchange their occupations whenever the independent hopping work reading is strictly downhill. The reservoir accepts that work. With \(m\) occupied positive modes there are exactly \(m\) negative holes. Precisely \(m^2\) pairs reduce \(m\) by one and every other accepted move keeps \(m\) fixed. The probability per attempt is \(m^2/\binom{2N^3}{2}\), regardless of the occupations inside either band. Thus the expected number of attempts to reach the sea is \(\binom{2N^3}{2}\sum_{j=1}^{m}j^{-2}\); absorption occurs almost surely. Reading that all negative modes are filled certifies the unique ground-state occupation. This preparation uses global hopping-mode readout, pair exchanges and a work reservoir. Interval Echo would have to supply the corresponding operations on any proposed native field sector. [Native critical exchange and its selection control](../../../../../model/experiments/interval_echo_einstein_station/critical_exchange/selection/README.md)

Fix a compactly supported occupation \(0\le f(p)\le1\). In the excited signal occupy a fraction \(f(p)\) of positive-energy modes and empty a fraction \(f(-p)\) of negative-energy modes. Both global preparations commute with the same hopping. Their covariance difference is

\[
 F_a(p)=f(p)P_+(p)-f(-p)P_-(p),\qquad F_a(-p)=-F_a(p).
\]

For both preparations publish the output modes \(b_x=\sqrt q\,a_x+\sqrt{1-q}\,e_x\), with \(e_x\) independent tracial fermion modes and \(q=1-2\delta\). This passive, sitewise interferometer has a unitary completion. It acts once in the common protocol and commutes with a change of mode basis applied to both copies. The two output covariances are \(C_0=\delta I+qP_-\) and \(C_1=C_0+qF_a\). They are stationary and satisfy \(\delta I\le C_i\le (1-\delta)I\). Thus the positive occupation floor is a property of the acquired signals before restriction. The ideal ground and the excitation occupy the signal port independently of detector resolution.

For a regional mode inclusion \(J\), observers obtain \(C_{i,B}=J^\dagger C_iJ\) by recording pair correlations, accessible by occupation counts after two-mode interferometry. The hopping stress follows from the separately specified action. In the continuum,
\[
 H=-i\boldsymbol{\sigma}\cdot\nabla,\qquad
 \rho=T^{00}=\int \operatorname{tr}[H(p)F(p)]\,\frac{d^3p}{(2\pi)^3}.
\]
The same finite input gives the independently measurable energy difference from the expectation of each nearest-neighbour hopping. The output readout multiplies this difference by \(q\). The limit has \(q\to1\).

## The continuum operator and the endpoints of the logarithm

Let \(C_B=1_B P_-1_B\) on \(L^2(B,\mathbb C^2)\), \(h(t)=t\log t+(1-t)\log(1-t)\), and \(k(t)=-h'(t)\). The gauge-invariant CAR modular generator on regional modes is \(k(C_B)\), up to the sign convention for field creation versus annihilation. With the filled-negative-energy convention and positive energy cost for a particle, the vacuum modular Hamiltonian of a ball of radius \(R\) acts on smooth interior spinors as
\[
 k(C_B)=2\pi\,\{w_R,H\}/2,\qquad
 w_R(x)=(R^2-|x|^2)/(2R).
\]
The conformal transformation of the vacuum double cone proves the independently normalized \(2\pi\). For a Weyl field the differentiated flow on its Cauchy data gives precisely the symmetrized product including the spinor derivative of the weight. Translation, dilation and the Lorentz representation of the same field give the other balls and frames. [La Piana and Morsella, Theorem 3.4 and equation (3.10)](https://doi.org/10.1007/s00220-025-05253-1) · [Casini, Huerta and Myers, section 2](https://arxiv.org/abs/1102.0440)

Dilate each rest ball to \(B_1\). The continuum vacuum covariance is invariant under this unitary rescaling. The rescaled covariance difference \(R^{-4}\Delta C_{B_R}\) has kernel
\[
 \frac{q}{R}\int e^{iRp\cdot(x-y)} F(p)\,\frac{d^3p}{(2\pi)^3}.
\]
Oddness cancels the constant and quadratic terms. Its limit \(A\) has kernel \(i\int p\cdot(x-y)F(p)\,d^3p/(2\pi)^3\). It has finite rank, spanned by constant and linear spinors restricted to \(B_1\). The remainder is \(O(R^2+\delta)\) **in trace norm**. To see this, factor each exponential kernel into the rank-one position vectors \(e^{iRp\cdot x}\) and \(e^{iRp\cdot y}\), for each of two spin components. Three derivatives in \(R\) give rank-one products of polynomials times the same vectors. Their trace norms are bounded by a constant times \(|p|^3\) on the unit ball. Integrate the Taylor remainder over the fixed bounded momentum support.

The true measured modular matrix is \(k_\delta(C_B)=k(\delta I+qC_B)\). Its norm is bounded by \(\log((1-\delta)/\delta)\). This controls the trace-norm remainder after taking \(\delta\) proportional to the largest shrinking radius in fixed inverse-length units. For the leading finite-rank term, all affine spinors \(v\) lie in the domain of \(k(C_B)\). Cut \(v\) off at depth \(\epsilon\) from the boundary. The cutoff derivative costs \(O(1/\epsilon)\) on the strip and \(w_1=O(\epsilon)\) there. Both the omitted spinor and its weighted derivative have \(L^2\) norm \(O(\sqrt\epsilon)\). The graph closure of the interior action therefore contains \(v\). Spectrally \(|k_\delta(t)|\le |k(t)|\) and \(k_\delta(t)\to k(t)\) for \(0<t<1\). Functional calculus and dominated convergence now give \(k_\delta(C_B)v\to k(C_B)v\) on each of the finitely many vectors in \(A\). This step is necessary for removing the detector floor.

For the same reason smooth plane waves restricted to the ball are in the modular domain. Evaluate \(\operatorname{Tr}[k(C_B)A]\) by acting on the first variable of its kernel and taking the diagonal. The derivative of \(w_1\) multiplies \(A(x,x)=0\); the remaining derivative gives \(\operatorname{tr}[H(p)F(p)]w_1(x)\). Since \(\int_{B_1}w_1\,d^3x=4\pi/15\),
\[
 R^{-4}\Delta K_{B_R}
 =R^{-4}\operatorname{Tr}[k_\delta(C_B)\Delta C_{B_R}]
 \longrightarrow \frac{8\pi^2}{15}\rho .
\]

## Relative entropy at finite excitation

For finite gauge-invariant quasifree states, even when their regional covariance matrices fail to commute,
\[
 D(\rho(C_1)\Vert\rho(C_0))
 =\operatorname{Tr}\big[h(C_1)-h(C_0)-h'(C_0)(C_1-C_0)\big]
 =\Delta K-\Delta S .
\]
Diagonalizing each covariance proves the formula from its quadratic Fock-space logarithm. The Hessian of \(\operatorname{Tr}h\) in an eigenbasis weights each squared matrix entry by the divided difference of \(h'\), including \(h''\) on coincident eigenvalues. All these weights are at most \(1/[\delta(1-\delta)]\) along the entire segment. Consequently
\[
 0\le D\le\frac{\|\Delta C_B\|_{\mathrm{HS}}^2}{2\delta(1-\delta)} .
\]
Oddness gives zero coincidence kernel. Its first derivative is uniformly bounded by the first momentum moment of the source. In lattice units each entry at \(x,y\in B_R\) is bounded by \(C a^3|x-y|\) by subtracting the zero-separation Fourier sum. For \(a\le R\), the sum of the squared entries is at most \(C'R^8\). The continuum double integral obeys the same bound. Fix an inverse length \(\kappa>0\) and take \(\delta=\kappa R_{\max}<1/2\) for the whole menu. Then \(R^8/\delta\le R^7/\kappa\) for every smaller radius. This proves that relative entropy cannot contribute another leading \(R^4\) coefficient.

## One Cauchy record for moving observers

Give an observer a direction \(n\) and rapidity \(\eta\). Use the one-particle generator \(B_{n,a,L}=\{g_{n,L},H_a\}/2\), where \(g_{n,L}=(n\cdot x)\chi_L(x)\) and the smooth torus taper equals one in the central cube of side \(L/2\). Set \(U=\exp(-i\eta B_{n,a,L})\). Tilted ball mode columns are \(U J_B\). Their covariance matrices are \(J_B^\dagger U^\dagger C_i UJ_B\). Different observers operate on the same global generators and the same preparations. The finite source can be rerun on copies to read incompatible observables. The independent tilted energy reading conjugates the weighted hopping observable by that same \(U\). The velocity convention is \(u=(\cosh\eta,-n\sinh\eta)\). A boost here uses the supplied Weyl signal and laboratory clock.

The continuum Noether boost on time-zero spinors has generator \(\{n\cdot x,H\}/2\). One can check this by differentiating the spinor boost and using its evolution equation: the derivative of the retimed field is \((n\cdot x)\partial_t\), and the spin term supplies the derivative of \(n\cdot x\) in the symmetrization. The continuum vacuum is invariant. Opposite signed mass-shell four-momenta remain opposites under the boost, and their spectral spin projectors agree after simultaneous energy and momentum reversal. In the new frame the symbol still has the form \(f_u(p)P_+(p)-f_u(-p)P_-(p)\), with \(f_u\) the positive-shell scalar occupation pulled back by the Lorentz transformation. In a Cauchy-space Fourier normalization, the Jacobian and unitary wavefunction factors cancel in conjugating this multiplication operator. Its momentum support remains bounded for finitely many rapidities. The previous odd-kernel and domain estimates apply with \(T_{uu}\) from the retimed signal in place of \(\rho\).

## Order of finite limits

For each fixed \(R,\delta>0\) the finite lattice can approximate those readings to arbitrary absolute accuracy. The following order justifies a single common diagonal.

1. At fixed torus length \(L\), send the cell spacing \(a\) to zero. Embed lattice spinors by normalized piecewise-constant cells. Ball multipliers converge strongly since a sphere has measure-zero boundary. On each fixed Fourier mode the negative-band spin projector converges to the Weyl projector, and the cell embedding converges in norm. Fourier polynomials are dense and the projections bounded, giving strong convergence on the entire Cauchy space. The extra lattice cones do not affect this argument. The excitation on this torus is a finite sum of rank-one plane waves in the fixed momentum window, giving trace-norm convergence after restriction.
2. At fixed \(L\) with smooth taper, the discrete boosts converge strongly for finite rapidity. On smooth torus data, centered differences converge in norm to every first derivative. For the continuum evolution generated by the smooth symmetric first-order operator \(\{g,H\}/2\), all required derivative norms are bounded on a compact rapidity interval. Differentiate the difference between sampled continuum evolution and discrete unitary evolution and integrate its consistency error (Duhamel's formula). The same proof with \(-\eta\) handles adjoints; density and unitary norms extend it to \(L^2\). This preserves the previous strong reference and finite-rank source convergence on tilted modes.
3. Send \(L\) to infinity. In the central region a compact test spinor transported by a continuum boost has support within a finite larger ball: the characteristic speed is bounded by \(|x|\) where the taper equals one. For a fixed ball, center and rapidity, a large enough box puts this whole domain inside the untapered chart. By finite propagation, shrinking-ball readouts then use the true Minkowski boost on their compact test data. The torus vacuum projectors converge strongly on localized spinors after zero extension of centered boxes. For smooth compact tests their matrix elements converge by Fourier Riemann sums; convergence of the projected norms follows from \(P_-^2=P_-\). Localization and unitarity give the same statement on the boosted subspaces. For the source, restrictions of the transformed low-energy plane waves converge as rank-one kernels in the Riemann sums. The boosts agree locally with the true boost on each such smooth plane wave by the same domain-of-dependence argument. Bounded momentum support gives convergence of the restricted sums in trace norm.
4. For these two limits hold the **physical** \(\delta\) fixed. On \([\delta,1-\delta]\) the log odds have a uniform polynomial approximation. Bounded strong functional calculus on the reference and trace-norm convergence on the restricted source therefore give convergence of \(\Delta K\). The Hilbert–Schmidt norm of the source also converges. For any fixed finite menu choose \(L\) large and then \(a\) small enough that all modular expectation errors are below \(R_{\min}^5\) and all squared Hilbert–Schmidt errors below \(R_{\min}^8\). The independent low-momentum action matrix elements also converge. For a boost, this last claim follows on each of the finitely many smooth torus source plane waves by the same smooth consistency estimate in one spatial derivative: a first finite difference of the evolution commutes past \(g\) with coefficients bounded by first and second differences of the smooth taper. A discrete energy estimate and Grönwall's inequality control the first-difference norm of the evolution error uniformly at fixed \(L,\eta\). Multiplication by the bounded weight then gives convergence of the weighted bond readings. Enlarge \(L\) first as necessary for the finite momentum Riemann sums.
5. Choose a sequence of menus with \(R_{\max}\to0\), and use the same \(\delta=\kappa R_{\max}\), \(L\), and \(a\) for every observer and radius in each menu. The preceding endpoint argument and finite entropy bound give the displayed assertions. The proof supplies existence of a diagonal. The small numerical boards below do not certify an error rate for a particular relationship between \(a,L\) and \(R\). An arbitrarily fine local grid at fixed nonzero detector floor would retain a different modular operator.

Translations and a finite number of different centers are allowed while their measurement domains lie in the untapered chart. This comparison assumes the Weyl phase, its ball algebras and signal transformations. Recovering these data from native Interval Echo publications, and demonstrating the independent fixed-proper-volume area response on the same physical preparation, remain separate physical questions.

## Freely evolving and charged wave packets

Choose finitely many orthonormal positive-energy wave packets and finitely many orthonormal negative-energy packets. Give their Fourier spinors smooth compact support away from zero momentum. Add occupations between zero and one in the first set and remove occupations between zero and one in the second. This defines physical covariances at fixed excitation amplitude. The single hopping evolves the common preparation freely. Each finite frame maps the packet profiles into smooth compact momentum profiles. On the lattice, sample their Fourier amplitudes, project onto the corresponding bands, and orthonormalize within each band. For a sufficiently fine torus the Gram matrices approach the identity. All the restrictions then have the trace-norm approximation used above, including after finite free evolution and retiming.

Write \(F(x,y)\) for the signal covariance difference in one frame at the event \(x=y=0\). The rescaled kernel now has a constant term \(qF(0,0)/R\), which can be nonzero. On the centered continuum ball inversion of the spatial argument, \(Jv(x)=v(-x)\) with unchanged spinor components, satisfies \(JC_BJ=I-C_B\) because \(P_-(-p)=P_+(p)\). Thus \(Jk_\delta(C_B)J=-k_\delta(C_B)\). The constant-kernel operator is even under this inversion. Its trace pairing with the measured modular matrix vanishes exactly at every \(\delta\). This is an operator symmetry of the vacuum restriction and works for packets carrying charge.

The next kernel term is \(q\sum_j(x_j X_j+y_j Y_j)\), with \(X_j=\partial_{x_j}F(0,0)\) and \(Y_j=\partial_{y_j}F(0,0)\). The residual has trace norm \(O(R)\) by the same finite-rank Taylor estimate applied to the packet factors. Hence its pairing with the bounded log is \(O(R\log(1/\delta))=o(1)\). The affine term uses the same endpoint-domain proof. Integrating \(\int_{B_1}(\partial_j w_1)x_k=-\delta_{jk}\int_{B_1}w_1\) gives the local action density
\[
 T_{uu}(0)=-\frac{i}{2}\sum_j\operatorname{tr}\big[\sigma_j(X_j-Y_j)\big],
 \qquad R^{-4}\Delta K\longrightarrow\frac{8\pi^2}{15}T_{uu}(0).
\]
The bounded smooth covariance kernel now gives \(\|\Delta C_B\|_{\rm HS}^2=O(R^6)\), so \(D=O(R^5)\) still drops out of the leading entropy response. The bounds hold for each named event, radius and frame on the same finite diagonal. Small fixed-radius smooth spinor probes also admit an operator check: first take the finite regulator to the continuum at fixed detector floor, then remove the floor on a fixed finite-dimensional span in the modular domain. Strong convergence on a basis gives uniform convergence on its unit sphere. The result concerns these smooth inputs; a norm estimate over all ultraviolet lattice modes would demand much more.

## An axial area test distinguishes the existing rod geometry

The published parallel-sheet apparatus has metric \(ds^2=-dt^2+dz^2+ds_\perp^2\), where the transverse curvature \(K\) depends on its registered sheet densities. In an orthonormal frame its Ricci tensor has diagonal \((0,K,K,0)\). Hence the two axial null vectors \(k_\pm=(1,0,0,\pm1)\) obey \({\rm Ric}(k_\pm,k_\pm)=0\) for every preparation. A cosmological term vanishes under both contractions. In the same conventions the covariant Weyl excitation stress is
\[
 T_{ab}=2\int f(p)\frac{p_a p_b}{|p|}\frac{d^3p}{(2\pi)^3},\qquad
 p_a=(-|p|,\boldsymbol p).
\]
The half sum of its axial null contractions is \(\rho+T_{zz}\ge\rho>0\) for every nonzero positive-energy particle/hole occupation of this class. For the tested dipole \(T_{zz}=\rho/3\). Thus even calibrating sheet number to match the rest energy exactly cannot reproduce the two directional responses at any positive coupling. Equivalently, the stress difference for this Weyl source has zero spacetime trace and the change of the sheet Einstein trace is \(-2\Delta K\). [Parallel-sheet action and rod geometry](../../../../../model/candidate_laws/interval_echo_partition_repair/published_strings/GEOMETRY.md)

There is a rest-and-moving-observer test using only equal-volume boundary readings. For a common excitation let \(\Delta A_u\) denote the area change at fixed proper volume \(4\pi R^3/3\); let \(u_0\) be at rest and \(u_\pm=(\cosh\eta,0,0,\pm\sinh\eta)\). A same-source Einstein response at coupling \(\kappa\) requires
\[
 \lim_{R\to0}-\frac{15}{4\pi R^4}
 \left[\frac{\Delta A_{u_+}+\Delta A_{u_-}}{2}-\Delta A_{u_0}\right]
 =\kappa(\rho+T_{zz})\sinh^2\eta .
\]
The difference of the two preparations removes a constant vacuum term. Averaging the opposite boosts removes the sign convention for momentum. A boost along the sheets is an exact isometry of their product metric. The three sheet areas therefore agree at every common proper volume and the left hand side vanishes exactly for that geometry. For the Weyl dipole at \(\eta=0.32\) the mean boosted stress is \(1.141258\) times the rest stress. Matching the rest densities at 0.041228 gives Weyl null readings 0.051620 and 0.058321, against the two zero sheet readings divided by the coupling. The numerical quadrature integrates the dipole's angular moments through degree three and its radial polynomials through degree eight. [Same-spectrum axial calculations](axial_response.json)

The independently calibrated free Weyl signal itself gives a second check. The canonical anticommutator is an operator identity at every finite resolution,
\[
 \{b_i(t),b_j^\dagger(0)\}=(e^{-itH_a})_{ij} I .
\]
Indeed \(b_i(t)=\sum_m(e^{-itH_a})_{im}b_m(0)\). In a preparation of covariance \(C\) the two ordered correlations are \(e^{-itH_a}(I-C)\) and \(e^{-itH_a}C\) with the appropriate common row and column. A detector comparing them can measure a source-dependent occupation even as their sum, the causal linear signal for positive time, remains fixed. On the finite board the excitation changes regional entropy by 0.280740 and changes one ordered-signal column by up to 0.013161. The sum reproduces the free signal in both preparations to \(6\times10^{-17}\). A rod and clock protocol that uses only this summed kernel receives identical input on both runs. A source-dependent area would require an independently calibrated geometric response in another common record. [Finite signal readings](axial_response.json)

The homogeneous Einstein–scalar radar experiment supplies equal-volume geodesic balls and a proper clock for six **scalar oscillator** modes. Its source evolution explicitly couples their stress to its metric. The spatially flat FLRW record-density construction reads a supplied conformal-factor profile and identifies the source law and physical clock as further inputs. The native critical exchange realizes a spatial interval mode and an energy analyzer. None of these protocols acquires both the Weyl covariance above and its independently shifting proper-ball boundaries on a single preparation. The three-area equality supplies a falsifiable target for any proposed identification before attempting the full tensor in additional frames. [Scalar field and radar source](../../../quantum_scalar_radar/PROTOCOL.md) · [Scalar equal-volume protocol](../../../quantum_scalar_radar/BALL_PROTOCOL.md) · [FLRW count-measure derivation](../../../../../paper/tex_fragments/SOURCE_NET_FLRW_RECORD_DENSITY.tex) · [Native interval selection](../../../../../model/experiments/interval_echo_einstein_station/critical_exchange/selection/README.md)
