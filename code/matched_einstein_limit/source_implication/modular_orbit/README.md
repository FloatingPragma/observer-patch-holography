# Modular time on a volume-matched collar perturbation

The 80-cell Interval Echo station admits a prepared perturbation with constant survey volume, parked energy and each marginal entropy. The total head entropy increases and its survey area decreases. On this perturbation the full head modular clock leaves the response of every candidate null charge constant. Restricting the clock and charge to the right coin gives a varying signal with a different time dependence from affine null translation. These results concern a particular proposed identification of a modular clock with a boost. The station uses a measured tetrahedral survey. The physical Einstein premise concerns proper balls and independently measured stress on those balls. [Station and preparation](../full_algebra/PROOF.md) · [Einstein tower, finite modular covariance and compressed responses](../../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

## One source and a signed collar signal

Use the faithful local stationary source \(\chi=\chi_e\oplus\chi_o\) of the six-unit survey. Its ordinary calibration \(\eta\) and volume are frozen. The left marginal includes the parity label. Write the right coin marginal as \(\chi_R\), the two left marginal blocks as \(\chi_{L,b}\), and the product reference as \(\pi=\bigoplus_b\chi_{L,b}\otimes\chi_R\). Let \(b=e,o\), \(j=1,2\), and define the following four directions on their respective blocks.
\[
 H_e=(\chi_e-I/5)/(1/10),\qquad H_o=(\chi_o-I/20)/(3/100),
 \qquad D_{bj}=H_b^j-\tfrac14\operatorname{Tr}(H_b^j)I .
\]
Every \(\rho(t)=\chi+\sum_{bj}t_{bj}D_{bj}\) commutes with \(\chi\) and preserves each block weight. Hence its parked energy \(16+\operatorname{Tr}\rho_o\) is constant. The constraint map is \((V,S_L,S_R)(\rho(t))\). Its Jacobian minor on \(e1,o1,o2\) lies between \(-3.463\times10^{-5}\) and \(-3.380\times10^{-5}\). The implicit function theorem gives an actual smooth curve through \(\chi\) with all three constraints fixed and \(t_{e2}\) as its parameter. The station’s positive preparation weights cover a neighborhood of \(\chi\), so the ideal material source lotteries implement the local curve. Along the curve, the marginal density matrices can change orientation. [Interval implementation](measure.py) · [Preparation radius](../full_algebra/PROOF.md)

Let \(D\) be the tangent for increasing \(t_{e2}\). All slopes below are per unit of this parameter. Intervals are rounded outwards from the certificate.

| Native station response | Certified interval |
|---|---:|
| Total entropy \(dS\), nats | \([0.0076742,\,0.0079773]\) |
| Tetrahedral survey area \(dA\) | \([-0.00035143,\,-0.00032600]\) |
| Right coin \(d\langle Z\rangle\) | \([0.0016451,\,0.0044043]\) |
| \(dV,\ dS_L,\ dS_R,\ dE_{\rm parked}\) | exactly zero |

There are two independent calculations of the area sign. The measured ranges give the interval in the table. Stationarity on the fixed-volume surface also gives \(dA=-dS/\eta\). Since \(I(L:R)=S_L+S_R-S\), all the entropy change equals a loss of mutual information. The frozen collar logarithm \(K_c=\log\pi-\log\chi\) has \(\operatorname{Tr}DK_c=dS>0\). This is a nonzero collar response at fixed ordinary energy. Determining \(T_{kk}\) from that response requires an independent mapping of charges, directions and physical scales. [Numerical values and enclosures](results.json)

For the enclosures, the Hilbert–Schmidt radius of the certified stationary state is \(2\times10^{-7}\). If the scaled block error is \(u\), the linear and squared block powers have respective errors \(u\) and \((2+u)u\), since the rational central matrices have operator norm below one. Subtracting the trace projection is a Hilbert–Schmidt contraction. Rational read effects, interval derivatives of the survey geometry and spectral log enclosures bound the constraint Jacobian. Cramer's rule then bounds \(D\). The rate calculations subtract multiples of the three exactly vanishing constraint rows before evaluating the intervals. The binary floating-point multipliers serve only to choose that linear combination. [Exact kernels and spectral proof](../full_algebra/native_precision.py) · [Certificate calculation](measure.py)

## Full head modular time

Consider a candidate boost \(\alpha_s(P)=\chi^{is}P\chi^{-is}\) and a candidate null translation charge \(P=P^\ast\). A boost dilating the affine coordinate scales the translation generator by \(e^{\pm2\pi s}\), with sign set by ray orientation. This concerns the translation charge extracted from a nested inclusion. The modular Hamiltonian of a single cap is a different generator. For a state contrast define \(F(s)=\operatorname{Tr}(D\alpha_s(P))\). Since \([D,\chi]=0\), cyclicity gives \(F(s)=F(0)\) for all real \(s\) and all \(P\). The statement also holds for finite differences along the curve. At the modular time with scale factor two the affine discrepancy equals \(|F(0)|\). Any charge which sees the measured collar contrast therefore fails this boost identification. A null charge with \(F(0)=0\) is consistent with this test. [Modular inclusion theorem](https://arxiv.org/abs/math/0412061) · [Affine group and sign convention](https://link.springer.com/article/10.1007/s00220-022-04432-8)

The same conclusion holds for any number of independent copies, arbitrary collective measurements and independent stationary spectators on their supports. Tensor products of the perturbed states commute with the reference product. A common unitary applied to reference, preparations and instruments preserves that statement. In the canonical Hilbert–Schmidt realization of a density matrix by its square root, \(\Delta_\chi^{is}(X)=\chi^{is}X\chi^{-is}\) fixes \(\sqrt{\rho(t)}\). Even a quadratic-form expectation on these canonical vectors is constant wherever defined. Preparations on a larger entangled system can select other extensions of the reduced states.

There is a quantitative necessity on alternative regional inputs. For \(U_s=\sigma^{is}\), a self-adjoint bounded \(P\), and \(\mathrm{osc}(P)=\lambda_\mathrm{max}(P)-\lambda_\mathrm{min}(P)\), duality gives
\[
 |F(s)-F(0)|\le\tfrac12\,\mathrm{osc}(P)\Vert U_s^\ast D U_s-D\Vert_1
 \le \tfrac12\,\mathrm{osc}(P)|s|\Vert[\log\sigma,D]\Vert_1 .
\]
The first bound subtracts the midpoint of the spectrum of \(P\). The second integrates the modular derivative. If the doubling residual is at most \(\epsilon|F(0)|\), the middle expression must be at least \((1-\epsilon)|F(0)|\). Comparing errors at the scale of the measured source contrast matters even when both the total signal and its absolute residual tend to zero.

## Right coin and a calibrated control

Restricting the perturbed head to its right coin changes the experiment. Write \(\chi_R=(I+\mathbf r\cdot\boldsymbol\sigma)/2\) and \(D_R=\mathbf d\cdot\boldsymbol\sigma/2\). A two-state entropy is a strictly monotone function of \(|\mathbf r|\) away from zero. The measured \(\mathbf r\ne0\), and constant right entropy implies \(\mathbf r\cdot\mathbf d=0\). The independent \(Z\) reading proves \(\mathbf d\ne0\). Modular conjugation on this coin rotates about \(\mathbf r\) with angular speed \(\omega=\log[(1+|\mathbf r|)/(1-|\mathbf r|)]\) in \([0.0531709,0.0531748]\). For every self-adjoint right analyzer,
\[
 F_R(s)=A\cos(\omega s)+B\sin(\omega s),\quad
 F_R(s)+F_R(-s)=2\cos(\omega s)F_R(0).
\]
Set \(s_0=\log(2)/(2\pi)\) and normalize errors to \(|F_R(0)|>0\). The desired values at \(-s_0,+s_0\) sum to \(5F_R(0)/2\). At least one error is at least \(5/4-\cos(\omega s_0)>0.25001720\). Using only the forward orientation, the recurrence at \(0,s_0,2s_0\) gives a lower bound \((5-4\cos(\omega s_0))/(1+2\cos(\omega s_0))>0.33336009\) on the larger of the two forward errors. Tensor powers with a first-order perturbation of one or more independent, identical coins obey the same recurrence. In each term all other density matrices supply stationary expectations.

One can tune a positive rank-one analyzer so that \(F_R'(0)=2\pi F_R(0)\) exactly. Write \(\mathbf n=\mathbf r/|\mathbf r|\), \(\mathbf m=\mathbf d/|\mathbf d|\), and choose the analyzer axis parallel to \(\omega\mathbf m+2\pi(\mathbf n\times\mathbf m)\). The normalized response becomes \(\cos(\omega s)+(2\pi/\omega)\sin(\omega s)\). At the rounded stationary source the native rail circuit evaluates the response ratios \(0.306840,1,1.693126,2.386194\) where the affine ratios are \(1/2,1,2,4\). The largest operator residual of the rail rotation is below \(3\times10^{-16}\). These are computed expectation values. The derivative calibration is an input of this control. [Rail circuit and finite source preparation](controls.py) · [Output](results.json)

## A finite-frequency bound for other cap choices

A weak finite modular orbit with at most \(m\) distinct frequencies gives samples \(F_j=\sum_{k=1}^m c_k z_k^j\) on an equally spaced modular-time grid, where \(|z_k|=1\). This follows by diagonalizing the regional reference. A faithful \(d\)-dimensional density matrix has at most \(d(d-1)+1\) distinct modular frequencies, including zero. The recurrence for \(p(z)=\prod_k(z-z_k)=\sum_{j=0}^m a_jz^j\) is \(\sum_j a_jF_j=0\). Padding the list of roots with ones allows \(m\) to be any uniform upper bound.

Suppose \(F_0\ne0\) and the desired samples are \(r^jF_0\), with \(r>1\). The residuals \(E_j=F_j-r^jF_0\) satisfy \(\sum_j a_jE_j=-p(r)F_0\) and \(E_0=0\). Since \(|p(r)|\ge(r-1)^m\), \(|a_0|=1\) and \(\sum_j|a_j|\le2^m\), one has
\[
 \max_{1\le j\le m}\frac{|E_j|}{|F_0|}\ \ge\ \frac{(r-1)^m}{2^m-1}.
\]
Choosing \(r=2^{1/m}\) keeps all samples inside one prescribed doubling window. The positive lower bound is uniform in the spectra, analyzer norms and source coefficients. An approximating family retaining nonzero normalized charge and converging at all those times must have unboundedly many active modular frequencies. This bound can be small. Growing frequency count alone gives no convergence theorem. The station-specific right-coin identity gives the larger bounds above.

For example, take a product reference of identical cells and a charge that is a sum of operators on at most \(k\) cells each. The perturbation may correlate all cells. A local charge term acquires only sums of \(k\) one-cell frequencies, including zeros for smaller supports. If one cell has at most \(L\) possible frequencies, every such charge has at most \(\binom{L+k-1}{k}\) frequencies regardless of cell count, overlap of supports, or the number of terms. For identical right coins the frequencies are \(j\omega\), \(-k\le j\le k\), so \(m=2k+1\) suffices. Consequently even fully correlated changes of the prepared state cannot establish the affine weak limit with a product reference and charges of uniformly bounded support size per term. Connected regional references can have a different modular spectrum.

A physical modular route would need a family of regional algebras on the same dynamics whose reduced source, modular evolution and independently read null stress have the required weak scaling at the retained order. A direct optical proof of physical null balance can use other premises. In particular QNEC bounds the second *shape* derivative along a physical null cut by a null stress. The entropy derivatives in the table vary a *prepared state* at fixed survey volume. Equality of the optical focusing and the ordinary-plus-collar stress on shrinking proper balls is a separate measurement identity. [Quantum null energy condition](https://arxiv.org/abs/1706.09432) · [Finite null-net audit](../../../geometry/NULL_NET_AUDIT.md)

Run from the repository root with `PYTHONPATH=code:. .venv/bin/python -m matched_einstein_limit.source_implication.modular_orbit.measure`. The `--write` switch regenerates the JSON output.

A different protected sector has an exact finite energy representation of its eight nonconstant angular entropy components across twelve overlapping stars. Its positive interior fit requires an excess on the rim, and a separately prepared equal-energy contrast rejects the monotone interior kernel. [Angular energy operators and signed matter records](angular_energy/README.md)

Complete relative-entropy remainders on the stationary angular source and the protected routed clock family determine which signed responses survive at small source strength. [Nested regions and the entropy-selected tag reference](relative_entropy/README.md)
