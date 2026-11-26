# The wall spectrum of a twice-read strand

## Propagating every recorded outcome

Let A and B be the two orientation amplitudes at the second reading when the intervening coin parity is even or odd. The parent measurement derives these paths from the native 99-state seam table. If the two labels match, the post-read amplitude in each final parity sector is (A+B)/4. If they differ it is (A−B)/4 up to a sector phase. Their squared norms determine eight joint answers per source. The final parities are orthogonal, so their probabilities add. After pooling all four words, each final parity has mass one half. This argument requires complete home records and the fixed two-read schedule. [Native two-read reduction](../FREQUENT.md)

For a right-side wall, the even final parity gives the intact strand. The odd parity evolves the nine right-side dimers independently of the five on the left. It retains the right reduced density produced before that split. On either segment the physical generator differs by a constant from
\[
H=-h\sum_i X_i+J\bigl(Z_0-Z_{n-1}+\sum_i Z_iZ_{i+1}\bigr).
\]
The linear terms are boundary bonds to two fixed spins with values +1 and −1. Set \(A_j=(\prod_{r<j}X_r)Z_j\) and \(B_j=(\prod_{r<j}X_r)Y_j\) on this extended chain. Then \(X_j=iA_jB_j\), \(Z_jZ_{j+1}=iB_jA_{j+1}\). The ordered active list starts with the left fixed spin's B operator and ends with the right fixed spin's A operator. Its real skew evolution generator K has superdiagonal entries alternating between 2J and −2h. Every even product acts inside the physical fixed-end sector. On arbitrary input amplitudes its expectation can be obtained by an explicit bit flip and sign on the native orientations. No Gaussian assumption is made about the postmeasurement state.

For a normalized state the antisymmetric covariance Γ has entries ⟨iγₐγᵦ⟩ off the diagonal. Its operator norm is at most one. The Gram matrix of the γ operators on the state and its conjugate are positive, giving both signs of the bound on iΓ. Arbitrary normalized mixtures satisfy the same bound. A wall click on bond j has probability \((1+[R\Gamma R^T]_{2j,2j+1})/2\), where \(R=\exp(tK)\). A subnormalized branch replaces 1 by its squared norm. The small-strand regression compares every coin/terminal answer at more than one receiver and time with direct full orientation propagation. The independent six-radian data on the 14-dimer strand supply another check.

## A cut-only spectral ceiling

Let U diagonalize iK and let the real signed frequencies be νₚ. They are the positive and negative singular values of the open bidiagonal matrix with diagonal 2J and lower diagonal −2h. Thus a wall record can contain νₚ−ν_q, including differences and sums of positive frequencies. Its response in a unit-sum Hann window w at ω is a superposition with multipliers \(W((\nu_p-\nu_q-\omega)\Delta)\), where \(W(x)=\sum_t w_t e^{itx}\). The sites are read independently on fresh copies. Both density profiles can display frequencies well below the constituent gap, even though each positive mode of the nominal material lies between 15 and 17.

For any two normalized input states the constant in their wall probabilities cancels and \(M=U^*(\Gamma_1-\Gamma_2)U/2\) satisfies \(\|M\|_F\le\sqrt d\), with d twice the number of dual bonds. The coefficient norm for a receiver (a,b) obeys
\[
|D(\omega)|\le\sqrt d\left(\sum_{pq}|U_{ap}|^2|U_{bq}|^2 |W_{pq}|^2\right)^{1/2}
\le\sqrt d\max_{pq}|W_{pq}|.
\]
This Hilbert–Schmidt bound holds for any preparation or retained reference at the second reading. It constrains an explanation in which the entire remote source difference subsequently follows the cut segment; arbitrary right-side imprint before the split is allowed. Testing both sectors separately over all candidate low difference frequencies gives one nominal rejection of the cut-only hypothesis and zero rejections of the intact-only hypothesis. This identifies a contribution requiring a coupled history, within the two known post-read sectors. The magnitude at the selected frequency is a finite-window coefficient. It need not identify one exact eigenpair.

### Contact and clock sensitivity

A supplementary enclosure uses the [contact-meter sensitivity box](../../../impulse/README.md) with J from 7.996 to 8.004 and h from 0.4997 to 0.5003. Write m=10 for the number of cut-segment bonds. Left singular modes are proportional to \(\sin((m-j)k)\) with \(\sin((m+1)k)=(h/J)\sin(mk)\). Right modes are their reversed partners. The roots stay real for \(0\le h/J<1\): the singular values lie between 2(J−h) and 2(J+h) by the triangle inequalities, with strict endpoints for the finite shift. Starting at the roots \(\ell\pi/(m+1)\) of the zero boundary ratio fixes their branches by continuity of the simple tridiagonal spectrum. Hence \(k\in[k_0,\pi-k_0]\) for \(k_0=(\pi-\arcsin(\max h/J))/(m+1)\). Summing \(\sin^2(jk)\) for j=1,…,m gives a normalization at least \(m/2-1/(2\sin k_0)\). Each eigenvector of iK uses the normalized left and right singular vectors with amplitudes divided by √2. Every \(|U_{ap}|^2\) is therefore at most 0.156733.

The 33 by 9 grid covers both coupling intervals. Each singular value varies by at most \(2(|\delta J|+|\delta h|)\) from a grid point. Twice this allowance encloses each frequency difference on the entire box. Multiplying by the two endpoints of the independently measured clock-spacing interval encloses its phase increment on the supplied uniform schedule. For N appointments and \(d_0=2\pi/N\) the normalized Hann window obeys
\[
|W(x)|=\frac{\sin^2(d_0/2)\,|\cos(x/2)\sin(Nx/2)|}
 {N|\sin(x/2)\sin((x-d_0)/2)\sin((x+d_0)/2)|}
\]
by summing three geometric series, with the continuous limits at removable singularities. Use \(|W|\le1\) there. Sine is 1-Lipschitz, so the enclosed phase interval bounds the numerator and each denominator term, giving an upper envelope \(W^*_{pq}\). The first Hilbert–Schmidt bound above is at most \(\sqrt{2m}\,0.156733\,\|W^*\|_F\). The [numeric enclosure](sensitivity.py) uses padding for numerical evaluation. The resulting 0.006861 ceiling is an engineering bound on the stated box. Its clock interval and count resolution have simultaneous sampling coverage. Arbitrary nonuniform clock drift or variation in the material along the strand lies outside this test.

## Count and wavelength resolution

The 4,096-point positive frequency scan uses a periodic Hann window with \(\sum_t w_t^2=3/(2N)\). At each site the contrast is the difference of two independently sampled Bernoulli means, each pooled over all head labels. For either real or imaginary Fourier component Hoeffding bounds the error by \(\sqrt{\|w\|_2^2\log(2M/\alpha)/S}\), where S is the shot count. The decoder bounds both components, every enrolled bin, law, site, and checked nominal sector frequency simultaneously. Detector totals are drawn by conditional binomial quantiles on the integer probability grid. The report includes their grid allowance. The total error allocation is 0.004 for Fourier contrasts, 0.004 for the all-appointment differences and 0.002 for the clock quadratures. The seven-by-two complex profile matrix has spectral-norm error at most √14 times the complex component radius. Its second singular value thus tests whether one fixed spatial profile suffices for the two bins.

Any superposition of reflected and transmitted waves at one real k obeys \(v_{j+2}+v_j=2\cos(k)v_{j+1}\). A simultaneous complex error ε allows a residual at most 4ε. The decoder intersects the intervals of real \(\cos(k)\) permitted by all five overlapping triples. This uses only held detector amplitudes. The finite-window rejection precludes assigning one wavelength to either measured profile by this recurrence. Neighboring unresolved lines or several spatial modes can produce such a profile. [Checked recurrence](SpatialPair.lean)
