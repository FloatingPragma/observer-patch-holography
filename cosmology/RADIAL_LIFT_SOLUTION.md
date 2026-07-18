# Radial-Lift Theorem Packet

The radial-lift problem has a precise mathematical boundary. Restricting a
homogeneous three-dimensional field to one shell gives an exact angular
covariance, while unrestricted recovery of the three-dimensional spectrum from
that shell is nonunique. The theorem packet proves the obstruction and gives
two routes to uniqueness.

## Exact shell and window projection

For a normalized radial window \(W\), screen normalization \(Z_q\), and
primordial curvature spectrum \(\Delta_\zeta^2\), the flat-branch projection is

\[
C_{\ell,W}^{q}
=4\pi Z_q^2\int_0^\infty d\ln k\,
\Delta_\zeta^2(k)|\Psi_{\ell,W}(k)|^2,
\qquad
\Psi_{\ell,W}(k)=\int W(dr)j_\ell(kr).
\]

A thin shell at \(R_\star\) has
\(\Psi_{\ell,W}(k)=j_\ell(kR_\star)\). This formula is a forward map from
three-dimensional source power to screen covariance.

## Complete one-shell obstruction

The full angular covariance on a shell of radius \(R\) is

\[
K_R(\mu)
=Z_q^2\xi_\zeta\!\left(R\sqrt{2-2\mu}\right)
=\sum_{\ell\ge0}\frac{2\ell+1}{4\pi}C_{\ell,R}^qP_\ell(\mu).
\]

The variable \(R\sqrt{2-2\mu}\) covers the chord interval
\(0\le s\le2R\). Every shell multipole therefore determines the
three-dimensional correlation function only on this interval.

The kernel is infinite-dimensional. For every nonzero smooth function \(g\)
supported beyond \(2R\), define

\[
h_g(k)=\frac{2k^3}{\pi}\int_0^\infty r^2g(r)j_0(kr)\,dr.
\]

Radial Fourier inversion gives \(\xi_{h_g}=g\), so \(h_g\) changes no shell
multipole. Choosing \(B=|h_g|+w\) for any strictly positive integrable \(w\)
produces two distinct positive spectra, \(B\) and \(B+h_g\), with the same
complete shell covariance. Positivity does not remove the radial ambiguity.

## Physical source-dilation route

The source construction uses the scale-labelled cofinal refinement orbit,
rather than the angular modes of one cut. At finite regulator \(r\), let
\(U_r\) embed this source orbit into the physical adiabatic mode space, let
\(R_{s,r}\) implement source refinement, and let \(D_{s,r}\) implement physical
log-wavenumber dilation. The required commuting square is

\[
D_{s,r}U_r=U_rR_{s,r},
\qquad
C_{\zeta,r}=U_rC_{{\rm src},r}U_r^*.
\]

It gives the exact residual transfer

\[
D_{s,r}^{-1}C_{\zeta,r}D_{s,r}-u_r(s)C_{\zeta,r}
=U_r\left(
R_{s,r}^{-1}C_{{\rm src},r}R_{s,r}
-u_r(s)C_{{\rm src},r}
\right)U_r^*.
\]

The source and physical residuals therefore have equal operator norm. Strong
convergence of the covariance, dilation, and inverse dilation operators, a
uniform covariance norm bound, convergence of \(u_r\), and a vanishing
finite-stage residual imply

\[
D_s^{-1}C_\zeta D_s=u(s)C_\zeta.
\]

If the source survival cocycle has generator \(\theta\), then
\(u(s)=e^{-\theta s}\). On the physical \(d\ln k\) mode space,
\(C_\zeta=M_{\Delta_\zeta^2}\), so conjugation forces

\[
\Delta_\zeta^2(k)=A_\zeta(k/k_\star)^{-\theta}
\]

almost everywhere. A continuous source spectrum satisfies this identity
pointwise.

## Radial-tomography route

Complete radial cross-covariances give the independent unrestricted solution:

\[
C_\ell(r,r')
=\frac2\pi\int_0^\infty k^2P_\zeta(k)
j_\ell(kr)j_\ell(kr')\,dk.
\]

The spherical Hankel transform diagonalizes this covariance operator,

\[
Q_\ell=\mathcal H_\ell^{-1}M_{P_\zeta}\mathcal H_\ell,
\]

and recovers \(P_\zeta\) almost everywhere. One complete cross-radius family
\(C_\ell(r,r_0)\), known for all \(r\), also suffices because the zeros of
\(j_\ell(kr_0)\) are discrete. Finite auto-window spectra do not supply this
tomography theorem.

## Thin-shell amplitude and finite windows

On the flat thin-shell source-dilation branch,

\[
\int_0^\infty\frac{dx}{x}x^{-\theta}j_\ell^2(x)
=\frac{\sqrt\pi}{4}
\frac{\Gamma(1+\theta/2)}{\Gamma(3/2+\theta/2)}
\frac{\Gamma(\ell-\theta/2)}{\Gamma(\ell+2+\theta/2)}
\]

for \(-2<\theta<2\ell\). For all retained multipoles
\(\ell\ge2\), the common strip is \(-2<\theta<4\), and

\[
C_{\ell,R_\star}^{q}
=A_q\frac{\Gamma(\ell-\theta/2)}
{\Gamma(\ell+2+\theta/2)},
\]

\[
A_\zeta
=\frac{A_q}{\pi^{3/2}Z_q^2(k_\star R_\star)^\theta}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
\]

Finite source thickness uses the exact window transform and a weighted Hilbert
space bound. The bound preserves the ultraviolet decay of the Bessel kernel and
is recorded alongside the exact quadrature result.

## Claim boundary

The mathematical radial packet is complete at the conditional theorem and
algorithm level. A physical primordial source requires one finite evidence
bundle containing the source embedding or cross-covariance tomography,
source-derived amplitude, physical scale and mode basis, window, null-space,
positivity, convergence, and non-fitting forward residual. The simulator
contract records these objects in a source artifact tree and keeps observational
comparisons in a separate diagnostic tree.

A finite prior or regularizer produces a conditional radial continuation.
Temperature and polarization spectra require the independent transfer and
likelihood system.
