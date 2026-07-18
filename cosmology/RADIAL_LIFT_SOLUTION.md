# The radial lift

The radial lift turns a pattern on an observer screen into a
three-dimensional primordial spectrum. The screen carries an angular spectrum,
which measures variation by angular scale. Cosmology codes consume
\(\Delta_\zeta^2(k)\), which measures primordial fluctuations by physical
wavelength.

One spherical screen cannot recover an arbitrary three-dimensional spectrum.
The screen measures correlations between pairs of points on the sphere, whose
separation never exceeds the diameter. Even a noiseless screen containing every
angular multipole determines the spatial correlation only for
\(0\le s\le2R\). Changes confined to larger separations are invisible.

This ambiguity is infinite-dimensional and survives positivity. Distinct
positive primordial spectra can produce exactly the same complete screen
covariance. A numerical inverse, positivity constraint, or regularizer
therefore selects one representative without proving that the source selected
it.

## Physical source dilation

OPH carries information about how the source screen is refined. A physical
radial theorem must show that this refinement acts as a change of wavelength in
the three-dimensional mode space.

At finite resolution, the source embedding \(U_r\), source refinement
\(R_{s,r}\), and physical dilation \(D_{s,r}\) satisfy

\[
D_{s,r}U_r=U_rR_{s,r},
\qquad
C_{\zeta,r}=U_rC_{{\rm src},r}U_r^*.
\]

This commuting square says that refinement followed by the physical embedding
agrees with physical embedding followed by wavelength dilation. It also
transfers the source covariance error exactly to the physical covariance.

Strong convergence of the finite covariance and dilation operators, a uniform
covariance norm bound, and a vanishing operator residual give the continuum
identity

\[
D_s^{-1}C_\zeta D_s=e^{-\theta s}C_\zeta.
\]

On the physical \(d\ln k\) mode space this identity forces

\[
\Delta_\zeta^2(k)
=A_\zeta\left(\frac{k}{k_\star}\right)^{-\theta}.
\]

Once the source law fixes \(\theta\), the screen has one amplitude left to
determine. Every retained multipole must infer the same amplitude, which turns
the unused multipoles into consistency checks.

## Radial tomography

The independent route uses correlations between different radii,

\[
C_\ell(r,r')
=\frac2\pi\int_0^\infty k^2P_\zeta(k)
j_\ell(kr)j_\ell(kr')\,dk.
\]

The spherical Hankel transform diagonalizes this covariance operator and
recovers \(P_\zeta(k)\) almost everywhere. A complete family
\(C_\ell(r,r_0)\), with \(r\) varying continuously, also suffices. A finite
collection of auto-spectra does not supply this tomography theorem.

## Exact amplitude conversion

For a flat thin shell,

\[
C_{\ell,R_\star}^{q}
=4\pi Z_q^2\int_0^\infty d\ln k\,
\Delta_\zeta^2(k)j_\ell^2(kR_\star).
\]

The source-dilation power law gives the exact gamma-ratio spectrum and

\[
A_\zeta
=\frac{A_q}{\pi^{3/2}Z_q^2(k_\star R_\star)^\theta}
\frac{\Gamma(3/2+\theta/2)}{\Gamma(1+\theta/2)}.
\]

At the source-native pivot with \(Z_q=1\), the conversion factor is about
\(0.16080676\). It is an arithmetic consequence of the symbolic relation. It
does not supply a source amplitude by itself.

A source with radial thickness uses the exact window transform. The finite
window theorem bounds the difference from a thin shell in the weighted Hilbert
norm of the projection, preserving the large-wavenumber decay of the Bessel
kernel.

## Finite computation

A finite radial basis gives a matrix equation \(C=Ap\). The singular values,
rank, right-null basis, resolution kernels, positivity active set, and prior
sensitivity expose what the angular data determine. For prior center \(p_0\)
and \(Q\succ0\), the minimum-prior representative is

\[
p_*=p_0+Q^{-1}A^{\mathsf T}
(AQ^{-1}A^{\mathsf T})^+(C-Ap_0).
\]

This representative is a conditional radial continuation. It is useful for
simulation and cannot satisfy a source-derived promotion gate.

The ordinary spherical-Bessel and gamma-function formulas apply to flat
spatial branches. Open and closed spatial branches require their
hyperspherical eigenfunctions and spectral measures. An unresolved curvature
branch blocks physical radial promotion.

## Claim boundary

The mathematical packet contains the exact forward projection, complete
one-shell obstruction, positive counterexamples, physical source-dilation
theorem, radial tomography, amplitude conversion, finite-window bound,
finite-basis continuation, and fail-closed promotion rule.

A physical primordial source requires one finite evidence bundle containing
the source embedding or radial tomography, source-derived amplitude, physical
scale and mode basis, window, null-space, positivity, convergence, and
non-fitting forward residual. Temperature and polarization spectra require the
independent transfer and likelihood system.
