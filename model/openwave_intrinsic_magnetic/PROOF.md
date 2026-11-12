# Stationary intrinsic current and its finite magnetic field

## Declared envelope architecture

Let \(\gamma_j,\beta\) be the complete sixteen-row matrices of
`native_mass_channels.folded_generators`, using the scalar mass channel.
They satisfy \(\{\gamma_i,\gamma_j\}=2\delta_{ij}\) and
\(\{\beta,\gamma_j\}=0\). The eight folded valley labels are retained.
The executor constructs an exact unitary transforming these matrices into
four identical four-row Dirac representations. This is a full unitary
change of basis. No internal row or spectral sign is discarded.

On an odd open cubic matter grid of spacing \(a\), declare a link
\((T_j)_{x,x+e_j}=e^{iaA_{xj}}\). Its adjoint and every other spatial row
are included. Define

\[
D_W(A)=\sum_j\gamma_j\frac{i(T_j-T_j^\dagger)}{2a}
+\beta M(A),\qquad
M(A)=m+\eta|x|^2+
\frac{r}{2a}\sum_j(2-T_j-T_j^\dagger).
\tag{1}
\]

The scalar mass well and open boundary confine the packet. The last term
is a declared gauge-covariant Wilson mass regulator. Its purpose is to
couple the spatial parity sectors of the central-difference envelope.
The values used are \(m=8\), \(\eta=0.01\), \(r=0.5\), with no measured
mass or magnetic moment as input. The charge convention is \(q=-1\)
when (1) is written with conventional \(p=-i\partial\) and
\(\alpha_{\rm physical}=-\gamma\). A positive \(A\) phase in (1)
then corresponds to \(p-qA\).

The matter kernel is the explicit construction choice

\[
K(A)=|D_W(A)|=\sqrt{D_W(A)^2}.
\tag{2}
\]

Its full number-conserving Fock extension has the empty vacuum at zero
energy. Both signed branches of \(D_W\) occur in (2). The construction
is distinct from the number-nonconserving and PW119 interacting sources;
no vacuum, particle identity or primitive-repair implementation is
transferred from those Hamiltonians.

The field is a real cochain on a larger periodic cube. It has mean zero
and satisfies \(D^*A=0\), where \(D\) and \(C\) are the actual finite
incidence and curl matrices. Its magnetic field is \(B=CA\). For a
normalized one-particle state, the stationary energy functional is

\[
\mathcal E(\psi,A)=\langle\psi,K(A)\psi\rangle
+\frac{a^3}{2\kappa}\|CA\|^2,\qquad \kappa>0.
\tag{3}
\]

The field couples to the expectation of the actual matter current.
This is a semiclassical transverse Maxwell architecture. Longitudinal
electrostatic energy and feedback are excluded. No screening or
neutralization mechanism, full coupled quantum electromagnetism, or
free self-bound particle is inferred. The external confining well and
its mechanical support are not dynamical source factors.

## Exact positivity and finite stationary existence

For every real link potential the Wilson term is nonnegative. Its
quadratic form is a sum of \(|u_x-e^{iaA_e}u_y|^2\), with the missing
boundary neighbors supplying positive Dirichlet terms. Thus \(M(A)\ge m\).
The kinetic term anticommutes with \(\beta\), so

\[
\operatorname{Re}\langle\beta v,D_Wv\rangle
=\langle v,M(A)v\rangle\ge m\|v\|^2.
\]

Cauchy–Schwarz gives \(\|D_Wv\|\ge m\|v\|\), hence
\(K(A)\ge m\) for every \(A\). There is no zero-crossing ambiguity in
the square root or its derivative. The empty Fock vacuum is the unique
matter ground at fixed \(A\). Together with the gauge-fixed field its
zero-field energy is zero.

On the finite mean-zero transverse space,
\(\|CA\|^2\ge\lambda_1\|A\|^2\), with
\(\lambda_1=4\sin^2(\pi/N)/a^2>0\). The normalized matter sphere is
compact. The energy (3) is continuous, and its field term is coercive.
A global one-particle minimizer therefore exists. For example,

\[
m\le\inf\mathcal E\le
m+\eta\max|x|^2+6r/a+3/a=:E_{\rm trial}^{\rm upper},
\]

and every minimizing sequence can be restricted to

\[
\|A\|^2\le
\frac{2\kappa(E_{\rm trial}^{\rm upper}-m)}{a^3\lambda_1}.
\]

Differentiability supplies the stationary matter and Maxwell equations.
This exact existence statement concerns some minimizer of (3). It does
not certify that the displayed numerical branch is that minimizer,
prove its nonzero dipole with interval bounds, or prove dynamical stability.

## Current from the complete positive kernel

Let \(X=\delta K\). Differentiating \(K^2=D_W^2\) gives the exact
Sylvester equation

\[
KX+XK=D_W\delta D_W+\delta D_W D_W.
\tag{4}
\]

In a complete \(D_W\) eigenbasis its unique solution is

\[
X_{ij}=\frac{\lambda_i+\lambda_j}
{|\lambda_i|+|\lambda_j|}(\delta D_W)_{ij}.
\tag{5}
\]

This formula handles both signs and all degeneracies. On the positive
spectral subspace \(P_+XP_+=P_+\delta D_WP_+\); on the negative subspace
\(P_-XP_-=-P_-\delta D_WP_-\). The executable checks (4)–(5) against
independent finite differences of the complete matrix absolute value,
including a rotated degenerate eigenspace.

For the selected positive-\(D_W\) stationary state, define
\(J_e=-\langle\psi,\partial K/\partial A_e\,\psi\rangle\). With
\(u=e^{iaA_e}\), its actual bond value is

\[
J_{xj}=\operatorname{Re}(\psi_x^\dagger\gamma_j u\psi_{x+e_j})
-r\operatorname{Im}(\psi_x^\dagger\beta u\psi_{x+e_j}).
\tag{6}
\]

The second term is the regulator current. Omitting it fails the independent
Hamiltonian derivative test. No Pauli coefficient or \(g\) parameter is
inserted. A simultaneous transformation
\(A\mapsto A+D\chi\), \(\psi_x\mapsto e^{-i\chi_x}\psi_x\) preserves
the spectrum, current and magnetic field. Changing the field while holding
the state fixed fails this test.

Gauge covariance and eigenstate stationarity imply \(D^*J=0\).
A constant vector potential is pure gauge on the open matter box, so
\(\sum_eJ_{ej}=0\) in each direction. The current has no periodic
harmonic component. The field equation is consequently

\[
C^*CA=\frac\kappa{a^3}J,
\qquad A=\frac\kappa{a^3}(D^*D)^+J
\tag{7}
\]

componentwise in transverse gauge. The Fourier inversion uses the exact
finite symbol \(4a^{-2}\sum_j\sin^2(k_ja/2)\). The current is recomputed
from the eigenstate at every iteration; it is not a prescribed source.

## Numerical stationary family and controls

The solver alternates the full bound-state eigenproblem and (7), with
relaxation. Every original spatial mode occurs in the matrix. The exact
sixteen-row unitary reduction supplies four identical complete blocks;
the prepared packet occupies one block and the other Fock modes are
retained. The containing periodic field cube and open matter boundary
are specified separately. The matter Dirichlet ghost boundary is fixed
at coordinate magnitude three in the displayed refinement family.

At \(\kappa=0.1\), the three matter grids have \((n,a)=(5,1),(7,0.75),
(9,0.6)\), with field side lengths \(19,25,31\) in grid sites. Their
magnetic coefficients are

| Grid | Coefficient using supplied \(m\) | Coefficient using actual bound energy |
| --- | ---: | ---: |
| 5 | 1.42774518 | 1.47516960 |
| 7 | 1.63861960 | 1.68416493 |
| 9 | 1.74608113 | 1.78875065 |

These are finite-grid results. The monotone increase does not establish
a continuum value or a precision \(g=2\) prediction. The containing
periods \(Na\) are retained explicitly; a larger-field-box calculation
separately tests boundary sensitivity.

For the coarse grid, \(\langle S_z\rangle\simeq0.49879995\),
\(\mu_z\simeq-0.04450995\), the positive matter eigenvalue is
\(8.26572912\), and the magnetic field energy is
\(1.34608\,10^{-6}\). The empty-vacuum-relative bound excitation plus
field energy is \(E_{\rm bound}\simeq8.26573047\). It is distinct from
the supplied mass coefficient eight and from the matter phase frequency.
The reported coefficients use

\[
\boldsymbol\mu=\frac12\sum_e x_{e,\rm midpoint}\times J_e,
\qquad
g_m=\frac{2m\mu_z}{q\langle S_z\rangle},\quad
g_{\rm bound}=\frac{2E_{\rm bound}\mu_z}{q\langle S_z\rangle},
\quad q=-1.
\tag{8}
\]

There is no factor \(a^3\) in the sum: \(J_e\) in (6) is the derivative
of the total matter energy with respect to a link potential. Its density
normalization in Maxwell's equation is \(J_e/a^3\).

The numerical eigenstate residuals are below \(10^{-10}\), the fixed-point
residuals below \(2\,10^{-11}\), and the Maxwell residuals below
\(10^{-9}\). Current divergence, net current and gauge divergence are
below \(10^{-10}\). These are floating-point residual diagnostics,
not interval error bounds. Reversing the selected spin reverses the
current and generated field. A perturbed initial field and a different
relaxation recover the same branch. Changing \(\kappa\) changes the
current itself, demonstrating feedback. Convergence of this iteration
is not a proof of physical stability.

The central-difference source with \(r=0\) supplies a useful distinct
control: its selected stationary bound doublet has current below
\(10^{-10}\). Thus an externally probed Clifford spin response alone
does not establish a stationary intrinsic dipole. The nonzero result
uses the explicitly stated Wilson envelope architecture.

## Dipole field and weak-field response

The generated field is compared with the infinite-space axial expression
\(B_z=\kappa\mu_z/(2\pi r^3)\), using the independently computed moment
in (8). No dipole amplitude is fitted to the field. At coarse-grid radii
four, five and six, the ratios are approximately 0.937, 0.963 and 0.961.
Finite source extent, spacing and periodic images are retained; an
infinite-space or continuum field identity is not asserted by these ratios.

An independent uniform weak-field probe uses
\(A_{\rm ext}=\tfrac12 B_{\rm ext}\hat z\times x\) on the matter links.
The matter state and self-field are reoptimized at both signs of the probe.
Stationarity of (3) implies
\(d\mathcal E/dB_{\rm ext}=-\mu_z\). The numerical derivative agrees
with the actual-current moment to about \(5.3\,10^{-10}\) on the coarse
grid. This tests the same bound construction and energy reference.
No coupled time-dependent Larmor experiment, radiative anomaly,
physical electron mass, or measured electromagnetic coupling is derived.
