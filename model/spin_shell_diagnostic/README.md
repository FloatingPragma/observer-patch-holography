# Exact diagnostics for the protected-spin mass shell

The protected spin sector has nonzero low-energy spectral weight and two
polarization channels. A sharp mass shell has **not** been established.
This package derives the exact current and three-spin force whose ground-state
correlations would establish it. It also supplies a quantitative shell test and
checks why the current moment inequalities do not settle that test.

The Hamiltonian is unchanged:

\[
 H=J\sum_{\langle xy\rangle}(\mathbf S_x\cdot\mathbf S_y-1/4),
 \quad J>0,\quad S=1/2,\quad x\in\mathbb Z^3.
\]

The state is the field-selected ordered ground from
[the ordered-spin proof](../spin_goldstone_sector/README.md). Its positive GNS
generator is K. Momenta below mean hidden momenta recentered at the two physical
valleys specified by [the scaling construction](../spin_particle_scaling/README.md).
The internal spin and the source space-group labels are retained; neither is
silently replaced by Lorentz spin.

## What the primary theorem supplies

Koma's [primary paper](https://arxiv.org/html/2105.04970v2), equations
(2.11)--(2.16) and Theorem 2.1, defines the excitation energy as the mean energy
of an energy-filtered packet. It proves positive lower and upper multiples of
its momentum scale. Its proof uses first moments, inverse moments and a Ward
identity. It does not prove that the width divided by the momentum tends to
zero. Our prior construction implements that mean-energy theorem with explicit
constants and a fixed broad window; it does not already contain a sharp pole.

We retain that window. Set m=1/1000, alpha=m/16 and beta=64/m^2. For a smooth
packet with p/2<=|k|<=p<=1, unit Fourier norm
integral |fhat(k)|^2 dk/(2pi)^3=1 (equivalently sum_x |f_x|^2=1),
and the prior fixed smooth filter G, write

\[
 \Psi_p=G(K/(Jp))A_f\Omega,\qquad w_p=\|\Psi_p\|^2,
 \qquad {m^2p\over96}\le w_p\le {p\over4}.
\]

Its energy support is alpha Jp<=K<=beta Jp. No filter that narrows around a
proposed dispersion surface is introduced here. A shell must capture a
nonvanishing amount of the original response in this fixed broad window.

## A shell test that handles radial wave packets

Let e(k)=sum_j(1-cos k_j) and s(k)=sqrt(2e(k)). On the joint support of the
normalized vector psi_p=Psi_p/sqrt(w_p), define

\[
 X_p={K\over J s(P)},\qquad
 d_p=\langle X_p\rangle\langle X_p^{-1}\rangle-1.
\]

Even-translation invariance gives the commuting spectral representation on
the small recentered momentum support; no global logarithm of lattice
translations is required. Since (2/pi)|k|<=s(k)<=|k|, we have
alpha<=X_p<=pi beta. If X and Y are independent random variables with this
spectral law, direct multiplication gives

\[
 d_p={1\over2}\mathbb E{(X-Y)^2\over XY},\qquad
 \alpha^2 d_p\le\operatorname{Var}(X_p)\le(\pi\beta)^2d_p.
\]

Thus d_p tends to zero exactly when this ratio distribution concentrates.
For one common c>0, valid for every admitted momentum envelope and direction,

\[
 \left\|\left({K\over Jp}-c{|P|\over p}\right)\psi_p\right\|
 \le \pi\beta\sqrt{d_p}+|\langle X_p\rangle-c|+{cp^2\over12}.
\]

The last term follows from
0<=|k|^2-2e(k)<=|k|^4/12 and |k|+s(k)>=|k|.
Chebyshev's inequality bounds the normalized spectral weight outside an
eta-neighborhood of the rescaled cone by the square of this expression
 divided by eta^2. Conversely convergence to that cone in this norm implies
both ratio-moment conditions, since s(k)/p>=1/pi.

This equivalence concerns concentration of the **entire selected response**.
It is sufficient, but is not necessary for the existence of any sharp branch:
a branch may coexist with a broad multiparticle continuum in the same window.
For example, half an atom plus half a uniform distribution on[1,2] has a
nonzero atom and a strictly positive ratio-moment defect.

An alternative sufficient route is to construct an independently justified
spectral component (a positive submeasure of the actual joint response, or a
reducing projection from its symmetry/dynamics) with weight at least zeta p
for one zeta>0. One can then establish concentration or a dispersion atom for
that component while retaining the continuum in the complement. The ratio
and residual tests below apply to the normalized component as well; for a
reducing projection Q, apply Q also to the current/force residual. Merely
choosing an artificially shrinking energy filter proves neither the required
nonvanishing weight nor the existence of such a branch.

More concretely, it suffices to prove that tubes of width delta_p tending to
zero around one fixed rescaled cone contain a fraction at least zeta>0 of the
original broad-window weight. For any fixed closed tube of width delta>0,
the shrinking tubes eventually lie inside it. The closed-set inequality for
weak convergence therefore gives limiting mass at least zeta in every such
fixed tube; continuity from above gives mass at least zeta on the cone itself.
This permits a continuum in the remaining weight. The fixed-fraction lower
bound is the substantive requirement; normalizing a narrow filter cannot
supply it.

The speed is v=Jc in lattice units. A separate fitted c at each momentum,
envelope or system size does not meet this condition. Raw energy variance
is not an adequate test for a radially broad packet: even a perfect cone
has different energies at different |k|. The ratio X removes that spread.
A sharp scaling cone would still not, by itself, prove an isolated pole at
every nonzero lattice momentum, a Lorentz field algebra or scattering.

## Exact current and force identities

Choose one orientation of each bond, let A_f^a=sum_x f_x S_x^a, and let
(Lf)_x=sum_(y~x)(f_x-f_y). Exact spin commutators give

\[
 [H,A_f^a]=iJ\sum_{\langle xy\rangle}(f_x-f_y)
                  (\mathbf S_x\mathbin\times\mathbf S_y)^a.
\]

For spin 1/2, h_xy=J(P_xy-1)/2 in terms of the swap. Therefore
ad_(h_xy)^2 S_x^a=J^2(S_x^a-S_y^a)/2. Separating repeated and adjacent bonds,

\[
 [H,[H,A_f^a]]={J^2\over2}A_{Lf}^a+F_f^a,
\]

\[
 F_f^a=-J^2\sum_{\langle xy\rangle}(f_x-f_y)
 \left\{\sum_{z\sim x,z\ne y}[(\mathbf S_x\times\mathbf S_z)\times\mathbf S_y]^a
 +\sum_{z\sim y,z\ne x}[\mathbf S_x\times(\mathbf S_y\times\mathbf S_z)]^a\right\}.
\]

These are operator identities, not replacements of spins by their means.
`algebra.py` checks arbitrary complex coefficients and all three axes on
single bonds, both orientations of a two-bond path, and a triangle. Disjoint
bonds commute, so these local checks exhaust the general identity. On the
cubic lattice L has Fourier symbol 2e(k), whence

\[
 K^2 A_f\Omega=J^2 e(P)A_f\Omega+F_f\Omega.
\]

A candidate speed must include the possible leading force renormalization.
For one c>0 define

\[
 R_p(c)=G(K/(Jp))\{F_f+(1-2c^2)J^2 e(P)A_f\}\Omega.
\]

It equals (K^2-c^2J^2s(P)^2)Psi_p exactly. On the broad window,

\[
 {\|R_p(c)\|\over(\beta+c)J^2p^2\sqrt{w_p}}
 \le {\|(K-cJs(P))\Psi_p\|\over Jp\sqrt{w_p}}
 \le {\|R_p(c)\|\over\alpha J^2p^2\sqrt{w_p}}.
\]

Thus a ground-specific filtered force bound
||R_p(c)||^2=o(J^4 p^5), uniformly over the admitted envelopes, is sufficient
and equivalent to concentration of this entire selected response on the
scaling shell. It is not necessary for a nonzero sharp component coexisting
with a continuum. The first-order form is an equivalent filtered current residual,

\[
 G(K/(Jp))\{[H,A_f]-cJ A_{sf}\}\Omega,
\]

whose squared norm must be o(J^2 p^3). It involves current-current four-spin
correlations and current-spin cross correlations. The second-order form
exposes the corresponding six-spin force correlations. Neither of these
vanishing bounds follows from the existing separate first/inverse moment
upper bounds. Setting F to zero would also prescribe an unjustified bare
speed J/sqrt(2).

## Why an operator-norm shortcut fails

There is an exact useful obstruction to a state-independent bound. Group the
force by unordered length-two paths x--y--z. Put u=f_x-f_y and v=f_z-f_y.
The path contribution divided by J^2 is

\[
 v S_x^a(\mathbf S_y\cdot\mathbf S_z)
 +u S_z^a(\mathbf S_x\cdot\mathbf S_y)
 -(u+v)S_y^a(\mathbf S_x\cdot\mathbf S_z).
\]

Here the superscript a denotes the fixed spin component. All terms
with three identical components cancel. The remaining six Pauli words are
orthogonal in normalized trace, with squared norms 1/64. The path squared
norm is consequently J^4(|u|^2+|v|^2+|u+v|^2)/32. Different three-site supports
are trace-orthogonal; on a simple triangle-free graph each path support has
one center. For degree D, the neighbor-pair sum at a site is

\[
 (2D-3)\sum_d|a_d|^2+\left|\sum_d a_d\right|^2.
\]

For a normalized plane wave on an even cubic torus of side at least four,
D=6, sum_d|exp(ik.d)-1|^2=4e(k), and sum_d(exp(ik.d)-1)=-2e(k). Hence

\[
 \tau(F_k^*F_k)={J^4\over8}e(k)(9+e(k)).
\]

The force is trace-orthogonal to every linear spin observable. For any real
lambda its linear subtraction instead gives

\[
 \tau\bigl(|F_k+(1-\lambda)J^2e(k)A_k|^2\bigr)
 =J^4\left\{{e(k)(9+e(k))\over8}
              +{(1-\lambda)^2e(k)^2\over4}\right\}.
\]

Operator norm dominates normalized trace norm, so an unfiltered
state-independent small-force estimate cannot establish the desired
remainder. This is **not** a no-shell theorem: the tracial state is not the
ordered ground, and its estimate contains no ground-state low-energy filter.
The needed cancellation may hold only in that ground and energy window.

## The moment bounds allow a continuous response

A direct positive-measure control also distinguishes the proved inequalities
from a shell. Fix k=(p,0,0), 0<p<=1, and any 0<m<=1/2. Let rho be uniform
probability on z in[1,2], with energy E=Jpz. Take

\[
 d\mu_p={p\over10}d\rho,\qquad
 d\nu_p={5m^2\over2p}d\rho,\qquad
 dC_p={im\over2}d\rho.
\]

The two-by-two measure with diagonal mu,nu and off-diagonal C is positive
of rank one at each z, and 2 Im integral C=m realizes the Ward commutator.
Since 11p^2/24<=e(k)<=p^2/2 and e(k+Q)=6-e(k), all the prior first-moment,
inverse-moment and weight inequalities hold. The checks use
56/81<log(2)<56/81+1/540<7/10, obtained from the first two terms and geometric
remainder of 2 atanh(1/3). Nevertheless this response has no atoms. Its
scaled energy variance is 1/12; even its best constant squared-energy
residual is 34/45, and its ratio-moment defect exceeds 1/27.

This is a control on those moment and Ward constraints, not a construction
of a competing Heisenberg state or a contradiction of OPH. It proves that
using only those constraints cannot establish a shell.

## Actual finite-source control and readout scope

`experiment.py` computes the full 256-state, 12-bond open-cube Hamiltonian
with J=1 and its exact-source current matrices. All eigenvectors and every
spectral weight enter its first, second and inverse moments. The three
nontrivial character classes have positive ratio defects, approximately
0.008684, 0.006647 and 0.031631; independent current moments agree with the
spectral calculation. This numerical control rejects a premature one-line
replacement on the finite cube. It is not an interval eigenvalue proof.
Its momenta do not satisfy p<=1, so it provides no evidence for or against
the infinite-volume low-momentum shell. No new native pulse execution is
claimed.

The ratio X and X^-1 are bounded continuous functions on the fixed joint
window. Consequently the existing source comparison transfers their
moments and d_p along any convergent scaling subsequence, provided the same
initial-state accuracy and normalization conditions are met. It does not
force d_p to vanish. Bounded smooth extensions off the compact window allow
the same finite-history Fourier readouts. The existing absence of an
effective thermodynamic state-preparation modulus is not removed here.

Reproduce the new local checks and record with

```sh
python -W error -m pytest model/spin_shell_diagnostic -q
python -W error -m model.spin_shell_diagnostic --write
python -W error -m model.spin_shell_diagnostic --check
```

The record pins only its new sources and explicit stable mathematical inputs;
it does not change or replay historical receipts. It sets
`mass_shell_established=false`. The specified ground-state current/force
cancellation, or equivalently common-speed ratio-moment saturation, would
establish concentration of the chosen broad-window response. Alternatively,
a separately justified component with proved nonvanishing weight may carry
a sharp branch alongside a continuum. Neither sufficient route has been
established for this ground state. Measured species and relativistic masses
require further identification.
