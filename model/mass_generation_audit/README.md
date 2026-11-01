# Interaction-generated charge gap in the unchanged Weyl chain

The existing repulsive contact interaction generates a positive thermodynamic
charge gap in the **half-filled canonical one-dimensional ground sector**. No
Dirac mass, measured particle mass, or new hopping coefficient is inserted.
At the existing source values `c=1/8,U=1`, the gap is

\[
\Delta_c=0.7715517157380375\ldots .
\]

This is a material join beyond `common_source/mass_response`: that earlier
alternative Hamiltonian assigned its hopping and repulsion from logarithmic
reference statistics. Here the original `periodic_momentum` Weyl/contact
operator itself is transformed exactly. The positive result concerns a
specified many-body reference in one spatial dimension. It is not a massive
three-dimensional particle, a globally gapped spin-and-charge vacuum, or a
demonstration that the complete repair process prepares or preserves that
reference.

## Exact unchanged-law join

Take a periodic source geometry `(L,1,1)`, `L=4m`, with the original link
`(ic/2) sigma_x` and contact `U n_up n_down`. The onsite SU(2) frame

\[
G_x=(-1)^x(i\sigma_x)^x,\qquad
G_x(i\sigma_x)G_{x+1}^\dagger=-I
\]

has period four. Its CAR lift fixes each local vacuum and doublon, so it
preserves the interaction, total number and all canonical energies exactly.
The transformed Hamiltonian is the standard repulsive Hubbard chain

\[
H=-t\sum_{x,s}(d^\dagger_{x,s}d_{x+1,s}+\mathrm{h.c.})
  +U\sum_x n_{x,\uparrow}n_{x,\downarrow},\qquad t=c/2=1/16.
\]

`chain.py` verifies the four local identities with Gaussian integers, binds
the original matrices on lengths 4,8,12,16, and independently checks the full
256-dimensional Fock-space identity on the four-site ring. The local period
and link proof establishes the whole `L=4m` family; matrices for every length
are not enumerated. At `U=0` its four exact zero modes and chiral pairing give
zero half-filled charge gap. Contact vanishes identically in the one-particle
sector, so this mechanism does not gap the previously measured Weyl modes.

## Imported theorem and exact numerical enclosure

The mathematical input is the thermodynamic repulsive Hubbard-chain charge
gap, rather than an unproved claim about a finite diagonalization. The
[Lieb–Wu analysis, §7](https://arxiv.org/html/cond-mat/0207529) gives
`mu_+ + mu_- = U` and the charge discontinuity

\[
\Delta_c=\lim_{L\to\infty}
 [E_L(L+1)+E_L(L-1)-2E_L(L)]
=U-4t+8t\int_0^\infty
 \frac{J_1(x)}{x(1+e^{Ux/(2t)})}\,dx.
\]

Expanding the Fermi denominator and using the Laplace transform of `J1/x`
gives, for `r=2t/U<1`,

\[
\Delta_c=U-4t+8t\sum_{k\ge1}
 {1/2\choose k}r^{2k-1}\eta(2k-1).
\]

The `k=1` sum is the alternating harmonic series `log 2`; all `k>=2`
double sums converge absolutely. The binomial series converges uniformly
for `r/m<=r<1`. After `K` terms its absolute tail is at most

\[
8t\left|{1/2\choose K+1}\right|
 \frac{r^{2K+1}}{1-r^2}.
\]

`gap.py` uses `K=12` and exact rational Euler sums for Dirichlet eta. Their
nonnegative finite-difference coefficients lie in `[0,1]` by the integral
of `exp(-x)(1-exp(-x))^n x^(s-1)/Gamma(s)`; the omitted Euler tail after
96 terms is at most `2^-96`. The resulting exact interval has width below
`1.67e-25`. Floating quadrature is only an independent test, not a premise
of the enclosure. Existing rational elementary-function utilities are
reused; the old alternative Hamiltonian's couplings are not consumed.

| Quantity in original model units | Value | Reference convention |
|---|---:|---|
| Charge discontinuity `Delta_c` | 0.7715517157380375… | invariant under common `-mu N` |
| Addition cost `E(L+1)-E(L)` | 0.8857758578690188… | original `H`, `mu=0` |
| Removal cost `E(L-1)-E(L)` | −0.1142241421309812… | original `H`, `mu=0` |
| Symmetric charge threshold `Delta_c/2` | 0.3857758578690188… | declared `K=H-(U/2)N` |

The negative original removal cost is essential: the half-filled reference
is not the grand-canonical ground of the original `mu=0` action. Specifying
particle number or a reservoir centered at `U/2` is an additional reference
condition. The spin sector remains gapless. Neither a gap for all observables
nor a Lorentz rest mass follows from this charge threshold alone. No source
word preparing the thermodynamic ground, repair-invariance theorem, or
measurement of this threshold is asserted here.

## Fixed-law continuum obstruction and explicit running condition

The positive-integral form of the same charge gap is recorded in
[Jeckelmann–Gebhard–Essler](https://arxiv.org/html/cond-mat/9911281):

\[
\Delta_c=\frac{16t^2}{U}\int_1^\infty
 \frac{\sqrt{y^2-1}}{\sinh(Ay)}\,dy,
\qquad A=2\pi t/U.
\]

Writing it as a positive odd-index sum of `K1` and using
`1<=sqrt(1+x)<=1+x/2` in the integral representation of `K1` proves the
global bounds

\[
\Delta_0\le\Delta_c\le
\Delta_0\frac{1+3/(8A)}{1-e^{-2A}},\qquad
\Delta_0=\frac8\pi\sqrt{Ut}\,e^{-A}.
\]

For the conventional fixed-velocity refinement `t=v/(2a)`, fixed `U>0`,
the gap therefore vanishes as `a^(-1/2) exp[-pi v/(aU)]`. The nonzero
lattice charge gap does not survive this fixed-interaction continuum limit.
`scaling.py` supplies rational enclosures through `a=1/256`, where the
upper bound is below `1e-40` for `v=1/8,U=1`.

A precise additional coupling-family condition can retain a charge threshold:
given a **chosen** `M>0,v>0`, set

\[
a(A)=\frac{4v e^{-A}}{M\sqrt{2\pi A}},\quad
t(A)=\frac{v}{2a(A)},\quad
U(A)=\frac{\pi v}{a(A)A}.
\]

Then `Delta_0/2=M` exactly and
`1<=Delta_c/(2M)<=(1+3/(8A))/(1-exp(-2A))->1`.
This derives a sufficient running-coupling condition, not a value of `M`
from the fixed source. Establishing the full massive relativistic scaling
theory, its source-prepared reference, and any three-dimensional continuation
requires further work. No SI calibration or observed-species fit is used.

## Verification and scope

Run `.venv/bin/python -W error -m pytest model/mass_generation_audit -q`
and `.venv/bin/python -W error -m model.mass_generation_audit --check`.
The cold receipt binds this proof, code, tests, imported source manifests
and the focused paper fragment. Its exact intervals are rational strings;
displayed decimals are informational. The imported thermodynamic theorem
is identified explicitly and is not presented as a new machine-checked
proof of Bethe ansatz completeness.
