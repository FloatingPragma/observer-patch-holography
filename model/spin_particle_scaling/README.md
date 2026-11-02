# Spatial quantum numbers and scaling of the ordered spin sector

This package resolves the source spatial quantum numbers of the positive
gapless spin packets proved in [FA92](../spin_goldstone_sector/README.md).
It constructs two orthogonal, equal-weight polarization channels and
nonvanishing normalized energy--momentum scaling measures. It also closes
the shrinking-weight problem in transferring these readouts to the actual
source: a single explicit volume/time/coupling diagonal makes the dynamical
errors small relative to the packet weight.

The result concerns the selected ordered ground of the fixed-spacing
antiferromagnet. It neither assumes nor proves a sharp particle pole,
relativistic helicity, a Lorentz-invariant field limit, or scattering.
The scale limit below is a limit of excitation probability measures;
it is not a construction of a continuum field algebra.

## Inputs and conventions

Use the unchanged source frame `d_x=R_x c_x`,
`R_x=(i sigma_1)^x1(i sigma_2)^x2(i sigma_3)^x3` from
[field_symmetry](../field_symmetry/README.md). The original Weyl-frame spin
is `s_x^a=c_x^* sigma_a c_x/2`; the hidden spin is
`S_x^a=d_x^* sigma_a d_x/2`. The identities here hold on the full onsite
four-dimensional CAR space, including the empty and doubly occupied states.
On the one-fermion-per-site carrier, [FA90](../protected_spin_limit/README.md)
gives `H=J sum_<xy>(S_x.S_y-1/4)`, `J>0`, and its all-record source comparison.
FA92 constructs a zero-field infinite-volume ground `omega` selected by a
positive staggered field along internal axis 1. It has staggered
magnetization `m_s>m=1/1000`.

The external mathematical input is precisely FA92's matched DLS/Koma
input. In particular [Koma, arXiv:2105.04970v2](https://arxiv.org/html/2105.04970v2),
equations (2.6)--(2.9), (3.4), (4.1), (5.3)--(5.12), supplies the
field-selected state, smooth packet limit and transverse susceptibility
bound. We use the explicit spectral-window constants already proved in
FA92; no additional assumption of a spin-wave dispersion is introduced.

## Exact source spatial action

Conjugation by the onsite frame gives

```
S_x^a = chi_a(x) s_x^a,     chi_a(x)=(-1)^(sum_j x_j-x_a).
```

For an actual proper cubic rotation `R` and integer translation `t`, write
the active position map as `x=R y+t`. The Weyl spin transforms as the
ordinary vector `s'(x)=R s(y)`. In the hidden frame the exact action is

```
S'(x)=Q_t R S(y),
Q_t=diag((-1)^(t2+t3),(-1)^(t1+t3),(-1)^(t1+t2)).
```

Thus the staggered order parameter transforms as
`n' = D_t R n`, with `D_t=diag((-1)^t1,(-1)^t2,(-1)^t3)`.
`algebra.py` checks every cubic rotation and translation parity, including
negative lattice coordinates, directly from the local frame identity.
The source Hamiltonian is covariant under this physical space group.
Finite repair maps outside the protected carrier are not asserted to have
this symmetry; their complete-record limit is the identity spin repair.

For `n=e1`, the unbroken physical translations are
`2 Z x Z x Z`. The unbroken affine operations obey `D_t R e1=e1`;
there are 32 such rotation/translation-parity pairs out of 192. Pure
site-centered rotations preserving the oriented axis form C4.
The original-spin order is `omega(s_x^1)=m_s(-1)^x1`, with the other
components zero. More generally the three order components have carriers
`X_a=pi e_a`, whereas a hidden near-zero component has carrier
`M_a=pi(1-e_a)`.

For the transverse components, reduction modulo the actual unbroken
translation lattice gives

| Hidden component | Original full-zone carrier | Residual-zone carrier |
| --- | --- | --- |
| `S^2` | `(pi,0,pi)` | `(0,0,pi)` |
| `S^3` | `(pi,pi,0)` | `(0,pi,0)` |

These are two distinct momentum neighborhoods, exchanged by the physical
quarter turn about axis 1. A rotation about a cube center has an additional
translation. Consequently this site-centered action does not replace or
contradict [FA95](../cubic_particle_channels/README.md), whose centered
finite cube has a different spatial little group action.

## Two nonzero orthogonal polarization channels

The selected state is invariant under hidden rotations about axis 1,
under all even-parity hidden translations, and under a one-step hidden
translation combined with a pi spin rotation about axis 2. These follow
already at positive field from invariance of the finite Hamiltonian and
its unique ground, and pass to the stated limits.

Let `A_a(f)=sum_x f_x S_x^a`, `a=2,3`, with Fourier support in
`p/2<=|k|<=p<=1`. This support and its translate by `Q=(pi,pi,pi)`
are disjoint. For any bounded positive function of excitation energy,
the two-by-two spectral matrix of these operators is a scalar multiple
of the identity. Here is a useful explicit proof. Even-parity translation
invariance only permits Fourier pairs `k'=k` or `k'=k+Q`. The support
condition removes the second possibility. Hidden axis-1 rotations give
equal diagonal entries. The combined one-step translation and axis-2
pi rotation preserves the state and energy, but negates the off-diagonal
entry at the surviving equal momentum; its translation phases cancel.
That entry is therefore zero. The proof holds for spectral projections
and hence for the full spectral measures, not only their first moments.

The circular operators `A_+=(A_2+i A_3)/sqrt(2)` and
`A_-=(A_2-i A_3)/sqrt(2)` carry residual internal U(1) charges `+1,-1`.
Their filtered vectors are orthogonal and have exactly the same nonzero
energy measure. A U(1)-invariant state alone would not imply equality:
the one-spin polarized countercontrol has one zero and one nonzero
ladder norm. The combined translation symmetry above is essential.

Choose a Fourier envelope invariant under a site-centered C4 turn.
That rotation sends `(A_2,A_3)` to `(A_3,-A_2)` and the circular channels
have characters `-i,+i`. They occupy the two-valley momentum star;
each circular channel is generally not a single physical momentum
eigenchannel. The real components give the separate translation channels.
These are exact cubic and internal quantum numbers, not relativistic
helicity. In particular this construction does not identify a photon.

## Nonvanishing normalized low-energy scaling measures

Fix a smooth radial function `phi`, supported in `1/2<|q|<1`, normalized
by `integral |phi(q)|^2 dq/(2pi)^3=1`. For example, normalize
`exp(-1/(|q|^2-1/4)-1/(1-|q|^2))` on that annulus and set it to zero
outside. Put `fhat_p(k)=p^(-3/2) phi(k/p)`.
Choose one smooth energy function `G`, equal to one on
`[m/8,32/m^2]`, bounded between zero and one, and supported in
`(m/16,64/m^2)`. The physical energy filter is `G(E/(Jp))`.
Write `Psi_(a,p)` for the resulting GNS vector and `w_p` for its norm squared.
The two-channel result and FA92 give, for either real or circular channel,

```
c p <= w_p <= p/4,       c=m^2/96,
energy support:  m J p/16 < E < 64 J p/m^2,
mean E <= 48 J p/m^2.
```

The upper weight bound also follows directly from FA92: the total
inverse moment is at most `1/(10J)` because `E(k+Q)>=5`, while the
first moment is at most `Jp^2/2`. Spectral Cauchy--Schwarz gives total
weight at most `p/sqrt(20)<p/4`; the filter is a contraction. Thus the
natural rescaling of the two-point spectral measure by `1/p` already
has a uniformly positive finite mass, without a fitted normalization.

For each real channel the commuting energy and unbroken physical
translations have a joint spectral measure. Recenter its momentum at
the corresponding residual-zone carrier and push forward by
`(E,k)->(E/(Jp),k/p)`, then divide by `w_p`. This is a probability measure
supported on the fixed compact set

```
[m/16,64/m^2] x {q:1/2<=|q|<=1}.
```

Every sequence `p->0` therefore has a weakly convergent subsequence
with total mass one, positive scaled energy and nonzero scaled momentum.
The measures rescaled by `1/p` instead have subsequential masses in
`[m^2/96,1/4]`, with a scalar two-polarization spectral matrix.
The two real channels give the same recentered measure, by internal
rotation symmetry. For the radial envelope it is C4 invariant in q.
Both orthogonal circular channels retain the same scalar energy measure.
This is a genuine two-channel scaling sector; it does not assert uniqueness
of the measure, support on one dispersion surface, surviving interactions
among normalized packets, or a continuum canonical field algebra.

## A source diagonal that preserves the shrinking weight

This section closes a further input beyond fixed-observable source
convergence. Put `p=1/n` and let `F` be the inverse Fourier transform of
`phi`. Define finite constants
`C_r=sup_x (1+||x||_infinity)^r |F(x)|` and
`D_r=integral (1+|u|)^r |g(u)| du`, where `g` is the inverse transform
of G. All are finite since the templates are smooth and compactly supported.
Since `f_p(x)=p^(3/2) F(px)`, cube integration gives the explicit bound

```
sum_x (1+p||x||_infinity)^(-8)
 <= (3/2)^8 p^(-3) integral_R3 (1+||y||_infinity)^(-8) dy
 = (6561/1120) p^(-3) < 6 p^(-3).
```

Indeed on the unit cube about x, `1+p||y||<=3(1+p||x||)/2`.
The integral is `24 integral_0^infinity r^2/(1+r)^8 dr=8/35`.
Consequently the filtered spin operator has norm at most
`K n^(3/2)`, `K=3 C_8 D_0`; the factor 1/2 is the onsite spin norm.
The time filter is `g_p(t)=Jp g(Jpt)`, so its L1 norm is D0,
independent of p. No uniform operator-norm assumption is made.

Truncate space to radius `n^2` and filter times to `n^2/J`.
The omitted operator norm is bounded by

```
K_tail n^(3/2)/(1+n)^8,
K_tail=3(C_16 D_0+C_8 D_8).
```

For quadratic readouts the corresponding error divided by `w_p>=c/n`
is `O(n^4/(1+n)^8)=O(n^-4)`, with constants fixed by these templates.

Take the source box radius `L_n=128 n^3`, with open-cube edge count
`E_n=6 L_n(2L_n+1)^2`, and keep FA90's couplings
`c_source=J E_n^4`, `U_source=J E_n^8`. A scaled characteristic-function
test on any fixed compact time/translation set uses, for large enough n,
total evolution time at most `3 n^2/J` and support radius at most `2n^2`.
The extra radius includes the translated packet. Rounding translations
to the unbroken lattice gives a vanishing `O(1/n)` phase error in the
rescaled momentum variable; the known valley phase is removed exactly.

FA90's boundary bound, applied to each local spin before summing the
absolute coefficients, is at most

```
B_n = 18 n^2(2L_n+1)^2 2^(-128 n^3+218 n^2).
```

Here `exp(48Jt)<=2^(72Jt)` uses `log 2>=2/3`. For `n>=2` the exponent
is at most `-19 n^3`. Thus the boundary error is also negligible after
the `O(n^4)` normalization cost. The source complete-record bound is

```
delta_n = 24/E_n^3 + 24 (Gamma/J) n^2/E_n^2 + 18 n^2/E_n = O(n^-7).
```

Weighted two-point/characteristic readouts incur at most
`K^2 n^3 delta_n`; division by the true weight costs one more n.
Hence the normalized source contribution is `O(n^4 delta_n)=O(n^-3)`.
One can justify the readout bound by expanding the two packet insertions
into onsite spin operators and the time integrals: each local bounded
correlation has the FA90 all-input, all-record comparison, and the sum
of the absolute coefficient products is at most `K^2 n^3`.
FA90's Stinespring estimate is uniform over untouched references, so
controlled local contractions or unitaries can be inserted and the free
segments telescoped. The first term `24/E_n^3` includes up to three free
segments; the remaining terms already use their total time. Equivalently
the insertions can be tested with controlled local unitaries. No collision
outcome is discarded.

Normalization is admitted only after the complete error epsilon is smaller
than `c/n`. If numerator and denominator errors are `epsilon_N,epsilon_D`,
the normalized characteristic error is bounded by
`(epsilon_N+epsilon_D)/(c/n-epsilon_D)`. This rejects a small absolute
error that is still large relative to the detected weight.

The finite initial spin state in this quantitative dynamical comparison
is the restriction of omega to the source box. Such density matrices
exist and embed in the one-fermion-per-site carrier. Their finite-cost
preparation, or an effective convergence rate of finite-field finite-volume
ground states to these restrictions, is **not** supplied by this theorem.
Computable finite-volume ground preparations may approximate each finite
restriction in the established iterated limits; an additional state error
must satisfy the same `o(n^-4)` accuracy requirement before normalization.
Thus the result is a controlled dynamical/readout diagonal, not an
unproved effective preparation theorem. Finite source correlations need
not themselves be stationary ground spectral measures; they converge to
the stated positive limiting characteristic functions.

## Verification and limits

`algebra.py` proves the exact local frame, affine space-group and C4
identities. `bounds.py` checks the rational lattice-sum constant, shrinking
weight normalization and explicit volume/coupling schedule. `experiment.py`
independently checks all energy clusters of the full 256-state positive-field
periodic-cube control: scalar transverse spectral matrix, circular charges,
orthogonality and actual spatial actions. This diagnostic retains its 24
parallel periodic bonds and does not prove an infinite-volume theorem.
Eigensolver residual and orthogonality checks reject corrupted proposals.
Raw floating tail weights are not serialized as exact claims; the receipt
records the checked eigenspace identities and significant positive response.

The receipt hashes a narrow, explicit list of mathematical inputs and
reports exact local snapshots of inherited receipts from fixed checkpoint
`ba312dbf0` separately. Their original hashes are checked before parsing;
runtime access to Git history is unnecessary, including on a shallow clone.
Live broad CI receipts are not direct inputs. It does not refresh
or silently repair previous evidence. Reproduce with

```sh
python -W error -m pytest model/spin_particle_scaling -q
python -m model.spin_particle_scaling --write
python -m model.spin_particle_scaling --check
```

The remaining particle inputs include a sharp pole or controlled scattering
description, a full relativistic interacting field scaling limit, physical
spin/helicity, measured species and empirical energy calibration. None is
assigned by the lattice C4 characters or by the nonzero probability limits.
