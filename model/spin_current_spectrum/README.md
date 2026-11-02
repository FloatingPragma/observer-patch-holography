# The spin-current part of the energy sum rule

The unchanged ordered antiferromagnetic source obeys an exact continuity
identity. It identifies the first-energy tail of a long-wavelength spin
excitation with an inverse-energy tail of a bond-current excitation. A
nonzero fraction of this current response is forced into energies
proportional to momentum by FA105. All fractional energy moments of order
strictly between zero and one are uniformly integrable. Neither fact proves
uniform integrability of the first energy moment or a sharp mass shell.

This is a reduction to a particular four-spin correlation function in the
actual ground state, together with exact finite controls. It supplies no
phenomenological magnon decomposition and no numerical particle mass.

## 1. Source and exact current

Use exactly the FA100/FA105 source

```
H_J = J sum_<xy> (Sx.Sy-1/4),       J>0,
```

on the cubic lattice, with each unoriented bond counted once. The ground
state is the same zero-field state selected by an infinitesimal staggered
field along internal axis 1. Write its positive ground generator as K and
its order as m>1/1000. The internal spin index is not a Lorentz-spin label.
Put Q=(pi,pi,pi), epsilon_x=(-1)^(x1+x2+x3), and

```
j_(x,i)^a = J (Sx cross S_(x+ei))^a.
```

For a summable test f, direct spin commutation gives

```
[H_J, sum_x f_x Sx^a] = i sum_(x,i) (f_x-f_(x+ei)) j_(x,i)^a.       (1)
```

Equivalently, the Heisenberg time derivative of Sx is minus the outgoing
current divergence. This is an operator identity on every finite graph;
it is not an expectation or a ground-state approximation. Finite positive
staggered field changes (1). For H_B=H_J-B sum epsilon_x Sx^1 its extra term is

```
-i B sum_x,b f_x epsilon_x epsilon_(1ab) Sx^b.                     (2)
```

In particular the axis-2 torque is -iB times the staggered axis-3 field.
The FA98 finite experiment at B/J=1/4 cannot be substituted in a zero-field
continuity formula. The field is taken to zero at fixed packet before the
long-wavelength limit, as in the inherited ordered-state construction.

For finite periodic boxes use the 1/sqrt(V) Fourier normalization with
exp(-ik.x), and put the current Fourier phase at the bond midpoint:

```
j_i^a(k) = V^(-1/2) sum_x exp[-ik.(x+ei/2)] j_(x,i)^a,
e(k) = sum_i (1-cos(k_i)),
j_L^a(k) = sum_i sin(k_i/2) j_i^a(k) / sqrt(sum_i sin(k_i/2)^2).
```

Then at nonzero lattice momentum (1) reads

```
[H_J,S^a(k)] = -sqrt(2e(k)) j_L^a(k).                            (3)
```

All signs in (1)--(3), and the field torque, are checked independently as
exact finite matrices in `algebra.py`. The thermodynamic statements mean
identities between momentum-smeared vectors and their positive joint
spectral measures; they do not assume that a point-momentum vector belongs
to the infinite-volume Hilbert space. Transverse two-point functions are
translation invariant: the selected state preserves even translations and
an odd translation combined with the appropriate pi spin rotation, and
both factors in the two-point function acquire the same sign. The
reflection-positive transverse structure bound makes the momentum marginal
absolutely continuous away from its exceptional momenta. Disintegrated
formulas below hold for almost every nonexceptional k, or equivalently after
integration against packet tests supported away from those momenta.

## 2. An exact spectral sum rule

Let nu_S(k,dE) and nu_j(k,dE) be the positive energy measures of the transverse
spin and longitudinal bond current in this state. Equation (3) and K Omega=0
imply the equality of positive measures

```
E^2 nu_S(k,dE) = 2e(k) nu_j(k,dE).                              (4)
```

The longitudinal current has no zero-energy component at nonzero k, since
its vector is in the range of K there. FA105 gives the exact first moment

```
integral E nu_S(k,dE) = J kappa e(k),
kappa = -omega(S0^1 Sei^1 + S0^3 Sei^3),
m^2/2 <= kappa <= 1/2.
```

The first-moment identity extends from FA105's real tests to complex
Fourier tests without replacing a mode by a symmetrized pair. Each selected
finite ground and its Hamiltonian are real in the S3 basis. For transverse
axis 3 both spin matrices are real; for axis 2 both are imaginary. Thus
`omega(Sx^a K Sy^a)` is real in either case. Its Hermitian quadratic form
is therefore symmetric, and real polarization determines the complex form.
The selected zero-field limit preserves this reality.

Consequently, including the entire positive spectrum,

```
integral E^(-1) nu_j(k,dE) = J kappa/2.                          (5)
```

This identity does not assert a frequency-independent response or a pole.
The general continuity/current sum-rule route appears also in Stringari,
[cond-mat/9311020, Eqs. (15)--(18)](https://arxiv.org/pdf/cond-mat/9311020).
We derive (1)--(5) directly; that paper's assumed magnon-plus-continuum split
is not an input. The inherited susceptibility bound is
[Koma, Eq. (3.4)](https://arxiv.org/html/2105.04970v2). It constrains a spin
inverse moment; it does not give a vanishing high-frequency current tail.

## 3. Precisely which current tail is missing

Take a real Schwartz G with smooth compact Fourier support away from zero;
all continuum Fourier integrals use d^3q/(2pi)^3. For sufficiently small p
there is no lattice aliasing. Use exactly the FA105 fields and clock

```
P_p(G)=-p sum_x G(px) Sx^2,          L_p=K/(Jp).
```

Let mu_p be the L_p-energy measure of P_p(G)Omega. Poisson summation and
(4) give the exact expression

```
T_p(R) := integral_(lambda>R) lambda mu_p(dlambda)
 = 2/(Jp^2) integral d^3q/(2pi)^3 e(pq) |Ghat(q)|^2
       integral_(E>JpR) E^(-1) nu_j(pq,dE).                     (6)
```

The full first moment is therefore

```
D_p = kappa/p^2 integral e(pq)|Ghat(q)|^2
    -> D = (kappa/2) integral |q|^2 |Ghat(q)|^2,
D_p <= D.                                                       (7)
```

Uniform integrability of these first moments is **equivalent** to
lim_(R->infinity) limsup_(p->0) T_p(R)=0. Thus its missing input is exactly
the weighted inverse-current tail on the right of (6), not an unspecified
choice of particle labels. The current-current expectation involves four
spin operators. The proved energy sum rule fixes its total (5),
not its distribution over frequency. Locality and the susceptibility
bounds used here give no vanishing estimate for (6).

Even this uniform-integrability condition would not alone prove a shell.
It would remove first-energy loss in the weak limit of the added FA105
field vectors. Their projection into the original history representation
may require control, and a broad positive limiting spectral measure
can have exactly the same moments. A sharp subbranch could also coexist
with other spectral weight; its existence does not require the entire
response to concentrate.

## 4. A nonzero soft-current fraction follows

There is a quantitative positive result about (6). Fix the same all-history
subsequence as FA105. Let ||G||_2=1 and supp Ghat lie in r/2<=|q|<=r. For
0<m0<=m, FA105 proves that the actual Stone-field energy measure contains
at least w=m0^2 r/48 in [a,b], where

```
a=m0^2 r/8,                  b=24r/m0^2.
```

Choose a nonnegative continuous compactly supported function h with h<=lambda,
h=lambda on [a,b], and support in (a/2,2b). The weak limit of the source
measures dominates the actual Stone-field measure, by the reducing
projection constructed in FA105. Thus

```
liminf_p integral_(a/2<lambda<2b) lambda mu_p(dlambda)
    >= a w = m0^4 r^2/384.                                      (8)
```

Since D_p tends to D<=kappa r^2/2, at least

```
m0^4/(192 kappa) >= m0^4/96                                   (9)
```

of the *whole* weighted inverse-current response in (5)--(7) remains in
energies Jp(a/2,2b) along this subsequence. This lower bound is measured
against the original full response, not against a renormalized narrow
filter. It is positive even though its conservative numerical value is
small. The interval has fixed nonzero width after scaling; (8) is not a
claim of concentration on a ray E=v|k|.

## 5. Fractional energy is tight; first energy is not established to be tight

For any positive spectral measure with first moment at most D and any
0<s<1,

```
integral_(lambda>R) lambda^s dmu <= D/R^(1-s).                    (10)
```

Equation (7) makes this uniform in p. The analogous FA105 Q_p bound uses
D=6kappa||F||_2^2. The measures have uniformly bounded total mass, so weak
convergence implies convergence of every such fractional moment, to the
measure in the enlarged representation. For actual bounded Weyl-history
vectors, whose spectral measures converge in the original FA100
representation and whose first moments are uniformly bounded by its
energy-kick estimate, (10) gives convergence of fractional energy directly
in that representation. No extra equality for the projected unbounded
Stone-field vectors is inferred.

There is no uniform R-decay at s=1 in (10). The exact countercontrol in
`bounds.py` has a uniform continuous part on [1,2] of mass 1/24 and an atom
at n of mass 1/(16n). Its first moment is exactly 1/8, half of which escapes
to infinity. Its total and inverse moments satisfy the conservative FA105
P bounds even for m0=1/2 and kappa=1/4 at unit momentum. Its weak limit is
atom-free. Assigning the current measure using (4) enforces the continuity
sum rule as well. This is a positive-measure control, **not** a constructed
Heisenberg ground state. It proves only that the quoted moments and Ward
identity are insufficient to supply the missing estimate.

## 6. Bare twist energy versus minimized energy

For any finite smooth Hamiltonian family with a simple ground level,
ordinary differentiation of its eigenvalue equation gives

```
E0''(0) = <H''(0)> - 2 <H'(0) Q (H-E0)^(-1) Q H'(0)>,             (11)
```

where Q projects off the ground vector. The second term is a nonnegative
inverse-current response when the parameter is a bond twist. FA105's
fixed-state twist form supplies the first term; it is not an evaluation of
the relaxed curvature (11). The exact two-site singlet control gives
<H''>=1/2 and twice the inverse-current moment=1/2: an open-bond twist is
pure gauge and its ground energy curvature is zero. This finite example
is not an assertion that the infinite ordered source has zero stiffness.
It prevents identifying the bare coefficient with a physical speed without
controlling the current response and the thermodynamic/twist limits.

## Replay and scope

```
python -W error -m pytest -q model/spin_current_spectrum
python -W error -m model.spin_current_spectrum --write
python -W error -m model.spin_current_spectrum --check
```

The packet recomputes exact finite algebra and rational bounds, and pins
this proof and the narrow inherited mathematical sources. It does not rerun
all prior producers. It leaves the actual current-tail estimate, sharp
subbranch residue/dispersion, interacting scattering limit, and empirical
particle identification unproved.
