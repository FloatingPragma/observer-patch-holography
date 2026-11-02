# Nonzero energy in the actual continuum spin fields

The regular continuum fields constructed in [FA100](../spin_canonical_field/README.md)
have genuine positive excitation spectrum, not merely a nonnegative implementing
generator. Both canonical Stone fields create nonzero positive-energy vectors.
For arbitrarily small momentum annuli their energy measures retain an explicit
positive weight in an energy annulus proportional to momentum. This proves
nonstatic, gapless dynamics inside the actual continuum field representation.
It does not prove a sharp mass shell or a nonzero connected scattering amplitude.

We also derive an exact microscopic twist-energy form with a strictly positive
coefficient. Its limit bounds continuum energy from above. Equality would require
an additional control on energy escaping the limiting representation; it is not
inferred from convergence of bounded field histories.

## 1. Unchanged source, state and clocks

The source is the same cubic spin-1/2 antiferromagnet

```
H_J = J sum_<xy> (Sx.Sy-1/4),        J>0,
```

with each unoriented nearest-neighbor bond counted once. The state omega is the
same infinitesimal-axis-1-staggered-field-selected zero-field ground. Put
`m=omega((-1)^x Sx^1)>1/1000`. The positive ground generator is K. The FA100
fields use `X=p*x`, `T=J*p*t` and the positive scaled generator `L_p=K/(Jp)`:

```
Q_p(F)=p^2 sum_x (-1)^x F(px) Sx^3,
P_p(G)=-p sum_x G(px) Sx^2.
```

Their regular limiting Weyl groups have Stone generators Q(F), P(G), with
`[Q(F),P(G)]=i m integral FG` in the Weyl sense. Their source-frame residual
valley is `(0,0,pi)` as in FA100. No spatial or Lorentz spin identification is
added here.

All Fourier integrals below use `d^3k/(2pi)^3`. Initially take real Schwartz
tests with smooth compact Fourier support avoiding zero. This class is preserved
by multiplication of the Fourier transform by `|k|` or `1/|k|`. Its functions
are limits of compact spatial tests in the FA100 continuity norms.

The inherited reflection-positive input is precisely
[Koma, arXiv:2105.04970v2, Eq. (3.4)](https://arxiv.org/html/2105.04970v2),
with the source normalization already matched in
[FA92](../spin_goldstone_sector/README.md):

```
chi_a(k) <= 1/(2J e(k+Q)),       e(k)=sum_j(1-cos k_j),
Q=(pi,pi,pi),                   a=2,3.
```

It first holds for the unique finite positive-field ground; transverse means
vanish. The finite bound is uniform in the field at fixed momentum away from
its singular point. Weak spectral convergence and the lower-semicontinuous
function `1/E`, assigned value infinity at E=0, show that the corresponding
fixed-p packet has no zero-energy component in the zero-field ground. Thus a
zero-energy atom is not silently removed by writing an inverse moment over
positive energies.

## 2. A strictly positive exact twist coefficient

The selected state is invariant under internal rotations about axis 1, cubic
coordinate permutations, and the appropriate even translations. It is also
real in the S3 product basis: each finite positive-field Hamiltonian is real,
and its unique ground can be chosen real before taking the stated limits.
Consequently the nearest-neighbor spin covariance is diagonal:

```
omega(Sx^a Sy^b)=diag(c_parallel,c_transverse,c_transverse)_(a,b).
```

The potentially antisymmetric transverse covariance vanishes because it contains
one imaginary S2 matrix in a real state. All bond directions have the same
coefficients. Define the actual, not fitted, number

```
kappa = -(c_parallel+c_transverse).
```

For either transverse axis and a real summable f, the exact double commutator
and ground identity give

```
omega(A_f K A_f)=(J*kappa/2) sum_<xy> (f_x-f_y)^2
                 =J*kappa integral e(k)|fhat(k)|^2,
A_f=sum_x f_x Sx^a.
```

For a staggered test replace e(k) by e(k+Q). This is equality before the
continuum limit. The transverse bond operator has norm 1/2, so `kappa<=1/2`.
To prove a positive lower bound, use a normalized real packet near k=0 in

```
m = |omega([A_2(f),A_3((-1)^x f)])|
  <= 2 sqrt(chi_2(f) * omega(A_3((-1)^x f) K A_3((-1)^x f))).
```

The first factor tends to at most `1/(12J)`, while the second tends to
`6J*kappa`. Hence

```
m^2/2 <= kappa <= 1/2.                                  (1)
```

In particular `kappa>1/2000000`. This is a source-specific bound on an exact
short-distance correlation coefficient; it is not a measured stiffness or
speed and does not give its numerical value.

For a real angle field theta about either transverse axis,

```
q(exp(i sum theta_x Sx^a)Omega)
    =J*kappa sum_<xy> [1-cos(theta_x-theta_y)].             (2)
```

The linear sine term vanishes by the same real-state covariance. Applying (2)
to the two scaled fields, and using Taylor's remainder plus Riemann sums,
gives the exact microscopic limits

```
q(exp(iP_p(G))Omega)/(Jp) -> (kappa/2) integral |grad G|^2,
q(exp(iQ_p(F))Omega)/(Jp) -> 6*kappa integral F^2.           (3)
```

They also add for the symmetric Weyl product `W_p(F,G)` of FA100. One direct
error estimate explains why no cross term is omitted. Let O_x be the SO(3)
rotation induced by its two onsite rotations, let C be the diagonal bond
covariance, and put `a_x=-pG(px)`, `b_x=p^2(-1)^xF(px)`. The bond energy is

```
-(J/2) Tr[C (O_x-O_y)^T (O_x-O_y)].
```

Linearizing O at the identity gives `L2*(a_x-a_y)+L3*(b_x-b_y)`.
The C-weighted cross trace is zero and each diagonal trace equals `-kappa`.
If `t=max_x(|a_x|+|b_x|)<=1`, the operator-norm linearization error on an
edge is at most `t*(|a_x-a_y|+|b_x-b_y|)`. Since `||C||_1<=3/4`, the absolute
energy error is at most

```
3 J t sum_<xy> [(a_x-a_y)^2+(b_x-b_y)^2]=O(Jp^2).
```

Thus the microscopic limit is the positive quadratic form

```
D(F,G)=(kappa/2) integral |grad G|^2+6*kappa integral F^2.  (4)
```

Bounded spectral convergence only implies
`||L^(1/2)W(F,G)Omega||^2<=D(F,G)` for the actual continuum generator L.
Equation (4) is not asserted to be equality in the continuum representation.

## 3. Bringing the Stone vectors into the same history representation

This is the step that distinguishes an actual field excitation from a separate
sequence of variational vectors. Fix any FA100 all-word subsequence. Add the
vectors `v_p=P_p(G)Omega` and `w_p=Q_p(F)Omega`, for a countable dense test set,
together with every bounded Weyl history applied to them, to the Gram matrices
used for compactness. Their norms are uniformly bounded by FA100's static
bounds. Their scaled first energy moments are bounded by the exact identities
above. The homogeneous form of the FA100 energy-kick inequality,

```
sqrt(q(U psi)) <= sqrt(q(psi))+delta*||psi||,
```

bounds the energy of every added finite history vector as well. Passing to a
further subsequence therefore gives an enlarged Hilbert space H_tilde with a
positive strongly continuous time group. The original FA100 history space H
is a closed reducing subspace for every bounded Weyl history and for this time
group. Write E_H for its orthogonal projection. This construction does not
change any of the already selected bounded-field correlation functions.

We claim `E_H v=P(G)Omega`. This needs weak differentiation, not an unproved
strong convergence of unbounded operators. For an original history vector
xi_p, the exact finite-p identity is

```
<xi_p,(exp(isP_p)-1)Omega>
    = i integral_0^s <exp(-irP_p)xi_p,v_p> dr.              (5)
```

The v_p norms are bounded. FA100's full-history amplitude continuity makes
the integrands asymptotically equicontinuous in r when tested against these
original history vectors. Approximate the integral by finite Riemann sums,
pass to the chosen Gram limit, and then refine the sums. Since H reduces the
limiting Weyl group, the result is

```
(W(s)-1)Omega = i integral_0^s W(r) E_H v dr               (6)
```

weakly against a dense subset of H, hence strongly. The integrand is strongly
continuous by FA100 regularity. Dividing by s and using Stone's theorem proves
that Omega is in the generator domain and `P(G)Omega=E_H v`. The same argument
identifies `Q(F)Omega=E_H w`. No second-moment uniform-integrability assumption
has been used. Extra norm in the enlarged vectors, if any, is not declared to
be part of the original field representation.

Because E_H commutes with time translations, the actual field spectral measures
are dominated by those of v and w in H_tilde. Lower semicontinuity transfers
positive first and inverse moments; for the inverse moment again use value
infinity at E=0. In particular both actual field vectors have no component in
the entire zero-energy subspace.

## 4. Energy bounds for the actual continuum fields

Let nu_P,G and nu_Q,F be the positive energy spectral measures of the actual
vectors P(G)Omega and Q(F)Omega for the continuum generator L. The preceding
argument, exact lattice first moments, and the susceptibility bound give

```
integral E dnu_P,G <= (kappa/2) integral |k|^2 |Ghat(k)|^2,
integral E^(-1) dnu_P,G <= (1/12) integral |Ghat(k)|^2,

integral E dnu_Q,F <= 6*kappa integral |Fhat(k)|^2,
integral E^(-1) dnu_Q,F <= integral |Fhat(k)|^2/|k|^2.      (7)
```

Both measures give zero mass to E=0. All integrals on the right use the stated
normalized Fourier measure. The static bounds remain

```
||P(G)Omega||^2 <= (1/4) integral |k| |Ghat(k)|^2,
||Q(F)Omega||^2 <= 3 integral |Fhat(k)|^2/|k|.
```

The Weyl CCR and these finite variances imply the usual uncertainty inequality
in quadratic-form sense; no product-domain assumption for QP is needed.
Choose `Fhat=|k|Ghat` in that inequality, or `Ghat=Fhat/|k|`, respectively.
This gives the strictly positive lower bounds

```
||P(G)Omega||^2 >= (m^2/12) integral |k| |Ghat(k)|^2,
||Q(F)Omega||^2 >= m^2 integral |Fhat(k)|^2/|k|.             (8)
```

Their means vanish by the inherited transverse spin symmetry, or by the same
weak differentiation applied to the zero microscopic means. Thus these are
excitation norms, not disconnected one-point contributions.

## 5. A nonzero continuum band at every small momentum scale

Take either test to have unit L2 norm and Fourier support in
`r/2<=|k|<=r`, with any fixed r>0. A positive measure with inverse moment chi
and first moment h obeys

```
nu([a,b]) >= nu(total)-a*chi-h/b,       0<a<b.
```

Using (1), (7) and (8), for any certified `0<m0<=m`, yields

```
nu_P,G([m0^2*r/8, 24*r/m0^2]) >= m0^2*r/48,
nu_Q,F([m0^2*r/16,12*r/m0^2]) >= m0^2/(2*r).               (9)
```

For example m0=1/1000 gives positive P weight at least `r/48000000`
in `[r/8000000,24000000*r]`. These are rigorous interval/weight bounds,
not fitted dispersion curves or individual energies.

Translation covariance of the Stone fields makes G -> P(G)Omega and
F -> Q(F)Omega intertwine envelope translations with the continuum spatial
unitary group. Their joint momentum support is therefore contained in the
chosen annulus (reflection of momentum does not change this radial statement).
The positive weight in (9) belongs simultaneously to that momentum annulus.
Letting r decrease to zero proves that the actual continuum field generator
has positive spectrum accumulating at zero. In particular its time action is
not identically the identity. No narrowing energy filter has been used to
manufacture concentration: the intervals have a fixed ratio of endpoints,
and their weights are measured in the original field vectors.

The same statements hold in every FA100 subsequential limit, since the
construction above can be carried out after fixing any such limit. They do
not prove uniqueness of that limit, a spectral atom, a relativistic pole,
Lorentz symmetry, or a nonzero scattering amplitude. Zero spectral threshold
is not an identification with a measured massless species.

## 6. What energy tightness is still missing

The first moments in (7) give tightness of the total spectral weight at large
energy; the inverse moments exclude loss into E=0. They do not give uniform
integrability of energy itself. The positive measures

```
nu_n=(1-1/n) delta_1+(1/n) delta_n
```

have unit mass, inverse moments at most one, and first moments tending to two.
They converge weakly to delta_1, whose first moment is one. This is a spectral
countercontrol, not an alternative Heisenberg ground state. It shows why the
exact microscopic form (4) does not yet identify the full continuum Hamiltonian.
Likewise continuous measures can obey all the displayed coarse moment bounds;
those bounds alone cannot prove a sharp shell. The stronger ground-specific
current/force estimate in [the shell diagnostic](../spin_shell_diagnostic/README.md)
remains a distinct requirement.

## 7. Verification and custody

The new code checks the exact two-site spin and rotation identities, rational
window margins, positive stiffness interval, and both spectral negative
controls. Finite checks audit the proof; they do not replace its infinite-volume
argument or claim a native continuum measurement. The receipt binds narrow,
current mathematical inputs and marks inherited producer replay as unperformed.
No earlier source, ground, clock, receipt or global claim ledger is changed.

```
python -W error -m pytest -q model/spin_canonical_energy
python -W error -m model.spin_canonical_energy --write
python -W error -m model.spin_canonical_energy --check
```
