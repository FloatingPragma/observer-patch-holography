# Controlled interacting two-body band-edge limit

Owners: [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142),
[#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143).

This closes a positive, explicitly recentered two-body limit of the same Weyl/contact
family. It does not close a many-body continuum limit, identify a physical particle,
or derive a relativistic massive/massless spectrum. All eight edge valleys are retained.
The earlier fixed-zero, `c=v/a` transparency theorem remains valid on its different domain.

## Declared family and comparison

On the infinite lattice at total momentum `K=0`, take the original pair Hamiltonian
`H_a=H_0(c_a)+U_a |d><d|`, `c_a=sqrt(3)/a^2`, `a>0`. Its free upper edge is
`B_a=2 sqrt(3)c_a=6/a^2`. The contact resolvent is the full-torus average

```
g_a(E) = < E/(E^2 - 4 c_a^2 sum_j sin(k_j)^2) >,
<f> = integral_[-pi,pi]^3 f(k) d^3k/(2pi)^3.
```

The comparison `J_a=B_a-H_a` explicitly shifts energy and reverses Hamiltonian time.
Consequently `S_J=S_H^*`; replacing this by `S_H` would give the wrong orientation.
No change to the finite source action or unannounced vacuum subtraction is made.
The parameters `a`, `c_a`, `U_a`, and the new positive inverse-length parameter `gamma`
are supplied family choices, not derived constants of nature. In particular the bare
coupling below is defined by a source integral, not by experimental fitting.

At each of the eight maxima `k_j=+/-pi/2+a p_j`, the free edge deficit tends to `|p|^2`.
Writing `Q=sum sin(a p_j)^2`, rationalization and
`|sin(x)^2-x^2|<=x^4/3` give the explicit bound

```
|B_a-2 c_a sqrt(3-Q)-|p|^2|
 <= (a^2/3) [sum_j p_j^4 + (sum_j p_j^2)^2].
```

This is a quadratic band edge, not a derivation of Lorentz symmetry.

## Full-torus threshold estimate

Set `C=sum cos(k_j)^2`, `W(eta)=<1/(C+eta)>`, `W0=W(0)`. The exact identity is

```
g_a(B_a+z) = (a^2/2)(1+a^2 z/6) W(a^2 z+a^4 z^2/12).
```

The map `q_j=cos(k_j)` includes every preimage and gives measure
`pi^-3 J(q)d^3q` on `[-1,1]^3`, where `J=product_j(1-q_j^2)^(-1/2)`.
Inside the ball `r<=R=1/2`, `0<=J-1<=r^2`: the product of `1-q_j^2`
is at least `1-r^2`, and `(1-r^2)^(-1/2)<=1+r^2` for `r^2<=1/4`.
The complementary measure is at most one and `C>=R^2` there. Jensen and the ball/tail
split give `2/3<=W0<5`; the ball contributes at most `4/(sqrt(3)pi^2)<1`
and the complement at most `4`. Thus `W0` is a finite, unambiguously source-defined number.

Let `A(r^2)` be the spherical average of `J`. Its absolutely convergent expansion
`A(t)=sum_n A_n t^n` has `A_0=1` and `0<=A_n<=1`: expand each factor with
coefficients `binom(2m,m)/4^m<=1`, then majorize the sum by the multinomial
expansion of `(u_1^2+u_2^2+u_3^2)^n=1`. This coefficient argument also bounds
the derivatives on compact subintervals and justifies termwise principal values below.

For exterior `eta>=0`, the constant-J ball contribution to `W0-W(eta)` is
`(4/pi^2)sqrt(eta) atan(R/sqrt(eta))`. Subtracting `(2/pi)sqrt(eta)` costs
at most `8 eta/pi^2`; the `J-1` correction costs `2 eta/pi^2`; the tail costs
`16 eta`. Therefore

```
|W(eta)-W0+(2/pi)sqrt(eta)| <= 18 eta.
```

For interior `eta=-t+i0`, `0<t<=1/8`, the constant-J ball's real correction is
bounded by `16t/pi^2`; its imaginary part is exactly `-(2/pi)sqrt(t)`.
The tail's real correction is bounded by `32t`. For `n>=1`, polynomial division gives

```
PV integral_0^R r^(2n)/(r^2-t) dr
 = sum_(j=0)^(n-1) t^j R^(2(n-1-j)+1)/(2(n-1-j)+1)
   - t^(n-1/2) atanh(sqrt(t)/R).
```

Using `A_n<=1`, the sum of absolute values is at most
`R/[(1-t)(1-R^2)] + t/[R(1-t/R^2)(1-t)] <= 4/3`.
The `J-1` real correction is therefore at most `16t/(3pi^2)`. Its imaginary
correction is at most `(2/pi)sqrt(t)t/(1-t)<=8t/21`. Altogether,

```
|W(-t+i0)-W0+i(2/pi)sqrt(t)|
 <= [32+16/9+16/27+8/21]t < 36t.
```

These bounds use the actual limiting boundary value, with no finite imaginary regulator.
Applying the exact Green identity, `W0<5`, and
`|(1+/-h/6)sqrt(1+/-h/12)-1|<=h/4` for `0<=h<=1/8`, yields

```
|(g_a(B_a+z)-a^2 W0/2)/a^3 + sqrt(z)/pi| <= 19a|z|,
```

for real `|a^2 z|<=1/8`, using the upper-half-plane square root
(`sqrt(-s+i0)=i sqrt(s)`). On the exterior the stronger constant 11 holds.
The branch and the limit at zero are explicit; no derivative estimate is inferred from
this value estimate.

## Nonidentity scattering with a uniform rate

Choose **exactly**

```
1/U_a = a^2 W0/2 - a^3 gamma/pi,  gamma>0.
```

For `a gamma<=1`, this is positive by `W0/2>=1/3>1/pi`.
At source energy `B_a-s`, `s>=0`, the denominator is

```
D_a(s)/a^3 = (-gamma+i sqrt(s))/pi + e_a(s), |e_a(s)|<=19as.
```

The rank-one stationary formula `S_H=conj(D_a)/D_a` consequently gives

```
S_H(B_a-s) -> (gamma+i sqrt(s))/(gamma-i sqrt(s)),
S_J(s)    -> (gamma-i sqrt(s))/(gamma+i sqrt(s)).
```

For `0<=s<=S`, `a^2 S<=1/8`, and `836aS<=7gamma`, the limiting denominator
has modulus at least `gamma/pi`, the perturbed denominator at least `gamma/(2pi)`,
and

```
sup_[0,S] |S_H-S_H,limit| <= (1672/7) a S/gamma.
```

The same estimate holds for `S_J`. It proves operator-norm convergence of the
scalar multipliers on each compact interval in the explicitly identified contact
absolutely continuous spectral representation. In particular it controls all unit
packets supported there. It does not identify full physical momentum-space scattering
operators or construct a continuum Hilbert-space embedding. On the orthogonal
complement of the cyclic contact subspace the rank-one perturbation is zero;
those free channels are retained, not removed by a one-valley approximation.
The comparison is extended continuously at threshold where `S=1`.

For `0<epsilon<=2`, the exact rational schedule

```
a <= min(1/gamma, 1/[8(1+S)], 7 epsilon gamma/[1672 S])
```

ensures the uniform error is at most `epsilon`. No gamma-zero threshold theorem
is asserted. A generic choice `U_a=kappa/a^2` without the displayed critical
countertuning instead has `D_a/a^2 -> 1/kappa-W0/2 !=0`, hence `S->1`.

## Bound pole retained

For `E>B_a`, `g_a(E)` is strictly decreasing from `a^2 W0/2` to zero:
it is the Stieltjes transform of a positive unit-mass spectral measure.
The chosen positive `1/U_a` lies strictly between these endpoints. Thus there
is exactly one upper bound pole `E_a=B_a+delta_a`. This existence proof does
not use the stronger sufficient condition `U>B_a` in the earlier finite fiber helper.
For `a gamma<=1/256`, the exterior estimate gives opposite signs at
`delta=gamma^2/4` and `4gamma^2`, so

```
gamma^2/4 <= delta_a <= 4gamma^2,
|delta_a-gamma^2| <= (5016/7) a gamma^3.
```

Indeed `|sqrt(delta_a)-gamma|<=19pi a delta_a`, and multiplication by
`sqrt(delta_a)+gamma<=3gamma` gives the displayed rate. The corresponding
`J_a` eigenvalue tends to `-gamma^2`. The original energy diverges as
`6/a^2+delta_a`; it is not a finite physical rest mass. We do not claim convergence
of eigenvectors, a pole residue, a norm-resolvent limit, or stability of a many-body vacuum.

## Verification and imported mathematics

`limit.py` evaluates exact rational bounds and separates numerical diagnostics.
The diagnostics integrate the original torus after analytic third-axis integration,
using the arcsine convolution density `2 K(s(2-s))/pi^2`. Tests check the Green
identity against an independent three-dimensional torus sum, the on-shell jump,
uniform estimates, all eight dispersions, pole brackets, rounding boundaries,
missing countertuning, and reversed scattering orientation.

Imported ingredients are elementary coarea, geometric-series and principal-value
calculus, the spectral theorem, and the already declared rank-one stationary
scattering formula. The model-specific claims are the full-valley coefficient,
explicit countertuning, and quantitative error bounds above. Lattice contact
renormalization and repulsively bound pairs are established mechanisms, not new
physical discoveries: compare [Koerber et al., arXiv:1912.04425](https://arxiv.org/abs/1912.04425)
and [Winkler et al., cond-mat/0605196](https://arxiv.org/abs/cond-mat/0605196).

Remaining: justify a physical vacuum and time/energy convention; control nonzero
total momentum, all particle sectors and interaction stability; derive appropriate
massive/massless relativistic sectors and units; supply the physical detector map.
