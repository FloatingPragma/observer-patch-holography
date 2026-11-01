# A controlled moving bound-pair continuum sector

This extends the [zero-total-momentum band-edge result](../band_edge_contact/README.md)
to `K=aP` with all three components of P variable. It closes a moving-pole and
wavepacket-propagation premise on the original infinite-lattice two-body action.
Owners: [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142)
and [#1143](https://github.com/FloatingPragma/observer-patch-holography/issues/1143).

The result is one composite-particle sector, with quadratic limiting dispersion.
It is not a many-body QFT, a relativistic mass shell, a physical particle identity,
or an on-shell scattering theorem at nonzero total momentum. All spin components,
edge valleys and orthogonal channels remain in the source Hilbert space.

## One global Hamiltonian comparison

Retain `H_a=dGamma(h_a)+U_a sum_x n_x,up n_x,down`, with
`h_a(k)=-c_a sigma.sin(k)`, `c_a=sqrt(3)/a^2`, and the same coupling as before:

```
1/U_a = a^2 W0/2 - a^3 gamma/pi,
W0 = <1/(sum_j cos(k_j)^2)>,     gamma>0.
```

The full-torus average W0 is a fixed source integral, with `2/3<=W0<5`.
The comparison is the **single global** operator `J_a=6/a^2-H_a` on the two-body
space. It shifts energy and reverses model time. In particular we do not subtract
a momentum-dependent band edge from the Hamiltonian: that would remove the
center-of-mass kinetic energy that this construction is meant to derive.
The scaling, gamma, source state and time/energy convention are supplied model
choices, not unexplained physical identifications or measured inputs.

## Exact free edge and four-spin Green function

For `|K_j|<=pi/2`, put `q=k-K/2`, `C_j=cos(K_j/2)`, `S_j=sin(K_j/2)`,
`d_j=C_j sin(q_j)`, `e_j=S_j cos(q_j)`, and `R^2=sum C_j^2`.
The two one-particle spin vectors are `d+e` and `-d+e`. Concavity of square root gives

```
|d+e|+|-d+e| <= 2 sqrt(d^2+e^2) <= 2R.
```

Every `q_j=+/-pi/2` attains equality. Hence the exact upper edge is
`B_a(K)=2c_a R`; all eight extrema are retained. Pauli algebra, or the Schur
complement of the singlet/triplet blocks, gives the exact singlet resolvent

```
G_K(q,E) = E(E^2-4c_a^2 e^2) /
           [E^4-4c_a^2(d^2+e^2)E^2+16c_a^4(d.e)^2].
```

For `E=B_a(K)+z`, `z>=0`, let `t=E/(2c_a)`, `eta=t^2-R^2`, and
`V(q)=sum_j cos(K_j) cos(q_j)^2`. Then the full contact Green function is

```
g_K(E) = t/(2c_a) < (1-e^2/t^2)/(V+eta+(d.e)^2/t^2) >.
```

At the edge the formula is interpreted by its integrable limiting value.
It does not use the restrictive sufficient condition `U>B_a(K)` of an older
finite-fiber numerical helper. The physical fermion exchange subspace is
preserved by the free resolvent and contains the constant spin singlet d.

## Uniform momentum perturbation estimate

Set `kappa=max_j|K_j|<=1/2`, `Q=sum cos(q_j)^2`, `h=a^2 z<=1/8`,
and `alpha=R/sqrt(3)`. Then

```
1-kappa^2/4 <= alpha <= 1,     R^2 >=45/16,
e^2 <= kappa^2 Q/4,           (d.e)^2/t^2 <= e^2,
eta = alpha h+h^2/12,        V = Q-2e^2.
```

Let `D0=Q+h+h^2/12` and `DK=V+eta+(d.e)^2/t^2`.
These exact inequalities imply

```
(1-kappa^2/2)D0 <= DK <= D0,
0 <= D0-DK <= (kappa^2/2)D0.
```

Writing `IK=(1-e^2/t^2)/DK` and `I0=1/D0`, one obtains

```
|IK-I0| <= (4 kappa^2)/(7D0) + 32 kappa^2/315.
```

The first term controls the denominator difference; the second uses
`e^2/(t^2 DK)<=kappa^2/(4R^2(1-kappa^2/2))`. Since `<I0><=W0<5`
and `alpha+h/6<=49/48`, the prefactor difference gives

```
|g_K(B_a(K)+z)-g_0(6/a^2+z)|
 <= a^2 kappa^2 [5/8+(49/96)(20/7+32/315)]
 < 3 a^2 kappa^2.
```

The previous exterior threshold bound has constant11. Therefore for `K=aP`,
`max|P_j|<=Pmax`, and the displayed domains,

```
| (g_K(B_a(K)+z)-a^2 W0/2)/a^3 + sqrt(z)/pi |
 <= a(11z+3Pmax^2).                                      (1)
```

This is a real exterior estimate, including the finite edge limit. It is not
used as though it were a bound for interior principal values or derivatives.

## Moving bound pole and uniform dispersion

Assume the exact rational condition

```
256 a (gamma^2+Pmax^2) <= gamma.                          (2)
```

It implies `a gamma<=1/256`, `a Pmax<=1/512`, positivity of U_a, and all
domains in (1) for `gamma^2/4<=z<=4gamma^2`. The pole denominator divided
by a^3 differs from `(sqrt(z)-gamma)/pi` by at most `a(11z+3Pmax^2)`.
It is negative at `z=gamma^2/4` and positive at `z=4gamma^2`:
the endpoint errors are bounded by `3gamma/256` and `44gamma/256`,
respectively, while the limiting magnitudes exceed `7gamma/44` and
`7gamma/22`. The contact Stieltjes transform is strictly decreasing above
the complete free spectrum. Thus there is exactly one upper bound pole
`E_a(K)=B_a(K)+delta_a(K)` for every P in the cube, and

```
gamma^2/4 <= delta_a(K) <= 4gamma^2,
|delta_a(K)-gamma^2| <= (22/7) a gamma (132gamma^2+9Pmax^2). (3)
```

The latter follows by multiplying
`|sqrt(delta)-gamma|<=pi a(11delta+3Pmax^2)` by
`sqrt(delta)+gamma<=3gamma`. This is a pole-offset theorem; it does not
claim a residue or eigenvector limit without a comparison map.

For the same global J_a, the bound-band energy is
`epsilon_a(P)=6/a^2-B_a(aP)-delta_a(aP)`. Rationalizing the edge and using
`|sin(x)^2-x^2|<=x^4/3` gives

```
|6/a^2-B_a(aP)-|P|^2/4|
 <= (a^2/48)[sum_j P_j^4+(sum_j P_j^2)^2]
 <= a^2 Pmax^4/4.
```

Consequently

```
sup_P |epsilon_a(P)-(|P|^2/4-gamma^2)| <= E_a,
E_a = a^2 Pmax^4/4+(22/7)a gamma(132gamma^2+9Pmax^2).       (4)
```

The source-derived limiting inertial coefficient is `1/4` in these model
coordinates. A single free upper edge has coefficient `1/2`; thus the
composite has twice that model inertial mass. No SI mass, relativistic rest
energy, physical vacuum or physical time calibration follows. The original
source pole still diverges; the reflection convention is essential.

## Explicit wavepacket map and controlled propagation

Let `R_K=(E_a(K)-H_0,K)^(-1)` and
`psi_a,K=R_K d/<d,R_K^2 d>^(1/2)`. The pole equation shows directly that
this normalized, exchange-antisymmetric vector is the bound eigenstate.
Its overlap with d is positive, fixing its phase. Above the edge the
formula is continuous in K, hence supplies a measurable fiber vector.

For continuum packets `f` supported in `[-Pmax,Pmax]^3`, define

```
(W_a f)(K,k) = a^(-3/2) f(K/a) psi_a,K(k),
```

and zero outside the corresponding K cube. Use `dP/(2pi)^3`, `dK/(2pi)^3`
and the normalized relative-torus measure consistently. Changing variables
`K=aP` proves `W_a*W_a=I`. This is an explicit bound-band isometry into the
full source space, not deletion of its other states. Its range reduces H_a
on this two-body direct integral. The band need not be isolated from the
full free spectrum at *different* total momenta.

Let `h_pair` multiply f by `|P|^2/4-gamma^2`. Equation (4) gives the
operator estimate, including arbitrary untouched reference factors,

```
|| exp(-it J_a) W_a - W_a exp(-it h_pair) ||
 <= min(2, |t| E_a).                                     (5)
```

The code gives an exact rational spacing schedule for any finite Pmax,
duration cap and positive tolerance. For `a<=1`, replace E_a by a times
`Pmax^4/4+(22/7)gamma(132gamma^2+9Pmax^2)` and combine with (2).
No arbitrary finite-volume K is claimed to be a reciprocal-grid label.
There is no many-pair interaction or many-body limit in (5).

## Exterior transition response and numerical controls

For real `0<=z<=Z`, `|sqrt(z)-gamma|>=d>0`, and
`44a(11Z+3Pmax^2)<=7d`, (1) also gives

```
|a^3 tau_a(B_a(K)+z)-pi/(sqrt(z)-gamma)|
 <= (968/49) a(11Z+3Pmax^2)/d^2.
```

This is the real, exterior renormalized rank-one transition coefficient,
away from its pole. It is not a quantum channel or an on-shell S matrix.

Numerical diagnostics integrate the *difference* from the independently
known K=0 Green function over all six pyramids of the relative-momentum
cube. The radial Jacobian removes its integrable origin singularity;
an algebraically subtracted numerator avoids large near-edge cancellation.
Two quadrature orders, direct four-spin resolvents, all eight edge extrema,
the forbidden momentum-dependent subtraction and the opposite time
orientation are independent controls. Floating agreement is diagnostic;
the analytic inequalities prove (1)--(5).

Remaining: a physical detector and vacuum map, nonzero-K on-shell scattering,
many-body stability and interactions, relativistic massive/massless sectors,
and the joint local stress/geometry/Einstein construction.
