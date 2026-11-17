# A matched regional Einstein limit

Use signature `(-,+,+,+)` and units `hbar=c=1`. Every tensor component below
is measured in one local orthonormal tetrad. The theorem is conditional on
regional source data; no such data are constructed by the finite examples.

## The common source interface

For each refinement n, point p and admitted state, the same source supplies
a metric, a symmetric matter stress readout T, and regional entropy/area
readouts on the ten observers specified below. Their source-state and
refinement indices agree. The metric is attached to the source action and
causal/proper-volume readout. The geometric tensor G is its actual Einstein
tensor, or an approximation with a controlled component error. The stress
is the same source that enters the matter action and propagation; replacing
it by independently chosen modular charges is outside the interface.

The required data are the following.

1. A common positive proper radius `ell_n` tending to zero. On each observer's
   local spatial geodesic slice, intrinsic geodesic balls and their comparison
   balls have equal proper volume. Each regional readout is available on
   that slice, using the metric's proper-clock normalization.
2. A bulk entropy contrast `s_(n,r)`, an area contrast `a_(n,r)`, a positive
   area coefficient `eta_n`, and one isotropic reference scalar `Lambda_n`.
   The coefficient and reference scalar are common to all ten observers at
   the point. Their values are independent of the choice of probe.
3. The three uniform bounds in equation (1), on those same readouts, with
   nonnegative errors `e_bulk,n`, `e_area,n`, `e_stationarity,n`.
4. Each error divided by `ell_n^4` tends to zero on this same refinement.
   The tensor readouts and `Lambda_n` converge. The area coefficient converges
   to one positive **spacetime-independent** constant `eta_*`.
5. For the smooth-field conclusion, the metric converges sufficiently in
   local charts to identify the tensor limit with `G[g]`, for example local
   C2 metric convergence with nondegenerate Lorentzian limit. Stress limits,
   tetrads and all source readouts use these same charts and state family.
   Independent Ward conservation, contracted Bianchi conservation and
   connectedness supply the constant-Lambda conclusion.

These are mathematical and physical hypotheses, not labels accepted by the
packet. Reading the simulator specification can identify the laws and
readout maps; it does not establish their regional estimates or convergence.
The finite reciprocal source has no common spatial refinement or regional
edge-area law meeting this interface. Its thermal coefficient depends on
two selected geometry values, temperature and cutoff, so its positivity
alone does not establish item 4.

## Geometric normalization and the three errors

At a smooth point p and future unit timelike vector u, take spacetime normal
coordinates with time axis u. The hypersurface `exp_p(u^perp)` has zero
second fundamental form at p: its coordinate embedding has zero second
derivatives there and the ambient Christoffel symbols vanish there.
This is a local construction for each u; it requires no globally
time-symmetric foliation. The source must admit the corresponding regional
readouts. The Gauss constraint gives `R^(3)(p)=2 G_ab u^a u^b`.

For an intrinsic three-dimensional geodesic ball of radius r, the normal
coordinate volume expansion gives

```
V(r)=4 pi r^3/3 * [1-R^(3) r^2/30+o(r^2)],
A(r)=4 pi r^2 * [1-R^(3) r^2/18+o(r^2)].
```

At fixed volume `4 pi ell^3/3`, solving for r gives
`r=ell+R^(3) ell^3/90+o(ell^3)`. Substitution gives
`A=4 pi ell^2-2 pi R^(3) ell^4/15+o(ell^4)`.
Subtract a common local vacuum reference with Einstein tensor
`G_ref=-Lambda g`, so its timelike contraction is Lambda. Thus the leading
fixed-volume area contrast is `-C_ell (G_uu-Lambda)`, where
`C_ell=4 pi ell^4/15`. On another prescribed slice the correction is
`R^(3)/2-G_uu=(K_ij K^ij-K^2)/2`; it must be included in the area error.
A fixed nonzero correction is order `ell^4`, so it does not disappear after
division by `ell^4`.

The radial diamond-kernel integral is
`integral_0^ell 4 pi r^2 (ell^2-r^2)/(2 ell) dr=C_ell`.
The geometrically normalized modular factor `2 pi` therefore gives the
bulk coefficient `2 pi C_ell=8 pi^2 ell^4/15`. These are the exact
normalizations in `EinsteinBranch/SmallBall.lean`. Identifying the source's
regional modular generator and its actual stress with this kernel is a
separate hypothesis. In particular it is not obtained from a temperature
or clock fitted after tensor reconstruction.

For every probe r the matched estimates are

```
|s_r-2 pi C_ell T_uu| <= e_bulk,
|a_r+C_ell (G_uu-Lambda)| <= e_area,
|s_r+eta a_r| <= e_stationarity.                         (1)
```

The last line is generalized-entropy balance on the specified fixed-volume
family. Regional edge/center terms and ultraviolet subtractions must either
produce `eta a_r` with the stated error or be included in that error. Global
thermal entropy production or stationarity of a Markov distribution does
not supply it. Retained records and their entropy cannot be omitted from
the source's specified entropy functional.

## Ten observers and an exact inverse

Let e0 be the unit time direction and e1,e2,e3 the orthonormal spatial
directions. Use

```
u0=e0,
ui+=(5/3)e0+(4/3)ei,    ui-=(5/3)e0-(4/3)ei,    i=1,2,3,
uij=3e0+2ei+2ej,                                      1<=i<j<=3.
```

Their metric squares are respectively `-1`, `-25/9+16/9=-1` and
`-9+4+4=-1`. Their time components are positive. For any symmetric tensor R
write its ten measurements as `r0`, `ri+`, `ri-`, `rij`. Expanding R(u,u)
and solving gives

```
R00=r0,
R0i=9(ri+-ri-)/80,
Rii=9(ri++ri-)/32-25 r0/16,
Rij=rij/8+7 r0/16-99(ri++rj+)/320+9(ri-+rj-)/320.        (2)
```

The absolute coefficient sums are `1`, `9/40`, `17/8`, `99/80` for the four
types of component. Consequently

```
max_(mu,nu)|R_mu,nu| <= (17/8) max_r |R(u_r,u_r)|.       (3)
```

The constant is sharp: choose charges `r0=-1`, `r1+=r1-=1`, with all
other charges zero. The invertible decoder gives a symmetric tensor with
`R11=17/8`. Removing any probe also destroys determination: the decoder
of the unit charge at the missing row is a nonzero tensor invisible to the
other nine. The packet generates all ten such witnesses and independently
inverts the rational design by Gaussian elimination. The Lean proof checks
both inverse compositions, all timelike norms, (3) and the sharpness witness.
This is a coordinate maximum-norm bound in the stated orthonormal frame;
it is not a Lorentz-invariant positive tensor norm.

## Quantitative Einstein residual

Define `kappa=2 pi/eta` and the single residual tensor
`R=G+Lambda g-kappa T`. For every unit timelike probe,
`R_uu=G_uu-Lambda-kappa T_uu`. The same s and a in all three lines of (1)
give the exact identity

```
eta C_ell R_uu
  = (s-2 pi C_ell T_uu)
      + eta [a+C_ell (G_uu-Lambda)] - (s+eta a).          (4)
```

The triangle inequality, `eta>0`, `ell>0` and (3) prove

```
max_(mu,nu)|G_mu,nu+Lambda g_mu,nu-(2 pi/eta)T_mu,nu|
 <= 255/(32 pi eta ell^4)
       * (e_bulk+eta e_area+e_stationarity).            (5)
```

Thus the tensor balance is derived from regional estimates, without a
separate assumption of null-tensor balance or a fitted tensor equation.
All quantities in (4) must be those supplied by the common source interface.
Three unrelated scalar witnesses cannot be substituted into this identity.

For quantitative convergence to limiting components in a common tetrad,
suppose `eta_n>=eta0>0`, `eta_*>=eta0` and `max|T_*|<=M_T`. Besides the
right side B_n of (5), the error in the limiting Einstein equation is at most

```
B_n + ||G_n-G_*||max + |Lambda_n-Lambda_*|
    + (2 pi/eta0)||T_n-T_*||max
    + (2 pi M_T/eta0^2)|eta_n-eta_*|.                   (6)
```

Here the tetrad metric is exactly diag(-1,1,1,1); any chart or tetrad
identification error is included in the component convergence estimates.
The normalized error identity

```
(e_bulk+eta e_area+e_stationarity)/(eta C_ell)
 =15/(4 pi) * [(e_bulk/ell^4)/eta
                + e_area/ell^4+(e_stationarity/ell^4)/eta]
```

shows that the common assumptions imply B_n tends to zero. Taking the limit
in (6) gives `G[g]+Lambda_* g=(2 pi/eta_*)T`. With
`G_N=1/(4 eta_*)`, the coupling is exactly `8 pi G_N`.
Pointwise assumptions give a pointwise equation. Uniform versions on compact
subsets give local uniform residual convergence. A finite list of profile
samples does not establish either assertion.

## Constant cosmological term and nonlinear scope

On a connected smooth component let the limiting fields satisfy the
contracted Bianchi identity, metric compatibility and the independent
matter Ward identity. Since kappa is spacetime-independent, covariant
divergence of the derived tensor equation gives

```
0 + partial_b Lambda = kappa*0,
```

so Lambda is constant on that component. A vacuum reference fixes its value;
the entropy calculation does not select that value. The Lean theorem also
proves the corresponding exact finite-difference statement using the
Ward/Bianchi interface of `Tensor.lean`. Its coordinate differences do not
claim to construct a continuum covariant derivative.

The constancy of eta is necessary for this inference. With variable
`kappa=2 pi/eta`, conservation instead gives
`partial_b Lambda=T_ab partial^a kappa`. For example, on a flat component
`G=0`, `T=c g` with constant c is conserved. Any positive variable kappa
and `Lambda=c kappa` satisfy the tensor equation, while Lambda varies.

The nonlinear conclusion uses (1) for the actual regional entropy/area
contrasts of each admitted state, with actual curvature and stress values.
This can be stated as a finite-deficit hypothesis, or supplied by an
additional integration theorem for a connected state family. Ordinary
first-order entropy stationarity at a single reference state does not give
that hypothesis. If s,a,G,T in (1) are only first variations at one vacuum,
equations (4)-(6) prove only the corresponding linearized Einstein relation.
An equality of first derivatives cannot be integrated without control of
the state path, changing metric, clocks, source identification and remainder.

## What the replay establishes

The finite replay proves exact frame inversion, its sharp bound, failure
under every probe deletion, the coefficient arithmetic and source custody.
Its rational profiles `ell_n=1/(n+1)`, `eta_n=2`, `e_n=ell_n^5` are supplied
algebraic controls. The Lean module proves universal real-variable residual
and scaling theorems, conditional on their displayed hypotheses. Neither
the profile samples nor the conditional proof establish that an executable
OPH source supplies the regional entropy, fixed-volume stationarity,
universal coefficient, local stress conservation or curvature refinement.
