# Lattice realization of the Weyl stress-response family

This package constructs finite-support lattice state perturbations and
quadratic stress detectors whose nine-by-nine response matrix converges,
with explicit error bounds, to the nondegenerate Weyl family in FA91.
It supplies a **free-sector lattice-to-continuum matter-response input** to
G5. The source specialization is exactly `U=0`, `c=v/a`, in its infinite
filled-band ground state. The proof does not infer convergence of a
derivative composite observable from convergence of field distributions.

The eight Weyl valleys are retained. We excite the `(1,1,1)` valley, whose
source spin frame has `h=+sigma.p`; the other seven limiting fields remain
in their vacuum. All finite-spacing operators act on the full lattice CAR
algebra. No momentum truncation, deletion of doublers, or finite-particle
cutoff is applied. The selected quadratic detector need not represent the
sum of stress operators on arbitrary multivalley states: its asserted
response is on this specified shell family, where the limiting variation
of the seven other sectors vanishes.

## Exact source and existing geometry target

The source, centered ground, eight-valley limit and spin-frame map are
those in `continuum_field_net`, `continuum_free_limit`, and FA87
`weyl_geometric_modular`. With `U=0` the centered and original Hamiltonians
coincide. The one-particle source symbol is
`-(v/a) sum sigma_j sin(k_j)`. Put `x0=v*t`, so the generator used below is
divided by `v`. Modulation by `(-1)^(n1+n2+n3)` makes its symbol
`h_a(p)=sum sigma_j sin(a p_j)/a`. The filled-band projector is `Pi_a`.
Its undefined values at the finitely many zeroes are irrelevant to the
L2 operator and to the unique free filled-sea state.

`EinsteinBranch/Tensor.lean` already proves conditional tensor completion;
we do not replace it. `SmallBall.lean` supplies the target `R^4`
normalization and requires common limit control. `RegisterSurface.lean`
also requires geometric null balance, area, stationarity and other supplies.
Here the new input is a concrete source-bound realization of the same
matter variation space, with a joint lattice/radius limit. It is not the
whole register, nor the full G5 geometric/collar hypothesis.

The continuum stress and rest-ball modular normalization are inherited
from FA87. Its primary modular import is La Piana and Morsella,
[The fermionic massless modular Hamiltonian, Theorems 3.4 and 4.1](https://arxiv.org/html/2405.07888v2),
Communications in Mathematical Physics 406, 81 (2025),
[DOI](https://doi.org/10.1007/s00220-025-05253-1).
The signed CAR Cauchy generator and opposite particle/hole convention are
exactly those recorded there by FA87. No new modular theorem is assumed.
The estimates below are proved directly.

## A gapless projection estimate without aliasing assumptions

Use centered cells and the contraction

```
(J_a f)(n)=a^(-3/2) integral_(a(n+[-1/2,1/2)^3)) f(x) dx.
```

Then `J_a^*J_a=P_a` is cellwise conditional expectation and
`||(P_a-I)f||<=a sum_j ||partial_j f||`. Define
`p_a,j f=[f(.+a e_j)-f(.-a e_j)]/(2 i a)` in the continuum as well as
the corresponding lattice operator. Cell averaging intertwines these
translations exactly, including arbitrary high momentum components.
For `h=sigma.p`,

```
R_a=h_a J_a-J_a h,
||R_a f|| <= (a^2/6) sum_j ||partial_j^3 f||.             (1)
```

Here and below the constant source modulation is understood. Let `Pi`
be the negative-band projector of `h`. For every `f in H^3 intersect L1`
and `epsilon>0`,

```
||(Pi_a J_a-J_a Pi)f||
 <= (a^2/epsilon) || |p|^3 f||
       + epsilon^(3/2) ||f||_1.                         (2)
```

To prove (2), first restrict the positive continuum band to `h>=epsilon`.
Set `X=Pi_a J_a`, `A=-h_a|negative>=0`, `B=h|positive>=epsilon`.
The exact Sylvester equation is `A X+X B=-Pi_a R_a`. On H3 data its
semigroup solution is
`Xf=-integral_0^infinity exp(-tA) Pi_a R_a exp(-tB)f dt`.
Differentiate `exp(-tA) X exp(-tB)f` to verify the formula; the boundary
at infinity vanishes because `X` is bounded and `B>=epsilon`.
The derivative seminorms commute with continuum spectral projections and
with `exp(-tB)`, giving `a^2 sum ||partial_j^3 f||/(6 epsilon)`.
The negative-continuum/positive-lattice component has the same bound.
Their sum is at most `a^2 |||p|^3 f||/epsilon`, since
`sum ||partial_j^3 f||<=sqrt(3)|||p|^3 f||` and `sqrt(3)/3<1`.
On `|p|<epsilon`, use `||Pi_a J_a-J_a Pi||<=2` and the unitary Fourier
convention to obtain
`||1_(|p|<epsilon)f||<=epsilon^(3/2)||f||_1/(sqrt(6) pi)`.
This proves (2), including its deliberately enlarged constants.

The same argument applies after each exact source spin rotation and
valley modulation. Equation (1), rather than a Brillouin-zone ansatz for
band-limited data, controls all aliases of the compact shell.

## Explicit shell Sobolev constants

Use exactly FA91's shell `g_n`, normalized as `f_n=g_n/||g_n||`, with
`L<=|x|<=2L`, `phi(s)=s^2(s-1)^4(2-s)^4`, and its nine listed directions.
The shell is C3 across both boundaries, with its first three derivatives
zero there. Its weak fourth derivatives are piecewise bounded and have
compact support, so it is H4. No H5 statement is required.

For a unit helicity vector the angular decomposition of
`sigma.omega [a+sum omega_j b_j]` has harmonic degrees 0, 1 and 2.
After removing their common `pi^3`, the respective squared angular norms
are `396900`, `176400`, and `99225/2`. The degree-zero component is
`sum sigma_j b_j/3`, the degree-one component is `sigma.omega a`, and the
remaining even component has degree two. Orthogonality and the shell
normalization give `||g||^2/L^3=(150088295/4249388)pi^3`.

For each degree define `L_l phi=phi''+2 phi'/s-l(l+1)phi/s^2`.
The radial integrands for moments m=1,2,3,4 are respectively

```
s^2 phi'^2+l(l+1)phi^2,
s^2 (L_l phi)^2,
s^2 (L_l phi)'^2+l(l+1)(L_l phi)^2,
s^2 (L_l^2 phi)^2.                                     (3)
```

Integration by parts has no boundary distribution, by the stated
vanishing jets. The code expands the weighted Laurent polynomials and
integrates every term exactly on [1,2]. The possible logarithmic
coefficients cancel; this cancellation is checked, not assumed. The four
exact normalized moments have common denominator `4288237` and numerators
`90678696`, `5320580800`, `529844446656`, `108752906109312`.
Their integer square-root ceilings are `(K1,K2,K3,K4)=(5,36,352,5036)` and certify
`|||p|^m f_n||<=Km/L^m`, uniformly in all nine directions. These large
constants are sufficient bounds; they are not numerical fitted rates.
Rotating the helicity direction rotates the shell and its spin frame,
so these isotropic Sobolev norms do not depend on the direction.

The shell volume is `28 pi L^3/3<49 L^3`. Thus
`||f||_1<=7 L^(3/2)` and `||partial_j f||_1<=7 K1 L^(1/2)`.
Put `q=a/L<=1` and `epsilon=a/L^2`. With `F=Pi f`, `F_a=Pi_a J_a f`,
equation (2) and the exact difference intertwiner give

```
||F_a-J_a F|| <= q K0,               K0=K3+7,
||p_a,j F_a-J_a p_j F|| <= (q/L) KD,
KD=K4+7 K1+K3/6.                                        (4)
```

For the second line apply (2) to `p_a,j f`, whose H3 seminorm is bounded
by `K4/L^4` and L1 norm by that of `partial_j f`, and add the central
difference error `a^2 K3/(6L^3)`. Also
`s_a=||J_a f||^2` satisfies `0<=1-s_a<=9 q^2 K1^2`.
For `q<=1/(6K1)` this implies `s_a>=3/4`.

## Actual lattice state perturbations and detectors

Normalize the compact lattice spinor as `f_a=J_a f/sqrt(s_a)`, with the
source modulation restored. Put
`B_a=c^*(f_a)c(f_a)-omega_a(c^*(f_a)c(f_a)) I`.
It is an even selfadjoint contraction of zero ground expectation. For
each of the nine shells take the same construction and set
`rho_a,theta(X)=omega_a((I+sum theta_j B_a,j)X)`.
On the algebra of sites disjoint from all shell cells these are positive,
normalized, normal states for `||theta||_1<1`, affine in theta. The global
positive vector extensions are generally nonlinear, as in FA91.

For a future unit vector `u=(5/4,3n/4)`, define in the selected frame

```
A_u=u0 I-u.sigma,        D_u=u0 h-u.p,
Q_u(w)=1/2 {w A_u,D_u},
Q_a,u(w)=1/2 {w_a A_u,D_a,u},       w_a(n)=w(a n).          (5)
```

Both factors in each anticommutator are selfadjoint. The lattice
`Q_a,u` is finite rank for compactly supported `w`, so its vacuum-normal-
ordered second quantization is an actual bounded finite-support CAR
observable. Restoring modulation makes it an observable of the original
source. Formula (5) is the symmetric stress contraction `T_mu nu u^mu u^nu`,
not an unsymmetrized product. At `u=(1,0)`, `2pi Q_u(w_R)` is FA87's
rest-ball modular form, where `w_R=(R^2-|x|^2)_+/(2R)`. The other nine
vectors give **lab-slice stress probes**; these operators are not asserted
to be boosted-ball modular generators or logarithms of lattice densities.

Take `rho=R/L<=1/8` and `q<=1/16`. The detector and its nearest-neighbor
sites are disjoint from the shell cells: it suffices that
`R+(1+sqrt(3))a<L`, which these constraints imply. Exact free Wick
contraction and the equal-time disjoint-support CAR therefore give

```
S_a,u,n(w)=d/dtheta_n rho_a,theta(:dGamma(Q_a,u(w)):)|0
         =-<F_a,Q_a,u(w)F_a>/s_a,
S_u,n(w)=-<F,Q_u(w)F>.                                  (6)
```

No evolved finite-speed support statement for the lattice is used. The
states and detectors in (6) are at one time. The CAR spectral measurement
of a bounded quadratic observable is specified; a native gate program
implementing arbitrary shell coefficients is not supplied by this package.

## Quantitative stress convergence and the joint limit

Let `W=||w||_infinity` and `Lip(w)` its global Lipschitz constant. The ball
weight is admissible despite its boundary corner: `W=R/2`, `Lip(w)=1`.
Using (4), `||A_u||<=2`, `||D_a,u F_a||,||D_u F||<=4K1/L`, and
`||D_a,u F_a-J_a D_u F||<=6q KD/L`, one obtains

```
|S_a,u,n(w)-S_u,n(w)|
 <= (W/L) q A +16 q Lip(w) K1,
A=16 K1 K0+24 KD+48 K2+144 K1^3.                         (7)
```

Here are the separate terms, so that the derivative and normalization
requirements can be audited. Replacing `F_a,D_a F_a` by `JF,JDF` costs
`(W/L)q(8K1 K0+12KD)`. Replacing the sampled weight by cell averaging
costs at most `8a Lip(w)K1/L`, since the cell radius is less than `a`.
Replacing `J^*J` by the identity on `DF` costs at most
`24WaK2/L^2`, by the H2 estimate. Finally divide by `s_a`, use
`1/s_a<=2`, and bound the extra normalization term by
`144(W/L)q K1^3`. This proves (7) with no pointwise projector estimate.

FA91 gives `delta_n T(0)=gamma_L ell_n ell_n^T`, with
`gamma_L=c_gamma/(pi^3 L^4)` and `c_gamma=4249388/150088295`.
The same spatial Taylor bound used there applies to these unboosted
contractions, with the conservative Hessian majorant
`M*=243(3H1 H2+H0 H3)=195689447424`, where
`(H0,H1,H2,H3)=(16,1024,131072,25165824)`.
Indeed the lab derivative directions have smaller coefficient sums than
the boosted ones; alternatively use `||u||_1<3` in the same product rule.
Odd spatial Taylor terms integrate to zero. Thus the matrix

```
M_a,ij = 2pi S_a,u_i,n_j(w_R) / [(8pi^2 R^4/15) gamma_L]
```

satisfies, entry by entry,

```
|M_a,ij-D_ij| <= C_lattice q/rho^4+C_ball rho^2,
D_ij=(5/4+(3/4)n_i.n_j)^2,
C_lattice=(60/c_gamma)(A/2+16K1),
C_ball=(3/14) M* 64/c_gamma.                             (8)
```

The lattice constant follows from (7), `rho<=1` and `pi<4`.
The continuum constant is the unchanged FA91 bound. Its exact matrix
has `||D^-1||_infinity=68/5`. Choose any

```
0<rho<=2^-28,
0<q<=min(1/16,1/(6K1),rho^7/(1+C_lattice)).               (9)
```

Then the entry error is at most `(C_ball+1)rho^2` and
`9(68/5)(C_ball+1)rho^2<1/4`, checked with rational arithmetic.
Every such finite-spacing response matrix is invertible and has inverse
norm at most `272/15`. For example `rho_m=2^-28/(m+1)` and any positive
`q_m` obeying (9) gives a common explicit refinement for all nine states
and detectors. This sufficient spacing is extremely small; no assertion
of optimality or practical laboratory resolution is made.

## What is closed, and what is not

The result supplies a concrete joint limit of compact source perturbations
and their finite lattice stress readouts, with derivative control,
normalization, all-valley source mapping and a uniform finite inverse.
It closes the free-sector source-realization gap left explicitly open in
FA91 for these averaged stress responses. Equal-time correlations of
other valley smearings have zero cross-valley limit by the oscillatory
cell estimate in `continuum_field_net`; the new projection estimate also
controls their source filled-sea contractions. Wick's rule then shows
that the shell perturbation changes only the selected continuum sector.

This is not a theorem about nonzero-U composite stress convergence,
regional lattice modular logarithms, boosted modular caps, or a regulator-
independent area/edge entropy. No Newton coefficient, geometric null
balance, generalized-entropy stationarity, or interacting/massive matter
sector is obtained. The original G1--G6 geometry theorem stack remains
the target. A stress detector approaching a modular form at rest does
not by itself produce a geometry equation.

Run the strict package tests and `python -m model.lattice_weyl_response
--write` / `--check`. The retained receipt binds current mathematical
inputs narrowly, checks the new exact controls, and separately reports
inherited broad-custody status. It does not rewrite frozen prior receipts.
