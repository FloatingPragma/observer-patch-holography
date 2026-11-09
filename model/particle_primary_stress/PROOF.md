# Local energy transport on the selected physical source

## Source and the reused finite theorem

Every member of the [primary contract](../particle_common_contract/PROOF.md)
has one open Wilson Hamiltonian H, sixty CAR modes per cell, three owned
outgoing link factors, the product word cutoff, the selected physical
polynomial density `rho=Q²/Tr(Q²)`, and the retained conjugate reference
with idle Hamiltonian zero. No action coefficient or energy zero is changed.
The clock is `alpha_t(A)=exp(itH) A exp(-itH)`.

The operator continuity and controlled-work theorem is
[the same-program energy theorem](../particle_gravity_response/PROOF.md#energy-and-work-from-one-hamiltonian).
The finite spatial and temporal sum identities are
[`CoarseEnergyWard.lean`](../../Lean/QFT/CoarseEnergyWard.lean).
This package supplies their actual primary cell partition, local supports,
source bounds and complete readout domain. It does not encode the primary
matrices as Lean values or repeat the finite sum proof.

## A local partition of the actual open action

Let `M_x=dGamma(3 beta_D/a)` be the onsite matter term, `E_(x,i)` the
complete electric energy of the owned outgoing link, `B_b` an interior
Hermitian Wilson bond and `P_p` an interior negative plaquette-trace
potential with its original energy zero. A plaquette based at x in axes
`i<j` has owner support `O_p={x,x+e_i,x+e_j}`. Its fourth corner owns none
of its four links. Define

```
h_x = M_x + sum_i E_(x,i)
      + (1/2) sum_(interior bonds incident to x) B_b
      + (1/3) sum_(interior plaquettes with x in O_p) P_p.       (1)
```

Every source term is allocated with total coefficient one, so `sum_x h_x=H`
exactly. Dangling outgoing links retain their electric terms; hopping and
plaquettes requiring exterior cells are absent. Ghost-endpoint Gauss
constraints remain imposed on the state. Regional energy means the sum
of this fixed allocation, including its shares of crossing interactions.
It is not a fresh isolated-box Hamiltonian or a density with interaction
energy removed.

The [primary locality proof](../particle_primary_locality/PROOF.md) proves
that this product link compression preserves each complete parent term's
support and cannot increase its norm. All terms in (1) are even and gauge
invariant before restriction to the physical subspace. The global Gauss
projection is not inserted as a local interaction. Since each h_x commutes
with that projection, the same operator identities restrict to physical
inputs without a factorization of the physical Hilbert space by region.

Every owner in `supp h_x` lies within Manhattan distance two of x. Even
operators on disjoint CAR/link supports commute, hence

```
J_xy = i[h_x,h_y] = J_xy† = -J_yx,
J_xy = 0 if |x-y|_1 > 4.                                  (2)
```

The radius four is a sufficient upper bound. The executable keeps every
pair whose actual energy supports overlap; overlap is not a claim that
the corresponding commutator is nonzero. A region retains all crossing
pairs, including those connecting cells more than one lattice step apart.
No nearest-face routing is inferred from the energy allocation.

## Norms and exact continuity at every finite level

The native beta matrix has thirty positive and thirty negative one-body
levels in the sixty-mode carrier. Thus `||M_x||=90/a` on the full Fock
space. The inherited full-Fock bond bound is `||B_b||<=30/a`, and
`||P_p||<=1/(7a)`. The native trace Casimir bound is `(29/5)ell²`; each
owned electric energy therefore has norm at most `609 ell²/(2a)`.
There are at most six incident bonds and nine owner plaquettes. Consequently

```
||h_x|| <= K_n = [180+(1827/2)ell²+3/7]/a,
||J_xy|| <= 2 K_n².                                      (3)
```

`cell_energy` returns sharper boundary-specific counts and the actual
term descriptors. `incidence_control` independently enumerates a small
open box, compares every allocated term with an independent action term
inventory, and checks the total share one for each term. At side three,
there are 27 matter terms, 81 electric terms, 54 bonds and 36 plaquettes.
The largest cell norm cap is `15315/14` at `a=ell=1`.
The bound (3) grows with cutoff and inverse spacing; it is not a uniform
continuum stress estimate.

The reused operator theorem gives

```
d alpha_t(h_x)/dt + sum_y alpha_t(J_xy) = 0,
H_R = sum_(x in R) h_x,
F_R = sum_(x in R,y outside R) J_xy,
alpha_T(H_R)-H_R + integral_0^T alpha_t(F_R) dt = 0.       (4)
```

Internal pairs cancel by antisymmetry. The whole source has zero outward
flux because all terms belong to that source; this does not omit flux
across a proper regional cut. Taking any joint state with an untouched
reference preserves (4). For the selected stationary source density,
`Tr rho F_R=0`; individual pair currents are not asserted to vanish.
The joint purification is generally not stationary under the source-only
clock. The operator law also applies after physical interventions, when
the prepared source marginal need not be stationary.

Setting exact endpoint increments and integrated pair fluxes from (4)
supplies `CoarseEnergyWard.coarse_continuity`. Arbitrary block maps and
successive time intervals then use `coarse_compose` and
`time_increments_telescope`. An endpoint current is not substituted for
an integrated current without an error bound.

## Complete records and a finite measurement error bound

For any nonzero Hermitian current J and a declared `M>=||J||`, `M>0`,
the effects

```
E_+=(I+J/M)/2,  E_-=(I-J/M)/2                             (5)
```

are positive, gauge invariant and complete. Labels `+M,-M` have mean
`Tr rho J`. The finite square-root instrument and its full pointer
isometry retain both outcomes and the original conjugate reference.
A zero current can use one identity outcome labeled zero. Equation (5)
is an unbiased two-outcome current instrument; its labels are not claimed
to be the spectral eigenvalues of J. The same construction applies to
regional energy and total outward flux. Full spatial matrix construction,
logarithm computation and apparatus execution are not claimed.

`regional_bounds` enumerates every potentially nonzero crossing pair for
the specified finite region without allocating the source Hilbert space.
Let its norm caps be `K_R>=||H_R||`, `M_R>=||F_R||`, and let `B>=||H||`
come from the selected preparation's actual spectral interval. A complete
output half-trace error at most delta, uniformly at every sampled time,
gives mean errors at most `2K_R delta` and `2M_R delta`. This comparison
includes preparation, evolution and detector errors and arbitrary retained
reference correlations; no outcome normalization is allowed.

For N equal subintervals of `[0,T]`, left-endpoint sampling of the flux
has quadrature error at most `T² B M_R/N`, since
`||d alpha_t(F_R)/dt||<=2 B M_R`. Thus the measured residual in (4) obeys

```
|measured endpoint increment + sampled integrated flux|
 <= 4 K_R delta + 2 T M_R delta + T² B M_R/N.            (6)
```

These endpoint and flux menus are alternative experiments on separate
copies of the same initial preparation. They do not reconstruct an
undisturbed trajectory from sequential measurements. If apparatus acts
during a single history, its disturbance and supplied work enter that
history's energy balance and complete instrument law.

The package returns this rational bound at the original source time.
Its large finite constants do not establish an efficient apparatus or a
refinement-uniform measurement of continuum stress.

## A coefficient response at fixed preparation

The declared control `H_f(t)=sum_x f_x(t)h_x`, with real differentiable
coefficients, has local energies `f_x h_x`. The existing work theorem gives

```
d< f_x h_x >/dt + sum_y f_x f_y <J_xy> = dot(f_x)<h_x>,
d<H_f>/dt = sum_x dot(f_x)<h_x>.                         (7)
```

At all f equal to one this is exactly the primary action. At a fixed
initial density, `partial H_f/partial f_x=h_x`; no perturbed ground state
is substituted. A varying f declares an external control and its work.
It does not define proper time or make f an observer-measured lapse.
The [constant-coframe action](../particle_geometric_response/PROOF.md)
has its own explicit deformation couplings and six strain operators.
Those additional couplings, and any identification with a proper metric,
do not follow from (1) or (7). This package asserts no coframe stress
identity by renaming a source coefficient.

## Regional entropy on the same state and its precise marginal

The [selected regional support consumer](../particle_spatial_modular_response/PROOF.md#3-physical-gauss-support-and-the-boundary-algebra)
gives `rho>=2^-N P_G` for a source-derived integer N. On a region of r
owner cells, `D_R<=2^b_R`, where `b_R=(60+27ell)r`. Its actual boundary
representation has supported partial-trace blocks `n_lambda/d_lambda`,
with integers `n_lambda>=1` and `d_lambda<=D_R`. Consequently

```
sigma_R >= m_R P_R,   m_R=2^-(N+b_R),
P_R = support Tr_out(P_G),
||K_R|| <= (N+b_R)log(2) <= N+b_R,
K_R = -log(sigma_R)|P_R.                               (8)
```

The floor applies only on P_R. Every physical boundary-center sector is
retained; the region is not projected to a separate Gauss singlet.
No ambient full-rank or limiting faithful-vacuum claim is made.

Take R to be the actual owner support of h_origin and
`U_R=exp(i theta h_origin)`, with `theta=1/1000`. Since `U_R tensor I`
commutes with P_G, U_R preserves P_R. Apply the complete instrument
`sqrt(1-s) I, sqrt(s) U_R`, retaining its classical coin and the conjugate
reference. The source regional marginal is

```
sigma_s=(1-s)sigma+s U_R sigma U_R†,
X=U_R sigma U_R†-sigma,  Tr X=0.
```

The [existing finite first-law theorem](../particle_spatial_modular_response/PROOF.md#4-modular-operator-and-finite-entropy-remainder)
on this fixed support supplies

```
S(sigma_s)-S(sigma) = s Tr(X K_R)-D(sigma_s || sigma),
0 <= D(sigma_s || sigma) <= s² ||X||_HS²/(2m_R) <= s²/m_R. (9)
```

Both endpoints and the entire mixture have floor m_R, and
`||X||_HS²<=2` for the difference of two densities. The receipt keeps
this potentially enormous constant factored; for
`s=2^-(N+b_R+4)`, it gives `D<=2^-(N+b_R+8)`.
There is no assumption that the response is nonzero. The logarithm and
entropy refer to the actual marginal of the selected polynomial state,
not an exact Gibbs replacement or `beta H_R`.

Equation (9) concerns the source marginal after a mathematical trace over
the coin and reference. Those registers remain in the complete observer
output. In particular the retained classical coin plus region has entropy
`H2(s)+S(sigma)`, whose derivative at zero is singular. Its entropy is not
the differentiable marginal functional in (9). Neither expression is
identified with Einstein generalized entropy or proper-area response.

## Executed independent quantum-link control

The exact control has two cells, a positive and a negative CAR species per
cell, and one link with electric labels `-1,0,1`. All 48 ambient states
are retained. The forward interaction is
`c_(0,+)† S c_(1,+)/3 + c_(0,-)† S† c_(1,-)/5`; add its adjoint.
Write N_x for the total occupation on cell x, E for the link electric
operator and B for this Hermitian bond. The action and partition are

```
H=N_0+2N_1+(3/2)E²+B,
h_0=N_0+(3/2)E²+B/2,   h_1=2N_1+B/2,
G_0=n_(0,+)-n_(0,-)-E,
G_1=n_(1,+)-n_(1,-)+E.
```

Both cell energies and their current commute with both Gauss generators.
Their common zero-charge support P has dimension six. The positive
amplitude `Q=P(10I-H)` gives a stationary density `rho=Q²/Tr(Q²)`;
`||H||<=257/30<10` ensures faithfulness on P. This is a separate quantum
control, not a compressed sector of the primary sixty-mode program.

The inherited exact partition checker verifies all operator currents,
weighted-source continuity and complete regional boundary flux. A local
phase `exp(i pi n_(0,+)/2)` acting on this density gives

```
<J_01>=-3600/124697,    d<h_0>/dt=3600/124697.
```

The current norm cap `M=6/5` produces both exact probabilities in (5):
`121697/249394` and `127697/249394`. They sum to one and recover that
mean with labels `+6/5,-6/5`. The stationary input has zero net current,
while its idle-purifier energy variance is strictly positive. All matrix
calculations include both species and every retained link label.

## The additional Einstein premises

The [matched regional theorem](../matched_einstein_limit/PROOF.md) requires
observer-defined proper geometry and curvature, a symmetric action-derived
Lorentz stress with covariant conservation, the same-state diamond modular
kernel, fixed-proper-volume area response and generalized-entropy balance,
and a common shrinking-ball limit with errors `o(radius^4)` and a positive
spacetime-independent area coefficient. The finite tensor decoder and
its conditional Einstein implication do not construct these inputs.

Equations (4) and (9) supply local source-energy conservation and a scoped
finite marginal entropy pairing on the selected preparation. They identify
neither proper geometry nor the other spacetime stress components. The
partial trace and its boundary center are specified by source ownership,
not by an independently measured geodesic ball. No gravitational coupling,
Einstein dynamics or continuum stress follows from this finite consumer.
