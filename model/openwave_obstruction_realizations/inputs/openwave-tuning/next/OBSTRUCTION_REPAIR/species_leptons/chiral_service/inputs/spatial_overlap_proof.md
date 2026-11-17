# Free spatial propagation and interacting quantum flux

This package checks two distinct finite controls: a three-dimensional free
spatial overlap Hamiltonian, and a quantum electric/loop/matter Hamiltonian
on a four-edge abelian graph. Their dynamics, preparations, cutoffs and
measurements are explicit. They do not form a single interacting chiral
continuum theory. The lattice spacing and electric-flux cutoff have no
derived identification with OPH capacity parameters or physical units.

## Source representation and conditional coordinate completion

The native integer charge blocks are `(q=6Y,dimension)=(-4,3),(-3,2),
(1,6),(2,3),(6,1)`. The sixty-state link is the vacuum plus the sum of
the five endpoint spaces `V_r* tensor V_r`. The existing
[quantum-link proof](../../particle_gauge_band/QUANTUM_LINK_PROOF.md)
derives its currents, finite coefficient products and conditional
compact-group realization. Those results are reused here.

For `G=S(U(3) x U(2))`, Schur orthogonality gives the isometry

```
J|0> = 1,
J|r;a,b> = sqrt(d_r) D^r_ab(g) in L²(G,Haar),
P_r M[D^r_ab] P_0 = J U_native,ab J* / sqrt(d_r).
```

The existing finite-menu operator is the rescaled compression
`sqrt(d_r) P_menu M[D^r_ab] P_menu`, whose vacuum amplitude is the
literal native raising amplitude. Multiplication has additional
representation transitions and is distinct from the raising operator.

The exterior representation is faithful: `Lambda²(U)=I` forces all five
eigenvalues to equal `z`, with `z²=z⁵=1`, hence `z=1`. Its matrix
coefficients and conjugates separate group points, so their unital algebra
is dense in `C(G)`. Given that coordinate completion and a normalized
endpoint-invariant cyclic blank, translation invariance selects Haar
integration. Selection of this completion by repair is not proved.

No finite-dimensional unitary charged shift obeys `U*EU=E+qI` for nonzero
`q`: taking the trace gives `q dim(H)=0`. This excludes that exact algebra,
not every finite-link continuum route. Representation cutoffs preserve
endpoint covariance while retaining boundary transport defects.

On a fixed finite cubic spatial lattice the coordinate completion admits
a conditional vectorlike Hamiltonian. The invariant positive metric
`g_ab=-tr(T_a T_b)` defines an electric Casimir Laplacian. A plaquette
trace potential is bounded. Spatial overlap has one-body norm at most
`2/a`; its lift to the full finite fermionic Fock space is also bounded.
The Wilson determinant is real analytic and nonzero at the flat
configuration. Its singular set therefore has Haar measure zero; a
bounded definition there gives the same multiplication operator almost
everywhere. Electric energy plus these bounded self-adjoint terms is
self-adjoint and gauge covariant. This supplies no uniform locality near
small Wilson gaps and no interacting continuum limit.

## A nilpotent one-bond test does not decide propagation

For `A=bI-B`, `B²=0`, the normalized native transport singular values are
`s_r=sqrt(d_r/6)`. The active spin/endpoint blocks satisfy
`D*D=2+4b/sqrt(4b²+s_r²)`. At `b=1`, the minimum singular value is
approximately `1.9464979789` and the 7,080 complementary states have
`D=2I`. At `b=-1`, that complement has `D=0`. At `b=0`, the required
Wilson gap fails. A sign change alone does not supply spatial transport.
The three-dimensional Wilson kernel has on-site coefficient `3-rho=2`
and three returning translations, giving `A(0)=-1`.

## Complete free spatial overlap control

Use source Hermitian gamma matrices, `beta=gamma4`, `alpha_i=i beta gamma_i`,
`k_i=a p_i`, and

```
u_i=1-cos k_i, B=sum_i u_i-1,
R²=B²+sum_i sin² k_i=1+2 sum_(i<j)u_i u_j >= 1,
V=(B+i sum_i gamma_i sin k_i)/R,
D=(1+V)/a, h=beta D.
```

Since `beta V beta=V*`, `h` is Hermitian. Its two positive and two
negative energies obey `E²=2(1+B/R)/a²`. Only `k=0` is zero; the other
Brillouin-zone corners have `E=2/a`. The modified axial charge obeys

```
q5=gamma5(1-V)/2, [h,q5]=0, q5²=1-a²h²/4.
```

It approaches `gamma5` near zero momentum and vanishes at cutoff corners.
The vectorlike pair is retained. This gauge-dependent charge generally
does not commute with electric energy. It is distinct from the Euclidean
Weyl projector in the adjacent electric package.

At box length `2pi`, momentum `(1,1,1)`, and `a=2pi/L`, the energies for
`L=32,64,128,256` are approximately `1.7448298797`, `1.7354582753`,
`1.7329158901`, `1.7322679027`, approaching `sqrt(3)`. The source free
bound `||h-alpha.p||<=a|p|²` applies when `a|p|<=1/4`. To avoid
cancellation for `B<0`, evaluate
`1+B/R=sum_i sin²(k_i)/(R(R-B))`.

The `L=8` control retains 512 momenta, four spin components and fifteen
internal components: 30,720 states. Its sixty zero states are the
spin/internal multiplicity at one momentum, not identified species.
The independent source check projects the full 324-dimensional source
Wilson/overlap matrix onto all 27 spatial momentum blocks at temporal
momentum zero.

For `E>0`, prepare `rho=P_negative/2` and measure `Pi_±=(I±beta)/2`
twice. All four outcomes remain in

```
p(a,b;t)=tr(Pi_b exp(-iht) Pi_a rho Pi_a exp(iht) Pi_b).
```

For `b_h=(1+B/R)/(aE)`, `Z=1-b_h²`, the recorded correlation is
`b_h²+Z cos(2Et)`. The four sequence values of `Z` are approximately
`.9706569017`, `.9927428157`, `.9981910158`, `.9995480921`. The
recorded transition frequency is `2E`. The positive-projector inelastic
and elastic weights are `Z/4` and `(1-b_h)²/4`, summing to its mean
`(1-b_h)/2`. A separate localized preparation with spin `I4/4` retains
all 64 position outcomes on the side-four lattice. Its probability away
from the origin at time `.5`, spacing one, is approximately `.09374198848`.

A fixed flat gauge configuration has Haar measure zero and is not a
normalizable fluctuating-gauge state. These free-fiber experiments do not
establish survival of the response in a source-selected interacting state.

## Complete abelian quantum-flux control

Restriction to the native integer `U(1)` subgroup gives its own
`L²(U(1))` realization, not a reducing subspace of `L²(G)`. Nonabelian
internal components become identical charge copies. For mobile charge
`q` at site `x` and static conjugate charge at site zero, the complete
Gauss basis is `|x,ell>`, with integer `ell` and

```
E_e(x,ell)=ell+q 1_(e<x),
G_v=q delta_(v,x)-q delta_(v,0)+E_v-E_(v-1)=0.
```

The forward hop `T` increments `x` and adds `q` to `ell` only on wrapping
to zero. The loop `W` increments `ell`. Thus `T⁴=W^q`, `[E_e,T]=q T Q_e`
and `[E_e,W]=W`, with complete position projector `Q_e`. The Hamiltonian is

```
H=sum_e E_e²/2-(T+T*)-(W+W*)/2.
```

The coefficients and preparations are declared controls. The loop is
the primitive charge-one loop, not a selected full-matter trace action.
Cutoff `|ell|<=K` retains every `4(2K+1)` state and compresses transitions
leaving it. Gauss and electric Ward identities remain exact; the boundary
transport unitarity defect has norm one.

All five charges and cutoffs `4,8,12,24,48` are executed. Complete
spectral arrays, including zero weights, can be exported. At `q=1,K=48`,
the ground-state position-zero inelastic weight is approximately
`.2324598386713`. Its first line above the declared `1e-10` residue
threshold has gap `1.681493978167` and weight `.1297411844112`.
At `K=24`, turning the loop coefficient from zero to one shifts the
ground energy by approximately `-.3768864726971`. This is a control
comparison, not a measured mass correction. The sixteen four-round
binary histories instead start in `|0,0>`.

## A cutoff bound for complete recorded experiments

The diagonal electric term is unbounded. The other interaction is bounded
by `B=2|J|+|kappa_B|`; an insertion changes loop flux by at most
`m=max(1,|q|)`. From initial support `|ell|<=K0<=K`, leaving the cutoff
requires `r=ceil((K-K0+1)/m)` insertions. The full and truncated
interaction-picture Dyson series agree at smaller orders, so

```
||psi(T)-psi_K(T)|| <= 2 sum_(n>=r)(BT)^n/n!.
```

Complete position instruments preserve flux support. Their unitary record
dilations obey the same estimate with total unitary time `T`, including
untouched references. Final measurement contracts trace distance, bounding
total variation of the entire recorded law by this bound, capped at one.
For all native charges, `m=6`, `K=48`, `K0=0`, `T=.5`, `B=3`, the
series is approximately `.0002486066359511`. It does not bound a ground-
state preparation without a separate initial-state approximation bound.

Quadratic electric growth and bounded interaction give compact resolvent.
Nested finite subspaces give variational convergence. A diagnostic lower
ground-energy bound minimizes the discarded electric diagonal, subtracts
`B`, then takes the lower eigenvalue of the two-by-two block bound with
off-diagonal norm `B`.

These inequalities are analytic statements. Matrix residuals, decimals and
evaluations of the bounds are floating-point diagnostics, not interval
certificates. The committed certificate checks fresh calculations against
explicit tolerances and stores stable admissions and source hashes.

The spatial-overlap/electric distinction agrees with
[Creutz, Horvath and Neuberger](https://arxiv.org/abs/hep-lat/0110009).
No common interacting chiral continuum, physical clock, energy scale,
source-selected vacuum or particle-mass identification is claimed.
