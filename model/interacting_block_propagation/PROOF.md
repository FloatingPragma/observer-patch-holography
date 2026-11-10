# A lattice of interacting blocks with a uniform gap and a moving band

## Declared source and inherited data

One block is the complete eight-cell Peter-Weyl source from
`atomic_spin_geometry`, `collective_spin_response` and
`interacting_ground_preparation`, with native hopping coefficient
`lambda=2^-34`. Its physical Hilbert space is the full joint kernel of
all twelve Gauss generators at every cell. All physical excited sectors
are retained. The larger kinematic carrier and the Gauss projectors are
retained as well; the gap theorem concerns the physical subspace.

Let h be the actual interacting block Hamiltonian. Its ground g is unique
with energy E0. Its first excited eigenspace has rank three and a common
excitation energy Delta. Let Pg and Pt be their exact spectral projectors
on the physical carrier, extended by zero on other Gauss sectors. The
native Fourier currents are

```
O_a=sum_x (-1)^x_a c_x^dag (sigma_a/2 tensor I15) c_x,
a=0,1,2.
```

Their full operator norm is 60. They are even CAR operators and commute
with every local gauge current. The actual three first-excited momentum
sectors are X_a. The source translation selection rule makes
`Pt O_a g` lie in its one-dimensional X_a sector. Cubic symmetry makes
their squared norms equal. The frozen collective-response theorem proves
the positive interval

```
r = ||Pt O_a g||^2/3600
    in [0.00204894419..., 0.00211090162...],
Delta in [8.20758426381e-22..., 8.21034719056e-22...].
```

The exact rational endpoints, rather than these displayed decimals, are
consumed by `model.py`. All input hashes of all three frozen
parents are validated. The native current, symmetry and full-occupation
operator checks are independently replayed.

Define the explicitly dressed native lowering operators

```
B_a=Pg O_a Pt/60.
```

Choose the phase of the exact normalized excited state e_a by
`e_a=Pt O_a g/(60 sqrt(r))`. Then

```
B_a=sqrt(r)|g><e_a|,
sum_a B_a^dag B_a=r Pt.
```

This is an exact identity on the full block carrier: B vanishes on every
other sector. Spectral dressing and the ensuing inter-block interactions
are declared architectural operations. Exact finite spectral projectors
can be specified by the complete finite block Hamiltonian and its isolated
clusters. No efficient native synthesis is asserted. This construction
is distinct from placing the original unmodified hopping between cells
of a uniform side-four source.

For a finite block graph of maximum degree z, use identical block
orientations and set

```
H_Lambda=sum_x h_x
  -kappa sum_<xy>,a (B_xa^dag B_ya+B_ya^dag B_xa).
```

The physical source is the tensor product of the full block Gauss
subspaces, with clock `exp(-itH_Lambda)`. Inter-block terms commute with
every block Gauss generator, preserve parity and raw L, and are supported
on two complete blocks. There are no inter-block charged gauge links in
this declared source. A global cubic rotation rotates both block positions
and the native vector index a, preserving the contracted interaction.

## Volume-uniform operator bound

Put `Q_x=I-Pg_x` and `J=kappa r`. Normalize the lowering operators as
`A_a=B_a/sqrt(r)`. For any two blocks and any vector v,

```
|<v,sum_a(A_xa^dag A_ya+A_ya^dag A_xa)v>|
 <= sum_a (||A_xa v||^2+||A_ya v||^2)
 = <v,(Pt_x+Pt_y)v>.
```

This proves both operator inequalities, including all higher block
sectors. The executable verifies their exact finite matrix structure with
the complete ground/triplet space and an untouched higher-sector control.
The exchange has norm one on the normalized two-block low space. The
inequality does not grow with the dimension of the untouched complement.

The frozen block theorem gives `h_x-E0>=Delta Q_x`. Consequently

```
H_Lambda-|Lambda|E0
 >= (Delta-zJ) sum_x Q_x.
```

The inter-block exchange annihilates the product ground
`Omega_Lambda=tensor_x g_x`. If `Delta-zJ>0`, it follows that this is the
unique finite-volume physical ground and the full gap is at least
`Delta-zJ`. No expansion in the total perturbation norm is used. In
particular, the bound is uniform in the number of blocks.

For the cubic block graph, z=6. Let Delta_- and r_+ be the inherited
rational lower gap and upper residue bounds. The declared source-only
coefficient is

```
kappa=Delta_-/(4 z r_+).
```

Then `J<=Delta_-/(4z)` and the full gap is at least
`3 Delta_-/4`. No measured value selects kappa. The exact value and its
source-unit bounds are in `results.json`.

The total number of excited blocks `sum_x Q_x` commutes with H. Higher
block sectors are retained and are untouched by the exchange at their
occupied site. In sectors with two or more excited blocks, the energy is
at least `2(Delta-zJ)`. A single higher block excitation has energy at
least the minimum of the frozen next-low-band threshold and the frozen
gap to the full complementary atomic band. The receipt verifies that
both lower bounds exceed `Delta_+ +zJ_+`. Thus the entire three-polarization
one-excitation band is isolated from the ground and every other physical
source state, uniformly in volume.

## Exact moving excitation and complete records

The vectors

```
|x,a>=e_a at block x, g at every other block
```

span an exactly reducing subspace of dimension `3|Lambda|`. Restricted
to it, the centered full Hamiltonian is

```
H_1=Delta I-J A_graph, independently for each a.
```

This is a consequence of the actual dressed native operator, rather than
a fitted hopping Hamiltonian. On an L-by-L-by-L periodic block lattice,
all Fourier vectors are exact eigenvectors, with

```
E_a(k)=Delta-2J[cos(k1)+cos(k2)+cos(k3)],
k_i=2*pi*n_i/L,  n_i=0,...,L-1.
```

The code checks every Fourier mode at L=3,4,6, corresponding to 216,
512 and 1728 physical source cells. Its diagonalizations are restricted
to this exactly proved reducing subspace. They do not claim numerical
execution of the full exponentially large block carrier.

The finite source correlation has a literal moving spectral pole:

```
<Omega| B_a(k,t) B_b(k)^dag |Omega>
 = delta_ab r exp(-it E_a(k)).
```

The residue is positive by the frozen theorem. The dressed operator has
zero spectral weight in the full outside complement by its definition;
the undressed native current retains the separately certified outside
weight. Block momentum k and the internal block momenta X_a are distinct
observables. The three internal momenta are not fermion generations.

A complete localized pump uses the dressed vector operator
`sum_a B_a^dag tensor |a> /sqrt(3)` at one block, with its Halmos-unitary
completion. On the exact product ground, success probability is r and
failure probability is `1-r`. The reference keeps the native vector axis
and inverse internal momentum. After success the full block-lattice
Hamiltonian moves the excitation. Position y and reference axis a have
unconditional probability

```
(r/3)|[exp(iJt A_graph)]_(y,0)|^2.
```

The common phase `exp(-iDelta t)` cancels in these probabilities. Every
position, every axis and the failure branch are retained. Their sum is
`1-r+r=1`. The code computes complete finite position distributions at
the source-defined clock reading `Jt=1/4`.

## Volume limit and propagation error

The finite-volume ground marginals are exactly compatible product states,
so their local ground-state limit exists without an approximation. They
have no inter-block ground entanglement; all interacting correlations
inside a block are retained. The gap statement concerns every finite
volume and the finite-excitation representation of this product limit.

For a localized one-excitation vector at distance R from a removed
boundary, the finite and infinite adjacency powers agree through order
R-1. Their norms are at most z^n. Taylor expansion of the two propagators
therefore gives the explicit norm bound

```
||U_finite(t)|0>-U_infinite(t)|0>||
 <= min(2, 2 sum_(n>=R) (zJ|t|)^n/n!).
```

It also bounds the total variation of any complete one-excitation
detector. A convenient exponential upper bound is
`min(2,2 exp(e zJ|t|-R))`. The receipt gives both bounds for several R.
This is a quantitative locality and volume-limit statement for the
actual moving excitation sector, with the same Hamiltonian clock.

General many-body locality theorems also have applicable finite-range,
uniform local-dimension hypotheses for this fixed-block architecture.
Their finite-velocity bounds do not imply Lorentz symmetry. No velocity
or spectral limit from another regulator is substituted here.

## Scalable finite preparation

The frozen constructive preparation uses the actual native local ground,
the cube singlet and a complete finite energy-agreement filter. Let p0
be its certified one-block trial overlap and epsilon its finite mixed
ground/excited filter leakage. Both are taken from the frozen receipt;
the filter length is unchanged.

If a source marginal has ground probability p, apply the same agreement
filter to two independently prepared copies. The proof remains valid for
mixed source marginals, since the ground/ground, mixed and excited/excited
spectral sectors are orthogonal and the filter preserves each. On success,
its one-copy excited-to-ground odds satisfy

```
r_next <= r^2+2 epsilon^2 r,  r=(1-p)/p.
```

Set `r0=(1-p0)/p0` and `c=r0+2epsilon^2<1`. Induction gives
`r_n<=r0 c^n`. This is a finite binary tree of the same declared filters,
using `2^n` source trials. All unused copies, spin references and failure
records remain in the output. It uses no interacting-ground-energy input.

For M independently prepared blocks, conditional on all tree successes,
the product ground fidelity is at least `(1+r0 c^n)^(-M)`. Hence its
trace distance from the exact product ground is at most

```
sqrt(M r0 c^n).
```

Choosing the least n with `M r0 c^n<=epsilon_global^2` gives an explicit
finite volume-dependent preparation. The receipt uses
`epsilon_global=2^-20` and retains the lower probability of all successes.
It does not claim efficiency or a volume-independent success rate.
After evolution by H_Lambda, the same trace-distance bound holds for all
real times because the target product ground is stationary. All complete
detector distributions differ by at most this error. A fixed one-block
accuracy is not silently treated as a global volume-independent accuracy.

## Spacing refinement and the nonrelativistic limit

Assign block spacing a and let p denote physical momentum, so k=ap. The
bottom energy is `m0=Delta-6J`. For bounded p,

```
E(ap)=m0+J a^2|p|^2+error,
|error|<=J a^4 sum_i p_i^4/12.
```

This follows from the elementary cosine Taylor remainder and holds on
the complete finite momentum menu in the stated bounded window. In the
fixed source clock, the kinetic coefficient vanishes as a approaches
zero. Keeping a nonzero coefficient by `J proportional a^-2` violates
the fixed-block gap control.

After subtracting m0 and using time `t=eta/(J a^2)`, the one-particle
generator tends to `|p|^2`, with operator error on a bounded momentum
window at most `a^2 sum_i p_i^4/12`. The phase error is bounded by this
number times `|eta|`. This is a controlled Schrödinger limit. Its
dispersion is not `sqrt(m^2+c^2|p|^2)`, and no Lorentz covariance or
relativistic spin assignment is proved. Matching the first quadratic
term of a massive relativistic dispersion would not establish that
stronger claim.

## Audit of theorem reuse

The eight-cell quartic remainder uses its own bipartite native parity.
The frozen proof explicitly excludes its extension to the even-displacement
paths of side four. A bound based on the total hopping norm would also
grow with volume. Neither estimate is used for inter-block gap control;
the exact number inequality supplies the uniform result.

Yarotsky's [quasiparticle theorem](https://arxiv.org/html/math-ph/0411042)
assumes a uniformly gapped nondegenerate product reference and a sufficiently
small finite-range perturbation, measured relative to that gap. Its stated
existence constants do not by themselves certify a particular source
coefficient. The present spectral dressing permits an explicit bound and
an exact band without importing such an unspecified threshold.

Nachtergaele and Sims's [locality and clustering theorem](https://arxiv.org/pdf/math-ph/0506030)
uses a uniform finite site dimension and a controlled interaction norm.
It supplies locality and gap-to-clustering implications, rather than a
relativistic dispersion law. Those hypotheses apply to fixed full blocks;
they would need independent uniform control if the block dimension,
interaction norm or spacing changed. The direct one-excitation tail
bound above states the actual approximation error used here.
