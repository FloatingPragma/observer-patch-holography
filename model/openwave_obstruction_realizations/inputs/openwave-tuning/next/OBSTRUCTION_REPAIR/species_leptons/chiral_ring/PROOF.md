# Full two-cell wall ring with Gauss-law electric dynamics

## Complete source

The carrier has two spatial cells, two regulator slices and both Dirac
coordinates at one fixed physical spin. All eight CAR modes and every
number sector are retained. The source uses h=1/8, Wilson r=1, and the
fifth and boundary mass terms of `chiral_service`. The two bare masses
are zero and the actual432 electron coefficient. The spatial links are
two distinct oriented links, including their returning counterparts.
There is no spatial plaquette in this one-dimensional ring, no magnetic
holonomy pinning term and no vacuum-energy subtraction.

Gauge fixing the interior tree link leaves one compact closing holonomy.
The matter Hamiltonian is `dGamma[A+C exp(i theta)+C* exp(−i theta)]`.
The literal A and C act on all eight modes. A contains the interior
hop; C is the closing hop. At fixed holonomy the full free spectrum is
the union of both momenta `k=(theta+2 pi j)/2`, j=0,1, with all wall bands.
No isolated low momentum replaces this complete sum.

The static references are −2 at each vertex, with an additional charge
−(N−4) at vertex zero in number sector N. If n0 is the matter number at
vertex zero, the charges are `q0=n0+2−N`, `q1=−q0`. Every Gauss state is

```
|occupation,ell>, E0=ell+q0, E1=ell,
H_electric=[(ell+q0)^2+ell^2]/4.
```

An interior hop updates q0. A closing hop also changes ell by one.
Both Gauss laws hold exactly on every retained occupation/flux state.
A fermion insertion changes the reference at zero and uses the fixed
tree dressing for a field at vertex one. Its electric flux change is
included by these formulas. Reference energy is zero by declaration.

Particle-hole conjugation uses Dirac X and reverses ell. It maps
occupation w to the Dirac-permuted complement, N to8−N and q0 to−q0.
The source satisfies `U A^T U*=−A`, `U C^T U*=−C*`; together with flux
reversal these prove an exact full-CAR intertwiner. The executable
checks every number sector, rather than assuming degeneracy from labels.

## Dynamical results and complete currents

Flux cutoffs K=4,8,16 retain every occupation. For K4/8 all finite
eigenstates are retained. For K16 the complete matrices remain in the
calculation, with sparse low eigenpairs and full-space linear resolvents;
the output explicitly does not certify their spectral ordering.
Every case has its lowest computed energy in N4. The complete K8
eigensystem and discarded-flux bounds also separate that sector from
every competing number sector in the infinite rotor.

At zero input mass the K16 results are approximately

```
E0=−48.42841728868,
<cos theta>=−0.83860017822,
charged static-reference transition=7.63866302645.
```

K8 agrees at the displayed precision. Its conservative infinite-flux
transition enclosure is [1.14085317460,13.11222132720], under the stated
finite arithmetic model. The K16 sparse values are diagnostics and do
not replace that certified-order enclosure.

The boundary field has four components: both Dirac components at both
spatial vertices, with all regulator bulk modes retained in its dynamics.
The complete particle and hole spectral matrices at K4/8 sum to I4.
At K16 both charged resolvents are solved on their entire1,848-dimensional
spaces, with residuals below2e−15. No omitted pole is set to zero.
The beta component of `−G(0)^−1` at zero input mass is approximately

```
[[4.04042026716, 0.11198709449],
 [0.11198709449, 4.16773366785]].
```

This differs substantially from the one-cell response. It is a dressed
boundary observable with a specified static reference, not a physical
lepton mass. The electron-input case has transition7.63866885034 and
changes this matrix by its actual source response; no fitted subtraction
or critical mass is applied.

As a diagnostic of the scalar wall hypothesis, assign the lower free
momentum its instantaneous `W_low=1−cos(|theta|/2)`. The event
`W_low<=1/2` is the exact holonomy interval `|theta|<=2pi/3`. Its
probability is evaluated from the complete rotor density matrix by
analytic Fourier integration, including all flux coherences. At m=0
it is approximately0.077691182. Thus the actual vacuum mostly occupies
holonomies near the critical W_low=1 region. This is a test of the
flat-background localization criterion; it is not a decomposition of
the interacting state into free particles.

## Cutoff and arithmetic bounds

The retained non-flux matrix A at zero input mass has exact characteristic
polynomial

```
x^8−512x^6+75776x^4−3670016x^2+26214400.
```

Eight rational isolating intervals give the lower edge of every CAR
number lift by summing its N lowest roots. The two boundary input
contacts change that edge by at most2|m|. The closing operator has norm
at most `8 sqrt(floor((n+1)^2/4))`, n=min(2,N,8−N), including spectators.
For discarded |ell|>=K+1 the exact electric lower bound minimizes
`[(ell+q0)^2+ell^2]/4` over every allowed q0 and the exterior integer
quadratic minima. Therefore

```
H_discarded >= d=E_electric,outside+a_N−2||C_N||.
```

The same two-block Schur lower root used by `chiral_vacuum` bounds each
full sector from its finite ground energy and the boundary coupling
norm. All coarse bounds, including those failing to separate sectors,
are stored. Complete K4/8 eigensystems use full residual/Gram checks
with the declared conservative standard-IEEE error envelopes. This is
not an interval-BLAS claim. K16 sparse ordering is marked uncertified.

The independent review imports no source module. It reconstructs the
one-body integer matrix, verifies its exact polynomial and all eight
Sturm root counts, checks the particle-hole identities, enumerates every
Gauss occupation/flux state and dressed creation update, and checks all
vacuum-record counts and particle/hole CAR sums. All finite-space
operator coefficients, occupations, bulk bands and references remain
part of the source.

## Spatial refinement boundary

The one-cell cutoff failure does not persist unchanged on this ring:
the boundary Green matrix and the charged reference gap both change.
The dynamical two-cell vacuum favours the region predicted by its
complete static filled-band sum. Two cells have no spatial momentum
deep inside the infrared when the twist is near pi.

The separate `wall_ring_control` calculation retains all free momenta
and wall bands for larger rings and bounds their infrared gaps. Its
spatial refinement recovers low free modes; it does not compute this
interacting Gauss-law vacuum. Consequently neither of these small-ring
charged gaps establishes an additive continuum mass or excludes a
protected infrared construction. A larger spatial interacting source,
its complete charged-family vacuum, the common pole clock and the
physical three-dimensional current remain distinct requirements.
