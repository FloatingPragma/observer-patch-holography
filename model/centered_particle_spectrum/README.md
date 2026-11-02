# Complete visible charged spectrum of the centered cube

The centered eight-site cube has **228 distinct positive energies in every
local one-spin addition or removal measure**, with every weight rigorously
positive. The signed retarded measure has 456 atoms. All eight-site residue
matrices are enclosed, and each has visible rank four per hidden spin.
The result is a complete local Lehmann resolution, including weights below
`10^-20`; no threshold discards a line.

At `c=1/8,U=1`, the certified reference and extreme visible energies are

| Quantity | Certified interval in source-clock energy units |
| --- | --- |
| Centered ground energy | `[-4.12405873337830,-4.12405873337828]` |
| Smallest visible excitation | `[0.36168271180597,0.36168271180599]` |
| Largest visible excitation | `[3.88643475495059,3.88643475495061]` |
| Conditional mean, either charge branch | `[0.49955163274809,0.49955163281375]` |

The first line has weight approximately `0.1329459819237663`. The smallest
weight is approximately `5.63e-21`, at energy `3.70680063381435`, and its
certified lower bound exceeds `4e-21`. Each charge branch has total weight
`1/2`. [lines.csv](lines.csv) contains every scalar energy and weight
interval; [results.json](results.json) also contains every spatial residue
entry, the exact error certificates and the conditional comparison.

These are finite excitation energies of the **new centered reference**.
They are not the 78 lines of the original zero-chemical-potential, `N=4`
ground. They are not relativistic rest masses. No dimensional calibration,
physical species assignment or native observer run is supplied here.
The visible energy spread is less than eleven. Consequently no choice of
three distinct lines, with one positive multiplicative unit factor, can
reproduce the fixed charged-lepton mass ratios if these gaps are proposed
as their rest energies. That conditional test excludes this finite gap
dictionary; it does not exclude unseen composites, a different source or
refinement, or OPH as a whole.

## Source, reference and imported ground theorem

The complete source is the existing open `2x2x2` Weyl/contact cube. Its exact
onsite hidden-spin frame gives

```
H = dGamma(A tensor I_2)/16 + sum_x n_xup n_xdown,
K = H-N/2,             A=A^T, A^2=3I_8.
```

All source matrices are extracted from the existing full 65,536-state
integer Hamiltonian. No hopping, interaction, observer clock or scalar
energy offset is fitted. The centering changes charged dynamics and uses
its own ground state.

The [centered-ground theorem](../centered_source_ground/README.md) proves
that every full ground has `N=8`; its proof uses the existing all-number
operator bound and a dimer variational state. Its explicitly imported
[Lieb theorem, Theorem 1](https://arxiv.org/html/cond-mat/9311033v2), gives
uniqueness and hidden spin zero for the connected, real, balanced bipartite
source with positive repulsion. Real hopping signs are permitted: the
review does not assume a constant positive edge hopping. We reuse this
proved ground selection, not a Ritz guess. `K` commutes with eta raising
and lowering. A nonzero eta ladder of its unique sharp-`N=8` ground would
be another ground in a different charge sector. Hence the ground has eta
zero as well.

The packet checks current source symmetries and the imported theorem's
finite source hypotheses directly. It separately reports historical
receipt pin discrepancies. Such discrepancies are not silently repaired,
and historical receipt authentication is not included in this packet's
`all_pass` scope. The old all-number operator theorem and Lieb theorem are
identified mathematical inputs, with their source/proof files bound.

## Exact spatial and spin/eta reductions

Each of the 48 cube automorphisms has a normalized real signed lift `R`
with `RA=AR` and site-zero gauge sign `+1`. The generated signed group has
96 elements. On even particle number its central sign acts trivially;
on odd particle number it can act as minus the identity.

The half-filled ground is invariant under all these spatial lifts. To
prove this, vary the repulsion continuously from zero to one in the fixed
`N=8` sector. At zero repulsion the nonzero kinetic band gap makes the
filled negative band unique. The spatial determinant on that occupied
band is squared across the two spins, so its character is
`det(R_negative)^2=1`. Lieb uniqueness holds for every positive repulsion.
A simple finite-dimensional ground projection is continuous; its character
under each finite-order symmetry cannot change. This proves the spatial
character at the target coupling without using numerical overlaps.

`invariant_orbits` constructs exact signed orbit columns `Q`. They have
disjoint supports, so `Q^T Q=D` is a positive diagonal integer matrix.
An orbit with inconsistent signs has no invariant vector. Every remaining
orbit is retained. The full spatial invariant balanced ground block has
dimension 120.

For the removal seed `c_0down G`, the hidden spin weight is `+1/2`, eta
weight is `-1/2`, and both `S_plus` and `eta_minus` annihilate it. It is
also fixed by the six-element origin stabilizer. The complete stabilizer
invariant `(Nup,Ndown)=(4,3)` block has dimension 656. There is no choice
of a particular small amplitude in this reduction.

Let `L=S_minus S_plus` or `eta_plus eta_minus`. Exact polynomial projectors
onto their kernels use the nonzero roots

```
ground weight: {2,6,12,20},       charged weight: {3,8,15}.
```

The code checks polynomial idempotence, annihilation by `L`, weighted
self-adjointness, commutation of the two projectors, and commutation with
`K`. Their product is the orthogonal projector onto the joint kernel.
Its rank is its exact integer trace. A column choice is accepted only
when an integer minor modulo 1009 proves that same rank. Thus the selected
columns span the *entire* projector range. The resulting dimensions are
32 for the ground and 228 for the local charged block.

The exact local CAR image of every retained ground-basis vector is checked
to belong to that charged block. Spin rotations give the other spin;
spatial transitivity gives every site. The exact signed bit-complement
particle-hole unitary intertwines the full balanced `N=7` and `N=9`
source matrices and all eight local CAR maps. It gives equal positive
addition and removal energies and weights. These reductions do not infer
Lorentz spin from the hidden total-spin label.

## Exact verification of the spectral proposal

`proposal.py` uses floating eigensolvers only to propose a complete finite
basis. The retained [witness.npz](witness.npz) stores its integer column
indices, vector coefficients with denominator `2^60`, and eigenvalue
proposals with denominator `2^56`. Normal cold checks read these fixed
integers and perform no eigensolver step.

Write the exact integer orbit basis of a reduced block as `X`, and its
source coefficient matrix as `T/16`; the code proves `K Q=Q T/16`.
For all proposed vector columns together, put

```
Y = Q X V / 2^60,       Lambda=diag(lambda_i),
R = K Y-Y Lambda,       delta=||Y^*Y-I||_infinity,
rho^2 = sum_i ||R_i||^2.
```

The Gram matrices, residual numerators and squared norms are evaluated
with checked integers and arbitrary-precision integer products. Square
roots receive outward rational enclosures using integer square roots.
Since `delta<1`, `Y` is an invertible map onto the complete reduced space.
The exact similarity satisfies

```
Y^(-1) K Y = Lambda + Y^(-1)R,
||Y^(-1)R|| <= rho/sqrt(1-delta) = epsilon.
```

The resolvent of the diagonal matrix proves that every eigenvalue lies in
the union of radius-`epsilon` discs. The same assertion holds along
`Lambda+t Y^(-1)R`, for `0<=t<=1`. Each isolated connected component of
those discs retains its initial eigenvalue count by root continuity (or
a finite contour projection). This is not an ordered Weyl inequality for
a nonunitary similarity. All 228 charged discs are disjoint, so each
contains exactly one simple reduced eigenvalue. The isolated lowest
ground disc also has rank one. Four auxiliary excited ground-block
vectors form an overlapping cluster near `-2`; no individual eigenvector
or ordered matching is claimed inside that cluster.

For an isolated column let `q_i=Y_i/||Y_i||` and let `r_i` be its exact
normalized residual bound. All other eigenvalues are at distance at least

```
gamma_i = min_(j != i)|lambda_i-lambda_j|-epsilon
```

from its proposed eigenvalue. The norm outside its spectral line is at
most `r_i/gamma_i`, so a phase of the unit eigenvector can be chosen with
error `eta_i<=2 r_i/gamma_i`. Its energy is within `r_i` of the proposal.
For the unit CAR map, the normalized overlap error is at most
`eta_ground+eta_i`. Exact rational overlap and normalization bounds then
give lower and upper *squared* amplitude bounds. Every lower bound is
strictly positive, including the weakest line. Completeness of the source
image gives unresolved local weight exactly zero.

## Spatial residues and states invisible to linear CAR probes

Let `psi_i` be the normalized origin-stabilizer-fixed charged eigenvector
and `w_i=|<psi_i,c_0down G>|^2`. The signed translations satisfy
`R_x^*R_y=kappa_xy R_(x xor y)` with checked signs `kappa_xy`. Therefore

```
B_i[x,y] = w_i kappa_xy <psi_i,R_(x xor y) psi_i>,
A_i      = eta B_i eta,       eta_x=(-1)^(x1+x2+x3),
G(z)     = sum_i [ A_i/(z-Delta_i) + B_i/(z+Delta_i) ] tensor I_spin.
```

Each translation expectation has error at most `2 eta_i`; rational
interval multiplication encloses every entry of `B_i`. The addition
matrix follows by the exact particle-hole sign rule. The known local
frame rotations return these matrices to the source Weyl frame; they
do not change the poles, local weights or residue ranks.

The orbit of `psi_i` under the finite signed group is irreducible. Any
nonzero irreducible component would contribute a nonzero stabilizer-fixed
vector at the same energy; simplicity of the complete fixed-space line
permits exactly one such component. The eight projective cosets give a
tight frame on this cyclic irreducible space. Central signs disappear
from its rank-one projectors. If its dimension is `d_i`, then

```
Tr(Gram_i^2)=8 sum_s |<psi_i,R_s psi_i>|^2=64/d_i.
```

The normalized approximate vector's eight squared moments have total
error at most `32 eta_i`. Their interval isolates one integer `d_i` among
1 through 8. Every line gives `d_i=4`. Consequently the all-site,
fixed-spin-weight cyclic dimension is 912, and all linear CAR probes
reach 1,824 states in each of the full `N=7` and `N=9` sectors. Each full
sector has 11,440 states; **9,616 states per sector are invisible to these
probes from this reference**. The signed two-spin response has cyclic
dimension 3,648. Other dark states can share a visible energy; rank four
is the visible spatial residue rank, not a claim about the full many-body
eigenspace multiplicity. The complete dark-sector spectrum and a global
neutral-gap estimate are not provided.

The exact local sum rules are

```
sum_i w_i=1/2,       M0=1,       M1=0,
M2=2 sum_i w_i Delta_i^2=3/256+1/4=67/256.
```

The rational energy/weight intervals enclose all these values. The
conditional mean energy is `2 sum_i w_i Delta_i`, using the exact branch
weight, rather than renormalizing an incomplete list.

## Replay and scope

```
python -W error -m pytest -q model/centered_particle_spectrum
python -W error -m model.centered_particle_spectrum --write
python -W error -m model.centered_particle_spectrum --check
```

A fresh candidate can be generated separately with `--propose`; it must
pass the same exact checks and is not substituted during a normal check.
No numerical solver success flag, threshold-zero weight or inherited
receipt `all_pass` flag replaces an exact check in this packet. This
construction establishes the finite source/reference response. A native
preparation and measurement of this changed ground, a justified
interacting relativistic refinement, dimensional calibration, physical
species correspondence and the shared gravity premises require their own
derivations.
