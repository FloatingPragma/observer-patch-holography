# Compact opposite-degree textures and isolated fermion bands

## Chosen finite-energy geometry

Use the same two-spin nearest-neighbor hopping as the radial construction,
with lattice spacing and hopping coefficient chosen to be one. Fix an integer
`a>=1` and put `R=a+1`. On lattice sites define

```
f(x,y,z)=(x,y,z^2-a^2)                  for max(|x|,|y|,|z|)<=R,
f(x,y,z)=(0,0,1)                      otherwise.
n=f/|f| when f!=0; n=0 at c_+=(0,0,a), c_-=(0,0,-a).
h=T+M n.sigma,                        M>=8.
```

This texture and its amplitude zero at the two cores are declared design.
They are not selected by bare repair. A continuous compact extension uses
`t(r)=min(1,max(0,R+1-r))` and
`f=(t x,t y,t(z^2-a^2)+1-t)`. Its only zeros are `c_+` and `c_-`:
for `t>0`, a zero requires `x=y=0`; in the transition region `|z|>=R>a`
the third component is positive. At the two zeros the Jacobian is
`diag(1,1,+2a)` or `diag(1,1,-2a)`, respectively.

The actual unit cubic shells around the cores have exact degrees `+1` and
`-1`. Every pair of vertices on each interpolation triangle has nonnegative
scalar product, so the normalized barycentric interpolation is defined.
The code counts signed rational regular-value preimages on every oriented
triangle using the complete shell routine. The outer shell of radius `a+2`
is constant and has degree zero. The local shell edges are disjoint for
`a>=2`; for `a=1` their shared face is retained and its energy is not counted
twice. The integer charge classifies these interpolated shell maps. It is
not a superselection operator on the fermionic Fock space.

For the ordered-spin source, fill each core direction by north and rotate the
selected ordered ground so that all its directions agree with this assignment.
The reference direction may be globally rotated from `e1` to north because the
source Hamiltonian is rotation invariant. At most `s^3` sites differ from the
uniform reference, where `s=2a+3`. Exactly `3s^2(s+1)` nearest-neighbor bonds
touch that cube. Each rotated-ground bond has excess in `[0,J]`, hence

```
0 <= Delta E_spin <= 3s^2(s+1) J.
```

The lower shell bound is `pi J kappa` at `a=1`, and `2pi J kappa` for disjoint
unit shells at `a>=2`. The spin-core fill and the zero fermion mass at the
cores are distinct declared parts of the design. No equality between an
ordered-spin expectation and the fermion mass amplitude is assumed.

There is also a uniform bound on the fermionic energy change relative to the
constant north texture. At a noncore site the traceless one-body difference
has eigenvalues `+/- M|n-north|`. Its second quantization has operator norm
`M|n-north|<=2M`. At each core it has norm `M`. Thus on every finite Fock space,

```
||dGamma(h_texture)-dGamma(h_uniform)|| <= 2M(s^3-1),
|E_vac(texture)-E_vac(uniform)| <= 2M(s^3-1).
```

These bounds are independent of the exterior cube radius. They remove the
infinite-tail energy problem for this compact pair. They do not establish
its formation law, equilibrium, positive binding energy, lifetime or the sum
of two independent constituent rest masses. Fermion occupation changes no
texture winding.

## Four isolated modes, with no assumed zero-mode index

Let `P` project onto the four spin coordinates at the two cores and `Q=1-P`.
The cores are not adjacent, so `PhP=0`. The commuting open-cube centered
differences give `||T||<=sqrt(3)<2`. The onsite potential has four zeros and
all other eigenvalues `+M` or `-M`. Weyl bounds therefore give exactly four
finite-volume eigenvalues in `[-2,2]`, with every complementary eigenvalue
outside `(-M+2,M-2)`.

Write `h=[[0,B*],[B,D]]`. On `Q`, the onsite potential is invertible and
`D` has no spectrum in `(-M+2,M-2)`. Direct nearest-neighbor multiplication
gives

```
||B||^2=b^2=7/4 (a=1),                b^2=3/2 (a>=2).
```

For `a=1`, the complete eigenvalues of `B*B` are `5/4,5/4,7/4,7/4` because
the two cores share one neighbor. At larger separation `B*B=(3/2)I_4`.
Every low eigenvalue obeys its Feshbach equation

```
E u = -B*(D-E)^(-1)B u.
```

Since `|E|<=2`, define `r=b^2/(M-4)` and obtain

```
|E|<=r,
||Q psi||^2/||P psi||^2 <= b^2/(M-2-r)^2,
||P psi||^2 >= 1/[1+b^2/(M-2-r)^2]       for ||psi||=1.
```

All bounds are uniform in exterior volume. They also apply on the infinite
lattice. There the inverse of `D-E` exists throughout the gap, and the four
eigenvalues of `-B*(D-E)^(-1)B-E I_4` are strictly decreasing functions of
`E`: the matrix derivative is `-I_4-B*(D-E)^(-2)B<=-I_4`. At `E=-2` it is
positive definite and at `E=2` negative definite. Each eigenvalue crosses zero
exactly once. Feshbach inversion excludes other gap spectrum. Thus the
infinite lattice has exactly four isolated modes in this interval, counted
with multiplicity; this is not a continuum-limit statement.

The leading zero-energy Schur term is also exact. With `V` the exterior
onsite potential, all single-core diagonal terms in `-B*V^(-1)B` cancel
between opposite neighbors. At `a=1` the common neighbor gives off-diagonal
blocks `-sigma3/(4M)`, so the leading four eigenvalues are `+/-1/(4M)`, each
twice. At `a>=2` no neighbor is shared and the whole leading term vanishes.
This explains the first separation-dependent cancellation without fitting a
spectral trend. Higher resolvent terms remain in the complete spectrum.
The tiny numerical splittings at `a=2` have no separate zero-exclusion
certificate; the exact four-dimensional band and localization bounds apply.

## An explicit exterior-volume error

Let the finite exterior radius be `L>=R+1` and compare it with any larger
cube or the infinite lattice. For real `|E|<=r`, expand the exterior resolvent
around its onsite part `V`:

```
(D-E)^(-1)=sum_{n>=0}[-(V-E)^(-1)T_Q]^n (V-E)^(-1).
```

Its ratio is bounded by `q=2/(M-r)<1`. A path starting at a neighbor of a
core lies initially inside radius `R`. For finite and infinite compressed
resolvents to differ, it must reach outside radius `L` and return. This
requires at least `2d` hopping factors, `d=L+1-R`. Hence the two four-row
self-energies differ in norm by at most

```
delta_L = [2b^2/(M-r-2)] [2/(M-r)]^(2d).
```

The monotone eigenvalue crossings have slope at most `-1`; Weyl's eigenvalue
bound then gives `|E_j(L)-E_j(infinity)|<=delta_L` for all four ordered
energies. The estimate concerns the same potential on nested lattices and
includes all exterior modes.

For `a=1,M=64,L=3`, an exact dyadic complete-basis certificate encloses all
686 eigenvalues of the literal normalized source matrix. It includes outward
rational square-root errors, exact integer Gram and residual calculations,
and a polar-decomposition/Weyl exhaustion bound. No floating residual is
used as proof. All four low intervals exclude zero. Subtracting `delta_3`
from their minimum absolute endpoint also excludes zero on the infinite
lattice. The numerical low energies are approximately
`+/-0.003905899328` and `+/-0.003906437151`; each mode has core probability
above `0.99957`. The filled-sea energy difference from the uniform source is
approximately `127.9700` at this volume. This energy includes the changed
background and differs from the isolated mode gap.

## Native currents and pairing on the same band

For any complete four-column band isometry `W`, set `A=W* h W` and let
`C=D_stagger i sigma2`. The literal hopping and every real Pauli texture obey
`hC=-C conjugate(h)`, with `C` antisymmetric and unitary. The symmetric
spectral window is invariant under this particle-hole map. Thus

```
hW=WA,  C conjugate(W)=WS,  S=W* C conjugate(W),
S^T=-S,  S*S=I,  A S=-S conjugate(A).
```

Extend by the actual native fifteen-state current representation and three
declared family copies. Every spatially uniform native current intertwines
through `W tensor I_45`. For any symmetric internal kernel `K`, choose the
same local pairing `Delta=C tensor conjugate(K)` as the radial construction.
The complete Nambu inclusion
`E=diag(W tensor I,conjugate(W) tensor I)` then has core

```
B_core = [[A tensor I, S tensor conjugate(K)],
          [(S tensor conjugate(K))*,-conjugate(A) tensor I]].
```

The full and core spectra consist of the respective energies of `h` and `A`
plus or minus the singular values of `K`. In particular, the nonzero compact
band energies shift the retained lepton coefficients. Dropping `A` gives a
false intertwiner. A complex symmetric-kernel check evaluates every source
row, and all twelve actual native generators satisfy the tensor inclusion.
The 360 core Nambu rows include spin, families and redundant holes; they are
not a count of physical species. Local gauge links and Gauss constraints
are outside this tensor-current statement.

A Pauli-z perturbation at the two cores leaves every surrounding shell datum
unchanged but changes the four energies. Its continuous-texture interpolation
may move a zero within a core cell. Thus these local shell charges do not force
zero spectral energy. Finite-support spin erasure paths and pair annihilation
are compatible with the energy bounds. Their exclusion would require an
additional restriction on admitted dynamics. No autonomous defect-selection
law, relativistic mass shell, interacting particle catalogue or Einstein
response follows from this construction.
