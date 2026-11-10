# Narrow charged-energy windows on the unchanged native source

## Source and exact local input

The source is the complete physical eight-cell Hamiltonian
h=H0+lambda T, lambda=2^-34, with its actual unique neutral ground g and
energy Eg. Every native Fock, scalar, link and Gauss sector remains in
the carrier. The full hopping norm bound is v=360*2^-34. A native e^c
removal is compensated by the exact two-spin one-particle sector of the
retained CAR reference, whose energy is 4. The complete reference Fock
carrier and every other outcome remain as in the charged-response source.

The exact local ground energy e is the least real root of

```
f(e)=[e(e-1)]^3-313[e(e-1)]^2+31708e(e-1)-1041468.
```

The complete visible charged target polynomial, rebuilt from the actual
CAR intertwiners and target matrices, is

```
q(x)=(x-1)(x^4-2x^3-212x^2+213x+11114)
           (x^4-2x^3-200x^2+201x+9858).
```

All nine real roots epsilon_j, their positive weights w_j and the exact
identity sum_j w_j=1/2 are supplied by the complete local resolvent. The
code refines the roots with rational Sturm isolation to width at most
10^-100, identifying their embeddings inside the parent's root intervals.
It retains the exact algebraic local weight bounds from the freshly
replayed resolvent. The root at x=1 is an exact rational singleton.

The physical atomic cube singlet u has each local doublet maximally mixed.
After the compensated contraction B its orthogonal components psi_j have
norm squared w_j and exact atomic source/reference energies
a_j=7e+epsilon_j+4. The frozen full-source bound is
||g-u||<=eta after phase alignment, with eta<0.000169320. Since ||B||<=1,
the actual and atomic pumped states differ by at most eta. No atomic
reference is substituted for g in the final probability.

## A quadratic neutral-ground shift bound

Let P be the complete atomic ground projector and Q=I-P on the physical
native source. The parent proves PTP=0 and
Q(H0-8e)Q>=g_at Q, where g_at>0.04117058. The complete perturbation obeys
||lambda T||<=v<g_at. For a unit vector z with p=||Pz|| and q=||Qz||,

```
<z,(h-8e)z> >= (g_at-v)q^2-2vpq
                >= -v^2/(g_at-v).
```

The last inequality completes a square and uses p<=1. Conversely every
atomic ground vector has expectation 8e because PTP=0. The variational
principle gives the full-source bound

```
8e-delta0 <= Eg <=8e,       delta0=v^2/(g_at-v).
```

This includes every physical source complement. It uses no ground-state
rank comparison and no charged-sector truncation. Its order-v^2 width
replaces the order-v neutral-energy uncertainty when reporting source
excitation intervals. The retained reference cost 4 is subtracted only
when explicitly reporting source-only excitation energy.

## Pythagorean transfer of the selected atomic weight

Let I be a rational interval centered on a certified isolator midpoint
for a_0, with radius r=cv. Let P_I be the **actual full h+H_ref** spectral
projector. Define d_0 as the distance of a_0 to the complement of I, and
d_j for j>0 as the distance of a_j to I, using positive rational lower
bounds. Every source state outside the probe module remains present in
P_I. The exact atomic eigenvectors satisfy

```
||(h+H_ref-a_j)psi_j|| <= v ||psi_j||.
```

The spectral theorem gives ||(I-P_I)psi_0||<=v sqrt(w_0)/d_0. Orthogonality
of P_I and I-P_I then yields

```
||P_I psi_0|| >= sqrt(w_0[1-(v/d_0)^2]).
```

This controls the selected norm directly. It is stronger than subtracting
v sqrt(w_0)/d_0 from sqrt(w_0). For every other visible atom,
||P_I psi_j||<=v sqrt(w_j)/d_j. Therefore, with b_j the atomic upper
weights and a_0 the atomic lower weight, define

```
ell_other=sum_(j>0) v sqrt(b_j)/d_j,
A_lower=sqrt(a_0[1-(v/d_0)^2])-ell_other-eta,
A_upper=sqrt(b_0)+ell_other+eta.
```

For the strictly positive A_lower in every recorded setting,

```
A_lower^2 <= ||P_I B g||^2 <= min(1/2,A_upper^2).
```

The code uses exact rational arithmetic and outward 128-bit dyadic square
roots. It rounds ell_other upward before forming the probability bounds.
It does not omit dark eigenstates, assert that I has rank two, or use an
eigensolver for the full charged source. Only the complete **visible atomic
support** is finite and explicitly enumerated in the residual estimate.

Every refined I is strictly contained in the original radius-1/1024 bin.
The actual per-spin Gram matrix is scalar under the proved finite native
cubic G1 symmetry, so the family-pump weight also equals the per-spin
unnormalized residue. Exact antiunitary symmetry of the actual source
gives total CAR success 1/2. The complete true-ground readout law is

```
pump failure: 1/2,
refined energy bin: w_sharp,
complete energy complement: 1/2-w_sharp.
```

All three weights are positive for the recorded windows. The energy
projector, complement and reference occupation commute with the actual
clock. Original physical two-copy preparation can be retained before
this instrument: its success q is followed by the declared P0 filter
onto both actual native grounds, with conditional success p. The complete
grouped law is then

```
1-q, q(1-p), qp/2, qp w_sharp, qp(1/2-w_sharp).
```

Both source copies, all spin references, filter counters, static charged
reference and outcome pointers remain present. P0 is a separate exact
spectral apparatus service with its own failure and cost. This supplies
no efficient native gate implementation. The source clock is the full
sum of the two unchanged h copies and the declared reference Hamiltonian;
no failed branch silently stops that clock.

## Resolution and retained residue

For each radius multiplier c in {2,4,8,16,32,64,128,256,512}, the receipt
gives the full source/reference bin, the source excitation interval above
the actual Eg, its width and its unnormalized weight. If [e_lo,e_hi]
isolates e and I=[I_lo,I_hi], the latter energy interval is exactly

```
[I_lo-4-8e_hi, I_hi-4-8e_lo+delta0].
```

At c=128 the lower residue exceeds 0.419706, hence exceeds the original
certified lower residue 0.4197041521592. Its energy width is about
5.36441804*10^-6, more than 364 times narrower than the original source
bin. At c=2 the lower residue exceeds 0.31477 and the source energy width
is about 8.38190424*10^-8. These are source-derived finite energies with
the retained reference cost made explicit, without a GeV calibration.

## Spectral selection and its exact limitation

Within a chosen refined bin, select E_star_sharp by maximum per-spin
residue, with a lower-energy tie break. Since the bin has positive weight
and the full source dimension is D, its selected residue is at least
w_sharp,lower/D. Its normalized origin CAR vectors form a G1 doublet
and are actual eigenvectors of the unchanged h. The exact native
antiunitary supplies the opposite source-charge partner with equal energy.
Equal complete reference costs for both partners require the enlarged
number-preserving dual reference of the native join; particle-hole on the
old 4N reference would change its energy and is not used.

This is an explicit spectral selection inside a nested bin. The original
E_star maximizes residue over the entire old bin. Neither a positive nested
weight nor a small outside weight proves that this original maximizer is
inside the narrower window when the inside can contain many eigenvalues.

An exact norm-perturbation counterexample demonstrates the logical issue.
Take v=1, H0=diag(0,...,0,3) on dimension17 and initial atomic vector e_1.
Let R rotate e_1 and the final basis vector with sine1/4 and cosine
sqrt(15)/4. Inside the first16 coordinates, conjugate the distinct energies
`(2j-15)/150`, j=0,...,15, by the normalized order16 Hadamard matrix.
Call that block A and set H=R diag(A,3) R^T. Then

```
||H-H0|| <= 3/4+1/10=17/20 <v,
```

because the rotation changes the rank-one high-energy projector by norm
1/4 and ||A||=1/10. The actual e_1 spectral weights are 15/256 at each
of the sixteen low energies and 1/16 at energy3. A radius2 window contains
15/16 of the total weight, but the largest individual residue is the
outside atom at3. Both lie inside a radius4 window. Rescaling gives the
same obstruction at arbitrary perturbation norm. This is a counterexample
to inference from the bounds, not a claim that this spectrum occurs in
the native source.

The obstruction also applies when almost all weight is retained. With
high energy 256v, rotation sine1/512, a normalized Hadamard matrix on
N=2^19 low states and low energies of norm at most v/10, the perturbation
norm is at most 3v/5. A radius128v window contains weight
1-1/262144, while the outside atom has weight1/262144 and each inside
atom has smaller weight (1-1/262144)/N. The outside energy256v is still
inside the original radius-1/1024 bin at the actual v used here. This
construction is specified by exact finite matrices without storing them.
It identifies a needed full-source spectral-rank or separation premise,
which the local visible resolvent and norm bounds alone do not supply.

Consequently an affine or free-Dirac construction using E_star_sharp must
declare that selection and its corresponding matrix-unit couplers as a
separate source choice. The frozen coarse selection and every existing
coupling remain unchanged. The finite G1 doublet alone identifies no
continuous physical spin, Standard Model species or relativistic mass.
