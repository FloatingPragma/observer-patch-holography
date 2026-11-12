# Complete finite-frequency three-point function

Use the literal 320-state lepton/Higgs Hamiltonian, six fermion fields and
two scalar fields of `model/yukawa_defect_vertex`. The primitive interaction
and scalar cost are one. Its ground projector P has rank two, energy -2
and gauge-singlet range. Prepare I_320/320 as the marginal of the canonical
maximally entangled source/dual-reference state. The complete source energy
measurement gives rho=P/2 on its ground record, with original success
probability 1/160; all other energies and the quantum reference are retained. The frequency
pairs (1,2), (2,3), (4,8) are fixed in source-time units without particle
data. Matrices have dimension at most 320; no spatial or family extension
is included.

## Full source resolvent

Write H in blank/occupied scalar blocks as [[0,W*],[W,I]]. The exact
eigenprojectors P_r of W*W have r=0,...,6. Let A_r be P_r embedded in the
blank block, X_r the off-diagonal W P_r plus its adjoint, and Z_r the
occupied W P_r W*. The occupied dark projector is
Q_dark=Q_occupied-sum_(r>0) Z_r/r, including all 198 dark dimensions.
For z outside the spectrum,

```
(z-H)^-1 = A_0/z + Q_dark/(z-1)
          + sum_(r=1)^6 [(z-1)A_r+X_r+z Z_r/r]/[z(z-1)-r].
```

The executable checks its full matrix inverse identity exactly at every
consumed imaginary frequency, with D=H+2. No cyclic compression or
intermediate charged state is omitted. A separate dense inverse checks
the first frequency numerically.

## Time orderings and amputation

For O=(c_i,c_j,phi_k*) and frequencies (w1,w2,-w1-w2), define the
Euclidean three-point transform by ordering the two positive time
intervals. For each permutation (a,b,c), its contribution is

```
sgn_F(a,b,c) Tr[rho O_a (i w_a-D)^-1 O_b (-i w_c-D)^-1 O_c].
```

The sign is negative precisely when the two fermions exchange order.
Indeed integrating exp[-(D-i w_a)u] and exp[-(D+i w_c)v] over positive
u,v gives the two resolvents and two cancelling minus signs. The charged
end operators have zero matrix elements within the singlet ground range,
so these integrals converge; the full matrices retain that range anyway.
Every one-point function vanishes exactly, hence all disconnected
three-point products vanish. This definition fixes the overall sign of
the amputated coefficient; no physical scattering convention is assumed.

The full propagators are

```
G_ij(z)=Tr[rho(c_i R(z)c_j* - c_j* R(-z)c_i)],
D_ab(z)=Tr[rho(phi_a R(z)phi_b* + phi_b* R(-z)phi_a)],
R(z)=(z-D)^-1.
```

These equal the complete addition/removal spectral formulas. In the
verified native frame G is diagonal and D is scalar on the doublet.
Each of the fourteen scalar spectral rows has equal addition and removal
matrices, so D(z)=D(-z). The executable checks this evenness at the actual
external frequencies, justifying the conjugate scalar leg convention.
The amputated tensor is T_ij,k/[G_ii(i w1) G_jj(i w2) D_kk(-i(w1+w2))].
All denominators are nonzero at the fixed frequencies. The scalar moment
normalization is Z_H=14/5, so replacing phi by phi/sqrt(Z_H) multiplies
this amputated tensor by sqrt(Z_H). Its dependence on the external
frequencies is calculated rather than replaced by the contact coefficient.
The calculation exhausts all 72 components at (1,2), with eight nonzero
entries, and checks exchanged-frequency fermion antisymmetry.

For example the (0,4,0) amputated component is 474925/2394532 at (1,2)
and 11640331925/26140075569 at (2,3). The source has no selected physical
mass shell, spacetime momentum or renormalized Standard Model coupling.
The finite amputated response is not identified with an observed Yukawa.

## A complete invariant observation

Sum the six ordered operator products over the four native oriented
pair/scalar components. Each field has operator norm one, verified by
its partial-isometry identity, and ||R(i w)||<=1/|w|. At (1,2), the sum
of these product bounds is eight. Thus their sum divided by eight is a
contraction B. Exact full-matrix commutators show [B,Q_a]=0 for all twelve
native currents. The four effects

```
(I+Re B)/4, (I-Re B)/4, (I+Im B)/4, (I-Im B)/4
```

are positive and sum to I. Their square-root instrument is a declared
finite apparatus preserving the complete source and its untouched
reference; it is not claimed to be a compiled sequence of native pulses.
The executable constructs all four square roots on the full 320-state
source, checks completeness, and independently matches its expectation
to the exact ordered three-point traces. All thirteen failed preparation
energies and four successful pointer outcomes give seventeen original
records with total probability one. No successful branch is renormalized.
