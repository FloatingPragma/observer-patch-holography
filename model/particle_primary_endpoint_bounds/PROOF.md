# A local electric bound at the interacting primary endpoint

## Actual source and conclusion

Use the [primary member](../particle_common_contract/PROOF.md) with its
full sixty-mode Wilson matter, all native link words, physical Gauss
carrier, selected polynomial state, retained reference and unchanged
source clock. Write Hhat=aH and E00=-90V as in the
[spectral-scale proof](../particle_primary_spectral_scale/PROOF.md).
Every bound below is evaluated at the actual hopping coefficient one,
electric coefficient 105 and magnetic coefficient 1/7. No second-order
energy coefficient is extrapolated to that endpoint.

The result is an exact rational lower bound on the full physical ground
energy. It is obtained before the Gauss restriction, so it does not rely
on a factorization of the physical carrier. A stronger variational upper
bound can be combined with it for the same member. This ground enclosure
is not an isolated excitation branch or particle rest energy. The
retained-penalty calculation below exposes the separate branch-isolation
inequality; it does not declare that inequality satisfied.

## A full two-site Fock bound with one quantum link

Let D_x=dGamma(beta_D)_x+30I be the nonnegative matter defect count.
For an interior bond xy, assign onsite coefficient m>0 at each endpoint
and electric coefficient e>0 to its link. Its shifted cluster is

```
K_xy=m(D_x+D_y)+(e/2) C_link+T_xy.                    (1)
```

The literal one-body hopping block is
`S_mu=(-beta_D+beta_D gamma_mu)/2`. It is a rank-two partial isometry.
Pointwise in the link coordinate, a native unitary change of basis at
one endpoint removes the matrix R(U) from the fermion matrix. This is
only a pointwise spectral calculation, not a conjugation that discards
the quantum electric operator. For each internal column the actual
eight-dimensional two-site Dirac characteristic polynomial is

```
[z^4-(2m²+1)z²+m^4]^2.
```

Its eigenvalues are `±(sqrt(m²+1/4)+1/2)` and
`±(sqrt(m²+1/4)-1/2)`, each twice. The full Fock minimum is the sum
of all negative one-body energies. Including all fifteen internal
columns and the shifted constant 60m gives the operator inequality

```
m(D_x+D_y)+T_xy >= -delta(m) I,
delta(m)=60(sqrt(m²+1/4)-m).                          (2)
```

The product cutoff is an orthogonal compression of the complete parent
term, so (2) survives it. This is a full-Fock bound, not a one-particle
Hamiltonian used in place of the many-body source.

Let P be the constant function projection on this link, with identity
on both Fock factors, and Q=I-P. The all-word Casimir proof gives
`C_link>=3Q/10`. The hopping has zero Haar mean, so P T_xy P=0.
The onsite defect term commutes with P. Thus

```
P K_xy P=m(D_x+D_y)P>=0,
Q K_xy Q >= Q[(e/2)C_link-delta(m)]Q.                 (3)
```

When `3e/20>delta(m)`, the Q block is positive. Resolvent order and
completion of the square in the P,Q blocks reduce a lower bound to

```
F=P T_xy Q[Q K_xy Q]^(-1)Q T_xy P
 <=P T_xy Q[(e/2)C_link-delta(m)]^(-1)Q T_xy P.        (4)
```

The inverse in (4) acts only on the actual single-entry link vectors
created from the constant function. They lie inside every selected
cutoff. Different native irreps and their conjugates are Haar
orthogonal. The denominator in native block r is consequently the
exact `e C_r/2-delta(m)`.

## The Haar contraction retains all occupations

Use the two active Dirac rows of the bond partial isometry as flavors
s=1,2. For an irrep of dimension d, put

```
A_ij=sum_s a_si* b_sj,
T_r=sum_ij R_r(U)_ij A_ij + adjoint.
```

All five native charges are nonzero and no pair of them is conjugate.
Schur orthogonality eliminates two forward coefficients and gives

```
P T_r² P=(1/d) sum_ij(A_ij A_ij*+A_ij* A_ij)
         =N_A+N_B-N_A N_B/d-(4/d) S_A.S_B.            (5)
```

Here N_A,N_B count the active flavors and S_A,S_B are their ordinary
U(2) flavor-spin generators, summed over the d internal columns. They
are auxiliary finite algebra labels, not physical Lorentz spin.
The inactive Dirac rows remain in the full Fock carrier.

At fixed occupations x,y in [0,2d], the spins obey
`j_x<=min(x,2d-x)/2` and likewise for y. Coupling the two spins gives
`S_A.S_B>=-min(j_x,j_y)(max(j_x,j_y)+1)`. This yields the sharp bound

```
0<=P T_r² P<=2(d+1) P.                               (6)
```

For completeness, write x=d+sA,y=d+tB with A>=B>=0, A,B<=d and
s,t in {-1,1}. Substituting the maximal spins bounds (5) by
`2d+2-A-B-2A/d+(1-st)AB/d`. For st=1 it is at most 2d+2; for st=-1,
use `2AB/d<=2B<=A+B`. Equality is attained at x=y=d and
j_x=j_y=d/2 in their coupled singlet. The exact finite CAR controls
check (5) and (6), including occupations outside the product vacuum.

Combining (4)-(6), completing the block square, and then using P<=I
gives the all-input one-link bound

```
K_xy>=-c(m,e) I,
c(m,e)=sum_r 2(d_r+1)/(e C_r/2-delta(m)),
provided 3e/20>delta(m).                              (7)
```

The rational calculator uses an outward upper bound on the square root
in delta(m). Its denominators are lower bounds and its c is an upper
bound, so the displayed ground floor is directed outwards. It rejects
the Schur step when the Q-block margin is not positive.

## Sum the clusters on the unchanged open source

Let d_max=3 for the side-two cube and d_max=6 for every larger selected
box. Choose r in [0,1), t=1-r, m=3t/d_max and e=105t. Retain rH0
explicitly. Each interior bond receives the cluster (1). A cell of
degree d_x receives md_x<=3t of its onsite defect cost; the remainder
is nonnegative. Each interior link is assigned to exactly one bond.
The unassigned electric terms, including the exterior links, are
nonnegative and are not traced out or deleted from the source.

The magnetic term has the inherited operator floor -P_plaq/7. With
the actual counts B=3L²(L-1), P_plaq=3L(L-1)², equation(7) proves

```
Hhat >= E00-B c(m,e)-P_plaq/7 + r(H0-E00).            (8)
```

This is a decomposition used to prove an inequality for Hhat(1).
The allocated m,e are not changed source parameters. All cross-bond
fermion and gauge interactions remain in (8); no independent-bond
state or invariant tensor product is assumed. Restricting the ambient
operator inequality to the actual Gauss carrier preserves it, with
arbitrary untouched references.

At r=0, the smallest primary box uses m=1 and has
`delta=60(sqrt(5)/2-1)<63/4`. Larger boxes use m=1/2 and also satisfy
`delta=30(sqrt(2)-1)<63/4`. The resulting lower bound is therefore
admitted at the actual coefficients for both size domains, uniformly
in the product word cutoff. Dividing (8) by the actual spacing a gives
the code's `physical_ground_lower` in the unchanged source energy units.
The source polynomial state and its canonical reference are not used as
trial states or changed by this spectral inequality.

## Retaining a defect penalty and testing branch isolation

The decoupled physical H0 has a unique ground and exact gap 6. The
min-max principle applied to (8) gives the ordered full physical
second eigenvalue, counting multiplicity,

```
E_1(H) >= [E00-B c(3t/d_max,105t)-P_plaq/7+6r]/a.    (9)
```

A proved physical ground upper bound U0 would isolate the ground
branch only if the right side of (9) exceeds U0. The program evaluates
the positivity margin and (9) at the stated r values; it does not turn
a positive one-link Q-block margin into a global spectral gap.

The function delta(m) decreases with m, so c(3t/d_max,105t) decreases
with t. Every admitted r therefore has c>=c(3/d_max,105). Consequently
the available second-eigenvalue lower bounds have the method ceiling

```
[E00-B c(3/d_max,105)-P_plaq/7+6]/a.                 (10)
```

For this ceiling the implementation uses a lower enclosure of the exact
c, obtained from a lower square-root enclosure. This reverses the
rounding direction used for a ground floor. Monotonicity is applied to
the exact function, not to a rounded square-root staircase.

This ceiling is an upper bound on the numbers supplied by this method,
not an upper bound on the actual second eigenvalue. If it lies below
U0, no choice of r in this retained-penalty family establishes branch
isolation with that U0. A tighter estimate must retain the defect-count
dependence of the effective quartic operator (5), control larger local
clusters, or supply another independent spectral estimate. The global
maximum 2(d+1) discards that occupation dependence even though it is
sharp on its stated all-input domain.

The source-authenticated join with the
[eleven-packet endpoint upper bound](../particle_primary_variational_readout/PROOF.md)
evaluates this test on the smallest selected member. With V=8, B=12,
P_plaq=6 and a=1, the local loss cap is approximately 1.61432705912 per
bond. The exact rational producers imply the outward enclosure

```
-740.229068 <= E0(H) <= -723.767819.
```

The unrounded certified interval has width approximately 16.46124762383.
The ceiling (10) minus the actual rational variational upper is less
than -10.46. Thus every admitted retained fraction r in this method
fails to isolate the ground branch using that variational upper. This
is a quantified failure of this bound, not evidence against an actual
gap. The larger-box m=1/2 allocation also passes the local Schur test,
with loss cap approximately 2.66658836812 per bond; it supplies the
same kind of full-source floor without a branch claim. Decimal values
in this paragraph are displays of the rational receipt, not inputs.

The result provides genuine endpoint ground information with an exact
source fingerprint. It supplies neither original probe residue nor a
momentum-resolved particle branch. Those conditions, together with the
chosen finite approximation errors, determine whether a finite-member
rest-energy comparison is admitted. No physical mass or GeV number is
inserted into the construction.
