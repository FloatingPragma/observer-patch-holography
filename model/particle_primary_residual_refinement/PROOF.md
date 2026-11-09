# Exact local moments for the physical Ritz residual

## Source and trial

Use the unchanged [primary source](../particle_common_contract/PROOF.md),
with aH=H0+V at interaction coefficient one. The
[eleven-vector trial](../particle_primary_variational_readout/PROOF.md)
has the same exact rational coefficients t_h,r and t_l,r throughout.
Its auxiliary anchor Omega fills the two negative-beta Dirac rows and
has constant links. The selected polynomial density and retained idle
reference are not replaced by this trial. Only the estimate on the
trial's full omitted-space residual is changed.

Work first in the parent Haar coordinate algebra before the product
word cutoff. Write V=J+M and V_t=J_t+M_t. J is all five hopping blocks;
J_t weights them by t_h,r. M is the real character multiplication
operator, weighted by t_l,r in M_t. These parent operators commute in
the combinations needed below: the native hopping blocks have disjoint
CAR modes and magnetic multiplication is scalar on their Fock factors.
The compressed finite operators need not commute.

## Full-Fock fourth moment through particle-hole sectors

At fixed link coordinates, let the one-body hopping matrix in its empty
and occupied Dirac subspaces be

```
h=[A B; B* D],       J=dGamma(h),       Tr(P_- h)=0.
```

Its action on Omega creates one particle-hole pair, with coefficient B.
Applying J again produces orthogonal zero-, one- and two-pair parts.
Put z=Tr(B*B), K=B*B and C=AB-BD. The three squared norms are
z², ||C||_HS² and 2(z²-Tr K²). Consequently

```
<J⁴>=3z²-2 Tr K²+||C||_HS².                         (1)
```

This is the full Fock action, including the two-pair sector. It does
not replace J by its one-particle restriction. For the native-weighted
J_t, the off-diagonal block in native block r is t_h,r B_r. A singular
value decomposition within each block diagonalizes both pair-creation
coefficient matrices. The coefficient of a two-pair state i<j in
J J_t Omega is b_i b_j(t_i+t_j). Its squared norm is
z z_t+q_h²-2 sum_r t_h,r² Tr K_r², where

```
z=15 B_count,  z_t=B_count sum_r d_r t_h,r²,
q_h=B_count sum_r d_r t_h,r.
```

The zero-pair norm is q_h² and the one-pair norm is
sum_r t_h,r² ||C_r||_HS². Thus

```
<J² J_t²>=z z_t+2q_h²
           +sum_r t_h,r² (||C_r||_HS²-2 Tr K_r²).     (2)
```

The literal finite CAR controls execute these sectors independently of
the identity. Shared native link coordinates do not create mixed Fock
covariances between different native blocks: their anchor factors are
independent and their second moments are constant in those coordinates.

## Exact two-step path counts

Let d_x be the interior degree, B_count the number of bonds,
S2=sum_x d_x², N2=S2-2B_count the number of directed nonbacktracking
two-step paths, and T2 the number of directed straight two-step paths.
The actual Wilson spin blocks are obtained from
S_mu=(-beta+beta gamma_mu)/2. Their occupied and empty diagonal blocks
are opposite adjacency matrices divided by two. For K_r=B_r*B_r,
backtracking gives (d_x/4)I on each occupied two-spinor. Every other
two-step path has squared spin matrix norm1/8. Distinct paths between
the same endpoints contain a link occurring only once in their cross
product. Haar integration kills that cross term. Hence

```
integral Tr K_r² = d_r (S2-B_count)/4.                (3)
```

For C_r, backtracks cancel. A straight path has squared spin norm1/2;
a right-angle path has norm1/4, from the Clifford trace relations.
The same Haar argument removes distinct-path cross terms. Therefore

```
integral ||C_r||_HS² = d_r (N2+T2)/4.
G(L)=(S2-T2)/4,
integral <J²J_t²>=z z_t+2q_h²-G(L) sum_r d_r t_h,r². (4)
```

For the open cube of side L>=2,

```
S2=36L³-66L²+24L,    T2=6L²(L-2),
G(L)=3L(L-1)(5L-4)/2.
```

At L=2, the per-native-dimension values in (3), the C norm and G are
15,12,18. Counts use only actual interior bonds. Ghost links introduce
no hopping paths. The argument uses unitarity and zero Haar mean of
each nontrivial native representation; it does not replace nonabelian
links by phases in the primary source. Small phase controls test the
path combinatorics separately.

## Exact magnetic and mixed moments

Let x_r=chi_r+conjugate chi_r and S=sum_r x_r. On one plaquette,
f=-S/210 and g=-sum_r t_l,r x_r/210. Orthogonality gives

```
u=<f²>=10/210², v=<g²>=2 sum_r t_l,r²/210²,
w=<fg>=2 sum_r t_l,r/210²,
c=<f²g²>=sum_rs t_l,r T_rs t_l,s/210⁴,
T_rs=integral S² x_r x_s.
```

The code computes the integer T matrix in two independent ways: the
actual exterior-weight Weyl constant term, and fourfold SU3/SU2 tensor
singlet multiplicities with exact integer hypercharge. They agree and
sum_rs T_rs=714. No continuum matter representation is substituted.

A product of four distinct elementary plaquettes has an edge in only
one loop. Each of the other three plaquettes can cover at most one of
its four edges. The same unique-edge argument removes multiplicity
patterns3+1 and2+1+1. For two distinct doubled plaquettes, integrating
one private link in each proves that their holonomies are independent
Haar variables even if the plaquettes share a link. For P plaquettes,

```
<M² M_t²>=P c+P(P-1)(u v+2w²).                      (5)
```

All odd total hopping moments vanish pointwise under the bipartite
sublattice sign operation. The second hopping moments z,z_t,q_h are
pointwise constants. Define z_l=P u, s_l=P v and q_l=P w. Equations
(4)-(5) give the exact parent mixed moment

```
M_exact=<V²V_t²>
 =z z_t+2q_h²-G(L)sum_r d_r t_h,r²
  +z s_l+z_t z_l+4q_h q_l+P c+P(P-1)(u v+2w²).       (6)
```

This replaces an extensive operator-norm estimate by explicit local
Haar and CAR contractions, with all connected and disconnected terms
retained.

## Full finite residual and its boundary

Let P_ell be the complete product word projection. Both Omega and
V_t Omega belong to every selected cutoff. Thus

```
||V_finite V_t Omega||²
 =||P_ell V_parent V_t Omega||² <= M_exact.           (7)
```

The vector V_parent V_t Omega has word degree at most two on each link.
For ell>=2, (7) is equality. At member zero ell=1, it remains a rigorous
upper bound. Gauge invariance places the full vector and the retained
projection in the physical carrier; the full outside space is counted.

Use the unchanged trial notation q=q_h+q_l,
s=z_t+s_l, g=P/22050, C_loop=-N/420 from its exact compressed matrix,
and rational Schur-root endpoint d_lo. Its Rayleigh shift is e. The
exact projection of V V_t Omega into the eleven-vector trial space has
squared norm q²+g||C_loop t_l||². H0 maps that space into itself, so

```
R_out² <= M_exact-q²-g||C_loop t_l||²,
||(Hhat-E00-e)psi||²
 <= [(q+e)²+(d_lo+e)²s+R_out²]/(1+s).                (8)
```

For ell>=2 the squared residual in (8) is exact for the declared
rational trial. The square-root display is rounded outward, and source
energy divides the dimensionless residual by the actual a. No trial
eigenvalue is promoted to a full physical eigenvalue. The spectral
theorem puts some full eigenvalue inside the displayed residual window
and bounds this auxiliary trial's weight outside radius R by r²/R².
It supplies no lower bound on the selected density's particle residue.
Ground isolation needs a separate excluded-complement spectral estimate;
a rest mass also needs the same source's observable branch and momentum
attachment. The ground upper bound, clock and preparation are unchanged.
