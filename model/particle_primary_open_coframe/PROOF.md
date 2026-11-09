# Boundary-completed strain of the open primary source

A declared constant-coframe coupling can reproduce the actual open Wilson
source exactly at its flat point. Its six spatial derivatives and complete
mean readouts then use the selected source, preparation, reference and clock.
This supplies an action-matched spatial response. The coframe is a declared
architecture input; no proper geometry, Lorentz stress or Einstein equation
is inferred from this construction.

## The boundary completion

Use the [primary contract](../particle_common_contract/PROOF.md), including
all sixty CAR modes per cell, product link-word compression, physical Gauss
support and the polynomial state. Work first in the parent coordinate
algebra. Each covariant shift T_i on the open matter box is a partial
isometry. Its initial-space defect is the positive projection
Pi_i=I-T_i†T_i, supported on one missing face. Let

```
D_i=(T_i-T_i†)/(2i),  K_i=I-T_i,
L_ij=(K_i†K_j+K_j†K_i)/2 + delta_ij Pi_i.
```

In particular L_ii=2I-T_i-T_i†. This adds precisely the term missing from
an unmodified periodic-coframe substitution. For every real positive
matrix C,

```
sum_ij C_ij L_ij = sum_ij C_ij K_i†K_j + sum_i C_ii Pi_i >= 0.
```

The first equality uses the symmetry of C. The first term is a positive
Gram form and the second is a sum of positive boundary projections. The
shifts, differences and the added projections transform by conjugation
under arbitrary local source gauge transformations.

Declare an orientation-preserving real coframe V, and set G=VᵀV and
v=det V. With the actual source Clifford matrices beta and
alpha_A=i beta gamma_A, define the complete matter one-body action

```
h(V)=sum_Ai alpha_A (V^(-T))_Ai D_i/a
                   + beta sum_ij (G^(-1))_ij L_ij/(2a).
```

At V=I this is exactly the primary onsite3 beta/a plus the positively
oriented hopping(-beta+beta gamma_i)T_i/(2a) and its adjoint. The boundary
completion is made before second quantization and product link compression.
It is not a modification of the flat source or its energy zero. The full
Fock space is retained. The identity holds for arbitrary parent native
link coordinates, rather than just the finite phase controls used by tests.

The added term relative to the naive periodic formula at the flat point
is dGamma(beta sum_i Pi_i)/(2a). Its physical norm is45L²/a, as established
in [the boundary audit](../particle_primary_gravity_premises/PROOF.md).
The independent control checks the signed difference, not just its norm.

## Electric and magnetic couplings on the actual open graph

Let E_i^r(x) be the left electric generators of the three outgoing links
owned by x, and set

```
E_ij(x)=sum_r (E_i^r(x) E_j^r(x)+E_j^r(x) E_i^r(x))/2,
H_E(V)= e/(2a) sum_xij (G/v)_ij E_ij(x),  e=105.
```

This is the same positive electric Gram prescription as the
[coframe theorem](../particle_geometric_response/PROOF.md). It preserves
every Gauss constraint. Ghost endpoint constraints force their dangling
links to constants; electric generators annihilate those factors. All
owned links and constraints are retained.

At each owned vertex use the three cyclically oriented plaquettes, when
all their vertices lie inside the box. If plane i is absent, put C_i=0;
otherwise let C_i=I-R(U_p), using the actual fifteen-row representation.
Define

```
Q_ij(x)=Re Tr(C_i† C_j)/30 - delta_ij 1_(plane i present),
H_B(V)= kappa/a sum_xij (G/v)_ij Q_ij(x),  kappa=1/7.
```

All loops have the same gauge base point. The diagonal of Q is exactly
-ReTr R(U_p)/15 when the plaquette is present, and zero otherwise.
Consequently H_B(I) is the literal open primary magnetic potential.
No torus plaquette or missing-face scalar is added. Cross terms are a
declared strain coupling; complete parent products are compressed only
after multiplication. The electric and magnetic terms are not second
quantized as matter one-body matrices.

Thus the full physical finite operator is

```
H(V)=P_word [dGamma(h(V))+H_E(V)+H_B(V)] P_word,
H(I)=H_primary.
```

Every term is Hermitian, even and gauge preserving. The same physical
Gauss carrier admits every V. The positive electric forms and bounded
remaining terms give finite matrices on each product cutoff. An internal
constant frame rotation V↦RV has its inherited Clifford spin lift; no
Lorentz boost, proper-metric identification or spatial-coordinate tensor
law beyond the specified coupling is asserted.

## Six derivatives and dilation

For real symmetric S and V(t)=I+tS, differentiate at zero. The coefficient
identities are

```
d(V^(-T))/dt=-S,   d(G^(-1))/dt=-2S,
d(G/v)/dt=2S-(Tr S)I.
```

The matter strain includes both parts of L:

```
h'_S= -sum_Ai alpha_A S_Ai D_i/a - beta sum_ij S_ij L_ij/a.
```

The boundary contribution relative to the uncompleted derivative is
-dGamma(beta sum_i S_ii Pi_i)/a. The electric and magnetic derivatives
replace G/v by2S-(Tr S)I in their displayed formulas. These are six actual
Hermitian observables K_S=∂H(V(t))/∂t on the unchanged flat source.

For isotropic strain S=I, let H_D and H_W denote the flat Dirac and the
completed Wilson parts. The exact scaling along V=lambda I is

```
H(lambda I)=lambda^(-1)(H_D+H_E+H_B)+lambda^(-2) H_W,
K_I=-H_D-2H_W-H_E-H_B=-H-H_W.
```

With the declared coordinate volume V_cells a³ lambda³, the corresponding
mechanical pressure mean is Tr(rho(H+H_W))/(3 V_cells a³). This defines a
response under the stated volume convention. It supplies no equation of
state or equality with proper geometric pressure. The polynomial
preparation is held fixed when taking this response; no equilibrium
readjustment or alternative Gibbs state is substituted.

## Full-carrier bounds and complete records

Let N=60 V_cells, ell be the source link degree, and
C_S=2S-(Tr S)I. Write |S|_1 for the sum of absolute matrix entries,
d_S=sum_i |S_ii|, c_d=sum_i |(C_S)_ii| and
c_o=sum_(i≠j)|(C_S)_ij|. The source norms obey

```
R_m = N/a (5 |S|_1+d_S),
R_e = V_cells e/(2a) (29 ell²/5) (c_d+c_o),
R_b = V_cells kappa/a (c_d+2c_o),
||K_S|| <= R_S=R_m+R_e+R_b.
```

Here ||D_i||≤1, ||K_i||≤2 and ||Pi_i||≤1 give ||L_ij||≤4 for
off-diagonal entries and at most5 on the diagonal. Full second
quantization is bounded by N times the one-body norm. The inherited
per-link Casimir cap is29 ell²/5; row-operator Cauchy-Schwarz bounds
||E_ij|| by that cap. The magnetic diagonal is bounded by1, and each
off-diagonal entry by2, including every absent-plane pattern. Compression
and restriction to the physical support contract these norms. The bound
is conservative and does not require constructing the enormous matrix.

For nonzero S, the effects E_±=(I±K_S/R_S)/2 and Kraus operators sqrt(E_±) give a complete
physical two-record instrument, with numerical labels±R_S. At S=0 both
effects are I/2 and both numerical labels are zero. Its original
mean is exactly Tr(rho K_S); both signs and every idle quantum reference
are retained. A complete-output half-trace error e changes that mean by
at most2R_S e. Exact positive-root services are part of the declared
apparatus; the large primary strain histogram is not executed. Measuring
this observable can disturb the input state.

## Finite-theorem and physical domains

The local energy Ward theorem and primary current construction apply at
the identical flat Hamiltonian. The finite modular-work comparison and
supported marginal first law also retain exactly their selected state
and clock. No continuity or entropy theorem is reinterpreted as an area
law. The matched Einstein theorem requires proper geometry, a physical
diamond modular kernel, the area coefficient and common limiting errors.
These objects are not produced by specifying H(V).

The coupling need not be unique. Adding (Tr G-3)A for any physical bounded
Hermitian observable A leaves H(I) and every fixed-flat experiment
unchanged but shifts K_S by2(Tr S)A. This is an explicit architectural
choice, not an obstruction requiring a uniqueness theorem. The
construction above fixes one coupling and derives its consequences; a
claim that proper geometry selects it requires the corresponding source
identification. Neither this spatial response nor a gravity proof is a
prerequisite for converting an independently identified rest energy to GeV.
