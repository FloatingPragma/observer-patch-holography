# The physical state and its boundary-charge records

## 1. The selected finite state

Take a member of the
[common source contract](../particle_common_contract/PROOF.md), with its
open spatial box, sixty fermion modes per cell, native fifteen-state gauge
representation, coordinate links and product word cutoff. The actual
Hamiltonian, physical projection and prepared state are

```
H = P_ell [105 sum_e C_e/(2a) + dGamma(beta_D D_W)
                - sum_p Re Tr_M15(U_p)/(105a)] P_ell,
P_G = product_x integral_G U_x(g) dg,
Q = P_G [T_m(beta(U I-H)/(2s))]^s,
rho = Q²/Tr Q².                                             (1)
```

The finite parameters, positive Taylor polynomial, spectral interval
and canonical basis are fixed by that member. Its conjugate reference
has idle Hamiltonian zero; its source clock is `exp(-itH)`.
The [Gauss preparation proof](../particle_gauss_ground_state/PROOF.md)
constructs Q by finite matrix operations and a deterministic Householder
map from the physical blank to `vec(Q)/||Q||_F`. The exact map is an
explicit preparation primitive. No acceptance conditioning defines rho.
Both source and reference lie in their separate Gauss kernels. The
source density commutes with H; its purification need not be stationary
under the source-only clock.

The constraints include every endpoint of every owned link. Outside
endpoints have zero external charge and no matter register. Empty matter
and constant link functions supply a nonzero physical blank. On the
physical space, every eigenvalue of Q is strictly positive: the argument
of T_m is positive semidefinite and its constant coefficient is one.
Thus `supp rho=P_G`, without requiring a gap or a unique ground state.

Write `D<=2^b` for the certified ambient dimension bound, `W=U-L`, and

```
N = b + ceil(2 beta W),     c=2^(-N).
rho >= c P_G,
rho_R=Tr_out rho >= c Tr_out(P_G).                         (2)
```

Indeed `1<=Q|phys<=exp(beta W/2)`, so `Tr Q²<=D exp(beta W)`;
`exp(-x)>=2^(-ceil(2x))` for x>=0. This is the same bound consumed by
[the regional modular proof](../particle_spatial_modular_response/PROOF.md).
It proves faithfulness on the actual regional support, rather than on
the whole ambient regional matrix algebra. The possibly tiny c is kept
as a factored dyadic; no exponential-size integer or matrix is allocated
by the member calculator.

## 2. Outer boundary and regional boundary are different constraints

Choose the region R to be the origin cell, its sixty fermion modes, and
its three owned outgoing links. Every selected member has side at least
two. These links end at the three interior vertices `e1,e2,e3`. Each
endpoint gauge transformation factors between the inside and outside
registers. Gauss invariance of a full physical vector requires its two
boundary representations to pair to a singlet. It does not require the
inside factor alone to be a singlet.

Let `Q_j` be the normalized integer hypercharge generator acting on the
right endpoint of the inside link from the origin to e_j. This is the
central U(1) generator of the gauge group at that endpoint. It commutes
with the full local gauge action, including its nonabelian factors.
The three Q_j commute. A physical vector has total endpoint charge zero,
so an inside outcome q_j is paired with outside charge -q_j. The source
state contains this correlation before any readout. It requires no
additional charged preparation reference. The separate outer ghost
vertices retain the zero-external-charge boundary condition in (1).

The cut convention assigns each link to its owner, as in the regional
restriction theorem. It does not split one link into extra edge degrees
of freedom. A different cut must carry its own factorization map.

## 3. Complete finite charge readout

The native irreducible charge and dimension pairs are

```
(q,d) = (-4,3), (-3,2), (1,6), (2,3), (6,1).             (3)
```

On a native coefficient `R_r(U)_(ab)`, right endpoint multiplication by
`exp(-itQ6)`, with integer-normalized `Q6=6Y`, acts with charge -q;
its complex conjugate has charge +q.
The constant word has charge zero. For the full per-link word space
V_ell, define

```
A={0, +/-1, +/-2, +/-3, +/-4, +/-6},
S_ell={a1+...+aell : a_j in A},   S_0={0}.               (4)
```

This is the exact endpoint charge spectrum. Every word of degree at
most ell has a charge in (4). Conversely, for each sum choose diagonal
native entries or conjugates with the selected charges and multiply
them, padding by the constant word. Each diagonal entry equals one at
the identity, so this product is not the zero function. Its Haar norm
is positive and its endpoint charge is the prescribed sum. Linear
dependencies in the Haar Gram quotient cannot remove an entire
eigenspace with this nonzero witness.

For r in S_ell let

```
P_r(Q_j)=product_(t in S_ell, t!=r) (Q_j-t I)/(r-t),
E_(r1,r2,r3)=P_r1(Q_1) P_r2(Q_2) P_r3(Q_3).             (5)
```

Functional calculus makes these orthogonal projectors. They sum to the
identity on the entire regional carrier. All `|S_ell|³` tuples are
retained, even when their weight in a physical state vanishes. The
neutral record apparatus

```
V psi = sum_r E_r psi tensor |r>_record,
V* V = I,
Pr(r)=Tr(rho_R E_r)                                    (6)
```

is a complete isometry. Enumerate the outcomes and let S_r be the cyclic
record shift taking blank label zero to label r. The explicit unitary
`sum_r E_r tensor S_r` realizes (6), since mutually orthogonal E_r sum
to identity. It commutes with every full source Gauss action
because the Q_j are central boundary charges. This is a physical
measurement of regional electric charge with compensating outside
correlations retained. The usual source clock acts before or after this
specified apparatus operation. A charge readout is not assumed to
commute with the interacting Hamiltonian or to have zero work cost.

The menu resolves the full boundary hypercharge tuple. Distinct
nonabelian irreducible representations may share a tuple, so this is a
coarse graining of the complete nonabelian boundary center. Its exact
resolution neither discards such representations nor computes their
full spatial multiplicities. Their finite computation can instead use
the retained Haar basis, internal Gauss kernels and boundary
intertwiner spaces. No multiplicity inferred from dimension is used
in (2) or (5).

`charge_projectors` checks exact functional calculus on a supplied
finite charge matrix, including Hermiticity, orthogonality, completeness
and the eigenvalue relation. The executed source control uses the
complete V_1 Haar coefficient space. Character orthogonality gives

```
dim V_1 = 1+2 sum_r d_r² = 119.                         (7)
```

The constant, every native matrix coefficient and every conjugate are
retained. Its eleven exact spectral projectors have ranks 1 and d_r²
for the corresponding nonzero charges. This 119-state link calculation
does not execute the full spatial density.

## 4. Actual nonzero charge sectors in the prepared spatial state

Use the elementary plaquette with vertices
`0,e1,e1+e2,e2`, and one native irreducible representation r. Its
character state on empty matter and otherwise constant links is

```
L_r = chi_r(U_(0,1) U_(e1,2)
                       U_(e2,1)^(-1) U_(0,2)^(-1)).     (8)
```

All four links are distinct. A closed trace is invariant under all
endpoint gauge transformations, and each link occurs at degree one
(native or conjugate). Thus `P_G L_r=L_r` and `P_ell L_r=L_r` for every
selected member. Product Haar integration reduces its norm to the
character norm on G. The native blocks are irreducible; the executable
control independently verifies that their twelve-generator commutants
have dimension one. That generator packet uses exterior-mask order and
central current `-i Q6/2`; the identified charge packet has a different
basis order. The control obtains invariant blocks from the generator
packet itself and checks the normalized charges and dimensions in (3).
It never applies one packet's indices to the other. Schur orthogonality
gives `||L_r||=1`.

The two links owned by the origin have boundary charges `(-q,+q,0)`.
The outside path has the opposite pair. Replacing r by its conjugate
reverses these signs. All ten nonzero tuples obtained from (3) and its
conjugates are distinct; the physical blank gives the eleventh tuple
`(0,0,0)`. These are actual states of (1), not an independently supplied
representation-pair experiment.

To compute their regional reductions, cyclically reorder the trace as
`chi_r(A B)`, with A the two-link inside path and B the outside path.
Haar multiplication preserves Haar measure. The matrix coefficients
`sqrt(d) R_r(A)_(ab)` and their outside counterparts are orthonormal.
Consequently (8) has d² Schmidt terms, each of coefficient 1/d:

```
Tr_out |L_r><L_r| = (1/d²) Pi_r,   rank Pi_r=d².         (9)
```

The inside boundary representation is the product of the appropriate
dual/native representations at e1 and e2. Its dimension d² describes a
regional charge span, not a particle multiplicity. Since L_r is a
normalized physical vector, `P_G>=|L_r><L_r|`. Equations (2) and (9)
therefore give actual selected-state statements:

```
Pr_rho(boundary tuple (-q,+q,0)) >= c > 0,
rho_R >= (c/d²) Pi_r.                                  (10)
```

The blank has the analogous lower bound c. All eleven sectors have
strictly positive weight in the exact selected polynomial state. The
bounds neither sum to one nor claim to list every occupied sector.
The complete normalization and all other outcomes are supplied by
(5)-(6). A low-temperature limit may lose some of these weights;
no volume- or refinement-uniform lower bound follows from (10).

## 5. Finite preparation and refinement have separate domains

For each selected level, (1)-(10) give an explicit physical state,
deterministic preparation map, outer boundary condition, actual regional
support and a complete boundary-charge readout on one fixed clock.
`common_contract_payload` reconstructs and authenticates that source
member before deriving these consequences; the bound receipt includes
its source, state, reference and clock identifiers. Mutation of any
member parameter or result is rejected by recomputation.

The construction is finite and effective as an exact matrix algorithm.
Generic state preparation and full representation decomposition can
have exponential volume cost. The full sixty-mode-per-cell spatial
matrices and their gate tape are not executed by this certificate.
This is compatible with resolving the finite physical-state preparation
requirement under the declared architecture.

The selected sequence changes volume, spacing, word cutoff and cooling
parameters. Its states are not proved to agree under a refinement
embedding or converge jointly. The fixed-spacing, fixed-cutoff local
ground subsequence and the separate fixed-volume link-cutoff theorem
have their stated domains. Neither can be interchanged with the other
limit without uniform estimates. A physical relativistic field, a
positive interacting four-dimensional mass gap and a particle dictionary
require those additional conclusions.
