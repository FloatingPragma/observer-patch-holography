# Local influence on the selected physical source

## 1. The exact source and the reused theorem

Take a member of the
[primary source contract](../particle_common_contract/PROOF.md), with
spacing `a=2^-n`, side `2^(2n+1)`, per-link degree `ell=n+1`, sixty CAR
modes per cell and the selected Gauss-constrained polynomial state. The
source clock is `exp(-itH)`; the conjugate reference has idle Hamiltonian
zero and is retained in each observer output. Its exact preparation may
be spatially entangled. The locality theorem below is uniform in the
input, so it imposes no product-state or thermal-state assumption.

The [finite Wilson cone](../particle_physical_causality/PROOF.md),
sections2-4, proves an all-state complete-record bound on the
coordinate/Fock parent. Its source ingredients are the exact native
Dirac/M15 bond spectrum and the even-interaction CAR commutator-chain
argument in [the field locality proof](../particle_field_causality/PROOF.md).
The present consumer verifies that the primary product-cutoff Hamiltonian
satisfies the same hypotheses directly. No strong-limit comparison of
the chosen prepared state with an uncut state is needed.

Each cell owns its three outgoing links. A hopping term is supported
on its two endpoint cells. A plaquette can be supported on its four
vertices; with link ownership its minimal support has three cells and
Manhattan diameter two. Terms extending outside the open box are absent.
Electric terms on all owned links, including links to ghost endpoints,
are onsite under this ownership convention.

## 2. Product compression preserves the local interaction bounds

Write the product word projection as a tensor product of single-link
projections, with identity on all fermion factors. For a parent term
Phi_Z supported on owner cells Z, its restriction to the cutoff carrier
is

```
Phi_Z,ell = P_Z,ell Phi_Z P_Z,ell tensor I_outside,
supp Phi_Z,ell subset Z,
||Phi_Z,ell|| <= ||Phi_Z||.                              (1)
```

The finite source forms each complete parent product before compression.
Thus (1) applies to its actual plaquette and bond operators. It would not
follow from compression by a total-degree projection over all links.
The present source explicitly uses the product projection.

For each native Wilson bond, its complete two-site one-body spectrum is
`+1/a` and `-1/a`, each with rank30, and zero with rank60. The full Fock
norm is therefore `30/a`, the sum of all positive one-body eigenvalues.
The coordinate representation is unitary pointwise; compression cannot
increase this norm. The normalized real plaquette trace has norm at most
one, so its actual coefficient gives `||Phi_p||<=1/(7a)`.

The onsite contribution contains the electric operator with coefficient
`105/(2a)` and the onsite Wilson matter term. The onsite interaction
picture is a product of cell unitaries, preserving the support, parity
and norm of every interaction. Although electric norms grow with ell and
1/a, their growth contributes no spatial interaction-chain weight.
This statement uses the exact finite-cutoff onsite operator; no electric
energy tail of the prepared state is discarded.

For spatial weight lambda=2, six incident bonds and at most twelve
incident plaquettes give

```
sup_x sum_(Z contains x) |Z| ||Phi_Z,ell|| 2^diam(Z)
 <= W_star/a,
W_star = 360*2 + 48*(1/7)*2² = 5232/7.                  (2)
```

This bound is uniform in volume and ell. The finite incidence control
independently enumerates every interior open-box bond and plaquette,
using their actual owner supports. Its smaller weighted sum satisfies
(2). The imported source bond calculation retains all fifteen internal
columns and every fermion occupation.

## 3. Physical restriction without a nonlocal Hamiltonian substitution

Let P_G be the product of local Haar constraints at all link endpoints.
The source Hamiltonian and each gauge-invariant instrument commute with
P_G. The state and both sides of the comparison are supported there.
The ambient source bound therefore applies to these same experiments
restricted to the physical carrier, since restriction cannot increase
an operator norm. The argument does not declare `P_G H P_G` to be a
sum of local operators on a newly factorized physical space. Locality
is proved on the ambient tensor/CAR carrier before using physical
support. The physical subspace itself need not factor by region.

The same distinction applies to gauge-invariant local apparatus. Its
ambient local isometry is specified first, then restricted to physical
inputs. Bare charged local operators without compensating references
are outside this instrument domain. The complete boundary-center
instrument used below is gauge invariant without a charged auxiliary
register.

## 4. Complete records, including the unchanged quantum reference

For disjoint owner supports X,Y at graph distance D, the reused
commutator chain gives, for even local contractions,

```
||[alpha_t(A_X),B_Y]||
 <= 2 |X| 2^(-D) [exp(2 W_star |t|/a)-1].              (3)
```

Take the same selected canonical purification in two experiments. One
applies a remote even gauge-invariant unitary at time zero and the other
applies identity. Both evolve with the same H and time T. Their observer
outputs comprise the local field output, every local apparatus record,
and the untouched conjugate reference. The remote field and any remote
private apparatus are traced out; directly delivering its record is a
different observation domain.

The complete-history proof allows arbitrary reference dimension and
input correlations. Adjoin its idle register to the observer cell X:
source interactions act trivially on that factor, so the same support
and norm estimates hold for all joint observer/reference effects.
For k intermediate local interventions followed by
a final local output on X, duality and telescoping give

```
D_half(outputs) <= min(1,2(2k+1)|X|
                      exp[-(log(2) d-2 W_star T)/a]),
d=aD.                                                   (4)
```

All local effects or their retained Stinespring dilations have norm at
most one. The number of outcomes does not multiply this bound. For the
single final observation used here, k=0 and the minus-one factor in (3)
can be retained:

```
D_half(outputs) <= min(1,2 |X| 2^(-D)
                            [exp(2 W_star T/a)-1]).      (5)
```

At T=0 the outputs agree exactly, including the reference, because the
remote unitary commutes with every operator of the observer/reference
output algebra. This is an all-state finite non-signalling identity.
At nonzero time the bound has tails and does not assert exact spacelike
commutativity.

## 5. An actual complete physical boundary measurement

Choose X to be the origin-owned cell, containing its sixty fermion modes
and three outgoing links. The
[physical boundary-state construction](../particle_physical_boundary_state/PROOF.md)
provides its complete endpoint hypercharge projectors E_q and neutral
pointer unitary `V=sum_q E_q tensor S_q`. Each E_q is an even local
projector of norm at most one; V is an even unitary on X and its pointer.
Because the endpoint hypercharge generators are central, this apparatus
commutes with all source Gauss constraints. All `|S_ell|³` pointer labels
are retained, including zero-charge and zero-probability outcomes.
Nonabelian irreducible labels are not inferred from these hypercharges.

The remote operation is `exp(i B_p)`, with
`B_p=P_ell ReTr_M15(U_p)P_ell/15`, on an actual interior plaquette in the
directions e2,e3. This is even and gauge invariant; `||B_p||<=1`.
The default anchor is `(side-1,0,0)`; its three owner cells have minimal
distance `D=side-1` from X. Each link and all four vertices are inside
the selected box. The finite geometry validator rejects overlapping
supports or a plaquette crossing the box boundary. These operations
preserve the selected physical carrier and use no global energy or
momentum projection as a purported local detector.

For fixed source duration `T=1/10000`, (5) gives nontrivial bounds for
the finite selected members. If `x=2 W_star T/a<=2`, a positive Taylor
sum with an explicit geometric tail bounds `exp(x)-1` using rational
arithmetic. Otherwise, set `mu=1/2<=log2` and use

```
g=mu*d-2 W_star T,
min(1,2 exp(-g/a)) <= 2^(-max(0,floor(g/a)-1))           (6)
```

when g>=0; a negative margin gives the trivial bound one. The short-time
calculation retains a factored `coefficient*2^(-D)` and never allocates
an integer of size exponential in the graph distance. At levels0,1,2,3,
the default complete-output half-trace bounds are respectively
`2^-2,2^-7,2^-30,2^-124`. These are conservative analytic bounds, not
executed full spatial histories or observed particle probabilities.

A separate comparison holds the assigned physical separation fixed.
At level6, `a=1/64`, choose the remote anchor `(64,0,0)`, so d=1.
The same fixed duration gives `g=3067/8750` and bound `2^-21`.
The code authenticates that actual finite member and its complete
boundary measurement; no same-state assertion between different levels
is used. Exact preparation and exact apparatus have no additional
implementation error here. Approximate implementations require their
own complete-output errors on both compared experiments.

## 6. The closed premise and its mathematical boundary

The selected state, reference, source action, product cutoff and complete
boundary measurement satisfy a finite local-record theorem on the same
clock. Their interaction weight obeys (2) at every selected level.
Consequently the conservative cone in assigned coordinates is uniform
under this source scaling. Its speed bound is
`2 W_star/log2`, with rational upper bound `20928/7`; it is not identified
with a physical light speed or the quadratic characteristic speed.

If independently identified joint state and observable records converge
under refinement, (4) makes their remote influence vanish outside this
cone. No such joint convergence is obtained by applying a state-uniform
bound separately at each level. This result supplies no Lorentz action,
physical mass shell, scattering limit, nontrivial interacting
four-dimensional continuum or Einstein dynamics. Its preparation and
apparatus specifications may have exponential implementation cost.

The typed consumer regenerates the entire primary member before deriving
its consequences. Changing the source, reference, clock, preparation,
geometry or claimed result fails authentication. Parent scientific
sources and their cutoff domains are unchanged.
