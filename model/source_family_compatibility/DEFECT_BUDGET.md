# Sparse occupation defects and complete repair records

Use the two-flavour full-Fock repair in
`model/defect_field_bridge/full_domain.py`, including its collision hold.
This is the declared full-domain extension of the binary repair, not a
derivation of that extension from the original binary rule. The coherent
propagation is the signed source Hubbard action. Both operations use the
same occupation basis. All coin, seam, law and raw records are retained.

## Exact one-attempt comparison on the whole Fock space

Let d_x=1_(n_x!=1), D=sum_x d_x, and let Q_e be the projection onto at
least one defective endpoint of seam e. The complete repair isometry V_e
agrees with the raw-zero identity comparison V0_e on 1-Q_e. Thus

    A_e=V_e-V0_e=A_e Q_e,    A_e* A_e <= 4 Q_e.             (1)

The actual branch maps give a stronger identity. Let R_e project onto
total endpoint occupation one, and C_e onto at least one doubly occupied endpoint.
They are disjoint. With E-law probability theta in [0,1],

    A_e* A_e = (1+theta) R_e + 2 C_e <= 2 Q_e.             (2)

Indeed raw zero is always an identity hold with positive sign. On R_e,
half the C coins and no E coins have raw zero. A collision writes raw
eight for either law. Both empty endpoints and both singly occupied
endpoints write raw zero for every coin. The other raw branches are
orthogonal to the identity comparison archive. This proves (2) on all
matrix units, including spin coherences and coherences between occupation
sectors. CAR signs in nonzero raw branches cannot alter this Gram identity.

On a graph with E edges, use one uniform global seam attempt. The archive
retains the seam label, so the squared difference is the average of (2).
Since sum_e Q_e <= sum_x deg(x)d_x <= d_max D, every normalized state rho,
with arbitrary reference, obeys

    ||(V-V0) purification(rho)||^2
           <= 2 d_max Tr(rho D)/E.                       (3)

This bound has no one-hole or no-doublon restriction.

## Complete histories with interacting propagation

Interleave m repairs with identical coherent source propagations, or
identical reference-stable channels with the same record-adapted policy.
Let rho_j^0 be the complete ideal source/archive state immediately before
attempt j, with every preceding repair replaced by V0. Write
b_j=Tr(rho_j^0 D), tracing all other registers. Telescoping the complete
dilations with the ideal prefix and actual suffix gives

    ||rho_actual-rho_ideal||_1
       <= min(2, 2 sum_(j=1)^m sqrt(2 d_max b_j/E)).        (4)

Each suffix is contractive, and each local term is bounded using the
ideal-prefix state in (3). For mixed inputs use a common purification;
the pure-state trace-norm difference is at most twice the vector norm.
Record commitment, FIFO export and common subsequent channels preserve
the bound. The ideal archive has raw zero but retains the same seam,
law and coin distribution. The actual archive retains all raw outcomes;
no postselection occurs. Exact branch-map statements do not certify a
hardware or native pulse implementation error, which must be added.

The bound uses the entire propagated state and therefore remains valid
when hopping creates doublon--hole pairs. A fixed initial bound on D
alone is insufficient; the expectations b_j must be controlled along
the ideal evolution. Equation (4) also gives the weaker coefficient
4 sum_j sqrt(d_max b_j/E) directly from (1).

## A Hubbard energy budget bounds every ideal defect expectation

Write the existing source exactly as

    K=T+U D_double-(U/2)N,
    H=K+U|V|/2=T+(U/2)D,    U>0,                         (5)

because D=2D_double-N+|V| on every Fock basis vector. Each signed edge
has hopping amplitude c/2 per flavour and operator norm c, so
||T||<=cE. Any state with Tr(rho H)<=C, C>=0, satisfies

    Tr(rho D) <= min(|V|, 2(C+cE)/U).                    (6)

This does not require U>=4cE. The full ground energy E0 of H obeys E0<=0,
by the Rayleigh quotient of any singly occupied basis vector. Thus a
positive ground-relative budget Tr(rho(H-E0))<=C implies the premise
of (6), without assuming a ground sector or uniqueness. Evolution by
this same H conserves its expectation for every duration. Energy-preserving
record-adapted controls also preserve the averaged budget. If other ideal
controls inject energy, their certified cumulative mean work must be
included in C at each attempt. Conservation is required only for the
ideal comparison; actual repairs need not conserve H.

A local estimate improves (6) when U is large compared with degree,
rather than graph volume. For P_e=1-Q_e, the actual hopping satisfies
P_e T_e P_e=0. Cauchy--Schwarz for a purification gives, with
q_e=Tr(rho Q_e),

    |Tr(rho T_e)| <= 2c sqrt(q_e(1-q_e))+c q_e.

For b=Tr(rho D), summing and using sum_e q_e<=d_max b yields

    Tr(rho H) >= a b-2c sqrt(E d_max b),
    a=U/2-c d_max.                                      (7)

When a>0, Young's inequality gives the rational bound

    b <= 2C/a + 4c^2 E d_max/a^2.                       (8)

The exact quadratic-root bound is
`b <= ((c sqrt(E d_max)+sqrt(c^2 E d_max+a C))/a)^2`.
The producer uses the minimum of (6), (8) and |V|. For U>=4c d_max,

    b/E <= 8C/(UE)+64c^2 d_max/U^2.                     (9)

Combining this with (4) controls a complete interacting history on any
duration for which the same energy budget holds. For C=0 the trace error
is at most `min(2,16 sqrt(2) m c d_max/U)`, and the rational upper bound
`min(2,24 m c d_max/U)` is valid. Thus U/(m c d_max) tending to infinity
is a sufficient repair-comparison condition for this energy class.
It does not select that trajectory or show that original odd charged
modes have finite energy on a spin clock. The separate full-ground
projection exclusion retains its own hypotheses and conclusion.

## Exact finite controls and physical scope

The producer reconstructs every signed Hubbard number block and every
full-domain repair column on two- and three-site graphs. It checks the
whole-Fock Gram identity (2), all padded raw slots, the defect identity
(5), and P_e T_e P_e=0. Full matrices retain coherent empty, single and
double occupation entries; no restricted repair channel is substituted.

An exact dimer control uses c=3/2, U=4, and an initially singly occupied
opposite-spin basis state. Its source H expectation is zero. Actual
Hubbard evolution for time pi creates a doublon--hole probability 9/50,
so b=9/25. The full one-repair dilation difference has squared norm 9/25,
independently of theta. Thus the comparison includes genuine leakage
from the singly occupied space while its energy budget is conserved.

These are same-source operational error bounds for a declared action,
state and schedule. They do not supply an autonomous state selection,
a physical duration per global repair attempt, a relativistic limit,
particle quantum numbers or a mass. Neither a Gibbs reference nor an
empirical energy scale enters the construction.
