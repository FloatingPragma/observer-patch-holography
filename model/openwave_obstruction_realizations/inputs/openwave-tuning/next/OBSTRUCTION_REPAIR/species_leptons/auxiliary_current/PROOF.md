# Full auxiliary dressing changes the axial current, but not enough

This is an independent current calculation on the hadron worker's frozen
local source:220 three-quark CAR states times25 auxiliary states, dimension
5500. The auxiliary consists of one scalar state and an SU3-adjoint physical
spin vector. Its Hamiltonian is

    H=ΔΠ1+g Σai Gai (|ai><0|+|0><ai|),
    Gai=Σflavour c†(Ta tensor σi/2)c,
    Δ=3κ, g=1/8, κ∈{1/2,1,2}.

Every coefficient was fixed by the owning source before the axial
comparison. This source is a new local architecture. It has not been joined
to the two-cell vacuum, spatial Dirac dynamics, or a physical weak clock.

## 1. Actual dressed states and exact current

For a normalized colour-singlet triquark v with spin/isospin1/2,
ΣGai² v=4v. The lower state is exactly

    D(v)=c|0>v − (s/2)Σai |ai>Gai v,
    s²=(1−Δ/sqrt(Δ²+16g²))/2, c²=1−s².

The auxiliary components are actual quark colour-adjoint states. They are
not replaced by a scalar energy shift. Total colour includes the auxiliary
adjoint generators and total angular momentum includes its spin-vector
generators. The run checks the full Hamiltonian residual, total Gauss
annihilation, and total spin projection on all four dressed nucleon states.

V+=Σc†u cd commutes with every Gai, so its dressed normalization is exactly1.
For Aa+=Σc†σa τ+ c, the one-body CAR commutator identity and the exact
colour/spin Casimirs give

    Σbi [Gbi,[Gbi,Aa+]] = (8/3) Aa+.

Between the two bare spin-half singlets, ΣG²=4 on both sides. Therefore

    <p|ΣG Aa+ G|n> = (4−4/3)<p|Aa+|n>.

The normalized auxiliary component has axial ratio(2/3)(5/3)=10/9. The
complete dressed current is consequently

    gV=1,
    gA=(5/3)(1−s²/3).

This is an exact full-state result. The script independently builds all
220-state CAR bilinears using Gaussian integers, reconstructs and verifies
an exact colour-singlet spin/isospin highest-weight quartet, and verifies
the bare5/3 and auxiliary10/9 matrix elements as exact fractions. It also
checks the full dressed2×2 vector/axial spin blocks numerically against
this formula.

For the frozen source members:

| κ | Auxiliary probability s² | gA |
| --- | ---: | ---: |
| 1/2 | 0.025658350975 | 1.652412027236 |
| 1 | 0.006803038084 | 1.662887201064 |
| 2 | 0.001727120878 | 1.665707155068 |

Every positive-gap member has0≤s²<1/2. Thus even allowing arbitrary finite
g at fixed positive Δ gives the source-family bound

    gA > 25/18 = 1.388888... .

The family cannot reach the earlier local axial comparison near1.2753.
This excludes this particular one-excitation positive-gap auxiliary as a
complete axial repair. It does not exclude a larger auxiliary Fock carrier,
a different interaction, or an actual relativistic quark source.

## 2. All current complements remain

The physical total-colour Gauss sector has80 states. They decompose into
20 lower coupled states,20 upper coupled states, and40 dark auxiliary
states at energy Δ. Both coupled branches are constructed, verified
orthonormal, and retained. Their spin1/2 and spin3/2 projectors are kept
separately. The full current maps remain on the5500-state domain; no
non-Gauss output is silently projected away.

For the complete contraction instrument K=J/3, the norm bound follows
from J being a sum of three one-quark currents of norm at most1. Every
input has six recorded outcomes: lower spin-half, lower spin-three-halves,
upper spin-half, upper spin-three-halves, dark complement, and current
failure I−K†K. All24 records (three κ, four currents, two initial spins)
sum to one. The dark part is independently checked to have energy Δ.

Exact CAR moments give the axial current norm squared11/3 on a bare
neutron and2 on its normalized auxiliary component. Hence

    ||A D(n)||² = 11/3 − (5/3)s²,
    low-spin-half weight = gA²,
    all complementary weight = 8/9 + (5/27)s² − (25/81)s^4.

Dividing these weights by9 gives the selected instrument probabilities;
the failure probability is16/27+(5/27)s². For the vector current the entire
selected weight1/9 is in the lower proton band and the failure is8/9.
No observed axial ratio is extracted by conditioning away a complement.

The full Hamiltonian checks have residual below2.3e−15; total Gauss and
spin checks are smaller. These are numerical execution checks of exact
source identities. The exact current formula, fractions, and positive-gap
bound do not depend on a floating eigensolver fit.

The local dressing mechanism is now explicit and independently checked.
Its physical axial discrepancy, interacting-vacuum join, native preparation,
Dirac/current transfer, and decay-clock identification remain unresolved.
