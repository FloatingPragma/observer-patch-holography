# Odd elementary response on the actual interacting ground

## Same source, charged intervention and static reference

The source is the full side-two Peter–Weyl construction, with thirty
native fermion modes and five scalar states at each of eight cells, and
twenty-four 119-state links. Its Hamiltonian is
`H_lambda=H0+lambda T`, with all onsite Yukawa, scalar and link
coefficients one and `lambda=2^-34`. Its source clock is the actual
`exp(-it H_lambda)`. The unique neutral physical ground g_lambda and
its controlled atomic cube trial u are those of
[collective_spin_response](../collective_spin_response/PROOF.md).
This calculation does not replace the source by a charged local model
or treat the atomic trial as its exact stationary ground.

At the origin adjoin a full thirty-mode CAR reference in the same native
gauge and spin representation as source creation operators. Its complete
Fock dimension is `2^30`; every occupation is retained. The reference
Hamiltonian is `mu N_ref`, with mu=4. There is no reference spatial
hopping. The full physical constraint at the origin includes the actual
reference gauge current; every other source Gauss constraint is retained.

For a native family r of internal dimension d_r, the map from reference
vacuum to its one-particle family sector is

`B_r=(1/sqrt(2d_r)) sum_(spin,a in r) r^dagger_(spin,a) c_(0,spin,a)`.

These are actual fermion-transfer words, with the global CAR grading.
On reference vacuum, orthogonality of the reference one-particle states
and the CAR give `B_r^dagger B_r=N_r/(2d_r)<=I`. Thus B_r is a
contraction on the complete source carrier. Its Halmos unitary acts on
the direct sum of reference vacuum and the selected one-particle family.
It acts as identity on the complete complementary reference space,
including every other family and all multiple occupations. The failure
map is `sqrt(I-B_r^dagger B_r)` with reference vacuum retained.

The source part of a success has odd fermion parity; the reference has
one fermion, so the complete instrument preserves total parity. Exact
sixty-mode one-particle commutators verify all twelve gauge generators,
all spin generators and raw lepton number for every one of the five
native family probes. Second quantization gives their full-CAR
commutators. The same reference orientation is required for removal;
the dual orientation gives a nonzero hypercharge commutator. The
reference stays quantum. No gauge-variant internal color or weak basis
is promoted to a classical measured label.

The spectral calculations below concern e^c removal and L removal.
They are elementary holes in a filled interacting ground, specified by
their native fields. They are not identified with an observed electron
or neutrino. A quark spectral response is not inferred from the checked
quark intervention algebra.

## Complete positive local spectral measures

The exact local physical ground has a spin-half doublet and the unique
140-dimensional maximal quark module. The native e^c and L operators
act on the lepton factor and therefore retain that quark module. The
complete charged target spaces and actual CAR intertwiners are those of
[charged_ground_exchange](../charged_ground_exchange/PROOF.md). For
removed spin up/down their dimensions are 17 and 7 for e^c, and 33 and
12 for L. All target images are retained.

Let e be the actual local ground energy, selected as the least root of
`[e(e-1)]^3-313[e(e-1)]^2+31708e(e-1)-1041468`. The rank-one
magnetic-up ground projector is evaluated exactly over Q(e) from the
square-free source polynomial p, using

`P_e=sum_i [sum_(j>i) p_j e^(j-i-1)] H_local^i/p'(e)`.

For each complete target generator H_s and actual removal map C_s,
the positive moments are

`m_k=sum_s tr(P_e C_s^dagger H_s^k C_s)`.

The target Gram matrices enter the adjoints. Both spin contributions
are added; the spin-difference functional used for exchange is not a
positive spectral measure. The exact zeroth moment is one. The spin
sum is a scalar, so the pure magnetic-up and mixed ground-doublet
references give the same measure. Gauge covariance gives the same
measure for every internal coordinate of a family. After the full pump
normalization, its atom weights are one-half of these positive weights,
independently of d_r. Their exact sum is 1/2, and failure is 1/2.

The complete target annihilators are

`q_ec(z)=(z-1)(z^4-2z^3-212z^2+213z+11114)`
`                 *(z^4-2z^3-200z^2+201z+9858)`

and

`q_L(z)=(z-1)(z^6-3z^5-303z^4+611z^3+30631z^2-30937z-1033304)`
` *(z^10-5z^9-490z^8+1990z^7+95682z^6-294032z^5-9310567z^4`
`   +19113519z^3+451524798z^2-461130896z-8729229568)`.

The moment sequence determines the proper resolvent n(z)/q(z).
Every atom at a root a of q has normalized pump weight
`n(a)/(2q'(a))`. Rational root isolation, algebraic-field coefficients
and outward interval arithmetic verify strict positivity at every one
of the nine e^c and seventeen L energies. The local complement has
probability zero and is retained. The lowest local energy differences
from the neutral ground and their pump weights are approximately

| Removal field | Lowest visible source-energy difference | Pump weight |
| --- | ---: | ---: |
| e^c | 0.691656680938395 | 0.419951524978232 |
| L | 0.096888497892281 | 0.321280694969404 |

The remaining positive weights sum to approximately 0.0800484750 and
0.1787193050. Their complete individual records, including very small
high-energy weights, are in the receipt.

## Comparison with the complete static-reference sector

Visible response energies need not equal a sector's ground. The code
therefore enumerates every raw lepton occupation and spin in the
maximal quark module at the compensating gauge charge. It imposes
both the actual weak and spin highest-weight kernels, checks the
independent character dimension subtraction, and computes every full
block characteristic polynomial. There are eight nonempty e^c blocks
and nine L blocks. No raw lepton charge is fixed for this ground comparison.

The reference e^c and L irreps are color singlets. Every other
color-admissible quark module has Casimir at most 42, from the complete
occupation inventory, irrespective of its weak/hypercharge label.
The existing full operator bound gives local energy greater than
`-1113/100` on all such modules, with every scalar and lepton state
included. Comparing that bound with the complete maximal-module blocks
proves positive compensated-local gaps above the actual neutral ground:
at least 0.1005852942 for e^c and 0.0968884978 for L. The L minimum
lies strictly below -11.13 and is therefore the exact full local
compensated ground. The exact e^c sector ground is not identified by
this exclusion bound, but it cannot lie below the neutral ground.

This local comparison extends to all link sectors of the lattice.
Write `l0=-1171/100`, `delta=E0-l0`, and use the proved unrestricted
local bound H_local>l0. For m nonblank links, at most 2m cells are
incident. If the reference origin is incident, every incident cell can
be bounded by l0 and every other cell by E0, giving
`H0>=8E0+m(1-2delta)`. If the origin is not incident, its stronger
compensated-local gap is added to that bound. For blank links its
compensated-local gap applies directly. Thus the complete one-reference
sector at zero reference energy has an atomic gap at least

`min(delta_reference,1-2delta)`

above 8E0. The receipt conservatively uses the frozen atomic gap g,
which is no larger than `1-2delta`. The full hopping norm is
`v=lambda*360`. Min-max comparison with the same neutral H_lambda
ground gives a one-reference-sector gap at least
`min(delta_reference,g)-2v`, exceeding 0.0411705465 for both families.
This includes every link set and every charged source state.

A full Fock reference with mu=0 has further occupation sectors. In
particular the filled reference is a gauge singlet because every native
generator has zero trace. It gives the same source-ground energy as
reference vacuum, so uniqueness over the entire reference Fock space
fails. A lower global energy in other zero-cost reference occupations
is not excluded here. This distinction is retained explicitly.

For the declared mu=4, all nonvacuum reference occupations obey the
unrestricted full-source estimate

`H_lambda+4N_ref >= 8 l0-v+4`.

The neutral source ground obeys `E_lambda<=8E0+v`. Their difference
is at least `4+8(l0-E0)-2v>0.1646823120`. Consequently the actual
neutral source ground tensor reference vacuum is the unique physical
ground of the entire declared source/reference Hamiltonian. Every
reference occupation is covered. The reference cost is a declared
apparatus energy, not a derived charged-particle mass.

## Finite full-H spectral windows

The atomic cube trial u has each local ground doublet maximally mixed.
After B_r it decomposes into orthogonal vectors psi_j of squared norm
w_j, at exact total energies
`a_j=7E0+epsilon_j+mu`. The actual ground g_lambda obeys the full
carrier bound `||g_lambda-u||<=eta_g` after phase alignment, so the
pumped vectors differ by at most eta_g.

Choose a rational energy bin I of radius 1/1024 about the center of a
certified rational isolator for the lowest a_j. The exact rational bin
endpoints, rather than a rounded substitute measurement, are recorded.
Let d_j be the distance of a_j to the complement of I for the selected
atom, and its distance to I for each other atom. All d_j have strictly
positive rational lower bounds from the complete local support.

Let P_I be the actual full-H spectral projector. Since
`||(H_lambda+H_ref-a_j)psi_j||<=v||psi_j||`, the spectral theorem gives
selected-atom leakage at most `v||psi_j||/d_j`, and leakage into I from
every other atom at most the same expression. Summing and using
`||psi_j||<=1` gives the amplitude bound

`eta=eta_g+v sum_j 1/d_j`.

This argument uses the complete full source Hamiltonian. It does not
assume isolation from charged states outside the local probe module,
retain only a truncated charged block, or assert unchanged rank of
P_I. Additional full-source eigenstates inside I are allowed. Both the
selected-bin amplitude and its complete outside amplitude differ from
the atomic comparison by at most eta. For an atomic probability in
[a,b] this gives

`[max(0,sqrt(a)-eta)^2,min(1,sqrt(b)+eta)^2]`.

Total pump success and failure use eta_g and their atomic probability
1/2. The complete conditional records at `lambda=2^-34` have the
following approximate bounds:

| Field | Selected bin | Complete outside | Pump failure |
| --- | --- | --- | --- |
| e^c | [0.4197041521,0.4201989707] | [0.0799404942,0.0801565287] | [0.4997605750,0.5002394823] |
| L | [0.3210643365,0.3214971263] | [0.1785579463,0.1788807366] | [0.4997605750,0.5002394823] |

The exact probabilities sum to one. Their interval endpoints are
marginal enclosures, not independent free parameters. The corresponding
source excitation bins, after explicitly subtracting reference energy
four, are approximately [0.6906800974,0.6926332644] and
[0.0959119144,0.0978650814]. Total source/reference excitation bins are
these intervals plus four. The source's absolute atomic energy and
scalar exchange shift are not silently discarded; the comparison uses
`E_lambda in [8E0-v,8E0+v]`.

The spectral bin, its complement and reference occupation commute with
the actual clock. Their probabilities therefore hold at every real
source time after the intervention. Every successful source is odd,
including the entire outside outcome. This proves a finite interacting
odd response, without identifying an isolated charged pole, momentum
dispersion, rest mass, empirical species or continuum limit.

## Preparation and all retained failures

The exact I/D filter comparison is explicitly kinematic pre-physical
preparation. Its successful Gauss/ground filter has probability 1/D.
Its initial failure has probability 1-1/D and generally does not
satisfy source-only Gauss. It is retained as such. Together with pump
failure, selected bin and complete outside, there are four positive
coarse records, whose exact probabilities sum to one. The full static
reference Fock space is retained; unpopulated occupations have zero
probability for this specified preparation.

The physical construction uses
[interacting_ground_preparation](../interacting_ground_preparation/PROOF.md).
It starts from physical source vacua, applies the gauge-preserving local
and cube controllers, and retains every Cayley-filter ancilla outcome.
All its success and failure branches obey source Gauss. Its conditional
source is within trace distance `2.355922487e-7` of the actual ground.
Expanding each detector probability interval by that error gives the
receipt's constructive-preparation bounds. The second source and all
preparation references are retained. The charged detector preserves
joint source/static-reference Gauss on every branch. Neither an
approximate stationary state nor a source-only physical kinematic
failure branch is asserted.

The preparation succeeds with probability at least 0.999999942661984.
Its remaining probability is retained failure. These constructive
records use explicitly declared finite algebraic synthesis and Cayley
services; no full-source efficient native circuit execution is claimed.
