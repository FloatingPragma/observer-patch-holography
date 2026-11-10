# Native quark removal on the interacting local reference

The Hamiltonian and scalar register are the full thirty-mode native local
source in `chiral_joint/local_reference`. The local reference is its actual
physical ground projector divided by two. Its quark factor is the unique
140-dimensional collective Spin(5) module with Casimir 54; its lepton factor
and scalar register are retained. The source has raw lepton charge one and
zero baryon charge. The quark factor is a native color and physical-spin
singlet. This calculation does not replace that reference by the empty Fock
state.

## Literal fermion map and normalization

The actual particle–hole automorphism exchanges creation and annihilation
on the native weak singlets. The transformed one-body Yukawa matrices have
six four-mode quark components. Signed permutation intertwiners, checked
against all four actual transition matrices, identify each component with
the representative `(0,12,22,23)`. The exterior-square Casimir-four subspace
has dimension five and positive metric `diag(2,1,1,1,1)`.

Let `v_i` denote that five-state basis and let `n=(n_1,...,n_5)` have total
degree `d`. Define its symmetric tensor embedding by

`J_d z^n = (prod_i n_i!)/sqrt(d!) sum_distinct v_(i_1) tensor ... tensor v_(i_d)`.

Its Gram matrix is `prod_i n_i! g_i^n_i`. The six factors are distinct actual
four-mode Fock components, ordered through a canonical fermion permutation.
Every source factor has even occupation two, so the annihilator on the
first component has no preceding odd-factor sign. The exterior powers of
the checked signed permutations supply the remaining five actual component
frames. Collective transitions intertwine this embedding exactly.

The Casimir-54 subspace of `Sym^6(5)` has dimension 140. Its embedded image
lies in the full native quark Fock space and is invariant under the same
ten generators. The full character calculation in `local_reference` gives
exactly one such irreducible in that Fock space. Thus this construction
fixes the actual source module and the actual CAR map, including its
normalization; it does not infer the probe matrix from matching energies.

Splitting the first factor gives the exact identity

`J_6 p = (1/sqrt(6)) sum_i v_i tensor J_5 partial_i p`.

Consequently an elementary transformed annihilator or creator has image

`D_alpha p = sum_i (c_alpha v_i) tensor partial_i p`, or
`D_alpha p = sum_i (c_alpha^dagger v_i) tensor partial_i p`.

The rational coordinates use target embedding `J_5/sqrt(6)`, so the target
metric is the five-fold symmetric metric tensor the odd-sector identity,
divided by six. The literal four-mode CAR matrices determine all signs.
The code checks the ten generator covariance identities, full image
inclusion and

`sum_alpha D_alpha^dagger D_alpha = 2 I_source`

for both removal and creation. The adjoints in this identity use the exact
positive metrics. The factor one sixth is required by the actual tensor
splitting, even though it does not affect the Hamiltonian eigenvalues.

## Complete invariant image

Native `Q` removal is transformed annihilation and has transformed quark
number eleven. Native `u^c` and `d^c` removal are transformed creation and
have transformed quark number thirteen. The full ambient target spaces are
`Sym^5(5) tensor exterior^1(4)` and
`Sym^5(5) tensor exterior^3(4)`, each of dimension 504.

For each of the four CAR maps, direct multiplication gives

`C_target D_alpha = (95/2) D_alpha`.

At every weak/hypercharge weight the combined image rank equals the full
kernel rank of `C_target-95/2`. The total is exactly 224. Every target
transition and actual weak/hypercharge current preserves this kernel.
Thus the entire Hamiltonian orbit of each elementary image is retained.
No assumption that a representation label implies nonzero probe support
enters this argument.

The actual transformed quark hypercharge is the collective traceless
current plus `Nq'/6-2`. It therefore adds `-1/6` on the annihilation target
and `+1/6` on the creation target. Omitting these terms selects the wrong
native charge space. Representative original removed modes are zero for
`u^c`, twelve for `d^c`, and twenty-three for the lower weak component of
`Q`. Their source charge changes are respectively

| Probe | Quark number after transformation | Highest weak weight | Hypercharge change |
|---|---:|---:|---:|
| `u^c` removal | 13 | 0 | `+2/3` |
| `d^c` removal | 13 | 0 | `-1/3` |
| `Q` removal | 11 | `+1/2` | `-1/6` |

The full native color and quark-spin currents commute with the local
Hamiltonian. The initial quark module is a singlet under both groups.
These exact symmetries relate the color and quark-spin partners of each
probe. The two weak components of `Q` are related by the actual weak
current. Lepton spin commutes with each quark probe and the Hamiltonian,
so one highest lepton-spin component supplies the same measure as the
mixed physical ground doublet. These multiplicities do not remove the
static reference needed to make an individual colored probe physical.

## Scalar and lepton factors

The target retains all 64 lepton Fock states and all five scalar states
before imposing the conserved transformed lepton number three and
lepton spin one half. The complete weak highest-weight kernel is imposed
at the actual target hypercharge. Kernel dimensions are checked by the
independent weight subtraction formula. The scalar excitation projector,
all quark transitions and all lepton transitions retain their actual
coefficients. Every Hamiltonian image and every elementary source image
is checked in the full weight space. The Gram identity `G H=H^dagger G`
holds before and after the kernels are imposed.

The invariant image theorem concerns the complete spectral support of
these elementary probes. It does not identify the ground of every
colored compensated sector, nor does it establish propagation or a
physical quark pole in the interacting many-cell source.

## Exact positive energy measures

Let `p` be the square-free annihilator of the complete sixteen-state source
magnetic block and let `e` be its actual least ground root. Its rank-one
spectral projector is evaluated as

`P_e = [p(H)/(H-e)]/p'(e)`.

Polynomial division occurs before matrix evaluation. The source trace of
this polynomial is checked to equal one. For each literal probe matrix
`A`, all target moments

`m_j = tr(P_e A^dagger H_target^j A)`

are computed in the exact number field of the degree-six neutral ground
root. The complete target annihilator `q` is checked directly on the full
target matrix. The rational resolvent numerator is reconstructed from its
moments. At every real root `epsilon` of `q`, the spectral weight is the
residue `numerator(epsilon)/q'(epsilon)`. Exact rational root intervals and
outward number-field arithmetic certify every visible weight strictly
positive. The exact numerator/denominator gcd identifies every zero
residue in the complete target spectrum; these roots are retained as
zero-probability records. No numerical eigenvector or numerical eigensolver
is used.

Each of `u^c` and `d^c` has thirteen lines. Their complete rational
resolvents are exactly equal in this source with all Yukawa coefficients
one. This computed degeneracy is not a physical up/down mass equality.
The `Q` response uses its separate twenty-five-root target annihilator;
its exact resolvent cancellation leaves thirteen visible roots and twelve
exact zero residues. Each
measure has exact total weight one half. The equality includes all dark
energy-one states through their aggregate spectral projector.

The gauge-covariant probe is

`B_r = sum_(spin,a in r) r_(spin,a)^dagger c_(origin,spin,a)/sqrt(2 d_r)`

on the vacuum of the full thirty-mode static reference. This reference
transforms in the same native representation as source creation operators;
the opposite orientation fails the hypercharge commutator. The existing
full sixty-mode operator certificate checks all twelve Gauss, all three
spin and raw lepton commutators for all five native families. The joint
probe is even, and the reference remains an unmeasured quantum system.
Its effect is `N_r/(2d_r) <= I`. A Halmos dilation supplies the failure
operator and acts identically on the complete orthogonal reference
complement. Every branch therefore retains the Gauss constraints.

The individual color, spin and weak components have equal measures by
the exact symmetries stated above. The normalized sum over orthogonal
reference components consequently has the same measure as the computed
one-component probe. There is no additional factor one half. Pump success
and failure each have probability one half in the actual local ground.
Outside the complete local support the probability is zero. The static
reference Hamiltonian is `4 N_ref`, so every successful line has reference
energy four and the failed branch has reference energy zero.

## Full eight-cell source and actual preparation

Let `h=H0+2^-34 T` be the unchanged eight-cell native Hamiltonian, including
all scalar, fermion and PW119 link states. Its actual unique physical
ground is `g`; its frozen atomic cube trial is `u`. The inherited bounds
are `||g-u|| <= eta` and `||2^-34 T|| <= v=360*2^-34`. Each local quark
spectral component of `B_r u` is an exact atomic vector `v_j` with

`H0 v_j=(7e+epsilon_j)v_j`, `||v_j||^2=w_j`.

This holds in the spin-entangled cube trial because the quark probe and
Hamiltonian commute with the separate lepton-spin action. The associated
local moment operator is scalar on the complete local ground doublet.
Thus no two-cell source or empty Fock reference is substituted.

The retained Hamiltonian for preparation and observation is

`H_total = h_active + h_partner + 4 N_ref`.

The actual physical two-source preparation, including all its agreement
filter records, is inherited from `interacting_ground_preparation`. Its
success probability is `q` and its conditional active-source trace distance
from `|g><g|` is at most `delta`. No exact ground projection is added. The
prepared conditional state is not asserted stationary. Its distance from
the stationary ground is preserved under every common autonomous wait.

Choose disjoint bins of radius `128v` around each isolated atomic total
energy `7e+epsilon_j+4`, with rational intervals accounting for the root
enclosures. Let `P_j` be the exact spectral projector of
`h_active+4N_ref` onto its bin. If `d_j` is the distance of its atomic root
interval to the bin complement, then

`||(I-P_j)v_j|| <= v sqrt(w_j)/d_j`.

For another atomic line `k`, with interval-to-bin distance `d_jk`,

`||P_j v_k|| <= v sqrt(w_k)/d_jk`.

Orthogonality of `P_j` and `I-P_j` gives

`||P_j v_j|| >= sqrt(w_j(1-(v/d_j)^2))`.

The lower bin amplitude is this last expression minus all other-line
leakages and `eta`, clipped at zero. Its square is the recorded lower
probability. The upper amplitude is `sqrt(w_j)` plus all other-line
leakages and `eta`. Every atomic line enters these sums, including tiny
residues. The interval arithmetic rounds outward before the complete
finite sum. Exact-ground pump success remains one half by the full native
antiunitary occupation identity from `interacting_charge_conjugation`.
This source-only symmetry preserves `h` and reverses the evolution
parameter. It fixes the origin and sends its occupation to its complement;
uniqueness of the physical source ground therefore gives occupation one
half in each mode. No charge-conjugation invariance of the retained
reference Hamiltonian `4N_ref` is asserted. The probability outside all listed bins
is bounded using their orthogonality and this exact total.

For the actual prepared source, each bin and the single outside-all-bins
event change by at most `delta`. Each unconditional probability lower
bound multiplies the prepared lower endpoint by the inherited lower
bound for `q`. Tiny positive local residues need not retain a strictly
positive eight-cell lower endpoint under these errors; the receipt
retains their intervals and the full outside record without a positivity
claim for them.

The complete law consists of preparation failure, pump failure, every
listed energy bin and the outside-all-bins event:

`1-q`, `q tr(K_fail rho K_fail^dagger)`,
`q tr(P_j B_r rho B_r^dagger)` for every `j`, and
`q tr((I-sum_j P_j) B_r rho B_r^dagger)`.

These events sum exactly to one. The apparatus retains the static
reference and all native preparation references and pointers. The energy
instrument is a declared exact full-source spectral service with
`degree(q)+1` energy labels. It is a finite service on the full carrier;
no efficient synthesis or execution of its enormous source matrix is
claimed. Pump and energy services each last one source unit. For a
unitary completion `U` and Hermitian logarithm `L` with `exp(-iL)=U`, the
interaction-picture control
`V(t)=exp(-it H_total) L exp(it H_total)` gives full evolution
`exp(-i H_total) U`. Its pointer records the interval-start observable.
Common autonomous evolution continues on every failure and no-operation
branch. Energy projectors commute with that evolution, so the displayed
law includes the pump/readout duration without a long-time perturbation
comparison. Bare-source work per completed service is bounded by
`2||H_total||`; controller workspace is a separate implementation resource.

The source-only excitation bins subtract the explicit reference cost
four and the actual neutral ground interval, including its certified
second-order atomic energy shift. The bins are source-clock energies,
not rest masses. Static color compensation at the origin does not supply
asymptotic colored states, propagation, confinement or Standard Model
matching.
