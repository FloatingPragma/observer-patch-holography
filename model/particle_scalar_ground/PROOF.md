# The complete finite scalar ground and its records

The Hamiltonian is the literal finite quantum-scalar and canonical
anticommutation-relation source of the
[interaction construction](../particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md):

\[
 H=W(Y)+W(Y)^\dagger+N_H,\qquad Y=\operatorname{diag}(1,2,3).
 \tag{1}
\]

The scalar has the vacuum and two occupied states. Its occupation is at
most one, so the carrier has dimension `3*2^18=786432`. The fermion-number
energy offset is zero, the scalar pair cost is one, and time is the actual
parameter in `exp(-itH)`. These are declared action coefficients, without
measured particle inputs. No spatial kinetic or link term is added to (1).

The [307-state calculation](../particle_observer_outputs/HIGGS_RECORDS.md)
gives the entire `Q*=N_F-2N_H=0` block, with ground energy −7. The
[ambient variational control](../particle_observer_outputs/BARE_EFFECTIVE_BOUNDARY.md)
proves that this is not the ground of (1). The argument below covers every
charge sector of the original carrier. It reuses the integer-witness
admission principle of the existing
[full-carrier ground calculation](../interacting_vacuum/global_ground.py)
and [ground-gap certificate](../centered_source_ground/gap.py): numerical
proposals carry no proof assumption.

## A complete cell decomposition

Writing scalar-empty and scalar-occupied inputs separately gives

\[
 H=\begin{pmatrix}0&A^\dagger\\A&I\end{pmatrix},\qquad
 A=\begin{pmatrix}P_1^\dagger\\P_2^\dagger\end{pmatrix},
 \qquad M=A^\dagger A.
 \tag{2}
\]

Here `P_h^dagger` are the actual two pair-creation components of W on the
full fermionic carrier. Every positive eigenvalue λ of M gives the two
energies `(1 ± sqrt(1+4lambda))/2`. Empty-scalar null directions have
energy zero and occupied-scalar directions in `ker A^dagger` have energy
one. Thus the full ground is determined by the largest eigenvalue of M.
This is a singular-value decomposition of the complete operator, not a
truncation to low fermion number.

For diagonal Y the modes form six disjoint three-mode cells. Each contains
one right-handed conjugate mode e and its two weak partners L1,L2. The
pair spin sign is absorbed into a local basis phase. Its conserved integer
`k=n_L1+n_L2-n_e` has four sectors, of dimensions `1,3,3,1`. In a signed
occupation basis the two pair maps have the following exact forms:

| Cell sector | Basis interpretation | Pair component a |
|---|---|---|
| k=−1 | e alone | zero |
| k=0 | vacuum, L1e, L2e | `|a><0|` |
| k=1 | filled cell, signed L2, signed L1 | `|0><a|` |
| k=2 | L1L2 | zero |

The source weights are `(1,1,2,2,3,3)`. Even CAR pair operators on disjoint
cells commute; the signed fermionic mode permutation supplies their exact
tensor embedding in the original occupation ordering. `local_frame` and
`original_column` specify every sign and original mode address. The
dimension sum is `(1+3+3+1)^6=2^18`, including all inactive cells.

An active cell has either of the three-dimensional representations in the
middle rows. The largest active matrix therefore has dimension `3^6=729`.
Every coefficient in its two pair maps is nonnegative, so every entry of
M is nonnegative. Permuting cells of equal weight gives exact unitaries:
the 64 active patterns reduce to 27 count patterns, without changing
multiplicity. Each M further splits into the connected components of its
actual nonzero-entry graph.

Inactive cells cannot attain the global maximum. Replace one by the
fundamental three-state cell and embed its original input in `|0>`. The
compression of the enlarged M is exactly `M_old+2w^2 I`: cross terms have
zero vacuum expectation. Thus the enlarged active sector has top
eigenvalue at least `lambda_old+2w^2`. Repeating this for every inactive
cell bounds its top eigenvalue by the full active maximum minus at least
two. This comparison covers all omitted occupations; it does not discard
their states or assign them zero detector weight.

## Exact global comparison

For each connected integer matrix C, a numerical routine proposes a
strictly positive integer vector v. Acceptance uses Python-integer
multiplication and the rational inequalities

\[
 \min_i (Cv)_i/v_i\leq\lambda_{\max}(C)
     \leq\max_i (Cv)_i/v_i.                         \tag{3}
\]

They follow from the nonnegative-matrix Perron theorem, or directly by
diagonal conjugation and the minimum/maximum row-sum bound. Neither an
eigensolver residual nor a numerical success flag is admitted. Every
connected component of every active count pattern is included.

The largest component has 53 original scalar-empty occupation states.
Its cell types are dual, dual, fundamental, fundamental, fundamental,
fundamental. Its top eigenvalue is the largest root of

\[
 p(x)=x^6-244x^5+23737x^4-1171722x^3
       +30638880x^2-394771968x+1881169920.            \tag{4}
\]

The code calculates the characteristic polynomial of that actual integer
matrix, verifies divisibility by (4), and uses exact Sturm root counts
to isolate the root in its rational interval. It also verifies that every
other root in the same component is below `703/10`; (3) places every
competing component below the same bound. The inactive-cell bound lies
below it too. The winning count pattern has no multiplicity from equal
weight permutations, and its graph is connected. Perron simplicity and
the strict comparison prove a unique full-carrier ground.

The admitted intervals imply

\[
 \lambda_{\max}\simeq71.09007333321,\qquad
 E_0={1-\sqrt{1+4\lambda_{\max}}\over2}
        \simeq-7.94630530665.
 \tag{5}
\]

The complete Hamiltonian has an excitation gap greater than `0.04690056`.
The bound follows monotonically from the two square roots at the certified
top-root lower bound and the competing upper bound `703/10`. Dark energies
zero and one are included and lie higher. Equation (5) is in the source's
energy units; it is not a particle rest mass.

## Gauge, spin and preparation

Every original occupation in the winning empty-scalar support has
`N_L=4,N_e=2,N_F=6`. Its scalar image has `N_L=5,N_e=3,N_H=1`.
Both therefore have `Q*=6` and hypercharge
`(-N_L+2N_e-N_H)/2=0`. The full source Ward identities commute H with
weak gauge SU(2) and the declared spin SU(2). A unique ground transforms
in a one-dimensional representation of either connected SU(2), which is
trivial. The ground is a weak-gauge and spin singlet. Color acts trivially
on this eighteen-mode lepton source. These are the actual finite source
actions; they do not supply relativistic boosts or a continuum vacuum.

Let u be the normalized positive eigenvector of the 53-state M, including
its signed embedding into the original CAR basis, and let
`v=Au/sqrt(lambda_max)`. The exact ground is

\[
 |G\rangle=\sqrt\alpha\,|0,u\rangle
             -\sqrt\beta\,|1,v\rangle,
 \quad \Delta=\sqrt{1+4\lambda_{\max}},
 \quad \alpha={1+1/\Delta\over2},\quad\beta=1-\alpha .
 \tag{6}
\]

The notation `|1,v>` includes both occupied scalar components and their
entire fermionic image. All coefficients are computable algebraic numbers.
The elementary seed with all six first-family fermion modes filled and
the scalar empty has the same conserved charge and is also a gauge/spin
singlet. A Householder reflection from this seed to (6), identity on the
complete orthogonal complement, is a full unitary preparation. Because
the seed and target are singlets, that preparation commutes with the
actual gauge and spin actions. It also preserves Q*.

The filled-family seed has mean H energy zero. Preparing (6) changes that
mean by E0, which is negative; the declared controller receives the released
energy. The Householder is a preparation service, not autonomous evolution
under H or an energy-conserving unitary on the isolated source. Its
controller/work degrees of freedom belong to any closed implementation.

The executed rational witness projects its proposed u onto the exact weak
and spin singlets before normalization. On `Jz=0`, the finite Casimir
polynomial `product_{j=1}^{jmax}(I-J^2/[j(j+1)])` is the singlet projector;
`jmax=2` for the four left weak modes and `jmax=3` for six spin modes.
The code applies the literal CAR raising/lowering operators and verifies
both raising operators annihilate the resulting vector. Its actual A
image gives a trial Rayleigh energy. The complete ground gap converts its
energy excess into an infidelity bound below `2.6*10^-11`. Hence the pure
preparation trace distance is below `5.1*10^-6`.

The exact or approximate full preparation can be a declared coherent
service, or can use the existing
[finite-instrument compiler](../spin_operational_position/PROOF.md) under
its operational premises. No large primitive tape is asserted to have
been executed. The original pair Hamiltonian preserves Q*, so evolution
from the empty Q*=0 state alone cannot prepare this Q*=6 ground. A choice
of preparation sector or a charge-changing source operation must be part
of that program; uniqueness of the architecture is not required.

## Complete records on the same Hamiltonian and clock

Use the complete scalar apparatus `P0=|0><0| tensor I` and
`P1=I-P0`, with ranks `2^18` and `2*2^18`. Every original carrier input
and external quantum reference is retained. On (6), the two singular
directions form an exact invariant pair, so the ground scalar expectation
is β and the connected ordered correlator is

\[
 \langle N_H(t)N_H(0)\rangle_G-\beta^2
       ={\lambda_{\max}\over\Delta^2}e^{-i\Delta t}. \tag{7}
\]

An actual first occupation measurement has backaction. Its full two-time
joint law is

\[
 J_{00}=\alpha(1-T),\quad J_{01}=\alpha T,\quad
 J_{10}=\beta T,\quad J_{11}=\beta(1-T),\qquad
 T={4\lambda_{\max}\over\Delta^2}\sin^2(\Delta t/2).
 \tag{8}
\]

All four outcomes sum to one. Other carrier states are retained by the
same projectors, even though their input weights in this ground experiment
are zero. The code encloses all four original probabilities at the fixed
ticks `t=1/32` and `2t=1/16` using rational sine bounds. The same actual
Hamiltonian step is repeated in the second experiment. Their unconditional
ratio satisfies `J01(2t)/(4J01(t))=cos^2(Delta t/2)`; the source bounds
certify `0<Delta t<pi` and positive denominator, so the measured ratio
determines the finite gap on that branch. No accepted-only histogram or
separate clock is substituted.

This occupation transition has `Delta` about 16.89261. It differs from
the complete ground-isolation gap bounded below by 0.04690056: this
particular apparatus couples the two singular directions in (6), and has
zero residue on the lower-energy transitions outside that pair.

For a preparation infidelity η and each declared compiler operator-norm
error budget, the complete two-tick/two-read law differs in total variation
by at most `sqrt(eta)+epsilon_prep+2epsilon_tick+2epsilon_read`. The same
bound holds with untouched references. This accounts for preparation and
every read, including failure labels in a compiled realization. Agreement
between observers can copy and interpret these orthogonal records through
the existing record instrument. It does not turn (7) into a relativistic
mass or determine a continuum geometry.

## General family coefficients and the coherent source

For any complex Y, take `Y=U Sigma V^dagger`. Replacing the left creation
modes by U and the right creation modes by the complex conjugate of V
turns the full pair operator into the diagonal singular coefficients.
This is an exact one-particle unitary and its full exterior-algebra lift,
not a projection. It commutes with the family-independent gauge/spin
currents and preserves the scalar occupation, reference domain and clock.
Equations (2) and the complete cell-sector reduction therefore hold for
arbitrary positive singular coefficients, repeated twice for the two
spins. For noninteger coefficients the positive-matrix bounds (3) accept
rigorous coefficient intervals in place of integer entries.

The actual coherent source at unit phase `z=1 (theta=0)` and pass amplitude
one has two equal singular coefficients and a third larger by
`34+24sqrt(2)`. Its family-unitary transport is covered by this theorem.
The specific root (4), uniqueness comparison, ground charge and numerical
record tables certify only the rational fixture in (1). They are not
assigned to that different coherent action. That action requires its own
weight-dependent ground comparison. No fitted particle masses or phase
selection is inferred from either form of the finite block theorem.
