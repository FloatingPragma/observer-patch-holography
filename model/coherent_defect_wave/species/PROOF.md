# Shared-link CAR extension and species boundary

## Domain and declaration

Let `L` be the ten link words of the parent's fixed-boundary one-regional-defect
sector, and let

\[
 h=\mathbb C^{12}_{\rm port}\otimes\mathbb C^5_{\rm internal}
       \otimes\mathbb C^2_{\rm repair}\otimes\mathbb C^{60}_{\rm address}.
\]

The parent clock is `U=C4 C3 C2 C1 C0 B`, where
`B=sum_g |g><g| tensor W_g`, `W_g=Q_g R`, and `C_e` applies the involution
`T_e` on the same link word precisely when the wave port is an endpoint of e.
The [parent proof](../PROOF.md) establishes the unitary factors, fixed boundary
and central gauge constraint.

Declare the full wave Fock algebra and Hilbert space

\[
 \mathcal F(h)=\bigoplus_{r=0}^{7200}\Lambda^r h,
 \qquad \mathcal H_F=\mathbb C^{10}_{\rm links}\otimes\mathcal F(h).
\]

This is one link register with several wave occupations. Taking the exterior
algebra of the parent's entire 72000-dimensional joint space instead would
duplicate the links; that different construction is not used. Sparse fixed-r
executions do not truncate the algebra. A literal truncation at r=2 would fail
full CAR at its top occupation boundary.

Let `N_e` count all occupied wave modes whose port is an endpoint of e.
For a fixed declared Boolean function `f` on nonnegative integers, set

\[
 \begin{split}
 B_F&=\sum_{g\in L}|g\rangle\langle g|\otimes\Gamma(W_g),\\
 C_{e,F}&=T_e\otimes P_{f(N_e)=1}
              +I\otimes P_{f(N_e)=0},\\
 U_F&=C_{4,F}C_{3,F}C_{2,F}C_{1,F}C_{0,F}B_F.
 \end{split}
\]

Here `Gamma(W)=direct_sum_r exterior_power_r(W)`. The implemented primary law
is `f(n)=n mod 2`; the explicit comparison is `f(n)=1[n>0]`. Both satisfy
`f(0)=0, f(1)=1`. Choice of f is a construction input. Neither measured masses
nor coupling fits enter it.

## Unitary and source identities

Each `Gamma(W_g)` is unitary because an exterior power preserves inner
products of wedges through their Gram determinant. The link projectors are
mutually orthogonal, hence `B_F` is unitary. Each `C_e,F` is an involution:
its two occupation projectors are orthogonal, `T_e^2=I`, and T acts on links
while N acts on wave modes. Reversing the ordered controls and applying
`Gamma(W_g*)` gives the exact inverse.

On occupation one, `Gamma(W_g)=W_g` and `f(N_e)` is exactly the parent's
endpoint projector. Thus the map

\[
 |g\rangle\otimes|p,j,r,a\rangle
 \longmapsto |g\rangle\otimes a^\dagger_{p,j,r,a}|0\rangle
\]

intertwines U and `U_F` for every state, preparation and one-wave detector.
The executable checks 180 basis inputs across both link-sensitive and
link-blind internal sectors, every internal index, two repair values, several
ports and addresses. These are finite regression witnesses for the universal
identity just proved, not a substitution for its proof.

For ordered occupation sets, left creation inserts a mode with sign
`(-1)^(number of occupied modes below it)` and annihilation removes it with
the same sign. If i=j, exactly one of creation-then-removal or
removal-then-creation acts as the identity. If i differs from j, exchanging
the order changes one insertion/removal sign. Therefore

\[
 \{a_i,a_j^\dagger\}=\delta_{ij}I,
 \qquad \{a_i^\dagger,a_j^\dagger\}=0,
 \qquad \{a_i,a_j\}=0.
\]

The argument includes arbitrary occupied spectators. Four-mode tests execute
all 16 basis states and all 16 pairs of mode indices. The statistics are
supplied by this declared exterior architecture. No rotation action or
spin-statistics theorem is inferred from CAR alone.

## Gauge and charge preservation

Write `J|g>=|g xor star(5)>` and let D flip the sign of wave modes at port5
in internal indices3 and4. The parent's constraint is `G=J tensor D`, with
`W_(g xor star) D=D W_g` and `T_e J=J T_e`. Its lift is

\[
 G_F=J\otimes\Gamma(D),\qquad G_F^2=I.
\]

The exterior functor carries the wave covariance identity to every degree.
D is diagonal in occupation and commutes with every endpoint-number
projector. Every control commutes with `G_F`, so `[U_F,G_F]=0`.
The complete flux-face projector sums both gauge representatives and also
commutes with `G_F`. Each detector branch preserves the positive constraint.

The number N and every total internal occupation `N_j` commute with the
clock. In particular the additive charge
`Q5=2(N0+N1+N2)-3(N3+N4)` is conserved. Every possible link transition is
the parent's T_e, so the fixed boundary bits, exterior syndrome and regional
odd flux are preserved. These identities concern the central finite
constraint. They do not supply full local `SU(3) x SU(2) x U(1)` Gauss laws.

## Exact two-wave observation

Prepare occupied modes `(port5,internal3,repair0,address0)` and
`(port5,internal4,repair0,address0)`. Their D sign is positive. The equal
superposition of link words1 and6553635, divided by square root2, lies in
the positive `G_F` sector. One complete clock step has 16200 nonzero exact
amplitudes and norm1. Its complete flux-face distribution is

| Carrier face index | Parity-controlled probability |
|---|---:|
| 0 | 7417/7500 |
| 2 | 449/202500 |
| 5 | 1/202500 |
| 16 | 112/50625 |
| 17 | 1343/202500 |

The probabilities sum to1. Every branch satisfies the positive constraint;
wave number2, Q5=-6 and regional charge1 hold. The presence-controlled law
has a different complete distribution, including face17 probability
`809977/810000`. Both laws restrict to the identical one-wave U. Therefore
the one-wave source does not select its multi-wave response without a
declared occupation control. The fixed parity law provides one executable
choice, with actual coherent shared-link interaction and exact readout.
This experiment is not a relativistic collision or a two-particle S matrix.

All ten empty-wave vectors `|g> tensor |0>` are fixed by U_F. Their positive
Gauss subspace has dimension5, one vector per gauge-paired face location.
They are matter vacua inside an odd-flux sector. They neither supply a
zero-flux vacuum nor select a unique physical vacuum. No pair creation,
annihilation or wave-number-changing term is present.

## Exact internal reduction

Every W_g and every control preserves j. The link signs are identical on
j=0,1,2 and identical on j=3,4. Thus, after ordering the factors,

\[
 U=I_3\otimes U_t\;\oplus\;I_2\otimes U_d,
 \qquad \dim U_t=\dim U_d=14400.
\]

The constraint has the same decomposition, with positive-sector dimension
7200 per internal copy. Consequently all matrix units within each internal
block commute with U and the corresponding full `U(3) x U(2)` action does
too. The endpoint-number controls preserve this symmetry in the Fock lift.

This is also the exact commutant among pure internal unitary transformations.
For a cross-block matrix unit, the commutator is that matrix unit times
`U_t-U_d`. The receipt exhibits a basis vector whose cross-block commutator
has squared norm `421/225`, so `U_t-U_d` is nonzero. Every cross-block entry
of a commuting internal matrix must vanish. The trace-free anti-Hermitian
part therefore has dimension `9+4-1=12` and is `s(u(3) plus u(2))`.
This identifies an actual global internal symmetry of this source; making
that symmetry local requires link variables and Gauss constraints realizing
its noncentral transformations.

For any spectral function F, including every spectral projector,
`F(U)=I3 tensor F(U_t) direct_sum I2 tensor F(U_d)`. Corresponding
excitations in the three triplet copies have equal spectral outputs. Naming
these copies electron, muon and tau cannot produce unequal corresponding
masses. Different preparations or apparatus can break the internal symmetry;
an arbitrary detector weighting is not a source spectral mass splitting.
Different eigenbranches can have different eigenphases, but their physical
species interpretation requires its own observable proof.

The parent uniform triplet fixed vectors are related by this same internal
symmetry. Selecting one such vector is a reference choice; it does not
generate a three-family hierarchy.

## The exact fifteen-state comparison

The auxiliary five-mode construction `Lambda^2 C5 plus Lambda^4 C5` has
15 states. With `q=6Y=-Q5`, an auxiliary mask containing n3 triplet modes
and n2 doublet modes has `q=-2n3+3n2`. Exterior color and weak weights give

| Auxiliary occupation counts | SU(3), SU(2), Y | Dimension | Standard Model comparison |
|---|---|---:|---|
| (2,0) | (3bar,1,-2/3) | 3 | up singlet conjugate |
| (1,1) | (3,2,1/6) | 6 | quark doublet |
| (0,2) | (1,1,1) | 1 | charged-lepton singlet conjugate |
| (3,1) | (1,2,-1/2) | 2 | lepton doublet |
| (2,2) | (3bar,1,1/3) | 3 | down singlet conjugate |

These are the gauge representations of one all-left-handed Standard Model
fermion family, with no right-handed neutrino. The conventional left doublets
and right singlets are described in the
[PDG electroweak review, section10.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf);
the table expresses right singlets as left-handed conjugates. This is a
representation comparison, with no measured mass input.

The receipt independently computes every individual color Cartan weight,
T3, Y and electric charge. It reproduces the representation in
[`standard_model_source_map`](../../standard_model_source_map/PROOF.md).
The auxiliary masks are neither the parent's five elementary j states nor
automatically physical fermionic creators. Indeed two disjoint degree-two
auxiliary creation products commute: moving two creators past two creators
gives sign `(-1)^4=+1`. A physical one-fermion space carrying these 15
internal states can be declared, as in the existing
[fermion attachment](../../standard_model_source_map/FERMION_ATTACHMENT.md),
but its dynamics and transfer from the present wave register are separate.

## Finite-theorem premise crosswalk

| Existing result | Reusable conclusion | What an application to this source requires |
|---|---|---|
| `Lean/Screen/A2HolonomyBridge.lean`, `compact_product_dimensions_of_projective_holonomy` | Conditional compact response-group dimensions from projective holonomy, fixed-space and classification data | Actual response-generated continuous group, full holonomy and the theorem's fixed-space/classification witnesses; central signs alone do not supply them |
| `Lean/Screen/GaugeSectorBracketCompletion.lean`, `mirror_pair_invariance_one_ray` | A kinetic coefficient ray for its declared field tuple and mixing generators | The same dynamical generators and kinetic strengths must be extracted from U_F; internal spectator symmetry alone is insufficient |
| `Lean/Screen/FermionSectorAssembly.lean`, `fermionSectorAssembly_receipt` | Gauge-charge-compatible finite algebraic pairings on its declared coefficient bundle | The module explicitly has no derivative, physical Spin action or covariant fermion kinetic operator; those must act on the joined source |
| `Lean/Geometry/CanonicalLorentzModule.lean`, `det_toMatrix`, `finrank_Herm2` | Hermitian2x2 matrices form a four-dimensional Lorentz quadratic module | A map from the actual clock, spatial readouts and correlations to that module, compatible with its transformations |
| `Lean/Geometry/MassShellKinematics.lean`, `lorentzQ_fourMomentum` | The declared mass and timelike frame obey the mass shell | A measured source energy/momentum branch, positive energy and a unit convention; a declared mass parameter cannot calculate its value |
| `Lean/QFT/FiniteUnitaryScatteringNoGo.lean`, `finite_unitary_powers_have_no_limit` | A nontrivial finite unitary power sequence cannot converge | Define finite-time scattering records, or provide a controlled large-domain limit; a finite recurrence is not an asymptotic particle S matrix |
| `Lean/QFT/InheritanceMatrix.lean` | Exact typed boundaries for a particular finite-net construction | A source-specific symmetry action, vacuum, sector and scattering carriers; the matrix gives neither an automatic inheritance theorem nor a universal obstruction |
| `model/standard_model_source_map/FERMION_ATTACHMENT.md` | Exterior lifting, CAR and transferred internal currents on a declared fermion one-particle space | The present shared-link lift supplies CAR for h; identifying h with the physical 15-state chiral field space requires an explicit intertwiner |

The related `repair-law/work/SIM/CONSOLIDATED.md` occupation enumeration and
`repair-law/audit/gauge/REPORT.md` gauge checks concern declared charge
representations and their own response laws. Their algebraic comparisons do
not transfer physical species labels or continuous gauge dynamics to U_F.
The central code action `diag(I3,-I2)` and a global commutant are compatible
with many continuous gauge extensions. They do not determine which extension
has physical local gauge fields.

The concrete source-to-particle requirements are: a common spatial clock and
reference carrying controlled dispersion; spatial rotation and chirality
actions on its excitations; observable local gauge currents and interactions;
a controlled many-wave scattering/readout regime; and rest-energy extraction
with physical units. These are mathematical attachment conditions on the
specified source. Architectural uniqueness is not required.

## Spatial Fock circuit on the same local links

For a finite declared cell set X, use one local link register per cell and
one global wave Fock algebra:

\[
 \mathcal H_X=(\mathbb C^{10})^{\otimes X}
                 \otimes\mathcal F\left(\bigoplus_{x\in X}h_x\right).
\]

Every local `U_F,x` is the occupation extension proved above, acting on its
own link and wave modes. On global occupation one, an occupied cell executes
the literal parent U and the empty cells execute the identity. This preserves
the source identity of the onsite law without making each cell permanently
contain one wave.

Choose the canonical representative `min(g,g xor star)` in each link orbit.
Let R_g act on a wave mode with sign -1 precisely when g is the other
representative and the mode has port 5 with internal index 3 or 4. For two
neighboring cells x,y define the one-body dressed swap

\[
 S_{g_x,g_y}= (R_{g_x}\oplus R_{g_y})\,S_{xy}\,
                    (R_{g_x}\oplus R_{g_y}),\quad
 V_{g_x,g_y}=aI+bS_{g_x,g_y},\quad
 a=(1-i)/2,\quad b=(1+i)/2.
\]

Here `S_xy` swaps the two one-body summands h_x and h_y. It does not move
either link register. Since the dressed swap is a self-adjoint involution,
V is unitary, has V squared equal to S, and has the conjugate-coefficient
inverse. Define the actual many-wave hopping gate by exterior quantization,

\[
 H_{xy}=\sum_{g_x,g_y}|g_x,g_y\rangle\langle g_x,g_y|
                      \otimes\Gamma(V_{g_x,g_y}).
\]

Under a gauge flip at x, `R_(g_x xor star)=D_x R_gx`. This gives
`V_(g_x xor star,g_y) D_x = D_x V_(g_x,g_y)`. Applying Gamma proves that
H commutes with `J_x tensor Gamma(D_x)`; the same holds at y. Both local
links are retained exactly. Each gate preserves total wave number and Q5.

Declare a clock consisting of one onsite layer and six disjoint nearest-
neighbor matching layers, the two parities in each of three cubic axes.
All gates are even in the CAR algebra. Gates on disjoint cells commute;
odd fields on disjoint cells anticommute, while even observables commute.
Conjugating an observable by one layer expands its support at most across
that layer's edges. After t ticks its coordinate radius is at most 2t cells
and its Manhattan radius at most 6t. These are exact finite circuit bounds.
They do not prove a Lorentz-invariant continuum limit.

The empty-wave vector is fixed by every Gamma(V), and the onsite clocks are
identity at empty occupation. Tensor products of positive-Gauss link-pair
states with this wave vacuum therefore give explicit common stationary
references. They have odd regional flux in each cell and are not selected
physical vacua. The two-cell interpreter executes these statements without
allocating the full Fock space.

A constrained preparation with one distinct mode in each cell produces
complete receiver occupation probabilities `(1/4,1/2,1/4)` for counts
`(0,1,2)`. All branches preserve both Gauss constraints. The code executes
the inverse exactly, verifies fixed local links during hopping, and evolves
a one-wave hopping output through the literal receiver U with its inverse.

The gate is distinct from the tensor-register partial SWAP used in
`../spatial/`. On two occupied equal modes, exterior hopping acts by
`det([[a,b],[b,a]])=-i`. Tensor partial SWAP acts by `a+b=1`. For two
different modes, exterior hopping also produces two-wave occupation of one
cell; a tensor construction with exactly one wave per cell cannot do that.
Consequently the tensor model's reference and dispersion do not transfer
without an explicit sector intertwiner. The spatial Fock circuit is a
declared architecture with the same literal onsite source law, and its own
reference and dispersion questions.

## Joint translation and clock readout on the spatial occupation law

Use two identical cells with the same link space, Fock law and local clock.
Define their full tick to be

\[
 \mathcal U=H_{01}U_{F,1}U_{F,0}.
\]

Define T to swap both complete cells: it swaps the two link words and
exchanges the two wave-mode blocks by the exterior lift of their one-body
permutation. In an ordered occupation basis, swapping the blocks contributes
the CAR sign `(-1)^(N_left*N_right)`. Thus T is a unitary involution.
The two onsite clocks commute because they are even and have disjoint
support. T conjugates each onsite clock to the other. Conjugating the
dressed hopping by T also exchanges its two link arguments and its two
one-body summands, leaving the hopping operator unchanged. Consequently

\[
 [T,\mathcal U]=0,\qquad T^2=I.
\]

T exchanges the two local Gauss constraints. It therefore preserves their
joint positive subspace; it need not commute with each separate Gauss
operator on the entire ambient space. The restricted translation projectors
are legitimate observables on that jointly constrained Hilbert space.
Their physical-space extension can be written with the joint Gauss
projection on both sides.

Take the empty-wave reference with the identical normalized local link
pair `( |1> + |6553635> )/sqrt(2)` at both cells. It is fixed by both
T and the full tick. Prepare one wave in the left cell in mode
`(port7, internal3, repair0, address0)` and apply the two positive-Gauss
preparations, retaining their complete success/failure records when viewed
as preparation instruments. The normalized successful state psi is
orthogonal to the empty-wave reference by wave number. This is an explicit
one-wave preparation, not spontaneous production from empty occupation by
the number-conserving idle clock.

Measure the complete translation projectors and the complete two-bin clock
instrument on that same psi:

\[
 P_\eta=\frac{I+\eta T}{2},\qquad
 K_\zeta=\frac{I+\zeta\mathcal U}{2},\qquad
 B_{\eta,\zeta}=K_\zeta P_\eta,
 \quad\eta,\zeta\in\{+1,-1\}.
\]

The P projectors sum to I on the physical space. Unitarity gives
`sum_zeta K_zeta* K_zeta=I`; commutation retains the momentum outcome.
All four branch probabilities sum to one, and all four output states
satisfy both Gauss constraints. The executable evaluates the literal
onsite pair and hopping, verifies the commutator on the prepared state,
and retains these complete records:

| Translation eigenvalue | Dimensionless cell momentum | Clock phase bin | Probability |
|---:|---:|---:|---:|
| +1 | 0 | 0 | 1/120 |
| +1 | 0 | pi | 59/120 |
| -1 | pi | 0 | 1/4 |
| -1 | pi | pi | 1/4 |

The momentum marginal is `(1/2,1/2)`. On a clock eigenphase phi, the
phase-instrument response is `(1+zeta*cos(phi))/2`. These broad finite
clock kernels are not exact energy eigenprojectors. The momentum labels
refer to the declared two-cell cyclic translation; no physical distance,
clock duration or three-dimensional dispersion is inferred. The executed
tick also has positive weight outside the original local flux-face orbit,
so the same calculation retains actual link transport. The receipt gives
that probability exactly.

The table is conditional on the normalized successful preparation. From
the definite raw link words and localized wave, the four Gauss outcome
pairs each have probability `1/4`. Retaining the three failed preparations
and continuing the positive-positive branch gives original-trial weights
`1/480, 59/480, 1/16, 1/16` for the displayed records, plus three failure
weights `1/4` each. Their total is one. Separate full-tick flux detection on
the normalized prepared excitation has moved-flux probability `1/450`.

This supplies a joint momentum-phase precursor on one spatial occupation
law, state, reference and clock. A three-dimensional relativistic dispersion,
energy branch and physical mass require further source-specific analysis.
