# Scalar-number and quantum-link-loop spectroscopy

## One source and three exact operators

The source, complete neutral Gauss ground and clock are those of
[`periodic_spectrum`](../periodic_spectrum/PROOF.md). Both five-state quantum
Higgs memories, both six-state oriented links and all twelve active lepton
CAR modes are retained. The source coefficients are unchanged:
`y_e=omega=kappa=1`, `y_u=y_d=0` and oriented hopping magnitude `1/4`.
The apparatus acts on the same ground used for the charged and neutral
composite records.

The local scalar operator is the actual number projector at site one,
`N_H=diag(0,1,1,1,1)`. For a retained internal block `B`, define

\[
 R_B=\sum_{a,b\in B}U^{01}_{ab}U^{10}_{ba},\qquad
 A_B=(R_B+R_B^\dagger)/d_B.
\]

Each `U_ab` is the original raising matrix unit `|a,b><blank|`. The singlet
block has `d_B=1`; the weak-doublet block has `d_B=2`. Both oriented links
remain quantum. On the complete 36-state two-link memory, `R_B` maps the
blank pair to a vector of squared norm `d_B^2` and vanishes elsewhere.
Consequently `A_B` has eigenvalues `-1,0,+1` with multiplicities `1,34,1`.
The code verifies `A_B^3=A_B`, `Tr(A_B^2)=2`, Hermiticity and all 24 local
current commutators, twelve at each endpoint. The scalar projector commutes
with all twelve scalar currents. These full-register identities extend to
every matter occupation.

Translation swaps the two oriented links. Relabeling the two contracted
indices leaves `R_B` unchanged, which the exact 36-state swap verifies.
The two loop operators therefore preserve the ground's `T=+1` phase. The
local scalar operator has both even-sector translation components.
Neither a two-link trace nor a scalar cutoff number is identified here with
a continuum gauge field or an observed elementary boson.

## Complete interventions

The apparatus records an equal-probability choice among the three settings.
For the scalar setting, use success `K=N_H` and failure `F=I-N_H`.
For either loop setting use

\[
 K=A_B/2,\qquad F=I+(\sqrt3/2-1)A_B^2.
\]

The exact spectra imply positive failure operators and `K^dagger K+F^dagger F=I`
on every input. The explicit two-outcome dilation
`U=[[F,-K^dagger],[K,F]]` is checked as a unitary on the complete local
operator and record registers. The apparatus is even and gauge neutral.
It is a declared finite intervention; the idle source Hamiltonian is
unchanged.

Each raw output is checked against the inherited exact Gauss embedding.
Success maps are projected into both `T=-1` and `T=+1` sectors by their exact
spectral projectors. Their Gram effects sum to the unresolved success
effect. Together with the full unresolved failure map, the source effect
is exactly the identity on the entire 251-dimensional input block. No
failure operator is defined by taking a square root after translation
compression.

## Energy, translation and retained alternatives

Every isolated energy projector is retained with its full multiplicity.
The loop instruments have exactly zero success map into `T=-1`; all energy
records in that branch are retained with exact zero probability. The
odd-parity translation phases `+i,-i` also have zero probability. The
`T=+1` ground projector is retained as an elastic event, separately from
positive energy transfers. No basis vector inside a degenerate eigenspace
is treated as a distinct physical energy outcome.

For each setting and each energy/phase projector `P`, the joint probability
is `Tr(P K rho_ground K^dagger)/3`. Rational approximate ground vectors and
projector bounds come from the same periodic source. Since `||K||<=1`,
`(delta_ground+delta_P)/3` bounds the joint probability error. All intervals
are rounded outward at denominator `2^27`. The receipt includes the full
success sum, the retained failure interval and the lowest certified
positive inelastic record of each setting when one exists.

The failure state is the unnormalized quantum state `F rho_ground F/3`;
energy and translation are not read on that branch. The ground filter is
also retained: success has probability `1/3424`, and failure has probability
`3423/3424` with unnormalized state `(I_neutral-P_ground)/3424`. Multiplying
the subsequent conditional law by `1/3424` gives the original experiment
law. Both complete probability enclosures contain one. Successful energy
branches are never renormalized.

## Relation to other finite boson controls

[`particle_candidate_bosons`](../../../particle_candidate_bosons/PROOF.md)
uses a different eight-vertex, twelve-edge Wilson-character source and a
rare empty-matter sector of its selected density. Its cube multiplets and
residue bounds are not transferred to this periodic quantum-link ground.
[`particle_candidate_scalar_gauge_response`](../../../particle_candidate_scalar_gauge_response/PROOF.md)
uses a fermion-pair radial contraction and a different primary Wilson
source. The scalar number here is the explicit five-state quantum scalar
operator of this source.

The present statements concern scalar and closed-link detector channels,
complete finite probabilities and source-unit energy transfers. They supply
no observed Higgs, photon, W, Z or gluon identification, intrinsic Lorentz
spin, colored spectrum, continuum particle pole or GeV mass. The lepton
reduction contains no colored propagating target sector.

## Reproduction

The receipt binds all package code and proofs, transitive local imports,
native representation inputs, the consumed periodic-spectrum receipt and
its complete input custody. The consumed receipt's pins must validate
before generation or replay. `--write` regenerates the receipt; `--check`
performs independent read-only replay and compares certified intervals and
custody exactly.

Projector errors use the exact nearest-interval form of the inherited
full-complement bound, with ordered endpoints and approximate-energy
containment verified on every target block. The helper is shared with the
neutral-composite calculation; its rational bounds are unchanged.
