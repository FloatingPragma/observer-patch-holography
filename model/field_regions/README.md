# Regional CAR maps on the continuing interacting field

CP-08 / [#1144](https://github.com/FloatingPragma/observer-patch-holography/issues/1144)
and CP-13 / [#1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137).
This component supplies explicit spatial algebra and record restriction maps
for the same [Weyl/contact field](../interacting_field/DERIVATION.md). Its
Hamiltonian is `H=dGamma[-(1/8) sigma.sin(k)] + sum_x n_x0 n_x1`.
The declared cubes, their coordinates, two CAR modes per site and Hamiltonian
model time are inputs. A growing cube is a volume limit at fixed lattice
spacing. It is not a spatial continuum or physical Lorentzian reconstruction.

The nearest results are the existing
[`CPRegionalNet`](../../Lean/QFT/CPRestrictionNet.lean), its trace-preserving
conditional-expectation and tower laws, and the full-outcome CAR/source
history estimates in the field derivation. The existing
[`no_scalar_hom_of_matrix`](../../Lean/QFT/NoncommutativeWitness.lean)
explains why the older
[`FiniteCausalObserverNet`](../../Lean/QFT/FiniteCausalObserverNet.lean)
with a scalar bottom cannot admit its homomorphic restriction from a full
matrix factor. These are reused results. The contribution here is a concrete
signed CAR attachment, retained local particle sectors, compatible observer
records, and finite-volume embeddings for this particular interacting field.
The general proof below is mathematical prose with independent numerical
and exact finite controls; it is not a Lean kernel proof.

## A full-Fock split before sector restriction

For an ordered finite mode set `L` and `K subset L`, define

```
S_K : F(L) -> F(K) tensor F(L\K),
S_K |I> = (-1)^(#{(q,p): q in I\K, p in I cap K, q<p})
           |I cap K> tensor |I\K>.
```

The sign follows by anticommuting the retained creation operators past the
exterior creation operators. This is a unitary permutation of full Fock
bases. Restrict a joint state with any untouched reference by conjugating
with `S_K tensor I` and tracing the exterior. In particular, a global
**N=2** input produces local **N=0,1,2** blocks. No block is discarded and
no surviving trace is renormalized. Keeping only the local N=2 block would
postselect that both particles remained inside the region.

For each exterior occupancy `J`, the explicit Kraus map is
`R_J = (I tensor <J|) S_K`, restricted to the input N=2 sector. Their
completeness is `sum_J R_J* R_J = I_N2`; therefore the map preserves trace
and is completely positive, including after adjoining any reference.
For nested retained sets `K subset V subset L`, the two signed splits differ
only by a permutation of the two discarded factors. Tracing both gives
`Tr_(L\K) = Tr_(V\K) Tr_(L\V)`. This proves restriction composition even
when the retained density is mixed and includes all three particle sectors.

The [implementation](regions.py) exposes:

```python
regional_basis(kept_modes)
restrict_pure(space, vector, kept_modes)
restrict_density(space, rho, kept_modes, reference_dimension=1)
restrict_regional(rho, ambient_modes, kept_modes, reference_dimension=1)
```

`space` is the existing `TwoFermions` or `WeylLattice`. The output basis is
vacuum, ascending singletons, then lexicographic pairs. Reference indices
are fastest; a pure input of shape `(space.dimension, ref)` is one joint
pure vector, not a classical mixture. Retained modes must be sorted,
distinct global indices. The four modes at the origin and its +x neighbor
have an eleven-dimensional regional output. The maps accept arbitrary finite
matrices linearly and do not repair positivity by clipping eigenvalues.

## The actual local algebra and conditional expectations

On N=2, use **number-preserving** local observables. Their full-Fock CAR
inclusion preserves this sector and restricts to a unital star homomorphism
`j_K`. Two such algebras on disjoint mode sets commute: their monomials are
even, and all cross-region CAR signs occur in pairs. Odd CAR operators do
not commute on disjoint regions. Arbitrary even operators such as pair
creation also need not preserve N=2. The implementation rejects a local
number-changing operator instead of treating `P_N2 X P_N2` as a
homomorphism. Indeed `P a0* a1* P=0` while
`P a0* a1* a1 a0 P` is a nonzero projector.

For an ambient matrix `X` on N=2, the trace-preserving conditional expectation
onto this local algebra is

```
E_K(X) = direct_sum_n [Tr_out(P_n S_K X S_K* P_n)/d_n] tensor I_(d_n),
d_n = choose(|L\K|, 2-n),
```

conjugated back by `S_K`; absent sectors are omitted. `P_n` selects local
particle number n. Explicit Kraus operators in that block are
`I_local tensor |a><b|/sqrt(d_n)`. They prove complete positivity, ordinary
trace preservation and unitality. The map fixes the regional algebra and
is the trace-orthogonal projection onto it. Since `j_K(A_K) subset j_V(A_V)`
for `K subset V`, the orthogonal-projection identity gives
`E_K E_V=E_K`, as well as idempotence. This is the conditional-expectation
contract in `CPRegionalNet`, attached to actual CAR regions. The state
restriction and ambient expectation have different codomains and are not
interchanged in the API or proof.

There is also a concrete coverage witness if the **declared** generating
family contains every two-mode region. On N=2, `n_p n_q` is the rank-one
projector onto `|p,q>`. A local hopping operator between two modes, sandwiched
by such projectors, gives every matrix unit for configurations differing by
one occupied mode. The graph of two-element subsets joined by such changes
is connected; multiplying along paths gives all matrix units. Thus these
regional algebras generate the full N=2 matrix algebra. The reusable source
hopping/FSWAP controls can address those pairs. This coverage is a statement
about that supplied control menu. The two-site detector used below is not
claimed to reconstruct the complete field algebra from its four records.

The ambient expectations are **not multiplicative**. For four modes,
retain modes 0,1 and let `X=c2* c3+c3* c2`. Exactly `E_K(X)=0`, but
`E_K(X^2)` is the nonzero projector onto the sector with one exterior
particle. This is a concrete control for the existing homomorphic-retraction
obstruction, not a further no-go. Nor do separated marginals uniquely determine
a joint state: `(|0,2>+|1,3>)/sqrt(2)` and its minus-phase partner have
identical marginals on regions {0,1} and {2,3}, but orthogonal joint states.
The existing two-slot correlation controls have the same logical boundary.

## Complete histories and the finite-volume join

The receipt executes the existing field source compiler on the full N=2
spaces of the first two cubes. There are two complete Ramsey number reads,
with conditional phase feedback and another evolution between them. The
field is not reset. All four unnormalized quantum-record outputs are
restricted to the same two-site region. Their traces reproduce the complete
record probabilities, and restriction through the two-site region to one
site agrees with direct restriction. A separate admitted entangled input
retains a complex two-dimensional reference throughout both reads. Its
preparation is supplied; only the subsequent field/instrument words are
executed in that test.

A local number-preserving instrument `K_b` satisfies, branch by branch,

```
Tr_out(j_K(K_b) rho j_K(K_b)*) = K_b Tr_out(rho) K_b*.
```

This follows from the signed split, not from labeling an operation local.
It holds with untouched references. Summing all branches gives local
no-signalling for a disjoint retained algebra. Conditioning on a selected
outcome does not give that unconditional theorem. The source tests retain
both outcomes and verify the square against the independently constructed
local CAR effect. Source floating contraction residuals are numerical
diagnostics, not outward error certificates.

For nested cubes define the one-mode inclusion by identical coordinate and
spin labels and let `J` be its exterior square. The maps are isometric,
compose, and preserve all CAR orientation signs. They satisfy exactly

```
J* H_large J = H_small.
```

They do **not** satisfy `H_large J=J H_small`: hopping can cross the smaller
cube's boundary. The receipt retains a nonzero intertwining defect, so an
exact Hamiltonian compression is never substituted for a dynamics theorem.
For the initial apparatus-supported states and common local instruments,
reuse the existing complete-history Taylor bound

```
d_half(output_R, output_infinite) <= 2 (B T)^R/R!,
B=6c+U=7/4,  T=sum_j |t_j|,
```

plus the separate ideal splitting errors. Apply the same regional trace
map and use trace-distance contractivity; the bound therefore controls
regional quantum states, retained records and arbitrary untouched
references together, without a factor for the number of outcomes. Spatial
inclusions are exact; finite-volume dynamics converge rather than
intertwine. Bounds refer to the ideal compiled operation; floating and
physical service errors are separate. This gives a usable fixed-spacing
local observation limit and does not establish a UV limit, strict relativistic
microcausality or topology of physical spacetime.

## Premise attachment and remaining work

PR-16 and PR-52 require a source-selected physical event/geometry dictionary,
two-way order/cone faithfulness, calibrated count/volume and rods/clocks,
manifold and dimension/topology tests, and a controlled unique Lorentzian
limit. None follows from choosing cubic CAR labels. PR-53 requires
physical propagation/frame/clock/readout identification. PR-70, PR-71 and
PR-72 concern the committed edge-direction action, its physical continuous
operator sector and cofinal gluing into that comparison. The present cubes
supply no map from that source edge action to this Hamiltonian. Their
premise-register dispositions are unchanged.

The finite mathematical result is concrete: the supplied interacting
field has signed regional CP maps, local algebra/expectation identities,
source-instrument naturality and a record-preserving fixed-spacing limit.
The outstanding physical selection and continuum premises are explicit.

Run from the repository root:

```bash
python -m pytest model/field_regions -q
python -m model.field_regions --check model/field_regions/results.json
```

The [campaign](campaign.json) fixes the scope before replay. The
[receipt](results.json) binds its source and nearest theorem files. `--write`
uses exclusive creation and never overwrites an existing receipt.
