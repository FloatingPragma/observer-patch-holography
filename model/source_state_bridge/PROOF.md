# Direct identification of the native state and a field model

The named native observer has a nonuniform stationary state and a compiled
modular clock. Their source specification is explicit. The results below
test direct identifications of those objects with the separate Higgs
benchmark. They do not exclude an implementation with additional quantum
resources, a different source family or a restricted observation menu.

## Exact state and ideal modular target

The actual default program in `native_program.interpreter.program` uses
port weights `(25,25,8,8,8,8,8,8,8,8,8,8)` on the binary total-two sector.
The [local repair proof](../local_repair_family/README.md) derives the unique
quantum stationary state for its fixed positive selector and endpoint
coins, including contraction on reference-entangled inputs. Its basis
probability is the product of the occupied port weights divided by their
sum. Enumerating the actual 66 configurations gives

| Occupation of the two distinguished ports | Eigenvalue numerator | Multiplicity |
| --- | --- | --- |
| two | 625 | 1 |
| one | 200 | 20 |
| zero | 64 | 45 |

The common denominator is `Z=625+20*200+45*64=7505`. Thus the state is
faithful and nontracial. Modulo a scalar, its modular Hamiltonian is

```math
K=-\log\rho_*=(\log Z-\log625)I+\delta(2I-n_0-n_1),
\qquad \delta=\log(25/8)>0.
```

It has levels `0,delta,2delta` modulo that scalar. On all 4356 matrix
units, the commutator `[K, .]` has gap multiplicities
`{-2delta:45,-delta:920,0:2426,delta:920,2delta:45}`. No floating logarithm
is needed for this census. In particular its normalized derivation obeys
`D(D²-I)(D²-4I)=0`, where `D=[K,.]/delta`.

Let `J:A_target -> M_66(C)` be an injective linear observation map satisfying
`D_native J=J D_target`. For `p(t)=t(t²-1)(t²-4)`, the identity
`0=p(D_native)J=Jp(D_target)` and injectivity give `p(D_target)=0`.
An injection in the opposite direction constrains only its image, not
the whole target algebra.
A common positive clock rescaling changes the gap unit, not the five-gap
pattern. This is a necessary condition on that exact equivariant map.
It is not a numerical exclusion of the Higgs frequencies, whose displayed
precision ladder is not an interval certificate.

The exact spectrum concerns the ideal stationary-reference modular
target. The record-inferred source words have their separate inference,
finite-preparation and compiler/export errors. They are not exactly this
unitary or a calibrated physical time. The state and potential are
consequences of the declared port weights; the weights themselves are
source parameters. `record_clock_member` explicitly supplies an H66-to512
address embedding and nine-qubit coherent work services. The original
carrier and that larger work register are different operational resources.

## A change of representation cannot purify the full state

Let `V` be an isometry from the native 66-dimensional carrier into any
larger space. Its image state `V rho_* V*` has the same 66 nonzero
eigenvalues. For a target density `tau` supported on a projector `Q` of
rank at most `r`, the variational characterization of trace distance gives

```math
D(V\rho_*V^*,\tau)\geq 1-\operatorname{Tr}(Q V\rho_*V^*)
\geq 1-\sum_{j=1}^{\min(r,66)}\lambda_j^\downarrow(\rho_*).
```

The last inequality is the finite variational bound for a rank-r
projection, obtained by diagonalizing the state: its diagonal projection
weights lie between zero and one and sum to r, so their maximum pairs
with the largest eigenvalues. For every pure target the bound is
`1-625/7505=1376/1501`. The top native eigenvector attains it. For a
finite-prepared input within trace distance `epsilon` of `rho_*`,
contractivity and the triangle inequality give the lower bound reduced
by `epsilon`, capped below by zero.
This uses the controlled source marginal; it does not replace the full
growing preparation archive by an independent state.

This compares the **full image state** with a pure target. An isometry
into field tensor memory can instead leave the field pure and place the
mixedness in memory: `V|i>=|Omega> tensor |i>` is an exact countercontrol.
Tracing that memory, or a reset/cooling channel, can produce a pure field.
Such a preparation needs its actual source operations, reference resources
and complete output accounting. The theorem excludes a mere full-state
relabeling; it does not exclude pure subsystem preparation. No rigorous
ground-uniqueness claim for the Higgs benchmark is needed: the comparison
applies to any prescribed pure vector.

## One quantum carrier with arbitrary classical records

Fix a logical target dimension `D` and a single communication/storage cut
where the entire retained quantum system has dimension `d`. Any finite
classical record may also cross that cut. There is no additional retained
quantum memory, transmitted environment, quantum bypass or pre-shared
entanglement across it. Let a complete encoder instrument have Kraus maps
`A_jl:C^D -> C^d` and classical flag `j`. A flag-dependent complete
decoder has Kraus maps `B_jk:C^d -> C^D`. Failures are retained as ordinary
complete outputs; there is no success conditioning.

The composite channel has Kraus maps `K_jkl=B_jk A_jl`, each of rank at
most d. For the normalized maximally entangled vector `Phi_D`, its overlap
after applying this channel to one half is

```math
F= D^{-2}\sum_{jkl}|\operatorname{Tr}K_{jkl}|^2
\leq d D^{-2}\sum_{jkl}\operatorname{Tr}(K_{jkl}^*K_{jkl})=d/D.
```

The inequality follows from
`|Tr K|<=||K||_*<=sqrt(rank K)||K||_F`; the equality uses completeness.
Testing the Bell projector proves

```math
\tfrac12\|\mathcal D\mathcal E-\mathrm{id}_D\|_\diamond
\geq\max(0,1-d/D).
```

The same reasoning is the standard Schmidt-number restriction on local
channels and classical communication; see
[Terhal and Horodecki, *A Schmidt number for density matrices*](https://arxiv.org/abs/quant-ph/9911117).
The source-specific content here is the bound for the actual carrier and
the complete target-sector dimensions. For `d=66`, the neutral `D=438`
and compensated charged `D=1596` sectors give respectively `62/73` and `255/266`.
The code reads these dimensions from the complete exact Higgs basis.
An exact reversible full-algebra representation is therefore impossible
through this specified single-carrier interface.

The cut condition matters. Sequentially extracting information into a
decoder's additional quantum memory violates it. So does using the named
observer's entire 512-dimensional coherent work register as the transmitted
resource; this is an existing source service, not a hypothetical loophole.
Its larger dimension alone does not construct the Higgs action or state.
The bound tests arbitrary logical inputs with references, and does not
exclude a comparison restricted to one prepared state and six probes.
Many native carriers, protected encodings, additional quantum memory and
controlled effective subalgebras require their own common-model maps.

## Executable checks and the unresolved source map

`bridge.py` reconstructs the stationary weights from the actual program,
checks all source modular gaps, and computes the rational state and
channel bounds. Exact small encoder/decoder controls retain classical
flags and saturate the channel bound when D is a multiple of d. Other
controls preserve every classical basis state while destroying coherence,
and explicitly move mixedness into residual quantum memory. These tests
check the hypotheses; the all-dimensional proof is the argument above.

The needed positive result is a specified source operation/observation map
joining its state and clock to the interacting action, with all additional
quantum resources, preparation errors and complete records counted. The
Higgs scalar, hopping, gauge terms and prepared ground do not follow from
the displayed product-reference logarithm. No claim excludes a different
encoding; no physical particle spectrum, unit calibration or Einstein
limit follows from these necessary conditions.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/source_state_bridge
```

The integrated receipt consumes `bridge.certificate()` and `SOURCE_FILES`.
