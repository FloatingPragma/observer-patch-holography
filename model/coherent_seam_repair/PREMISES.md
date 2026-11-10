# Executed chart agreement and finite-theorem premises

## One construction and its accepted comparisons

The [coherent seam program](PROOF.md) has one declared carrier, code, walk
clock, preparation service and apparatus grammar. Its gauge-coupled state
space has dimension 7200. A clock event is one application of `W_G=Q_G R`;
an optional apparatus event is separately recorded. The fixed five-level
code splits into charge eigenspaces of `Q5=diag(2,2,2,-3,-3)`.

An accepted chart comparison uses one finite-cover frame word at each named
source port. [agreement.py](agreement.py) compiles those words from the same
120-lift controller as the source transport. Its apparatus wrapper applies
the inverse frame, executes the owned basis detector or quarter-turn, then
restores the frame. This wrapper is an explicit architecture service. It
implements the transported query and does not assert an arbitrary effect
compiler. A record label means a query in that declared frame; the frame
word is part of the experiment specification.

The [exact receipt](agreement_results.json) binds the executable import
closure and dynamically loaded carrier bytes. It checks all sixty oriented
local operator identities, an independently executed transformed-link walk,
and every branch of a finite adaptive filter/phase/walk/detector experiment.
It also checks a strict wrong-apparatus control: a bare pointer has
probability `1/120` after a nonmonomial code rechart, whereas the correctly
transported query retains probability `1/60`.

## Finite agreement theorem

Let `V_a` be an implemented unitary frame at port `a`. Write `S` for their
direct sum, acting trivially on parity and seam address. Transform every
oriented link by

`G'_ba = V_b G_ba V_a*`.

The seam isometries obey

`m_G' = (V_a direct_sum V_b) m_G V_a*`,
`p_G' = (V_a direct_sum V_b) p_G V_a*`.

The reversible local exchange is therefore conjugated by the same frame.
The address reflector commutes with `S`, so

`Q_G' = S Q_G S*`, `W_G' = S W_G S*`.

The callable `ChartedObserver.repair` executes the transformed links
directly; it does not define the new walk by conjugating the old walk.
The certificate compares these independent executions. Each admitted frame
preserves the charge blocks, so it also commutes with `Q5`.

For either outcome of the complete preparation filter,

`K_+(W)=(2I+W)/3`, `K_-(W)=sqrt(2)(I-W)/3`,

polynomial covariance gives `K_±(W_G')=S K_±(W_G) S*`.
Consequently the actual reduced heralded reference transforms as

`rho_G' = S rho_G S*`,

where `rho_G=[(5/4)I+(W_G+W_G*)/2]/2620`. Its normalization, positive lower
bound `1/10480` and stationarity with its own walk are preserved. A changed
chart does not mean leaving the input density matrix unchanged.

The declared detector and phase grammars transform as

`{P_i,I-P_i} -> {S P_i S*,S(I-P_i)S*}`,
`I+(i-1)P_i -> S[I+(i-1)P_i]S*`.

The frame wrapper executes exactly these operators. For any finite adaptive
protocol made from this grammar, use the same setting-selection rule on
the same interpreted previous records. If `B_h` is the product of branch
operators for a complete history `h`, induction on its length proves

`B'_h = S B_h S*`.

Therefore for every input density matrix `rho`,

`B'_h (S rho S*) B'_h* = S (B_h rho B_h*) S*`.

Taking the trace proves equality of the entire history law, including
preparation failures and detector complements, and identifies the
conditional output states for every nonzero-probability history. Summing
histories proves agreement for every coarse-grained event. Completeness
follows inductively from the completeness of each instrument. The argument
also covers a finite retry protocol when each restart transports the same
preparation and preserves its recorded failure; discarding failures and
comparing conditional samples alone is a different experiment.

All identities tensor with the identity on an untouched reference `R`.
They therefore hold for correlated inputs `rho_HR` under `S tensor I_R`.
The retained preparation archive is not silently reset or asserted to be
stationary. If the archive is read, its joint state and record meaning enter
the experiment. The stated stationary reference concerns the reduced walk
algebra only.

The concrete adaptive experiment retains all four outcomes:

| Preparation | Detector | Exact probability in either chart |
|---|---|---|
| Success | Yes | `43537/3888000` |
| Success | No | `3700463/3888000` |
| Failure | Yes | `1181/1944000` |
| Failure | No | `70819/1944000` |

Their sum is one. This is an executed finite accepted-chart agreement
construction. It is not a theorem about every possible observer federation.

## Which finite premises this construction supplies

The canonical hypotheses are the
[premise registry](../../claims/axiom_registry.yaml) and shared
[A1–A3 statement](../../paper/tex_fragments/PREMISE_BASIS.tex).
Each application requires its own source-bound witnesses; counting formal
theorems does not supply their hypotheses.

| Interface | Source witness | Exact scope or obligation |
|---|---|---|
| Finite carrier | The pinned twelve-port, thirty-seam, twenty-face oriented icosahedron | A single carrier, not a directed refining spherical federation. |
| A1 central port records | Classical port addresses can label experiments | The coherent position projectors are not central in the full walk matrix algebra. Treating them as central changes the algebra and requires a separate interface. |
| A1 complete infinitesimal response | The declared geometric code response has an exact rank-twelve Gram and Lie image `s(u(3)+u(2))` | Completeness for every quotient-visible response, naturality under federation/refinement and the quotient separating optional apparatus controls require proofs. The larger apparatus algebra cannot be discarded by naming it internal. |
| A2 finite chart meaning | The theorem above supplies complete history agreement, transported apparatus and arbitrary untouched-reference covariance | The accepted domain is the implemented finite frame groupoid. Seam restrictions, triple overlaps and refinement of a growing federation require their own commuting diagrams. |
| A2 holonomy witness | The source-generated 120-cover controller has typed reversible arrows, recorded inverses and carrier actions on the same code | The [holonomy bridge](../../Lean/Screen/A2HolonomyBridge.lean) permits such a cover; it does not require an `A5 -> PU(5)` section. Its response-generated flow, naturality and endogeneity hypotheses retain their explicit domains. A machine formal identification with [SL2F5PortCover](../../Lean/Screen/SL2F5PortCover.lean) is separate from exact executable checks. |
| A3 reference covariance | The complete filter prepares a faithful reduced `rho_G`, and the same filter and reference transform jointly with `W_G` | The declared `alpha=1/2` filter is not an information-projection dynamics. An A1-generated constraint grammar, separating observer cover, natural positive weights and realization of the prescribed optimizer require constructions. Choosing the target state as its own reference does not derive those interfaces. |
| Finite quantum records | The implemented Kraus instruments give complete positive effects and normalized history probabilities in the declared complex-state semantics | This validates that quantum construction. It does not derive its probability semantics from classical repair counts. |
| Operational Born bridge | Finite events have additive probabilities and effects; chart-equivalent experiments agree | [OperationalBornBridge](../../Lean/EventAlgebra/OperationalBornBridge.lean) also requires implementation independence and joint refinement for all coexistent effects of its declared domain. Scalar-walk tomography spans Hermitian operators but does not implement every positive effect; the gauge walk splits into charge blocks. |
| Nontrivial readable modular action | The faithful reference is nontracial and fails to commute with a readable register projector | Its modular generator is a function of `W_G+W_G*`. Equality with the walk clock, geometric cap normalization and a relativistic time interpretation are separate conditions. |
| Gauge response and charge | Static links affect a complete detector law; `Q5` is conserved under idle evolution and admitted frames | Charge conservation and observable background holonomy do not establish dynamical gauge fields, physical spin/statistics, particle poles, masses or a continuum limit. |

The distinction between a spanning set of tests and an operational effect
domain is substantive. The finite theorem
[`nonlinearBinaryWeight_not_affine`](../../Lean/EventAlgebra/FiniteEffectClosureBoundary.lean)
exhibits a continuous, normalized binary-projector rule that is not Born.
This package avoids that inference by declaring its quantum instrument
semantics and proving the consequences on the apparatus it executes.

## Reproduction

```sh
python -m model.coherent_seam_repair.agreement --check model/coherent_seam_repair/agreement_results.json
python -m pytest model/coherent_seam_repair/test_agreement.py -q
```

The certificate uses exact arithmetic. Its finite executions support the
displayed examples; the operator proof supplies the universal finite-history
and untouched-reference quantifiers within the declared grammar.
