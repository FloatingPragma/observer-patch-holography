# Probabilities on the complete native program grammar

For each prescribed preparation, the same positive normalized state gives
the complete public record probabilities for every finite native decision tree.
This includes written collars, free contact histories, record feedback, early
stops, classical random streams and labelled record unions. The [Lean induction](../../../../Lean/EventAlgebra/NativePublicEventGrammar.lean)
proves the result for arbitrary finite dimension and depth. The
[executable source bridge](grammar.py) checks the native 99-state operators,
complete records and independent event signatures against the normal candidate
entry points. The quantum branch response is a declared construction choice;
this result supplies no derivation of that response from scalar consensus.

## The full finite grammar

Fix the literal hardware allocation of any registered finite board, including
its ternary material ports, private collar memories, reference registers,
optional fuel and enrolled reader/controller registers. A program is either a
halt with a retained label, or a finite primitive instrument followed by one
program for each labelled outcome. Each primitive is a finite family K_i with
`sum_i K_i* K_i = I`. There is no bound on the finite program's depth and no
restriction to an inventory, protected code, selected effect table, real
subalgebra or blank collar preparation. An individual valid program and its
allocated hardware are finite. The theorem quantifies over all such programs.

The primitive interface is a mathematical superset of the native menu; it
asserts no native reachability of an arbitrary supplied Kraus family. Native
instances enter by their actual maps: the Hermitian full-state generator gives
unitary free evolution; the echo is the specified 99-state unitary; native
references and routed permutations are unitary; the complete sharp material
or enrolled pointer reading is a projective family; and independently supplied
classical random streams give scalar branches `sqrt(p_i) I`. A deterministic
finite controller chooses the continuation from its preceding records.
Reversible apparatus operations are permutations on their entire declared
register alphabet. Unobserved environmental/Kraus labels are auxiliary indices,
with their public record map applied after summing the complete branches.
No private collar reading is introduced.

The actual source maps, entry-point calls and semantic boundaries are in the
[receipt](results.json). In particular, `candidate.run` delegates to
`replay.main`; `replay.generate` executes the material clock and cohort source.
The default clock uses native fair draws, sharp initialized material registers
and the coherent hand's squared-amplitude terminal read. The full bulk program
uses declared independent 0.04 tag offers and 0.16 matter offers, then native
sharp-sector transitions. Their exact complete eight-word initialization law
is retained. Quantum channels and this initial classical law belong to the
construction. Fixed pseudorandom seeds reproduce sampled arrays; they do not
prove independence of physical randomness.

## State, effects and equivalence

For each complete history h, multiply its primitive maps chronologically to
get L_h. Induction on the decision tree proves

`sum_h L_h* L_h = I`.

A normalized positive initial law gives a diagonal initial density. Coherent
source preparation can be included as part of the same program. More general
supplied positive normalized input densities and arbitrary quantum references
are also admitted by the theorem. The runtime branch rule is the specified
quantum update and trace: `p_h=Tr(L_h rho L_h*)`. Cyclicity derives
`p_h=Tr(rho L_h* L_h)`. Summing histories of any public event S gives

`E_S=sum_(h in S) L_h* L_h`,

`p(S)=Tr(rho E_S)`.

The event operators are positive and bounded by identity, the complete event
has value one, and disjoint labelled unions add. These facts hold for every
generated public effect of every finite program. Labelling, forgetting private
indices, appending a conditional unitary, refining a record with a complete
coin and finite stopping all preserve the relevant identities. This establishes
a positive normalized state on the complete generated effect system, rather
than testing only a finite selection of settings.

Reading equivalence has a state-independent certificate. Equality of event
operators is equivalent to equality on a complete Hermitian operator signature.
The Lean theorem `operator_signature_ext` proves the symmetric/antisymmetric
matrix-unit signature is injective. The executable independently reconstructs
a d-dimensional Hermitian effect from d² normalized pure probe values: its
d diagonal probes, and the `(e_i+e_j)/sqrt(2)` and
`(e_i+i e_j)/sqrt(2)` probes for each i<j. This tests every matrix coordinate,
not equality on the one target input. Equal certified effects then give equal
probabilities for every admitted source and reference state. The probe theorem
is a mathematical equivalence certificate; arbitrary native preparation of
all probes on the entire written-memory hardware is not asserted.

The theorem's probability rule is explicitly primitive quantum semantics.
It does not infer the same law from context-wise normalization. A complete
three-outcome scalar instrument with effects `(I/4,I/4,I/2)` and the separate
binary instrument `(I/2,I/2)` have the same union effect I/2. A sampler which
squares the primitive weights and renormalizes each context gives respective
union probabilities 1/3 and 1/2. Both contexts remain positive and normalized.
The exact discrepancy 1/6 is checked in Lean and rational arithmetic. Thus
record normalization alone does not establish implementation independence.

## Native execution and limits of an entry point

The 99-state regression includes coherent density entries across charge
blocks, arbitrary dirty collars, a native free contact, echo, complete material
reading, conditional native evolution and retained refusal/early-stop labels.
A 1:3 record refinement and a native postrecord echo have the same nine
material-event operators as the original program, while their terminal maps
and refined records differ. The independent d² signature has 9,801 coordinates.
The program propagates every native input column, rather than fitting one
source probability table. A separate regression runs a complete clock trip
through the candidate entry point. For each owner's published elapsed count,
the 99-state effect calculation reproduces the terminal instrument's two
first-rail probabilities. Both owners recover the recorded outcome bit.

The native clock makes at most twelve address attempts at each of at most
eleven subsequent walk steps. An out-of-range attempt is part of the history;
a timeout releases the primary instead of resampling an accepted cohort. For
q=3, the exact retained timeout probability is `(1/4)^12`. Record-map coarsening
can hide the timeout reason while retaining its probability. For a five-bit
address field the one-service upper bound is 672 classical fair draws. Ideal
real uniform thresholds and finite 53-bit uniform thresholds differ by at most
2^-53 per scalar draw, so a finite N-draw execution has coupling discrepancy at
most `N*2^-53`, with operator arithmetic error separate.

Invalid allocation or ownership, an unprepared sharp read, numerical invariant
failure and exhausted finite elapsed registers raise host exceptions in the
normal entry point. They are not physical observer outcomes. The theorem
applies to complete valid native executions and explicitly represented
refusals/stops. It does not turn those exceptions into an emitted native record
or certify that every user-supplied journey respects capacity. The default
registered bounded schedules and the source's finite control menu have their
stated preconditions. No success-conditioned sample is used to assert their
unconditional normalization.

This proves the candidate FT-BORN criterion for generated public effects under
the declared quantum state, instrument and sampling semantics. The PR-03
operational route to the finite Busch–Gleason theorem also requires a native
common refinement for every coexistent pair in its target matrix algebra.
The unrestricted instrument constructor in Lean supplies no compiler for those
pairs on arbitrary written collars. Recovery of a unique ambient state requires
public effects that distinguish states.

## Replay

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.public_event_grammar.grammar
OPENBLAS_NUM_THREADS=1 .venv/bin/python -B -O -W error -m model.experiments.interval_echo_consensus.public_event_grammar.grammar
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/public_event_grammar
cd Lean
lake build EventAlgebra.NativePublicEventGrammar
```

`--write` is exclusive first-record creation. Discrete records compare exactly. The deposited observation source and its
reviewed successor with optional attendance, scheduling and history recording
are pinned separately. The replay driver also has a reviewed revision that protects
files during publication. Both drivers run the same measurement generator. Each
accepted combination of producer, observation source and driver has a recorded
digest for all three. Other source bindings must match exactly. Floating replay uses absolute tolerance 1e-11.
