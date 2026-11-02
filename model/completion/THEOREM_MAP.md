# The finite campaign's theorem-to-model map

28 October 2026. CP-13 / [#1137](https://github.com/FloatingPragma/observer-patch-holography/issues/1137),
with a CP-02 reference construction. The [machine-readable map](theorem_map.json)
binds eleven reviewed results to twelve existing source files, their inspected
statement/proof anchors and SHA-256 hashes. These are scoped attachments to the
existing [premise register](../../tracking/premise_register.json) and
[generated premise map](../PREMISES.md), not a replacement scientific ledger.
The
[accepted finite map](finite_attachment_map.json) resolves
CP-13's integration contract. CP-02 and the physical refinement/dictionary
joins remain open; the earlier campaign's restricted service domain is
preserved below.

The [finite campaign](campaign.json) fixes the original total-two carrier,
tagged C/E mixtures, identity holds, positive repair/exploration floors and
bounded classical control. It excludes new coherent operations. This review
inspects actual existing theorem statements and the native model maps. It
neither recompiles the whole Lean library nor claims to review every theorem.
The campaign is a bounded component investigation; its boundary is narrower
than the full completion program.

| Result | Actual model attachment | Scope and owning gate |
| --- | --- | --- |
| Conditional resampling and stationarity | Binary66 configurations; uniform reference; a fixed seam's exterior occupations and total endpoint load define the fibre. | The native diagonal transition is uniform resampling. This is a classical marginal statement, not idempotence of the complete C quantum channel. CP-01/02; PR-07. |
| Stationarity of mixtures | `mixture_stationary` takes fixed weights `a : C -> R`, independent of source state. | Fixed mixtures preserve the reference. Record feedback can correlate source and scheduler; the exact countercontrol below prevents extending this conclusion to all source marginals. CP-01/02; PR-07. |
| Native quantum mixing | Existing connected carrier, fixed independent edge weights, and records discarded. | Preserve its exact uniform C/E gap and binary reference; do not apply those constants to an adaptive source/controller channel. CP-02/03. |
| Complete history recovery | Each recorded native path is a weighted partial injection; retained currents determine translations. Tagged mixtures and identity holds preserve the form. | Full tape plus final source supports the stated population recovery. A bounded memory summary alone does not inherit the full-tape theorem; export and halt are separate outputs. CP-01. |
| Native coherence boundary | Every allowed branch is a partial injection or identity, with classical conditioning and diagonal ancillary preparation. | Diagonal inputs stay diagonal. Existing coherent inputs can still retain quantum output; no assertion that the whole channel is classical. CP-01/02/06; PR-65. |
| Normal endpoint uniqueness | The Lean theorem requires observation preservation and completeness; endpoint existence also requires weak normalization. | Continuing native load motion is not a terminating normal-form solver. The distinct fixed-history constraint solver retains its own result. CP-01/13; PR-01. |
| Gibbs/modular identity | The retained benchmark chooses a faithful Gibbs state of its supplied Hamiltonian on the full work algebra, with `t=-beta s`. | Preserve that positive construction. The new finite reference below supplies neither those quantum controls nor a physical clock. CP-02; PR-23/15. |
| Native defect invariants | The [new exact branch-algebra certificate](defect_algebra.md) checks all complete C/E branches on H78 and H66. | Scalar commutant excludes nontrivial bare-source Hermitian charges commuting with every branch and its adjoint, and nontrivial reducing sectors. This is not an exclusion of arbitrary expectation-fixed observables, histories, metastability or refinement. CP-03. |
| Actively corrected native memory | Existing encoders and branch-conditioned coherent recovery satisfy their code identities. | Keep the achieved capability as an external control for this campaign. Its extra controls are excluded here; it is not autonomous defect protection. CP-01/03. |
| Einstein tensor completion | Actual theorem arguments include common tensors, null matching, Ward/Bianchi divergence identities and connectivity. | These physical maps are missing here. Markov stationarity does not supply generalized-entropy stationarity at fixed proper volume. CP-09/10; PR-23/25. |
| Original q=3 scalar history | The 8,192-level onsite recipe and its mathematical/resource bounds remain as originally stated. | The large native word and 106-qubit field/pointer state remain unexecuted. This campaign neither executes nor retires that broader obligation. CP-00/01/13. |

The closest existing finite stationarity result is
[`mixture_stationary`](../../Lean/Thermodynamics/FiniteConditionalRepair.lean).
Its proof moves the same coefficient `a c` outside a sum over source states.
A controller correlated with that source does not supply this identification.
The new contribution is an exact native counterexample to the attempted
transfer, followed by construction of the correct joint reference for the
same declared controller. No general stationarity theorem is being replaced.

## The adaptive reference countercontrol

Start with the actual uniform binary source and a blank classical bit `m=0`.
At each native attempt, choose a seam from

```text
q(e|m) = (1/2)/30 + (1/2) 1[e=e_m],
e_0=(0,1), e_1=(0,5).
```

Execute the usual fair-coin C or E repair and set `m'=1[current>0]`.
Both laws expose this current and have identical diagonal transitions, so
every tagged convex mixture has the same source/bit marginal. The first
attempt leaves the source marginal uniform while correlating it with the
bit. The second attempt gives **30 nonuniform source probabilities**, with
maximum absolute deviation **1/960** and total variation **13/1320** from
`I_66/66`. All quantities are exact rational numbers; there is no sampled
trajectory or numerical rank threshold.

The complete distribution and source custody are in
[adaptive_reference_results.json](adaptive_reference_results.json).
Keeping the preferred seam fixed, or setting exploration to one, preserves
the uniform source marginal. Identity holds only slow the effect in the
separate controller that leaves its bit unchanged on a hold. The production
observer writes a hold record with zero current instead. The shared
countercontrol and stationary construction both use repair probability one,
so this distinction is absent from their actual comparison.

This refutes a particular stationarity transfer. It does **not** show that
the adaptive process lacks a reference state.

## A reference for the actual joint process

For the same controller, form its complete classical kernel on
`(x,m) in binary66 x {0,1}`. Every transition has denominator 120. Solve
`K^T pi=pi` with `sum pi=1` over the rationals and check the original residual
exactly, including the equation replaced by normalization.

[stationary_reference.py](stationary_reference.py) constructs a **unique
classical stationary distribution**, positive on **129 of the 132 joint
states**. Exact transition-support analysis finds one closed recurrent class
of size 129 and three transient singleton classes. A retained self-loop of
probability `14/15` makes the recurrent class aperiodic. Thus the finite
Markov chain converges from every initial distribution, including the
declared uniform-source/blank-bit preparation. This result supplies no
finite attempt budget for a desired preparation accuracy.

The source marginal is positive on **all 66 binary states and nonuniform**.
The full exact distribution, support and residual certificate are in
[stationary_reference_results.json](stationary_reference_results.json).
The independent test replays this distribution through the native event
executor and recovers it exactly, rather than only rechecking the matrix
used to solve the equations.

The [full-record observer attachment](observer_reference_results.json) lifts
this bit reference through the actual CP-01 instrument and verifies exact
stationarity again. For repair probability one and one retained FIFO slot,
the C endpoint has 128 positive full-record blocks and E has 192; each has
4,218 positive source/record diagonal entries. The projection back to the
source/current-sign bit equals the solved reference exactly. This is a
reference for the finite source and latest-record marginal; exported older
records continue to accumulate, so the full growing history is not declared
stationary. Each record block is faithful only on its positive source
support. The [observer guide](OBSERVER.md) retains its complete record and
resource contract.

Because native branches preserve diagonal inputs, this also constructs an
invariant diagonal density of the quantum source/classical-bit channel.
Classical uniqueness does not by itself prove uniqueness of all its quantum
fixed states. The joint density is not faithful on the entire 132-dimensional
space. It is faithful on its 129-dimensional support; the source marginal is
faithful on H66. Neither assertion makes it faithful on H78 or on an
unrestricted source/memory algebra.

The algebra choice matters. Compressing to the support does not implement
access to all its matrix units. The modular automorphism flow is trivial on the finite
commutative record algebra. Nonuniform source weights would define
nontrivial phases on a larger matrix algebra, but the admitted native
operations have not implemented those phases or a calibrated detector for
them. Moreover, the source marginal belongs to a correlated joint stationary
state, not an autonomous closed source dynamics. Regional modular clocks,
their comparisons and the common physical energy normalization remain CP-02
obligations.

## Reproduction and retained limits

From the repository root:

```bash
python -m model.completion.adaptive_reference --write
python -m model.completion.stationary_reference --write
python -m model.completion.adaptive_reference --check
python -m model.completion.stationary_reference --check
python -m unittest model.completion.test_adaptive_reference model.completion.test_stationary_reference model.completion.test_theorem_map
```

The attempt count is a native event parameter. These calculations use no
laboratory time or rest-energy unit and add no physical comparison. Their
source binding and structural checks certify which computations and theorem
statements were reviewed, not the truth of an unproved physical attachment.
Positive modular-clock, conditional 3+1-field, active-memory and Einstein
kernel evidence remains on its existing domains. No frozen packet is edited.
