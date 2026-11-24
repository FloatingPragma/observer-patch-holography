# Review of the proposed state and decoder

**Author-side checks verified. Full-credit filing on HOLD.** OpenWave's row-30 cell for M10 is untested at its verified main revision. The source, registration and acceptance patch are local proposals. Source publication, a maintainer decision adopting this candidate, a public acquisition and the cell decision are outstanding. [Public cell](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L467) · [Proposed package](PACKAGE.json)

| Artifact | Identity |
| --- | --- |
| Frozen RER source | `92370140edfbdb0ac2be647b2b8166d821e5f1a2` |
| Proposed OpenWave tree, including author-side exposure disclosure | commit `cb16815014dfabfd96edaaa84fbc7cfbfcc0dfb1` on `d4fc3e4046c472f4ad38468df8dddfa1f97bac97` |
| Candidate declaration SHA-256 | `17486cf931c62cb4f5e2669dc062b69b4da5de1a4f8f6f39779f657ce0ac908a` |
| Terms SHA-256 | `4e88e5d7116f225d3b56a83d831212dfe5b9ee9d10d4443c422213d64b1ffccc` |
| Source manifest SHA-256 | `a6be1385ba2df9783f3449d51052102fe6b07464f17694c9fcf60b7f914c3d80` |

**Candidate admission.** The verifier checks an affirmative decision on this row and declaration, its source pointer, the sole canonical pin through each first-parent snapshot, published RER bytes and the actual OpenWave program and inputs. The synthetic history tests include a temporary pin change followed by restoration. The two real status paths returned `unaccepted` (5) before an official draw. The local patch includes the known rejection and subsequent author validation as a separate finding. [Binding tests](../../row_30_collective/confirmation/test_binding.py) · [Port tests](../../row_30_collective/confirmation/test_run.py) · [Prior-results disclosure](KNOWN_RESULTS.md)

**Acquisition integrity.** The exposed coupling mutation, a missing count, a corrupted count and interrupted acquisition exercise the independent outcomes. The original v2 validation replay returned `rejected` (2), with all 130 arrays present. A separate reviewer checked the shapes, types, budgets, digests and terminal differences of all 170 fresh v3 arrays and rehashed the 130 v2 arrays. All three v3 cohorts replayed from their saved seeds. An independent checkout imported from the source bundle reproduced the v3 validation verdict. A default-thread inspection returned `invalid` (4) because one 80-cell recurrence diagnostic differed by about 7.14 × 10⁻⁸. A fresh one-thread replay reproduced the 170 arrays and returned `passed` (0). The [replay instructions](README.md#reproduce-from-the-frozen-checkout) give the numerical-library settings. [Acquisition tests](../../row_30_collective/confirmation/test_acquire.py) · [V2 validation](v2/validation.json) · [V3 validation](validation.json)

**Spatial discrimination.** The fresh v3 minimum excess of the launch-label residual above its own threshold exceeded 791,000. The measured wave numbers passed the same held receiver data on each board, and all 22 measured frequencies passed the unchanged dispersion gates. The result concerns a supplied coherent, one-dimensional material. The lower-band continuum tolerance does not uniquely select a continuum law. [Spatial tests](../../row_30_collective/confirmation/test_spatial.py) · [All three cohorts and interpretation](README.md)

The frozen collective experiment suite passed 54 tests. The scorecard suite passed 23 tests. The proposed patch applies to its official base and the OpenWave roadmap linter passes. The staged state decision and the validation data confer no official row award.
