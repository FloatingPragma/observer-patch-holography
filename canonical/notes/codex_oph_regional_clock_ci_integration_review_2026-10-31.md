# Regional clock CI integration: reviewed provenance, not blanket rehashing

Date: 2026-10-31. Worktree and branch remain the dedicated regional-clock
agent worktree and `agent/k3-oph-regional-clock-implementation-0-review-2026-10-31`.
Historical construction base: `1f14865fe560f5788a990a1b25884d6f91d12858`.
Previous validation base: `4d485074bc6a6e80f722dc9744634166dea592b7`.
Current merged validation base: `a444d48606ab3611f02ac0e0e29a883f0cac56fc`.

## Verdict by claim

| Claim | Verdict | Evidence |
|---|---|---|
| TM-11 construction conclusion survives the reader-guide change | GO | The full reviewed q=3 paragraph is byte-identical to the last validated source |
| Imported theorem-map provenance resolves to the reviewed live source | GO | Every source hash and anchor checked; four downstream pins checked separately |
| The former `7acd...` pin identifies a retrievable historical source | HOLD | No match in the 191 available historical file snapshots; not treated as evidence |
| The divided-step results extend to h=0 | REQUEST CHANGES to that claim | Division by h and a one-hop h=0 counterexample; original nonzero hypotheses retained |
| Regional implementation and its resource/error claims changed | NO | Regenerated receipt has identical scientific payload |
| Inventory accounts for the new active README | GO | Exactly one added surface, no removals or allowlist changes |
| CP-02 or a physical clock/field join is accepted by this repair | HOLD | No acceptance or physical premises changed |

## Source review before accepting the updated custody link

The CI log supplied by the owner expected construction-source hash
`7acd841f74df47285baaa99a72067bbcf8e2f5a4c8abeac52421aa63d879742e`.
Searching all available path history, including merge history, examined
191 snapshots and 139 distinct file hashes. None matched that hash.
At `1db86257d`, the commit introducing that pin into the theorem map,
the actual committed CONSTRUCTION.md already hashed to
`257af3dba03e0bab0b13e2990846d8509cddb81aa12c7d80bf67fedc25686e45`.
LF/CRLF variants of the known sources do not explain the stale pin.
This is an unresolved historical custody pointer, not proof of a different
mathematical statement. No missing historical bytes were invented.

The last validated construction source, at both the original construction
base and the previous validation base, hashes to
`bf4c2130fbd7f77292f9247e83bc7e1695b6f70ecb654a666d742de642764888`.
Comparing that actual source with the current source gives one local prose
change around the record-law classification: E adds seam-parity readout,
C0 is a distinct coherent current-only instrument, and classically forgetting
C's coin leaves C rather than producing C0. This correction does not alter
the q=3 preparation paragraph, its resource/error limitations or its
explicit statement that the enormous word and 106-qubit state were not
executed.

Only TM-11 consumes the `construction` source entry. Its reviewed anchor,
complete supporting paragraph and unexecuted/outside-campaign disposition
remain valid. The mapping is historical to its own finite campaign; it does
not silently reverse the later V2 q=3 adoption decision. The other ten
mapping conclusions and every reviewed source anchor are preserved by the
imported fix. No source field is updated locally by an indiscriminate scan.

Upstream `a444d4860` already supplies the corrected single source pin and its
consumers. The theorem map's live hash is
`12db32c565d3233da6879b96c532745474e0d9410d7c18bf10c5614f0fb936ab`.
It was independently checked at all four consumer locations:

* claims/completion_program.json: stage progress evidence and next-wave evidence;
* code/particles/runs/status/postdiction_ledger.json: both rendered counterparts.

The canonical premise-register navigation and the human completion-program
link also remain valid. These imported files are unchanged relative to the
new upstream. `build_completion_program --check` and the distinct canonical
model/theorem-map generator check both pass; neither generator is used to
rewrite unrelated history.

## h != 0: retained mathematical restriction and changed presentation

Previously, the flagship's finite-results table claimed a conserved hopping
charge and explicitly required a declared nonzero step. The rewritten
flagship removes that table entry. Its remaining supplied native-field
operations do not assert this divided-step hopping-charge theorem. The
relevant included current finite-construction material likewise does not
make that assertion.

The restriction itself remains essential. In
Lean/Geometry/PortChargeMinimalCoupling.lean,
`hoppingCurrent = -((q/h) * hopIndicator)`, and `hopping_continuity` still
requires `(hh : h != 0)`. Its proof uses cancellation with that hypothesis.
The corresponding charge-fixed and transported field theorems retain their
nonzero-step hypotheses. Independently, at h=0 and a unit charge crossing
one edge, the load change is exactly (-1,1) while h times the current
boundary is zero; the continuity equation fails. This is not merely a
stylistic condition that can be dropped.

The claim registry retains the nonzero-step assumptions and zero-step
falsifiers. The detailed observers paper still states both the declared
nonzero step and the declared nonzero field step. The imported surface-test
correction therefore scopes only the flagship phrase check to an actual
hopping-charge claim there. Its canonical claim checks, Lean-signature
checks and both detailed-paper checks remain unconditional. No hypothesis,
falsifier or mathematical source is weakened by this integration.

## Owned regeneration and diff review

The regional receipt was regenerated only through
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.regional_clock.verify --write`
and then replayed with `--check`.
Comparing every field with the previous receipt shows only:

* `validation_upstream` changes from the previous validation base to the current base;
* three source hashes change: the imported V2 contract, the regional README
  documenting the base, and the regional contract documenting that base.

All target identities, regions, confidence qualifications, literal word
hashes, exact error bounds, floating diagnostics and resource budgets are
unchanged. The original construction base remains recorded. No original
training, preparation or heavy K3 pool was rerun.

The active-surface inventory was regenerated only with
`python3 tools/check_axiom_consistency.py --inventory`. Its count changes
from 775 to 776, adding exactly `model/regional_clock/README.md` as an active,
non-allowlisted `model/**/*.md` surface. No old entry or exclusion changes.
The inventory writer's regression test also writes this same generated
file, without changing its payload.

The merged completion workflow now splits upstream checks across jobs.
Every upstream job and its ordered steps were independently compared with
the merged workflow. The only additions remain the two regional path
triggers and two regional test/replay commands, now in `program-native`.
No check is removed or duplicated.

## Verification and NEXT

The targeted pass checks the regional suite, theorem-map custody, inventory
guard, post-r2029 public hypothesis surfaces, replay workflow coverage and
the updated pytest shard scheduler. Its only initial failure was the
inventory writer running in the filesystem sandbox's read-only worktree;
the authorized rerun of that writer's complete 13-test suite passes.
The remaining 418 targeted tests passed without altering their assertions.
The regional receipt replay, inventory consistency, completion-program
generator and canonical model/theorem-map generator checks pass.

The previously failing full `program-contract` command is replayed with
the original test selection, warnings as errors, and one BLAS/OMP thread.
The original mandatory post-r2029 command group is also replayed with its
owner's warning policy. An exploratory stricter invocation of that group
encounters a collection deprecation in the locally installed pytest for an
existing itertools.product parametrization; no test source or warning
filter is changed to hide it. The strict targeted hypothesis tests pass.

The old head's GitHub Actions logs were inspected directly, rather than
inferring all failures from the supplied excerpt:

| Historical failed job | Actual failed check | Current repair / replay |
|---|---|---|
| Completion program / program-contract | theorem-map construction SHA mismatch | Imported reviewed source pin and all consumers; full job command replay |
| Completion program / program-native | regional receipt differs from replay | Owned regional regeneration and successful complete replay |
| Mandatory shard 1 | committed inventory omits the new regional README | Exact one-entry inventory regeneration; 13 guard tests and live check pass |
| Mandatory shard 9 | obsolete flagship hopping-charge phrase assertion | Reviewed scoped upstream fix; complete 176-test post-r2029 command group passes |
| Observer regression shards 1, Linux/Windows | same theorem-map SHA mismatch | Same reviewed provenance repair; strict custody test passes |
| Observer regression shards 2, Linux/Windows | `KeyError: program` in two workflow coverage tests | Imported all-job workflow traversal; strict replay/sharding tests pass |

Historical runs are
[Completion program](https://github.com/Arithmon/observer-patch-holography/actions/runs/41206140626),
[Mandatory Suite](https://github.com/Arithmon/observer-patch-holography/actions/runs/41206140524)
and [Runnable observer](https://github.com/Arithmon/observer-patch-holography/actions/runs/41206140532).
Successful checks on the old head were not mistaken for validation of the
new tree, and no failed job is rerun against that obsolete SHA.

The 176-test mandatory command retains its owner's normal warning policy.
Pytest emits one collection deprecation for the existing iterator
parametrization; it is visible and unsuppressed. The actual mandatory owner
command does not use `-W error`; the additional exploratory invocation had
made that warning fatal. Replaying the owner command preserves its policy
and all h!=0 assertions. Local replay is on Linux; no local Windows result
is asserted.

NEXT remains source-service and physical normalization review under CP-02.
This CI repair establishes current provenance and test custody; it does not
accept any new physical identification or normalization.

Final result: the complete previously failing `program-contract` command
passes **422 tests with warnings treated as errors**, in 916.72 seconds.
Its longest setup is the existing all-66 source-word certificate (502.22
seconds); it was executed once. All four remaining job commands also pass:
acceptance audit, V2 check, independent program review and finite attachment
map check. The supplied log's theorem-map failure is therefore resolved on
the actual merged validation tree, rather than bypassed by test selection.

The strict targeted selection validates 419 distinct cases: 418 pass in the
initial sandboxed invocation, and the sole filesystem-blocked inventory
writer passes in its authorized 13-case suite rerun. The 176-case public
surface group passes under its unchanged mandatory-suite warning policy.
Final diff/whitespace, source-hash and generator consistency checks pass.
Only the regional receipt and the active-surface inventory are locally
regenerated; all imported theorem-map and consumer bytes match upstream.
