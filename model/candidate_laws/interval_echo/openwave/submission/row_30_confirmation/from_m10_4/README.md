# OW-30 registration from the eleven-contact M10.4 pair

**Ready for author and maintainer review.** At the reviewed OpenWave revision, the M10.4 source is `09f062f77f7f540f2783f128ae6314c5d8b818ea`, its pair retains the row-4 partial grade, and row 30 is untested. The authors’ proposed M10.5 source is the published v3 revision `92370140edfbdb0ac2be647b2b8166d821e5f1a2`. It commissions a coherent strand of loaded dimers on the same seam law. Clocks and separated spatial receivers measure frequencies and wavelengths. Both the older rejected two-board spatial comparison and the three-board v3 validation are disclosed in the request. [Current public candidate](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#L82-L95) · [Official row](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L467) · [Original frozen package](../README.md)

| Review stage | What a reviewer can check |
| --- | --- |
| [01-request.patch](01-request.patch) | Adds an author request, roadmap backlog entry, previously disclosed results, launcher and exact frozen inputs. The official task path, source pointer and decision are absent. The sole pin still selects M10.4. This is the patch for a request for review. |
| [02-maintainer-decision.patch](02-maintainer-decision.patch) | A separate draft for a maintainer who affirmatively adopts this source under these exact terms. The operative task, affirmative decision, sole canonical pin and source pointer must first land together. The earlier M10.4 source, review, terms and grade remain in its history. This private rehearsal carries no official decision. |

The platform head used for both patches is `c531265eb6085725c6268e156804499c810ba5e4`. At that revision there are official tasks through M10.4. A separate author-side single-load proposal also requests M10.5. The two protocols claim the same task pathname. Coordinate numbering with the other authors and the maintainers before filing. If that number goes to the other protocol, this frozen executable would require an explicitly versioned registration for another number. Its terms, scientific manifest, binding code and selector domain include M10.5 or v3 identities. Merely renaming the patch would leave the checker pointed at the old task.

The frozen checker requires the November 21 OpenWave base to be an ancestor on the official first-parent chain. The intervening adoption of M10.4 satisfies this requirement. From the first landing of the exact official task onward it checks every snapshot for one consistent source pin, decision and pointer. This permits a review request at a separate path before the decision. The verification fixture rejects terms posted one commit early even if the decision is supplied later. Its separate wrong-pin control also rejects. [Offline verification](verify.py) · [Verified identities and outcomes](receipt.json)

The selectors are deterministic functions of the first qualifying official Git commit. Whoever constructs that commit can preview it and can vary its bytes. This binds the selection to the resulting first-parent history. Unpredictability or impartial selection would require separate evidence. The two patches were rehearsed without evaluating any public selector. A landing party can disclose any preview or comparison of prospective commits. In a public replication the authenticated OpenWave launcher recomputes the selectors from official history. Low-level helpers accept caller-provided seeds and cannot establish that provenance. [Frozen launcher and selector](https://github.com/FloatingPragma/observer-patch-holography/blob/92370140edfbdb0ac2be647b2b8166d821e5f1a2/model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation/run.py)

## Verify the review packet

Use a complete OpenWave main history containing the reviewed base, a frozen RER checkout at the proposed source with its native ancestor objects, and a working copy of this packet. Git metadata and blobs must already be present for a disconnected check. In a blob-filtered clone, Git may fetch missing objects on first use. This command builds temporary indexes and unreachable fixture commits while leaving the reviewer’s files and branches in place.

~~~sh
python3 -I -B "$packet/verify.py" --platform "$openwave" --source "$frozen"
~~~

The [receipt](receipt.json) identifies the exact resulting request and decision trees. The checker compares the entire execution payload, binding files and author disclosure with the legacy proposal, down to the Git blobs. It exercises the actual source checker on the current history, a request-only history, a private affirmative fixture and the two malformed histories. For a platform or source checkout tracking current official main, the frozen public status commands remain authoritative.

## Reproduce the validation on a private review branch

The following commands assume `openwave` is a separate clean checkout at the exact base above and `python` selects the versions in the frozen requirements. Apply only the request for this local replay. Use a fresh directory for the validation records.

~~~sh
git -C "$openwave" apply --index --check "$packet/01-request.patch"
git -C "$openwave" apply --index "$packet/01-request.patch"
runner="$openwave/openwave/xperiments/m10_oph/research/scripts/m10_5/run.py"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
"$python" -I -B "$runner" check-sources --source "$frozen"
"$python" -I -B "$runner" acquire --mode validation --source "$frozen" --release "$release"
"$python" -I -B "$runner" replay --mode validation --source "$frozen" --release "$release"
cmp "$release/manifest.json" "$packet/../validation.json"
~~~

The private in-platform replay passed on all three previously exposed validation inputs. It regenerated 46, 62 and 62 count arrays and reproduced the complete archival manifest byte for byte. On the real official head, both the frozen source and platform status calls returned `unaccepted` (5), as did the public acquisition entry before any release or draw. The science score therefore remains the existing local validation. A public test is pending admission. In a future acquisition retain the per-board verdicts. The aggregate verdict prioritizes invalid and incomplete records ahead of a measured rejection, so a missing board can coexist with a completed failed board. [Archived validation](../validation.json) · [Verdict implementation](https://github.com/FloatingPragma/observer-patch-holography/blob/92370140edfbdb0ac2be647b2b8166d821e5f1a2/model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation/acquire.py)
