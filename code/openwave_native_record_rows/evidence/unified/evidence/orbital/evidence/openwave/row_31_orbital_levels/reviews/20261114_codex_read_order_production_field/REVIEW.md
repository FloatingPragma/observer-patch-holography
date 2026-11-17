# Same-group review: corrected native production field and clock

Reviewer Codex MFR, rows21–25 lane. This review retains an initial failed
decoder check and the corrected retry separately. It is same-group review,
not G6 or a hydrogen-row award. After the correction below, no blocking
issue remains in the finite production source/field/clock claim.

## A real failure was found and corrected

The initial executor assigned a NumPy column view `u=proj[:,index]`, then
normalized `u` in place after appending `proj` to its projector list. That
operation changed one column of the supposed orthogonal projector. My
independent reconstruction failed the recorded band projections, with a
maximum readout discrepancy1.8948645374959483. It independently reproduces
the bad decoder exactly, identifying the cause rather than treating the
disagreement as numerical tolerance. `VERIFICATION.json` preserves this
failed initial review. The producer retained all original files under
`failed_initial` and posted a preregistration amendment before the retry.

The corrected executor copies the column before normalization and checks
projector symmetry, trace, idempotence, resolution and an independent
polynomial construction. `verify_correction.py` compares every corrected
band readout with independently reconstructed original-source projectors;
the maximum error is exactly0 in this computation. Every native event,
initial location, occupation count, clock, censor, constant-sector and curl
record is identical between original and corrected attempts. The population
receipt is unchanged except its archive hash. Thus the original complete
native replay transfers by exact data identity, while the previously invalid
band-readout claim is replaced only by the corrected, independently checked
one. Neither the original failure nor its source costs are erased.

## Complete source-record checks

The independent reader imports neither producer nor native code. It checks
all82,008 local endpoint/spectator census records, all262,144 global outer
macro records and all2,613 adaptive inner records. Every preparation, edge,
coin, retry, target/escape predicate, proposal shortcut and feedback choice
is regenerated from the declared seeds. All9 held proposals are actual
native zero-load events. The768 producer full-vector checks were reviewed
as code and accounted for separately, not claimed as my own literal calls.
All347,533 producer calls are therefore accounted for.

No hidden source location or external reset enters the reconstructed
controller. The known current location determines whether an internal
proposal can matter. Relocation uses current-confirmed native moves from
the retained previous inner endpoint. Every timeout would be retained;
none occurred. The matched disabled-feedback protocol retains the inner
work and executes its declared outer comparison.

Exact rational orbit recurrences independently reproduce the capped probe
probabilities. Preparation can depend on history, so the exact operator is
the almost-sure-preparation/capped-probe comparison. The fully capped
implementation has per-active-attempt error at most12*2^-128 and the stated
history bound. The proof correctly keeps these two processes distinct.
The observed absence of failure does not remove their theoretical difference.

The active proposal counts are adapted, not a fixed number of independent
resets. The chosen radius is a conservative Hoeffding-Azuma bound after a
union over active-menu lengths and channels. The exact-preparation
conditional mean uses escape-by-cap/2, with the preparation allowance added.
All observed/expected move totals satisfy it. No precision field covariance
measurement is inferred from these small active-event counts.

## Same-source field and finite clock scope

After correction, all internal bands, constant fields, curls, occupation
counts, first-hit endpoints and actual ground-mode record readouts are
independently reproduced. Source clock tick counts are[5,2,2,3]; control
counts are[3,1,1,2]. Their final open intervals are respectively
[5251,18363,1203,8281] and[6686,23319,8944,204]. All remain in the archive.

The ideal finite spatial matrices and all27 levels in each positive internal
band were recomputed independently. The population tables refer to the
source-derived exact-preparation model at block length100,000,000; they are
not empirical estimates from the eight32768-macro tapes. The occupation
contact is retained. Internal feedback annihilates the constant sector,
so the massless curl source is unchanged. For each positive spatial level s,
the period3 first-hit equation and Green ratio give1/(1+s), not the
infinite-slice exponential. Every next-plane hazard remains1/(2W), since
feedback changes only internal proposals. Variable inner-event costs are
retained separately from the macro count.

The source small member uses a=1, inner boundary2 and outer period3. It is
an actually run resource setting of the parameterized controller, not a
point falsely claimed to execute the asymptotic Coulomb diagonal. The
separate shell-limit proof is not reviewed or accepted by this finite
execution review. Physical nuclear charge, particle identity, irreversible
scattering and platform acceptance are not inferred from these records.
