# Identifying the complete current-resolving repair law

Nine source-executed scalar measurements identify all nine real coordinates
of the existing [record-law family](../record_laws/README.md). This determines
the **coarse definite-current instrument on all 78 total-two states**, with
arbitrary spectator and reference coherences. It extends the [C/E experiment](../observer/source_discrimination/README.md): an admissible
collision law agrees with E on the entire binary sector, yet this experiment
separates it. A reference-uniform bound propagates measurement uncertainty
through finite feedback histories.

The executable result comprises 36 compiled probes and ten complete
two-event histories, five with untouched quantum references. The histories
retain raw records, aborts, quantum outputs, a current-dependent source
intervention and both repairs on the same state. They also run through the
existing [canonical observer runtime](../primitive_observer/README.md),
including serialized checkpoint validation and continuation. No new observer
runtime or checkpoint format is introduced.

The repository calls the contrasted laws **C and E**. Neither the existing
requirements nor these synthetic executions establish which law describes
nature. The result supplies a complete distinguishing experiment within the
classified family, rather than presupposing that those two candidates exhaust
the possibilities. It does not identify arbitrary fine record refinements or
coherent dilations from coarse current observations.

## One observer contract

| Component | Executed specification |
| --- | --- |
| State | Twelve-port total-two carrier; the existing nine-qubit code represents its 78 configurations. All 512 raw amplitudes, including unused codes, survive source execution. |
| Repair | One addressed seam, one of four explicitly compiled instruments below, fresh fair coin and blank records per call. |
| Preparation | A coherent pair of native configurations, with phase 0 or pi/2. Reference tests entangle that pair with an untouched two-dimensional reference using an emitted source word. |
| Read | A source-compiled rank-one plus projector on the current-translated pair and its entire coherent complement. Every current, both pointers and all aborts remain. |
| Public records | `(current,pointer)` with current -1,0,1 and pointer -1,+1, plus abort. Coin, contrast and auxiliary records remain in the raw execution but are inaccessible to the common controller. |
| Feedback | After public `(0,+1)`, apply emitted Z on source bit zero before the second repair. No branch receives a second preparation. |
| Clock and fields | Each attempted repair advances the addressed observer's event clock, including zero-current repairs. The integrated signed current and pointer correlation are computed jointly from the same complete history. Aborts halt further requests. These are dimensionless event/count quantities. |
| Continuation | Serialize the first event's canonical history; import by independent prefix replay; resume the second event. Compare every final quantum block with direct execution. |
| Refinement | Narrowing admissible probability tubes narrows the complete-history bound. Enlarging the finite call cap spends that bound linearly. This is statistical/operational refinement, not a claim of spatial or continuum refinement. |

The coherent proper-code load/export operations, blank preparations,
addressing and classical read services are declared apparatus. The source
controls are the existing six-level generators; the implementation executes
their pulse words. Exact target operators are used only as comparison oracles.
The apparatus is not a derivation of an autonomous physical clock or a QFT.

## Complete identification theorem

Write the existing admissible table as `(k,g,h,m,n)`, with complex k,g,h,n
and real m. Its positive matrices `G_j` have the fixed classical diagonals.
The classified current instrument is

```text
Phi_j(|x><y|) = G_j(key(x),key(y)) |x+jD><y+jD|.
```

Seam locality, the classical update, number-phase covariance, endpoint-swap
covariance and the committed definite-current record contract are precisely
the premises reused here; the two covariances are declared premises of the
construction, not consequences of the classical rule. Their classification supplies the full map once
the table is known. It does not restrict the comparison to populations.

Choose the existing representative configuration for each seam key below;
the other ports hold the remaining load. Prepare
`(|x> + exp(i theta)|y>)/sqrt(2)` and read the plus projector on
`(|x+jD> + |y+jD>)/sqrt(2)` after the selected current.

| Coordinate | Input seam keys | Selected current j | Baseline b | Phases theta |
| --- | --- | ---: | ---: | --- |
| k | 00,11 | 0 | 1/2 | 0, pi/2 |
| g | 00,01 | 0 | 3/8 | 0, pi/2 |
| h | 11,01 | 0 | 3/8 | 0, pi/2 |
| m | 10,01 | 0 | 1/4 | 0 |
| n | 01,02 | -1 | 3/8 | 0, pi/2 |

**Theorem 1.** These nine unconditional joint probabilities determine the
entire classified current instrument, and the exact inverse rejects tables
outside its CP domain.

*Proof.* The two diagonal terms contribute
`b=(P(j|x)+P(j|y))/4`; the cross terms contribute
`Re(z exp(-i theta))/2`, where z is the named multiplier. Consequently
`Re z=2(p_0-b)` and `Im z=2(p_pi/2-b)`. The nine-coordinate response matrix
is `(1/2) I_9`, with determinant `1/512`. The reconstructed table is checked
against every principal-minor CP condition from the classification. The
formula then fixes every native matrix unit, hence the full complex-linear
map and all its reference extensions. No division by a successful-current
probability occurs. The zero table has an admissible coordinate cube of
radius `1/16`: absolute row sums give positive diagonal margin at least
`3/16`. Thus the family has real dimension nine, and fewer than nine scalar
affine statistics cannot identify an open subset. This last statement is
not a lower bound on joint multi-outcome settings or trials. QED.

An exact point inverse is appropriate for exact laws. Noisy measurements can
produce a non-CP point estimate even when the underlying law is admissible.
The evidence interface retains probability tubes; it never clips or projects
such an estimate to manufacture a law.

## From calibration to complete histories

For two admissible tables let `d_z` bound the modulus of their multiplier
difference. Real-plus-imaginary absolute values suffice; the executable exact
bound also majorizes sqrt(2) by 3/2 when needed.

**Theorem 2.** The half diamond distance of their common-current instruments
is at most

```text
epsilon = min(1, max(d_k+2d_g, d_k+2d_h, d_g+d_h+d_m+d_n, d_n)/2).
```

*Proof.* The seam-key projectors `P_a` are mutually orthogonal, including
arbitrary spectator multiplicities. For any positive normalized state with
an untouched reference, set `p_a=tr(P_a rho)`; then
`||P_a rho P_b||_1 <= sqrt(p_a p_b)`. Current translation is an isometry
within its supported input keys. The equal classical diagonals cancel.
The triangle inequality for the direct-sum current output and
`2 sqrt(p_a p_b) <= p_a+p_b` bound its half trace norm by half the largest
absolute coefficient row sum. Those four rows give the displayed expression.
The channel difference is Hermiticity preserving, so optimizing these states
and arbitrary references gives its diamond distance. QED.

If each of the nine probability differences is at most delta, a complex
multiplier has real-plus-imaginary difference at most `4 delta`, and m at
most `2 delta`. Thus `epsilon <= min(1,7 delta)`. If two admissible tables
both fit the same radius-r tubes, their mutual probability discrepancy is
at most **2r**, giving `epsilon <= min(1,14r)`.

For identical admissible intervening instruments and public-record feedback,
telescoping channel differences and trace-norm contractivity give

```text
complete_output_distance <= min(1, preparation_error + N epsilon + terminal_error).
```

Here N is a worst-branch repair cap; all stopping and timeout outcomes remain.
An observed or average call count cannot replace N without another theorem.
The result covers arbitrary untouched references and any derived current,
clock or field statistic on this common output. An event probability changes
by at most this distance; a bounded observable with spectral range width W
changes in expectation by at most W times it. Comparing two different fine
coin/contrast interfaces, or intervening coherently on their environments,
is outside this theorem. A shared seam law must actually apply at each
compared call; heterogeneous seams require their own calibrations/budgets.

For the retained radius `1/1000`, two-call history uncertainty is at most
`7/250`. This is an exact probability-tube statement. It is separate from
the floating source execution residuals and from physical apparatus accuracy.

## Executed alternatives and a decisive blind spot

| Emitted candidate | k | g=h | m | n |
| --- | ---: | ---: | ---: | ---: |
| C | 1 | 1/2 | 0 | 1/2 |
| E | 1 | 0 | 0 | 0 |
| C_phase_quarter | 1 | -i/2 | 0 | 1/2 |
| E_collision_half | 1 | 0 | 0 | 1/4 |

The quarter-phase candidate executes C followed, on zero current, by an
emitted pi/2 phase on the waiting single-load states. The collision candidate
executes C and always copies seam parity at zero current. At moving current
it uses an additional source fair coin and copies parity only on one choice.
Its moving branch is therefore the half C/half E mixture; its zero-current
branch is E. All auxiliary outcomes are retained. These are actual compiled
operations, not table values injected into source execution.

The collision candidate and E agree on **every binary-shell input operator**.
There is no 02 or 20 key there, so n is invisible. Any history confined to
that shell and using the common interface remains identical, even with
references and adaptive controls. In the collision probe their unconditional
`(-1,+)` probabilities are respectively **1/2 and 3/8**. This proves that a
successful C/E comparison restricted to the binary sector cannot identify the
general repair law. The probe here removes that obstruction within the stated
family. Additional exact controls vary h, m and n independently near the
interior and exercise the imaginary quadratures.

## Resources, numerical reductions and evidence

`results.json` binds the implementation and source receipts. It retains all
seven probabilities for every probe, all 49 public probabilities per history,
raw outcome counts, complete-state comparison diagnostics and resource
envelopes. The producer compares every quantum output block for all 36
single-event probes as well as the ten continuing histories; matching only
the selected probability is insufficient. Independent endpoint controls
also check every emitted repair column on all 78 states and complex
reference inputs outside the prepared probe chart.
Raw branch costs include every requested pulse (including zero
angle), proper-code load/export, coin, record and pointer service, feedback
and preparation. Resource maxima range over all branches, including zero
probability outcomes; they are not expected operating costs. Full tapes and
amplitudes are regenerated rather than stored as a large artifact.

The raw execution has 512 amplitudes and retains small numerical leakage.
The canonical runtime uses eight active work codes, including both padding
codes. Before restricting it, the bridge measures the complete stacked-Kraus
tail and checks both full and reduced trace-preservation residuals. For full
stack V and its embedded reduction W,
`||V-W|| (||V||+||W||)/2` bounds the represented map difference in exact
arithmetic. For a raw final vector `(a,b)`, reduction changes its unnormalized
state by at most `||a||||b||+||b||^2/2` in half trace norm. These quantities
are evaluated separately: the final-vector term is added to the complete
output comparisons, and the stacked-map diagnostics are guarded separately.
No state is normalized or polar-corrected. The numerical evaluations use a fixed `2e-9` guard and
are **not certified machine-roundoff bounds**. Zero failure outcomes remain
in the validated raw menu. Anti-Hermitian discrepancies cannot disappear
through symmetrization: comparisons use the singular-value trace norm.

Reference preparations use the runtime's supplied-input interface, which
deliberately charges no preparation. The adapter adds the independently
executed joint preparation cost explicitly. Each raw cost is checked against
the canonical componentwise service envelope. The reference is not reset,
measured or used for feedback.

`calibration.analyze_counts` accepts all nine ordered settings with seven
integer bins each, including aborts, summing to the attempted-trial total.
Its exact schema and an example row are in [calibration.py](calibration.py).
An immutable `calibration.sampling_design` declares the nine attempted-trial
counts and sampling radii separately from the observed counts. Each observed
total must match its declared horizon. The horizons and radii must be chosen
before, and independently of, all analyzed outcomes; matching the design
does not establish that chronology. Apparatus/service errors must also be
supplied explicitly.
Under this fixed design and independent attempts within each setting, the
simultaneous sampling
failure bound is `min(1, sum_i 1/(4 m_i r_i^2))` by the variance bound and
Chebyshev plus a union bound. It is conservative and does not require
independence between settings. Counts do not establish these premises.
Outcome-dependent stopping or choosing a radius after inspecting outcomes
does not satisfy this guarantee. A retained bounded Bernoulli stopping
counterexample exceeds the incorrectly reused fixed-sample failure bound.
The separate continuing binary C/E classifier does not supply collision
preparation or this independent-trial guarantee.

The analysis reports compatible/excluded locked candidates and the diameter
bound for any two admissible tables in the tubes. An empty candidate list
does not rule out the whole family; full-family feasibility then is undecided. No physical observations are fabricated in this packet.

## Reproduction and remaining obligations

```sh
python -m model.repair_identification --check model/repair_identification/results.json
python -m pytest -q -W error model/repair_identification
python -m model.verify --hashes-only
```

The finite component is complete: full-family identification, executable
distinguishing operations, canonical continuation and a composed error budget.
This does not settle physical law selection. Independent observations or a
separate selection principle must constrain the repair table; apparatus
accuracy and physical clock/field correspondences need evidence. Fine record
access and coherent continuation require their own operational specification.
The original large scalar benchmark and continuum/relativistic recovery are
not closed by these smaller finite histories.
