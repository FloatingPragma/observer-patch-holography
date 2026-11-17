# RA1: independent component review of the row-1 repair package of the auditor's root

Reviewer: Fable-B / `reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B`, 2026-11-14, 16:33 to 17:35 (+07). The reviewer wrote
none of the reviewed code. `audit.py` imports no candidate and no producer code.

Candidate: `openwave-tuning/next/OBSTRUCTION_REPAIR/row1/`, frozen by `MANIFEST.json`, SHA-256
`85f9eefcf4c2292c1a2f47f522149cf1e28852705376d3abc2a65024243c7f3a` (created 16:39:48). Two components:

- Part A: the adaptive phase decoder `adaptive_phase.py` with `verify_adaptive.py` (sync entry OR002).
- Part B: the neutral-start record-only pointer join `neutral_record.py` with `verify_neutral_record.py`
  (candidate id `native-neutral-pair-record-only-v1`, sync entry OR003).

This is a component review. It gives a verdict per statement (proved, executed, not established) and no verdict
on the literal criterion of row 1. No status is awarded. Official platform status of the row is unchanged.

## Bytes reviewed

First copy 16:33:37 (`consumed/candidate/`, before a manifest existed), second copy 16:40:50
(`consumed/candidate_snapshot2/`). The review is of the second copy. Its ten files equal the ten pins of the
manifest (10 of 10). Against the first copy: `adaptive_phase.py`, `verify_adaptive.py`, `ADAPTIVE_RESULT.json`,
`RECORDED_REPLAY.json`, `NEUTRAL_REPLAY.json`, `neutral_record.py` and `NEUTRAL_RECORD_EXAMPLES.json` (hashed live
at 16:33:37, `e28c52cb...`) are unchanged. `verify_neutral_record.py` and `NEUTRAL_RECORD_RESULT.json` changed
before the freeze: the six controls accumulate with "and" over all cases instead of "or", and the executions of
the wrong-pointer control are counted (34). `README.md` and `MANIFEST.json` are new. Only the frozen versions of
the two changed files were read. Fourteen of the fifteen source pins of the manifest equal the live files; the
fifteenth is the root's own review result, read after the verdict (last section).

The public package is unchanged: `PUBLIC_MANIFEST.json` `d293ba59...`, `SNAPSHOT_MANIFEST.json` `86c990ef...`, and
127 of 127 snapshot files equal their pins. Both components are additive files outside the package tree.

| File of the manifest | SHA-256 |
| --- | --- |
| `adaptive_phase.py` | `3d371995e460af04db9fa2335ae62615a14d72bedaacdecd3c1be75d4d48ecc5` |
| `verify_adaptive.py` | `4f92c69c2cb1a1e3fd6ce46b078c0d8caeb6e0c9f68e3594064f8e304d13b80f` |
| `ADAPTIVE_RESULT.json` | `999c6627a053ae249c4a4c19c5d80e9e7658c188f2affda086f97fb15fe751c4` |
| `neutral_record.py` | `6fd90489c748ab87d7f44a878f640d716c6fef8c9429e1380952df16b925da60` |
| `verify_neutral_record.py` | `b1ed86ce8d028ca81133a29f05879c182deed1618c0f1a7984453f5c0710134b` |
| `NEUTRAL_RECORD_RESULT.json` | `36ba2528b9c86b347dc5eceeffa01c123745e38143b898a884d3490a14d88aaa` |
| `NEUTRAL_RECORD_EXAMPLES.json` | `e28c52cb8eaad8bdcaa0db77016160dd02b16753d667bc27045300bd8abae164` |
| `README.md` | `503baab12f23a2b0a682c9921dc89f553a486f08a156d4797ca6c72347d218da` |
| `RECORDED_REPLAY.json` | `71549bb87a9234a90cfdb52cca5ee0c90cb7c68a7ee4dea525706ab688992ee4` |
| `NEUTRAL_REPLAY.json` | `935d9f2efcc6e76af510cc67b51b561a537757ca60e52f2a2439dfb446730daa` |

Every other file read is in `consumed/` with its hash in `MANIFEST.json` of this directory.

## Verdicts

| Statement | Verdict |
| --- | --- |
| A1 increment below pi/4 at n_e = 8 abs(F_e) + 8, for every F_e | proved, and executed |
| A2 the unwrapped sum is the exact degree | proved, and executed |
| A3 the density rule is a resolution rule: a field with another winding is read as the field's | executed on a decoder with a field input; not established for `read_phase` as shipped, which has no field input |
| A4 rejections: fractional input | executed |
| A4 rejections: a field that jumps by more than the bound | not established |
| A4 rejections: abs(F_e) beyond the resource cap | executed |
| A5 the present filing texts admit this decoder | not established: six passages state the opposite; list below |
| B1 the inputs of the join hold no per-port reference and no terminal load | executed, with two riders |
| B2 translation covariance and the single-pair invariant give a scalar reference and a charge in {-1, 0, +1} | proved (exhaustive on the finite sector) |
| B3 the counts 23,940 and 288; failures kept with their weights | executed |
| B4 what remains supplied | classified below; no calibration |
| B5 the objection "relative to a remembered start" | moved and reduced, not removed |

## Part A: the adaptive phase decoder

### The encoding

Public bytes: `topology.phase_degree`, lines 124 to 147 of
`code/openwave_charge_quantization/workspace/openwave-tuning/row_completion/ROW_01/topology.py`
(`9f30570a...`, pinned in `SNAPSHOT_MANIFEST.json`), and `curve` in `code/openwave_row1_native/native_charge_winding.py`,
lines 99 to 112 (`8526fbf0...`). On the dual edge e = (a, b) of the icosahedral carrier (a, b among the twenty
dual vertices), with integer current F_e:

    z_e(t) = exp(i phi_e(t)),   phi_e(t) = g_a + (g_b - g_a) t + 2 pi F_e t + A_e sin(2 pi t),   0 <= t <= 1,
    g_v = 0.31 sin(v + 1),   A_e = 0.23 cos(a + b + 1)        (deformed);   g = A = 0   (undeformed).

The admissible fields of the declared encoding are these two per integer tape. The field is a function of the
tape. The package routine samples 64 subdivisions per edge and refuses abs(F_e) > 3 (line 130). The task file
says 32; 32 is the count of the law controls in the filing's G4 text. The ancestor routine already takes its
count from the tape: `max(16, ceil(8 D / pi))` with D the derivative bound.

### A1. The increment bound: proved

The increment between consecutive samples at n subdivisions is

    Delta_k = (g_b - g_a + 2 pi F_e) / n + 2 A_e sin(pi/n) cos(pi (2k+1) / n).

n = 8 abs(F_e) + 8 is a multiple of 8, so both signs of cos(pi/n) occur among the k, and the largest increment on
the edge is exactly

    max_k abs(Delta_k) = abs(g_b - g_a + 2 pi F_e) / n + abs(A_e) sin(2 pi / n).

With abs(g_b - g_a) <= 0.62, abs(A_e) <= 0.23 and sin x <= x this is at most
(2 pi abs(F_e) + 0.62 + 0.46 pi) / (8 abs(F_e) + 8), the candidate's bound, and

    pi/4 - bound = (1.54 pi - 0.62) / (8 (abs(F_e) + 1)) > 0     for every F_e,

because the inequality is 0.62 + 0.46 pi < 2 pi and does not contain F_e. `audit.py` checks that inequality in
exact rationals with 333/106 < pi < 355/113. The bound is analytic, for all F_e, on the mathematical
interpolant. It also covers every field of the same form with abs(g) <= 0.31 and abs(A) <= 0.23.

Exact worst case. On the envelope the largest increment is
E(m) = (0.62 + 2 pi m)/(8m + 8) + 0.23 sin(2 pi/(8m + 8)) at abs(F_e) = m. E increases with m (checked for
m = 0 to 20000) and its supremum is pi/4 = 0.785398, not attained. On the declared geometry
(max abs(g_b - g_a) = 0.605122, max abs(A_e) = 0.229991) the worst edge is dual edge 10 with positive current:

| abs(F_e) | n_e | declared geometry | envelope E | candidate bound | pi/4 minus bound |
| --- | --- | --- | --- | --- | --- |
| 0 | 8 | 0.215517 | 0.240135 | 0.258142 | 0.527257 |
| 3 | 32 | 0.647092 | 0.653294 | 0.653584 | 0.131814 |
| 4 | 40 | 0.674836 | 0.679798 | 0.679947 | 0.105451 |
| 128 | 1032 | 0.781118 | 0.781311 | 0.781311 | 0.004087 |
| 12494 | 99960 | 0.785354 | 0.785356 | 0.785356 | 0.000042 |

The closed form equals the direct evaluation of the phase on 30 edges and 13 currents to 1.0e-13. Executed: on
192 own reads the largest sampled increment is 0.780991. In binary64 the read at abs(F_e) = 12494, the largest
that the default budget allows on one dual cell, has largest increment 0.785347 and error 0.0: rounding is far
below the margin 4.2e-5 at that current, and below the distance 3 pi/4 to the value pi where a principal
increment stops being the true one.

### A2. Exactness of the degree: proved and executed

Proof. Each increment is below pi, so the principal value of the sampled ratio is the true increment. The sum
on an edge telescopes to g_b - g_a + 2 pi F_e. Around the oriented pentagon of a dual cell the g terms cancel
(own check of `GEOMETRY.json`: every cell is a closed oriented pentagon whose signs are the incidence of its
port, and every dual edge bounds two cells with opposite signs), so the degree of the cell at port v is
sum of s F_e = (B F)_v, and of a region the sum over its ports.

Executed with an own decoder (own contour assembly, gluing and closure checks, `atan2` increments):

- Witness `witness_four_turns.json` (`6d4023c2...`), replayed with an own endpoint rule: the tape ends with -4
  turns on seam 1, max abs(F_e) = 4, B F = G - n, degree of region {5} equal to 1, final loads one unit at
  port 11. The own adaptive read is 1 (float 0.9999999999999999), resolutions {1: 40, 5: 8, 18: 8, 21: 8, 22: 32},
  101 samples, largest increment 0.6693661085056494. The amendment's record: degree 1.0000000000000004, the same
  resolutions, 101 samples, 0.669366108505649.
- Own sample, seed 2026111401: 96 tapes with signed currents up to abs(F_e) = 128, random regions, both
  deformation flags: 192 of 192 reads equal sum_R (B F)_v, largest error 5.68e-14, 156 reads with nonzero
  degree, degrees from -567 to 629.
- Comparison with the amendment. `verify_adaptive.py`, run from a copy, reproduces `ADAPTIVE_RESULT.json` byte
  for byte. On the 192 own reads the copy of `read_phase` returns the status, degree, sample count and largest
  increment of the own decoder in 192 of 192; the two float degrees differ by at most 1.1e-13.

A fixed 64 subdivisions, with the cap lifted, also read the witness exactly (largest increment 0.418410). The
refusal at four turns is a declared cap and no limit of the resolution: 64 subdivisions certify abs(F_e) <= 15
under the public guard pi/2 (72 of 72 own tapes with abs(F_e) <= 15 read exactly), and the registered histories
reach at most 4 (`b/H1_filing_audit/AUDIT.md`, part C.1).

### A3. What the density is taken from

The rule n_e = 8 abs(F_e) + 8 is a resolution rule. It needs the magnitudes of the current record and nothing
of the field, and an observer who holds the record applies it before the first sample. It does not make the
gate circular. Mutated-field test on the own decoder, which takes the field as an argument and the tape only
for the resolution:

| Mutation of the field (resolution from the tape F) | Read |
| --- | --- |
| M1: field of the reversed tape, -F | the field's degree in 78 of 78, the tape's in 0 |
| M2: one more turn on a boundary edge of the region | the field's degree in 96 of 96, the tape's in 0 |
| M3: 3 abs(F_e) + 3 turns on that edge (increments between pi/2 and pi) | rejected with no degree in 96 of 96 |
| M4: n_e more turns on that edge | the tape's degree in 96 of 96, the field's in 0 |

M4 is the limit of every sampled reading: n_e more turns give the same samples, so a sampled decoder certifies a
degree only under a bound on the derivative, and here that bound comes from the tape. The field's degree is
read whenever the field has at most 2 abs(F_e) + 1 turns on each edge.

The candidate itself cannot be given the test. `read_phase(geometry, current, region)` has no field input: it
builds the samples from `current[e]` (lines 38 to 40) and decodes them. Its output equals sum of s F_e for every
integer tape by A2. No control of `verify_adaptive.py` separates "reads the declared field" from "returns the
sum of the tape": the 168 reads are compared with the sum of the same cochain. The degree clause therefore
stays what the filing text already calls it, an identity for every integer tape, and the amendment changes the
domain on which the sampling is certified and nothing in the evidence. On a one-substitution mutant of the
candidate copy (`replay/adaptive_phase_field_mutant.py`: the phase built from a second current, the resolution
from the argument) the candidate's own arithmetic gives the same table: M1 78 of 78 and M2 96 of 96 read as the
field's degree, M3 96 of 96 `numerical_failure` with no degree, M4 96 of 96 read as the tape's degree.

What the task document has to say so that a maintainer cannot read it as circular:

1. The number of samples on an edge is computed from the received current record before any sample is taken;
   it is a resolution, and the bound of A1 is the reason it suffices.
2. The field is the observer's encoding of the same record. The degree clause is an identity for every
   integer tape, with or without the rule. The rule earns nothing for the cell.
3. One registered control: resolution from the tape, field from the reversed tape; the read is the field's
   degree. A decoder that returns the tape's value fails it.
4. The limit: the read is the degree of any field of the declared form with at most 2 abs(F_e) + 1 turns per
   edge; a field with n_e more turns on an edge gives the same samples.

### A4. Rejections and resource stops

- Fractional input: executed. 0.5, `Fraction(1, 2)`, 4.0, `True`, a string and `None` raise `ValueError`
  ("one actual integer current per seam is required") on the copy of the candidate; no degree is returned. The
  reason is the text of an exception, while a resource stop is a returned record.
- A field that jumps by more than the bound: not established. The function has no field input. Its guard
  stands at pi/2, not at the certified pi/4, and no integer tape reaches it, so no control of the package
  executes it. On the mutant: 3 abs(F_e) + 3 turns give `numerical_failure` with no degree in 96 of 96, under a
  reason text that names finite precision; with one more turn, 94 of 96 reads have an increment at or above
  pi/4 and are returned as `read` (the value is still the field's); n_e more turns are returned as `read` with
  the tape's degree in 96 of 96.
- A field that does not glue: the amendment has dropped the gluing and closure assertions of the public
  routine (lines 141 and 143) and checks no integrality. With one dual edge given another end vertex it returns
  status `read` and the degree -0.0338; the own decoder reports an open curve (defect 0.212).
- Resource cap: executed. At the default budget of 100000 samples, abs(F_e) = 12494 on one dual cell needs
  99997 samples and is read exactly; 12495 needs 100005 and returns `resource_stop` with a reason and no
  degree; 10^100 stops before any float conversion; a zero budget stops; a negative or boolean budget, 29
  currents, a repeated port and port 12 raise. The budget is an argument of the call, with a default.

### A5. Scope: what the filing has to state if it registers this decoder

Rule: "A changed law, memory, state or decoder is a new candidate, with its empirical ancestry kept visible ...
Moving the pin is a pull request that says what changed" (`m10_theory_canonical.md`, lines 20 to 24). The present
texts (`b/H1_filing_audit/fixed_patches/1_gates.patch`, `b/G1_row1_strengthen/branch_gates_v2/1_gates_v2.patch`)
register the package routine and say the opposite of the amendment in six places. A filing with this decoder
states:

1. Candidate and pin. The decoder is part of the candidate: its file, commit and SHA-256 in the candidate
   table and in § 1.1 ("Added decoder"), in the author's public repository outside the package tree, so that
   the tree `029da7a4...` and both manifests stay as filed. If the first candidate is already pinned, this is a
   new candidate id and a pin move of its own. Ancestry: written on 2026-11-14 after an audit found a history
   of the registered kind with four turns; the 168 reads (seed 711034) and the witness are development outcomes.
2. G4. "Sample counts are fixed before the run on the domain abs(F_e) <= 3: 64 per edge in the candidate"
   becomes the rule n_e = 8 abs(F_e) + 8 with its one-line bound, the budget as a number (100000 samples, on
   one dual cell abs(F_e) <= 12494) and the guard. The fail condition "a tape outside the domain accepted"
   becomes "a tape beyond the budget given a degree".
3. The no-search rule. "No ... sample density ... changes after this filing" names the rule and its three
   constants in place of a count.
4. K6. "A tape of 64 turns on one seam ... rejected by G4's domain check" is false under the amendment (64
   turns are read with 521 samples on that seam). The control becomes a tape beyond the budget.
5. The sampling-domain paragraph, the row "The driver on a witness outside the sampling domain" and the
   required outcome of replay command 6 ("the sampler refuses it"): the witness is read, degree 1, 101 samples.
   `witness_four_turns.json` carries `"expect_G4": "outside_domain_refused"` and `heldout_driver.py` implements
   the refusal; both change, with their SHA-256 in the arms row and a new run of the 171 recorded examples.
6. The gluing clause of G4, the one clause the text says can fail. The amendment evaluates no gluing: a
   non-integer current is refused by a type test before the first sample. The filing names the routine that
   carries the clause (the package routine, or the amendment after repair D-A2), or K1's "no curve closes" has
   no carrier.
7. The four sentences of A3.

The present v2 text needs none of this: it records the refusal of the witness and reports the exact contour
degree. Registering the amendment adds a sampled read on the whole registered domain and no evidence for the
cell.

### Defects of part A, with the smallest repair

| Id | Defect | Material | Smallest repair |
| --- | --- | --- | --- |
| D-A1 | `read_phase` builds the field from the tape; no control can fail for a decoder that returns the tape's sum | for the sentence "reads the field"; not for the bound | give the function a field argument (a callable, the declared field as default) and add control M1 to `verify_adaptive.py`; about ten lines |
| D-A2 | no gluing, closure or integrality check; an unglued field is returned as `read` with a non-integer | if the decoder replaces the package routine in G4 | assemble the contour per dual cell as the package routine does, return `open_curve`, and return the rounded integer only when the sum is within 1e-9 of it |
| D-A3 | the guard is pi/2, the certified bound pi/4; the guard's reason names finite precision; no control reaches it | no | state the guard as the public pi/2 and name the outcome by its cause, or guard at pi/4; add control M3 |
| D-A4 | the budget is a caller's argument; fractional input is an exception while a resource stop is a record | for a filing | fix the number in the task document; return a rejection record |

Supplied in part A, all architecture and no calibration: the circle encoding with its constants 0.31 and 0.23
(declared reading of the tape), the rule's constants 8 and 8, the budget, the guard. Nothing is found by the
decoder that the integer tape does not already hold.

## Part B: the neutral-start record-only pointer join

The construction: every port starts at one filling q; one unit is moved across a named seam by a recorded
intervention; then native events under the unchanged endpoint rule, the closing sweep of the thirty seams with
contrast reads, counters modulo six, and the pointer of the public package commanded with the decoded degree.

### B1. What the observer holds at admission: executed, with two riders

Inputs of `admit(filling, tape, region, pointer, retries)`, by reading the code and by a trace of the copy on
logging dictionaries (six tapes; keys read: `records` and `sweep_order` of the source tape; `kind`, `seam`,
`coin`, `current`, `contrast` of its records; `records`, `success`, `reason`, `retained_signed_cycle_steps`,
`winding`, `fine_history_probability` of the pointer tape and its rebuilt `final_loads`; `stage`, `seam`, `coin`,
`current`, `requested_cycle_direction` of its records):

| Input | Content |
| --- | --- |
| `filling` | one integer q >= 1 |
| `tape['records'][0]` | the preparation: seam and sign of the moved unit |
| `tape['records'][1:]` | per native event: seam, coin, current; on the last thirty also the contrast |
| `tape['sweep_order']` | a permutation of the thirty seams |
| `region` | a set of ports |
| `pointer` | records (stage, seam, coin, current, requested direction), success, reason, signed steps, winding, weight; its `final_loads` is removed before any read and rebuilt from the records for each of the twelve possible columns |
| `retries` | the retry bound, default 2 |
| fixed | the seam list, the endpoint rule as a predictor of each current, the pointer's policy, the table 0, 1, 5 |

No argument is a per-port array and none is a terminal load.

Mutations, own join (own observer, own pointer state machine) and the copy of the candidate:

- Terminal loads changed, records unchanged. Altered terminal snapshots of source and pointer, given to the
  copy of the candidate's `admit`: the identical verdict on 288 of 288 recorded tapes. A real transfer of one
  unit after the last record: the records are equal and the verdict is equal in 400 of 400 own histories, and
  the admitted charge then differs from the charge of the terminal loads in 202 of them. The reading is a
  property of the records; "nothing moves after the last record" is a premise.
- A transfer with no record before the sweep: rejected in 400 of 400 by the own join and in 400 of 400 by the
  candidate (the closing contrasts disagree with the reconstruction). This control can fail and does not.
- One field of one source record changed, 7,860 mutations on the 8 recorded source tapes, each with its three
  regions: 7,032 rejected; 150 admitted with the same charge and another weight (a prefix event with no current removed
  or doubled); 678 with the identical verdict, each of them a history that the rule itself produces from the
  same preparation (the coin of an event on an even difference, 615; the seam of a prefix event with no
  current, 63); none admitted with another charge. Every one of the 696 other preparation seams and the
  reversed sign are rejected.
- One field of one pointer record changed, 8,896 mutations on 96 tapes: 7,288 rejected, 1,608 identical, all
  of them the coin of an attempt on a seam without the load, another fine history of the same weight.
- A completed pointer tape of another command: rejected in 59 of 59.
- The candidate's admission and the own one give the same verdict on 13,432 of 13,432 mutated tapes (7,860
  source mutations and 5,572 pointer mutations), and on 400 of 400 own histories; the candidate's source
  executor and the public pointer executor produce the reviewer's tapes in 400 of 400.

Rider 1. The value of the filling does not enter. Every other positive integer gives the identical verdict on
288 of 288 recorded tapes (own join with 1, 2, 9 and 10^6; the candidate's `admit` with q + 1 and 10^6). The
premise in force is that all twelve
ports hold the same load before the pair is moved, so that the total is twelve times an integer.

Rider 2. The observer rebuilds all twelve loads at every step, from `[q] * 12`, the preparation record and the
tape, and predicts every current. It holds a complete per-port start; the start is supplied by the uniformity
premise and one record, not by an array. What makes the reading independent of a remembered start is a fact
the candidate does not execute: the thirty closing records alone (contrasts and currents) fix the terminal
loads up to their mean, and the regional charge is the deviation of the region's load from that mean. An own
decoder that reads only those thirty records gives the admitted charge on 288 of 288 recorded tapes and 400 of
400 own histories, with no filling, no preparation record and no prefix.

### B2. The proposition: proved

Proposition. Let j(n_a, n_b, c) = [d - (2c - 1)(abs(d) mod 2)] / 2 with d = n_a - n_b be the endpoint rule on a
seam (a, b) with coin c. (i) j depends on d and c only, so a uniform shift of all loads changes no current.
(ii) Let S_q be the 133 configurations q + delta with delta = 0, or one +1 and one -1 on two different ports. Every
event on every seam with either coin maps S_q to S_q with a current in {-1, 0, +1}: on a seam with d = 0 nothing
moves; with abs(d) = 1 the unit (excess or deficit) crosses or holds, one coin each; with abs(d) = 2 the adjacent
pair annihilates under both coins. (iii) For n in S_q and every region R, sum_R (n_v - q) is -1, 0 or +1, and the
counter modulo six of the record decodes it.

Proof: (i) is the formula. (ii) and (iii) are finite and are enumerated by own code: 7,980 of 7,980 branches
closed, with unit current and shift covariance, at each of the fillings 1, 7, 100000, 2, 3 and 12 (at filling
1: 5,460 branches with no difference, 1,200 moves, 1,200 holds, 120 annihilations); 544,768 of 544,768 cuts
bounded and decoded. The 60 prepared states reach exactly the 133 configurations. The uniform state is
absorbing, and since j is nonzero for some coin whenever d is nonzero, the uniform fillings are the only
absorbing states. Controls: under a half-current law 2,400 branches leave the integer sector; under j = d, 60
seam visits carry a current of two; two prepared pairs give a regional charge of two, outside the decoder.

### B3. The counts: executed

- 23,940 = 3 fillings x 133 configurations x 30 seams x 2 coins, recomputed in full with an own rule.
- 288 = 8 source histories x 3 regions x 12 pointer columns. Replayed in full with own code: source tapes
  288 of 288, terminal loads 288 of 288, charge from the records equal to the charge of the terminal loads
  288 of 288, verdicts, compatible columns and weights 288 of 288, pointer tapes regenerated from their coins
  288 of 288. `verify_neutral_record.py`, run from a copy, reproduces the 288 tapes byte for byte.
- Failures. 106 of the 288 tapes are timeouts (78 in routing, 28 in measurement). Each is admitted with no
  readout and with the weight 2^-N of its complete record, N the number of source and pointer coin attempts.
  Nothing is conditioned on success. Against the exact pointer law at retry bound 2 (own state recursion,
  180 of 180 laws sum to one and equal the closed formula at retry bounds 1, 2 and 3) the 288 tapes hold 182
  successes where 182.02 are expected (standard deviation 6.75), 78 routing timeouts for 76.80 and 28
  measurement timeouts for 28.90.
- What the 288 tapes cover. Eight source histories. 25 of their 28 prefix events carry no current. The pair
  is present at the end in six and annihilated in two. A region with nonzero charge occurs in two histories
  (cases 0 and 1, 48 tapes); the complete nonzero readouts are 11 (six of -1, five of +1). The comment in the
  script that each history carries every regional value does not hold for six of the eight.
- No complete law of the joined instrument is stored. The README says so. The stored complete laws of the
  package are at retry bound 1 (`AUDIT.md`, D.1 row 6); the tapes run at 2. Own exact numbers at retry bound 2:
  unconditional success between 0.459 and 0.934 for a commanded zero and between 0.082 and 0.166 for one turn.
- The sweep disturbs what it reads. With no event before the sweep the pair survives it with probability
  between 0 and 227/256 = 0.887 over the 60 preparations, mean 0.572; four preparations (seams 7 and 9, both
  signs) annihilate with certainty. The reported charge is that of the latch after the sweep.

### B4. What remains supplied

| Piece | Class | Note |
| --- | --- | --- |
| The uniform start (the "scalar filling") | preparation, architecture | its value does not enter; the premise is uniformity and a total of twelve times an integer |
| The prepared pair | preparation by an external recorded intervention | the rule creates no pair; the uniform state is absorbing |
| The pointer carrier, its controller, route, cycle and retry bound 2 | architecture (apparatus) | displays the observer's decoded degree; no second measurement of the source |
| The integer alphabet and the endpoint rule | architecture (law), unchanged | |
| The closing sweep and its faithful channel | architecture and premise | a replaced closing tape is not authenticated |
| The scheduler of prefix events | architecture | |
| The dictionary charge = load minus filling = minus degree | declared reading | |
| The counters modulo six and the region | architecture (memory), observer's choice | |

No number is set from a target: no calibration. No strength enters.

### B5. The objection

Paragraph for the filing's risk file:

> The record-only neutral arm (`native-neutral-pair-record-only-v1`, manifest `85f9eefc...`) moves the objection
> and removes half of it. The reference is no longer a remembered per-port start. Every port starts at one
> filling, a recorded intervention moves one unit across a named seam, and the charge of a region is its load
> minus the filling. An independent review (`reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B/`) executes three facts.
> The value of the filling does not enter the reading: every positive integer gives the same verdict on 288 of
> 288 tapes. The thirty closing records alone give the regional charge on 688 of 688 tapes, so the reading is
> the deviation of the present loads from their mean and needs no memory of the start. The -1 is a deficit
> that moves and annihilates under the same rule, not the place a unit has left. What stays supplied: neutral
> means uniform, with a total of twelve times an integer (the uniform fillings are the only states the rule
> leaves unchanged); the pair is put in by hand, because the rule creates none; the charges of a pair sum to
> zero and lie in {-1, 0, +1}; the closing channel is faithful and nothing moves after the last record; the
> pointer displays the observer's own register. The closing sweep disturbs what it reads: with no earlier
> event the pair survives it with probability between 0 and 0.887 by preparation. Under a strict reading the
> caveat changes from "relative to a retained preparation" to "relative to the uniform filling, with a
> prepared pair and zero net charge", and the supplied preparation remains a reason to withhold validated.

### Defects of part B, with the smallest repair

| Id | Defect | Material | Smallest repair |
| --- | --- | --- | --- |
| D-B1 | "one true scalar neutral filling" is listed as an input and a retained premise; its value does not enter | wording of a filing | state the premise as uniformity with a total of twelve times an integer |
| D-B2 | the fact that answers the objection, the charge from the closing records alone, is not executed by the package | for the claim "no remembered start" | add the closing-only decoder as a gate; about twenty lines (`closing_only_charge` in `audit.py`) |
| D-B3 | coverage: nonzero regional charge in 2 of 8 histories, 11 complete nonzero readouts, 25 of 28 prefix events without current; the script's comment says otherwise | for the evidence, not for the proposition | choose the regions from the reconstruction (excess alone, deficit alone, both, neither) and the prefix seams at the pair; report the populated cases |
| D-B4 | no complete law at the executed setting; the disturbance by the sweep is not stated | for a filing | store the 36 exact pointer laws at retry bound 2 and the exact sector law of the sweep |
| D-B5 | the two "hidden" controls pass for any code that does not read a key named `final_loads` | no | add the transfer with no record, before the sweep (rejected) and after it (same verdict) |

None of the five contradicts the announced claim: no per-port array and no terminal load enter admission; the
scalar filling and the prepared pair remain.

## Reading for the platform's maintainer

This paragraph is a reading and no failed gate. Under the stricter note on OPH's column, part A earns nothing:
the sampled degree is a float round trip of the integer tape through a declared encoding, with or without the
rule. Part B is an observable computed on the carrier's own process with the preparation, the dictionary and
the apparatus supplied, which the column's sentence makes a case for partial with the caveat in the cell; the
integer is the declared alphabet, the pair is prepared and the net charge is zero. Nothing in either part
follows from installed standard physics, and nothing in either part derives the source.

## Reproduce

    cd openwave-push/reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B
    PY=reverse-engineering-reality/.venv/bin/python   (3.12.6; OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1)
    $PY audit.py own                          # RUN_own.log, RESULT.json: 60 of 60 own checks
    (cd replay; $PY verify_adaptive.py; $PY verify_neutral_record.py; $PY compare_candidate.py)   # RUN_replay.log
    $PY audit.py final                        # RUN.log: embeds replay/COMPARE.json in RESULT.json: 76 of 76 checks

`THRESHOLDS.md` was written at 16:49:48, before the first run. `RUN_own_stopped_1708.log` is a first run stopped
by the reviewer: its pointer laws at retry bound 16 exceeded the CPU budget; the bound was replaced by 3.

## Other parties, read after the verdict

The verdict above was written at 17:29:18 (`VERDICT_BEFORE_OTHER_REVIEWS.txt` holds the hashes of `REVIEW.md`,
`REVIEW.json` and `RESULT.json` at that time). Then read: the root's opposing evidence
`openwave-tuning/next/OBSTRUCTION_REPAIR/hadrons/row1_review/RESULT.json` (`00880e87...`, equal to the fifteenth
source pin of the manifest) with its `review.py` (`f7003384...`), both in `consumed/other_reviews/`. Under
`row_completion/reviews/` no review of this candidate exists; `ROW_01_ROOT` and `ROW_01_DYN17` are of the earlier
row-1 candidate.

Agreements:

- The bound. Its script reduces it to the same inequality (0.62 < 1.54, from pi > 1). Its 120 reads up to
  abs(F_e) = 512 have increments below pi/4 and error at most 4.5e-13; the 192 own reads and the read at 12494
  agree.
- The resource stop, the weights 2^-N on compatible columns, the ignored `final_loads` keys, and the equality
  of the charge from the records with the charge of the terminal loads (its 480 reads at fillings 1 and 19,
  with 47 prefix events and a shuffled sweep order; here 288 recorded and 400 own histories).
- The scope: the event records, the alphabet and the field dictionary are declared, and no physical charge
  identification follows.

Differences:

- Its script imports `adaptive_phase` and `neutral_record` and calls them. Its phase reads compare the
  candidate's output with the sum of the same cochain, the identity of A3. It has no decoder with a field
  input and no mutated field. This review adds both, and with them D-A1 and D-A2. Its verdict is one flag,
  `component_review_pass: true`; here two statements of part A are not established (A3b, A4b).
- It passes `sample_budget=200000`: the cap is whatever the caller gives (D-A4). It executes no fractional
  input, no unglued field and no guard.
- It lists the "scalar preparation filling" among what remains declared. Here the value of the filling is shown
  not to enter (288 of 288 on the candidate's own `admit`), the premise is uniformity, and the thirty closing
  records alone give the charge (688 of 688).
- It checks weights on 42 of the 288 tapes and no law. Here all 288 are replayed by own code and compared with
  the exact pointer law at retry bound 2 (182 successes, 182.02 expected), and the coverage of the tapes
  (nonzero charge in 2 of 8 histories, 11 nonzero readouts) and the disturbance by the sweep are counted.
- It does not address what a filing has to state (A5) or the objection (B5).
