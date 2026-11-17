# RA2: independent component review of the row-1 interface, version 2

Reviewer: Fable-B / `reviews/ROW_01_INTERFACE_V2_FABLE_B`, 2026-11-14, 17:59 to 18:22 (+07). A delta review against
this reviewer's first review (`reviews/ROW_01_ADAPTIVE_DECODER_FABLE_B/`, manifest
`066b2a351c015de54fefacbd23ad21ce534f7ede5a221ba5ca4439743102c2dd`). The reviewer wrote none of the reviewed
code. `audit.py` imports no code of the package; it loads the reviewer's own first `audit.py` for the functions
written there.

Package: `openwave-tuning/next/OBSTRUCTION_REPAIR/row1_v2/`, `MANIFEST.json` SHA-256
`3c954ae36058bd6f685976e376490f9220601e5584529d19206e3ef7e43492a2` (created 17:53:33). Copied 17:59:09 into
`consumed/package/`: 11 of 11 pins equal the copies, the live files still equal the pins at the time of the
audit (18:09), 18 of 18 source and review pins equal the live files. The first review is byte-identical: its
manifest has the recorded hash and 95 of 95 files it lists have their hashes.

This is a component review. It gives a verdict per statement and none on the literal criterion of row 1. No
status is awarded. The lane's own description (sync index, "Response to B042") was read after the bytes.

| File | SHA-256 |
| --- | --- |
| `phase_reader.py` | `152de9ba4b06414d52a28b0e7095741e3c07bd6edfc4cf4e5d454cbd815433be` |
| `closing_reader.py` | `0a9f3fb66737343a2d805b8e0267c8cd5f5c8b6b5b2f0a8d55470f716b59f96a` |
| `verify.py` | `5c947dc04f1e9665d3d7a5bd65fb55ae0cd74bc59e1f5e7772dd75b0fa10a973` |
| `PROOF.md` | `7adbcb22ab548634c30bf2bfbf9dea0b0bc22492b95f8e73cd3663b9c9714231` |
| `README.md` | `02c977f58f536ef9de3296d4a46564a8ce10341486f0e13f4b52c3b773b63c7b` |
| `RESULT.json` and `run.log` | `c6149d5efbc1d0acd1fc50ef7c1a18481a2ace42a396d3d8dc61b70ccaa84f2b` |
| `EXACT_LAWS.json` | `acf809b0ddd6e5278fd8836c6d3028ff2fb123a58a51dfff40899cd3e73318c5` |
| `EXAMPLES.json` | `8997749bbbafaac214249a39bcfc7f9a32ba76d445545cf374e3c18db678b076` |
| `FINAL_CODE_RECHECK.json` | `0604607c2b247a2b57a02657b04a103ffc745f2db1b2d03641c84b7a214b1ed3` |
| `recheck.log` | `ac547a3e54bfdd0bf7a585c5a9f14d09b26b870d9387167b9adb5d62d113d64c` |

## Verdicts

| Statement | Verdict |
| --- | --- |
| V1 the phase reader takes a supplied field; the tape fixes the resolution only | executed |
| V2 the reading is certified on a stated field class | proved, and executed |
| V3a that class is a property of the declared encoding of the recorded tape | proved |
| V3b that class holds for a field from another producer | not established: an added premise on the source |
| V4 the closing reader takes thirty E records and a region and nothing else | executed |
| V5 an excluded input changes no verdict; one changed record changes the verdict or is rejected | executed |
| V6 the charge on own histories is that of the terminal loads | executed |
| V7 all 2^30 coin histories per start by exact dynamic programming | executed, as a count for one sweep order |
| V8 the joined instrument is complete | executed (assembled here; the file's flag is a constant) |
| V9 the remembered-start objection is removed for this observer | executed, with what remains in section 3 |
| V10 no unrecorded later motion is detected | executed; wider than "later" |
| V11 the present filing texts admit version 2 | not established: passages in section 4 |

Checks: 24 of 24 own, 15 of 15 on copies of the package, and `replay/WITNESS_V2.json`.

## The earlier defects

| Id | Earlier defect | Now | Own test |
| --- | --- | --- | --- |
| D-A1 | the reader built its field from the tape | repaired | `read_phase` takes `field` as a required argument and never builds one. With the resolution from the tape and own fields: reversed tape, one more turn and +-(2 abs(F_e) + 1) turns are read as the field, 48 of 48 each on the own decoder and 40 of 40 each on the package copy |
| D-A2 | no gluing, closure or integrality check | repaired | package copy: an unglued field is `open_curve`, a field off the unit circle `invalid_field`, a fractional sign or an unclosed cell `rejected_geometry`; the degree is an integer. Rider: a closed cell with reversed orientation is accepted and read with the opposite sign; the geometry is a caller's argument and has to be pinned by the task |
| D-A3 | guard at pi/2 against a certified pi/4, misnamed, never executed | repaired | the guard is the edge of the registered class and its outcome is `increment_guard`; own fields between 2 abs(F_e) + 2 and 6 abs(F_e) + 5 turns: 238 of 240 rejected, 2 read as the field, 0 wrong |
| D-A4 | the budget a caller's argument; fractional input an exception | repaired | package copy: 100001 and 10^9 are `rejected_input`, 12494 turns read, 12495 and 10^100 stop before the field is called; a non-integer record entry returns the record `rejected_input`. Rider: an exception of the supplied field outside five types (a `ZeroDivisionError`) leaves the reader as an exception |
| D-B1 | the filling listed as an input | repaired | the closing reader has no filling. The same relative state placed at fillings 0, 7 and 10^6 gives the same thirty records, 400 of 400 |
| D-B2 | the closing-only gate absent | repaired | `observe_closing` is that gate. Own reader and package copy agree on 2,664 of 2,664 reads; charge equal to that of the terminal loads on 2,664 of 2,664 own reads |
| D-B3 | thin coverage of charged cases | repaired in part | 133 starts, 266 sweeps, 2,394 joined effects. Counted on the replay: 186 of the 2,394 effects are on a charged region (93 of +1, 93 of -1) and 17 are complete nonzero readouts (10 of +1, 7 of -1); 2 discovery, 371 routing and 144 measurement timeouts. `RESULT.json` reports no count of charged cases, and the regions are taken by index, not from the state |
| D-B4 | no complete law at the executed setting; the disturbance not stated | repaired in part | the 133 sweep laws and the 36 pointer laws are stored and equal the own recursions (133 of 133, 36 of 36). They are for the order 7..29, 0..6; the 266 executed sweeps use shuffled orders. The joint flag in the file is a constant (assembled here: 4,788 of 4,788 sums equal one). The survival of the pair is in the file and in no sentence: 0 to 227/256 for the 60 adjacent starts, 1169/4096 to 477/512 for the 72 others |
| D-B5 | the hidden-snapshot controls cannot fail | repaired in part | the source snapshot is no input any more. Motion after the last record is stated as not detected, by one arithmetic line that calls no reader (`verify.py` 201 to 207). Own test: 400 of 400 such transfers change no record, and the admitted charge then differs from the terminal one in 55 |

## 1. The phase reader

Read from the bytes: `read_phase(geometry, resolution_record, region, *, field, sample_budget)`. The record gives
n_e = 8 abs(F_e) + 8 and nothing else. The samples are `field(e, k/n_e)`. Checked before a degree is returned:
integer record, budget at most 100000, closed oriented cells with signs +-1, finite unit samples, gluing at
every corner and at the closing point (1e-9), every principal increment below pi/2, an integral sum per cell
(1e-9). Every return carries `regularity_verified_from_samples: False`.

The class on which the reading is certified, stated exactly:

- No wrong degree. If on every sampled edge the phase lift moves by less than 3 pi/2 between consecutive sample
  points (in particular if sup abs(phi'_e) < 3 pi n_e / 2), the reader returns the field's degree or a rejection.
  Proof: an increment below pi/2 in magnitude is its own principal value, and the sum telescopes; an increment
  between pi/2 and 3 pi/2 has a principal value of magnitude at least pi/2 and meets the guard.
- Always read. If sup abs(phi'_e) < pi n_e / 2, the package's registered class, no increment meets the guard and
  the degree is returned. The declared family with abs(F'_e) <= 2 abs(F_e) + 1 lies in it: the largest increment
  is at most [2 pi (2 abs(F_e) + 1) + 0.62 + 0.46 pi] / (8 abs(F_e) + 8) = pi/2 - (1.54 pi - 0.62)/(8 (abs(F_e) + 1)),
  the inequality of the first review; on the declared geometry 1.000915 at abs(F_e) = 0 and 1.570752 at 12494,
  against pi/2 = 1.570796.
- Outside. A field whose lift moves by a nonzero multiple of 2 pi, give or take less than pi/2, between two
  sample points is returned as `read` with a wrong degree.

Executed on own fields with the own decoder (48 tapes, seed 2026111403, abs(F_e) up to 64; one boundary edge
changed) and on the package copy (40 of these tapes: 520 of 520 reads with the status and the integer of the own decoder):

| Field (resolution from the tape F) | Own decoder |
| --- | --- |
| same tape; reversed tape; one more turn; +-(2 abs(F_e) + 1) turns | read as the field, 48 of 48 each |
| 2 abs(F_e) + 2 turns | 46 rejected, 2 read as the field |
| 3 abs(F_e) + 3, 4 abs(F_e) + 4, 5 abs(F_e) + 5, 6 abs(F_e) + 5 turns | rejected, 48 of 48 each |
| n_e more turns; n_e fewer turns | read as the tape, 48 of 48 each |
| n_e + 1 more turns | read as neither the field nor the tape, 48 of 48 |
| one extra turn inside a single sample interval (own field; the samples are those of the tape's field) | read as the tape, 48 of 48; read as the field at 16 n_e on that edge, 48 of 48 |

The package keeps 72 aliases of the kind "n_e more turns" as its counterexample. The last two rows show that the
limit is wider than that case: a wrong degree need not be the record's, and the extra turn need not be spread
over the edge.

Is the class a property of the declared encoding? For the field that the registered construction has, yes; for
any other field, no. The task document calls the field "the observer's reading of the tape". That field is
`declared_field` of the recorded tape, F' = F, and it lies in the class by proof (largest increment below pi/4,
first review A1). Then the reader cannot return anything but the sum of the tape, and the degree clause stays
an identity for every integer tape. A field from any other producer lies in the class only if that producer
bounds its phase between the sample points, by a relation to the record (at most 2 abs(F_e) + 1 turns per edge).
That is an added premise on the source, and no finite sample set checks it; the package says so. Version 2
repairs the reader. It supplies no field source: nothing in the candidate produces a field apart from the tape.

What the row-1 task document has to say:

1. The reader takes the field as an argument; the record fixes the number of samples and nothing else.
2. In the registered gate the field is the declared encoding of the same tape. It is in the certified class by
   proof, and the degree clause is an identity for every integer tape. The rule earns nothing for the cell.
3. The reading of a field is certified when its phase moves by less than pi/2 between consecutive samples
   (always read), and no wrong degree is returned below 3 pi/2. A field with a full turn hidden between two
   samples is read wrongly; the read is conditional and says so.
4. The controls of the reader, with their counts: reversed field, extra turn, guard, and the retained aliases.
5. The geometry is pinned by SHA-256: the reader checks closure and signs of the cells, not their orientation.

## 2. The closing reader

Read from the bytes: `observe_closing(records, region)`. Thirty records of kind `E`, each with seam, coin,
current and contrast; each seam once; current and contrast in the alphabet of one pair and equal to the endpoint
rule; the thirty difference equations solved on the connected seam graph and checked on every cycle; the mean
removed; the start required to be integral and to hold one pair or none; the sweep replayed with the contrast
and the sector checked at every step. The charge is the sum of the end state over the region.
`admit_closing(records, region, pointer, retries)` adds the pointer of the public package, commanded with minus
the charge, for all twelve columns, with the weight 2^-(30+m).

- Inputs. The arguments of the package copy are the records and the region (`records`, `region`). Records dressed with a
  filling, terminal loads and a preparation give the same read, and the keys read are `kind`, `seam`, `coin`, `current`, `contrast`.
- Excluded inputs changed, own histories. Filling, preparation and prefix: the thirty records of a full history
  (filling, prepared pair, prefix of up to twelve events) equal the records of a sweep from the same relative
  state placed directly at fillings 0, 7 and 10^6 with no preparation and no prefix, 400 of 400. A terminal
  snapshot is no argument; a transfer after the last record changes no record, 400 of 400.
- One record changed: 23,030 mutations on 70 own tapes. 19,454 are rejected. 3,576 give the identical verdict,
  and each of them is a sweep that the rule itself produces from the same start: the coin of an event on an
  even difference (1,711 of 2,100) and two neighbouring records that commute, exchanged (1,865 of 2,030). None
  gives another charge. Current +-1, contrast +-1, another seam, another kind, a missing contrast, a deleted and
  a doubled record are rejected in 2,100 of 2,100 each. The package copy gives the own verdict on 23,030 of 23,030.
- The charge on own histories: 400 histories and 266 sweeps from the 133 relative states with shuffled orders,
  four regions each: 2,664 of 2,664 reads equal the charge of the terminal loads, with the measured start and
  end equal to the true ones (372 reads of +1, 744 of -1, 1,548 of 0). The package copy: 2,664 of 2,664, equal to the own reader.
- Outside the sector (two pairs, one excess with a total of 12 q + 1, a unit of two, a current of 10^50): no
  charge, own reader and package copy.
- The joined effect: on 150 own tapes with own pointer tapes the package copy gives the admission, charge, readout, compatible columns and weight 2^-(30+m) of the own join (150 of 150); an altered pointer terminal load changes nothing (150 of 150); a completed pointer tape of another command is rejected (21 of 21).
- "All 2^30 coin histories per start by exact dynamic programming": the 133 stored laws are counts of coin
  words by recursion over the thirty seams in the order 7..29, 0..6, each summing to 2^30. The own recursion
  gives the same 133 laws. Brute force on a reduced instance: all 1,024 coin words of the first ten seams for
  each of the 133 starts, and all 16,384 words of the first fourteen seams for three starts, give the
  histograms of the recursion (133 of 133, 3 of 3). The neutral start and a pair on the first swept seam end
  neutral in all 2^30 words, as they must. The 36 pointer laws equal the own state recursion. The sentence
  holds as a count; no history is enumerated at depth thirty.
- `verify.py`, run from a copy, reproduces `RESULT.json`, `EXAMPLES.json` and `EXACT_LAWS.json` byte for byte.

## 3. What it does not do, and the remembered-start objection

The package's words: it "does not derive the integer E channel from electromagnetism or detect unrecorded later
motion". Both hold. The second is wider than its wording. A port is read as of its last visit in the sweep, not
as of the thirtieth record: of 400 transfers of one unit with no record between two records of the sweep, 353
are rejected, 44 are read with the charge of the terminal loads and 3 are read with another charge. A transfer
before the first record is read as the state it leaves (109) or rejected because that state is outside the
sector (291), never as another charge. The premise in force is: from the first closing record on, every motion
is one of the thirty recorded events.

The objection of the filing audit (`b/F1_filing/RISKS_ROW1.md`, point 2): "The instrument needs the true
preparation G ... a moved unit reads +1 where it is and -1 where it was. That is a displacement record against a
memory. A charge would be read from the present state through a cut."

For this observer the objection is removed in its first two sentences. The observer holds no preparation, no
filling, no record of the pair and no tape before the sweep; the reference is the mean of the present loads,
measured by the contrasts; that the state holds one pair or none on an integral mean is measured too, and a
state outside is rejected. The -1 is a deficit present on the carrier. What stays, and is not this objection:
the charge is read from the present state through all thirty seams, not through a cut; it is the state after a
sweep that moves and annihilates the pair (survival 0 to 227/256 from adjacent starts); the charges of a pair
sum to zero and a single excess is rejected, so no net charge is read; the pair is prepared outside the
instrument, because the rule creates none; the closing channel is taken as faithful and complete; the integer
alphabet and the dictionary are declared. The caveat of a cell would change from "relative to a retained
preparation" to "relative to the mean of the present loads; the pair is prepared; zero net charge".

## 4. What the row-1 filing would have to change to register version 2

The package says of itself that the observer policy is changed and that this is a new candidate. Passages, by
file and line of the present texts (`consumed/filing/`):

`b/FILING_ROW1_READY.md` (`5ac492f4...`)

1. Lines 9 to 16, the table of files: two more files to publish before the run, `phase_reader.py` and
   `closing_reader.py`, with their SHA-256; the hash of the task document in step 2 changes with the edits below.
2. Lines 32 to 39, item 4: it names `OBSTRUCTION_REPAIR/row1/` and says that registering needs two repairs and
   six text changes. It becomes: version 2 (`3c954ae3...`), the repairs made, this review, and the list below.
3. Line 45: "the charge is relative to a retained preparation or a known filling".

`b/H1_filing_audit/fixed_patches/0_pin_move.patch` (`04f0cde7...`), section 1.1

4. Line 42, Candidate: "a declared observer apparatus that holds the preparation and receives the event records
   ... and thirty closing contrast reads", and the candidate id ending in `record-only-v1`.
5. Line 45, Added memory: "The committed current tape F (one integer per seam, zero at preparation)". The
   closing observer keeps no tape from the preparation.
6. Line 46, Added decoder: the closing reader and the phase reader join the sweep and the pointer.
7. Line 49, Empirical ancestry: the two interfaces were written on 2026-11-14 after an independent review of
   the first version; development seeds 411719 and 411720.

`b/G1_row1_strengthen/branch_gates_v2/1_gates_v2.patch` (`93a2cbd7...`), the task document

8. Line 9, the index row, and lines 62 to 63, the author's sentence on the two arms.
9. Line 71, the candidate id; line 74, the arms directory with the SHA-256 of every file.
10. Line 130, A field: produced by the declared encoding from the tape and handed to a reader that builds none;
    the class of section 1.
11. Line 134, The reference; line 135, The observer's memory and decoder; line 136, The closing channel (the
    whole reading now rests on it: add the premise of section 3).
12. Line 137, The neutral start: "a reader of the loads at each reading ... It is no record-only observer".
13. Line 157, Complete pointer laws: "No complete law is on record at the settings of the 171 examples (retry
    bound 2)". Version 2 stores 36 pointer laws at retry bound 2 and 133 sweep laws for one order.
14. Line 159, the witness row; line 277, the required outcome of command 6; lines 337 to 342, the sampling
    domain: the four-turn tape is inside the budget and is read, no longer refused (package copy: degree 1, 101 samples; `replay/WITNESS_V2.json`).
15. Line 212, G4: the fixed counts of 64 and 32 on abs(F_e) <= 3 become the rule, the budget and the guard. Its
    clause "a non-integer current leaves the curve open" needs a producer: `declared_field` refuses a non-integer
    current (`ValueError`), so control K1 has to hand the half-turn field to the reader itself (from another producer it is `open_curve`, defect 2.0).
16. Line 213, G5: the reconstruction from the received tape becomes the read from the thirty closing records.
17. Line 216, G8: "the excess and the deficit being read from the loads at that reading".
18. Line 219, the no-search rule: "sample density".
19. Line 238, K6: 64 turns are read (521 samples on that seam); the control becomes a tape beyond the budget.
20. Line 242, K10: controls of the closing reader (one record changed, a state outside the sector).
21. Lines 276 to 277, the replay commands; lines 346 to 353, the drivers.
22. Lines 363 to 369, the first stated limit: "it reads the loads at every reading, so it is no record-only
    observer ... a record-only form of the neutral start is not part of this task".

`b/F1_filing/RISKS_ROW1.md` (`f9b0bf16...`), point 2: the rows "What the packet answers", "Strict reading" and
"Smallest construction that closes it".

Nothing has to change if version 2 is not registered: the present task document registers the refusal of the
four-turn tape and the load-reading neutral arm and is consistent as it is.

## New findings on version 2, with the smallest repair

| Id | Finding | Material | Smallest repair |
| --- | --- | --- | --- |
| N1 | a port is read as of its last visit; motion with no record inside the sweep is read with another charge in 3 of 400 | for the sentence "reads the state at the last recorded event" | state the premise of section 3; or add a second pass over the thirty seams and require the two reads to agree, and state the added disturbance |
| N2 | the stored sweep laws are for the order 7..29, 0..6; the reader accepts any order and the 266 executed sweeps are shuffled | for a filing | name the registered order in the task document, or have the reader check it |
| N3 | `all_joint_effects_sum_to_identity` is written as a constant | no | compute it (done here from the stored laws) |
| N4 | `RESULT.json` reports no count of charged reads or of complete nonzero readouts | for the evidence | report the tally (section on D-B3) |
| N5 | the geometry's orientation is not checked by the reader | no | pin `GEOMETRY.json` by SHA-256 in the task |
| N6 | the control for motion after the last record is one arithmetic line | no | run the transfer through the source and compare records and charge |

Supplied, unchanged in class from the first review: the circle encoding and its constants, the rule's constants
and the budget (architecture, decoder); the single-pair sector as the domain of the instrument, the sweep order,
the pointer with its controller and retry bound 2 (architecture, apparatus); the integer alphabet and the
endpoint rule (law); the faithful and complete closing channel (premise); the prepared pair (preparation,
outside the instrument). No number is set from a target.

## Reading for the platform's maintainer

This paragraph is a reading and no failed gate. The phase reading stays a round trip of the integer tape
through a declared encoding; version 2 makes the reader able to fail on a field it is handed and adds no field
source. The closing observer computes its charge on the carrier's own records with no remembered start; the
integer is the declared alphabet, the pair is prepared, the net charge is zero, the channel is taken as faithful.
Under the column's sentence that is a case for partial with the caveats in the cell, as before.

## The lane's description, read after the bytes

"Response to B042" in the sync index (`consumed/sync/`, `ab54f4f0...`). Agreements: every count it gives is
reproduced by the replay (72 reversed reads, 57 nonzero comparisons, 72 extra turns, 72 guards, 72 aliases, 266
sweeps, 798 reads, 2,394 effects, 36 laws); "the tape fixes resolution only"; "finite samples do NOT
authenticate off-grid regularity"; "no filling, preparation, prefix or terminal snapshot". Differences: "this
removes the remembered-start dependency from this observer" holds for the inputs, with the remainder of section
3; "detect unrecorded later motion" is wider than "later" (N1); "exact dynamic programming retains all 2^30 coin
histories per start" holds for one sweep order (N2).

## Reproduce

    cd openwave-push/reviews/ROW_01_INTERFACE_V2_FABLE_B
    PY=reverse-engineering-reality/.venv/bin/python   (3.12.6; thread variables 1, PYTHONDONTWRITEBYTECODE=1)
    $PY audit.py own                        # RUN_own.log, RESULT.json: 24 of 24 own checks
    (cd replay; $PY compare_v2.py)          # RUN_replay.log, COMPARE.json: copies of the package against own code
    $PY audit.py final                      # RUN.log: embeds replay/COMPARE.json in RESULT.json

`THRESHOLDS.md` was written at 18:06:43, before the first run.
