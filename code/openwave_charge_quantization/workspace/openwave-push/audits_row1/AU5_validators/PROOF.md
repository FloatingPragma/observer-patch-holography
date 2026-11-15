# The closing contrast sweep: statement, proof, hypotheses

Frozen callable: `closing_sweep.closing_sweep(G, retained_tape, source_state, order, coins)`.
The argument is the one of ledger entry A033 (`openwave-tuning/sync/AUDITOR.md`), written out with its hypotheses.

## 1. Setting

Ports V (twelve), seams E (thirty), each with a fixed orientation a → b. B is the incidence matrix (+1 at a, −1 at
b). A native event on a seam with coin c reads the contrast d = n_a − n_b, transfers the current

    j(d, c) = (d − (2c − 1)(|d| mod 2)) // 2,        n_a −= j,   n_b += j,

and emits the record (seam, coin, contrast, current). The current is a function of the contrast and the coin.

The observer knows the preparation G and holds a retained tape of current records with per-seam sums F̂. The
true per-seam sums are F. Loads of the source when the sweep starts (latch L0): n = G − B F. The observer's
reconstruction: n̂ = G − B F̂. Discrepancy: δ = n − n̂ = B (F̂ − F). Its entries sum to zero, because every column
of B sums to zero.

The sweep appends events on seams e_1, …, e_m with coins c_1, …, c_m. At step k the source performs the native
event on e_k = (a_k, b_k) and emits (c_k, d_k, j_k). The observer compares d_k with the contrast of its
reconstruction, d̂_k = n̂_{a_k} − n̂_{b_k}, and then advances the reconstruction by the received current j_k.

## 2. Statement

**Hypotheses.**

- H1 (connected graph, every seam swept). The swept seams connect all ports. The thirty seams of the carrier do.
- H2 (the sweep's own records are received faithfully). The contrast, coin and current of each of the m sweep
  events reach the observer unchanged.
- H3 (known G). The reconstruction is G − B F̂ with the true preparation G.
- H4 (total load conserved). The source changes only by seam currents, so n = G − B F for an integer F and the
  totals of n and n̂ are equal.

**Theorem.** Under H1 to H4 the sweep shows no contrast mismatch if and only if n̂ = n at L0. In that case
n̂ = n after every sweep event, in particular at L1.

**Proof.**

1. *Invariance.* A sweep event changes n by −j_k at a_k and +j_k at b_k. The observer applies the same received
   current to n̂ (H2). So δ = n − n̂ is the same before and after every sweep event. (If the observer instead
   advanced n̂ by its own predicted current, the two currents would agree until the first mismatch, because the
   current is a function of contrast and coin; δ is then invariant until the first mismatch. The first mismatch
   is the same in both conventions.)
2. *What one comparison says.* d_k − d̂_k = (n_{a_k} − n_{b_k}) − (n̂_{a_k} − n̂_{b_k}) = δ_{a_k} − δ_{b_k}. The
   record of step k matches if and only if δ takes the same value on the two ports of e_k.
3. *All records match.* Then δ_a = δ_b across every swept seam; by H1 δ is constant on V; by H3 and H4 its
   entries sum to zero; twelve times the constant is zero; δ = 0.
4. *Conversely* δ = 0 gives d_k = d̂_k at every step by (2), and δ stays 0 by (1). ∎

## 3. Consequences

- **C1. The mismatching seams are fixed by δ alone.** They are the swept seams (a, b) with δ_a ≠ δ_b. They do
  not depend on the order of the sweep or on any coin. Detection of δ ≠ 0 is certain.
- **C2. One lost moving record.** A record with current j = ±1 on seam (a, b) is lost: δ = −j (e_a − e_b).
  The mismatching seams are (a, b) and the four further seams at a and at b: 9 of 30 on this carrier. The first
  mismatch is the first swept seam that touches a or b.
- **C3. Signed reconstructions.** A damaged tape can give a reconstruction with a negative entry (source e_w,
  reconstruction e_w + e_u − e_v). The proof uses only that the current is a function of contrast and coin.
- **C4. The latches.** The sweep is thirty native events on the source. The loads at L1 are the loads at L0
  advanced by the thirty sweep currents. A pass certifies the reconstruction at L0, and with the sweep currents
  appended, at L1. A regional charge decoded from the retained counters refers to L0, a state the source has
  left; counters advanced by the thirty sweep currents describe L1. Both states and all thirty records are
  returned by the callable.

## 4. What the sweep does not detect

- **Lost records whose currents cancel.** If B (F̂ − F) = 0 (a hop and its return on one seam both lost; lost
  unit currents round a closed loop), then δ = 0: no mismatch. The reconstruction of the loads is correct; every
  regional charge is correct. The tape differs from the true tape by a closed current.
- **A lost record plus a corrupted sweep record.** If H2 fails, a pass certifies nothing. Replacing each
  received contrast by the reconstruction's own prediction produces a pass with δ ≠ 0. The sweep replaces the
  premise "every source record was received" by the premise "the thirty sweep records were received".
- **A discrepancy that is constant on every component of the swept graph.** If the swept seams do not connect
  all ports (H1 fails), a pass certifies only that δ is constant on each component. Sweeping the 21 seams that
  avoid both ports of the lost seam shows no mismatch. A spanning tree of 11 seams is enough for H1.
- **A constant discrepancy.** δ = c on all ports is invisible to every contrast. H3 and H4 exclude it; lost seam
  records cannot produce it.

## 5. Status of the instrument

The sweep is an appended instrument: thirty further native events on the source and thirty contrast records.
It is not a retrospective completeness property of the current-only tape. From current records alone a tape
with a lost record cannot be told from a complete one whenever the shortened tape is itself a native history
from G: such claims pass the strict tape validator of this directory (Section 6). Without the sweep the loss is
visible only against the terminal load, which the controller does not hold.

## 6. Executed verification (`run.py`, `RESULTS.json`)

**One load.** Every history of one to three events from every start port over the sixty (seam, coin) letters:
2,635,920 histories with 655,260 moving records, hence 655,260 damaged tapes (one moving record removed). They
fall into 660 distinct pairs (source loads, reconstruction); the 7,260 damaged tapes of histories with up to two
events were passed literally to `reconstruct`. For each pair and each of the twelve orders (identity, reverse,
seeded orders 1 to 10) every one of the 2^30 coin words was covered by a transfer recursion through `sweep_step`
that merges words with equal (source, reconstruction, mismatch positions): 7,920 sweeps, 8,504,035,246,080 coin
words, **none without a mismatch**. Every sweep has exactly nine mismatching seams, the same nine for every coin
word; the first mismatch lies between positions 0 and 21; δ is unchanged in every leaf.

*Restriction on the coin words: none.* An even contrast gives the same contrast and current for both coins;
those two words are carried together. The recursion was compared with a literal enumeration without merging: 51
sweeps, 3,852 words of the decision coins (the coins of events with an odd contrast) executed through the frozen
callable, each standing for the 2^(30 − decisions) words that share it; the weights sum to 2^30 per sweep and
the verdicts agree.

*No false alarm.* The twelve end states of complete tapes, twelve orders, all 2^30 coin words (154,618,822,656
words): no mismatch. These states are the ends of all 2,635,932 complete histories (the 2,635,920 above and the
twelve empty ones).

**Two loads, reduced set.** The four preparations with units on {0,0}, {0,1}, {0,2}, {0,3} (doubly loaded,
adjacent, distance two, antipodal), every history of one or two events: 14,640 histories, 4,578 damaged tapes
(all passed literally to `reconstruct`), 389 distinct pairs; 4,668 sweeps, 5,012,226,834,432 coin words, none
without a mismatch, nine mismatching seams in every sweep. No false alarm on the 61 end states and twelve orders
(785,979,015,168 words; 14,644 complete histories).

**The frozen callable on explicit coin words.** `closing_sweep` was called on one real damaged tape per distinct
pair and order: 7,920 calls (one load) and 4,668 calls (two loads). Verdict, mismatch positions (C1),
reconstruction at L0, and the equality of source minus reconstruction at L0 and L1 agree in all calls.

**Not detected (Section 4), executed.** 10,560 tapes with two lost records whose currents cancel (one load, up
to three events): reconstruction equal to the source and verdict `consistent` on all 144 state and order pairs.
Sweep records replaced by the reconstruction's own prediction: 1,320 sweeps, all pass although the loads differ.
One sweep record with its current changed by one: 1,320 sweeps, all flagged `sweep_record_not_native`. Sweep of
the 21 seams that avoid both ports of the lost seam: 1,320 sweeps, no mismatch, verdict
`inconclusive_sweep_not_spanning`. Sweep of an 11-seam spanning tree: 1,320 sweeps, all detected.

**Change of the source by the sweep.** One load: the source ends the sweep where it started with probability
between 1/32 and 25/128, depending on port and order (1/32 is the probability that no sweep event moves the
unit). Two loads: between 0 and 37/512 (0 for a doubly loaded source: the first incident sweep event divides it).

**With the tape validator.** 1,680 claims with one dropped moving source record, recomputed from the shortened
tape, pass `strict_tape_validator.validate`. The sweep on the actual source detects every one of them in every
order: 20,160 sweeps, nine mismatching seams each.

**The coordinator's sampled check** (`bridgework/COORD_checks/completeness_sweep.py`: 1,983 of 1,983, 2,925 of
2,925, 3,389 of 3,389 lost records detected; no false alarm in 4,000 each) is the early-stopping form of the same
comparison. By C1 its outcome is certain for every history, order and coin word. Three loads were not enumerated
here.

Not executed: any configuration of total load four.
