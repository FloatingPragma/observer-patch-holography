# Row 11 status: ready for full-row review

Criterion (SHA256 49297c0d…): "exactly 3 charged-lepton minima; ratios 1 : 207 : 3477". Candidate
`row11-generation-register-into-higgs-wilson-charged-source-fable-v1`. No OpenWave validation is reported; the
official M10 cell is untested. The producer does not set `full_criterion_met`; see the question below.

## Clause table

| Clause | Gate (fixed 15:33:27, before the runs) | Value | Result | Declared, with class |
| --- | --- | --- | --- | --- |
| "exactly 3 charged-lepton minima" | G1a: every face stabilizer of order 3, regular on three corners | 20 of 20; 60 rotations | pass | carrier: architecture |
| | G1b: commutant of the third-turn | dimension 3; Hermitian real dimension 3 | pass | |
| | G1c: distinct positive levels at (1/2, 2/9) and at 2²⁰ rounds | 3 (2.379438, 0.040350, 0.580212) | pass | the third-turn register P1: architecture; selected among the carrier's {2, 3, 5} |
| | G1d: one minimum per coefficient on the Higgs-Wilson symbol (side-four cube; gradient and Hessian exact) | [1, 1, 1], three distinct rest values, total 3; frozen census [1, 1, 1] | pass | Higgs-Wilson term, coefficient one: architecture; padding label excluded by declaration |
| "ratios 1 : 207 : 3477" | G2a: \|r_μ − 207\| ≤ 1/2 and \|r_τ − 3477\| ≤ 1/2 | limit 206.770316 and 3477.4728 (0.2297, 0.4728); 2²⁰ rounds 206.769530 and 3477.4589 | pass | analyzer (balance 1/2); exclusive read, activity 1/3, unit kick, many rounds (phase 2/9): see classification |
| | G2b: own code against the repository strings and the frozen midpoints | 1.9e-25 relative; 1.3e-41 absolute | pass | |
| | G2c (reported, not gated): measured values, CODATA 2022 and PDG 2024 | muon +9.83 ppm = 442 standard deviations; tau +0.61 standard deviations of the tau mass | tau within 3, muon outside | electron mass: calibration, scale only |

Controls, all as required. Count: phase 0 gives 2 levels; no charged part 1; a square face 4 (the carrier has none
and no rotation of order 4); seam pairs with a fifth-turn or a half-turn give registers of 5 or 2 states; the
padding label counted gives 4 minima. Ratios: phase 0 (1, 67.94); plain record (1.895, 2.574); one round (19.06,
278.1); native stationary count (88.53, 1628.1); phase + 0.001 (218.85, 3662.6); one load per cell (8.287, 245.1).
A pair fitted to the measured ratios passes G2a and is rejected as not a menu value.

## Classification of the declared pieces (group rule)

| Piece | Fixes | Class | Ancestry on file | Menu number |
| --- | --- | --- | --- | --- |
| carrier, coherent seam law | the split agree : repair + : repair − = 2 : 1 : 1 | architecture | none | |
| two seams with a third-turn, register F, coupling, oriented record | the count three | architecture | the count was known | 3 options {2, 3, 5}; 120 of 420 seam pairs; 60 of 120 at one port |
| analyzer flags | balance 1/2, Q = 2/3 | question below | "target-informed balance premise" (Koide paper); FZ-10, 2026-07-28: "first abstracted from the measured triple" | 3 of 14 analyzers |
| exclusive read on two independent objects | count 2/9 | question below | phase architecture "selected after exposure to the charged spectrum", 2026-07-13; ledger: "retrospective phase 2/9" | 4 of 57 reads |
| activity 1/3 | count 2/9 | preparation | declared in P2 | 1 of 45 activities |
| one radian per unit count, limit of many rounds | phase 2/9 radian | question below | "an additional rule" (particle paper) | 1 of 6 coefficients; 0.2053 at one round |
| squared reading, order assignment | ratios from levels | architecture (convention) | Koide's form | not enumerated |
| Higgs-Wilson coefficient one | one minimum per coefficient | architecture | "No muon or tau target is used" | |
| electron mass | MeV values | calibration | CODATA 2022 | enters no ratio |

All three kinds together: 82 of 35,910 labelled choices (1 of 6,018 distinct pairs) give both ratios within 1 %,
all 82 the pair (1/2, 2/9), which ranks first; the next pair is at d = 0.136. At the declared activity: 12 of 798.
Unit-fraction menu: 1 of 654. Source-listed menu: 1 of 88. Within 3 standard deviations on both ratios: none.
With a wider bound on the activity (denominators to 20 or 30; fixed after the run) a second pair lies 0.3 % away;
within 100 ppm the chosen pair stays alone and first.

## The question for the reviewers

Are the balance 1/2 and the phase 2/9 architecture, or "a number set from the target it is compared with"?

Producer's answer: calibration by discrete selection. The repository's own ledger dates both numbers before the
register (2026-07-13 and 2026-07-28 against 2026-11-11): the balance "first abstracted from the measured triple",
the phase architecture "selected after exposure to the charged spectrum". The
menu shows that they are a selection and not a tuning: one pair of 6,018, reached by 82 of 35,910 choices, the
nearest pair to the one solved from the measured ratios, with the next pair 14.6 % away in the muon ratio. In the
limit of many rounds the construction has no continuous parameter. What the row has is a retrodiction with a
discrete selection: the count three proved on a declared register, two ratios executed from two selected numbers,
the tau ratio within 0.61 standard deviations, the muon ratio 442 standard deviations off.

Independent of the target: three levels of a third-turn register; one minimum per coefficient; Q = 2/3 if and only
if the balance is 1/√2; balance 1/2 at every activity for the three analyzers that read the event's 2 : 1 : 1
split; the electron anchor.

## Verdict

Ready for full-row review. No computation is missing for the decision. Two rulings decide the count.
1. Ancestry. If calibration (the producer's answer), clause 2 is accommodated and the row is not a derived full
   criterion. The remaining obligation is then exact: a statement from the source, with no mass on file, that
   fixes (a) the third-turn register, (b) the analyzer among the 14, (c) the exclusive read among the 57 with the
   activity 1/3 and one radian per count, or a ratio published before its measurement. The smallest construction
   that would produce (b) and (c) is a consensus theorem of the kind section 7 of the register package's PROOF.md
   gives for the pair on one port (class-3 flux: one consensus state), extended to the analyzer and the kick. No
   such statement is on file (this directory's PROOF.md, section 8).
2. Source. The three minima are on the declared Higgs-Wilson charged source. The frozen connection does not claim
   the transfer to the shared routed field source, and neither does this candidate. If the reviewers require the
   row on the routed source, the remaining obligation is that transfer.
Whatever the rulings: a finite round count is an integer knob (N = 405,309 reproduces the measured muon ratio, step
5e-9); the frozen 2²⁰ was declared with no target and leaves +6.03 ppm.

## Reading under the platform's stricter note

The one credited cell of this row (M5, partial) has three minima measured in-platform and its masses "accommodated
by the article's anchor": electron and muon fix a circulant, "Koide imposed", the tau held out at +0.43σ. Both
quoted phrases stand in that one cell of PRECEDENTS.md; M4 and M6 are negative and no column is validated. By parity,
this candidate has three minima by a declared register, no measured mass in the arithmetic of the ratios, and two
numbers that the repository records as target-informed; it holds out nothing, and its muon ratio misses by 442
standard deviations where the M5 circulant takes the muon as input. A strict maintainer reads it as the present
M10 cell does: "A fit by the author's own parameter count: at least five choices against two observables, the
register designed with the masses known", untested until run in-platform, and at most partial with the M5 tag
after such a run. The menu sharpens the cell's count: three kinds of discrete choice with 45, 14 and 57 options
plus the kick coefficient and the round count, against two observables, with a preregistered null of 82 in 35,910
for a blind draw; the platform's M8 precedent treats such a null as the test of a "within ×3" claim. It does not
change the word.
