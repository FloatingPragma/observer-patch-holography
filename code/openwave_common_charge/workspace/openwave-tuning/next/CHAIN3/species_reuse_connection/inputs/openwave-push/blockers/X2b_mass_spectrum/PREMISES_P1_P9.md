# The nine premises of the lepton register program, classified

Worker X2b, 2026-11-14. Source of the premises: `reverse-engineering-reality/model/lepton_register_round/PROOF.md`,
section "Premises" ("The exact checks above use these declared constructions; the program does not derive
them."), and `certificate.py`, dictionary `PREMISES`. Executed evidence: `run.py`, gates R1 to R5 in
`logs/run.log`. **No OpenWave validation is reported. Row 11 and every other M10 cell is untested on the platform.**

## The two numbers that decide both ratios

The register program gives m_mu/m_e = 206.770315972729 and m_tau/m_e = 3477.47283710460 from two numbers and one
reading:

- the balance (b/a)^2 = 1/2 of the circulant C = a I + b (z R + z^-1 R^-1);
- the phase theta = 2/9 (z = e^{i theta});
- masses proportional to the squared eigenvalues, electron, muon, tau = characters w^1, w^2, w^0.

The electron is the anchor: both ratios are quoted against the electron character, and every value in MeV uses
the measured m_e = 0.51099895069 MeV (paper, Table `tab:lepton-masses`: "input").

Without the phase the same weights give 1 : 1 : 34 + 24 sqrt 2 = 1 : 1 : 67.94 (gate R5; the repository's
`PHASE_AND_PREMISES.md`, equation 8). The phase alone separates electron from muon. The hierarchy is the square of
a near-cancellation: the electron level 1 + sqrt 2 cos(2/9 + 2 pi/3) = 0.0403 vanishes at theta = pi/12 = 0.2618,
and 2/9 lies 0.0396 below that edge. d ln(m_mu/m_e)/d theta = 56.14.

## Classes

| Class | Meaning |
| --- | --- |
| C | carrier: follows from the twelve ports, thirty seams, twenty faces and sixty rotations |
| N | native law: a record probability of the committed rule (uniform ordered seam, fair coin), executed here on the literal function `model/observer/native.py::native` |
| K | coherent completion: holds on the public coherent seam program (`Observer.repair`, `step`, `inverse_step`), a unitary dilation of N. The native law does not fix it: paper, Proposition "Classical records select no coherent phase" |
| D | declared choice on which neither of the two numbers depends |
| T | declared choice on which the balance or the phase depends, shown by an executed control. The repository lists the balance and the phase themselves as architectural choices made with the masses known (quotes below) |

## The nine premises

| Id | Premise (package wording, shortened) | C | N | K | D | T |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | two seam objects C on (2,3) and E on (2,10), each with load, record, two-orientation address; neutral register F with the third-turn R of the face (4,5,11); coupling T; clock K_theta = exp(i theta J); oriented record V | the two half-turns compose to a third-turn that fixes the faces (4,5,11) and (6,7,8); it acts on the ports in four orbits of three | none | the object event is `Observer.repair` at the two addresses of the seam; the 432-state tick is unitary | the register F (a three-state system that is no port, record or address of the carrier); the coupling; the face (one of twenty, with a sense) | the oriented record V: with V = I the compression is a I + b (z R - z^-1 R^-1), not Hermitian, no mass reading (gate R2.plain_record_is_not_hermitian). theta is a free parameter of the theorem |
| P2 | preparation with activity 1/3, load at the screen tail; analyzer Phi = Fourier zero over three flags | activity 1/3 = 20 of 60 addresses for two loads on non-adjacent ports, 4 of 12 on the hexagon register | given an addressed load: agree 1/2, repair 1/2, orientation uniform, so the three letters weigh 2 : 1 : 1 | none | the load at the tail | balance: the analyzer. Its neutral flag is uniform over both orientations of its active object (support of eight basis states: one letter for "agrees at + or at -"), each charged flag holds one orientation (support of four: two letters for "repairs at +" and "repairs at -"), with equal weights and square roots of probabilities as amplitudes. Controls: neutral flag on one orientation gives (b/a)^2 = 1; weights sqrt2 : 1 : 1 give 1/4; neutral flag idle/idle gives 1/8 (gates R2). Phase: in the 432-state register the exclusive count is 3 alpha/4 - alpha^2/4, equal to 2/9 only at activity 1/3 (gate R3) |
| P3 | hexagon objects with two loads and own twelve-address registers; standing state S = (D_2 + D_10)/sqrt 2 | none | a double load moves with certainty when addressed, for both coins: P(+) = P(-) = 1/n on a carrier of n ports | none | the standing superposition; the twelve-address register | phase: two loads on one port. With one load per cell (the content of the public program) the counts are 1/24 and 1/12 and the exclusive count is 17/144, ratios 8.29 and 245 (gate R3.one_load_per_cell_reads_17_144). The paper lists "the operation that brings two loads onto one port" as open |
| P4 | weak round: kick on the exclusive read of fresh records, inverse tick, return test, N rounds per unit clock time | none | none (the committed law has no inverse tick) | `inverse_step(step(x)) = x` | the round | phase: the exclusive read. The same two counts read otherwise give 17/72 (or), 1/72 (and), 1/12, 1/6, 1/4 (one register); ratios 530/8334, 1.26/76.5, 4.22/155, 25.8/579, 2710/39882 (gate R4). The limit N to infinity: N = 1 gives 0.2053 and m_mu/m_e = 92.29; N = 302.28 gives the measured ratio (gate R5) |
| P5 | declared clock H = sum (I - Re U_x) + sum (pi/4) dGamma(I - S), coefficient one of the kick | none | none | the consensus block of the tick is I - L/60, the native mean law; the edge gate equals the clock at tau in 1 + 4Z | the clock and its coefficients 1 and pi/4 | phase: the coefficient one of the kick. It is the conversion of an event probability into radians. With coefficient kappa the phase is kappa x 2/9; kappa = 1 - 7.9e-7 would cancel the muon residual |
| P6 | oriented toggle T(eta) = exp(i eta J (x) Z_F Z_I) | none | none | the commutator of the toggle with the tick is the kick with h = 2 eta | eta ("sets the lepton scale, is open") | phase: the product Z_F Z_I of two flag parities is the exclusive read; one parity alone reads 1/12 or 1/6 |
| P7 | two cells with independent sixty-address registers, two loads each | n = 12; each port is first port of five addresses and second port of five | P(R+) = 1/12, P(R-) = 1/12, P(written) = 1/6 at every port (gate R3) | none | two carriers that tick together | phase: independence supplies the joint term 1/72. One register gives 1/4 (gate R3.one_register_gives_1_4). `lepton_cycle/PROOF.md`: "Fixed marginals 1/12, 1/6 allow p_FH anywhere in [0, 1/12], giving phases [1/12, 1/4]" |
| P8 | class-k flux as the link phase w^k on one hexagon seam | none | none | consensus states by class: class 0 two separated and one pair, classes 1, 2, 4, 5 none, class 3 the pair only | the link phase; lepton = class 3 | none: neither number depends on it. It is the package's argument for the content of P3: under the coherent law with a class-3 link phase the two-load consensus state is unique, is fixed by the tick, and its blank-record part is the uniform double on the six hexagon ports (weight 1/3) with tick weights 2/3, 1/6, 1/12, 1/12. The committed native law reads no link phase; its stationary two-load state gives the count 445/2178 |
| P9 | return test = projector on the uniform address state; fresh preparation after a "no"; with n partners the kick is (1/n) sum Qx_l | none | none | none | the test; the fresh storage | phase: without the reset the count is 2/9, 7/9, 7/9, 2/9 in the four test sectors; partners read at once give 73/216 and 563/1296; the normalisation 1/n keeps 2/9 (package section 8) |

## What is native in the program

Executed on the literal rule `native(loads, seam, coin)` with the uniform ordered seam and the fair coin
(gates R3 of `run.py`; sampled again with a seed in `verify.py`, V3).

- At an addressed load the record (coin, current) takes four equally likely values, (0,+1), (1,0), (1,-1), (0,0).
  The current letter alone weighs 0 : +1 : -1 = 2 : 1 : 1. The balance (b/a)^2 = 1/2 is this ratio read through
  the analyzer of P2, which uses the current letter and merges the two holds. It holds on every carrier tested
  (tetrahedron, octahedron, cube, icosahedron, dodecahedron) and at every activity.
- Two loads on non-adjacent ports of one carrier read P(+) = 1/12, P(-) = 1/12, P(hold) = 1/6. Two loads on one
  port read P(+) = 1/n, P(-) = 1/n, P(written) = 2/n on a carrier of n ports, 1/12, 1/12, 1/6 at n = 12. Two
  loads on adjacent ports read 1/15, 1/15, 2/15. One load reads 1/24, 1/24, 1/12.
- The exclusive count of two independent cells is 3/n - 4/n^2 for the double load: 2/9 at n = 12, 1/2, 7/18,
  5/16, 7/50 at n = 4, 6, 8, 20. The carrier enters through n only. One load on a six-port carrier gives the
  same pair (1/12, 1/6).
- The native law does not hold the register's content. A double load survives one event with probability 5/6
  and is never formed from two separate loads (the 66 two-port configurations are closed under the rule). The
  stationary two-load state is uniform on those 66 configurations and reads P(+) = 5/66, P(hold) = 5/33: the
  balance stays 1/2, the exclusive count of two independent cells is 445/2178 = 0.20432, and the register
  reading (balance 1/2, count as phase) then gives 88.53 and 1628.1 in place of 206.77 and 3477.47. The value 2/9 needs the loads put back on
  non-adjacent ports, or on one port, before every round: the reset of P9 and the content of P3.

## The phase 2/9: derived from what, chosen when

Derived from: 2/9 = p + q - 2 p q with p = 1/12 and q = 1/6, the exclusive read of "one sign in the defect's
cell" and "a written record in a second cell", for two loads on one port of a twelve-port carrier and two
independent address draws. As an angle it is that probability at one radian per event and clock unit.

The repository's own statements:

- Paper (`deriving_the_particle_zoo_from_observer_consistency.tex`, after the Proposition "Classical records
  select no coherent phase"): "setting a probability equal to an angle in radians is an additional rule."
- `model/charged_lepton_response/SOURCE_PHASE_SELECTION.md`: "A measured event frequency, such as 2/9 in a
  declared register, is a probability. Converting it into radians requires an explicit phase-per-event service
  and clock; a frequency alone is not a unitary phase."
- `model/lepton_cycle/physical_comparison.json`, `exposure`: `architectural_choices` lists "theta=2/9",
  "balance=sqrt(2)", "equal analyzer weights", "independent fresh source pairs", "squared-eigenvalue rate
  interpretation", "electron/muon/tau mode assignment"; `prior_empirical_exposure`: "Particle masses and Koide
  relations were known during architecture development; no blind parameter-selection claim."
- `model/lepton_mass_compiler/results.json`: `programmed_phase_two_ninths_is_passive_source_prediction: false`.
- `docs/POSTDICTION_LEDGER.md`, rows muon and tau: "charged analyzer, retrospective phase 2/9".
- Paper, Table `tab:lepton-register`: "No measured mass enters the evaluation of the ratios; the architecture
  and comparison are not blind."
- Paper, section "The charged mass operator", on the earlier response law that held the phase origin 2/9: its
  charged shape "uses a three-phase Koide-type parametrization" (references Koide 1983 and Brannen 2006); "The
  architecture has target-informed provenance".
- Paper, section "The charged-lepton register": "The balanced circulant was first written down from the measured
  masses; here its balance and its phase are the register's own read counts."
- On the balance: Koide paper (`extra/koide_identity_from_positive_c3_face_circulants.tex`): "The packet
  realizes the target-informed balance premise through an asymmetric block choice, which is where that premise
  enters"; frozen row FZ-10: "the balance was first abstracted from the measured triple".

Dates from the repository history: the certificate that the data "do not select the Stage-5 phase 2/9"
(`code/particles/leptons/derive_charged_z3_phase_holonomy_no_go.py`: "Pure A5/C3 finite geometry cannot emit
delta=2/9 as a character holonomy") was added on 2026-07-12; the response law with the phase origin 2/9
(`derive_charged_mcpr_completion_conditional.py`: "historically selected after exposure to the charged spectrum")
on 2026-07-13; FZ-10 was frozen on 2026-07-28; the register package was first committed on 2026-11-11.

Classification of the phase: the count 2/9 is C + N under the declared reads (P3, P7, P4). Its use as the angle
of the circulant is T: the value was the phase origin of the charged circulant in the repository four months
before the register, and the repository files the comparison as a conditional retrospective postdiction.

## Which premise carries the target

| Number | Native content | Premises that turn it into the number | Control that removes it |
| --- | --- | --- | --- |
| balance (b/a)^2 = 1/2, Koide quotient 2/3, hence m_tau from m_e and m_mu | P(+)/P(written) = 1/2 | P2 (analyzer: one neutral letter for both orientations against one charged letter per orientation, equal weights, square roots as amplitudes) and the oriented record V of P1 | one orientation in the neutral letter: 1; weights sqrt2 : 1 : 1: 1/4; V = I: not Hermitian |
| phase theta = 2/9, hence each of the two ratios | p = 1/12, q = 1/6 | P5 (coefficient one: probability to radians), P4 (exclusive read, N to infinity), P7 (independent cells), P3 (two loads on one port), P9 (reset, partner normalisation), P6 (product of the two parities); in the 432-state form also the activity 1/3 of P2 | any other read, load content, activity, round count or kick coefficient |

No mass ratio of the program is an output with no target as input and no retrospectively chosen phase. The
combination that needs the balance only (the Koide relation) rests on P2 and V. Each single ratio needs the phase
as well. The output that uses neither is the set of three native read probabilities; with the analyzer and no
phase it reads 1 : 1 : 67.94, which the measured triple excludes.
