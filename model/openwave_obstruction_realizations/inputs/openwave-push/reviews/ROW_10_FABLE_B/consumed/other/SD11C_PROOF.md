# SD11C: a law of the source that keeps the pair

Worker SD11C, second Fable session, 2026-11-14. Row 11, criterion "exactly 3 charged-lepton minima; ratios
1 : 207 : 3477". No OpenWave validation is reported; the official M10 cell is untested. Nothing here is a status.

Sources. A: the coherent seam law (public completion, tick "repair after the address reflection", two-load walks
and link phase of `model/lepton_register_round`). N: the committed native rule `native(loads, seam, coin)`. Theory
repository HEAD fd461f1d, read only; every imported file has the same SHA256 before and after each run. Predecessor:
`b/SD11_lepton_selection`, manifest 5655c5b7 (copy in `consumed/sd11`). Gates: `PREREGISTRATION.md` (17:28:07),
`PREREGISTRATION_ADDENDUM.md` (17:34:02), both before `c2`, `c3`, `c4`; `PREREGISTRATION_ADDENDUM_2.md` (17:50:59)
before `c5`. Scripts and logs: `c1_flux_class.py`, `c2_kept_state.py`, `c3_join.py`, `c5_sixty_orbit.py`,
`c4_members.py`, `logs/`. 27 executed checks: 25 hold, and the two outcome gates "kept count against 2/9" report
another value. A first version of `c1` with a slower loop for gate G-1d was stopped before it wrote output and
replaced by the block form of the same check; a first run of `c3` reached the time limit before any output and was
rerun with the same gates in rational and modular arithmetic.

## 0. Result

| Item | Finding | Number |
| --- | --- | --- |
| 1 flux class | class 3 for the lepton doublet is a quoted statement of the source and follows from its definitions; the link sign -1 is the sign of the public central link on a weak level; an odd flux at the register's face and the weak level of the loads are declared | 3 of 6 flux classes give the class-3 table for a weak-level load |
| 2 kept state | the class-3 law has one kept state in the two sectors; its count is another value | 13/162 (hexagon object, R-emit), 5/9 (sixty-address carrier, R-rec); declared 2/9 |
| 3 join | one law with two regimes; read at every event it is the native rule and the phase cancels | 445/2178 for every class |
| 4 restoring process | forward ticks conserve the agreement with the kept state; on the hexagon object the blank pair returns after six ticks with probability one; on the sixty-address carrier it has no exact return | hexagon: counts per tick 2/9, 260/729, 77/729, 20/729, 77/729, 0, mean 298/2187; sixty addresses: long-time mean 0.410951 |

## 1. The flux class

### 1.1 Statements of the source

Quotes are verbatim up to line breaks; display formulas are written in line and, in S5, TeX markup is removed. The
quoted files are copied in `consumed/rer`.

- S1, `paper/tex_fragments/DEFECT_PARTICLE_CLASSES.tex:7`: "A face holds a flux class $k\in\mathbb Z_6=\mathbb
  Z_3\times\mathbb Z_2$, the product of a colour-centre charge and a weak-centre parity."
- S2, same file, Theorem "Centre classes of the declared code" (lines 25 to 37): "a flux of class $f$ acts by the
  character $e^{2\pi i f\cdot 6Y/6}$. The class seen by a pattern is $k\equiv 6Y \pmod 6$, $k\bmod 3=n_3\bmod 3$
  (triality), $k\bmod 2=n_2\bmod 2$ (weak centre)."
- S3, same file, table rows at lines 62 and 63: "$3$ & $1$ & $-1$ & lepton doublet $L$, its conjugate, the Higgs
  doublet & frozen" and "$0$ & $1$ & $+1$ & $e_R,\nu_R$; photon, gluons, $W$, $Z$ & no defect". The same rows are
  in `model/defect_particle_classes/README.md:16-17`.
- S4, lines 69 to 73: "a quantum of class $k$ taken around a flux of class $f$ acquires the Aharonov--Bohm phase
  $e^{2\pi ifk/6}$"; "The electrically charged $e_R$ sits in class zero and sees none of these fluxes."
- S5, `extra/fine_structure_constant_derivation.tex:388-389`: "The central element diag(1,1,1,-1,-1) of the code
  group is the class-three flux of the protected defects and the $2\pi$ rotation of the carrier's cover"; line
  398: "the shortest generators of the five nontrivial central classes give" 12/5, 3, 4, 3, 12/5.
- S6, `model/coherent_seam_repair/program.py:36`: "One explicit program; fibre=5 admits a declared static central
  link."; line 103: `link = -1 if self.central_edge == tuple(sorted((a, b))) and internal >= 3 else 1`;
  `model/coherent_seam_repair/README.md:44-45`: "executes the declared central-link background. The code does not
  evolve that background or derive a many-particle interaction."
- S7, `model/lepton_register_round/certificate.py:69-70` (P8): "a class-k flux through the hexagon is the link
  phase exp(i pi k/3) on the hexagon seam (2, 3), entering as the gauge conjugation of the completion at both
  addresses of that seam".
- S8, `DEFECT_PARTICLE_CLASSES.tex:13`, 211, 311: "A face-flux class supplies no particle decoder, generation,
  chirality or mass."; "a class-three defect alone carries no generation label"; "A class-three flux keeps one
  two-load content, the pair on one port".

### 1.2 Statements with executed checks (`c1_flux_class.py`, `c2_kept_state.py`)

Proposition 1.1 (the class of the lepton doublet). On the 32 patterns of the code C^3_(-1/3) + C^2_(1/2), k = (3 n_2
- 2 n_3) mod 6 has the multiplicities (4, 6, 6, 4, 6, 6); the pattern with one weak level and no colour level has
k = 3; 3 is the one element of order two of Z6; the doublet L with 6Y = -3 has k = 3 (gate G-1a). The pattern with
two weak levels (the charged weak singlet, 6Y = 6) and the right-handed charged lepton (6Y = -6) have k = 0
(control K-1a). The five plaquette counts of S5 are reproduced from the shortest traceless generators of
exp(2 pi i k Y): 12/5, 3, 4, 3, 12/5, with the class-3 generator 2 pi i diag(0, 0, 0, -1/2, 1/2) (gate G-1b). The
class label of the defect register, the centre class of the gauge-structure paper and the "five nontrivial central
classes" of the law for v are one group, the centre Z6 of the code group.

Proposition 1.2 (the link phase). By S2 and S4 a load of class q carried around a flux of class f acquires
exp(2 pi i f q/6). The link phase w^k of P8, w = exp(i pi/3), is this phase for f q = k mod 6; P8 with k equal to
the flux class is the case of a load of class 1. A load on one weak level (q = 3) acquires (-1)^f, the exponents
(0, 3, 0, 3, 0, 3) for f = 0..5; a load on one colour level (q = 4) acquires the exponents (0, 4, 2, 0, 4, 2). Five
of the 36 pairs (f, q) give the exponent 3: (1, 3), (3, 1), (3, 3), (3, 5), (5, 3) (gate G-1c). For a weak-level
load the link phase therefore has two values, +1 and -1, and the value of class 3 is exp(2 pi i 9/6) = -1 = w^3.

Proposition 1.3 (the public program holds this sign). `Observer(fibre=5, central_edge=(2, 3)).repair` on one load
equals, on the internal levels 3 and 4, the fibre-1 repair conjugated by (-1)^(load on port 3) at both addresses of
the seam (2, 3), on 2880 of 2880 basis states, and equals the fibre-1 repair on the levels 0, 1, 2, on 4320 of
4320 (gate G-1d). This conjugation is the one of `classes.repair_with_flux` at k = 3. Controls (K-1d): the
conjugation applied to a colour level agrees on 4296 of 4320 states, and the plain repair on a weak level on 2864
of 2880.

Proposition 1.4 (the table by flux and load). The package's table is reproduced: kept states (separated, pair) =
(2, 1) for k = 0, (0, 1) for k = 3, (0, 0) for k = 1, 2, 4, 5 (gate G-2a). Through k = f q mod 6: for a weak-level
load the flux classes with exactly one kept state are f = 1, 3, 5, and f = 0, 2, 4 give three; for a colour-level
load no flux class gives exactly one (gate G-1e).

Proposition 1.5 (the class of the pair). The class is additive, so two weak-level loads on one port have the class
6 = 0 mod 6 and the hypercharge of the pattern with two weak levels, the charged weak singlet of S3. Their link
phase is the square, +1, around every flux (G-1c: the double exponents are 0 for every f). Around an odd flux the
separated doublet loads have no kept state and the kept state is built on the pair that has the class of the
charged singlet. This is a reading of Propositions 1.2 and 1.4 and adds no check of its own.

### 1.3 Decision

- The class of the lepton: a statement of the source (S3) and a consequence of its definitions (Proposition 1.1).
  It is the class of the doublet. The right-handed charged lepton is in class 0 (S3, S4).
- The link phase: fixed. For a weak-level load it is (-1)^f; the sign -1 of `classes.py` at k = 3 is the sign of
  the public central link (Proposition 1.3). The value exp(2 pi i k/6) of P8 for a general k is the phase of a
  class-1 load.
- Declared, with menus: (a) an odd flux at the face of the register. S6 calls the central link declared and
  static, and the program does not evolve it. Menu: six classes; for a weak-level load two link values, three
  classes each. (b) The level of the two loads. The two-load walks of the package have fibre 1 and carry no
  internal level. Menu: weak level (class 3) or colour level (class 4).
- Statement that would have to hold for a derived class: a face that carries the charged-lepton register holds a
  flux of odd class and the loads of the register are weak-level quanta, as an outcome of a law that evolves the
  link background.

## 2. The kept state (`c2_kept_state.py`, `c3_join.py`)

Proposition 2.1 (multiplicity). On the hexagon object, of the 84 (configuration, record) states 36 are fixed by
every event in every class (their record is one the completion does not read there); the two sectors hold the
other 48. In the two sectors the class-3 law has one kept state. On the sixty-address carrier with the sign -1 on
a string of five seams from the face (4, 5, 11) to the face (6, 7, 8), a single load collects -1 around exactly
these two faces, the separated sector (132 states) has no kept state (the matrix 120 (D - 1) is invertible modulo
2^61 - 1, hence over Q), and the pair sector (72 states) has one (gate G-3b). Two different strings give the same
result, and so does the sign on the one seam (2, 3), the form of the public central link (gate G-3c).

Proposition 2.2 (the state). Hexagon object: 18 components of equal modulus, 6 blank doubles and 12 emission
states; blank weight 1/3, record 2 weight 1/3, record 3 weight 1/3; fixed by the tick (exact in Q(w)).
Sixty-address carrier: 72 components of equal modulus, 12 blank doubles and 60 emission states; blank weight 1/6,
record 2 weight 5/12, record 3 weight 5/12; fixed by the event at each of the 60 addresses. In both the written
weight is v/(v + 1) for the valence v of the port graph (2 and 5): each double has 2v emission states and shares
each with one neighbour.

Proposition 2.3 (the count that the law keeps).

| Instrument pair | First tick of a blank double | Kept state | Kept count |
| --- | --- | --- | --- |
| H: hexagon object, read R-emit | forward 1/12, in-dimer 1/6, count 2/9 | forward 1/36, in-dimer 1/18 | 13/162 |
| S: sixty-address carrier, read R-rec | record 2 1/12, written 1/6, count 2/9 | record 2 5/12, written 5/6 | 5/9 |
| hexagon object, read R-rec (judges nothing) | 1/6, 1/3, count 7/18 | 1/3, 2/3 | 5/9 |

Outcome by the preregistered rule: the law selects another value, in pair H and in pair S (the two outcome gates).
The declared count needs a written weight 1/6. The first-tick written weight is 2/n for n ports, and the kept
written weight is v/(v + 1), at least 1/2 for every valence. The statement that would have to hold for 2/9 on the
kept state is "the blank pair is itself a fixed vector of the tick"; the tick moves it with weight 1/3 (hexagon
object) and 1/6 (sixty addresses).

Proposition 2.4 (what the flux does). The pair sector of class 3 has the letters of the pair sector of class 0,
term by term (gate G-2e): the sign on the flux seam multiplies both emission states of that seam and a double
sees its square. The flux removes the kept states of separated loads: class 0 keeps two on the hexagon object,
with an R-rec count anywhere between 0 and 1, and one on the sixty-address carrier (controls K-2a). The flux
selects the pair against separated loads. It does not move the pair's letters.

Proposition 2.5 (all fixed vectors of the tick). On the 216 states of the pair sector with its 12 addresses the
class-3 tick is unitary within 3e-16 and its fixed space has dimension 8: one vector with the uniform address, the
kept state, and seven with zero address mean (gate G-2d). The package's definition counts the first.

Proposition 2.6 (no reset, from the blank part). The uniform blank double has the squared overlap 1/3 with the
kept state at every tick (13 ticks, exact). It returns exactly after six ticks. Weight of the blank pair at ticks
1 to 6: 2/3, 2/9, 1/9, 2/9, 2/3, 1. R-emit counts: 2/9, 260/729, 77/729, 20/729, 77/729, 0, mean 298/2187. R-rec
counts on the hexagon: 7/18, 91/162, 44/81, 91/162, 7/18, 0, mean 11/27 (gate G-2c).

## 3. The join (`c3_join.py`)

Proposition 3.1 (event for event). One event of the coherent law at the ordered seam (a, b), from a blank record
and with d = load(a) - load(b):
- d = 0: nothing changes, record 0. Native rule: current 0.
- |d| = 1: the load is found at either end with weight 1/2 and the record is 1 with weight 1/2. Native rule: the
  coin moves the load with weight 1/2.
- |d| = 2: one load moves with amplitude 1; the record is 2 for the current +1 and 3 for the current -1, and the
  configuration is the native rule's, for both coins.
The link phase of class k multiplies the amplitude of an image by w^(k (n' - n)), n and n' the loads on the marked
end of the seam before and after. It has modulus one and each image of a basis state has one amplitude, so every
event probability is unchanged. Executed for the six classes on the 78 configurations and 60 addresses (28,080
coherent events against 9,360 native events): kernel rows equal 78 of 78 and letter rows equal 78 of 78 in every
class, doubles event for event (gate G-3a). Control K-3: a completion in which a double never emits agrees on 66 of
78 rows.

Proposition 3.2. The native rule with a class-k link phase is the native rule: the uniform distribution on the 66
two-port configurations is stationary, a double survives an event with weight 5/6, the letters are 5/66 and 5/33
and the count is 445/2178 for every k (gate G-3a, control K-2b).

Decision. One law: the coherent tick, with the link phase inside the completion. It has two regimes that differ by
one piece of the reading. Regime N: record and address are read and renewed at every event; this is the native
rule, the phase cancels, the double decays and the count is 445/2178. Regime A: no read between ticks; the phase
acts through the interference of events of different ticks, the class-3 kept state exists and counts 13/162 or
5/9. The join that SD11 uses is an identity for one event from a blank record. It carries no kept state from one
regime to the other. The piece that stays declared is the rate of reads: every event (N), never (A), or once per
round after a restoration (the package).

## 4. The restoring process (`c2_kept_state.py`)

Proposition 4.1 (agreement is conserved). The tick W is unitary and W v = v for the kept state v, so
<v|W^t psi> = <v|psi> for every state and every t. Executed: 1/3 for the blank pair at 13 ticks, 1/18 for each of
the six one-port doubles at 4 ticks (gate G-4a). Forward ticks do not raise the agreement with the kept state; a
process that raises it has a step outside the unitary law.

Three processes, executed.

| Process | Pieces outside the unitary tick | A round ends in the pair | Count over the rounds |
| --- | --- | --- | --- |
| P-fwd: six forward ticks from the blank pair (hexagon object, class 3) | none | blank pair with probability 1 after six ticks; 2/3, 2/9, 1/9, 2/9, 2/3 in between | read at every tick: mean 298/2187; read at the first tick of each six: 2/9 |
| P-test: tick, read, test of agreement with the kept state, preparation of the kept state on "no" | the test and the map on "no" | kept state with probability 1; "yes" has probability 1 without a kick; the blank pair has weight 1/3 (hexagon), 1/6 (sixty addresses) | 13/162 (pair H), 5/9 (pair S) |
| P-pkg: tick, inverse tick, return test on the blank pair, reset (the declared round) | inverse tick, test, reset | blank pair after every round | 2/9 (control K-4) |

P-test with the declared kick exp(-i h Qx), h = 1/N, on two cells: the "yes" amplitude is 1 - P + P exp(-i/N) with
P the kept count. Phase per unit clock time at N = 1, 10, 1000: 0.069997, 0.080144, 0.080247 for pair H (limit
13/162 = 0.080247) and 0.560626, 0.555601, 0.555556 for the R-rec count 5/9 (gate G-4b). The uniqueness of the
class-3 kept state gives the "no" branch one target; in class 0 the target is a space of dimension three. After
one tick from the blank pair this test says "yes" with probability 1/3.

Declared in P-test: the test as an operation (`program.py:160-164` files the return detector under "the optional
apparatus, not idle repair") and the map on "no" (menu of three: continue the ticks, prepare the kept state,
prepare the blank pair with the uniform address as in P9).

The revival of P-fwd is a process of the law that returns the pair with nothing added: the uniform blank double
lies in the block of the discriminant with the eigenvalues 1 and 1/2, so the tick has the phases 0 and +-pi/3 on
it. With one read per six ticks, at the first tick after the return, every read is a first-tick read and counts
2/9 with no inverse tick, test or reset. The rule "one read per six ticks, in phase with the return" is then the
declared piece; the six phases of that rule give 2/9, 260/729, 77/729, 20/729, 77/729, 0.

Proposition 4.2 (the sixty-address carrier, `c5_sixty_orbit.py`). On the pair sector with its 60 addresses (4320
states) the event is a signed permutation. The uniform blank double is the sum of the kept state, with weight 1/6,
and of an eigenvector of the discriminant with the eigenvalue 4/5, with weight 5/6 (exact, gate G-5a). The orbit
is psi_t = sqrt(1/6) n_1 + sqrt(5/6) (sin(t theta) b - sin((t - 1) theta) a)/sin(theta) with cos(theta) = 4/5; it
agrees with the direct tick over 2000 ticks within 3e-14 (gate G-5b). arccos(4/5) is no rational multiple of pi
(Niven's theorem), so the blank pair has no exact return; the closest approach within 2000 ticks is the distance
0.0009 at tick 869. R-rec counts at ticks 1 to 6: 0.222222, 0.491980, 0.551675, 0.540359, 0.525070, 0.545790; weight
of the blank pair: 0.8333, 0.5156, 0.3540, 0.3988, 0.4435, 0.3793. Long-time means over the equidistributed angle:
record 2 weight 0.2160494 (35/162), written weight 0.4320988 (35/81), blank weight 0.5679012 (46/81), count
0.410951074531; the Cesaro mean of the direct ticks is 0.411237 over 2000 ticks (gate G-5c). On this carrier the
first-tick count 2/9 is read once and no rate of reads returns to it exactly.

## 5. Members (`c4_members.py`)

Balance 1/2, phase equal to the count, exclusive read, unit kick and many rounds as declared in SD11. The last two
columns are the measured ratios divided by the member's.

| Member | Count | m_mu/m_e | m_tau/m_e | Levels positive | Measured / member |
| --- | --- | --- | --- | --- | --- |
| kept state, pair H | 13/162 | 3.98703 | 149.629 | yes | 51.86, 23.24 |
| kept state, pair S | 5/9 | 17.9475 | 79.6279 | no, one level is -0.2467 | 11.52, 43.67 |
| no reset from the blank pair, mean of six ticks, pair H | 298/2187 | 12.2518 | 325.535 | yes | 16.88, 10.68 |
| no reset from the blank pair, long-time mean, pair S | 0.410951 | 37.4142 | 278.956 | no | 5.526, 12.47 |
| native rule with any class phase (regime N; control) | 445/2178 | 88.5270 | 1628.08 | yes | 2.336, 2.136 |
| declared register, or one read per six ticks on the hexagon object | 2/9 | 206.770 | 3477.47 | yes | 0.99999, 0.99997 |

For the phases 5/9 and 0.410951 the levels leave the positive chamber (its edge is pi/12 = 0.2618) and the Koide
quotient 2/3 of SD11's Proposition A does not apply. Every member other than the declared control row misses both
measured ratios by a factor above 2.

## 6. What remains declared

1. An odd flux at the face of the register (six classes; three give the sign -1 for a weak-level load).
2. The weak level of the two loads (two alternatives).
3. The rate of reads: at every event, never, once per six ticks in phase with the return, or once per round after a
   restoration. The declared count 2/9 holds for the last two.
4. In a round with a test: the test and the map on "no".
5. Unchanged from SD11: the exclusive read, the unit kick, the limit of many rounds, the reading rule of the
   balance, the source operation that puts two loads on one port.

## 7. Limits

- The reads are the two of the preregistration. A read outside them is not counted.
- R-emit uses the address register; on the kept state the address is uniform, and the read classifies a written
  component only at the two addresses of its seam.
- The count is read on the kept state with the uniform address. The seven other fixed vectors of the tick are not
  read.
- On the sixty-address carrier the classes 1, 2, 4, 5 are run for the one-event kernel only; the class-0 count of
  separated kept states is a floating-point count (1) with a modular upper bound (1).
- The long-time mean on the sixty-address carrier is a floating-point number. 3595/8748 is its nearest fraction
  with a denominator below 100000 and no proven closed form.
- The six-tick return is shown for the uniform blank double of the hexagon object. Other starts of that object
  overlap blocks with other eigenvalues and are not run.
- No independent review of this directory exists.

## 8. Deviations from the preregistration

- G-2c names a spectral mean with the tolerance 1e-10. The orbit returns exactly after six ticks, so the mean is the
  exact mean over the period; the return is also confirmed on the 216-state matrix within 2e-15.
- G-3b names exact rationals for the kept states of the sixty-address carrier. The executed form is: the nullity of
  120 (D - 1) modulo 2^61 - 1 as an exact upper bound (0 for the separated sector, 1 for the pair sector), and a
  candidate vector verified exactly at each of the 60 addresses. The 108 states outside the two sectors are not
  counted there.
- G-1d is executed in block form: 24 calls of the public repair per program return the images of all basis states,
  because the repair acts inside each block (address, internal level).
- `c5` and its gate G-5 are added after the first results, by `PREREGISTRATION_ADDENDUM_2.md`, to supply the count
  over the ticks on the sixty-address carrier.
