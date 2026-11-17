# Row 10 full-row review: "3 flavors + oscillation; the PMNS angles"

Reviewer: Fable/ROW_10_FABLE_B (second Fable session). Producer: first Fable session. Producer and reviewer are two
sessions of one party; the register's owner decides whether this review counts as one of the two. `audit.py` imports
no candidate code. No other review of this candidate exists (section 8).

## 1. What is reviewed

- Candidate `row10-carrier-residual-symmetry-fable-v1`, directory `openwave-push/rowfull/ROW_10_symmetry/`,
  `MANIFEST.json` SHA256 `4855979585e0fbc523e080dfaeee00dc097730f1988adcaec5d564e125c7ea3d` (equal to the announced
  hash). Pins checked at 18:23 +07: 37 of 37 files of the directory; 33 of 36 consumed files. The three that
  differ were rewritten after the seal: `reverse-engineering-reality/flagship/from_observer_consensus_to_standard_physics.tex`
  (working tree of the repository) and two ledgers of another party (`row_completion/README.md`,
  `INTERNAL_VALIDATION.json`). None enters a number of this review.
- Criterion: `openwave-tuning/row_completion/ROW_10/TASK.md` (SHA256
  `277590eae05029086698d605a753945255efce94ac05b57e3f8f81cb3471253f`) holds "Literal criterion: 3 flavors +
  oscillation; the PMNS angles." and "Derive three neutral flavors and their mixing/oscillation probabilities from one
  fixed source; audit the physical PMNS angle comparison and ancestry. Reject target-fitted mixing matrices and
  degenerate one-frequency surrogates." SHA256 of the criterion string without the period:
  `d476f5ead64371c4e44b4f6cf90912478632756704c328c50e2a51d943ac4c28`, as the candidate states.
- The candidate's own status: not ready; the row stops on the clause "the PMNS angles"; `full_criterion_met` is
  not set by the producer.
- Construction under review: the rotation group of the carrier acting on the family triplet (3 or its conjugate);
  a residual subgroup G_e for the charged states (the generation register's face third-turn) and a hypothetical
  residual subgroup G_ν for the neutral states; the matrix of squared moduli between the two sets of lines.

## 2. Point (1): the enumeration, with the reviewer's code (`audit.py`, `RUN.log`, `RESULT.json`; 3.7 s of CPU)

Own route: the 120 automorphisms of the carrier's seam graph by backtracking (the code of the row-25 addendum), an
own embedding of the ports, the 60 elements of determinant +1 as matrices, the conjugate triplet by √5 → −√5 on
the identified entries, all subgroups by closure of every pair of elements, lines by character projectors, squared
moduli Tr(P Q), patterns up to row and column permutations. The candidate's route is exact arithmetic in a
cyclotomic field; the two share nothing.

| Count | Here | Candidate |
| --- | --- | --- |
| rotations; element orders 1, 2, 3, 5 | 60; 1, 15, 20, 24 | the same |
| subgroups; by order 1, 2, 3, 4, 5, 6, 10, 12, 60 | 59; 1, 15, 10, 5, 6, 10, 6, 5, 1 | the same |
| commutant dimension by order | 9, 5, 3, 3, 3, 2, 2, 1, 1 | the same |
| subgroups that fix three lines | 21 (10 of order 3, 6 of order 5, 5 Klein) | the same |
| ordered pairs on each triplet | 441 | 441 |
| patterns with G_e of order 3 (menu H); adding order 5 (X); adding Klein (A) | 7; 9; 12, on both triplets | the same |
| patterns found here and absent from the candidate | 0 | |
| fixed columns from a half-turn G_ν | 8, of which 4 with G_e of order 3 | the same |
| menu H: labelled matrices, with undefined angles, angle triples | 85, 2, 79 | the same |
| menus X, A: matrices and triples | 121 and 114; 169 and 162 | the same |
| menu R (rows by the register, both orientations) | 34 triples | 34 |
| distinct values of θ12, θ13, θ23 on menu H | 31, 21, 25 | the same |
| distinct entries of all patterns and columns; entries between 0.0203 and 0.0239 | 33; none | 33; none |

Entry by entry: each of the twelve candidate patterns is found on both triplets with every entry equal at the
comparison's resolution (5 × 10⁻¹⁰), the candidate's closed forms evaluated numerically included; the number of
ordered pairs agrees for every pair of subgroup types (for example P04: 30 pairs of third-turns; P06: 60; P07: 20
third-turn and Klein pairs and 20 the other way). The identifications agree: transposes P02 with P10 and P08 with
P11, all others self-transpose; the conjugation √5 → −√5 exchanges P03 and P05 and fixes the rest; |J| is 0 for
every pattern except P07, where it is 0.0962250 = 1/(6√3); the pair generates a tetrahedral subgroup for P06 and
P07, the subgroup itself for P01 and the whole group otherwise. The eight fixed columns agree with their pair counts on both triplets. The 79
angle triples of menu H equal the candidate's as a set (none missing either way).

One count of the first run differed and was this review's error: the per-angle counts came out 36 and 32 because
equal values were compared at a fixed rounding; with clustering within 10⁻⁷ they are 31 and 25. The first log is
kept (`RUN_first_run_per_angle_count_by_fixed_rounding.log`).

**The list is complete for what it enumerates**: full residual subgroups of the rotation group on both triplets,
and half-turns as G_ν. Not enumerated, as the candidate says: the double cover, generalized CP transformations,
groups that are not the carrier's.

## 3. Point (2): the selection

| Quantity (menu H) | Here | Candidate |
| --- | --- | --- |
| triples within 3σ on all three angles; inside the three 3σ ranges; within 30 % on each | 0; 0; 0 (also on X and A) | the same |
| nearest by χ² | P04: sin² = 0.225884, 0.016211, 0.774116; χ² = 480.6045 | the same |
| its pulls; factors against the reference | −7.47, −10.24, +17.89; 0.733, 0.732, 1.647 | −7.5, −10.2, +17.9 |
| reactor angle alone | nearest of 21 values 0.016211, −10.24σ; none within 3σ | the same |
| random unitary with χ² at or below 480.6 | 0.008506 ± 0.000046 (4 × 10⁶ draws; by QR on 4 × 10⁵: 0.008395) | 0.008617 |
| expected number among 79 random draws; probability of at least one | 0.67; 0.49 | 0.68 |

Reference values as the candidate states them (NuFIT 6.0, normal ordering: 0.308, 0.02215, 0.470 with their
errors); not checked against the publication here.

Which statement selects the pair.
- **G_e, derived from declared architecture.** The register's operator is a polynomial in the face third-turn,
  so its lines are that third-turn's lines for every balance and phase. The register's numbers (row 11) enter the
  row order only.
- **G_ν, not selected.** Seven searches of `model/`, `paper/`, `flagship/`, `claims/` (residual symmetry, Klein
  group, tri-bimaximal or trimaximal, golden-ratio mixing, neutrino with subgroup, neutral with stabilizer, PMNS
  with subgroup or residual) find no statement that names a subgroup or a set of lines of the carrier for the
  neutral states: the hits are Klein-Gordon, a code stabilizer of the virtual-exchange package, a spin-statistics
  artifact, a statement about the group itself, a bibliography entry, and uses of "residual symmetry" in three
  code files, three historical completion records and the prior-work list, none about a neutral basis. The neutral package's own proof calls
  its coefficient "a symmetric, explicitly chosen member of the allowed family-tensor space". That coefficient is
  rebuilt here from its stated architecture with the reviewer's port vectors (`audit_neutral_coefficient.py`): no
  rotation other than the identity commutes with Y or leaves M_ν unchanged (smallest change 0.0080 in the
  Frobenius norm on both triplets; the candidate quotes 0.0044 and 0.0035 in a norm this review did not match; all
  are non-zero). The candidate's sentence "nothing in the repository selects G_ν" holds as far as these searches
  and this rebuild reach; the searches are not exhaustive.
- **One source statement bears on it.** The source puts the lepton doublet as a whole in flux class 3 at a face
  (the quotes S1 to S3 of `b/SD11C_kept_pair/PROOF.md`). Read as a statement about the neutral component, it puts
  the neutral states on the charged states' third-turn: G_ν = G_e, pattern P01, no mixing (sin²θ12 = 0, θ13 = 0:
  −28σ, −38σ). This is the reviewer's reading; the source does not state a neutral basis.
- **The Klein statement, derived and confirmed.** A real symmetric invariant operator of a cyclic subgroup has a
  degenerate pair (smallest gap at most 8.9 × 10⁻¹⁶ for all 16 cyclic subgroups, at least 0.285 for the 5 Klein
  groups; a complex Hermitian one of a third-turn: 0.836). With a symmetric neutral coefficient the menu shrinks to
  P02 (0.12732, 0, 0.5; factors 0.413, 0, 1.064; phase undefined) and P07 (0.5, 0.3333, 0.5; |J| = 1/(6√3), phase
  ±90°).
- **The nearest pattern is a selection with the measured angles in view**, and the candidate says so: one of 79
  triples (6.30 bits), one of 34 with the register's rows (5.09 bits), phase 180° (J = 0) against 212°. A list
  from which one member is picked afterwards is a menu; the candidate's statement about that is exact, and its
  comparison with random unitaries is reproduced.

## 4. Point (3): oscillation probabilities

The candidate evaluates the standard vacuum formula on a grid. Here the three-state equation
i dψ/dx = 2 U diag(0, r, 1) U† ψ is integrated by a fourth-order Runge-Kutta step from x = 0 to π/2 (4,000 steps) on
P04 with its best assignment and compared with the closed form (difference 4.3 × 10⁻¹⁵, unitarity 5.6 × 10⁻¹⁵).

| At Φ31 = π/2 | P(e → e) | P(μ → e) | P(μ → μ) | P(μ → τ) | Frequencies | Slow period in fast periods |
| --- | --- | --- | --- | --- | --- | --- |
| levels of the repository, r = 0.0046379 | 0.93617 | 0.04938 | 0.27370 | 0.67691 | 3 | 215.6 |
| reference splittings, r = 0.029805 (input) | 0.93476 | 0.04938 | 0.27512 | 0.67550 | 3 | 33.6 |
| candidate (first row) | 0.936 | 0.049 | 0.274 | 0.677 | 3 | 215.6 |

Controls: three equal levels: appearance 1.2 × 10⁻³¹; two equal levels: one frequency, appearance 0.677 (the
one-frequency surrogate the task rejects); G_ν = G_e: appearance 0.

What this is: an evaluation of the installed vacuum evolution on a member of the menu, with a level law taken
from another package of the repository and joined by the producer. The splittings: the three levels are rebuilt
here from the package's stated architecture (five port projectors on the ring of neighbours of port 0, weights
12⁴, 12³, 12², 12, 1 over 22621): 0.0020642, 0.0636598, 0.9342760 on the triplet and 0.0047330, 0.0608109,
0.9344561 on its conjugate, equal to the quoted values; their ratio 0.0046379 (0.0042094) is 0.156 (0.141) of the
reference ratio. The weights are declared architecture of that package (an event probability 1/12). The second
mass set is entered; the source has no absolute scale. Any three unequal levels with a matrix that is not a permutation oscillate; the clause is carried in that
sense and by no record on one fixed source.

## 5. Point (4): classification of every input

- **Derived** (group theory of the declared carrier, reproduced here): the 60 rotations, 59 subgroups, 21
  basis-fixing subgroups, 12 patterns, 8 columns, their multiplicities and J; the restriction of a symmetric
  neutral coefficient to Klein groups; G_e as the face third-turn of the declared register and the absence of a
  charged correction; the levels of the repository's neutral coefficient and their ratio 0.004638, given its
  declared weights (rebuilt here).
- **Selected** (with the measured angles in view, labelled by the candidate): the nearest pattern P04 and its
  assignment (6.30 bits; 5.09 with the register's rows); the half-turn lane with two fitted numbers (a
  calibration, reported and not offered). Carried both ways and not selected: the triplet (3 or its conjugate) and
  the orientation of the third-turn.
- **Declared**: the carrier and its triplet as the family space; the identification of the register with the
  third-turn's lines; the hypothesis that a subgroup fixes the neutral basis; the standard vacuum evolution; the
  join of the register with the neutral package. **Entered**: the reference angles and splittings (distance
  columns and one mass set). The register's balance and phase have row 11's ancestry and enter the row order only.

## 6. Point (5): verdict

`full_criterion_met`: **false**. This agrees with the producer's own status.

| Clause | Verdict | On which construction |
| --- | --- | --- |
| "3 flavors" | carried in a weaker sense: three lines of an abelian residual subgroup on the declared triplet; the count is the triplet's dimension; no neutral state of the source is exhibited | this construction; the repository's package by reference (not rebuilt here) |
| "oscillation" | carried in a weaker sense: installed vacuum evolution on a menu member with joined levels; three frequencies; controls fail as they should | this construction's formula; the repository's package by reference |
| "the PMNS angles" | **not carried**; the negative is exact on the enumerated class: 0 of 79 triples (0 of 162) within 3σ, reactor angle at best −10.2σ, nearest pattern as near as the nearest of 79 random unitaries | none |

**Is the candidate's negative statement exact.** Yes, on its stated class: the menu is complete and its entries
are confirmed one by one; the selection statement holds as far as the searches reach; the nearest pattern is as
near as a random unitary is expected to be among 79 draws (0.67 expected, probability 0.49).

**Exact remaining obligation.** One object of the source: a neutral coefficient on the triplet that the source
produces (the one on file is "explicitly chosen" and has no residual rotation, section 3). The enumeration fixes its smallest admissible form: symmetric, invariant under one half-turn that
lies with the register's third-turn in one tetrahedral subgroup. That fixes the column (1/3, 1/3, 1/3) of the
second state (60 pairs give this column, confirmed), hence sin²θ12 = 1/(3 cos²θ13) = 0.3409 (+2.74σ); the
coefficient's block on the plane across the half-turn's axis then gives the remaining angle and the phase, and its
levels give the splitting ratio (0.0298 is needed; 0.00464 is on file) and, with a scale from the source, the
splittings. Until such a coefficient is derived, a rule for θ13 and δ is two fitted numbers.

## 7. Reading for the platform's maintainer (a reading, no failed gate)

Untested. The residual-symmetry hypothesis is declared and returns a menu; the calibrated member's six
oscillation parameters are inputs and earn nothing; the uncalibrated member's mixing is close to the identity, as
the cell says. The one credited cell of the row (M5, partial) rests on a parameter-free tri-bimaximal reading
(θ12 35.26° against 33.7°) with the dynamics not run. The parameter-free readings of this carrier are farther on
θ12 (28.4° for P04, 20.9° for P02), so the angle clause gives no case for partial; what this construction adds to
the cell is a second sentence: the carrier's symmetry alone gives no pattern within 3σ.

## 8. Other reviews, register entry, precedent

No review of row 10 by any party exists under `openwave-tuning/row_completion/reviews/` or `openwave-push/reviews/`
at 18:22 +07. Register entry 10 (`INTERNAL_VALIDATION.json`): literal gate "3 flavors + oscillation; the PMNS
angles", targeted gap as quoted in section 1, `internal_validated` false, `work_state`
"external_producer_running". `platform/PRECEDENTS.md` lines 375, 1028, 1079: no validated precedent; M5 partial
with "Oscillation dynamics not yet run"; M10 untested with the cell text used in section 7. Not read within the
budget: the candidate's scripts and library, `JOIN.md` beyond its tables, the repository's neutral package beyond
the head of its proof. The candidate's joins of that coefficient with the register (its angles 0.034, 0.480,
0.930) and the package's own mixing are not repeated here.
