# Row 10 status: not ready; the row stops on the clause "the PMNS angles"

Criterion (SHA256 d476f5ea…): "3 flavors + oscillation; the PMNS angles". Candidate construction
`row10-carrier-residual-symmetry-fable-v1` (this directory). No OpenWave validation is reported; the official M10
cell is untested. The producer does not set `full_criterion_met`.

## Clause table for the row (which construction carries which clause)

| Clause | This construction: the carrier's rotation group and its residual subgroups | The repository's neutrino packages (by reference: `openwave-push/steps/S14b_neutral_weak/SPECIES_SHEET.md`, row 10, SHA256 d82e7457…) | Carried by |
| --- | --- | --- | --- |
| "3 flavors" | pass, structural: three lines fixed by an abelian residual subgroup on the three-dimensional representation (3 or 3′), which is the repository's family triplet. Gates G1a, G1b, G1c pass (60 rotations, 59 subgroups, 21 basis-fixing; 882 pairs exact by two routes; invariant operators to 1.1 × 10⁻¹⁴). Controls: dimensions 1, 4, 5 are not three; order 2 fixes one line; non-abelian subgroups force equal levels; a real symmetric operator with a cyclic group has two equal levels. One control (non-invariant term) does not fail as required on the permutation pattern | executed on `model/openwave_neutrino_oscillation`: three unequal levels, three distinct frequencies; the count follows from the three families of the declared carrier | either |
| "oscillation" | pass on the executed vacuum evolution: 6 runs (3 candidates, 2 mass sets): largest appearance 0.68 to 1.00, unitarity to 1.8 × 10⁻¹⁵, three frequencies per flavour, closed form to 4 × 10⁻¹⁴. Controls fail as required: equal masses 7.8 × 10⁻³²; two equal masses: one frequency, rejected; G_ν = G_e: 0. Level law from the repository: r = 0.004638 against 0.029805 (−84.4 %, −31.6 σ); no absolute scale in the source | executed with a control that can fail (zero port phase: zero appearance); largest cross-family appearance 0.0069 | the repository's package (a record on a source); this construction adds the controls and the level-law miss |
| "the PMNS angles" | **not met.** Menu H (the register's third-turn against every basis-fixing subgroup): 7 patterns, 79 angle triples; 0 within 3 σ on all three, 0 inside the published 3σ ranges, 0 within 30 % on each angle. Nearest (χ² = 480.6): third-turns of two edge-adjacent faces, sin² = ((11−3√5)/19, (7−3√5)/18, (8+3√5)/19) = (0.2259, 0.0162, 0.7741): −7.5 σ, −10.2 σ, +17.9 σ; −26.7 %, −26.8 %, +64.7 %. Reactor angle alone: no menu value within 3 σ. Nothing in the repository selects G_ν; the nearest is a selection of 6.3 bits (5.1 with the register's rows) and is as near as the nearest of 79 random unitaries (N p = 0.68). Charged-side correction: exactly none. The repository's own neutral frame on the register's rows: sin² = (0.034, 0.480, 0.930) | fails on the uncalibrated member (1.89 × 10⁻⁸, 1.76 × 10⁻⁵, 2.14 × 10⁻³; reproduced here to 2.5 × 10⁻¹²); the six oscillation parameters are inputs on the calibrated member | none |

## Verdict

Not ready for full-row review as a candidate for the count. Two clauses are carried; the third is a negative with
its numbers on every source-derived construction examined:

| Construction | sin²θ12, sin²θ13, sin²θ23 | Pulls | Status of the pattern |
| --- | --- | --- | --- |
| nearest residual pattern (P04) | 0.2259, 0.0162, 0.7741 | −7.5, −10.2, +17.9 | selected from a menu of 79 (6.3 bits); not a derivation |
| residual pattern with a symmetric neutral coefficient (Klein only; P02) | 0.1273, 0, 0.5 | −16.4, −38.2, +1.8 | menu of 19 |
| the repository's neutral frame on its own charged frame | 1.9 × 10⁻⁸, 1.8 × 10⁻⁵, 2.1 × 10⁻³ | −28.0, −38.2, −36.0 | derived in the repository |
| the repository's neutral frame on the register's rows | 0.034, 0.480, 0.930 | −24.9, +818, +27.1 | constructed here (join); menu of 12 |
| splitting ratio from the repository's levels | 0.004638 against 0.029805 | −31.6 | derived in the repository |

## The exact remaining obligation

A source statement, made with no oscillation datum on file, that gives the neutral states a basis on the carrier.
The enumeration fixes what it cannot be and what the smallest admissible form is:
1. not a full residual subgroup of the rotation group, on either triplet, with any assignment (exhaustive: 441 pairs);
2. not the neutral coefficient on file (no residual rotation; near-identity mixing on its own frame; 0.48 on
   \|U_e3\|² on the register's rows);
3. the smallest form that survives: a half-turn residual symmetry lying with the register's third-turn in one
   tetrahedral subgroup (fixes the column (1/3, 1/3, 1/3) of the second state; sin²θ12 = 1/(3cos²θ13) = 0.3409, +2.7 σ),
   together with a rule from the source for the remaining angle and phase. The half-turns the register is built
   from do not qualify: they fix (1/2, 1/2, 0) on (e, μ, τ), χ² ≥ 8448. The source has neither such a half-turn
   structure on the neutral side nor the rule;
4. a neutral level law with ratio 0.0298 (on file: 0.00464, a factor 6.4 low; the register's spectrum carried over
   gives 0.00354; its levels read linearly 0.0592) and an absolute scale (none on file).
Items 3 and 4, once stated, run through `s2_angles.py` and `s4_oscillation.py` unchanged.

## Declared pieces and their class

| Piece | Class |
| --- | --- |
| carrier, its rotation group, the triplet (3 or 3′) | architecture; the triplet and the orientation of the third-turn are not fixed by the source: both carried |
| charged basis = lines of the face third-turn; rows by the register's levels | architecture with row 11's ancestry (balance 1/2, phase 2/9); no number of it enters a pattern |
| residual subgroup for the neutral states | hypothesis under test: not selected by the source: a menu |
| standard vacuum evolution | declared (installed dispersion); not a record on a common source of the group |
| levels of Y as neutral masses | derived in the repository from declared architecture |
| reference angles and splittings | input, distance columns and one labelled mass set only |

Nothing is fitted to a measured angle, phase, splitting or ratio.

## Reading under the platform's stricter note

The one credited cell of this row (M5, partial) rests on a parameter-free tri-bimaximal reading (θ12 35.26° against
33.7°) with θ13 from a breaking and "Oscillation dynamics not run". The corresponding parameter-free readings of this
model's own carrier are worse on θ12 (20.9° for the Klein pattern with θ13 = 0; 28.4° for the nearest pattern with
θ13 = 7.3° and θ23 = 61.6°); the tri-bimaximal pattern belongs to the cube's group and does not occur in the
icosahedral menu (control: cube, nearest χ² 1466.0 at (1/3, 0, 1/2)). The dynamics is run here. A strict maintainer
reads the row as now: untested, with the cell's sentence "the uncalibrated source gives a mixing matrix close to the
identity" confirmed and a second sentence available: the carrier's symmetry alone gives no pattern within 3 σ.

Files moved from `rowfull/ROW_10/`: none (nothing was written there).
