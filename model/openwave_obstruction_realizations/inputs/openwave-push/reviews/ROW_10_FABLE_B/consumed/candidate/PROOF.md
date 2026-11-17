# Row 10, construction `row10-carrier-residual-symmetry-fable-v1`: statements and derivations

Criterion: "3 flavors + oscillation; the PMNS angles". Every number below is in `RESULT.json` and in the log of the
script named beside it. Rules fixed in advance: `PREREGISTRATION.md` (16:17:50, before any script),
`PREREGISTRATION_ADDENDUM.md` (16:40:20, before `s3` and `s4`). Items fixed after a run are marked.

## 1. Setting

I: the rotation group of the icosahedral carrier, 60 matrices over Q(√5) (`s1`: closure at 60; g gᵀ = 1 and
det g = 1 for all; orders 1, 2, 3, 5 in numbers 1, 15, 20, 24; 59 subgroups, by order 1, 2, 3, 4, 5, 6, 10, 12, 60 in
numbers 1, 15, 10, 5, 6, 10, 6, 5, 1). Irreducible dimensions 1, 3, 3, 4, 5 (characters orthonormal, exact); the two
three-dimensional ones are 3 and its conjugate 3′ (√5 → −√5). The theory repository's family triplet is this
space: its port projectors are v vᵀ/|v|² for the 12 port vectors, with a `conjugate` option
(`model/common_source/native_family_record_memory.py`, `_projectors`).

The 10 subgroups of order 3 are the stabilizers of the 20 faces and act regularly on the three corners; the 6 of
order 5 are the vertex stabilizers; the 15 of order 2 the edge stabilizers (`s1`).

## 2. Which subgroups fix a basis

**Statement 1.** On 3 and on 3′ the commutant of a subgroup H has dimension (1/|H|) Σ_h |χ(h)|² = 9, 5, 3, 3, 3, 2, 2,
1, 1 for the orders 1, 2, 3, 4, 5, 6, 10, 12, 60 (`s1`, exact). A subgroup fixes three lines if and only if it is
abelian with three pairwise distinct characters: the 10 of order 3, the 6 of order 5, the 5 Klein groups. Order 2
fixes one line and leaves a plane. A non-abelian subgroup forces equal levels (two for orders 6 and 10, three for
orders 12 and 60): a degenerate pair gives one frequency, the surrogate the gap rejects.

**Statement 2 (no free number).** Let A, B be two of the 21 basis-fixing subgroups, P_ψ = (1/|A|) Σ_a conj(ψ(a)) ρ(a) and
Q_ψ′ the projectors of their characters. Every A-invariant Hermitian operator is Σ μ_ψ P_ψ (the commutant is spanned
by the P_ψ); with three distinct μ its eigen-lines are the P_ψ whatever the μ. The matrix of squared moduli between
the two eigenbases is M_ψψ′ = Tr(P_ψ Q_ψ′) = (1/|A||B|) Σ_{a,b} conj(ψ(a)) conj(ψ′(b)) χ(ab): a function of the pair of
subgroups alone. Executed: projector products and the character sum agree exactly on 882 pairs (441 on each
triplet); random invariant Hermitian operators reproduce the moduli to 1.1 × 10⁻¹⁴ with smallest level gap 0.024;
commutant dimension 3 in every case. Control: a non-invariant term of norm 10⁻³ moves the moduli by 5.2 × 10⁻⁵ to
1.2 × 10⁻² on 11 patterns (`s5`); on P01 the change is of second order (3.5 × 10⁻⁸ in `s1`, 2.6 × 10⁻⁷ and 2.6 × 10⁻⁵ at 10⁻³
and 10⁻² in `s5`), below the preregistered 10⁻⁶: the control does not fail as required on P01.

**Statement 3 (real operators; fixed after the runs, `s5`).** For a cyclic subgroup of order 3 or 5 the projectors of
the conjugate pair satisfy P_ψ̄ = P_ψᵀ. A real symmetric invariant operator, and a symmetric (Majorana-type)
invariant coefficient, therefore has equal levels on the pair. Executed: smallest level gap of a random real
symmetric invariant operator ≤ 8.9 × 10⁻¹⁶ for all 16 cyclic subgroups, ≥ 1.1 for the 5 Klein groups. Consequences:
a cyclic G_ν needs a complex Hermitian neutral operator (a phase, as the register's 2/9 on the charged side); a
symmetric neutral coefficient (the repository's Mν = V* Y V† is symmetric) leaves the Klein groups only. The
repository states the first half for the charged families (claim OPH-W5-STABILISER-POTENTIAL-BOUNDARY: "the C3-
and C5-fixed loci are one-dimensional and every nonzero fixed point has a doubly degenerate spectrum").

**Statement 4 (the two triplets).** The pattern of a pair on 3′ is the image of its pattern on 3 under √5 → −√5.
An outer automorphism of I exchanges 3 and 3′, so the menus coincide. Executed: the same 7, 9, 12 patterns with the
same multiplicities on both; the conjugation exchanges P03 and P05 and fixes the others.

**Statement 5 (CP invariant).** For two cyclic subgroups, complex conjugation exchanges the conjugate pair of rows
and of columns; two transpositions leave J unchanged while conjugation reverses it: J = 0. For two Klein groups the
overlap matrix is real. Exact J² from the moduli (J² = M₁₁M₁₂M₂₁M₂₂ − R², R = (M₁₃M₂₃ − M₁₁M₂₁ − M₁₂M₂₂)/2): zero for every
pattern except P07, where J² = 1/108.

## 3. The patterns (`PATTERNS.md`, `patterns.json`)

Menu H (G_e of order 3; G_ν Klein, order 5 or order 3): 210 ordered pairs, 7 patterns, one orbit of pairs each.
s = √(30 + 6√5)/2.

| Id | Pair | Ordered pairs | ⟨G_e, G_ν⟩ | Rows of squared moduli |
| --- | --- | --- | --- | --- |
| P01 | G_ν = G_e | 10 | order 3 | permutation matrix |
| P02 | third-turn, Klein | 30 | I | (0, (3−√5)/6, (3+√5)/6); twice (1/2, (3+√5)/12, (3−√5)/12) |
| P03 | third-turn, fifth-turn about a corner of the face | 30 | I | (0.01054, 0.18426, 0.80520); (0.18426, 0.63148, 0.18426); (0.80520, 0.18426, 0.01054); entries 1/3 − √5/15, 1/3 + 2√5/15, 1/3 + √5/30 ± (1/12 + √5/60)s |
| P04 | third-turns of two faces sharing an edge (axes 41.81° apart) | 30 | I | ((7−3√5)/18, 2/9, (7+3√5)/18); (2/9, 5/9, 2/9); ((7+3√5)/18, 2/9, (7−3√5)/18) |
| P05 | third-turn, fifth-turn about a far vertex | 30 | I | (0.03519, 0.48240, 0.48240); (0.48240, 0.16500, 0.35259); (0.48240, 0.35259, 0.16500); entries 1/3 − 2√5/15, 1/3 + √5/15, 1/3 − √5/30 ± (1/12 − √5/20)s |
| P06 | third-turns of two faces sharing a corner only (70.53°) | 60 | tetrahedral | 1/9 on the diagonal, 4/9 elsewhere |
| P07 | third-turn, Klein, in one tetrahedral subgroup | 20 | tetrahedral | every entry 1/3; \|J\| = 1/(6√3) |

Menu X adds G_e of order 5: P08 (fifth-turn, Klein; 30 pairs): (0, (5−√5)/10, (5+√5)/10), twice (1/2, (5+√5)/20,
(5−√5)/20); P09 (two fifth-turns; 30): (3/10 − √5/10, 2/5, 3/10 + √5/10), (2/5, 1/5, 2/5), reversed. Menu A adds the
Klein groups as G_e: the transposes of P02 and P08, and P12 (two Klein groups; 20): cyclic (3/8 − √5/8, 1/4, 3/8 + √5/8).
Order 2 as G_ν with G_e of order 3: 4 fixed columns: (1/3, 1/3, 1/3) 60 pairs; (0, 1/2, 1/2) 30; (1/4 − √5/12 twice,
1/2 + √5/6) 30; (1/2 − √5/6, 1/4 + √5/12 twice) 30; one angle and one phase stay undetermined.

Cross-check against the repository (fixed after the runs, `s5`): 6 + 10 + 15 axes, 16 rows of the pair table, smallest
non-zero acute angle 20.905157°, sine 0.356822, arctan(1/φ) between a five-fold and a two-fold axis: the numbers of
the repository's claim OPH-A5-REAL-AXIS-CABIBBO-NOGO (its quark-side menu of the same axes), recomputed by own code.

## 4. The charged side: assignment and correction (`s3`)

**Statement 6.** The register's operator on the triplet is E = 1 + β·mod·(e^{iθ}S + e^{−iθ}S†), S the face third-turn,
β² = 1/2, θ = 2/9 and mod = 1 in the limit, θ = 0.22222222222220767, mod = 0.99999991758370445 at 2²⁰ rounds. E is a
polynomial in S, so its eigen-lines are the three lines of S for every θ, mod and round count; its levels
1 + 2β·mod·cos(θ + 2πk/3) are distinct, so the lines are attached to masses: e on character 1, μ on 2, τ on 0 (the
invariant line). Executed on both triplets, limit and 2²⁰ rounds: commutator with S ≤ 2.2 × 10⁻¹⁶; 1 − overlap of
each eigenvector with its line ≤ 1.1 × 10⁻¹⁵ (both orientations of S); levels against the formula ≤ 8.9 × 10⁻¹⁶.
**The charged side supplies no correction to any pattern.** It supplies the row assignment: menu R (rows by the
register in two orientations, columns free) has 34 distinct angle triples.

The third-turn of the register is the product of the two seam half-turns of the repository's seams (2, 3) and
(2, 10): order 3, fixing the faces (4, 5, 11) and (6, 7, 8) (`s3`; the register's face is (4, 5, 11)).

## 5. Angles and distances (`s2`, `ANGLES.md`; a labelled comparison)

Reference: NuFIT 6.0, table 1, normal ordering, IC24 with SK atmospheric data: 0.308 (+0.012, −0.011), 0.02215
(+0.00056, −0.00058), 0.470 (+0.017, −0.013); in degrees 33.71, 8.56, 43.28.

| Pattern | Best assignment: sin²θ12, sin²θ13, sin²θ23 | Pulls | χ² | Miss in percent | D (degrees) |
| --- | --- | --- | --- | --- | --- |
| P04 | (11−3√5)/19 = 0.22588, (7−3√5)/18 = 0.01621, (8+3√5)/19 = 0.77412 | −7.5, −10.2, +17.9 | 480.6 | −26.7, −26.8, +64.7 | 18.34 |
| P05 | 0.50000, 0.03519, 0.50000 | +16.0, +23.3, +1.8 | 801.4 | +62.3, +58.9, +6.4 | 11.29 |
| P03 | 0.18623, 0.01054, 0.81377 | −11.1, −20.0, +20.2 | 932.1 | −39.5, −52.4, +73.1 | 21.15 |
| P02 | 0.12732, 0, 0.50000 | −16.4, −38.2, +1.8 | 1731.4 | −58.7, −100, +6.4 | 12.80 |
| P01 | 0, 0, 1 | −28.0, −38.2, +31.2 | 3214.4 | −100, −100, +112.8 | 46.72 |
| P06 | 0.50000, 0.11111, 0.50000 | +16.0, +158.9, +1.8 | 25495.3 | +62.3, +401.6, +6.4 | 11.29 |
| P07 | 0.50000, 0.33333, 0.50000 | +16.0, +555.7, +1.8 | 309044.4 | +62.3, +1404.9, +6.4 | 26.71 |

Menu sizes: H 7 patterns, 85 labelled moduli matrices, 79 angle triples (6.30 bits); X 9, 121, 114 (6.83); A 12, 169,
162 (7.34); R 34 triples (5.09). In every menu no triple lies within 3 standard deviations on all three angles, none
inside the three published 3σ ranges, none within 10 % or within 30 % on each angle (on sin² and on the angles).
Nearest by χ²: P04 as above in H, X, A and R (in R: orientation +, δ = 180°, J = 0). Nearest by the largest angular
miss: 11.29° in H (sin²θ12 = 1/2; P05 and P06 tie) and in R (P06), 7.49° (P09) in X and A.

Each angle alone on menu H: θ12 31 values, nearest 0.31878 (+0.9 σ); θ23 25 values, nearest 1/2 (+1.76 σ); θ13 21
values, nearest (7−3√5)/18 = 0.01621 (−10.2 σ, −26.8 %), none within 3 σ. The reactor angle decides.

Null (2 × 10⁶ Haar unitaries, seed 20261114; mean sin²θ13 = 0.33333): a random matrix has χ² ≤ 480.6 with probability
0.0086 and D ≤ 11.29° with probability 0.0297. For the 79 triples of H: N p = 0.68 and 2.35. The nearest member of the
menu is as near as the nearest of 79 random draws is expected to be.

Controls. K3a: the matrix built from the reference angles passes the tolerance and lies 0.091 (largest entry) from
the nearest of the 169 labelled matrices of menu A: a supplied target is rejected by membership. K3b: the same
enumeration in floating point reproduces the icosahedral menu sizes (7, 79, 12, 162; χ² of the nearest equal to
1.0 × 10⁻⁴) and gives for the tetrahedron 3 patterns, 10 triples, nearest χ² 3214.4, for the cube 5 patterns (G_e of
order 3), 38 triples, nearest χ² 1466.0 at (1/3, 0, 1/2), the tri-bimaximal pattern of the row's M5 precedent.

## 6. Selection (`SELECTION.md`, `s3`, `s5`)

Nothing in `model/`, `paper/`, `flagship/`, `claims/` names a subgroup of the carrier's rotation group, or a set of its
lines, as the symmetry or basis of the neutral states (class S3). What the repository has:
- its uncalibrated neutrino package (`model/openwave_neutrino_oscillation`): neutral coefficient Mν = V* Y V†, Y the
  family coefficient of `native_repair_mass` (five port lines on the orbit [1, 5, 11, 10, 7] of a fifth-turn with
  fixed ports 0 and 3, weights 20736, 1728, 144, 12, 1 over 22621), V = I + (−i − 1)P₂ a quarter-turn pulse on the line
  of port 2. Rebuilt here from the port frame: equal to the stored exact Y in all 9 entries; levels equal to the
  stored masses to 2 × 10⁻¹⁷; mixing parameters equal to the stored ones to 2.5 × 10⁻¹². No non-trivial rotation
  commutes with Y (smallest commutator 0.0044) or leaves Mν invariant (0.0035); the neutral frame is 0.25 (F) and
  0.28 (V F) away from the nearest subgroup frame. It corresponds to no enumerated subgroup;
- the ledger's statement (claim OPH-OBSERVER-MATTER-INTERACTION-WITNESS): "Physical source identity, defect
  assignment, mass and neutrino mass-eigenstate assignment are null; the physical particle map is empty.";
- the icosahedral residual-symmetry literature credited as prior work of the particle paper, used for the quark-side
  axis no-go and the stabilizer boundary, not for a neutral basis.
Hence the nearest pattern is a selection of 6.30 bits from menu H (5.09 with the register's rows), and it misses.

The repository's neutral frame, carried to the end:

| Charged frame | Neutral frame | sin²θ12, sin²θ13, sin²θ23 | Pulls | Class |
| --- | --- | --- | --- | --- |
| eigenframe of Y (the package) | V F | 1.89 × 10⁻⁸, 1.764 × 10⁻⁵, 2.138 × 10⁻³ | −28.0, −38.2, −36.0 | derived in the repository; reproduced here |
| register (rows by its levels) | V F (J1) | 0.03380, 0.48021, 0.93037 (orientation +); 0.03874, 0.48360, 0.92992 (−) | −24.9, +818, +27.1 | constructed here: join of two repository architectures on the same triplet |
| register | F (J0) | 0.02493, 0.48794, 0.95291 | −25.7, +832, +28.4 | constructed here |

Over both family coefficients of the repository, both orientations and both triplets: 16 labelled joins, 12 distinct
triples, none within 30 % on each angle; nearest χ² 73,121 (sin² = 0.4714, 0.1730, 0.2091). The heaviest line of Y is
the line of port 1, the same line as port 2 (antipodes), at 79.19° from the register's axis: its squared components
on (τ, e, μ) are (0.0352, 0.4824, 0.4824); with τ on the invariant line this puts 0.48 on \|U_e3\|².

## 7. Neutral masses (`s3`)

Reference ratio r = Δm²21/Δm²31 = 0.029805 (+0.000789, −0.000796); 3σ range 0.026843 to 0.032844.

| Level law | Relative masses | r | Miss | Pull | Class |
| --- | --- | --- | --- | --- | --- |
| levels of Y (the package's masses) | 0.0020642, 0.0636598, 0.9342760 | 0.004638 | −84.4 % (factor 0.156) | −31.6 | derived in the repository |
| the same on the conjugate triplet | 0.0047330, 0.0608109, 0.9344561 | 0.004209 | −85.9 % | −32.2 | constructed here |
| levels of the family record memory's coefficient | 0.00026483, 0.0516737, 0.9480614 | 0.002971 | −90.0 % | −33.7 | repository object read as masses here |
| register spectrum carried over, masses λ² | 1, 206.77, 3477.47 | 0.003535 | −88.1 % | −33.0 | constructed here; premise: the package's "equality of charged/neutral singular values" |
| register levels read linearly | 0.04035, 0.58021, 2.37944 | 0.059189 | +98.6 % | +37.2 | constructed here |

A neutral circulant "with its own balance and phase" drawn from the register's menu of this side's row 11 (6,018
pairs): 260 pairs (4.3 %) give r inside the 3σ range with masses λ², 169 (2.8 %) with masses \|λ\|; a matching pair
would be a selection of 4.5 bits. None is selected. The source states no absolute neutrino scale ("supplies neither
an absolute neutrino GeV scale nor measured mixing angles or mass splittings").

## 8. Oscillation (`s4`)

P(α → β) = \|Σ_k U_βk conj(U_αk) e^{−2iΦ_k}\|², Φ = (0, rΦ₃₁, Φ₃₁), Φ₃₁ = Δm²31 L/4E, on 70,001 points up to Φ₃₁ = 700.
Equivalent form from the squared moduli and J alone: δ_αβ − 4 Σ_{i<j} Re_ij sin²(Φ_j − Φ_i) + 2 Σ_{i<j} Im_ij sin 2(Φ_j − Φ_i),
Re_ij = M_αi M_αj (α = β) or (M_αk M_βk − M_αi M_βi − M_αj M_βj)/2, Im_ij = J ε_αβ ε_ij (sign fixed on the reference matrix).

| Candidate | Masses | Largest appearance | Unitarity defect | Frequencies per initial flavour | Form deviation |
| --- | --- | --- | --- | --- | --- |
| P04, nearest of H and R | levels of Y (r = 0.004638) | 0.6769 | 1.8 × 10⁻¹⁵ | 3, 3, 3 | 3.9 × 10⁻¹⁴ |
| P04 | reference splittings, inputs (r = 0.029805) | 0.6769 | 1.8 × 10⁻¹⁵ | 3, 3, 3 | 4.0 × 10⁻¹⁴ |
| J0, J1 on the register's rows | both | 1.0000, 0.9999 | ≤ 1.6 × 10⁻¹⁵ | 3, 3, 3 | ≤ 9.1 × 10⁻¹⁵ |

On P04 at the first fast maximum (Φ₃₁ = π/2): P(e → e) = 0.936, P(μ → e) = 0.049, P(μ → μ) = 0.274, P(μ → τ) = 0.677. The slow
period is 215.6 fast periods with the repository's levels and 33.6 with the reference splittings. Controls: three
equal masses: largest appearance 7.8 × 10⁻³²; two equal masses: one frequency per flavour (gate needs 3) with
appearance 0.677, the one-frequency surrogate, rejected; G_ν = G_e with three distinct masses: appearance 0.

## 9. One fixed column (two numbers free)

With G_e of order 3: 4 columns, 28 assignments; with the two free numbers fitted, one passes: (1/3, 1/3, 1/3) as the
column of state 2 (G_ν a half-turn lying with the third-turn in one tetrahedral subgroup): sin²θ12 = 1/(3cos²θ13) =
0.3409 (+2.74 σ), cos δ = 0.278 at the reference θ13 and θ23, χ² = 7.5. The others: χ² 224.0, 331.3, 1461.6. The
half-turns the register is built from (the two seam half-turns and their conjugate; they generate a group of order 6
with the third-turn, each inverts it, so each axis is perpendicular to the register's axis) fix the column
(1/2, 1/2, 0) on (e, μ, τ): smallest χ² 8448. Two fitted numbers are a calibration; this lane is reported, not offered.

## 10. What the source lacks (the exact remaining obligation)

1. A neutral structure on the carrier with a residual symmetry. The only neutral coefficient on file has none
   (section 6), and joined to the register it puts 0.48 on \|U_e3\|².
2. If it is to be a full residual subgroup: excluded here for every subgroup and both triplets (79 triples of H, 162
   of A; reactor angle at best −10.2 σ). The charged side cannot repair it (Statement 6).
3. If it is to be a half-turn with the trimaximal column: a half-turn in one tetrahedral subgroup with the register's
   third-turn (the seam half-turns are not: they give (1/2, 1/2, 0)), and a rule from the source for the remaining
   angle and phase, equivalent to deriving θ13 and δ; and a neutral operator with a phase or a symmetric coefficient
   compatible with Statement 3.
4. A neutral level law with r = 0.0298: the one on file gives 0.00464; and an absolute scale.
Smallest construction that would produce clause 3: items 3 and 4 as source statements made without an oscillation
datum, run through `s2` and `s4` unchanged.

## 11. Ancestry of every number

| Number | Class |
| --- | --- |
| the 60 rotations, 59 subgroups, 12 patterns, 8 columns, J, all multiplicities | derived (group theory of the carrier; `s1`) |
| rows: τ on the invariant line, e and μ on the pair; zero charged correction | derived from the register (β² = 1/2, θ = 2/9, 2²⁰ rounds: the register's own numbers, whose ancestry is row 11's question) |
| orientation of the third-turn; choice of triplet | not fixed by the source: both carried |
| Y, its levels 0.0020642, 0.0636598, 0.9342760, r = 0.004638, V, the package's mixing | derived in the repository from declared architecture (five-state memory, event probability 1/12, port 2, phase −i); reproduced here |
| J0, J1 and their angles; menu R; the register-spectrum level laws; the menu scan | constructed here from the repository's objects, each with its premise and menu |
| 0.308, 0.02215, 0.470, 212°, 7.49 × 10⁻⁵, 2.513 × 10⁻³ eV² and their errors | input: reference fit; distance columns and mass set (b) only |
| Haar seed 20261114, random operators | arbitrary, recorded |

Nothing is fitted to a measured angle, phase, splitting or ratio anywhere in this directory.

## 12. Literature (from memory, not from a repository file)

Patterns from a Klein G_ν in the icosahedral group are in the literature of discrete flavour symmetry: P08 is the
golden-ratio pattern (sin²θ12 = (5−√5)/10) of Everett and Stuart (2009, cited by the repository) and of Feruglio and
Paris (2011); P02 and P12 appear in the classification of de Adelhart Toorop, Feruglio and Hagedorn (2011); the
general result that no full residual symmetry of a small group reproduces all three angles is that of Fonseca and
Grimus (2014). A reviewer should check these attributions.
