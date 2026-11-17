# Row 10: for the two reviewers

Construction `row10-carrier-residual-symmetry-fable-v1`. Everything is in this directory; `MANIFEST.json` pins every
file and every consumed file. No script imports the theory repository or the OpenWave clone; the repository is read
as text and data (`s3`), the clone is not read at all. The producer's verdict is negative on the angle clause.

## Recompute first (interpreter and thread variables as in the brief; about 23 CPU seconds in total)

1. `python s1_enumerate.py` (5 s): the 60 rotations over Q(√5), the 59 subgroups, the 882 pairs by two exact routes,
   the 12 patterns with multiplicities (`PATTERNS.md`), the 8 fixed columns. Exact arithmetic is `r10lib.py`
   (Q(ζ₁₅) with integer coefficients; 490 lines).
2. `python s2_angles.py` (13 s): angles, pulls, menu sizes (79, 114, 162 triples), the null, the membership control, the
   tetrahedron and the cube (`ANGLES.md`).
3. `python s3_source_join.py` (2 s): the register's rows and menu R; the charged circulant at 2/9 and at 2²⁰ rounds;
   the repository's family coefficient rebuilt and compared exactly; the join; the neutral level laws (`JOIN.md`).
4. `python s4_oscillation.py` (1 s), `python s5_supplement.py` (2 s), `python collect.py`, `python seal.py`.

## Where to attack

1. **The enumeration rule.** Fixed at 16:17:50 before any script. It identifies patterns up to row and column
   permutations, counts ordered pairs, and admits every column order. Check that no rule you find natural adds a
   pattern: other representations have 1, 4 or 5 states; non-abelian subgroups force equal levels; order 2 fixes one
   column. Not enumerated: the double cover and spinorial representations; generalized CP transformations combined
   with a half-turn; products with a group that is not the carrier's.
2. **The one number that decides.** sin²θ13 has 21 values on menu H; the nearest to 0.02215 is (7−3√5)/18 = 0.016211
   (−10.2 σ). An entry between 0.0203 and 0.0239 would put the reactor angle within reach; none exists among the 33
   distinct exact entries of all patterns and columns (`patterns.json`; smallest gap between two entries 0.0057).
3. **Statement 3** (a cyclic residual group with a real symmetric or Majorana-type neutral operator gives two equal
   levels). It removes P03 to P06 for a symmetric neutral coefficient and leaves the Klein patterns. It was fixed
   after the runs; the menus of the tables do not use it (they are the larger ones).
4. **The identification of the register with the third-turn's lines.** The corner space of one face is invariant
   under the face stabilizer only; corner ↦ position vector sends its Fourier modes to the third-turn's eigen-lines
   on 3 (20 of 20 faces, exact), and to those on 3′ with the conjugate vectors. Whether the generation register lives
   on this triplet is a declared identification; the repository's family triplet is this space.
5. **The row rule** (τ on the invariant line): it follows from the register's levels with character 0 heaviest. If
   row 11's numbers are replaced, the lines stay and the assignment may change; the menus H, X, A do not use it.
6. **The join J0, J1** is this producer's construction: the repository keeps the register and the neutrino package
   apart ("This source is distinct from the analyzer-effect charged-lepton mass architecture"). Premise: one triplet,
   one port labelling. Its result is negative either way.
7. **The first run of `s3`** used the wrong family coefficient (that of `native_family_record_memory`); the exact
   comparison with the stored Y failed and exposed it. The run is kept in `logs/first_runs/`; the second run carries
   both coefficients and the checks pass (Y in 9 exact entries; masses to 2 × 10⁻¹⁷; mixing to 2.5 × 10⁻¹²).
8. **Reference values.** Central values and mass ranges are in the repository's frozen file; the angle errors and
   3σ ranges were read from the arXiv HTML table through a fetch tool on 2026-11-14 and agree with the producer's
   memory. Check them against table 1 of arXiv:2410.05380v2.
9. **Controls that did not behave as written**: K1d on P01 (second-order response); the fourth-route check in `s2`
   printed False from a tolerance set too tight (sizes agree; χ² differs by 1.0 × 10⁻⁴). Both are in `RESULT.json`
   under `deviations_from_the_preregistration`, with the vertex-set correction.
10. **Literature attributions** in PROOF.md section 12 are from memory.

## What the producer asks you to decide

Whether the negative on clause 3 is accepted as the row's present state on this construction, and whether the
remaining obligation of STATUS.md (a half-turn residual structure for the neutral states in one tetrahedral subgroup
with the register's third-turn, a source rule for one angle and one phase, a level law with ratio 0.0298) is the
right next object for the row's producers.
