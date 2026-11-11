# The integer lift of the defect register

Every statement carries one label. EXACT: an integer identity, checked by
the program on the listed states. IMPORTED: a count of the universe runs of
the public simulator, copied with its source and checked for arithmetic
only. ARGUED: a reading. DECLARED: a definition. Each exact statement is a
named check in [results.json](results.json); the prefix of the name is the
section number.

## 1. Definition

DECLARED. The register of the runs is a declared layer on the exact loads
of law C. A repair on the seam `(a, b)` with loads `x_a, x_b` moves the
integer current `j = old x_a - new x_a` and proposes the twist
`delta = j mod 6`. The face that traverses the seam from `a` to `b` changes
its class by `+delta`, the face that traverses it from `b` to `a` by
`-delta`. The twist is applied when
`cost[new] - cost[old] <= kappa j (d - j)` with `d = x_a - x_b`, the cost
table `(0, 5, 8, 10, 8, 5)` and `kappa = 5`.

The lift keeps `j` unreduced. Every seam carries an integer link `n_s`,
zero at the start, and every applied twist adds `j` to it. The value of a
face is

\[
m_f = \sum_{s:\ f \text{ traverses } s \text{ from } a \text{ to } b} n_s
    \;-\; \sum_{s:\ f \text{ traverses } s \text{ from } b \text{ to } a} n_s .
\]

An applied twist adds `j` to one face and subtracts `j` from the other. A
defect face has a nonzero class. A flat face has class 0. The value of a
flat face is a multiple of 6 and is called hidden when it is nonzero. The
minimal value of a class is `0, +1, +2, ±3, -2, -1` for `k = 0, ..., 5`.

## 2. Closed carriers

EXACT. The carrier is the icosahedron with 12 ports, 30 seams and 20
oriented faces. Every seam is traversed once in each direction, by two
different faces of the same carrier. Three ports of every carrier are
glued to neighbours, and every glue seam is traversed once in each
direction by two different screen plaquettes. The program checks this for
the federations of 20, 80 and 320 carriers, with Euler characteristic 2
for the carrier and for the screen.

Consequence (Stokes). For any integer link field the twenty face values of
a carrier sum to zero, and the plaquette values of the screen sum to zero.
An integer word of a carrier therefore has total zero. With the reading of
section 5 the total is `6Y = 0`.

## 3. The exact gates

EXACT on the reference run. The reference implementation executes law C,
the register and the lift in Python integers on the federation of 80
carriers for 200 sweeps, with the preparation and the schedule of the
public simulator. It applies 4786 twists. The program checks:

1. the value modulo 6 equals the class on every face, and the link modulo
   6 equals the register on every seam, after every sweep;
2. the value is the signed sum of the integer links around the face;
3. every carrier total is zero and the screen total is zero, after every
   sweep;
4. on the consistent cuts of depth 1, 2, 4, 16, 64 and 512 the applied
   twists are those of the full run with depth at most the cut, the value
   is the sum over those twists, and items 1 to 3 hold;
5. every applied twist has `|j| <= 3`;
6. a second lift that draws a random representative of `j mod 6` for every
   applied twist has the same classes and the same zero totals, and other
   values.

The reference run equals the public engine at the same size: loads,
register, classes and every applied twist with its sweep, position, seam,
current and depth. The digests are stored in
[imported_census.json](imported_census.json).

IMPORTED. The same gates on the universe runs, with the lift rebuilt from
the twist list of each run:

| Level | Carriers | Applied twists | Stored cuts | Face comparisons | Mismatches | Carrier states with total zero |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 5120 | 374750 | 282 | 29599284 | 0 | 1443840 of 1443840 |
| 5 | 20480 | 1462772 | 215 | 90266030 | 0 | 4403200 of 4403200 |
| 6 | 81920 | 5900072 | 181 | 303964522 | 0 | 14827520 of 14827520 |
| 7 | 327680 | 23149983 | 146 | 23920932 | 0 | 47841280 of 47841280 |

At levels 4 to 6 every stored cut lists both tiers. At level 7 the stored
cuts list the screen, and the end state agrees with the stored final
checkpoint on all 6717442 faces and 10321920 seams by hash. The identities
rebuilt from the twist list equal the stored catalogue in every column.

## 4. Kinematics of the lift

EXACT, by enumeration over the loads `0..5`, both coins, and the cost rule.

- K1. `|j| <= 3`. A contrast `|d| = 1, 2, 3, 4, 5` releases the potential
  `j (d - j) = 0, 1, 2, 4, 6`. In the balanced band the loads take two
  adjacent values and `|j| <= 1`.
- K2. A pair is born only at `|d| = 3` for classes 1 and 5, at `|d| = 4`
  or `5` for classes 2 and 4, and at `|d| = 5` for class 3. A birth on two
  faces of value 0 gives `(j, -j)`, the minimal values.
- K3. A face of value `m` that becomes flat ends at 0 when `|m| <= 2`, at
  0 or `±6` when `|m| = 3`, and at `±6` when `|m| = 4` or `5`. For every
  `|m| >= 4` the face keeps a nonzero multiple of 6. A defect of size 4 or
  more therefore does not move with its value.
- K4. In the band a class-1 defect next to a flat face only moves, by
  `j = -1` on its plus side, and a class-5 defect only by `j = +1`. A
  value `+1` or `-1` crosses faces of value 0 unchanged.
- K5. In the band a face of class 2, 3 or 4 with a flat neighbour does not
  change. Two adjacent defect faces can react in the band, for example
  `3 + 3 → 2 + 4` and `2 + 4 → 1 + 5`.
- K6. The rule of section 1 takes the two classes, `j` and `d`. It never
  takes a value. Two faces of equal class and different value have the
  same future, and no step of a run depends on the lift.

The reference run confirms K2 and K3 on its own twists: 144 births on two
faces of value 0, all minimal; 4493 emptied faces of size at most 2, all
ending at 0; 58 emptied faces of size at least 4, all keeping a multiple
of 6.

## 5. Readings

ARGUED. Three readings of the value were fixed before the census.

**R1, hypercharge.** `m = 6Y` of one multiplet. The values of `6Y` on the
32 occupation patterns of the five-level code are
`3 n_2 - 2 n_3 ∈ {0, ±1, ±2, ±3, ±4, ±6}` (EXACT). Up to size 6 the reading
excludes only `±5`.

| Class | Value | `6Y` of |
| --- | --- | --- |
| 1 | `+1` | `Q = (u_L, d_L)` |
| 5 | `-1` | anti-`Q` |
| 4 | `-2`, `+4` | `d_R`, `u_R` |
| 2 | `+2`, `-4` | anti-`d_R`, anti-`u_R` |
| 3 | `-3`, `+3` | `L` and anti-`H`; anti-`L` and `H` |
| 0 | `0`, `-6`, `+6` | `nu_R` and the gauge fields; `e_R`; anti-`e_R` |

With total zero and all flat faces at 0 the words are (EXACT):
`{4,4,4}`: `(+4,-2,-2)`. `{2,2,2}`: `(+2,+2,-4)`. `{3,3}`: `(+3,-3)`.
`{1,5}`: `(+1,-1)`. `{2,4}`: `(+2,-2)`, `(-4,+4)`.
`{1,2,3}`: `(+1,+2,-3)`, `(+1,-4,+3)`. `{3,4,5}`: `(+3,-2,-1)`,
`(-3,+4,-1)`. `{2,2,4,4}`: three words. `{2,3,3,4}`: four words.

Scope of R1. The value separates `u_R` from `d_R`. It does not separate
`u_L` from `d_L`: both have `6Y = +1`. It does not separate `L` from
anti-`H`. `T3` stays outside the lift. The orientation of the faces is a
convention, so the mirror reading `m = -6Y` with every label conjugated
fits the same counts.

Statistics. A face is a flux. A word names the flux that quanta of those
hypercharges see. Every pure flux content is a boson. The word
`(+4,-2,-2)` is the flux that a `u_R d_R d_R` content sees. It is not a
neutron.

**R2, minimal.** Every face takes the value of least size. The minimal
values of `{4,4,4}` sum to `-6` and those of `{2,2,2}` to `+6` (EXACT), so
section 2 excludes R2 on these carriers. Its surviving form R2' keeps the
defect faces minimal and puts the remainder on a flat face:
`(-2,-2,-2|+6)` and `(+2,+2,+2|-6)`.

**Null models.** N1: every nonzero current is added with no acceptance
rule. N2: the applied twists of the run, each with a random representative
of its class.

## 6. Imported census

IMPORTED. End state of each run, faces of the carriers.

| | Level 5 | Level 6 | Level 7 |
| --- | --- | --- | --- |
| Sweeps | 600 | 500 | 400 |
| Defect faces | 12859 | 51843 | 215988 |
| At a value of R1 | 12828 (99.76 %) | 51674 (99.67 %) | 215285 (99.67 %) |
| At the minimal value | 11994 (93.27 %) | 48154 (92.88 %) | 201117 (93.11 %) |
| Size 4 | 834 | 3520 | 14168 |
| Size 5 | 0 | 8 | 57 |
| Size 7 or more | 31 | 161 | 646 |
| Flat faces | 396741 | 1586557 | 6337612 |
| Flat faces at `±6` | 4031 | 16445 | 66299 |
| Flat faces at `±12` | 7 | 9 | 39 |

Values per class at the end of level 6:

| Class | Faces by value |
| --- | --- |
| 1 | `+1`: 681, `-5`: 5, `+7`: 3 |
| 5 | `-1`: 681, `+5`: 3, `-7`: 4 |
| 4 | `-2`: 20088 (91.66 %), `+4`: 1750 (7.99 %), `-8`: 73, `+10`: 4 |
| 2 | `+2`: 20070 (91.68 %), `-4`: 1770 (8.09 %), `+8`: 48, `-10`: 4 |
| 3 | `+3`: 3360, `-3`: 3274, `+9`: 13, `-9`: 12 |
| 0 | `+6`: 8235, `-6`: 8210, `+12`: 1, `-12`: 8 |

The screen plaquettes of level 6 hold 1891 defect faces: 83.3 percent
minimal, 15.7 percent of size 4, and class 4 at `+4` on 23.7 percent. 1936
of 39071 flat plaquettes hold a nonzero value.

Integer words of the carriers at the end:

| Composition | Word | Level 5 | Level 6 | Level 7 |
| --- | --- | --- | --- | --- |
| `{4,4,4}` | carriers | 66 | 334 | 1325 |
| | `(+4,-2,-2)` | 40 | 191 | 707 |
| | `(-2,-2,-2\|+6)` | 12 | 70 | 321 |
| | other words | 14 | 73 | 297 |
| `{2,2,2}` | carriers | 84 | 335 | 1376 |
| | `(+2,+2,-4)` | 47 | 179 | 779 |
| | `(+2,+2,+2\|-6)` | 21 | 78 | 294 |
| | other words | 16 | 78 | 303 |
| `{2,4}` | carriers | 4988 | 19513 | 78483 |
| | `(+2,-2)` | 4037 | 15650 | 62706 |
| | `(-4,+4)` | 29 | 126 | 548 |
| `{3,3}` | carriers | 768 | 3073 | 12438 |
| | signs opposite | 540 | 2270 | 8982 |
| | signs equal, with `∓6` on a flat face | 223 | 778 | 3362 |

Every listed word has total zero. On the stored cut of level 7 where
`{1,2,3}` is most numerous (depth 1024, 7278 carriers) the words are
`(+1,+2,-3)` on 4685, `(+1,-4,+3)` on 325, `(-5,+2,+3)` on 61, and other
words on 2207, most of them with a hidden value.

Applied twists of level 6 by kind, with the share that leaves the flat
faces unchanged: births 156400, of which 156394 on two faces of value 0;
moves 5352909 (97.93 %); annihilations 87419 (98.42 %); fusions 132826
(94.85 %); splits 48598 (99.71 %); exchanges 121920 (100 %).

Time. On the first stored cut of every run all defect faces are minimal
and no flat face holds a value: 129966 faces at level 6. Size 5 peaks at
depth 192 with 1.05 percent of the defect faces and is 8 faces at the end.
Over all stored cuts of all four levels the share at a value of R1 stays
above 98.5 percent. On the cuts of depth 2 and 3 every `{4,4,4}` carrier
holds `(+4,-2,-2)`: 54 of 54 at level 5, 181 of 181 at level 6 and 867 of
867 at level 7. The other words appear with the hidden values.

Mirror pairs, pooled over the four levels, in Poisson standard errors:
`+4` on class 4 against `-4` on class 2: 9326 against 9418, `-0.67`.
`+3` against `-3` on class 3: 17607 against 17777, `-0.90`. `+6` against
`-6` on flat faces: 43997 against 43834, `+0.55`.

## 7. Tests of the readings

IMPORTED counts against the readings of section 5.

| Question | R1 | R2' | N2 | Level 6 | Level 7 |
| --- | --- | --- | --- | --- | --- |
| Defect faces outside the values of R1 | 0 | 0 | large | 0.33 % | 0.33 % |
| Class-4 faces at `+4` | allowed | 0 | common | 7.99 % | 7.97 % |
| `{4,4,4}` with `(+4,-2,-2)` | all | none | mixed | 57.2 ± 2.7 % | 53.4 ± 1.4 % |
| `{4,4,4}` with `(-2,-2,-2\|+6)` | none | all | mixed | 21.0 ± 2.2 % | 24.2 ± 1.2 % |
| Flat faces with a nonzero value | none | on `{4,4,4}`, `{2,2,2}` only | common | 1.04 % | 1.05 % |

The same quantities under N2 at level 6: 37.0 percent of the defect faces
outside the values of R1, 32.1 percent minimal, 47.7 percent of the flat
faces at a nonzero value. Under N1 the value modulo 6 equals the stored
class on 16.65 percent of the faces, which is chance.

Results.

1. R1 fails as an exact statement. Values of size 8, 9 and 10 remain at
   the end, values of size 5 and 7 occur during the run, and flat faces
   reach `±12`.
2. R2 fails. Section 2 excludes it on `{4,4,4}` and `{2,2,2}`, and 6.8
   percent of the defect faces have size 4.
3. Neither word of `{4,4,4}` is the rule. At levels 5 to 7 the word of R1
   holds on 53 to 61 percent of these carriers, the word of R2' on 18 to
   24 percent.
4. The run is far from N2: the lift is a narrow function of the path.
   Sections 4 and 6 give the reason. Births are minimal by K2, and band
   moves carry `±1` unchanged by K4.

## 8. What the lift carries

EXACT. One integer per face beyond the class: the winding
`(m - minimal value) / 6`. For class 3 the two minimal values `+3` and
`-3` differ by one winding. By K6 the run never reads the winding.

IMPORTED. A value of size 4 or more does not travel with its defect. At
level 6, 15330 class-4 identities are born at `+4` and 129205 at `-2`.
The first group makes 5 moves in total, and each of them changes the
value, as K3 requires: `+4` leaves `+6` behind and continues as `-2`. The
second group makes 1443 moves and keeps `-2` in all of them. A class-1
defect takes the value `+1` plus the hidden value of the face it sits on:
of 2675515 class-1 moves, 2619953 keep the value, and the rest step on or
off a flat face that holds a hidden value.

IMPORTED. The reactions of class 3 have two lifts each. At level 7 the
fusion `1 + 3 → 4` runs as `(+1) + (-3) → (-2)` in 16289 events and as
`(+1) + (+3) → (+4)` in 7979. The split `3 → 1 + 2` runs as
`(+3) → (+1) + (+2)` in 7519 events and as `(-3) → (+1) + (-4)` in 8789.
Under R1 these are the down-type and the up-type lift of one vertex.

ARGUED. The lift is a record of how often the applied current wound around
a face. As a statement about the set of values, R1 holds on 99.7 percent
of the defect faces. The value is no conserved label of a member inside a
class: the law is blind to it, it changes when a defect moves, and one
composition carries several words.

## 9. Loads on the ports of class-3 faces

IMPORTED. The public parity control of the coherent defect program leaves
a defect unmoved exactly when both interior edges of its face see an even
number of loads. The same condition in the run's loads is equal parity of
the three port loads of a face. A double is a face with loads
`(u, u, u + 2)`.

In the band every load is 2 or 3 and no face holds a double, for every
class (EXACT). On the 150 band cuts of level 6 the three loads have equal
parity on 25.75 percent of the resting class-3 faces and on 25.75 percent
of the flat faces, and on 26.74 percent of the faces of classes 1 and 5.
Before the band, a double sits on 1.13 percent of the resting class-3
faces, on 1.36 percent of the flat faces and on 2.24 percent of the
class-3 faces that change before the next cut. The runs show no double
load at a resting class-3 defect. The register of the runs and the parity
control are different rules, so this count does not test the coherent
program.

## Premises and scope

- P1. The register, its cost table and `kappa` are the declared layer of
  the public simulator. The lift adds no rule to it.
- P2. The runs are quenches from the preparation `exact06` with one
  schedule per level. Counts are those of these runs.
- P3. The reading R1 uses the class table of
  `model/defect_particle_classes` and no measured value.
- P4. The run counts are imported. The program checks their arithmetic and
  does not recompute them. The reference implementation covers the gates
  and the kinematics at 80 carriers.

Masses, couplings, generations and `T3` are outside the scope of this
package.
