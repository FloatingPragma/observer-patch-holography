# Defect classes are centre classes

All arithmetic below is exact: integers, rationals, and the field
`Q(sqrt 5)` for the icosahedral cover and the carrier. Each statement is a
named check in [results.json](results.json); the prefix of the name is the
section number.

## 1. The class rule

The internal space of one quantum is the exterior algebra of
`V = C^3 ⊕ C^2` with hypercharge `diag(-1/3,-1/3,-1/3,1/2,1/2)`. A pattern
occupies `n3` colour levels and `n2` weak levels, so

\[
6Y = 3n_2 - 2n_3,\qquad Q = T_3 + Y,\qquad
k = 6Y \bmod 6 = 3\,(n_2 \bmod 2) - 2\,(n_3 \bmod 3) \bmod 6 .
\]

On all 32 patterns `k mod 3 = n3 mod 3` (colour triality) and
`k mod 2 = n2 mod 2` (weak centre parity). The class multiplicities are
`(4,6,6,4,6,6)` for `k = 0,...,5`.

## 2. The Z6 kernel

The element `(omega I3, -I2, exp(2 pi i Y))`, `omega = exp(2 pi i / 3)`,
acts on a pattern by the phase `n3/3 + n2/2 + Y` turns. This is an integer
on all 32 patterns, so the element acts as the identity. The same holds on
every Standard Model field with triality `t` and weak spin `j`: the phase
`t/3 + j + Y` is an integer, and the field's class gives its triality
modulo 3 and its value of `2j` modulo 2.

## 3. The class table

The 32 patterns, grouped by colour representation `Λ^{n3} C^3`, weak
representation `Λ^{n2} C^2` and hypercharge, are exactly one generation
`Q, u_R, d_R, L, e_R, nu_R` and its conjugates, component by component.

| `k` | Patterns | Fermion multiplets | Left-handed names | Bosons |
| --- | --- | --- | --- | --- |
| 1 | 6 | `Q` | `Q` | |
| 5 | 6 | anti-`Q` | | |
| 4 | 6 | `u_R`, `d_R` | | |
| 2 | 6 | anti-`u_R`, anti-`d_R` | `u^c`, `d^c` | |
| 3 | 4 | `L`, anti-`L` | `L` | `H`, anti-`H` |
| 0 | 4 | `e_R`, `nu_R` and conjugates | `e^c`, `nu^c` | `G`, `W`, `B` |

In the left-handed convention the generation `(Q, u^c, d^c, L, e^c)` with
`nu^c` has classes `(1, 2, 2, 3, 0, 0)`, and the Higgs doublet
`(1,2)_{1/2}` has class 3. The adjoint gauge fields have trivial centre
action and class 0.

The class is the centre class. It does not carry the generation, the weak
isospin `T3` or the rest mass, and a colour singlet and a colour octet share
a class.

## 4. Aharonov-Bohm visibility

A quantum of class `k` carried around a flux of class `f` acquires
`exp(2 pi i f k / 6)`. The phase is nontrivial for

| Flux `f` | Classes that see it |
| --- | --- |
| 2, 4 (colour) | 1, 2, 4, 5: the quark classes |
| 3 (weak) | 1, 3, 5: the left-handed doublets `Q`, `L`, `H` |
| 1, 5 | 1, 2, 3, 4, 5 |

## 5. Reactions

The register has 16 reaction types:
`1+5→0`, `1+1→2`, `5+5→4`, `1+4→5`, `1+2→3`, `4+5→3`, `1+3→4`, `3+5→2`,
`2+2→4`, `2+3→5`, `2→1+1`, `3→1+2`, `1→2+5`, `4→2+2`, `4→1+3`, `5→2+3`.
Together with the pair channels `2+4→0` and `3+3→0` each conserves the class
modulo 6. Their closure under conjugation gives 12 fusions, 12 splits and 3
pair channels.

For each channel `a+b→c` the program lists every triple of fields `A, B, C`
of classes `a, b, c` such that `A ⊗ B ⊗ C*` contains an
`SU(3) × SU(2)` singlet and `Y_A + Y_B = Y_C`. Singlets are counted as
the exact common kernel of the Chevalley generators on the tensor product.
Each realization is typed by its Lorentz content:

| Channels | Vertex with the same classes |
| --- | --- |
| `1+5→0`, `2+4→0` | fermion gauge current: quark and antiquark into a gluon, `W` or `B` |
| `1+3→4`, `3+5→2`, `2+3→5`, `4→1+3`, `5→2+3` | Yukawa: `Q` with the Higgs into `u_R` or `d_R`, and the conjugates |
| `1+2→3`, `4+5→3`, `3→1+2` | Yukawa: `Q` with `d^c` or `u^c` into the Higgs class |
| `3+3→0` | Yukawa `L` with the Higgs into `e_R` or `nu_R`; also `L` anti-`L` and `H` anti-`H` gauge currents |
| `1+1→2`, `5+5→4`, `1+4→5`, `2+2→4`, `2→1+1`, `1→2+5`, `4→2+2` | three quarks in a colour-epsilon singlet only, with no renormalizable vertex |

The four Yukawa operators `Q^† H d_R`, `Q^† H̃ u_R`, `L^† H e_R`,
`L^† H̃ nu_R` are gauge singlets with hypercharge sum zero; their classes
are `1+3→4` (quarks) and `3+3→0` (leptons). Every gauge-invariant triple
among the listed fields conserves the class. Read as a reaction of single
quanta, `1+1→2` changes fermion parity; read as a reaction of flux defects,
which section 7 grades even, it keeps it.

## 6. Closed carriers

The carrier is the icosahedron with 12 ports, 30 seams and 20 faces. With
faces oriented outward, every seam bounds two faces with opposite
orientation, so the face fluxes of any assignment of link values sum to zero:
the classes on a closed carrier add to 0 modulo 6. The contents
`{2,4}`, `{1,5}`, `{3,3}`, `{2,2,4,4}`, `{4,4,4}`, `{2,2,2}`, `{2,3,3,4}`,
`{1,2,3}`, `{3,4,5}` are class-neutral. Among the 21 neutral multisets of
nonzero classes with at most four members, every one holds an even number
of odd classes, so the grading `k mod 2` is trivial on closed carriers.

## 7. Spin and statistics

The binary icosahedral group `2I` is generated in exact arithmetic by
`(1+i+j+k)/2` and `(phi + phi^{-1} i + j)/2`; it has 120 unit quaternions.
It acts on the code by `rho(q) = R(q) ⊕ S(sigma q)`: the rotation `R(q)` of
the carrier on the three vector levels and the Galois-conjugate quaternion
`sigma q` (`sqrt 5 → -sqrt 5`) on the two spinor levels. The program checks:

1. `R(q)` is special orthogonal for every `q`; `sigma` maps `2I` onto a
   group of order 120; `rho` is faithful.
2. The centre of `2I` is `{1, -1}`, so `z = -1` is the unique nontrivial
   central element, and `rho(z) = diag(1,1,1,-1,-1) = exp(i pi 6Y)`.
3. The central involutions of `S(U(3) × U(2))` are `1` and
   `diag(1,1,1,-1,-1)`.
4. On every pattern `z = (-1)^{n2} = (-1)^{6Y}`, which is the holonomy of a
   class-3 flux.
5. No element of `S(U(3) × U(2))` acts as `(-1)^{degree}`: it would be
   `-I5`, of determinant `-1`.
6. Rotations permute the 20 faces, and exactly `±1` fix every face, so `z`
   acts as `+1` on every placement space of defects.
7. `z` is constant on each block `(n3, n2)` of patterns, so its two-sided
   action on link operators inside a block is `+1`.
8. With the port sectors of premise P5, `z` acts as `(-1)^k` at flux `k`,
   and the sectors 2 (vector) and 7 (spinor) reproduce `rho(z)`.

Gradings are solved over GF(2) on 32 pattern variables, 2 register
variables, 5 defect classes and 12 link blocks. The vertex constraints
(64 currents, 32 register toggles, 120 hoppings, 4 Yukawa vertices, 37 CPT
rows, 27 register channels) admit 16 gradings. On physical states
(premise P2 and section 6) they reduce to four, generated by `(-1)^{number of quanta}` and `(-1)^{coloured quanta}`.
The defect gradings kept by every channel are the trivial one and `k mod 2`.

The rotation by `2 pi` assigns `-1` to each quantum and `+1` to the
register, the defects and the links (premises P1 and P2, items 6 and 7).
This assignment is one of the 16 gradings. Under it all 32 patterns are
odd, the register, defects and links are even, a quark with its antiquark
is even, three quanta on a carrier are odd, and every pure flux content of
section 6 is even. The grading by `6Y mod 2` and the grading by `k mod 2`
are gradings of the vertices and differ from the rotation; the degree
parity is no grading.

## 8. Imported run counts

`results.json` copies counts of the level-10 universe run of the public
simulator (repository, commit and file hashes recorded under
`imported_census.source`): carriers by content, pair births and
annihilations per conjugate pair, identities per class. The package
imports these counts and does not recompute them. The program checks that every listed content is
class-neutral and gives `{4,4,4}` against `{2,2,2}` as `-0.75` Poisson
standard errors.

## 9. Central flux costs

The class-`k` element `(exp(2 pi i Y))^k` has eigenphases `-2k/3` (three
colour levels) and `k` (two doublet levels), in units of `pi` modulo 2. A
traceless generator `H` with that exponential has eigenvalues `i pi x_j`,
with `x_j` a lift of each phase and `sum x_j = 0`. The program enumerates
the lifts with shifts up to 3 and takes the minimum of
`tr_5 H^dag H = pi^2 sum x_j^2`. A lift with some `|x_j| >= 5` costs at least
`25 pi^2`, above every minimum found, so the window is complete.

| Class | Shortest `x` | `tr_5 H^dag H / pi^2` | `8 pi^2 / tr_5 H^dag H` |
| --- | --- | --- | --- |
| 1, 5 | `(-2/3, -2/3, -2/3, 1, 1)` and conjugate | 10/3 | 12/5 |
| 2, 4 | `(-4/3, 2/3, 2/3, 0, 0)` and conjugate | 8/3 | 3 |
| 3 | `(0, 0, 0, -1, 1)` | 2 | 4 |

The class-3 generator is `-2 pi i T3`, of the same length as `2 pi i T3`, and
its exponential is the central
involution `diag(1,1,1,-1,-1)` of section 7. For `alpha = g^2 / (4 pi)` the
Gaussian plaquette action of one class-`k` flux quantum, `tr_5 H^dag H / g^2`,
equals `2 pi / (beta_k alpha)` with `beta_k = 8 pi^2 / tr_5 H^dag H`. The
weak-centre value `beta_3 = 4` is the coefficient of the electroweak
transmutation exponent of the fine-structure derivation.

## Premises

The exact checks above use these stated premises; the program does not
derive them.

- P1. Each quantum carries a spin doublet on which the rotation by `2 pi`
  acts as `-1`.
- P2. On a Gauss-invariant state the code, register and link-end factors of
  `z` multiply to `+1`, so the rotation by `2 pi` equals
  `(-1)^{number of quanta}`.
- P3. The local vertices are those listed in section 7.
- P4. Occupations are lifted antisymmetrically. The tick alone admits both
  lifts; a filled reference forces the antisymmetric one, and this package
  does not compute that reference.
- P5. Vector and spinor levels carry port loads 2 and `-3`, seen at flux
  `11 × load mod 20`.
- P6. A carrier is a closed icosahedral screen whose register classes are
  face fluxes of seam links.

Gauge coupling values, masses and generations are outside the scope of
this package.
