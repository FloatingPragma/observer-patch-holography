# The executable construction of A-contrast-memory-spherical-v1: hardware and software

The construction on one sheet and panel by panel. Hardware is the carrier with its registers. Software is the
law that acts on it, with the instruments, the clock and the records of an observer.

A-contrast-memory-spherical-v1 is the executable model of the theory: one hardware, and the laws under test on it. Law 1, the
contrast-memory law, is executed by the package [contrast_memory_spherical_v1](../contrast_memory_spherical_v1/README.md), whose
[specification](../contrast_memory_spherical_v1/SPEC.md) is the binding text for it. Law 2, the scalar channel, is executed by
the package [scalar_channel](../scalar_channel/README.md) on one cell, and its cage by the subpackage
[cage](../scalar_channel/cage/README.md). Both laws are candidates under test, and the search for the correct
native repair law is work in progress. The two laws have different state spaces, and their results are not
combined into one claim.

[![The executable construction on one sheet: the idea, the cell, the registers, the federation, the two laws, a cage of the second law, one run](executable_spec.svg)](executable_spec.svg)

**Full size.** Click the sheet, or open [executable_spec.svg](executable_spec.svg) in a browser tab and zoom. From
a local checkout, [viewer.html](viewer.html) shows the sheet full screen with zoom and pan.

| File | Contents |
| --- | --- |
| [CONSTRUCTION_MANUAL.md](CONSTRUCTION_MANUAL.md) | the construction manual: every part, how to build it, why it is designed this way, every equation with what it says, and each section again in plain words |
| [executable_spec.svg](executable_spec.svg) | the whole sheet, 2400 by 4672 |
| [fig_idea.svg](fig_idea.svg), [fig_cell.svg](fig_cell.svg), [fig_registers.svg](fig_registers.svg), [fig_federation.svg](fig_federation.svg), [fig_law1.svg](fig_law1.svg), [fig_law2.svg](fig_law2.svg), [fig_cage.svg](fig_cage.svg), [fig_run.svg](fig_run.svg) | the eight panels as separate figures |
| [viewer.html](viewer.html) | a full-screen viewer for a local checkout |
| [make_diagram.py](make_diagram.py) | draws every figure from the carrier files |
| [make_paper_figures.py](make_paper_figures.py) | draws the four figures of the executable model in the papers, `paper/figures/model_*.pdf` |
| [inputs/CELL.npz](inputs/CELL.npz) | one cell with its seams and wedges, in the labelling of the manual |

The figures show the design. The packages named in
[section 11 of the manual](CONSTRUCTION_MANUAL.md#11-what-this-repository-executes) hold what is executed, with
their receipts.

## The idea: a Universe that computes its own history

![The idea](fig_idea.svg)

- **The loop.** The theory's closure hypothesis: the Universe is a timeless, self-referential structure. Its
  observers recover the specification of their own world from completed records and build hardware whose
  execution realizes that same history, its builders and their hardware included. With `U` the complete history,
  `D(U)` the account of its architecture recoverable inside it and `R` a realization of that account built
  inside it, the condition is `R(D(U)) ≅ U`. No host, no outside time and no global clock are supplied.
- **A brain that thinks.** The hardware is a network of very simple observers. Each holds a few records and
  compares them with its neighbours. Thinking is the settling of their disagreements, and what survives every
  comparison is a public fact.
- **Why the problem can be solved.** A consistent history is a fixed point: every record agrees with every record
  it overlaps. A search through all candidate histories grows exponentially with the number of records. The
  observers do not search. Each repairs what it shares with one neighbour, and an accepted repair lowers a global
  measure of disagreement by a fixed step, so the number of repairs is bounded by the disagreement at the start.
  For the native repair of the theory this is a theorem: an accepted repair moves one unit across a seam whose
  loads differ by two or more and lowers the sum of squared loads by at least 2. For the two laws described here
  a bound on the work is the target, and it is work in progress.
- **The size.** The tests run on one cell and on 20, 80 and 320 cells. The theory counts about `8 × 10^122` cells
  of equal area on the screen of the Universe we inhabit, a display on the measured side; the identification of
  a carrier cell with a screen cell is a separate theorem of the theory.

## Hardware 1: the cell

![The cell](fig_cell.svg)

- An echosahedron: twelve ports at the corners of an icosahedron, thirty seams along its edges, twenty faces.
  Every port has five seams. The port opposite to port `p` is port `11 - p`.
- The thirty seams, in the labels of the standalone cell file: `0-1 0-2 0-3 0-4 0-6 1-2 1-3 1-5 1-7 2-4 2-5 2-8 3-6
  3-7 3-9 4-6 4-8 4-10 5-7 5-8 5-11 6-9 6-10 7-9 7-11 8-10 8-11 9-10 9-11 10-11`.
- Federation files use a different port numbering. One map that keeps every seam sends standalone ports 0 through
  11 to federation local ports `0,1,5,7,11,9,10,8,4,6,2,3`; under it the antipode of federation port 0 is 3. Use
  the incidence of the carrier being run to address an observer. The map is checked on all 420 cells
  ([test_contrast_replay.py](../test_contrast_replay.py)).
- A wedge is a pair of seams at one port, written `(a, b, c)` for the seams `a-b` and `b-c` with `a < c`. A port has
  ten wedges and a cell has 120.
- The symmetry group has 60 rotations, 120 elements with the inversion.
- Exact numbers of the cell: the port Laplacian `L = 5 - A` has the levels 0, `5 - sqrt5` three times, 6 five
  times, `5 + sqrt5` three times. The Gram matrix of the port directions has the entries 1, `1/sqrt5`, `-1/sqrt5`,
  -1 at seam distances 0 to 3, with `G^2 = 4G` and rank three.

## Hardware 2: the registers

![What sits on a port, a seam and a wedge](fig_registers.svg)

| Place | Law 1 | Law 2 |
| --- | --- | --- |
| port | a load 0, 1 or 2; on a port with five seams one reference register at zero | a fuel mode at level 0, a parity reference, two pointers, a record memory |
| seam | a private memory of eleven states: blank, and ten records `(c, x)` with `c` in {0, 1} and `x` in {-2, ..., 2} | one mode that holds one object or none, at level `+m` |
| wedge | | one mode that holds one object or none, at level `-m` |

On one cell law 2 has 30 seam modes, 120 wedge modes, 12 fuels and 12 references: 174 modes of one algebra, written
on 432 two-state registers (120 transit bits, 120 bits in the codes of the ports, 192 apparatus bits).

## Hardware 3: the federation

![The federation](fig_federation.svg)

| Cells | Ports | Seams | inside cells | between cells | Wedges | Ports with six seams |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 240 | 630 | 600 | 30 | 2,700 | 60 |
| 80 | 960 | 2,520 | 2,400 | 120 | 10,800 | 240 |
| 320 | 3,840 | 10,080 | 9,600 | 480 | 43,200 | 960 |

- The cells sit on the triangles of the geodesic tower: the icosahedron, with every triangle cut into four at each
  refinement. A cell has three neighbours.
- Neighbouring cells share one external seam, from one port of each. Two external seams never share a port, so
  a port has five or six seams.
- The gluing rule is a declared convention of the builder: the seam towards a neighbour leaves through the port
  whose direction, in the local frame of the cell, has the largest inner product with the direction of that
  neighbour; collisions are repaired by an exact assignment. At 320 cells it gives 216 cells whose three external
  ports lie pairwise two seams apart, and 104 with ports at distances 1, 2 and 3.
- Files: [../contrast_memory_spherical_v1/carriers/](../contrast_memory_spherical_v1/carriers/) `FEDERATION_LEVEL0.npz`, `LEVEL1.npz`,
  `LEVEL2.npz`.

## Software 1: law 1, loads with a seam memory

![Law 1](fig_law1.svg)

Candidate `A-contrast-memory-spherical-v1`, the law of the [executable model](../contrast_memory_spherical_v1/README.md).

    H = (1/2) sum_e (1 - SWAP_e) + kappa C
    C = sum_e (n_a - n_b)^2 + sum_v (6 - d_v) n_v^2 = 6 sum_v n_v^2 - 2 sum_e n_a n_b

- **One seam as a program.** On the two loads and the memory of a seam, 99 states: `h_e = U (G0 + kappa G1) U`
  with `U = [[0, W^T], [W, 1 - W W^T]]`, `W |a, b> = (|0, a - b> |y_0> + |1, a - b> |y_1>) / sqrt2`,
  `y_c = (floor(s/2) + c p, floor(s/2) + (1 - c) p)` for `s = a + b` and `p = s mod 2`, `G0 = (1 - R)/2` with
  `R |c, x> = |c, -x>`, and `G1 = x^2`. With the memory blank this is `(1 - SWAP_e)/2 + kappa (n_a - n_b)^2`.
- **Costs.** Two loads on one seam 10 `kappa`, two loads apart 12 `kappa`, two loads on one port 24 `kappa`.
- **Theorem.** For `kappa > 3`, from every start on a seam and at every time, the probability of "no pair on a
  seam" is at most `40 / (2 kappa - 6)^2`, on every federation of complete cells
  ([proof](../contrast_memory_spherical_v1/PROOF.md)).
- **Kept exactly.** The number of loads and the numbers of ports with one and with two loads; for every seam,
  whether its memory is blank and the size `|x|` of a record it holds.
- **As events.** Every seam acts once per round by `exp(-i tau h_e)`, and each port observer counts its own
  rounds.
- **Two readings.** The reader keeps one record per seam. The presence reading keeps one answer for the whole
  cell, "a pair on some seam", through a reversible tree over the thirty readers. Both find a pair with the same
  probability and leave different states behind, and the records of a later reading tell them apart
  ([presence.py](../contrast_memory_spherical_v1/presence.py)).
- **Supplied.** The coupling `kappa`, the emitter with its fuel, the reader, the reference registers, the routing
  of the federation, the tree of the presence reading.

## Software 2: law 2, the scalar channel with a filled reference

![Law 2](fig_law2.svg)

Candidates `A-scalar-CAR-pair-transit-v2`, its form on ordinary registers `A-native-scalar-overlap-v3`, which
[scalar_channel](../scalar_channel/README.md) executes on one cell, the sibling with round seam terms, and
`A-native-caged-exciton-v1`, the same free law with two more operators in the observer's menu, which
[cage](../scalar_channel/cage/README.md) executes.

    D = [[m, B], [B^T, -m]]        D^2 = m^2 + diag(B B^T, B^T B)        K = B B^T
    B[(a,b), (a,b,c)] = +1/sqrt2        B[(b,c), (a,b,c)] = -1/sqrt2
    H_0 = sum_e n_e - sum_w n_w + (1/sqrt2) sum_w [(c*_ab - c*_bc) c_w + h.c.]        (m = 1)

- **Levels of one cell.** `K` has the levels 0, 1.382 three times, 3 five times, 3.618 three times and 5 eighteen
  times. The frequencies are 1, 1.5434, 2, 2.1490, 2.4495 and their negatives, with 91 wedge modes at `-1`.
- **Reference.** All 120 levels below zero are filled. A particle is one object more and lies at least `m` above.
- **The 150 modes.** 23 carry waves over the twelve ports (`K J = J L / 2` with the port Laplacian `L`), 36 are the
  flat branch (eighteen seam profiles at the level 5 with their wedge images, the cages of the next panel), 91 are
  circulations at `-m`.
- **On registers.** `H = H_0 + 1024 sum_b (1 - P_b) + 512 sum_l (P_S(l) - C_l) + 1024 sum_u (1 - |0><0|)_u`, with 259
  independent loop constraints among the 240 triangles of a star and the 20 faces. They select the objects with
  an even total number.
- **On the sphere.** Round seam terms: `D = [[m, G B T], [T B^T G, -m]]`, with a factor per seam from the area of
  its cell and an amplitude on the wedges of every external seam.
- **Supplied.** `m`, the pulses, the ramps, the apparatus bits, the round terms.

## Software 3: a cage in the flat branch of law 2

![A cage](fig_cage.svg)

Candidate `A-native-caged-exciton-v1`, executed by [cage](../scalar_channel/cage/README.md): the free law of law 2,
an object in its flat branch, and two operators that write and read it.

    e  = (1/2) ( |0-2> - |2-4> + |4-6> - |6-0> )        e' = B^T e / sqrt5
    2K = diag(d_a + d_b) - Q^T Q        Q e = 0   =>   K e = 5 e
    D (e, e') = (e, e') [[m, sqrt5], [sqrt5, -m]]        levels +E and -E,  E = sqrt(m^2 + 5) = sqrt6

- **The ring.** Four ports in a ring, here 0, 2, 4, 6, with the signs `+ - + -` on its four seams. `Q` is the
  port-by-seam incidence without signs and `d_a` the number of seams at port `a`. At every port of the ring a plus
  meets a minus, and with five seams at each of the four ports the profile has the level 5 of `K`.
- **What the law keeps.** The pair `(e, e')` on 4 seams and 28 wedges, coupled to nothing else. A cell has 30
  rings of rank 18, the whole flat branch. No profile on fewer than four seams exists.
- **The object.** In the reference the lower mode is filled. `X^+ = a_+^+ a_-` lifts its object into the upper mode:
  `[H_0, X^+] = 2E X^+`, a neutral excitation at `2E = 4.898979`. Its projector `P_X = n_+ (1 - n_-)` commutes with
  `H_0`.
- **The swing.** With the object on `e`, `P_e(t) = 1 - (5/6) sin^2(sqrt6 t)`: a constant and one line at `2E`, from
  rest. The pulse `H - 2E n(e)` for the time `pi / (2 sqrt(2E (E - m)))` writes that state and does the work 3.456.
  The pulse `H - (2E^2/m) n(e)` for `pi m / (2 sqrt5 E)` writes a sharp excitation, which does not swing.
- **Read from counts of one cell.** The line of the seam record at 4.898982 (standard error 0.000002) of the
  observer's own reference; the level difference of the pair, read separately with a fuel dial, at 4.898980
  (0.000008). The excitation at five waits: in the state 0.99995 after the sharp pulse, 0.70593 after the pulse of
  the swing and 0.000004 with no pulse, and the detector reads 0.9999 to 1.0002, 0.7059 to 0.7060 and 0.000004.
- **Under a changed law** with seam hopping round every port (same alphabet, routing and symmetry) the response has
  no single line, and the cage holds 0.79, 0.29, 0.21 and 0.14 of its weight at 1, 4, 16 and 64.
- **On a federation.** Rings on ports with five seams are exact: 92, 312 and 1,272 of them on 20, 80 and 320 cells.
  Of the eighteen profiles of a cell 6 or 7 keep the level.
- **On registers.** `n(e) = Q_e / 8` with `Q_e^2 = 8 Q_e`: 17 words of weight at most 15. The converter
  `c(e)^+ f + f^+ c(e)`: 8 words of weight at most 24. Both commute with the 260 loops and touch the registers of
  the ports 0, 2, 4, 6; one guard over the codes of those four ports stands in front of each.
- **Supplied.** The two operators with their guard, the fuel dial, the work and the phase of the pulse.

## Software 4: observers, instruments, clocks

![One run](fig_run.svg)

1. **Prepare.** Every loop is measured, its sign kept and corrected; the transits are filled; twelve star ramps of
   256 steps with the profile `3 s^2 - 2 s^3` bring the cell to its reference.
2. **Write.** A fuel at level 0 is coupled to a local target with strength 64 for a pulse of area `pi/2` or `pi/4`
   while the law keeps running.
3. **Run.** The law acts alone. The observer counts cycles of its own circulation
   `chi_b = (|a,b,c> + |c,b,f> - |a,b,f>) / sqrt3`, which sits at the level `-m` exactly.
4. **Read.** A second pulse with one of four phase settings; with the coupling off, the fuel is copied into two
   pointers. For a cage the fuel is read twice in one trial: first for the observer's phase, then, used as found,
   for the capture of the profile.
5. **Agree.** Identical complete reports of two observers are admitted. Every other outcome is kept as an outcome.
6. **Record.** Three copies of the admitted record; the order and the count of its own records are the observer's
   clock.

The reading: the probability of "fuel full" with the fuel prepared full, minus the same with the fuel prepared
empty, is the squared amplitude of one added object, whatever the background is. Four settings give one complex
number per count, and its spectrum over the counts shows the frequencies in units of the observer's own reference.

**Observers first.** Ports, seams and cells are the hardware of a computation, and they are no space. What is
compared with nature is what the observers record and agree on. No observer reads the state of the program, its
evolution parameter or a common clock.

## Regenerate

```sh
python model/construction_manual/make_diagram.py          # writes the sheet and the eight figures
python model/construction_manual/make_paper_figures.py    # writes the figures of the papers, with rsvg-convert and gs
```

The script needs NumPy and NetworkX from `model/requirements.txt`. It reads the cell from `inputs/`, the
federations from `../contrast_memory_spherical_v1/carriers/` and the numbers of the cage panel from
`../scalar_channel/cage/results.json`. The figures of the papers use the same geometry; their stability plot is
computed from the law on the federation of 20 cells, and their PDFs are rendered by `rsvg-convert` and normalized
by Ghostscript, as the diagrams of the book are.
