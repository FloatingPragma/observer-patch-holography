# A-contrast-memory-spherical-v1 executable specification

The A-contrast-memory-spherical-v1 model is the executable construction of the contrast-memory law.
The law is a candidate native repair law. This document fixes the carrier, the
state, the local law with its memory, the law on loads, the sectors, the emitter, the reader and its records, the
events and schedules, the resources and the supplied parts. The code in this folder implements every clause;
[PROOF.md](PROOF.md) holds the binding theorem and its proof, and [RESULTS.json](RESULTS.json) the computed values.

## 1. Carrier

**Cell.** One cell (the echosahedron) has 12 ports, 30 seams and 20 faces with the incidence of the
icosahedron. Its seams, in the labelling every cell uses, are `CELL_SEAMS` in [carrier.py](../../executable_model/carrier.py). The table
is isomorphic to the icosahedron built from the coordinates $(0, \pm 1, \pm\varphi)$ and their cyclic shifts, every
port has five seams, every seam lies on two faces, the automorphism group has 120 elements, and the adjacency
spectrum is $5, \sqrt5^{(3)}, -1^{(5)}, -\sqrt5^{(3)}$.

**Federation.** A federation is a set of complete cells joined by external seams. Port $12k + j$ is port $j$ of
cell $k$, and each cell keeps its 30 seams. The external seams form a partial matching on ports: every cell has
three external ports, so a port has five or six seams. An internal seam has exactly two common neighbours of its
ends; an external seam has none; a nonadjacent pair of ports has at most two common neighbours. The model uses
three federations, stored in [carriers/](../../executable_model/carriers/):

| Level | Cells | Ports | Seams | External seams | Ports with 5 / 6 seams |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 20 | 240 | 630 | 30 | 180 / 60 |
| 1 | 80 | 960 | 2520 | 120 | 720 / 240 |
| 2 | 320 | 3840 | 10080 | 480 | 2880 / 960 |

The routing of the external seams (which port of which cell meets which) is a supplied convention.

**The carrier is hardware.** Ports, seams, cells and federations are the hardware of the computation and no
physical space, and no clause of this specification asks that they form a three-dimensional lattice. Physical
space, time, fields and particles are effective descriptions that the observers inside the computation use for
the records they agree on (section 7).

## 2. State

Every port $v$ holds a load $n_v \in \{0, 1, 2\}$. Every seam $e$ holds a private memory of eleven states: one blank
$|\varnothing\rangle$ and ten labels $|c, d\rangle$ with coin $c \in \{0, 1\}$ and contrast $d \in \{-2, \dots, 2\}$. The
memory is coherent working memory shared by the two port observers of the seam. It is no public record, and no
observer reads it. A port with five seams carries one private reference register (a load register held at 0) with
its own eleven-state memory. The local space of one seam is the product of the nine load states $|n_a, n_b\rangle$
and the eleven memory states, 99 states in all.

## 3. Local law with memory

**Native event.** On loads $(a, b)$ put $s = a + b$ and $p = s \bmod 2$. The output for coin $c$ is

$$y_c = \big(\lfloor s/2 \rfloor + c\,p,\ \lfloor s/2 \rfloor + (1 - c)\,p\big).$$

The event conserves the load of the seam.

**Contrast write.** The isometry $W$ from the nine load states into labels times loads is

$$W\,|a, b\rangle = \tfrac{1}{\sqrt2} \sum_{c = 0}^{1} |c,\ a - b\rangle\,|y_c\rangle .$$

It is defined on all nine inputs, doubled loads included, and $W^* W = I$.

**Echo.** On blank plus record sectors,

$$U = \begin{pmatrix} 0 & W^* \\ W & I - W W^* \end{pmatrix}, \qquad U = U^*, \quad U^2 = I,$$

a Hermitian unitary on all 99 states, dirty memory included. Write, phase and reverse write form the echo.

**Memory generator.** On the labels $R\,|c, d\rangle = |c, -d\rangle$ and $D\,|c, d\rangle = d\,|c, d\rangle$. With coupling
$\kappa$,

$$G_\kappa = 0_{\varnothing} \oplus \Big[\tfrac12 (I - R) + \kappa D^2\Big] \otimes I, \qquad h_\kappa = U\, G_\kappa\, U^* .$$

**Exact reduction.** With $T$ the swap of the two loads, $C = (n_a - n_b)^2$ and $J$ the embedding of the blank memory,

$$(R \otimes I)\,W = W T, \qquad (D^2 \otimes I)\,W = W C, \qquad h_\kappa J = J\Big[\tfrac12 (I - T) + \kappa\,(n_a - n_b)^2\Big].$$

These identities hold as polynomial identities in $\kappa$, hence for every real $\kappa$; [law.py](law.py) checks
them with exact SymPy arithmetic on all 99 states. Every local term preserves the blank code, so the sum of the
local terms over overlapping seams preserves the product of blank memories, and on that code the law acts on the
loads alone.

**What the local term keeps.** The local space of a seam splits into three parts that the local term maps into
themselves: the blank code, on which it is the seam term; the written code $\operatorname{ran} W$, on which it
vanishes; and the other records, on which it is the generator. On the records the squared contrast $d^2$ commutes
with the projector on the written code and with the generator, and so with the local term. Every seam therefore
carries a label that no event changes: its memory is blank, or it holds a record of size $|d| = 0$, $1$ or $2$.
The write itself changes the memory and does not commute with $d^2$.

**Degree term.** On the pair formed by a port with five seams and its reference register the generator is the
$D^2$ part alone. Its reduction is $(n_v - n_{\rm ref})^2$, which keeps $n_{\rm ref} = 0$ and equals $n_v^2$ there. With it
every port counts six seams. The exchange part on this pair would move load into the reference and is no part of
the law.

**Binary loads.** On loads 0 and 1, with $n = (I - Z)/2$, the seam term is the ferromagnetic XXZ interaction

$$h_e = \tfrac14 \big[(1 + 2\kappa)(I - Z_a Z_b) - X_a X_b - Y_a Y_b\big].$$

The binding it produces is of the droplet kind known from that literature (Nachtergaele and Starr,
arXiv:math-ph/0009002; Fischbacher and Stolz, arXiv:1712.10276). The construction claims no new interaction; its
content is the memory realization, the constants on this carrier and the tests on records.

## 4. The law on loads

$$H_\kappa = \tfrac12 \sum_{e} (I - \mathrm{SWAP}_e) + \kappa\, C, \qquad
C = \sum_e (n_a - n_b)^2 + \sum_v (6 - d_v)\, n_v^2 = 6 \sum_v n_v^2 - 2 \sum_e n_a n_b .$$

- $\kappa$ is one dimensionless number, the agreement energy in units of the exchange. Tested values are
  0, 0.5, 1, 2, 4, 8, 16, 64 and 128; the headline values are 8 and 128.
- A seam whose two ports both hold a load (an agreed seam) lowers $C$ by 2. Two loads cost 10 as a pair on a seam,
  12 as separate loads and 24 as one doubled load.
- The law conserves the number of loads and the multiset of loads.
- On one load $H = L/2 + 6\kappa$ with $L$ the graph Laplacian. On one cell the levels of $L/2$ are 0 once,
  $(5 - \sqrt5)/2$ three times, 3 five times and $(5 + \sqrt5)/2$ three times. The triplet $(5 - \sqrt5)/2$ is the
  rank-three response from which the observers' effective three-dimensional description is built in the theory.

## 5. Sectors

The law splits into sectors of fixed load number and fixed multiset of loads. On $N$ loads 0 and 1,
$C = 6N - 2\,(\text{agreed seams})$, and the code writes the sector Hamiltonians as $H_\kappa - 6N\kappa$, so that $N$
free loads at rest lie at 0 ([sectors.py](sectors.py)):

| Sector | Dimension | Used for |
| --- | --- | --- |
| two loads on distinct ports | $V(V-1)/2$ | the pair, the theorem, the events |
| one doubled load | $V$ | reducing, $L/2 + 12\kappa$ in these units |
| all two-load states | $V(V+1)/2$ | the reader, the theorem's complement |
| three loads on distinct ports | $\binom{V}{3}$ | the face object |
| $N$ loads with doubled loads (one cell) | 12, 78, 352, 1221 for $N = 1..4$ | bindings and fragmentation |

**Changed laws** (controls with the same alphabet, routing and symmetry): the exchange alone ($\kappa = 0$) and the
agreement energy reversed ($-\kappa$).

## 6. Emitter

One fuel state $|F_e\rangle$ at a seam $e = (a, b)$ holds the two units of load of a pair. Its energy is the diagonal of
the law at the pair, $E_f = 10\kappa + (d_a + d_b - 2)/2$. During the pulse the law acts throughout, with

$$H_{\rm create} = H_\kappa + E_f\,|F_e\rangle\langle F_e| + g\,\big(|F_e\rangle\langle e| + |e\rangle\langle F_e|\big),
\qquad \tau = \frac{\pi}{2g},$$

and the converter is then off. The emitted state lies within $\sqrt{d_a + d_b - 2}/(2g)$ in norm of the pair on $e$
(PROOF.md, section 4). Fuel left in the emitter is a complete outcome. The emitter seam is $(0, 1)$ on one cell and on
20 cells and $(2, 3)$ on 80 and 320 cells; one end of each carries an external seam on every federation. The
coupling is $g = 64$.

## 7. Reader, records and tests

**Observers.** Every seam has two port observers, one at each end. Each copies the predicate
$B_e = \mathbf 1\{n_a = n_b = 1\}$, "both ends hold one load", into its private pointer; a courier carries one copy; the
comparison admits a public record only when both copies are valid and agree. The reader block of one seam has
eleven bits (two pointers with validity, the courier with validity, two decision bits, three public bits); its
four stages are reversible permutations on all 2048 words, and decisions distinguish agreement, mismatch and
unavailable ([reader.py](reader.py)).

**Outcome set.** On the sector of two loads the seam projectors are orthogonal. The complete outcome set is one
agreed record per seam and one outcome "no pair"; separate loads, doubled loads and fuel all read "no pair". On 20
cells that is 630 records and "no pair". On three loads the face reading has one record per face and "no agreed
face".

**Readings.** Each reading uses a fresh run; reading is never used to hold an object in place. The private seam
memory stays unread, and seam records are committed at the reading. The finite readout with the law acting
during the four stages ($\tau_{\rm read} = 1/16384$, 161,792 states on one cell) differs from the ideal reader at the
start of the window by at most $12\,\tau_{\rm read}$ in total variation on every federation.

**What the tests compare.** The state of the simulation does not contain the physical Universe or its spacetime.
Every test of the model is a statement about the probabilities of admissible observer records: here the agreed
seam records and the face records. The state vector, the levels, the evolution parameter and positions in a
sequence of events are tools of the calculation (audit). Free evolution for any duration supplies the state for
the first seam reading. A later reading takes the output of the preceding instrument as its preparation.
A difference between state encodings is no failure of the
model, and a resemblance between them is no success. "Read on other seams" is the observable meaning of motion.

**Particle stability, the observable.** A pair of loads written on one seam by the emitter is read as a pair on
some seam with probability at least one minus the emitter bound at the first complete reading after free evolution; it is read on other
seams than the first; nothing holds it but the law; and the law with the agreement energy removed loses it. Audit
on the spectrum adds that the lowest level of two loads lies below two free loads and is a pair on a seam with
probability above one half. What counts against the claim: "no pair" above the bound at that first reading; no level below
two free loads; a changed law with the same alphabet, routing and symmetry that gives the same answer on the whole
property. [stability.py](stability.py) implements the test; `test_model.py` checks that each part fails under a
changed law, start or reader.

## 8. Events and schedules

**Continuous parameter.** Free evolution is $e^{-itH_\kappa}$, computed by sparse propagation. The theorem is proved
for this form.

**Events.** An event is one seam acting once,

$$E_e(\tau) = \exp\Big(-i\tau\Big[\tfrac12 (I - \mathrm{SWAP}_e) - 2\kappa\, n_a n_b\Big]\Big),$$

on loads 0 and 1. The sum of the event generators over all seams equals $H_\kappa - 6\kappa N$ on $N$ loads 0 and 1. A port observer takes part in the events of its own seams and
counts them; that count is its clock, and no common clock exists. Events on seams that share no port commute, so a
program fixes the records only through the order of events that each port observer sees. Programs
([events.py](events.py), fixed seeds): rounds (every seam once per round in a fresh random order), rounds with two
rates (seams of one half of the carrier twice per round), independent seams (every event draws its seam at
random) and independent seams with two rates.

**Layers.** Layers of disjoint seams (a five-colouring of the 30 seams of the cell, used in every cell, and the
external seams as a sixth layer) run as a symmetric pulse: layers $1..m-1$ for $\tau/2$, layer $m$ for $\tau$, and
back down for $\tau/2$.

**Committed contrast.** A further member of the memory family commits the contrast of the acting seam (which end
holds a load) to a public record after every event.

## 9. The unrecovered seam

A seam whose memory is not blank is unrecovered. By section 3 blankness and the magnitude of the contrast are conserved: a recovered start
stays recovered, and an unrecovered seam is a supplied preparation. With no load anywhere the ten records of a
seam have the levels $0$, $\kappa$, $\kappa + 1$, $4\kappa$ and $4\kappa + 1$, each twice: $0$ for the records of
zero contrast, $\kappa$ and $\kappa + 1$ for the symmetric and the antisymmetric record of contrast $\pm 1$, and
$4\kappa$ and $4\kappa + 1$ for contrast $\pm 2$. A seam that holds a record of size one with empty ports holds no
load and lies $\kappa$ above the blank seam.

**A load beside it.** The local term vanishes on the coherent written code and equals the memory generator on its orthogonal complement in the labelled space. A product of a definite coin label and definite endpoint loads can overlap both subspaces. For example the matrix element from coin 0, contrast +1, loads (0,1) to coin 1, contrast +1, loads (1,0) is −κ/2 − 1/4. A written seam can therefore transfer a load while its contrast magnitude stays fixed. At contrast zero the entire labelled block has zero generator. Beside a record of size one a load on a port
of the seam lies below a load elsewhere, by about $2\kappa$ where the record fits the load and by about $\kappa$
where it does not, and it sits on a port of the seam. On one cell at $\kappa = 8$ a load elsewhere starts at
56.38, two levels lie at 41.94 and six between 49.81 and 50.94 ([collar.py](collar.py)).

**Reading.** A local fuel at level 0 is converted into a contrast record by a pulse of area $\pi/4$ with coupling
64 and the law acting; after $q$ cycles a second pulse with one of four phase settings reads "fuel present". From
$Z(q) = [p_0 - p_\pi - i\,(p_{\pi/2} - p_{3\pi/2})]/4$ the observer reads the phase rate per cycle of the symmetric
and of the antisymmetric record, $r_s$ and $r_a$. Their difference is the cycle in exchange units, and
$r_s/(r_a - r_s)$ is the coupling. The length of a cycle is not given to the observer.

The papers of the theory call the private memory of a seam its collar, and this state an unrecovered collar.

## 10. The presence reading

The reader of section 7 keeps one record per seam. The presence reading keeps one answer for the whole cell. A
tree of reversible OR gates over the thirty reader blocks copies the answer "some seam reports a pair" into a
public register, and the tree and the readers run backwards, so that every seam word returns to blank. With $P_e$
the projector on the pair on seam $e$, $P = \sum_e P_e$ and $Q = I - P$, a positive answer acts on a state $X$ of
the loads as

$$X \mapsto P X P \quad \text{(presence)},\qquad X \mapsto \sum_e P_e X P_e \quad \text{(reader with seam records)},$$

and a negative answer acts as $X \mapsto Q X Q$ for both. The first answer has the same law under both readings,
and the laws of later records differ. One read allocates 720 bits. [presence.py](presence.py) executes both
readings and the laws of two readings in one history, section 8 of [PROOF.md](PROOF.md) gives the instrument with
its error bounds, and [presence_protocol.json](presence_protocol.json) fixes the campaign.

## 11. Resources

Level 0, literal allocation: 240 load qutrits, 630 seam memories of eleven states, 180 reference qutrits with 180
further eleven-state memories, one two-state fuel and 630 reader blocks of eleven bits. Without controller and
feedback probes,

$$D_{\rm alloc} = 3^{240}\; 11^{630}\; 33^{180}\; 2\; 2^{11 \cdot 630}.$$

The sector of two loads has $V(V+1)/2$ states (78 on a cell, 28,920 on 240 ports with doubled loads; 66 and 28,680
without), so the tests are polynomial in the number of ports. The construction is no solver for general
constraint problems and claims no polynomial bound for one.

## 12. Supplied

- the coupling $\kappa$;
- the coherent generator $G_\kappa$ on the private seam memory;
- the reference registers of the degree term;
- the emitter with its fuel, coupling $g$ and duration $\pi/(2g)$;
- the evolution parameter of continuous-time runs and the independently specified schedules and angles of finite-event tests;
- the routing of the external seams of the federation;
- the reader's stage durations;
- the unrecovered seam as a preparation, with its fuel and converter;
- the tree of the presence reading, its wiring and its pulse schedule.
