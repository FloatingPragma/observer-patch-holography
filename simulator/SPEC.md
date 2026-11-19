# Executable specification of the reflected-routing family

This document fixes the carrier, the state, the update with its parameters,
the members of the family and their standing, the exact methods of
evolution, the observers' readings and the supplied parts. The code in this
folder implements every clause, and [test_simulator.py](test_simulator.py)
checks each method against the literal update. The
[journal](../docs/research/laws/JOURNAL.md) records the experiments that set
the standing of each member.

## 1. Carrier

A carrier is a torus with four axes of sizes $N_1, N_2, N_3$ and $P$. The
first three are spatial addresses $r = (x_1, x_2, x_3)$. The fourth is a
private cyclic register of $P$ positions held by every address. Every pair of
an address and a private position holds eight slots, numbered 0 to 7 and
paired as (0, 1) on axis 1, (2, 3) on axis 2, (4, 5) on axis 3 and (6, 7) on
the private axis. Slot $s$ belongs to axis $a(s) = \lfloor s/2 \rfloor$ and
points forward, $\sigma_s = +1$, for even $s$ and backward, $\sigma_s = -1$,
for odd $s$.

The carrier is the hardware of the computation. Space, time, fields and
particles are effective descriptions that observers give of the records they
agree on.

## 2. State

One excitation is a unit vector $\psi_s(r, p)$ in
$\mathbb C^{N_1 N_2 N_3 P \cdot 8}$. Two excitations are a unit vector in the
tensor product of two copies. The update commutes with the exchange of the
two copies, so an exchange-symmetric or antisymmetric preparation keeps its
symmetry.

## 3. Update

One update applies the coin at every position and then routes every slot.
For two excitations the encounter term acts first.

**Coin.** With $|e\rangle = 8^{-1/2}(1, \dots, 1)$ and the equal-slot phase
$\lambda$, $|\lambda| = 1$,

$$C_\lambda = (1 + \lambda)\,|e\rangle\langle e| - I .$$

For $\lambda = 1$ the coin reflects the equal-slot component against its
orthogonal complement.

**Routing.** The amplitude of slot $s$ moves one position along axis $a(s)$
in direction $\sigma_s$ and arrives in the partner slot $s \oplus 1$:

$$(R\psi)_{s \oplus 1}(r, p) = \psi_s\big((r, p) - \sigma_s\,\hat e_{a(s)}\big).$$

**Encounter.** For two excitations,
$E_\chi$ multiplies the amplitude by $e^{i\chi}$ on every configuration in
which the two spatial addresses coincide, whatever their private positions.
This form of the encounter term is a declared choice of this simulator.

$$U_{\lambda} = R\,C_\lambda, \qquad U^{(2)}_{\lambda,\chi} = (U_\lambda \otimes U_\lambda)\,E_\chi .$$

**Access.** A member's access rule says which operations inhabitants can
trigger after a fresh choice. Field-only access gives them local
preparations, local detectors and their own clocks; a message between
laboratories is a free-wave launch. Equipped access adds an addressed relay
that carries a chosen value one address per update.

## 4. Members

| Key | $\lambda$ | $\chi$ | Access | Standing |
| --- | --- | --- | --- | --- |
| `reflected-routing` | 1 | 0 | field-only | candidate: passes steps 1 to 3 |
| `reflected-routing-encounter` | 1 | any nonzero value | field-only | undecided: equal to the candidate on one excitation |
| `equipped-relay` | 1 | 0 | equipped | excluded in step 2 as a law with one causal speed |
| `phase-deformed` | $e^{0.05 i}$ | 0 | field-only | disfavoured in step 1 |
| `phase-i` | $i$ | 0 | field-only | registered control of step 1 |

The encounter member is stored with $\chi = \pi/2$; every value of $\chi$
belongs to the family. `python -m simulator laws` lists the members.

## 5. Exact methods

**Literal update.** [law.py](law.py) applies $C_\lambda$ and $R$ to the
full array of shape $(N_1, N_2, N_3, P, 8)$, or the pair update to the
array of shape $(N_1, N_2, N_3, P, 8)^{2}$.

**Private sector.** The update commutes with translations of the private
register. A state $\psi_s(r, p) = e^{2\pi i m p/P}\chi_s(r)$ keeps its
private order $m$, and routing along the private axis multiplies $\chi$ by
$e^{-i\sigma_s\theta}$ with $\theta = 2\pi m/P$. [sector.py](sector.py)
executes this reduction; it equals the literal update to rounding.

**Plane waves.** On $e^{ik\cdot x}$ the update of the candidate is an
$8 \times 8$ unitary whose mobile band has the eigenvalue $e^{-iE}$ with

$$\cos E = \tfrac14\big(\cos k_1 + \cos k_2 + \cos k_3 + \cos k_4\big).$$

A prepared packet is a finite sum of mobile-band waves, and its state at
update $n$ is the same sum with every wave multiplied by $e^{-inE}$.
[spectral.py](spectral.py) evolves clock doublets, optical pulses and gapless
references this way on carriers of 512 positions per spatial axis and 4096
private positions. Its marginals and four-aperture readings equal those of
the literal update on a small carrier to rounding.

## 6. Readings and observers

**Field.** The scalar field at an address is the equal-slot projection
$\varphi(r) = \langle e|\chi(r)\rangle$. A field released at rest,
$\chi_0 = f\,|e\rangle$, gives $\varphi_n = \sum_K \hat f_K \cos(nE_K)\,e^{iK\cdot r}$
because $\langle e|W(k)^n e\rangle = \cos(nE)$. On a chart of side one with
$h = 1/N$ and private angular wave number $Q$, refinement gives the
Klein–Gordon equation $\varphi_{tt} = \tfrac14\nabla^2\varphi - \tfrac{Q^2}{4}\varphi$
with error of order $h^2$ at fixed chart time. The run of
`python -m simulator run continuum` measures that error at several
resolutions.

**Clocks.** The beat clock is a doublet of private orders $m$ and $m+1$; its
hand is the four-aperture reading
$\arg\big[(p_0 - p_2) + i(p_1 - p_3)\big]$ at the private positions
$0, P/4, P/2, 3P/4$ of one address. The optical clock is a pulse travelling
around the third spatial axis; its reading is the number of turns of its
envelope. A laboratory's own clock counts the updates of its address.

**Observer view.** The replay app shows what observer A at the centre of
the chart can hold: the reading of an address at distance $r$ reaches A at
the signal speed $c = 1/2$, so the picture A holds at its clock time $t$
shows each address at $t - r/c$. The app computes this delay from the
recorded field; the continuum run does not simulate the delivery signals.
The chart description, every address at the same chart time, is an audit
view.

## 7. Run files

A run is one gzip-compressed JSON document with format
`oph-simulator-run/1`: the law with its parameters and standing, the carrier,
the preparation, the recorded readings, the comparison it makes, and the
software and machine that computed it ([runs.py](runs.py)). A run computed on
any machine replays in the app.

## 8. Supplied

- the carrier sizes and the private sector;
- the preparations: field patterns released at rest, clock doublets, optical
  pulses and reference waves with their momentum windows;
- the form of the encounter term and the statistics of two excitations;
- the readers, their apertures and the observer positions;
- the archived free reference speed 0.499350 used in the clock forecasts;
- in the observer view, delivery of every reading at the signal speed.
