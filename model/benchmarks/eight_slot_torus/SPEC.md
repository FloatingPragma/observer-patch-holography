# Executable specification of the torus law

This document fixes the carrier, the state, the torus law with its parameters,
the controls used in its experiments, the exact methods of evolution, the
observers' readings and the supplied parts. The code in this folder
implements every clause, and [test_simulator.py](test_simulator.py) checks
each method against the literal update. The [benchmark results](README.md) describe the measured clock and propagation records on these supplied axes. The long-wave dispersion follows from the Fourier coin-and-shift block of this specified architecture.

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
For two constituents the capture vertex acts first, then the contact update.

**Coin.** With $|e\rangle = 8^{-1/2}(1, \dots, 1)$ and the equal-slot phase
$\lambda$, $|\lambda| = 1$,

$$C_\lambda = (1 + \lambda)\,|e\rangle\langle e| - I .$$

For $\lambda = 1$ the coin reflects the equal-slot component against its
orthogonal complement.

**Routing.** The amplitude of slot $s$ moves one position along axis $a(s)$
in direction $\sigma_s$ and arrives in the partner slot $s \oplus 1$:

$$(R\psi)_{s \oplus 1}(r, p) = \psi_s\big((r, p) - \sigma_s\,\hat e_{a(s)}\big).$$

$$U_{\lambda} = R\,C_\lambda .$$

The eight-slot torus law has $\lambda = 1$ for every constituent.

**Contact.** The law acts on two distinguishable constituents, one of
which carries a conserved rail label. Off contact each constituent applies
$C_1$. Where both occupy the same address and private position, the eight
matched slot pairs $(s, s)$ receive the contact coin

$$D_\varphi = e^{i\varphi/2}\big[(1 + e^{-i\varphi})\,|e\rangle\langle e| - I\big]$$

and the phase $e^{-i\epsilon}$, while the 56 unmatched slot pairs there are
left unchanged. Both constituents then route. Rail $r$ has the contact phase
$\varphi_r$ and the extra phase $\epsilon_r$. The matched coincident states
form an invariant subspace, and so does its complement: a pair started
outside never enters it. On the subspace the pair is one walker with the coin
$D_\varphi$, whose mobile band obeys $\cos E = a(k)\cos(\varphi/2)$, $a$ the
mean of the four cosines.

**Capture.** The law gives every address a periodic wire of
neutral modes, longer than the run. Before the contact update, where the two
constituents meet with constituent A in slot 5 and B in slot 4 (a head-on
collision along the third axis) and the start of that address's wire is
empty, the vertex rotates the collision state into the matched equal-slot
pair plus one marker at the start of the wire:

$$|5, 4; 0\rangle \mapsto \cos\theta\,|5, 4; 0\rangle + \sin\theta\sum_s e_s\,|s, s; 1\rangle,
\qquad \theta = \pi/4,$$

and acts as the inverse rotation on the orthogonal partner. Every marker
moves one link along its wire per update. Both rails see the same vertex.

**Opening.** A newborn pair carries an age $j$, born at zero with its marker
and advanced by one at every update; the age register has 8192 states. On
the bound sector the update is the contact coin followed by the partial hop

$$H_\theta = e^{-i\theta}\big[\cos\theta\,I + i\sin\theta\,F\big],
\qquad \theta = \tfrac{\pi}{2}\,\frac{1 - \cos\big(\pi (j+1)/128\big)}{2}\ \text{for}\ j < 128,$$

between the two ends of every routed edge, $F$ the exchange of those ends
with the flip of the routing label. From age 128 on, $\theta = \pi/2$ and
$H_\theta = F$, the full hop. The opening has an exact inverse.

**Access.** Inhabitants have local preparations, local detectors and their
own clocks; a message between laboratories is a free-wave launch.

## 4. The law and the controls

`python -m model.benchmarks.eight_slot_torus laws` lists them. The eight-slot torus law is stored as `torus`.
Every other entry is a control used in the experiments, and changes one part
of the law:

| Key | What it changes |
| --- | --- |
| `control-encounter` | a scalar phase $\chi = \pi/2$ when two constituents share an address, in place of the contact coin |
| `control-relay` | an addressed relay available to inhabitants, one address per update |
| `control-phase` | equal-slot phase $e^{0.05 i}$ |
| `control-phase-i` | equal-slot phase $i$ |
| `control-spectator` | equal contact phases $\pi/16 + \pi/4096$ with extra phases $\mp\pi/8192$: a clock that counts services |
| `control-abrupt` | the full hop from age zero, with no opening |

## 5. Exact methods

**Literal update.** [law.py](law.py) applies $C_\lambda$ and $R$ to the
full array of shape $(N_1, N_2, N_3, P, 8)$, or the pair update to the
array of shape $(N_1, N_2, N_3, P, 8)^{2}$.

**Private sector.** The update commutes with translations of the private
register. A state $\psi_s(r, p) = e^{2\pi i m p/P}\chi_s(r)$ keeps its
private order $m$, and routing along the private axis multiplies $\chi$ by
$e^{-i\sigma_s\theta}$ with $\theta = 2\pi m/P$. [sector.py](sector.py)
executes this reduction; it equals the literal update to rounding.

**Plane waves.** On $e^{ik\cdot x}$ the update of one constituent is an
$8 \times 8$ unitary whose mobile band has the eigenvalue $e^{-iE}$ with

$$\cos E = \tfrac14\big(\cos k_1 + \cos k_2 + \cos k_3 + \cos k_4\big).$$

A prepared packet is a finite sum of mobile-band waves, and its state at
update $n$ is the same sum with every wave multiplied by $e^{-inE}$.
[spectral.py](spectral.py) evolves clock doublets, optical pulses and gapless
references this way on carriers of 512 positions per spatial axis and 4096
private positions. Its marginals and four-aperture readings equal those of
the literal update on a small carrier to rounding.

**Bound sector and capture branches.** [law.py](law.py) executes the pair
laws literally on two constituents (`step_contact`) and on the bound sector
as one walker (`step_bound`); the two agree to rounding. A marker never
returns within a run, so each capture event, labelled by its address and
update, starts a branch orthogonal to every other, and the branch lives in
the bound sector. [formation.py](formation.py) evolves the unmarked
constituents and every branch this way. The test suite compares its capture
probability with a literal update that carries the marker wires as explicit
occupation bits.

## 6. Readings and observers

**Field.** The scalar field at an address is the equal-slot projection
$\varphi(r) = \langle e|\chi(r)\rangle$. A field released at rest,
$\chi_0 = f\,|e\rangle$, gives $\varphi_n = \sum_K \hat f_K \cos(nE_K)\,e^{iK\cdot r}$
because $\langle e|W(k)^n e\rangle = \cos(nE)$. On a chart of side one with
$h = 1/N$ and private angular wave number $Q$, refinement gives the
Klein–Gordon equation $\varphi_{tt} = \tfrac14\nabla^2\varphi - \tfrac{Q^2}{4}\varphi$
with error of order $h^2$ at fixed chart time. The run of
`python -m model.benchmarks.eight_slot_torus run continuum` measures that error at several
resolutions.

**Clocks.** The beat clock is a doublet of private orders $m$ and $m+1$; its
hand is the four-aperture reading
$\arg\big[(p_0 - p_2) + i(p_1 - p_3)\big]$ at the private positions
$0, P/4, P/2, 3P/4$ of one address. The optical clock is a pulse travelling
around the third spatial axis; its reading is the number of turns of its
envelope. The rail clock of a bound pair is read by the four-outcome
analyzer $p_a = \tfrac14\sum_s |r_{0,s} + a\,r_{1,s}|^2$, $a \in \{1, i, -1, -i\}$,
on the two rail amplitudes at one address, with the hand
$\arg\big[(p_1 - p_{-1}) + i(p_i - p_{-i})\big]$. A laboratory's own clock
counts the updates of its address.

**Observer view.** The replay app shows what observer A at the centre of
the chart can hold: the reading of an address at distance $r$ reaches A at
the signal speed $c = 1/2$, so the picture A holds at its clock time $t$
shows each address at $t - r/c$. The app computes this delay from the
recorded field; the continuum run does not simulate the delivery signals.
The chart description, every address at the same chart time, is an audit
view.

## 7. Run files

A run is one gzip-compressed JSON document with format
`eight-slot-torus-run/1`: the law with its parameters and standing, the carrier,
the preparation, the recorded readings, the comparison it makes, and the
software and machine that computed it ([runs.py](runs.py)). A run computed on
any machine replays in the app.

## 8. Supplied

- the carrier sizes and the private sector;
- the preparations: field patterns released at rest, clock doublets, optical
  pulses and reference waves with their momentum windows;
- the statistics of two constituents and the rail preparation;
- the readers, their apertures and the observer positions;
- the archived free reference speed 0.499350 used in the clock forecasts;
- in the observer view, delivery of every reading at the signal speed.
