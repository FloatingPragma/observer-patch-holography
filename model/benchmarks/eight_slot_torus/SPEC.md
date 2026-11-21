# Executable specification of the torus law

This document fixes the carrier, the state, the torus law with its parameters,
the controls used in its experiments, the exact methods of evolution, the
observers' readings and the supplied parts. The code in this folder
implements every clause, and [test_simulator.py](test_simulator.py) checks
each method against the literal update. The [law](../model/README.md) page
states what is proven and measured.

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

### What a pair detector retains after a capture

Fix rail zero, private order zero and a periodic strip or cube of volume $V$.
Write $S$ for the equal-slot vector. At each translated launch site prepare
the sharp inward approach to the $(5,4)$ vertex, with amplitude
$f(x)=e^{ibk\cdot x}(1+e^{i(k\cdot x+s)})/\sqrt{2V}$ for a nonconstant
grating. One free routing brings the separated constituents to contact.
The vertex writes with probability $1/2$. In its positive branch the joint
state is $\sum_x f(x)|x,S\rangle|m_x\rangle$, where the outgoing wires
$m_x$ are orthonormal. Tracing the wires gives diagonal site weights
$|f(x)|^2$. Every material-only read is therefore invariant under a phase
mask on the incident sites, for any subsequent material unitary. A sharp
site mixture has exactly these same material records.

There is a resource bound on undoing this loss. An accepted wire effect $E$
which restores the undistorted material state for **every** superposition
must satisfy $\langle m_y|E|m_x\rangle=c$ for all $x,y$, up to a fixed
correctable phase mask. This follows by testing each site and every
two-site superposition. Its nonzero eigenvalue on their span is $Vc\le1$.
Projection onto the uniform wire state attains $c=1/V$ conditional on
capture. This bound concerns one accepted outcome with one fixed
correction. A full measurement of the wires in a conjugate basis with
outcome-dependent feedforward is a different apparatus.

One equipped receiver counts all slots at a single aperture. Four incident
phases give differences of opposite click counts divided by their sums.
The real difference uses $0,\pi$, the imaginary difference $3\pi/2,\pi/2$.
For a coherent age-zero input this reads $\langle u_n(bk),u_n((b+1)k)\rangle$.
For a capture it reads $\sum_x G_n(x)e^{-ik\cdot x}$, where $G_n$ is the
position probability from a point write. Free space coordinates enter the
calculation as indices of the apparatus. The receiver needs only counts,
phase labels and its optical booking. Post-run enrollment and the
local hit give the three joint outcomes. A coherent comparison supplies the
matched grating with fair enrollment. Disabling the vertex gives all refusals.

### The full age period and the coherent envelope

Put $L=8192$, $O=128$, $\mu=\pi/32$ and
$a=(1+\sum_{i=1}^3\cos k_i)/4$. For reflected routing $F$ the vectors
$S$ and $(FS-aS)/\sqrt{1-a^2}$ span an invariant space. At $k=0$ the
second ray is understood along the incident direction. In this basis

$$
 F=\begin{pmatrix}a&\sqrt{1-a^2}\\\sqrt{1-a^2}&-a\end{pmatrix},\qquad
 C=\operatorname{diag}(e^{-i\mu},-e^{i\mu}).
$$

The gate at opening age $j$ is
$M_j=e^{-i\theta_j}(\cos\theta_j I+i\sin\theta_j F)C$.
Write $P_r$ for the first $r$ gates, using $FC$ for ages $O$ through
$L-1$. At $n=\ell L+r$ the wave is exactly $u_n=P_rT^\ell S$ with
$T=(FC)^{L-O}M_{O-1}\cdots M_0$. Every wrap repeats the opening.
The cosine opening gives $\sum_{j<O}(\pi/2-\theta_j)=127\pi/4$.
The two rest eigenvalues are $z_0=e^{-iL\mu}$ and $z_1=-i e^{iL\mu}$.

For fixed finite Fourier menus $k=hq$ on a unit chart, the registered
rail has $z_0=1$, $z_1=-i$. Near rest the off-diagonal term of $F$ is
$|k|/2+O(|k|^3)$ and $F_{00}=1-|k|^2/8+O(|k|^4)$. With $T_j$ the
coefficients at orders $|k|^j$, the isolated bright eigenvalue has
quadratic coefficient $A=(T_2)_{00}+(T_1)_{01}(T_1)_{10}/(z_0-z_1)$,
by the two-row eigenvalue equation. Conjugating by $\operatorname{diag}(1,-1)$
reverses the signed radial variable, so the eigenvalue is even along each
ray. The quadratic coefficient is isotropic. The overlap error in the
continued bright ray is $O(|k|)$. After removing the common rest phase,
the coherent grating thus freezes at $n\sim t/h$ and has Schrödinger
multiplier $\exp(-itD|q|^2)$ at $n\sim t/h^2$, uniformly for bounded
$t$ on the finite menu, with $D=-\operatorname{Im}(A/z_0)/L=1.253192337$.
The single permanently mature coin would give $\cot(\mu)/8=1.269146298$.

For a finite-mode error budget let $\lambda$ and $\Pi$ be the continued
bright eigenvalue and projector. Put $\nu=D|k|^2$. The vector error from
the plane-wave forecast is at most
$\|P_r S-e^{-ir(\mu+\nu)}S\|+2\|(I-\Pi)S\|+
\ell|\lambda-z_0 e^{-iL\nu}|$. Splitting off the complementary
projector and telescoping the scalar powers proves the bound. Summing the
two vector errors bounds a coherent grating contrast. The replay reports
the finite bounds and observed forecast errors at both diffusive windows.
[Period operator and bounds](formed_clocks/fringe_walk.py)

### Which change produces a scalar wave

The rest eigenvalues coincide precisely at $e^{2iL\mu}=i$. The smallest
positive resonance is $\mu_*=\pi/(4L)$. At $\mu=\mu_*+hm$ expand
$T/z_0(\mu_*)=I-ihL H(q)+O(h^2)$ for $k=hq$. Its two diagonal entries
in $H$ are $m,-m$, by differentiating the rest eigenvalues. The linear
routing term is purely off-diagonal with magnitude $c_*|q|$,
$c_*=|(T_1)_{10}|/L$. Unitarity makes $H$ Hermitian. Thus
$H^2=(m^2+c_*^2|q|^2)I$. At complete periods the product converges to
$e^{-itH}$ for $hn\to t$ by telescoping unitary factors. The error on
a bounded time and finite Fourier menu is $O(h+|hn-t|)$ with fixed $L$.
The scalar projection has initial value $f$ and derivative $-im f$.
Its second derivative is $c_*^2\Delta f-m^2 f$ at time zero, and the same
wave equation holds throughout. An intensity reader in all slots includes
the longitudinal component, which the calculation retains.

Each opening gate can contribute at most $\sin(\theta_j)/2$ to the
first derivative. Each mature gate can contribute at most $1/2$.
Since the rest gates are unitary, summing yields
$c_*\le (L-O+\sum_{j<O}\sin\theta_j)/(2L)=0.4969226637<1/2$
at **every** resonance. At $\mu_*$ the mature contributions alone sum
as a geometric progression with successive phase $2\mu_*$. Subtracting
the opening bound gives $c_*\ge
(\sin((L-O)\mu_*)/\sin\mu_*-O)/(2L)=0.4367876296$. The full product
gives $0.4448918621$. The detuning $m$ is supplied. Sending the bare
phase to zero as $h$ leaves the two rest eigenvalues separated and gives
a frozen ballistic grating. A slow envelope at fixed rest phase by itself
cannot discriminate a massive relativistic theory from its low-momentum
limit. The two-sector test compares the measured phase rates of two
waves with identical grating and clock. [Resonance](formed_clocks/fringe_resonance.py)
· [Common aperture](formed_clocks/fringe_reference.py)

### Ballistic propagation from the marked birth

A point write samples the entire Brillouin zone. For a strip grating of
one turn, the marked contrast at $n/N\to t$ has real limit

$$
 C(t)=\int_0^\pi \cos(2\pi t\,v(k))\,\frac{dk}{\pi},\qquad
 v(k)=\frac{1}{L}\,\partial_k\omega(k).
$$

Here $e^{-i\omega(k)}$ is either eigenvalue of the determinant-normalized
axial period up to its constant phase. To obtain the formula, Fourier
transform the position density. Its first harmonic pairs adjacent
frequencies $k,k+2\pi/N$. In a smooth part of the zone, projectors on
different bands at adjacent momenta have product $O(1/N)$. Same-band
products tend to their weights on $S$ and the phase difference tends
to $2\pi t\,\omega'(k)/L$. The same tail $P_r$ on the two columns
cancels in this limit because $r<L$, the slots vary smoothly away from
rest and each $P_r$ is unitary. The trace discriminant is real analytic
and nonzero at rest, so isolated crossings and the interval endpoints
can be covered by arbitrarily small intervals of bounded integrand.
Riemann sums on the rest converge. The constant determinant gives
opposite velocities, the weights sum to one and spatial reflection makes
the full contrast real. This yields the displayed cosine average.

A single undamped line with unit visibility obeys
$C(2t)=2C(t)^2-1$. For the limiting mixture the excess equals
$2\operatorname{Var}(\cos(2\pi t v))$. Observers can test this identity
at finite size from the local repeat bookings alone. It concerns the
specified marked preparation and reader. Multiple velocities or a continuum
field with a distribution of occupied modes can produce such an excess.
[Marked flight and quadrature](formed_clocks/fringe_flight.py)

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
