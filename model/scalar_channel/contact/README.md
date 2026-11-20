# A local interaction between contrast memory and scalar occupation

A collar changes a scalar fuel reading. Filling the scalar source changes the
phase read from the collar. A mobile load pair and a scalar quantum also give
correlated readings at two local pointers. All three experiments use the same
Hamiltonian on one cell. The largest distance between the mobile joint record
law and the product of its marginals is 0.03634 for positive contact and
0.03127 for negative contact. At zero contact it vanishes exactly.

## Law and source map

The experimental member `contrast-square-scalar-contact-v1` takes the tensor
product of the two register families in the [construction manual](../../construction_manual/README.md).
Each port carries its load, each seam its eleven-state memory, and the scalar
even algebra uses the 432 link registers. In the manual's port labels set

$$H(g)=H_{\rm contrast}(8)\otimes I+I\otimes d\Gamma(D)
       +g\sum_s A_s\otimes n_s .$$

Here $A_s=U_s(0\oplus D_{\rm memory}^2)U_s^*$
is the coefficient of $\kappa$ in the complete
[99-state seam term](../../contrast_memory_spherical_v1/law.py), with identities on loads in
the memory generator. $n_s=(I-B_s)/2$ is the [encoded seam occupation](../compiler.py).
Its parity word commutes with the scalar loops. Each product uses the same
seam on both sides. Reversing the seam exchanges its ports and sends
$(c,d)$ to $(1-c,-d)$, preserving the echoed square. The
[port chart](../../contrast_replay.py) maps the two carrier labellings.
The degree compensation of the contrast source is unaltered. Energy is in
exchange units, with scalar rest level 1 and hop $1/\sqrt2$. That relative
normalization, the contact operator, the apparatus and the signs tested
($g=0,+1,-1$) are design choices.

On the empty scalar register, with no object on any seam or wedge, the full
contrast source is exact. On empty loads and blank memories every $A_s$
vanishes and the full scalar source is exact. The prepared reference of the
scalar channel is not empty: it holds occupation 0.2635 on every seam, so in
that reference the contact acts on every load.
These embeddings also hold with untouched quantum references. In particular
the one-load cell transport carries the port-Laplacian fingerprints
$5-\sqrt5$ and 6 after the native factor of one half and the constant
reference shift. The free scalar cage has line $2\sqrt6$. Every value of
$g$ has the same restrictions to these two sectors. Even their complete
record laws, for instruments staying in those sectors, cannot determine $g$.

## Two directions of the interaction

With empty loads the blank and the written states
$(|0,+1\rangle\pm|0,-1\rangle)/\sqrt2$ are eigenstates of $A_s$ with
values 0 and 1 respectively for blank and written. Both written signs are
orthogonal to the write image at empty loads. Their exchange energies differ
by one. With every scalar source mode vacant or every mode filled, hopping
annihilates the preparation and the written arm gains $g n$ for filling
$n=0,1$. The difference between the two symmetric-arm rates divided by the
vacant antisymmetric-minus-symmetric rate therefore gives $g$. Two finite
collar pulses and four phase settings read these rates. With the declared
phase branch, the largest phase per query is $10/16<\pi$. The observed
ratios are 0.00007, 0.99981 and −0.99973. The same two-arm restriction
intertwines inclusion of any spectator cells prepared in the same filling,
up to a common scalar phase. This particular local comparison is invariant
at every size for *every* $g$. It gives no selection of a running coupling.

For the reverse experiment the source begins with all 120 wedge modes filled
and every seam mode empty. Hopping is off at this point of the specified
ramp. A finite strength-64 pulse on seam $(0,2)$ writes the symmetric collar
and its fuel is read. Its transfer probability is
$(64/\Omega)^2\sin^2(\pi\Omega/128)$, $\Omega=\sqrt{64^2+4^2}$.
The two records select invariant sectors of the same law. A missed write
gives scalar generator $D$ and a successful one gives $D+g|s\rangle\langle s|$.
The contact acts throughout all 192 ramp substeps, every ring pulse, wait
and capture. The occupied frame includes all 150 source modes. Each of the
60 cage settings retains the write fuel and both copies of the scalar fuel.
The largest changes in the marginal scalar reading relative to zero contact
are 0.64345 and 0.64862. The sharp and undriven preparations last 0.28679
and 0.58947 exchange-time units. Their contrast therefore includes a change
in preparation time.

The coherent pulse and the undriven control both last 0.58947. Subtract their
full-fuel frequencies at the same wait and analyzer, then subtract that
difference for zero contact. The largest absolute differences are 0.621752
and 0.628802 for positive and negative contact. This secondary reading uses
the same retained counts. For one million simulated trials per setting,
Hoeffding's inequality and a union bound cover all $K=180$ full-fuel
frequencies simultaneously. Each difference uses four frequencies, giving
an absolute sampling allowance of
$4\sqrt{\log(2K/10^{-8})/(2\times10^6)}<0.01395$ with coverage at least
$1-10^{-8}$, including selection of the maximum. The allowance describes
sampling from the declared law. These are readings of fuel pointers.

## Exactly which cages survive

Take the [reducing pair of a degree-five ring](../cage/PROOF.md) with unit
seam profile $e$ and wedge partner $e'$. If stationary collars produce
seam potentials $w_s$, the term sends $e'$ to zero. The squared norm sent
from $e$ outside the pair is

$$\|(I-P)Ve\|^2=\tfrac14\sum_{s\in r}w_s^2-
                    \big(\tfrac14\sum_{s\in r}w_s\big)^2.$$

This is the variance on the four seams, since $|e_s|^2=1/4$ there. The
pair survives exactly when its four potentials agree. A single marked seam
gives $3g^2/16$ on each ring through it and zero on all other rings. This
local identity also holds on every federation ring satisfying the degree
hypothesis. It concerns leaving the two modes, without implying physical
motion. On one cell the equalities on the 30 rings have rank 29 on 30 seam
weights. Preserving every ring requires a uniform potential. Requiring it
for every placement of a single written collar forces $g=0$.

There is a quantitative cost to omitting the term. Starting from $e$, the
probability outside the pair after a short free interval $t$ has leading
term $3g^2t^2/16$ for an on-ring contact. On the binary mobile sector with
one scalar quantum, the difference of the full and additive generators has
norm at most $|g|$. Duhamel's formula bounds the total variation of any
terminal pointer law by $\min(1, |g|T)$ for the same input and an identical
apparatus sequence of total action time $T$. Neither estimate bounds an
unspecified face or rhombus term.

Precisely 17 of the 18 modes at each of $\pm\sqrt6$ survive a single seam
of nonzero strength. Let $P_E$ project onto either old eigenspace. A ring
through $s$ has a nonzero seam component there, so $P_E|s\rangle\ne0$.
Project the perturbed eigenvalue equation at $E$ with $P_E$. This forces
$\langle s|v\rangle=0$, after which the original equation holds. The
remaining space is the kernel of one nonzero functional on the old space.

## A mobile pair under the same law

On blank memory the [contrast identity](../../contrast_memory_spherical_v1/PROOF.md)
gives $A_{ab}=(n_a-n_b)^2$. Use the 78 load configurations, one residual
emitter fuel state and all 150 one-quantum scalar modes (11,850 amplitudes).
An isolated occupied apparatus mode supplies the even-parity companion.
The scalar starts on $(0,1)$. The resonant load emitter on $(0,2)$ has
strength 64. The contact is zero in its fuel branch, where all loads vanish.
At six waits two initially blank bits undergo simultaneous strength-64
controlled $X$ pulses of duration $\pi/128$, controlled by both end loads
being one on $(0,2)$ and by scalar occupation on $(0,1)$ respectively.
The source runs throughout. These local binary readers are supplied apparatus.
The record consists of the two terminal bits. Orthogonal pointer-$X$ sectors
give four source propagators whose signed sums yield the full record law.
At $g=0$ the initial state is a product. Each pulse, wait and terminal
projection factors between the load, its apparatus and the scalar with its
own apparatus. Hence the two terminal bits are independent at every wait.

There is also an all-time bound before the terminal read. In the binary
two-load sector, on a single scalar quantum, the contact is diagonal with
spectrum between $\min(0,g)$ and $\max(0,g)$. Doubled loads are invariant
and the emitter cannot reach them. The scalar width is $2\sqrt6$ by the
incidence identity. Thus the native cell's adjacent and separated diagonal
blocks retain ordered gap at least
$\gamma=2\kappa-1-3\sqrt6-|g|$, and their off-diagonal norm is 2.
[The block separation theorem](../../contrast_memory_spherical_v1/PROOF.md) bounds
the separated amplitude after an ideal adjacent start by $4/\gamma$.
During emission compare to a perfect two-level rotation tensored with the
free scalar propagator. In the emitted arm the contact residual has norm at
most $|g|$ and the eight outgoing swaps of size one half occupy orthogonal
load labels distinct from that arm. Integrating the sine of the rotation
bounds the preparation error by $\sqrt{2+g^2}/64$. Consequently

$$\Pr(\text{separated, doubled or fuel})\le
   \min\{1,(4/\gamma+\sqrt{2+g^2}/64)^2\}\quad(\gamma>0).$$

For either signed contact this is below 0.395. The largest simulated weight
at the six waits is below 0.0122. These weights audit the state-domain bound.
The separately computed two-bit probabilities are the observable prediction.
Multiple active scalar quanta require a different bandwidth estimate.

## A material clock and a travelling pulse on the same cell

Six stationary empty-load collars occupy seams 0, 5, 11, 16, 22 and 29 in the
scalar carrier's seam ordering. Both written exchange eigenstates have contrast
square one. At station $i$ the complete written doublet therefore adds $n_i$ to
the scalar generator, including when a local clock is precessing in that doublet.
The material splitting is one exchange unit. The seven collar preparations,
the read schedule, the source energy unit, $\kappa=8$, and the homogeneous
contact $g=+1$ are supplied. Every active scalar seam and wedge is retained.
One occupied inert companion fixes the even parity of the code. Preparing an
excitation at a named station and measuring its occupation are apparatus services.
[Source map](joint_source.py) · [Experiment specification](radar_protocol.json)

Let $P=|m\rangle\langle m|$ be the occupation of the seventh seam, $m=2$.
On the blank and symmetric written arms there the scalar generators are
$H_0=D+\sum_i n_i$ and $H_1=H_0+P$, after removing scalar-independent
energies. In the full clock evolution the written arm also accumulates its
native energy 8. Suppression and reversal of the marked contact keep all six
station potentials fixed.

**One local intervention changes both records.** Prepare the marked collar in
equal arms and a scalar in one of seven local modes. After two exchange-time
units, apply the scalar occupation term $64P$ for $\pi/64$ or wait for the
same interval. The source remains on throughout. Read a collar quadrature and
the scalar station on the same trial, retaining a seventh scalar outcome
for a detection elsewhere. This tests a pulse that a station can apply with
the Hamiltonian fixed between interventions. It costs a supplied phase drive.
For a scalar launched at the mark, six units after the drive the minus-Y
positive-pointer probability changed by $-0.18631$ and the probability of a
scalar at seam 0 by $-0.05435$. Subtracting the same pulse-minus-sham differences
with the marked contact suppressed gave $-0.18631$ and $-0.06169$.
Simultaneous binomial intervals over the entire menu, including selection
of this setting, give absolute lower bounds 0.18595 and 0.05422 for both
comparisons at family error $10^{-6}$. The $\pi/64$ drive was fixed before
the simulated counts. These numbers measure changes in pointer probabilities.
[Local intervention and reader](intervention.py)

**Which excitations can reveal the contact?** The cyclic subspace
$K=\mathrm{span}\{H_0^j|m\rangle:j\geq0\}$ has dimension 37. On its
113-dimensional orthogonal complement, $e^{-itH_0}$ and $e^{-itH_1}$ coincide
for every time. Conversely, equality for a given input at every time gives
$P e^{-itH_0}\psi=0$ by differentiation. All its derivatives vanish at
zero, hence $\psi\perp K$. The executable proves the dimension by exact
integer orthogonalization in a similarity chart whose wedge metric is 2.
If $k_a$ is the first order with $c_a=\langle m|H_0^{k_a}|a\rangle\ne0$,
the first changed amplitude at $b$ has coefficient
$c_b c_a/(k_a+k_b+1)!$ in powers of $-it$. In the difference of powers
expand each word containing $P$. The first surviving word has one $P$,
with $k_a$ factors on its right and $k_b$ on its left.
The bare-phase-corrected marked-collar coherence obeys
$Z'_a=-i\langle a|e^{itH_0}|m\rangle\langle m|e^{-itH_1}|a\rangle$.
Its first imaginary response has order $2k_a+1$. The six station orders
$k_a$ are 2, 4, 4, 4, 4 and 6. This gives an all-time obstruction for
invisible excitations and a source prediction for each visible station.
[Exact calculation](contact_response.py)

**Selecting between five record laws.** A separate finite strength-64
preparation and readout retains joint records with two fuel outcomes and all
150 scalar modes. It tests both written signs. A sparse propagation and
an independent spectral sum agree within $1.9\times10^{-13}$ per outcome.
Five stations at waits 0.5 and 2 select the positive marked contact over a
suppressed or reversed mark, a product of independent records, and a
dephased mark. With no refit the largest total-variation discrepancies for
the selected forecast are 0.00382 at station six, 0.00405 at wait 6 on the
first five, and 0.00479 at the marked reference. The simultaneous allowance
is 0.01445. The samples were generated by this contact Hamiltonian.
The selection measures the power of the apparatus to distinguish these
five hypotheses. [Joint laws, alternative hypotheses and scoring](joint_clock.py)

## Six-station radar and receiving clocks

In each ordered pair the receiving detector checks a local occupation at
half-unit intervals. Its first hit starts a return from that station. Each
leg has twelve opportunities. Both timeouts enter the accounting. For a
fixed written or blank mark the rank-one hit makes the two first-hit
amplitudes multiply. The renewal recursion predicts each amplitude from
unmonitored propagators. In 31,457,280 trials per mark preparation the
blank and written menus yielded 6,517 and 6,445 completed trips. Each
reported range is half the conditional mean successful duration in
cadence units. It is available only alongside its success probability.
[Radar instruments, holdouts and uncertainty](radar.py)

A coherent mark requires the conditional material state from the outward
leg. A return click swaps the scalar into a local isolated capture mode.
For a timeout the residual scalar continues freely to the common deadline.
For each of 157 histories $h$, let $v_0(h),v_1(h)$ be its two subnormalized
scalar vectors at that deadline, including the free collar phase. Its X
probabilities are $\|v_0\pm v_1\|^2/4$, with $\pm i$ giving the other
quadrature. Replacing the mark independently at turnaround changes the
scalar history distribution by at most $6.73\times10^{-7}$ in the 30
tested directions. The joint terminal instrument distinguishes initial
coherence from dephasing by up to 0.48921 in total variation. The first
number quantifies how little the scalar-only history reveals about this
distinction. [Coherent return and reset control](coherent_radar.py)

Separate receiving-clock cohorts read the written doublet only when the
signal arrives. The two signed pointers use plus-X and minus-Y, with zero
on a timeout. For a slot displaced by $s$, the complex mean is
$e^{i\omega s\Delta}\sum_n p_n e^{i\omega n\Delta}$, where $p_n$ is the
unconditional first-hit probability and $\Delta=1/2$. Independent
sender-clock records supply the phase rotation between slots. No carrier
coordinate or arrival timestamp enters this comparison. The unit-rate
receivers passed every tested comparison. Detuning only the receiver by
25% excluded the sender prediction on 22 of 30 paths from slot 0 to 16
for each mark preparation. The much shorter interval from 16 to 17
resolved no such exclusion. The clock measures an ensemble phase modulo a
turn. The check cadence, common launch slots and preparation of fresh local
clocks remain supplied. [Station clock acquisition and reader](station_clock.py)

**No three-dimensional Euclidean placement fits these ranges to 5%.**
For any five of the six stations, ten squared distances determine the Gram
matrix of four displacements. A three-dimensional placement requires that
matrix to have rank at most three. Vary each of the ten distances by the
allowed relative error. Its squared interval gives a box with 1,024 corners.
The Gram matrix is affine in the squared distances. If its four leading
principal minors are positive at every corner, the whole box consists of
positive definite matrices. One such subset excludes the six-station fit.

The executable encloses the unmonitored propagators by integer moments and
a degree-240 Taylor sum. The physical generator has norm below 8, so the
unitary integral remainder up to time 6 is at most $48^{241}/241!$.
Outward-rounded renewal amplitudes give success and conditional time
intervals, including denominators bounded away from zero. The determinant
check uses integer arithmetic on outward boxes. At 5%, all six five-station
subsets exclude the fit with the mark blank and three do with it written.
The 7.5% test yields no global exclusion. The original finite counts lacked
the precision for the 5% decision. A separately frozen cohort of $2^{34}$
trials per ordered pair and preparation excludes the fit on five and three
subsets, respectively, at simultaneous confidence $1-10^{-6}$. It also
excludes a single multiplier relating all ranges across the two preparations.
The cohort samples multinomial histograms of simulated records. These bounds
apply to this detector, deadline and conditional mean. [Range enclosures](certified_ranges.py)
· [Gram certificates](metric_certificate.py) · [Acquisition and controls](radar_checks.py)

**Which existing objects survive the intervention?** A four-seam cage with
$k$ of its four seams contacted has squared leakage $k(4-k)/16$ out of
its two-dimensional subspace. On the six-station background, 12 of the 30
rings remain invariant with the mark blank and 10 with it written. Ten
survive both. Each one misses the mark, so the marked memory produces
identical scalar evolution on its two arms for every state in that cage.
An exact line and indefinite persistence on those ten objects alone cannot
diagnose the reciprocal interaction. [Ring classification](radar_checks.py)

**A returned click can come from the sender's own pulse.** A separate
challenge encodes one bit in the launch time. The receiver records the first
arrival slot, captures the scalar, then releases it in the reply slot for
that value. In the silent control it keeps the scalar. If the receiver
missed, the outgoing wave continues into the sender's detection window
under both controls. The complete histories include wrong slots, both
timeouts and these false echoes. Across sixty directed menus the
unconditional false-echo probabilities range from 0.3691 to 0.7426.
Correct delivery through a receiver capture ranges from 0.0000155 to
0.000697. The million-trial menus resolve no positive lower bound for
receiver-tag information after conditioning on the source bit and capture.
Even on the best tested path, fewer than 0.224% of sender clicks followed
a receiver capture. [Time-bin law and silent control](tagged_radar.py)

**A receiver can still transmit a private bit.** The written memory
$|0,+1\rangle$ has equal overlap with the two exchange signs. A local sign
read selects whether the receiver keeps or releases a captured scalar for a
whole block. Both signs leave the endpoint contrast square at one. At the
sender, two labelled calibration blocks give the silent and active click
frequencies. Eight further blocks carry fresh receiver signs. The sender
gets only its ten click totals and the common block length. A midpoint
decoder abstains if the calibration intervals, enlarged for future sampling
error, overlap. False echoes and all timeouts contribute to these totals.

The sender recovered all 480 held-out bits over the 30 ordered paths and two
mark preparations, with zero abstentions. At the weakest pair, activity
changes the sender click probability by only 0.00003196. Each calibration
or bit uses $2^{42}$ independent trials. The archive samples their full
multinomial counts. Across all 600 blocks this represents about
$2.64\times10^{15}$ trials. A simultaneous Hoeffding bound including
acquisition bias limits the probability of any erroneous, non-abstaining
bit to $10^{-6}$ under the fixed independent-trial law. A 43-bit local
counter suffices for a block. Synchronized slots, fresh preparations, a
receiver-controlled latch, a deadline reset and storage for the secret
sign during the block are supplied. The block communicates a receiver
decision. It cannot identify which single click came from a reflection.
[Reader, hidden signs and resource bounds](remote_message.py)

The 31 physical criteria require transfers to their actual preparations,
apparatus and observer records. The scalar parity code, an externally
timed detector and an ideal captured mode delimit this finite calculation.
There is no measured particle mass, speed or gravitational coupling in its
selection, and it assigns no particle identity.

## Reproduce

From the repository root with [the pinned dependencies](../../requirements.txt),
set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` and run

```sh
python -B -O -W error -m model.scalar_channel.contact
python -B -O -W error -m model.scalar_channel.contact.pair
python -B -O -W error -m model.scalar_channel.contact.material_radar
python -B -O -W error -m model.scalar_channel.contact.radar_checks
python -B -O -W error -m model.scalar_channel.contact.message_replay
python -B -m pytest -q -W error model/scalar_channel/contact
```

The [collar protocol](protocol.json) and [mobile protocol](pair_protocol.json)
specify the acquisition and falsifiers. The first replay regenerates every
count with the specified seed, checks the sampling errors and reads the
couplings from counts.
The second recomputes each entry of the [joint laws](pair_records.json).
The seed replay compares the material radar counts with the original
[acquisition](radar_counts.npz). The clock archive reader authenticates the
[530,208 integer sampling weights](clock_weights.npz) against the digest in the
[reviewed acquisition reader](https://github.com/FloatingPragma/observer-patch-holography/blob/44e45533078dfc592451f33507cbf3fca1ef2d11/model/scalar_channel/contact/material_radar.py#L25-L33).
It compares all six families with the normalized regenerated source and
scores the saved counts against both distributions. A row has $K$ bins on
a lattice of $L=10^{12}$ units. Rounding gives each ordinary bin an allowance
of $1/(2L)$. One largest rounded bin absorbs the remainder with allowance
$K/(2L)$, which also bounds the row total variation. The separate
floating-point allowances are $1/(2L)$ per bin and $1/L$ in total
variation. The digest commits to the integer weights exactly.

Run the archive reader from the repository root with the two thread settings
above. Add `--details` for the per-family discrepancies and held-out folds.

```sh
python -B -O -W error -m model.scalar_channel.contact.material_radar --read-archive
```

This command reproduces the archive on x86-64 Ubuntu 24.04.4 with CPython
3.12.3 and the dependencies in `model/requirements.txt`. On Windows 11
x86-64 with CPython 3.13.2, NumPy 2.5.1 and SciPy 1.18.0, the default
NumPy dispatch reproduces every integer weight. With
`NPY_DISABLE_CPU_FEATURES=X86_V3`, 141 re-rounded entries differ by one
or two lattice units. Both runs pass the archive checks. The other two
radar programs also accept `--read-archive`.

For the large message blocks the seed replay uses the original
[integer sampling distribution](message_weights.npz). It verifies the
distribution against the regenerated tagged source before drawing any
counts. A change of one grid unit can affect several trials in a block.
The tests compare the first experiment with full native memory and forward
scalar evolution, and the second with propagation on both explicit pointer bits.
The cell calculations have finite register and menu sizes. Neither instrument
assigns an observed particle, an external time unit or a continuum field.
