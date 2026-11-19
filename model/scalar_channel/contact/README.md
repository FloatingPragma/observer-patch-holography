# A local interaction between contrast memory and scalar occupation

A collar changes a scalar fuel reading. Filling the scalar source changes the
phase read from the collar. A mobile load pair and a scalar quantum also give
correlated readings at two local pointers. All three experiments use the same
Hamiltonian on one cell. The largest distance between the mobile joint record
law and the product of its marginals is 0.03634 for positive contact, 0.03127
for negative contact and below $2\times10^{-14}$ for zero contact.

## Law and source map

The experimental member `contrast-square-scalar-contact-v1` takes the tensor
product of the two register families in the [construction manual](../../construction_manual/README.md).
Each port carries its load, each seam its eleven-state memory, and the scalar
even algebra uses the 432 link registers. In the manual's port labels set

$$H(g)=H_{\rm contrast}(8)\otimes I+I\otimes d\Gamma(D)
       +g\sum_s A_s\otimes n_s .$$

Here $A_s=U_s(0\oplus D_{\rm memory}^2)U_s^*$
is the coefficient of $\kappa$ in the complete
[99-state seam term](../../polynomial_sv/law.py), with identities on loads in
the memory generator. $n_s=(I-B_s)/2$ is the [encoded seam occupation](../compiler.py).
Its parity word commutes with the scalar loops. Each product uses the same
seam on both sides. Reversing the seam exchanges its ports and sends
$(c,d)$ to $(1-c,-d)$, preserving the echoed square. The
[port chart](../../contrast_replay.py) maps the two carrier labellings.
The degree compensation of the contrast source is unaltered. Energy is in
exchange units, with scalar rest level 1 and hop $1/\sqrt2$. That relative
normalization, the contact operator, the apparatus and the signs tested
($g=0,+1,-1$) are design choices.

On scalar vacuum the full contrast source is exact. On empty loads and
blank memories every $A_s$ vanishes and the full scalar source is exact.
These embeddings also hold with untouched quantum references. In particular
the one-load cell transport carries the port-Laplacian fingerprints
$5-\sqrt5$ and 6 after the native factor of one half and the constant
reference shift. The free scalar cage has line $2\sqrt6$. Every value of
$g$ has the same restrictions to these two sectors. Even their complete
record laws, for instruments staying in those sectors, cannot determine $g$.
The $1/6:1/12:1/12$, $\sqrt2$, $2/9$ and $3/8$ targets have no matched
preparation and reader in this experiment. The [centre class and plaquette
flux costs](../../defect_particle_classes/README.md) use a five-level code
and a protected link flux. Neither is a register in this contact. These
targets receive no score. The experiments use no measured masses,
fine-structure or gravitational constant. A physical residual requires
a unit map and a particle or field interpretation of a public record.

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
are 0.64345 and 0.64862. Subtracting the undriven reading from the sharp
preparation at matched waits and analyzers gives changes up to 0.60490 and
0.48675. Each setting uses one million simulated trials. These are readings
of fuel pointers. The calculator's eigenvectors and coordinates are private.

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

On blank memory the [contrast identity](../../polynomial_sv/PROOF.md)
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

There is also an all-time bound before the terminal read. In the binary
two-load sector, on a single scalar quantum, the contact is diagonal with
spectrum between $\min(0,g)$ and $\max(0,g)$. Doubled loads are invariant
and the emitter cannot reach them. The scalar width is $2\sqrt6$ by the
incidence identity. Thus the native cell's adjacent and separated diagonal
blocks retain ordered gap at least
$\gamma=2\kappa-1-3\sqrt6-|g|$, and their off-diagonal norm is 2.
[The block separation theorem](../../polynomial_sv/PROOF.md) bounds
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

## Reproduce

From the repository root with [the pinned dependencies](../../requirements.txt),
set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` and run

```sh
python -B -O -W error -m model.scalar_channel.contact
python -B -O -W error -m model.scalar_channel.contact.pair
python -B -m pytest -q -W error model/scalar_channel/contact
```

The [collar protocol](protocol.json) and [mobile protocol](pair_protocol.json)
specify the acquisition and falsifiers. The first replay checks the complete
count receipt against regenerated laws and reads the couplings from counts.
The second recomputes each entry of the [joint laws](pair_records.json).
The tests compare the first experiment with full native memory and forward
scalar evolution, and the second with propagation on both explicit pointer bits.
The cell calculations have finite register and menu sizes. Neither instrument
assigns an observed particle, an external time unit or a continuum field.
