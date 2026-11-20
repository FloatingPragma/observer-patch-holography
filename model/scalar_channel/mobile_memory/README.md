# Mobile charges, bound loads, and a written cage link

The [speed and clock tests](relativity/README.md) check this generator and a separately labeled
incidence transport variant on the same carrier. The [calibrated incidence experiment](calibrated/README.md)
combines charge and cage transport, mass tags, a native swap, and an onsite contact.

An external seam can carry a record that later changes which scalar detector fires. The
mobile-memory candidate keeps the contrast-bound pair and the scalar cage on the same carrier and
couples them through the native contrast square. It also carries two signs of charge in a mobile
lowest band. A signed charge-current read on the seam can control the same cage link or steer a
later conversion of opposite charges into a neutral pair.

## Declared source

The carrier, native load alphabet and contrast memory come from the A-contrast-memory-spherical-v1
candidate, and the scalar seam and wedge modes come from the scalar-channel construction. At each
port the mobile-memory candidate also carries the two charged fermion modes and the neutral
oscillator of the radiative pair converter. The scalar cage excitation is a particle and hole
in the seam and wedge modes. The neutral oscillator quanta belong to the pair converter. At each
external seam there is an auxiliary orthogonal current record with values $-2,-1,0,+1,+2$, with
zero as the blank value.
[Port carrier](../../contrast_memory_spherical_v1/carrier.py) · [Native contrast source](../../contrast_memory_spherical_v1/law.py)
· [Scalar cages](../cage/construction.py) · [Pair converter](../radiative_repair/source.py)

Use exchange units with $J=1$. Each charge sign has

$$ Q = \mu I + L_{\mathrm{within}} + \eta L_{\mathrm{between}},\qquad \mu=1,\quad\eta=1/8. $$

The two Laplacians are sums of positive two-port difference squares. The term $\mu I$ is an
onsite occupation cost. The neutral oscillator and on-port converter retain the radiative generator.
The native contrast Hamiltonian uses $\kappa=128$. On each cell choose the first labeled degree-five
ring from the carrier enumeration. Write $X_c=f_{c,-}^{\dagger}f_{c,+}$ for its particle–hole lowering
operator in the two exact cage modes. Each external seam $s$ between cells $c,d$ contributes

$$
g(A_s+R_s^+)\otimes(X_c-X_d)^\dagger(X_c-X_d),\qquad g=1,
$$

where $A_s$ is the native echoed contrast square and $R_s^+$ is the projector onto a positive
auxiliary current record. The term is indexed by the seam and supported in its two incident cells,
on the selected rings and their attached wedges. The controls, carrier labels, selected rings,
orthogonal pointer alphabet, local pulse schedule and all coefficients are supplied. The charge,
binding and readout experiments use the reducing preparations of this generator with the indicated
inactive registers in their references. The filled scalar reference is annihilated by each $X_c$.
With one particle–hole excitation, each scalar link is precisely the graph Laplacian of that link. The generator therefore
preserves the whole cage during transfer. [Charge generator](law.py) · [Full Fock check](cage_transfer.py)

## Binding with the cage link enabled

For binary two-load configurations on blank native memory, the actual one-seam operator restricts
to $A_s=(n_a-n_b)^2$. The external seams form a port matching, so at most two links are active in
each configuration. On the one-cage subspace their combined positive contact has norm at most
$4g$. With the auxiliary current register blank the interaction is diagonal in the load addresses.
It leaves the off-diagonal block between adjacent and separated load pairs unchanged. Subtract the
common free cage energy. The adjacent block rises by at most $4g$, and the separated block can
only rise. The ordered-block binding theorem then gives, whenever $2\kappa-6-4g>0$,

$$
\Pr(\text{no adjacent load pair at a read})\le
\min\left\{1,\frac{40}{(2\kappa-6-4g)^2}\right\}.
$$

The bound applies at every free wait with a fixed enabled subset of external links, including
the full set. The replay also reads a joint load-pair and cage outcome from the 20-cell carrier,
including separated pairs. It uses one enabled link and the full external interaction with the
same preparation. With empty loads, the stationary blank and written native collar arms
give exact link weights zero and one. A sharp cage sent across one enabled link arrives intact with
probability $\sin^2(gt)$. The finite fuel read includes the no-sender background under its source-on
capture pulse. [Binding theorem](../../contrast_memory_spherical_v1/PROOF.md#3-persistence-at-every-value-of-the-evolution-parameter)
· [Joint record replay](contrast_binding.py)

## Motion in the lowest charged band

The cell-uniform space has one mode per cell. The second energy of one complete cell is
$\mu+\Delta$, with $\Delta=5-\sqrt{5}$. The external matching has operator norm at most two.
Min–max confines the first $C$ levels of a $C$-cell federation to $[\mu,\mu+2\eta]$, with the
remaining levels at or above $\mu+\Delta$. The uniform state on a connected federation is the
unique lowest eigenstate. The rest of its lowest band has strictly greater energy. On cell-uniform
states the leading intercell compression is $\eta L_{\mathrm{cell}}/12$.

Prepare the uniform orbital in cell zero and apply the declared spectral filter for this band. The
acceptance is at least $1-(2\eta/\Delta)^2$ by projecting $(Q-\mu I)$ onto the higher band. Read
every cell at the waits 0, 8, 24, 48 and 96 and retain every rejected filter trial. The parked control
uses the same filtered preparation and external diagonal costs. It removes only the intercell
hops, so each cell's total occupation is constant. [Mobility replay](mobility.py)

## A current that changes a later record

Orient the selected seam $(i,j)$ from the sender towards the receiver. The signed-current instrument
diagonalizes $iq[Q_{ij},n_j]$ for charge $q$. Dividing by the positive seam coefficient gives
outcomes from $-1$ to $+1$ for a single charge and $-2$ to $+2$ for the pair. All equal-eigenvalue
subspaces are kept coherent, including the neutral part of the zero-current sector. The write is a
conditional translation on a cyclic orthogonal pointer.

For the one-charge test, use the same accepted low-band packet and free wait 24, with either charge
sign. A positive written value enables the neighboring cage link. The shuffled arm draws an
independent value with the same frequency within each charge sign. Both the current and receiver
fuel marginals agree between those arms, while the full joint record distinguishes them. The
replay gives the complete current × two-bit fuel records, the filter failure, and the enabled
and blocked no-sender controls.

For the two-charge test, balanced preparations have opposite imaginary seam phases and equal
position probabilities. After the current write an addressed seam pulse reads the stored sign to
route both charges to the target port. For seam weight $w$ the two pulse arms use the difference
square for times $\pi/(4w)$ and $3\pi/(4w)$, with other terms parked during the pulse. Release under
the same charge and conversion generator for $t=1/128$ in inverse exchange units, then read all
charged and neutral outcomes. Keeping the pulse arm fixed or shuffling its control gives the matched comparison. If $\lVert H\rVert\le M$ and
$\rho=(Mt)^3/6$, the target conversion amplitude has imaginary part $t+O(t^3)$ with remainder at
most $\rho$. Conversion from the other endpoint into the target read has zero coefficients through
order two and modulus at most $\rho$. The balanced feedback advantage is at least
$((t-\rho)^2-\rho^2)/2$ for this wait. [Current instrument](current.py)

With the charge parked during the conditional cage pulse, the address projectors commute with the
controlled action. The same writer can then run backwards and return its private pointer to blank.
An optional orthogonal copy retains a public record. A finite operator check compares every input
column with both write–control–unwrite maps. [Record retention check](memory.py)

## Co-active charges, loads, and a cage

Read the signed current of a filtered charge on the addressed seam. In the same trial, prepare
adjacent native loads and a cage on the sender, then let all external links act during the common
free interval. A positive current record adds one copy of the addressed link, with the native
contrast weights acting throughout. For a fixed record the charge and the load–cage state evolve
in disjoint factors, while the charge-current pointer controls the load–cage generator. The
terminal read retains the actual current, the final charge cell, load-pair presence, and the cage
cell. A shuffled control draws the enabled-link bit independently with the same frequency within
each charge sign. The read keeps the rejected filter trials in the denominator. The lowest-band
preparation applies before the current read, which can populate other charged bands.

For each configuration the native links have norm at most $4g$, and the extra addressed link
has norm at most $2g$. The interaction stays diagonal in the load addresses during the free
interval, so the same ordered-block argument bounds the probability of separated loads by
$40/(2\kappa-6-6g)^2$ when the gap is positive. At wait 0.5 the full records distinguish the
actual and shuffled writes by total variation 0.017444 and 0.047289 for the two charge signs.
The separated-load probability in the positive-pointer arm is $5.51874\times10^{-5}$, within
the all-time bound $6.71863\times10^{-4}$. A direct 48-state dense check compares the
charge–load–cage operator and all joint read bins. [Co-active replay](coactive.py)

An independent short-wait certificate uses the joint read of the *original* adjacent pair and
the addressed receiving cell. For this one load configuration, at most ten native swaps have
absolute matrix element $1/2$. The squared coupling out of the configuration is at most $5/2$,
so two applications of Duhamel's formula bound the return-block error by $5t^2/4$.
The cage generator compressed to that configuration has norm at most $6g$. Its receiver
amplitude starts at $igt$ for the blank pointer and $2igt$ with the positive pointer, with
remainder bounded in both arms by $\epsilon=(18g^2+5/4)t^2$. At $g=1,t=1/128$,
$(2t-\epsilon)^2-(t+\epsilon)^2$ is strictly positive. The nine allowed swaps for the
addressed pair sharpen the probability increase bound to $1077/8388608$. The replay checks
the slope in the assembled generator and reads probabilities $6.10278\times10^{-5}$ and
$2.44095\times10^{-4}$ with all links enabled.

## The addressed link on port registers

For the addressed external seam, embed the two selected cage pairs in a two-cell register
chart. Local mode rotations along chart links send the two cage modes of each cell to its
incident port's fuel and reference modes. The chart checks those profiles against the
full federation and checks the rotation and central words against its cycle constraints.
Conjugating the central link pulse by the routing rotations implements the same quartic on
the original modes.

The central generator expands into two Pauli batches on the two endpoint observers'
registers. One contains the two local occupation costs. The other exchanges the cage
excitation. The words in each batch commute pairwise, and the two batch sums commute
by an integer Pauli-polynomial check. The rotations therefore implement the central
pulse exactly, including the identity term's phase when the pulse is controlled.
A separate check compares the generator, occupation operators and pulse against the
literal even fermion sector of the bridge and its parity anchor. The chart and the
routing sequence are supplied for an isolated addressed pulse. For seam $(10,20)$ the
two batches have seven and eight words, and the full pulse agrees on all 32 literal
even-sector columns to within $10^{-14}$.
[Register compilation](register_link.py)

## Refinement of the cage links

In the stored 20-to-80 and 80-to-320 carriers, the parent of fine cell $j$ is $\lfloor j/4\rfloor$.
Each parent contains a central child joined to three corners. Each coarse seam has two fine
crossings, and the three sets of corners incident to the three neighboring parents omit different
corners. With all links written in the one-cage sector, remove the common free cage energy and
write $L_c,L_f$ for the coarse and fine cell Laplacians in link units.

For the isometry $P$ sending a coarse cell to the uniform superposition of its four children,
each fine crossing contributes one quarter of its corresponding coarse difference square.
It follows that $P^*L_fP=L_c/2$. The leakage $E=(I-PP^*)L_fP$ satisfies $E^*E=L_c/2$ as well.
To see the latter, let $d_1,d_2,d_3$ be the coarse amplitude differences from a parent to its
neighbors and $s=d_1+d_2+d_3$. The leakage on the center is $-s/4$ and on the corner missing
neighbor $k$ it is $(s-2d_k)/4$. Their squared moduli sum to $\sum_k|d_k|^2/4$ in that parent.
Each coarse seam is counted twice over all parents.

Composing the two checked isometries gives $P_{02}^*L_{320}P_{02}=L_{20}/4$.
For the composed leakage, the term orthogonal to the 80-cell lift and the lifted first-step
leakage are orthogonal. Their Gram matrices contribute $L_{20}/4$ and $L_{20}/8$, respectively.
The total unscaled leakage Gram matrix is $3L_{20}/8$.

There is a stronger local obstruction. Suppose a subspace is a direct sum of subspaces inside
individual parents and is invariant under $L_f$. Each off-diagonal block between neighboring
parents and its adjoint preserve the respective subspaces. Their product projects onto the
two crossing corners. Two such corner projections within a parent isolate one corner.
The internal star then takes any such corner into the span of the center and takes the center
into the span of all three corners. Invariance therefore forces a nonzero parent subspace to
contain all four child coordinates. A crossing propagates this condition to the adjacent parent.
Connectivity leaves only the zero and full subspaces. The same argument uses support projections
of the block products for nonzero positive seam weights. A proper exact parent-local
intertwiner for this one-cage generator is therefore excluded.

The finite detector prepares a cage uniformly within one parent and reads the receiving parent
across a specified seam. The fine detector has a short-time leading probability $t^2/2$ from two
distinct receiving children. The energy-matched coarse walk generated by $L_c/2$ gives
$t^2/4$. The replay checks the full distributions and the zero-link control on the actual
carriers. At $t=1/32$ the matched receiver-probability difference is about $2.43327\times10^{-4}$
on both steps, with the rational lower bound $17598039/107374182400$ from the sine and cosine
remainders. [Refinement certificate](refinement.py)

## Records on 20 cells

The charged, conversion and cage runs share the addressed external seam $(10,20)$. The probabilities
include the failed filter or complementary terminal outcomes in each trial denominator.

| Preparation | Observer read | Probability |
| --- | --- | ---: |
| Filtered lowest-band charge | Cell 1 at wait 0 | 0.000060 |
| Same preparation | Cell 1 at wait 48 | 0.157988 |
| Adjacent loads with a sharp cage and all links active | Both adjacent loads and a cage outside the sender at wait 0.5 | 0.218504 |
| Same preparation with zero contact | Both adjacent loads and a cage outside the sender | 0 |
| Balanced charge pair and recorded feedback | Neutral pair at the target at wait 1/128 | $6.10065\times10^{-5}$ |
| Same preparations and a fixed or shuffled pulse arm | Neutral pair at the target at wait 1/128 | $3.05033\times10^{-5}$ |

In the full coupled load-and-cage run the probability of separated loads is
$5.49024\times10^{-5}$ at the read, with the all-time bound $6.60983\times10^{-4}$.
A written link gives a receiver fuel-read probability of 0.703690 for a sharp sender and
0.296021 for the no-sender reference. In the one-charge experiment, the joint records with stored
feedback and shuffled records have total variation distances 0.015125 and 0.041004 for the negative
and positive charge respectively.

## Scope of settling

For a fixed proper target subspace with complement projector $Q_{\rm out}$ and any finite unitary
$U_t$, $\lVert Q_{\rm out}U_t\rVert=1$. Uniform settling of arbitrary inputs into that subspace
requires a specified open channel or a restriction on inputs. The unwrite above starts from the
declared blank-pointer frame and preserves unitarity on all inputs.

## Reproduction

The replay and checks use the dependencies in `model/requirements.txt`.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.scalar_channel.mobile_memory
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m pytest -q model/scalar_channel/mobile_memory
```

The quantities here are finite detector records under supplied preparations and controls.
Autonomous scheduling of the current reads, the link controls and the concurrent source terms
requires a schedule on the native registers. Relativistic dispersion, proper-time normalization,
coupling selection and a dynamical physical refinement require additional attachments to the
carrier and observer instruments.
