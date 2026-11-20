# Proofs for the seam program

## Realization by A-contrast-memory-spherical-v1 seam events

Mark port $m$ of a connected carrier and restrict the controls to seams avoiding $m$. Encode the moving slot
at port $j$ as the two-load state $|\{m,j\}\rangle$. On this subspace a seam has at most one of the two loads.
Its event acts by $\exp[-i\tau(I-T_{ab})/2]$, where $T_{ab}$ interchanges the endpoint basis states. All
other ports are spectators. The subspace is reducing. For $\tau=\pi$ this event exchanges two basis states
exactly. Conjugating the final edge of a path of length $\ell$ by the first $\ell-1$ swaps implements the
desired exchange between its endpoints in $2\ell-1$ events. Every intermediate port is restored, even if
it carried an unknown coherent amplitude. This construction needs paths avoiding the marker.
[Seam event](../contrast_memory_spherical_v1/SPEC.md#8-events-and-schedules) · [Compiler](compiler.py)

On one eight-port block label the slot basis by the three-bit words. For each nonzero mask $a$, the four
disjoint pairs $i,i\mathbin{\mathrm{xor}}a$ generate $(I-X_a)/2$, where $X_a$ translates by that mask.
All seven translations commute, and $\sum_{a\ne0}X_a=8P-I$ with $P$ the equal-slot projector. Acting
on each pair once with angle $\pi/4$ gives

$$
 \prod_{a\ne0}\exp[-i(\pi/8)(I-X_a)] = \exp[-i\pi(I-P)] = 2P-I.
$$

Applying this program on each block and then exchanging each positive outgoing slot with the negative
slot at its adjacent logical address implements $SC$ exactly. Arbitrary slot states and an untouched
quantum reference are covered. This is an identity at completed program boundaries. The paths, clock
counts and incidence of the physical seams between those boundaries belong to the chosen embedding.

There are two native expressions to check. On blank memory the original 99-state law reduces to
$(I-T_{ab})/2+\kappa(n_a-n_b)^2$. On binary loads the difference from the scheduled generator of
section 8 of the specification is $\kappa(n_a+n_b)$, which commutes with both generators. At the
declared $\kappa=8$, both $\kappa\pi/4$ and $\kappa\pi$ are multiples of $2\pi$. The two event unitaries
therefore agree, including relative phases between occupied and empty seams. This argument holds on
the entire binary pair sector, including preparation error. Blank seam memories return to blank.
Other couplings and pulse lengths require their own phase check. [Blank-memory identity](../contrast_memory_spherical_v1/law.py)

The finite emitter starts with fuel on the marker and source seam. Its error relative to the ideal
pair is bounded in norm by $\sqrt{d_m+d_j-2}/(2g)$. Both ideal and actual states then undergo the
same seam events. The error bound carries to the total variation of any terminal observer distribution
by unitary invariance and contraction of trace distance under measurement. For a selected slot, swap
its port onto the emitter seam, run the two endpoint observers and read their public decision. The
complete seam measurement also reports every other seam and the no-pair outcome. Residual fuel
contributes to no pair. A finite-duration read under the hypotheses of the native reader adds at most
$12\tau_{\rm read}$ to this bound. The execution stores ideal instantaneous readouts after the finite
emitter. Distinct selected-slot settings use fresh preparations.
[Emitter](../contrast_memory_spherical_v1/sectors.py) · [Reader](../contrast_memory_spherical_v1/reader.py)

**Preparation obstruction for the routing pulses.** On a connected control component of $N$ nonmarker
ports, each section 8 partial exchange fixes the uniform vector. So does an unshifted pulse with
$\kappa\tau\in2\pi\mathbb Z$, which includes all routing pulses above. Starting with a sharp marked pair, its overlap with that
vector is $1/\sqrt N$ under every word of these gates. Any target orthogonal to the uniform vector
is therefore separated from the reachable pure states by at least $1/\sqrt N$ in trace distance.
This follows by projecting the reachable state onto the orthogonal complement of the uniform vector
and bounding its overlap with the target by $\sqrt{1-1/N}$. Untagged nonzero-momentum packets
have zero total amplitude. Embedding them with zero amplitude at unused ports
preserves that property. Even the finite emitter retains a lower bound of $1/\sqrt N$ minus its
preparation error, by the triangle inequality.

The unshifted law also admits nonresonant durations. On the single cell, hold the marker at port 0
and start the mobile load at port 1. A pulse of duration $\pi/2$ on $(1,5)$ gives amplitudes
$c=(1-i)/2$ and $s=(1+i)/2$. A pulse of duration $3\pi/16$ on $(4,5)$ multiplies the sum on that
edge by $e^{-8i(3\pi/16)}=i$. The total amplitude is then $c+is=0$, with unit norm and no intervening
measurement. Thus the uniform-mode bound concerns the specified routing pulses. The two-pulse example
prepares one zero-sum state. The following protocol prepares a specified arbitrary packet.

## Preparing finite wave states from native fuel

For $M\ge2$ slot amplitudes provide $2M$ distinct nonmarker ports in a connected component. Use the
antisymmetric rail basis $|j_L\rangle=(|j_+\rangle-|j_-\rangle)/\sqrt2$. A four-port equal-slot reflection
uses three XOR matchings of angle $\pi/2$. On an initially sharp pair it gives amplitudes
$(-1,1,1,1)/2$. Query the final two ports for the marked pair, separately, bringing each to the native
seam reader and restoring the queried port. A yes ends the trial. Two agreed no records leave the first
logical basis vector with probability $1/2$ and a known overall sign. The probabilities of the two
early exits are $1/4$ each. Each absence measurement preserves coherence between all sites giving no.

For arbitrary angles retain the unshifted blank-memory contrast law, including its commuting endpoint
occupation phase. Apply one exchange between the positive ports of two rails and the same exchange
between their negative ports. On this pair of logical amplitudes the resulting matrix has diagonal
$e^{-i(\kappa+1/2)t}\cos(t/2)$ and off-diagonal $i e^{-i(\kappa+1/2)t}\sin(t/2)$. It acts identically
on all other logical amplitudes. An exchange within a rail acts on its antisymmetric basis state by
$e^{-i(\kappa+1)t}$, giving an arbitrary relative phase by a nonnegative duration. Intermediate path swaps
at angle $\pi$ are exact at $\kappa=8$.

Given a normalized complex target $\psi$, split the initial amplitude successively into labels
$j=1,\ldots,M-1$. At split $j$ choose $t/2$ as the angle with opposite side $|\psi_j|$ and adjacent side
$\sqrt{|\psi_0|^2+\sum_{k>j}|\psi_k|^2}$. If both sides vanish the donor is empty and may be left alone.
This gives the prescribed magnitudes. Track the common phase of each split and the extra factor $i$ on
the newly occupied rail. A final within-rail phase makes each amplitude equal to $\psi_j$ up to one global
phase. Pulse angles depend on the declared preparation, including its phases. They use no moving-clock
record, measured mass or other physical comparison.

Exclude the $M$ negative ports in a fixed order with the same seam reader, restoring the routes. Each
rejection has unconditional probability $|\psi_j|^2/4$. Following all no records, the **original single-rail
encoding** holds $\psi$ with unconditional probability $1/4$. Each rejection has its own agreed terminal
record. This proves that every specified finite packet admits a
heralded implementation on a carrier of sufficient size. The acquisition executes one complex four-mode
packet on 72 slots. For $M$ outputs the source uses at most $6+2(M-1)+M$ nonlocal exchanges before
routing, and $M+2$ binary queries on its longest branch. Each exchange on a path of length $\ell$ uses
$2\ell-1$ seam events. An aperture query swaps its port to the reader and restores the path after the read.
The [receipt](preparation_records/result.json) charges those edges and all four reader stages.

The finite fuel emitter replaces the ideal sharp pair throughout this complete decision tree. Contractivity
bounds the total variation on all stops and the final seam records by the initial emitter norm error.
On this binary sector the unshifted single-edge generator has norm at most $\kappa+1=9$. If pulse $j$
misses its duration by $\delta_j$, its operator error is at most $9|\delta_j|$. Telescoping and retaining
all classical branches adds at most $9\sum_j|\delta_j|$ on the longest path with these error budgets.
For a binary read followed by further coherent operations, include the diagonal phase acquired during
the read. In the fixed-number energy convention the free law is $H=K-2\kappa B$, where
$0\le B\le I$ is the adjacent-pair projector and $0\le K\le12 I$ by the maximal pair-graph degree.
At $\kappa=8$ the spectrum lies in $[-16,12]$. The decoupled fuel's energy is a diagonal entry in
this interval. Comparing the entire driven reader with the instantaneous reader by Duhamel therefore
costs at most $16\tau$ in joint state error for a read lasting $\tau$. The sharper $12\tau$ terminal
record bound can discard diagonal phases at the final measurement. Summing the joint bounds budgets the distribution
before conditioning on successful preparation. This budget treats the scheduler, finite fresh pointer
words and conditional stop/continue service as apparatus. The native clock time of that program and
a physical identification of the traveling state require separate observations.

## Reading phase agreement at the home seam

Let $K_\psi$ be the ordered product on the accepted branch of the ideal preparation. It contains the
unitary event words and the no projectors of the binary seam readers. In the marked-pair code,
$K_\psi|a\rangle=-|\psi\rangle/2$ up to an overall phase, where $|a\rangle$ is the input pair on
the fuel seam. Execute the instructions in reverse order, conjugating every unitary and repeating
each no projector with a fresh nominal reader. The surviving branch then applies $K_\psi^\dagger$.
The binary unshifted generators at $\kappa=8$ have eigenvalues $0,8,9$. Running an angle $-t$ modulo
$2\pi$ implements the adjoint with positive-duration pulses, including their phases on states outside
the marked code. Each reversed no projection has the same accepting operator as its forward read.
A yes stops that trial and writes its own outcome; the read is a fresh measurement.

For an emitted ideal pair with programmed source $\phi$ and receiver $\psi$, the unconditional
probability of returning to the home seam is consequently

$$
 |\langle a|K_\psi^\dagger K_\phi|a\rangle|^2=\frac{|\langle\psi|\phi\rangle|^2}{16}.
$$

All input and analyzer phase choices use the same seam predicate. Identical omitted phases therefore
restore the large match probability. The finite experiment includes unconverted fuel, disagreements
between the programs, every early termination at either end and the final complete seam measurement.
The finite-emitter norm bound applies to this entire instrument, independent of the number of exact
unitary events. Timing and finite read-window errors accumulate according to the separate budget above.
[Phase-sensitive receiver](receiver.py) · [Complete readout](receiver_records/result.json)
