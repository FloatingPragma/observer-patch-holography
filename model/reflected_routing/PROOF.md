# When does the scalar cone limit a message?

## Free propagation

Bernhard Mueller's eight-slot update is $U=SC$. The local coin is $C=2P-I$ with $P=|s\rangle\langle s|$ and
$s=(1,\ldots,1)/\sqrt8$. In a Fourier fiber the flip-flop shift $S(k)$ moves each outgoing amplitude by one
signed step and exchanges the two slots on its axis. Thus $S^2=I$. Three coordinates serve as spatial dials
and the fourth is private. The proof uses lattice units on the infinite lift. A periodic measurement uses it
only before a periodic image can reach the detector. [Source gate](https://github.com/muellerberndt/openwave/blob/260de06352eac5425d78fbbae636625544d439ab/simulation/experiments/observer_geometry/reunion/registered/ingredient_sources/radar_shell/replay.py)

For $v=Ss$ and $a=\langle s,Ss\rangle=\tfrac14\sum_j\cos k_j$, one has $Us=v$ and $Uv=2av-s$.
The dispersive frequencies obey $\cos E=a$. Their group velocity has components $\sin k_j/(4\sin E)$.
Since $\sum_j\cos^2 k_j\ge4a^2$, its norm is at most $1/2$, approached by long waves. The orthogonal
complement of $s,Ss$ evolves by $-S$ with frequencies zero or $\pi$. Limits cover coincident bands.
A bound on group velocity alone does not control finite-time communication.

## All-slot exponential bound

Let a source occupy $\hat n\cdot x\le0$ and a detector occupy $\hat n\cdot x\ge d$, for a unit vector $\hat n$.
Weight each position amplitude by $\exp(\theta\cdot x)$, where $\theta=\eta\hat n$ and $\eta>0$.
The shifted Fourier operator contains $S(k+i\theta)$, whose norm is $D=\exp(\max_j|\theta_j|)$ and whose square is the identity. Define

$$
 A(\theta)=\operatorname{arcosh}\left(\tfrac14\sum_j\cosh\theta_j\right).
$$

The complex $a$ lies in the ellipse with semiaxes $\cosh A,\sinh A$. In the real planar direction of angle
$\phi$ the support function of the $j$-th cosine ellipse is $\sqrt{\cosh^2\theta_j-\sin^2\phi}$.
The function $x\mapsto\sqrt{x^2-\sin^2\phi}$ is concave on $x\ge1$. Averaging proves the inclusion.
The two roots of $z+z^{-1}=2a$ therefore have modulus between $e^{-A}$ and $e^A$.

The row space spanned by $\langle s|$ and $\langle s|S$ is invariant. With second-kind Chebyshev polynomials
$V_j$, the row $g_j=\langle s|U^j$ equals $V_j(a)\langle s|-V_{j-1}(a)\langle s|S$, where $V_{-1}=0$.
The Laurent sum gives $|V_j(a)|\le(j+1)e^{jA}$, including repeated roots. Expand $U^t$ about $(-S)^t$
using $U+S=2SP$. Its $j$-th additional term is $2(-S)^{t-1-j}S|s\rangle g_j$. The factor left of
the ket has norm at most $2D$. Plancherel and removal of the position weight give

$$
 \|\Pi_{\rm det} U^t\Pi_{\rm src}\|
 \le D[1+(1+D)t(t+1)]\exp(t A(\theta)-\eta d).
$$

This covers every initial slot polarization, both flat bands and an untouched quantum reference.
Since $A(\eta\hat n)=\eta/2+O(\eta^3)$, fixed $d/t>1/2$ gives exponential suppression. For $B$ noninteracting
signal particles initially confined to the source half-space, the expected receiver occupation is at most $B$
times the square of this bound, by applying the estimate to each one-particle marginal. The probability of
detecting any particle is no larger. A polynomial number of independent occupied-rail trials cannot overcome
this suppression. Pre-existing remote carriers, intermediate writes and active receivers that transmit during
the flight require their own analysis.

## Sharp cost of a faster passive signal

For an equal-slot state localized at the origin, the squared norm of the weighted wave is
$M_t(\theta)=\mathbb E[\exp(2\theta\cdot X_t)]$ under its complete position measurement.
The norm bound before removal of the weight gives $\limsup t^{-1}\log M_t\le2A(\theta)$.
At $k=0$ in the shifted Fourier integral the equal-slot return is $T_t(\cosh A)=\cosh(tA)$.
For $\theta\ne0$, in an open neighborhood the larger root has modulus exceeding $\exp(A-\epsilon)>1$.
Its reciprocal cannot cancel it in $T_t=(z^t+z^{-t})/2$ at large $t$. Integrating the squared return over this
fixed positive volume proves the reverse limit. This argument also works for the equal-slot detector alone
and for the axial, transverse-uniform preparation, using its one-dimensional Fourier integral.

On an axis $A(\eta)=2\operatorname{asinh}[\sinh(\eta/2)/2]$. For $1/2<u<1$ take the unique $\eta>0$
with $A'(\eta)=u$. Under the measure tilted by $\exp(2\eta X_t)/M_t(\eta)$, $X_t/t$ concentrates on $u$.
Apply exponential Markov bounds on either side using $M_t(\eta\pm h)/M_t(\eta)$.
Differentiability makes both bounds decay exponentially for sufficiently small $h>0$. Removing the tilt on
$|X_t/t-u|<\delta$ gives the lower bound for that window. For the half-line take a window just above $u$,
then shrink it towards $u$. The upper bound is the original exponential Markov inequality. The same proof
normalizes the tilted equal-slot subprobability by its own moment. Consequently

$$
 -\lim_{t\to\infty}t^{-1}\log\Pr[X_t\ge ut]
 =2(\eta u-A)=4[u\operatorname{asinh} b-\operatorname{asinh}(b/2)],
 \qquad b=\sqrt{(4u^2-1)/(1-u^2)} .
$$

At $u=1$ only the straight path arrives. The equal-slot launch has probability $(1/8)16^{-(t-1)}$ there
over all slots. A directed launch $C|{+x}\rangle$ has probability $16^{-(t-1)}$. Both give rate $\log16$.
Beyond that front the probability is zero. Within the cone the best uniform exponential bound is trivial.
For the transverse-uniform equal-slot experiment put $P_t=2\,4^t T_t((6+z+z^{-1})/8)$.
Its integer Laurent coefficients follow $P_{t+1}=(6+z+z^{-1})P_t-16P_{t-1}$ from $P_0=2$ and
$P_1=6+z+z^{-1}$. Squared coefficients divided by $4\,16^t$ give the equal-slot position law on every site.
Other slot outcomes carry its missing normalization. [Exact finite tail numerators](records/result.json)

## A message with local regeneration

On an axial segment of length $r$, a directed launch reaches the first allowed arrival with probability
$p=16^{-(r-1)}$. Each forward step after the first uses the off-diagonal coin entry $1/4$. On an oblique
free route all shortest paths ending in the same slot add coherently. The forecast counts permutations by
first and last step. Acquisition scatters the actual eight amplitudes at each reachable site. Pruning outside
the monotone box is exact at the shortest arrival. Each receiver reads both rails of each bit.

With $B$ copies per bit the segment delivers a bit with probability $1-(1-p)^B$. For a $b$-bit frame
and $m$ equal segments each way, the ideal authenticated return probability is $[1-(1-p)^B]^{2mb}$.
With one service tick per regeneration or remote reply the round trip takes $2mr+2m-1$ ticks. Initial
preparation and terminal read are boundary services. A sufficient budget for failure at most $\epsilon$
is $B\ge\log(2mb/\epsilon)/[-\log(1-p)]$ for $p<1$. For fixed $r>1$ the required parallel copies grow
logarithmically with distance. The asymptotic axial speed is $r/(r+1)>1/2$. Taking $r$ of order
$\log_{16}d$ and using order $d\log d$ copies per bit per segment approaches speed one with polynomial
total resources. This statement charges the declared parallel service. The circuit depth of reading and
regenerating an arbitrarily large bank in serial hardware is outside it.

The stations have reserved, initially empty lanes for each trial, hop, bit and direction. A source uses
fresh preparations for each occupied rail. A relay prepares its next frame only from a complete local read.
Other branches and failed preparations need no coherent reset and are included in the counts. The passive
comparison uses the same parallel budget per launch and no relay. The addressed comparison carries the
payload explicitly by the original local courier. The wave runs use the courier only to calibrate the
chart and deliver terminal reports. The message between wave stations travels through detections.

The sender chooses a fresh nonce after the forecasts are sealed. Distinct request and acknowledgment tags
bind the paired endpoint key, trial identifier and payload. The home station checks the returned nonce.
Keys and sender entropy are disclosed for offline replay. This is a trusted laboratory transaction with
explicit setup. The scalar deadline uses a flat Euclidean interpretation of the three reciprocal chart
ranges and the independently measured wave speed.

## Two pulse clocks and their finite resources

The original clock combines adjacent private Fourier orders. Here each prepared mode has private order zero.
An independently enrolled pulse circulates in a *public transverse direction*. At real momentum let
$w=(Ss-as)/\sin E$. In the basis $s,w$ the update is a planar rotation. The preparation
$(s+iw)/\sqrt2$ has eigenvalue $e^{-iE}$. The two-dimensional packet superposes five Fourier orders on
each axis with the registered binomial weights. It is uniform in the other two coordinates. After one
native rotation per service, local counters measure all position bins, summing over slots. Combining
those reports gives the first circular position moment. Its phase measures the pulse's displacement.

Separate preparations supply the reference speed $c$, longitudinal speed $v_x$, and transverse rest
lap rate. They predict the moving-to-rest ratio $\sqrt{1-(v_x/c)^2}$. For an infinitesimally narrow,
long-wave pulse this follows from the magnitude $1/2$ of the continuum group velocity. A finite-width
grating need not attain that limit by increasing the grid size alone. The complete finite-gate forecast
uses the actual prepared modes and dispersion. Reported passes use the registered tolerances.

For $M$ equally spaced counter thresholds the cumulative probability error at any cut is at most $1/(2M)$.
Summation by parts around the circle bounds the error in its complex moment by $\pi/M$. The angular
uncertainty at a nonzero measured moment $z$ is at most $\arcsin[\pi/(M|z|)]$. The moment couples adjacent
orders, whose frequencies differ by at most $\pi/L$ by the group-velocity bound. The sum of their absolute
overlap weights is at most one by Cauchy–Schwarz. Between enrolled samples spaced by $L/8$, its path stays
within $\pi/8$ of its starting value. At the registered visibility threshold 0.45 this excludes an
unobserved full turn. A uniform transverse translation erases the clock hand while preserving its longitudinal marginal.

The local challenge adds three reversible tag bits. On crossing the transverse seam at the edge of the
launch-centered chart, the tag changes by the signed crossing. Inverse Fourier transformation of the
prepared wave, injection into tag zero and a Fourier transform on the eightfold lift preserve its norm.
At $y=0$ the eight lifted addresses separate the tag outcomes. An aperture fixes $x,y,z$ and the private
address and accepts every slot. Both uniform spectator coordinates contribute the probability factor
$1/L^2$. The experiment uses $2^{16}L^4$ midpoint preparations per aperture and wait. It includes every
miss. For waits below $6L$ an initially zero-tag amplitude ending at $y=0$ and tag one has precisely one
positive net cut crossing. A path with winding $-7$ or $9$ would exceed the nearest-neighbor support.
The initial preparation has finite width. A tagged click identifies a winding of that amplitude and
a sampled return peak, without defining a point worldline.

The sealed local forecast propagates the two dispersive components by a rotation and the orthogonal
complement by $(-S)^t$, including both flat bands excited by writing the tag. Acquisition independently
iterates the full eight-slot gate. Suppressing the write leaves the different windings coherent in
the original torus. Summing measured tags after the write instead sums their probabilities.

## Realization by polynomial-sv seam events

Mark port $m$ of a connected carrier and restrict the controls to seams avoiding $m$. Encode the moving slot
at port $j$ as the two-load state $|\{m,j\}\rangle$. On this subspace a seam has at most one of the two loads.
Its event acts by $\exp[-i\tau(I-T_{ab})/2]$, where $T_{ab}$ interchanges the endpoint basis states. All
other ports are spectators. The subspace is reducing. For $\tau=\pi$ this event exchanges two basis states
exactly. Conjugating the final edge of a path of length $\ell$ by the first $\ell-1$ swaps implements the
desired exchange between its endpoints in $2\ell-1$ events. Every intermediate port is restored, even if
it carried an unknown coherent amplitude. This construction needs paths avoiding the marker.
[Seam event](../polynomial_sv/SPEC.md#8-events-and-schedules) · [Compiler](seam_program/compiler.py)

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
Other couplings and pulse lengths require their own phase check. [Blank-memory identity](../polynomial_sv/law.py)

The finite emitter starts with fuel on the marker and source seam. Its error relative to the ideal
pair is bounded in norm by $\sqrt{d_m+d_j-2}/(2g)$. Both ideal and actual states then undergo the
same seam events. The error bound carries to the total variation of any terminal observer distribution
by unitary invariance and contraction of trace distance under measurement. For a selected slot, swap
its port onto the emitter seam, run the two endpoint observers and read their public decision. The
complete seam measurement also reports every other seam and the no-pair outcome. Residual fuel
contributes to no pair. A finite-duration read under the hypotheses of the native reader adds at most
$12\tau_{\rm read}$ to this bound. The execution stores ideal instantaneous readouts after the finite
emitter. Distinct selected-slot settings use fresh preparations.
[Emitter](../polynomial_sv/sectors.py) · [Reader](../polynomial_sv/reader.py)

**Preparation obstruction for the routing pulses.** On a connected control component of $N$ nonmarker
ports, each section 8 partial exchange fixes the uniform vector. So does an unshifted pulse with
$\kappa\tau\in2\pi\mathbb Z$, which includes all routing pulses above. Starting with a sharp marked pair, its overlap with that
vector is $1/\sqrt N$ under every word of these gates. Any target orthogonal to the uniform vector
is therefore separated from the reachable pure states by at least $1/\sqrt N$ in trace distance.
This follows by projecting the reachable state onto the orthogonal complement of the uniform vector
and bounding its overlap with the target by $\sqrt{1-1/N}$. The untagged, nonzero-momentum packets of
the clock measurements have zero total amplitude. Embedding them with zero amplitude at unused ports
preserves that property. Even the finite emitter retains a lower bound of $1/\sqrt N$ minus its
preparation error, by the triangle inequality.

The unshifted law also admits nonresonant durations. On the single cell, hold the marker at port 0
and start the mobile load at port 1. A pulse of duration $\pi/2$ on $(1,5)$ gives amplitudes
$c=(1-i)/2$ and $s=(1+i)/2$. A pulse of duration $3\pi/16$ on $(4,5)$ multiplies the sum on that
edge by $e^{-8i(3\pi/16)}=i$. The total amplitude is then $c+is=0$, with unit norm and no intervening
measurement. Thus the uniform-mode bound concerns the specified routing pulses. The two-pulse example
prepares one zero-sum state. The following protocol prepares a specified arbitrary packet.

## Preparing the clock's finite wave states from native fuel

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
record. This proves that every specified finite packet, including the finite clock preparations, admits a
heralded implementation on a carrier of sufficient size. The acquisition executes one complex four-mode
packet on 72 slots. For $M$ outputs the source uses at most $6+2(M-1)+M$ nonlocal exchanges before
routing, and $M+2$ binary queries on its longest branch. Each exchange on a path of length $\ell$ uses
$2\ell-1$ seam events. An aperture query swaps its port to the reader and restores the path after the read.
The [receipt](seam_program/preparation_records/result.json) charges those edges and all four reader stages.

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
[Phase-sensitive receiver](seam_program/receiver.py) · [Complete readout](seam_program/receiver_records/result.json)
