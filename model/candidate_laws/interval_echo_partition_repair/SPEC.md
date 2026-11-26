# Interval Echo with partition repair

The law specifies what joint records an observer can obtain, how an ordinary echo changes those records, and what preparation selects their equilibrium. The observer can read the local endpoint of an aperture, the local volume counter and any enrolled payload. No transition computes an entropy. The area coefficient is determined by the number of possible completions of a local record. A second, explicitly geometric version below uses the same counting principle on surveyed proper balls.

## Local dynamics and observers

Start with the finite Interval Echo carrier. Keep its addressed free-flight unitary. At every enrolled interface add a register with legal ordered endpoint words
\[
(0,0),\ (1,1),\ (1,2),\ (2,2),\ (2,1).
\]
The endpoint order is a material address carried with its ports under a relabelling. Each material cell also has a volume bit. On an interface repair, read the word in this basis, then use two fair bits to move one place backward, stay, stay, or move one place forward around the list. A volume service reads its bit. The chances of changing it are $1/4$ for $0\to1$ and $1/2$ for $1\to0$. Both services discard the temporary command. Measuring the sharp word removes coherences between its alternatives. [Finite transition and regional reader](regional.py)

On a native contrast service, write the two port loads into the original eleven-state memory using the 99-dimensional echo unitary $U$. Whenever the written record has nonzero signed contrast, transpose aperture words $(0,0)$ and $(1,1)$. Unwrite with $U^\dagger$. The operation extends linearly to arbitrary entangled inputs. For a blank native memory and closed aperture, precisely six of the nine ordered load preparations open it. This gate couples the existing material to the aperture and is reversible. [Native service and local reservoir](candidate.py) · [Checks with an entangled reference](test_candidate.py)

For equilibrium preparation enroll native payload registers of fixed size, including port loads and seam memories. A local reservoir opportunity replaces one payload factor by its normalized identity, preserving the reduced density matrix of the other factors. It can be implemented by uniform random local shifts and phases when the command is discarded. This changes the preparation ensemble and can exchange load with the reservoir. An addressed free-flight experiment can run between reservoir services. An observer enrolling extra nontracial fields must include their entropy and dynamics in the calculation.

One finite schedule is a sweep over each free-flight and contrast gate (execute or skip with equal probability), followed by one native reservoir service on every payload factor and one sharp repair on every interface and volume bit. All command bits are independent. The fixed state is the product of the tracial payload, uniform legal pairs and volume probabilities $(2/3,1/3)$. Each service preserves it. Two sweeps with all optional gates skipped have positive probability. On that event all payload factors have been replaced, coherences in the added records have been removed, and every classical transition on each local factor has strictly positive two-step probability between every pair of states. On this finite board the resulting channel therefore contains a positive multiple of replacement by the product reference. The residual channel is trace preserving after normalization. Iterating contracts trace distance, proving uniqueness and attraction. This is a statement about the prepared record distribution. It gives no relaxation time uniform in board size.

## Exact stationarity and a necessary test

A region $R$ sees both letters of an internal pair, the inside letter of each cut pair, and the volume bits of its cells. The regional law is the marginal of one global preparation for every overlapping or disconnected selection. Let $A_R=a$ times the number of open cut endpoints and $V_R=v$ times the sum of the enrolled volume bits. These two readings are made on a word. Entropy is measured from the distribution of the entire local word. With $c$ cut pairs, $i$ interior pairs and $n$ cells the reference probability of a visible sharp word $w$ is
\[
 P_{*,R}(w)=5^{-i-c}3^{-n}2^{n+A_R(w)/a-V_R(w)/v}.
\]
A fully seen pair contributes $1/5$. A single endpoint contributes $1/5$ for symbol zero or $2/5$ for either open symbol. Including the normalized identity on fixed enrolled payload factors adds a region-dependent constant to the local negative log density. On the accessible support,
\[
 K_{*,R}=-\log P_{*,R}=C_R-\eta A_R+\lambda V_R,
 \qquad \eta=\log(2)/a,\quad \lambda=\log(2)/v.
\]
For any regional distribution or density matrix $Q_R$ on this support the definition of relative entropy now gives
\[
 S(Q_R)+\eta\langle A_R\rangle_Q-\lambda\langle V_R\rangle_Q
 = C_R-D(Q_R\Vert P_{*,R}).
\]
Consequently the reference strictly maximizes generalized entropy on the fixed-mean-volume constraint. Every differentiable perturbation through it has zero first derivative of the relative entropy, even when remote records or different faces are correlated. Thus $\delta S_R+\eta\,\delta\langle A_R\rangle=0$ at fixed measured volume. The quantum statement uses the full accessible regional density matrix on this same support. The simulated population tests use sharp records.

There is also a converse that can be tested on another law. For a faithful regional reference and all normalization- and volume-preserving tangents, exact stationarity is equivalent to $K_{*,R}+\eta A_R$ belonging to the span of the identity and $V_R$. Indeed $dS=\mathrm{Tr}(dQ\,K_{*,R})$, and the orthogonal complement of the allowed tangents is precisely that span. If the dynamics permits only a smaller tangent space, test orthogonality only on accessible variations. A scalar entropy bound or agreement on a single tangent cannot establish the operator identity. A failure can be detected by shifting probability between two readings with equal area and volume and different reference weights. The implemented bias control gives $dS=\log(2/3)$ in exactly that experiment. [All finite transfers](REGIONAL.md)

## How to force the required accuracy for surveyed proper area

The aperture unit $a$ above has operational meaning on the registered instrument. Equating its sum with the boundary area of an arbitrary smooth reconstructed ball would require evidence. The cubical comparison uses a different reader that sees a shared bit when both incident cubes are enrolled. Independent rod readings give $\eta=\log(2)/(2a)$ for that apparatus. Its fixed cubical mesh has a directional surface bias on smooth balls. [Cubical calibration and holdout](surveyed_cells/README.md)

For an explicit corrected law on *smooth* balls, enroll the two documented ultrastatic survey geometries and a fair, readable common geometry tag $g$. In each of ten specified rest frames use ranges in that metric to identify a ball of volume $4\pi\ell^3/3$. Its independently computed proper area is $A_{g}$. Each distinct enrolled boundary cuts one shared pair of material endpoints. The other pairs are entirely inside or outside, with the same membership in both geometries. Store a fixed number $N$ of equiprobable legal pair words per pair in either geometry. The **added geometric coupling** sets the number of partners of each visible endpoint to
\[
 d_g=\left\lfloor M\exp(\eta A_g)\right\rfloor,\qquad
 M=\lceil\ell^{-5}\rceil .
\]
Here area and radius are expressed in the fixed survey units and $\eta=1$ in the published comparison. Choose $N$ divisible by both degrees and at least their squares. Given a uniform rank among $N$ words, take its quotient $i$ and remainder $k$ modulo $d_g$. Store $i$ at one end and $i+k$ modulo $N/d_g$ at the other. Exactly $d_g$ words complete either endpoint. Uniform pair-rank refresh and fair tag refresh select the reference for this declared apparatus. The tag is shared across the entire finite board. Communicating or generating it from native material is an additional dynamical problem.

Conditional on the common tag the ball sees $N/d_g$ possibilities at its one cut. Whole pairs contribute a tag-independent number of possibilities. Thus for an arbitrary normalized source tangent, including variations within each tag or correlations between pairs, $\delta S=-\sum_g\delta p_g\log d_g$. Its independently surveyed area slope is $\delta A=\sum_g\delta p_g A_g$, with $\delta V=0$. Since $|\log M+\eta A_g-\log d_g|\le 2/M$, the same source supplies
\[
 |\delta S+\eta\,\delta A|\le 2\ell^5\,\|\delta p_{\rm global}\|_1=o(\ell^4)
\]
for uniformly bounded tangents. This holds for every enrolled frame, on successively finer finite boards under the specified placement. Coarsely measured input areas contribute their own measurement error. The numerical order comparison does not certify rigorous integration intervals. [Geometry, compatibility on overlaps, signed responses and controls](proper_balls/README.md)

This rule makes a falsifiable constitutive choice. It predicts the ratio of compatible endpoint counts from a survey made without reading them. Reversing the two degrees holds the surveyed area fixed and reverses the entropy response. The construction supplies the metric and the geometry-dependent degrees. It does not derive a metric from a freely formed native mass. No identification of the calibrated regional stress of ordinary matter or of the collar with the same entropy slope has been established. The second derivative also contains relative entropy. For the fair geometry mixture, moving the tag probability from one half to five eighths incurs a finite relative entropy of at least 0.03158 nats at every radius. The symmetric secant cancels this even term. A straight-line count in the comoving representation also gives the wrong sign for the proper-area response of the tilted curved balls. [Finite-displacement bounds](REGIONAL.md) · [Gravity implication and outstanding physical inputs](../../../code/matched_einstein_limit/PROOF.md)
