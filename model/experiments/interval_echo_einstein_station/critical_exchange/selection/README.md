# Which matter reference does an Interval Echo observer select?

An observer on the 246-station guide can cool **every half-filled input** to the same correlated reference by accepting only exchanges that lower independently measured hopping energy. The observer can start with a sharp alternating load at the stations. The procedure uses the active law's native contacts, a stored comparison and a destination for the discarded preparation information. It uses no regional entropy, area, or desired population histogram as input. On release into the full federation, one exterior receiver detects leakage from this guide. [Frozen choices](PROTOCOL.md) · [Observer protocol](measure.py) · [Counts and selection trajectories](results.json) · [Action and release control](../README.md)

## Independent work calibration

Between adjacent stations the addressed forward and inverse contact sequence realizes the hopping action with matrix element $-1/2$. In the odd, 123-load sector the closure of the ring is periodic and its one-orbital energies are $\epsilon_k=-\cos k$. A native sequence of adjacent rotations maps the momentum orbitals to distinct readout addresses and back. A sharp occupation read at each resulting address measures a complete word and preserves the fixed total load. A swap of two such addresses is an allowed transfer after the word has been read. The rest of the federation is parked during this guided service.

Apply the symmetric/antisymmetric analyzer on a single internal bond of independent one-orbital preparations. Over all launches the three probabilities are $(1-\epsilon_k)/246$, $(1+\epsilon_k)/246$, and $1-2/246$. Their difference among the hits gives the signed energy and the quiet count verifies the aperture. Conditional on the hit count $n_k$, the sign has range two. Hoeffding's inequality and the union bound give all 246 energy intervals simultaneous coverage of at least 99% with half-width $\sqrt{2\log(492/0.01)/n_k}$. These are simulated counts. The electronic protocol accepts a transfer if the entire destination interval lies below the source interval. The frozen banks leave a gap greater than 0.02214 exchange units between the two half-filled shells. The analytic gap is $2\sin(\pi/246)=0.02554$.

## Exact selection from every input

After the first word read, let $m$ be the number of particles in the upper shell. There are equally many holes in the lower shell. Each proposal chooses an unordered pair uniformly from the $M=\binom{246}{2}=30135$ pairs. Exactly $m^2$ pairs connect one of these particles to one of these holes. They all pass the independent work test. Every other accepted comparison preserves $m$. Thus $m$ falls by one with probability $m^2/M$ at each attempt until it hits zero. The sum of the geometric mean waiting times is
\[
 \mathbb E[T\mid m]=M\sum_{j=1}^{m}j^{-2}.
\]
The largest expectation is 49,326.1 proposals. The bound $\mathbb E[m_{s+1}\mid m_s]\le(1-1/M)m_s$ gives $\Pr(T>s)\le123\exp(-s/M)$, which is below $1.385\times10^{-5}$ at the finite cutoff $s=482160$. Comparisons after absorption have no effect. At absorption, the 123 occupied negative modes constitute a *single* Slater word. Rotating the addresses back gives precisely the state read on all the overlapping intervals.

The proof includes non-Gaussian inputs and inputs entangled with an external reference. The orbital read produces a classical word correlated with that reference. Every branch cools to the same output. Summing all branches returns the original reference marginal tensor the cold material state. At finite cutoff the trace distance to this product is at most the unfinished probability. The controller retains the discarded information. The same pair count applies to any finite number-conserving one-particle action with a separated filling and an accessible mode analyzer.

| Input read before the comparisons | Upper particles initially | Pair proposals to absorption | Hopping work released |
|---|---:|---:|---:|
| Inverted filling | 123 | 11,828 | $+156.6127$ |
| Uniform half-filled mixture | 68 | 21,119 | $+85.1642$ |
| Sharp alternating station loads | 70 | 128,334 | $+86.3768$ |
| Cold filling with the work comparison reversed | 0 | 38,020 | $-156.6127$ |

For alternating sharp loads each momentum doublet $(k,k+\pi)$ contains one fermion in an equal coherent superposition. The orbital read therefore gives independent fair choices in the 123 doublets. The assigned hopping work uses the native action, and separate lower and upper work estimates in the archive use only the signed calibration counts. Under the reversed comparison all 123 particles finish above the gap. Parking the comparator leaves every initial orbital occupation invariant under the addressed hopping dynamics.

## A basin identified from energy density

There is an operational selection criterion independent of this cooling algorithm. On a ring of length $L$ at the same filling, let $\Delta E=\langle H\rangle-E_{\min}$. Let $d_k$ denote the occupation above the band edge, or the hole occupation below it. For every state, $\Delta E=\sum_k|\epsilon_k|\langle d_k\rangle$. For $0<g<1$,
\[
 \bar d=\frac1L\sum_k\langle d_k\rangle
 \le \frac{\Delta E}{Lg}+\frac{2}{\pi}\arcsin(g)+\frac4L.
\]
The last two terms count the equally spaced modes with $|\cos k|<g$. Averaging the origin of an observation uniformly projects its one-particle covariance onto the diagonal in momentum. Every separation then has covariance error at most $\bar d$ compared with the cold reference.

The criterion extends to the complete record of any fixed finite set of sites. Write each site annihilator as a particle annihilator above the edge plus a hole creator below it. In the translation-averaged preparation, the expectation of each of these particle or hole number operators at a fixed site is bounded by $\bar d$. Expand a local word in these operators, use their canonical anticommutation relations to order creators to the left, and subtract its reference expectation. The scalar contractions cancel. Each remaining term either starts with a creator or ends with an annihilator. Cauchy–Schwarz bounds its expectation by $\sqrt{\bar d}$ since the mode norms are at most one. For a fixed word the finite number and coefficients of these terms are bounded independently of $L$. Thus vanishing measured excess energy per site as $L$ grows selects all fixed finite, translation-averaged records of the cold state. One can take $g\to0$ after $\Delta E/L\to0$. A particular origin in a localized hot inclusion requires its own test. Arbitrarily long physical rings require their own carrier realization.

## Relation to the gravity balance

On the resulting reference, the independent hopping reader and the interval entropy reproduce both signs of a local energy change on neighboring overlapping intervals. This identifies a selectable *matter* source from records. The higher-dimensional Einstein test also needs the area read at fixed proper volume on the same source. [Interval measurements, finite spectral discrepancy and native all-open experiment](../README.md)

Even a collection of independently cooled ideal conformal strands would have a measurable dimensional discrepancy. In a rest ball of radius $R$ each parallel strand cuts a chord of half-length $r=\sqrt{R^2-s^2}$ at transverse offset $s$. In a homogeneous excitation of strand energy density $t$, the interval first law gives $(4\pi/3)r^2t$. Integrating over the transverse disk gives $(2\pi^2/3)R^4t$ at unit strand density. The three-dimensional ball requires $(8\pi^2/15)R^4t$. For observer velocity $u$, put $b=\sqrt{1+u_x^2+u_y^2}$ transverse to the strands and let $t$ denote the strand energy in the unit time direction projected onto its sheet. A longitudinal Lorentz change of coordinates makes the observer's time vector $(b,u_x,u_y,0)$. On its rest plane the squared transverse distance to a strand at $X$ is $|X|^2-(u_\perp\cdot X)^2/b^2$. The number of strands per observed transverse area therefore increases by $b$, while the action energy contracted twice with $u$ increases by $b^2$. After a single rest calibration the entropy-to-bulk ratio is $1/b$. A transverse boost with $u_x=4/3$ gives $3/5$. This discrepancy assumes independent strand algebras and homogeneous stress of positive energy. It supplies a reason to test transverse correlations directly. [Conformal ball generator](https://arxiv.org/abs/1102.0440) · [Three-dimensional correlated comparison](../../../../../code/matched_einstein_limit/source_implication/modular_orbit/weyl_balls/PROOF.md)

Run from the repository root with one numerical thread.

~~~sh
python -m model.experiments.interval_echo_einstein_station.critical_exchange.selection.measure
python -m pytest -q model/experiments/interval_echo_einstein_station/critical_exchange/selection
~~~
