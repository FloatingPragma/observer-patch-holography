# An entropy maximum on the complete two-head state

On the registered 80-cell radar experiment, a faithful preparation is the unique global maximum of stored-head entropy plus calibrated tetrahedral area at the measured comparison volume. The twelve directed ranges come from independent active-law probe effects. The gain comes from a separate ordinary exchange standard. The maximum covers all 31 parity-block coordinates, and pinching extends it to arbitrary coherences on the same eight-dimensional support. A positive mixture of native product preparations reaches the state. [Proof](PROOF.md) · [Certified bounds](certificate.json) · [State, source weights and held measurements](results.json)

Every other state at the same survey volume loses at least 0.14 times its squared trace-norm distance in the objective, measured in nats. A second faithful native preparation at that exact volume has a strictly positive first derivative between 7.47 and 7.49 nats per unit odd fraction. The source lotteries select which state the observers survey. A feedback calculation starts from the tracial state and uses the ordinary calibration and volume target to set those lotteries. For any number of disjoint copies, the product preparation uniquely maximizes the joint objective under both separate and pooled survey readings. The range cycle ratios at that equilibrium differ from one even after adjusting each station's clock scale. The same state and gain also fail stationarity when the clock timeout changes from six to four units.

On the same connected board, unconditioned, uniformly launched short probes relax the head coins to sector traces. One of these stationary states has the calibrated volume and a nonzero fixed-volume entropy-plus-area slope. Reusing the variational preparation for one full capped radar cycle decreases its left-coin reading on every charge and route. [Connected equilibration](PROOF.md#equilibration-on-one-connected-board) · [Complete radar readback](reuse/README.md)

### How much record capacity does the survey equilibrium use?

One joint coin observable certifies that the selected survey equilibrium is
**more than 0.11033 nats below its entropy ceiling at its own reversal
fraction**. Let $p$ be that fraction and $z=\langle Z_L Z_R\rangle$. Each
reversal sector contains four coin states. Its two possible $ZZ$ outcomes
each have multiplicity two. Write $h$ for binary entropy in nats. The sharp
maximum entropy given both readouts is

\[
S_{\max}(p,z)=h(p)+h((1+z)/2)+\log 2.
\]

To see this, dephase the two commuting readouts. Dephasing can only increase
entropy. Each resulting conditional coin block has entropy at most $\log 2$,
and the joint entropy of the readouts is at most the sum of their binary
entropies. Equality is attained by independent readouts with uniformly
populated conditional blocks. The native preparation can realize this
maximum for $0\le p\le1/2$: draw the left $Z$ bit fairly, draw an independent
product sign of mean $z$ to set the right bit, and load both signed rails.
Apply the unread left-sign write on an independent fraction $2p$ of the
trials. Omitting the $ZZ$ constraint gives the fixed-parity ceiling
$h(p)+\log 4$, attained by the sector traces selected by the uniform probe bath.
[Preparation and native pulses](PROOF.md#native-preparation-and-held-responses)

The existing native enclosure for the **ideal selected six-unit survey
reference** gives $0.4608947611<z<0.4608959612$. It includes uncertainty
in its rational density center and the certified distance to the true
reference. The bath predicts $z=0$ in every parity mixture. For $|z|<1$,

\[
\log2-h((1+z)/2)=\sum_{k=1}^{\infty}
  \frac{z^{2k}}{(2k)(2k-1)} .
\]

This follows by integrating the power series for $\operatorname{atanh}z$
from zero. The first five terms at $z=0.4608947$ already exceed $0.11033$
by exact rational comparison. This is a population bound on the ideal
reference. For an approximate laboratory realization, its own lower
confidence bound on $|z|$ can be substituted in the same inequality.
The quantity bounded here is ordinary stored-head entropy. The separate
fixed-volume theorem maximizes $S+\eta A$. [Native correlation and rational certificate](pulse_bath.json) · [Verifier](pulse_bath.py)

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.station
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.full_algebra.certificate
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/full_algebra
```

From the radar timestamps alone, 31 pre-cutoff return probabilities identify the full two-head state, and calibrated clock weights predict the independently measured signed force on both carriers. The twelve means and both force summaries together miss an entropy-increasing change. [Complete clock measurement and proof](PROOF.md#recovering-the-complete-source-from-the-original-clocks)
