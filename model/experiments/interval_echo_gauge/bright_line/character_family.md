# What can one colour determinant predict about the pair phase?

On the antipodal \(\pi\) pulse the one-load colour determinant is exactly \(i\). On the same pulse, an observer's return counts put the real part of the pair bright-line phase in **[−0.65555, −0.62265]** (outward rounding, simultaneous 99% coverage). Every integer power of \(i\) has real part \(-1\), \(0\), or \(1\). This one pulse therefore excludes the entire hypothesis that the bright phase of every line-preserving control is a continuous single-valued character of the colour block. [Native bands and phase detector](README.md) · [Count and eigenvalue check](character_check.py)

Indeed, every such character on \(U(3)\) factors through its determinant: it maps commutators to one and the commutator subgroup is \(SU(3)\). A continuous character of the remaining circle has the form \(z\mapsto z^k\), \(k\in\mathbb Z\), by lifting its angle from the real line. Thus if the phase is fixed by a character it must be \(\eta=(\det U)^k\) for one integer on the whole control family. A purely empirical interval around \(i\) would allow large powers to wrap around the circle. Here the determinant follows exactly from the trace \(1/4\) on the first one-load triplet for each of two antipodal controls and \(\det\exp(i\pi H)=\exp(i\pi\mathop{\rm tr}H)\). The native pair bright eigenvalue lies strictly between \(1/2\) and 1, giving a second, spectral derivation of the exclusion. The count check rebuilds the three pair return means from all seam and quiet bins, applies simultaneous binomial limits, and includes the finite phase-anchor error.

## Which subgroup results if a controller imposes a character?

Keep the two rotation triplets fixed. Write \(U\in U(3)\) for the one-load colour action, \(W\in U(2)\) for the action on the pair plane perpendicular to the bright line, and \(\eta\) for the bright phase. The full joint current group gives \(\eta\det W=(\det U)^2\). A controller that confines its words to \(\eta=(\det U)^k\) therefore acts through

$$G_k=\{(U,W)\in U(3)\times U(2):\det W=(\det U)^m\},\qquad m=2-k.$$

For even \(m=2l\), every element is uniquely \((U,(\det U)^l B)\), with \(B\in SU(2)\). This gives \(U(3)\times SU(2)\), whose primitive product cover has kernel of order three. For odd \(m\), the map

$$(A,B,z)\longmapsto (z^2 A,z^{3m}B),\qquad A\in SU(3),\ B\in SU(2),\ z\in U(1),$$

is onto: for a given \(U\), choose \(z^6=\det U\) and remove those scalar factors. Its kernel consists of the six triples \((z^{-2}I_3,z^{-3m}I_2,z)\) with \(z^6=1\). A primitive sixth root acts as \(-I_2\) in the weak factor. The case \(k=1\) thus gives the familiar diagonal order-six kernel. Many other odd choices give that same abstract kernel and change the physical weights of the two triplets. The kernel order alone cannot determine the phase rule.

These are graph subgroups selected inside the available current image. The exclusion concerns a phase prescribed on *every* line-preserving control by the full colour block alone. A restricted pulse menu, a separate variable carrying phase, or a relation on a covering group has different hypotheses. No measurement here selects a physical Standard Model group.

Run `.venv/bin/python -B -W error -m model.experiments.interval_echo_gauge.bright_line.character_check` from the repository root.
