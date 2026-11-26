# A signed test on two marked ports

**Prospective physical measurement.** The two-dimensional component and its one-dimensional phase reference admit four stratum-labeled assignments in the [PR-12 ballot](records.json). A native current with independently fixed port labels and polarities would let an observer separate two kinds of assignment. The numbers below are forecasts. The [P1 apparatus replay](port_pair/README.md) checks signed receiver counts with programmed P1 generators and an inverted port-1 control. Testing a naturally occurring current requires identifying its signed coupling to each marked port independently of that compiler.

## Which angles do the ports make?

In the frozen port numbering let `A` be the adjacency matrix of the lexicographically first 60-element orbit of ordered port pairs. Write $P_+$ and $P_-$ for its rank-three orthogonal projectors at eigenvalues $+\sqrt{5}$ and $-\sqrt{5}$. For either projector $P$, the direction induced by port $i$ on an active spinor is $n_i=P e_i/\|P e_i\|$, up to a common orthogonal frame and generator sign. The active compact factor acts by the spinor lift of its adjoint action. Its invariant norm is a scalar multiple of the port-sector norm, so the common scale cancels. The pinned projectors give $P_{00}=P_{11}=1/4$ and $P_{01}=\pm\sqrt{5}/20$. Thus $c=n_1\cdot n_0=\pm1/\sqrt{5}$. Interchanging the name of the two three-dimensional sectors without reattaching the marked ports loses this sign. [Exact B14 projectors](../../../../code/b14_jacobi/port_dual_metric_selector.py) · [Executable predictions](test_marked_ports.py)

| Matching representation on doublet plus singlet | Active sector | $c$ |
|---|---|---|
| P1, sampled plane point $(1,0)$ | plus | $+1/\sqrt{5}$ |
| P2, doublet on the first factor, if that factor is attached to plus | plus | $+1/\sqrt{5}$ |
| F2, trivial on $A_2$ | minus | $-1/\sqrt{5}$ |
| G2, trivial on $A_2$ | plus | $+1/\sqrt{5}$ |

On F2, the plus and five-dimensional sectors close together into $A_2$, leaving the minus triple as $A_1$. On G2 their roles interchange. On each plane axis the nonzero self-bracket selects the active triple. P1 also contains the other axis. A P2 apparatus with the doublet attached to minus gives the minus forecast. The abstract factor names alone leave that attachment free. The equal-weight oriented-face projection under the declared squared-distance repair picks the G point with $e=\sqrt{5}/20$ and $ab>0$, hence G2. Using this repair to select a source requires evidence that the material implements it. [Compact-family equations](../../../../code/b14_jacobi/classify.py) · [Projection parameters](../../../../code/b14_jacobi/oriented_face_bracket_selector.py) · [Declared measure and repair](../../../../code/b14_jacobi/port_dual_metric_selector.certificate.json)

## What the receiver would count

First mark ports 0 and 1 in the same material frame at two stations. Establish their current polarities independently of the doublet fringe and fix the gain and sign of each half-turn by separate continuous, signed pulse calibrations against the reference. A projective measurement on the doublet alone cannot distinguish $+\pi$ from $-\pi$. Check that the reference line stays fixed, and calibrate its propagation phase with the no-pulse source.

Prepare $(|s\rangle+|d\rangle)/\sqrt2$, where $s$ is the singlet and $d$ the positive eigenvector of the port-0 Pauli axis. Apply two $\pi$ rotations, first port 0, then port 1. With $U_i=-i\,n_i\cdot\sigma$, their product is
$$U_1 U_0=-c I-i(n_1\times n_0)\cdot\sigma.$$
The cross product is perpendicular to the port-0 eigenaxis, so $\langle d|U_1U_0|d\rangle=-c$. The terminal analyzer resolves the bright and dark vectors $(s\pm d)/\sqrt2$ and the orthogonal doublet address. The ideal, unit-transfer probabilities are

| $c$ | Bright $(1-c)^2/4$ | Dark $(1+c)^2/4$ | Orthogonal $(1-c^2)/2$ |
|---|---:|---:|---:|
| $+1/\sqrt{5}$ | 0.0763932 | 0.5236068 | 0.4000000 |
| $-1/\sqrt{5}$ | 0.5236068 | 0.0763932 | 0.4000000 |

Bright minus dark is $-c$ even for an arbitrary input doublet vector, since the cross term is pure imaginary. The separate three probabilities above use the specified eigenvector. Record every start and keep refusal or absent remote arrival as a fourth outcome. Independently check transfer efficiency when comparing absolute contrasts. Reversing just the second pulse swaps bright and dark. Exchanging the order of the two positive pulses keeps the real contrast. Removing both pulses puts all ideal counts in bright. Two full $2\pi$ turns return the doublet to identity.

For identical normalized port generators up to one common doublet frame, every ordered pulse product and every equivalently mounted terminal measurement has the same probabilities. Operators in the representation kernel act trivially on both arms. Repetitions on this specimen cannot distinguish algebras that differ only inside that kernel. A second detector carrying the other factor, or an independent measurement of brackets of the twelve currents, supplies a different observable. For the classified compact sources their derived dimensions are 3 on P1, 6 on P2 and 11 on F2 or G2. The sector that acts on the five distinguishes the last two. This protocol concerns internal port currents. A physical charged-lepton assignment also requires identification of the same freely traveling charge, spatial spin response and rest energy.

A separately calibrated native current now resolves both marked-port cosines on the isolated cell. A one- or two-turn pentagon launch gives opposite signed circulation at mark 1, with the sign verified by subsequent endpoint transfer. Moving the antipodal comparison from the receiver to the source reproduces the same isolated readings with an ordinary bond meter. Opening the federation changes the profiles. The coupling of either measured orbital current to a freely traveling charged doublet remains open. [Native sources, bare receivers and controls](ring_current/README.md)

The following command checks the forecasts from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m pytest -q model/experiments/interval_echo_gauge/matter_class/test_marked_ports.py
~~~
