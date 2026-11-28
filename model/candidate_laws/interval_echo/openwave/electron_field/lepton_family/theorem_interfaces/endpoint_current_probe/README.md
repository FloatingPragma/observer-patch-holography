# A neutral-current test on the charged Wilson endpoint

A signed interferometer on the eight-cell Wilson source gives the **same**
phase for pure hypercharge and for $T_3+Y$ on both charged spinor endpoints.
The two settings give different records on independently prepared weak
doublet controls on the same cube. The distinction uses a calibrated
auxiliary weak frame and a common receiver. The unchanged spinor sequence
by itself determines the hypercharge coefficient and leaves the weak
coefficient free. [Method](METHOD.md) · [Exact receipt](records.json)

The interferometer compares the original link-dressed transfer with a copy
conjugated by the receiver's current pulse. A signed charge counter marks
the transfer at both ends without reading which branch of the interferometer
ran. Summing the two phase ports recovers the original charge, energy and
idle-reference records. The original input is the full polynomial preparation
from the immutable Wilson source. It is the same source whose charged
endpoints carry the calibrated Dirac spin rotations. [Charged endpoint](../charged_endpoints/README.md)

The current settings give the following probabilities for the positive
sine port at exposure $\pi/2$. Each row is conditioned on its indicated
transfer, with a separate refusal record. A particle or hole here names the
occupied or emptied beta mode at the receiver. The signed native increment
for the doublet particle has the opposite orientation to the spinor particle.

| Accepted transfer at receiver | Pure $Y$ | $T_3+Y$ |
| --- | ---: | ---: |
| Spinor, native $+6$ | 1 | 1 |
| Spinor, native $-6$ | 0 | 0 |
| Independent upper doublet particle, native $-3$ | 0.146447 | 0.5 |
| Held lower doublet particle, native $-3$ | 0.146447 | 0 |
| Independent upper doublet hole, native $+3$ | 0.853553 | 0.5 |
| Held lower doublet hole, native $+3$ | 0.853553 | 1 |

With a uniform original receiver address, each spinor sign occurs on more
than **9.93083175%** of original starts, including more than **3.35939%** of
original starts in the archived transfer-energy band. In each independent
doublet bank the local stock accepts more than **78.320864%** of original
starts on the same polynomial preparation. This gives a sine-bank
distinction per original start greater than **0.2769** on the upper control
and **0.1146** on the held lower control, even with refusals included. These
are rigorous weight bounds and exact finite port calculations. [Protocol](PROTOCOL.md)

An observer who loses access to the independent frame receives identical
source records from the two control preparations. A common signed singlet
read and an upper-doublet null identify $T_3+Y$ within the calibrated
two-generator plane, and predict the held lower read. The physical
electromagnetic identification requires showing independently that the
source's own neutral propagation and symmetry-breaking response select
that same direction. [The reconstruction theorem](../EndpointNeutralDirection.lean)

From the repository root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.theorem_interfaces.endpoint_current_probe.measure
~~~
