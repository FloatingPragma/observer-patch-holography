# What two freely moving loads record at a three-arm junction

**IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.**

Prepare the same signed arrival at 1045, supply a second identical load at 1028 and open six material bonds together. At a local wait of $3\pi/8$, observers recover both original charge addresses and all six initial tree heads on **25.4885%** of enrolled native preparations. A fully spin-polarized fermion comparator with the same transfers, hard-core exclusion and port energies gives **2.56476%**. Cutting the spare arm makes all address probabilities identical in the two laws at every wait. A spin-singlet comparator reproduces the native statistics on the open junction exactly. These comparisons establish which spin preparation this particular measurement can distinguish. [Protocol and selection](DESIGN.md) · [All address counts](records.json) · [Matched routed exchange](../junction/README.md)

## What the counters register

The three physical arms from 1045 terminate at 1010, at 1028 via 1044, and at 1032 via 1052 and 1047. Both incoming 1045 face words already have the required internal tree heads on the supplied full-cell completion. The last bridge of the initial guide supplies its recorded direction. The other external heads and the second charge are supplied. During local flight every selected material bond has rate four, each port retains its intact-board potential $8(6-d)$, and every other bond and ring flip is shuttered. The two marker faces on the routed comparison are (1045,1049,1053) and (1021,1027,1028). Both lie off the six moving head bits and retain their original word differences. [Complete preparation and marker pairing](../junction/README.md) · [Source read at 1045](../../../fan_path/local_launch/WITNESS.md)

Two material counters report the unordered pair of occupied ports. A second bank also records the six active head bits. Removing an edge of the tree cuts off a set of ports. Every active edge internal to that set has exactly one head inside it. Conservation of incident heads minus signed load on the set fixes the remaining cut head from the total load inside. Thus each reachable address pair determines the complete active head word. In particular, a return to the initial pair restores those heads irrespective of how the two origins moved. All other written heads are spectators. Charge conjugation complements the six bits and gives the same position probabilities.

| At local wait $3\pi/8$ | Native | Polarized fermion | Distinguishable origin tags |
| --- | ---: | ---: | ---: |
| Six material contacts open | 25.4885% | 2.56476% | 14.0266% |
| Transfers to an occupied owner cell shuttered | 16.2088% | 9.08647% | 12.6477% |
| Spare arm cut | 2.24307% | 2.24307% | 2.24307% |
| Spare arm cut and occupied-owner transfers shuttered | 44.3958% | 44.3958% | 44.3958% |

Each percentage conditions on both inputs. With the second preparation and completion supplied with certainty, the first line's native and polarized simultaneous 99% intervals per **original** launch are **[0.40359%, 0.40404%]** and **[0.04058%, 0.04073%]** from the simulated banks. With the owner restriction they are **[0.25665%, 0.25701%]** and **[0.14384%, 0.14412%]**. Every bank includes all other address pairs and the initial refusal. The additional owner shutter removes a complete otherwise-enabled transfer, including its Laplacian diagonal. The ungated apparatus admits multiple ports in one owner subject to hard-core exclusion at each port.

## Why the return distinguishes exchange

In calibration copies let $A$ and $B$ be the complex amplitudes to recover the initial two addresses with origins in their original and exchanged order. Both contributions end with exactly the same head record. The native indistinguishable-material probability is $|A+B|^2$. The polarized comparator gives $|A-B|^2$. Keeping distinct origin labels gives $|A|^2+|B|^2$. Reading an origin in every trial would change the interference experiment. At $3\pi/8$ the distinguishable-origin calibration gives **3.35707%** for a swapped return with all contacts open, and **0.305439%** with the owner restriction. The two unlabelled laws at the later appointment π even reverse their ordering in the open case. [Individual waits and all refusal counts](records.json)

The native and polarized generators have the same diagonal. The canonical fermionic hopping sign changes when the other occupied orbital lies strictly between the origin and destination in the chosen ordering. Changing that orbital ordering conjugates the signed generator and the sharp endpoint readings by site-basis phases. An independent calculation of the distinguishable-origin amplitudes reproduces both signs. With the spare arm cut, order along the remaining physical path is conserved. A diagonal gauge then identifies every address probability for arbitrary waiting time. With the branch open, the shortest origin-exchanging return needs eight material hops, or twelve under the owner restriction. Integer matrix powers give the first nonzero terms in $P_{\rm polarized}-P_{\rm native}$ as $-16384 u^8/315$ and $-131072 u^{12}/467775$, respectively. [Exact finite generator, signed sectors and path-cut check](measure.py)

An additional opposite-spin Fock construction enrolls two spin orbitals per physical address and prohibits double occupancy at a port. Direct integer matrix products show that spin-independent hopping intertwines the spin singlet with the native spatial generator, and the $m=0$ triplet with the polarized generator. This holds at every wait and for every enrolled address distribution with the corresponding coherent spatial preparation. The field-marker operators act on fixed seam words throughout this release. A conditional return can preserve the original reference correlation at its old site even when a separate calibration shows the material origins exchanged. The two-reference measurement in the routed comparison resolves this distinction. Physical electron identity still needs a moving internal degree of freedom with the appropriate rotation and exchange readings on the same carrier. [Routed marker measurements](../junction/README.md) · [Five full-cell poses](../face_turn/FRAME.md)

The two origin assignments give an independently constructed double cover for the phase comparison. Every floating propagator is checked against a matrix exponential, with spectral-frame orthogonality error below $10^{-13}$ and residual below $10^{-12}$. Duhamel's bound at $u\le\pi$ fits within the $2\times10^{-10}$ probability guard. The count intervals include the independently bounded incoming face probability, the probability-grid rounding and a union over all reported return banks. [Independent rational return check](exact_return.py)

Run from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.electron_field.dressed_exchange.intact_path.junction_flight
.venv/bin/python -B -W error -m "$pkg.measure"
.venv/bin/python -B -W error -m "$pkg.exact_return"
~~~
