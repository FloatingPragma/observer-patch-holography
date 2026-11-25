# The same local charge ladder acts on all three winding characters

A single set of material contacts beside the receiving hub converts an arriving charge into a **regionally neutral vacancy–excess pair while transmitting its qutrit**. This works for both signs with conjugate preparation of the guard and donor. The originally transported vacancy or excess reads the three winding phases after conversion. The second member of the pair comes from a locally stocked donor. The measured charge transfer and winding preservation are properties of this driven instrument. [Local geometry and stock](apparatus.py) · [Native pulse and acquisition](measure.py) · [Complete records](results.json)

## What the receiver measures

After 32 units of ordinary free flight, a sharp signed arrival at port 160 selects the same carrier as in the [three-character transmission experiment](../README.md). The target 167 is an unused fifth neighbour of the hub, outside the three winding arms. The target, donor 166 and guard 158 form a local workstation with the hub. The receiver uses the established selective ladder on these addresses. Every original direction head remains untouched while its ternary endpoints can serve other contacts. The apparatus isolates all unaddressed contacts for this intervention. It opens its three incident equal-load shutters only after the signed arrival has been parked on blank spoke 160–165. A single exterior service opening supplies the held stock. [Registered ladder instrument](../../../../../interval_echo_gauge/arrival_ladder/README.md)

Two separate balanced $(2,0)$ pairs write even magnitude-two collars on 160–167 and 158–167. The ordinary private reference at 167 supplies the preparation phases. Two further compensated pairs provide the live donor and guard. All stock routes are physical blank-contact transpositions, which restore every intermediate port. The original captured defect returns to the hub before conversion. The full occupied and compensating board has zero net excess.

The operator reads the actual triple **(hub, target, donor)**. A vacancy starts as $(0,1,2)$ and converts to $(0,2,1)$. A positive defect with the reflected stock starts as $(2,1,0)$ and converts to $(2,0,1)$. Thus the hub plus target gains one unit in the first preparation and loses one in the second, arriving at zero excess in both cases. The donor supplies the opposite change. Each of the three loads is independently relayed to a calibrated signed-echo read against a zero reference. The complete count bank includes all 27 ternary words, missed arrivals and source refusals.

## Why the winding amplitudes survive

Call the two words $c$ and $n$. The temporary echo of the hub against the guard concentrates their extreme contrast on unit loads. Each stored magnitude-two collar costs zero on $(1,1)$ and 32 otherwise. The two prescribed dwells and the inverse echo give opposite phases $e^{i\phi/2}$ and $e^{-i\phi/2}$ for the two target occupations. A half blank exchange on the target and donor before those phases, followed by its inverse afterward, gives

$$
U(\phi)\big|_{\{c,n\}}=
\begin{pmatrix}
\cos(\phi/2)&-\sin(\phi/2)\\
\sin(\phi/2)&\cos(\phi/2)
\end{pmatrix}.
$$

The code computes the temporary echo on all 99 native states and checks restoration of its blank and reference. Every applied contact avoids the three unknown winding heads. The neutral-branch map therefore acts as $\sin(\phi/2)I$ on all eight receiving head words, including the full three-character code and any reference entangled with it. At $\phi=\pi$ it transfers each winding and their coherent superpositions with the same phase. The restriction commutes with the cyclic character read; the actual workstation still uses a particular spoke.

After reading the donor, its port and the guard port are returned to unit filling. Target 167 stays parked with load 2 or 0. The original vacancy or excess is still at the hub. It can carry out its own winding interferometer, or the common extractor can park its qutrit on 164, 162, 156. The ten active terminal ports together with target 167 contain one vacancy and one excess on the neutral branch. An additional probe used for a nondestructive unwound check can be returned to its stock. This larger aperture should be distinguished from the two-port region used for the conversion charge difference.

| Incident and material stock | Conversion conditional on arrival | Lower qutrit fidelity on neutral branch, simultaneous 99% |
| --- | ---: | ---: |
| Vacancy, positive stock, half setting | 1/2 | 0.69749 |
| Vacancy, positive stock, full setting | 1 | 0.73350 |
| Excess, negative stock, half setting | 1/2 | 0.69875 |
| Excess, negative stock, full setting | 1 | 0.73312 |

Each source and phase analyzer is applied in both the sharp-arm and three-character Fourier bases. Each setting has $2^{38}$ original simulated starts with both instrument apertures and every failure retained. Twice the same half pulse, with the donor left unread between them, gives full conversion. Recording the intermediate donor word changes the final conversion probability to $1/2$. Fresh signed-load banks implement both comparisons. The wrong incident sign on a fixed stock, the opposite guard alone, a depleted donor or a filled hub gives no conversion to the specified neutral pair. Fresh neutral stock reverses with unit probability. At half the registered contrast and frozen full-setting times the conversion probability is $1/2$. [Control and phase banks](results.json)

The measured regional charge transfer is a candidate matter current for this apparatus. The opposite-sign target was supplied locally and the run reads the pair at the conversion site. Independent free propagation of that neutral pair, a chiral spacetime response, an independently identified gauge current, and an autonomous flavor-changing decay require additional experiments. All three windings passed the same ladder. Their free-energy evidence is in the [parent experiment](../README.md).

From the repository root run the following commands.

~~~sh
p=model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge.neutral_ladder
d=model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/neutral_ladder
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m "$p.measure"
python -m pytest -q "$d"
~~~
