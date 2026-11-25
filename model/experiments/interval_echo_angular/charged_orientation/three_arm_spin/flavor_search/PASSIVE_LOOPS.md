# What can one charge record on a closed sequence of native contacts?

Consider a vacancy or excess on a unit-filled host. Route it once around a simple physical cycle by applying one addressed native contact at a time. Each transition has a sharp source and destination port. The successful branch records the prescribed arrival at **every** step. An observer who only checks the final port can collect other paths. The two experiments have different operators.

For the prescribed route, any head state or code which returns reversibly to itself on successive successful laps is supported entirely on blank traversed seams. Its return amplitude equals the ordinary blank-route amplitude. The result allows arbitrary coherent superpositions of all head words, entanglement with the rest of the material, and a retained quantum reference. The fixed-current winding calculation gives the restriction to one initialized, all-written sector. [Fixed-current winding proof](../../winding_memory/PROOF.md)

## Single-contact calculation

Let $B_e$ be the destination-from-source block of the registered local Hamiltonian, or of an isolated contact pulse $\exp(-it_e h_e)$. This block acts on the eleven head states. In the head basis of definite contrast magnitude and reversal parity, it has the form

$$B_e=b_eP_e+N_e,\qquad P_eN_e=N_eP_e=N_e^2=0,$$

where $P_e$ projects onto the blank. Both magnitude-one parities have one directed coin transition. For a vacancy crossing the oriented contact forward, they take coin 0 to coin 1. Traversing in reverse or replacing the vacancy by an excess interchanges the input and output coins. The generator has amplitude $-4$ in the even sector and $-9/2$ in the odd sector. For a blank contact $b_e=-1/2$. Magnitudes zero and two have zero oriented transfer.

Exact integer compression of $8h_e$ checks all 99 native coordinates. The blank source and its destination form an isolated, symmetric dimer. Each enabled written input belongs to one such dimer, while every blocked input in this basis is isolated from the destination. Consequently the same support holds for **every** isolated pulse duration. For an enabled dimer with diagonal $d$ and real off-diagonal $v$, the transfer amplitude is $-i\exp(-it_ed)\sin(t_ev)$. The code checks both orientations and charge signs at five appointments against the local exponentials. [Integer source and replay](passive_loops.py)

## Two laps determine every exactly restored code

Distinct contacts act on distinct head factors. Sharp arrivals fix the load word at each step. With $P$ the product of the blank projectors and $a=\prod_e b_e$, the ordered success operator on heads is

$$W=\bigotimes_{e\in C} B_e\otimes I_{\mathrm{rest}},\qquad
 WP=PW=aP,\qquad W^2=a^2P.$$

Tensoring with the identity on an arbitrary retained reference preserves these identities. Suppose $V$ encodes a code and one successful lap satisfies $WV=VT$ for an invertible code map $T$. This includes a unitary times a nonzero success amplitude. Since $P W^2=W^2$ and $W^2V=VT^2$, surjectivity of $T^2$ gives $PV=V$. It follows that $VT=aV$. If $a=0$, every two-lap amplitude is zero and there is no nonzero restored code. Otherwise the loop acts as the identity relative to the matched all-blank path. A bare odd product of negative Hamiltonian hops has a kinetic sign shared by that reference. [Lean proof of the code conclusion](PassiveLoop.lean)

On receiver ports 159, 164 and 165 the full integer three-contact operator has 27 nonzero entries on 1,331 head words for each charge sign. Its square has just the blank-to-blank entry. With the actual two blank heads and one live even head, the first-lap rank is one and the second-lap rank is zero. This includes every head superposition in that receiver preparation. A single return can nevertheless have nonzero overlap with the initial state. Discarding which-path records or projecting the heads again changes the experiment. [Reproducible receipt](passive_loops.json)

## Can the same comparison detect a loop phase?

On the all-blank loop a full native pulse at $t=\pi$ transfers a vacancy with amplitude $+1$ on each link. Supply the independent edge phases of the icosahedral magnetic fixture in just one interferometer arm. The matched ordinary route in the reference arm gives a relative phase $2\pi q/20$ around each of the 20 faces. For supplied flux $q=0,1,10$, the respective exact-model bright-port probabilities at the in-phase detector are $1$, $0.9755282581$, and $0$. The quadrature detector gives $0.5$, $0.6545084972$, and $0.5$. These calculated probabilities have unit route acceptance. Their phase comes from the supplied magnetic connection. [Magnetic fixture](../../flux_doublet/measure.py) · [Control receipt](passive_loops.json)

Free evolution can superpose paths with backtracking and repeated contacts. An active protocol can restore a direction bit or act jointly on the charge and detector. A partial trace can preserve logical information even when the full head word changes. The theorem concerns the single prescribed route with its complete material state restored. Interference between freely propagated paths, spin under a physical spatial rotation, and transport by a multiparticle excitation require their respective instruments.

From the repository root, use one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.flavor_search.passive_loops
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/flavor_search/test_passive_loops.py
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/flavor_search/PassiveLoop.lean)
```
