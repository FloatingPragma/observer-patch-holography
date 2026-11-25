# Fixed-priority repair with intervening native flight

On the 20-cell Interval Echo sphere, a sharp double at port 0 writes seam (0,1) if the installed priority sweep runs immediately. Leave the same board open for one native time unit before that sweep and a different seam is written with probability **0.72496213874**. If flight follows the sweep, the address stays (0,1). At both endpoints every later native flight and strict repair preserves the occupied address. The comparison uses the same clock exposure and priority assignments. Scheduling changes a permanent local record. A random guarded repair still terminates almost surely throughout this sector. [Integer checks and complete address counts](flight_confluence.json)

## The instruments and the common repair relation

Take the complete native federation, with each port initially at zero except a single double and every collar blank. Let \(|x\rangle\) denote a double at port \(x\). Sharps, mixtures, coherent superpositions and entanglement with an untouched reference are allowed in this domain. Use the installed priority at each port. \(F(x)\) is its highest-priority incident seam. Distinct incident priorities make this unique. The repaired endpoint from a sharp \(|x\rangle\) has two unit loads, one on each endpoint of \(F(x)\), and a written magnitude-two collar there.

The observer asks which physical seam is occupied. The two signs at one address give the same answer to this question. Let
\[
 \Pi_e=\sum_{x:F(x)=e}|x\rangle\langle x|.
\]
The addressed sweep on the doubled sector measures this family of projectors. Each strict ready check on blank inputs asks whether the double is at one of its endpoints. The legal priority order writes the first ready address. Two columns with the same winning seam have opposite endpoint origins and orthogonal written sign states. Their probabilities add when reading the occupied address. The same proof gives the reduced joint output with an arbitrary spectator. [Coherent sweep](coherent_sweep.py) · [Native write and sign columns](strict_pair.py)

On the blank doubled sector a native exchange moves the double as a whole. On every incident seam the squared-contrast energy is \(4\kappa\). The reference contacts give \(4\kappa(6-d_x)\), so the total contrast energy on each sharp doubled word is exactly \(24\kappa\), including at degree-six ports. The complete Hamiltonian is \(L/2+24\kappa I\), with \(L\) the unweighted native graph Laplacian. The stored two-load convention subtracts \(12\kappa I\). The probability calculation removes the remaining scalar too. On the local write columns, label reversal interchanges the two input orientations and squared contrast acts by four. The implementation verifies both with integer entries over the entire native written space, as well as the resulting 99-state numerical operator. [Degree-reference law](../../candidate_laws/contrast_memory_spherical_v1/SPEC.md) · [Native doubled Hamiltonian](../../candidate_laws/contrast_memory_spherical_v1/sectors.py)

After a commit all loads are binary. Unmarked seams preserve their blank binary sector during flight. At the one marked seam the two magnitude-two parity states with every binary endpoint input are invariant under its native terms, with no transition to blank. The individual sign may evolve. On binary inputs every later strict guard is quiet or sees an occupied collar. The mask consisting of one occupied address therefore persists under any continuation by these flights and guards, irrespective of subsequent schedule. Deliberate erasure or undoing of committed records would define a different relation. [Exact strict-sector columns](strict_pair.py) · [Sharp consistency and descent](ASYNCHRONOUS.md)

## The nonjoining fork

Compare complete sequences consisting of one priority sweep and one free interval \(t\), in opposite orders. Put \(U_t=e^{-itL/2}\). For launch \(|x\rangle\), the total variation distance of the two final occupied-address distributions is exactly
\[
 d_x(t)=1-\|\Pi_{F(x)} U_t|x\rangle\|^2 .
\]
The immediate sweep always writes \(F(x)\). A branch in the other sequence with a different chosen seam cannot be joined to this branch by any continuation. The same argument on unconditioned outputs uses their distinct, conserved address distributions.

A lower bound holds at every port on every such degree-five/six federation. The fiber of \(F(x)\) contains at most the two endpoints of that seam. At least \(m_x\ge d_x^{\rm port}-1\ge4\) neighbors have a different winner. For \(G=L/2-3I\) one has \(\|G\|\le3\), using the degree bound on the positive Laplacian. With \(Q=I-\Pi_{F(x)}\),
\(\|QG|x\rangle\|=\sqrt{m_x}/2\). The imaginary part of \(Q e^{-itG}|x\rangle\) has norm \(\|Q\sin(tG)|x\rangle\|\), since the basis and \(G\) are real. The spectral theorem gives \(\|\sin(tG)-tG\|\le(3|t|)^3/6\). Hence
\[
 d_x(t)\ \ge\ t^2\left(\max\{0,\sqrt{m_x}/2-9t^2/2\}\right)^2,\qquad
 d_x(1/4)\ge529/16384 .
\]
This bound survives increasing the number of cells. At a fixed time the result for a specific instrument also follows from the integer powers of \(L-6I\). The degree-22 Taylor approximation has operator error at most \((3t)^{23}/23!\). Only the one or two return amplitudes in the original fiber enter; projecting the approximation changes its squared norm by at most \(2\epsilon+\epsilon^2\). The table gives conservative enclosures, rounded outwards.

| Intact sphere | Launch | Different address at \(t=1/4\) | Different address at \(t=1\) |
|---|---|---:|---:|
| 20 cells | Port 0 | [0.075411614236, 0.075411614237] | [0.724962138737, 0.724962138752] |
| 20 cells | Range across all 240 initial ports | [0.060541286522, 0.075451397902] | [0.611751754288, 0.732000499924] |
| 80 cells | Port 0 | [0.075411614230, 0.075411614231] | [0.724961664764, 0.724961664780] |
| 320 cells | Port 0 | [0.075411614230, 0.075411614231] | [0.724961664764, 0.724961664780] |

On the 20-cell port-0 launch, separate simulated acquisitions gave 79,464 and 760,355 different addresses in \(2^{20}\) trials at the two times. An even simpler fork uses the isolated native contact (0,1) for a \(\pi\) pulse. Before repair it moves the double exactly from 0 to 1. The fixed sweep then marks (1,82). After repair the same pulse leaves the occupied (0,1) record intact. The total variation distance is one. Removing exchange during flight while retaining contrast, the degree completion, routing and the strict instrument makes the difference zero. Removing only the contrast coefficient leaves the blank double walk unchanged up to its common phase. [Exact controls and samples](flight_confluence.py)

## Which schedules can keep the address result?

For a single pre-sweep unitary \(U\), equality of the address distributions for all initial states is equivalent to \(U^\dagger\Pi_e U=\Pi_e\) for every \(e\), since these are the effects measured by the two schedules. This equality also preserves the joint classical address and arbitrary untouched spectator. For every real time in an interval this condition is equivalent to the generator commuting with each fiber projector. It gives an exact sufficient domain for controlled flight inside an address fiber.

More generally, suppose a decoder of the marked address outputs \(o\) with probability \(f_o(e)\). For a sharp doubled launch \(x\), the reordering condition for small \(t\) is
\[
 f_o(F(x))=\sum_y |(U_t)_{yx}|^2 f_o(F(y)).
\]
In the order-\(t^2\) coefficient each neighboring vertex has jump weight \(1/4\). Probability conservation supplies the opposite diagonal coefficient. Thus \(g_o(x)=f_o(F(x))\) is graph-harmonic. The maximum principle on a finite connected graph makes every \(g_o\) constant. Conversely any address decoder with identical output laws on all reachable winning seams survives this comparison. No distinction between reachable addresses, including a randomized one, can be retained by equating these schedules on this quotient. This statement concerns information read from the winning address. Different conserved observables can define other quotients.

## Finite descent and almost-sure termination with flight

On the declared sector, a rank of one for a single unrepaired double and zero for the binary state with its unique occupied seam is preserved by all flights and quiet checks. An accepted strict write decreases it by exactly one, also decreasing the squared loads by two. The two domains are distinguishable by the local ready and occupied flags in a complete survey. A complete, uninterrupted legal sweep catches any unrepaired input with probability one. This gives termination and completeness for that survey with an arbitrary preceding interval of native flight.

There is also a scheduler that allows flight between successive local attempts. At a single atomic attempt an independent fair choice offers one of the \(M\) seams. Its local controller tests its own readiness and the higher-priority competing requests, with the flight parked for that atomic guard. On the doubled input space all readiness projectors are diagonal. The effect of an accepted offer at \(e\) is exactly \(\Pi_e\), as each competing request shares the double only if it outranks \(e\) at its occupied endpoint. The code verifies all 630 guarded effects on the 20-cell sphere directly from the controller’s local comparison lists.

An arbitrary history-dependent flight preceding the fresh choice preserves the doubled sector. Conditional on each prior survival branch, with any coherent state \(\rho\), the next acceptance probability is
\[
  M^{-1}\sum_e \operatorname{Tr}(\Pi_e\rho)=1/M.
\]
The residual quiet branches stay in the doubled sector. Induction gives exact survival probability \((1-1/M)^n\) after \(n\) attempts, mean waiting time \(M\) attempts, and almost-sure termination. On the 630-seam sphere the three finite first-offer acquisitions at different source phases and dwells gave 1,617, 1,663 and 1,661 hits in \(2^{20}\) trials each. All branches either commit once into a persistent binary configuration or continue waiting. Termination of this scheduler does not enforce equality of the resulting physical addresses. [Guard construction, waiting law and outcomes](flight_confluence.json)

The protected coherent sweep therefore retains its existing schedule theorem. Enlarging its relation to allow the same free interval on either side of a commit is refuted on the persistent addressed quotient, even for sharp initial data with one double. The global finite normal-form criterion also asks for its claimed consistency predicate, descent and all forks on the declared physical domain. [Canonical FT-CONFLUENCE criterion](../../evaluation/EVIDENCE.md#ft-confluence)

Run `OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/test_flight_confluence.py` from the repository root.
