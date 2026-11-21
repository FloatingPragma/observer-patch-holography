# Strict repair and retained phase on an intercell wedge

A complementary service settles sharp discrepancies across the entire sphere with local priority checks and asynchronous requests. The proof covers every sharp load and signed collar configuration in that service. The measurements include every branch on a five-port path and nine complete-screen preparations. [Asynchronous sharp consensus](ASYNCHRONOUS.md)

A two-cell observer can settle two competing seam repairs at a common endpoint and retain a readable relative phase within the winning address. The controller has a stored priority and a service window containing both requests. It probes the preferred seam and commits immediately if it is ready. Only a quiet probe opens the competing seam. Reading both readiness flags in advance erases a phase that this priority service preserves. [Phase experiment and detector](phase.py) · [Phase probabilities](phase.json)

## Which operation counts as a repair

On a blank collar the native echo maps both (2,0) and (1,1) to a state with endpoint loads (1,1). For (2,0), the squared loads fall from 4 to 2. For (1,1), they start and end at 2. The pair pointer fires with probability one in both cases. Extending the original sharp-double acceptance test to an initially balanced adjacent pair would therefore count a transaction without strict descent. The written collar already distinguishes the two inputs by their contrasts. [Native echo](../../candidate_laws/contrast_memory_spherical_v1/law.py) · [Full input table and controls](results.json)

The strict instrument requires a blank collar and tests the **absolute** written contrast of two along with the pair pointer. Its positive subspace contains both signed contrasts. Keeping those signs coherent permits the local phase read. On all six initial two-unit basis states, each accepted write lowers the squared loads by exactly two. Every resulting state is stable under another strict attempt at either address.

## Order and phase controls

One internal and one external seam share a port. The test uses the first intercell bridge on each common federation and a disjoint physical priority rail. The three embeddings have the same local tensor operation. At fixed priority the complete endpoint channels for reversed requests differ by less than 10⁻¹⁰ in the Frobenius norm of their Choi matrices. Equality on all six input columns and all their coherences covers arbitrary superpositions and mixtures in the declared sector. The classical accepted address is part of the output. Under first-ready arbitration a sharp double at the shared port gives disjoint endpoints in the two orders, at trace distance one. [Protocol](PROTOCOL.md) · [Channels and material layouts](results.json)

With priority A, prepare one double in a superposition of the left and middle ports. Both alternatives should be accepted on A. Reading B’s readiness reveals which alternative occurred. Equal and opposite source phases then give identical material states even if both classical readiness results are kept. In the priority window the observer leaves B alone on a positive A. Undoing the winning echo, applying a half contact and reading its first rail against a vacancy distinguishes relative phases +i and −i with probabilities 0 and 1. The two-probe control gives 1/2 and 1/2. Reflected preparations at priority B give the same boundaries. These are exact finite-instrument probabilities to numerical tolerance. [Follow-up protocol](FOLLOWUP.md) · [Readout](phase.json)

The active observer occupies two cells. Other contacts are idle during service. The controller buffers both requests, with the stored priority, blank collars, finite window and phase-reference apparatus supplied. The endpoint excludes chronological attempt logs. The calculation certifies one intercell incidence, one total-load sector and a prepared priority. General convergence, arbitrary arrival with missing requests, and confluence of the clock’s transactions require separate arguments.

## Reproduce

From the RER repository root, run:

```sh
python -m model.experiments.interval_echo_consensus
python -m model.experiments.interval_echo_consensus.phase
python -m pytest -q model/experiments/interval_echo_consensus
```
