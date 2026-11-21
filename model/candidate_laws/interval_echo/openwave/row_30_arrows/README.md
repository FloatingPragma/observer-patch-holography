# Transport with a record on every seam

A single unit released at port zero can reach an observer in another cell while every physical seam retains a contrast record. The outcome depends strongly on the coin written on each seam before launch. At exposure one, a patch in cell 6 registered **34,107 arrivals in 262,144 simulated trials** under the second coin ordering, **595** under the first, and **zero** in the parked control. The physical sphere, 630 active seams, launch port and detector were identical. The initial written signs were in an even superposition on every seam. Each cell's blank degree references permitted an ordinary local contrast read after flight. These clocks are supplied exposure appointments. [Protocol](PROTOCOL.md) · [All ten detector outcomes at each appointment](results.json)

| Prepared coin order | Reachable microscopic states | Distinct ports | Distinct cells | Predicted patch detection at exposure 1 |
| --- | ---: | ---: | ---: | ---: |
| Ascending port IDs | 32 | 12 | 3 | 0.2346% |
| Ascending cells, reversed local labels after port zero | 2,730 | 54 | 7 | 12.9426% |
| Descending port IDs | 1 | 1 | 1 | 0% |

The receiver uses the nine degree-five ports of the cell across native gateway (1,82). At each port the observer applies the registered echo to its blank zero-load reference, copies the absolute-contrast-one outcome, and reverses the echo. One unit gives at most one reported address. Every trial, including quiet, contributes to the denominator. At exposure 1/4 both moving preparations also crossed the detection threshold. Every reachable configuration of all three runs was retained. [Native column checks, path enumeration and propagation](experiment.py)

### Why the background changes transport

For endpoints (a<b), a coin of one places an arrow head at (a) and zero places it at (b). In the sign-even, absolute-contrast-one sector the native edge term equals κ times the identity minus κ times the projector onto the equal superposition of (unit at (a), coin one) and (unit at (b), coin zero). A hop leaves an arrow head and flips that coin. Thus the difference between the port occupation and its number of incident arrow heads is conserved at each port. The experiment verifies the six vacuum and one-unit columns of the full 99-state operator including leakage, and includes every blank degree-reference term in the subsequent flight. [Registered interaction](../../../contrast_memory_spherical_v1/SPEC.md)

For these acyclic preparations the flipped edges always form one increasing path from the source. The unit can extend this path at the occupied port or retrace its final edge. Distinct routes to the same port remain distinguishable in the physical coins. Counting paths in the initial ordering gives the total state count independently of the explicit evolution. In particular every change of endpoint changes the coin word. Tracing the coin registers therefore gives **zero off-diagonal coherence between distinct ports for all times** on these sharp launches. The joint state of the unit and coins evolves coherently and can retrace a path.

### A patch can read the joint phase

The receiver applied a quarter contact or its conjugate on one seam after otherwise identical free flights. The native contact mixes the unit at one endpoint with the unit at the other and flips the seam's coin. The patch then read the target port through its degree reference. At exposure 1/4 on seam (82,78) under the second ordering, the conjugate readings gave **57,693 and 16,972 clicks** in 1,048,576 simulated trials per setting. Their probability difference measures the imaginary part of the joint load-and-coin overlap. A random native reference phase before the mixer erased that overlap without altering the preceding port occupation distribution. The observed contrast was **+0.03883** in the coherent arm and **−0.00024** after dephasing. On the same seam the coherent contrast changed sign at exposure one, reaching **−0.02817**. [Phase protocol](PHASE_PROTOCOL.md) · [All counts and operator checks](phase.json)

Eight of twelve link and time combinations crossed the 0.004233 family-wise margin in the second ordering. No setting crossed it in the ascending or parked orientations, or in any of the 36 dephased controls. All three arms used the same receiver window and timing convention. The readout resolves a local current quadrature across an external receiver's seam. Obtaining a matter-wave dispersion requires several preparations with independently measured spatial wavelengths and a frequency measurement on that same branch. [Joint phase reader](phase.py)

From the RER repository root, replay the complete receipt with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_arrows.experiment
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_arrows.phase
```
