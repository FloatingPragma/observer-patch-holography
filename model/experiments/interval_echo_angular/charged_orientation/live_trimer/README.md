# Can a freely moving charged composite carry a quantum state?

**Two vacancies and one excess carry orbital coherence from one cell to another.** The source encodes a qubit in which of two ports holds the excess. Three written memories participate in the flight and retain their physical direction coins. The receiver ignores all three. At the strongest sampled exposure, 131072 native units, **42.1153%** of maximally mixed launches place all three defects in the receiving code. The two sharp logical inputs reach it with probabilities **84.1372%** and **0.09334%**. A six-state process measurement detects transmitted entanglement after retaining all misses. This is a state-dependent quantum channel for a material pattern of net charge $-1$. [Protocol](PROTOCOL.md) · [Full simulated counts](results.json) · [Independent native evolution](test_live_trimer.py)

The two triangles occupy ports $(1,0,5)$ and $(82,74,78)$ on the twenty-cell federation. Their gateways are the endpoints of the actual seam $1$–$82$. Start one vacancy at the source gateway, one at a source tip, and the excess at the other tip. The receiving triangle starts uniformly filled. On each triangle the two tip words differ by an ordinary exchange of loads 0 and 2 on one blank edge. An isolated half pulse on that edge, and a private degree-reference phase on one tip, prepare and read the three Pauli bases. The corresponding local matrix elements are checked on the full 99-state seam. At each unit-contrast head the native echo on the even or odd superposition of (1,0) and (0,1), followed by a zero-coin read, supplies the selected parity with probability one half. The flight banks begin after preparing the heads and loads. Addressed load initialization, timing and controlled phase comparison are source and detector services. [Source, geometry and sector](measure.py) · [Local readout](detectors.py)

During flight every registered contact operates. Native zero-contrast writes on the exterior confine the load changes to these two triangles. The active interior contacts coexist with the original degree references. The matched receiver uses one unit-contrast edge on the source, one on the receiver, and one across the boundary. Each of their direction coins starts at zero and can change under free motion. The bank includes all final coin words. All three defects arriving together is a terminal coincidence. The experiment does not resolve whether a bound trimer stayed intact throughout the flight or measure a mass shell.

## Which degree of freedom travels?

| Measurement | Observer result |
| --- | --- |
| Matched receiver at exposure 131072 | The smallest eigenvalue of the reconstructed partial transpose is $-0.0122564$, with simultaneous 99% radius 0.001376. The native channel predicts $-0.0122629$. A negative value certifies entanglement with the input reference. |
| Receiver with two equivalent tips | Port interchange fixes the initial state and the Hamiltonian. Every accepted destination orbital is the equal superposition of the two tips, giving a rank-one output for arbitrary source rail superpositions. |
| Zero-contrast intercell seam | Complete arrival has exactly zero probability at every exposure. |
| Coherent input on the old launch head | Opposite phases between its two contrast parities give identical destination statistics even when all three defects arrive. At exposure 131072 the two banks record 28,225,270 and 28,226,646 arrivals of the same sharp tip pattern per $2^{26}$ launches. |
| Body turn of the arriving pattern | Two physical blank-edge exchanges cycle all three material positions on the receiving face. Three repetitions return the complete pattern, giving 28,261,399 bright, zero dark and 38,847,465 incomplete arrivals per $2^{26}$ mixed launches in the phase comparison. |

The reader reconstructs the unnormalized Choi matrix of the destination channel by preparing the six sharp and coherent tip states. It includes both signed outcomes of each three-axis read and the failed arrival in every bank. The sampling radius covers all settings and all sampled times, including the selected long exposure. A separate complex matrix exponential followed by seventeen squarings tests this long-time prediction on the full connected component. The threefold body turn cycles the loaded triangle on fixed heads. Its positive full-turn phase characterizes that material rotation. A spin assignment requires physical rotations of the complete moving dressing and exchange of two independently prepared copies.

## What do the native memories permit at the bottom of the spectrum?

On **every finite native board**, choose each seam's blank or written state. For a written seam specify its contrast magnitude $m$ and the parity $s$ that exchanges the labels $d$ and $-d$. Both are constants of free evolution. The direction coin and the material loads can evolve. In that basis the complete finite Hamiltonian has nonpositive off-diagonal entries. For every connected component of its configuration graph, the lowest eigenvalue is simple and its eigenvector can be chosen strictly positive. Consequently, two rotations that preserve this ground line have positive order commutator; the line cannot carry an intrinsic two-state spinor. This covers live written memories and arbitrary ternary material compositions. Degeneracies across distinct components and excited-state multiplets require separate tests. [Integer reduction and Perron–Frobenius proof](PROOF.md)

From the repository root, replay the native records and their independent checks with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.experiments.interval_echo_angular.charged_orientation.live_trimer.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/live_trimer
```
