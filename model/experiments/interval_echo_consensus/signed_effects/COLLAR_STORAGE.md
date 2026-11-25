# Storing a material address bit during free motion

A pointer prepared in either orientation of the neutral dimer, (2,0) or (0,2), can leave its classical bit in the seam's private memory. The two stored values are separately conserved under free evolution of the entire federation, including exchanges of material on neighbouring seams. Four disjoint pointer seams can therefore retain a complete four-bit address. The encoding uses a protected pulse program. Decoding requires protected access and a sharp terminal read. [Native acquisition and all nine endpoint loads](collar_storage.py) · [Terminal records](collar_storage.json)

## Deposit and conservation

In local port order (4,5), use the degree-five reference of port 5 for a dwell of π/64, open the pointer contact for π/2, and repeat the reference dwell. On the ordered input pair (2,0), (0,2), this acts as the Hadamard matrix with common phase exp(−iπ/4). Apply the native echo with blank seam memory. Both outputs have loads (1,1). Their memory states are equal-weight superpositions of the two extreme contrasts with a symmetric coin. The relative sign between contrasts is positive for input bit zero and negative for bit one. The four stations use cells 1, 2, 4 and 5 of the 20-cell federation. Each pointer has a physical degree-five reference.

For an endpoint load word q, write P for the projector onto q=(1,1). The full native seam Hamiltonian acts on these two memory lines as 32(I−P) and 33(I−P), respectively. The echo sends the balanced endpoint columns into the blank sector of the conjugated generator, where their energies vanish. At the other eight loads, the extreme contrasts are outside the write image and have the stated reversal energies. Each memory line times the entire material space is an invariant subspace of a Hermitian matrix. Other seam terms act on different memories, even when they share a material port. Degree corrections act on material. The two bit projectors therefore commute with every term of the free board Hamiltonian. Products of the projectors on the four disjoint memories conserve all joint address probabilities for every free duration. [Seam reduction](../../interval_echo_metric_response/measure.py) · [Test with a freely exchanging neighbour](test_collar_storage.py)

This statement concerns classical record probabilities. For example, prepare equal coherent amplitudes on the two memory bits and on the two material words (1,1) and (2,0), which have the same total load. One π exposure of the written contact gives a reduced memory purity of 1/2. Recovering the classical address does not establish storage of a coherent address superposition without environmental entanglement.

## Retrieval and controls

Supply (1,1) at the written endpoints for the inverse echo and inverse basis pulses. Two blank neighbouring contacts permit an exchange with a fresh pair of unit loads, if supplied at retrieval. The exchanges preserve whatever material word they displace in the scratch registers. A sharp read on the decoded material pair gives the original bit. Acquisition tests the procedure by exchanging each of the nine scratch inputs onto a deposited record, exposing the written contact and exchanging the unit loads back. The terminal measurement resolves all nine material words with blank memory and one occupied-memory rejection. The fixed inputs, isolated write and read services, and sharp measurement are apparatus.

Each of the 216 banks retains 4,096 simulated starts. With π exposure and either input bit the native encoded bit is recovered in all 73,728 starts across the nine disturbance settings. A direct echo without the basis pulses exchanges the bit with probability one on each of the eight disturbed endpoint words other than (1,1). Halving the exchange term preserves the alphabet, incidence and contrast interaction. With the same basis pulses and calibrated full swaps and inverse decoder, that law gives error probabilities 1/4 after π and 1/2 after 2π for each of those eight words. The balanced word recovers the bit in every arm. Coherent preparation, free persistence of the classical record and the eventual measurement have distinct physical requirements.

Run from the RER root with the repository's numerical environment.

~~~sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.collar_storage
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/test_collar_storage.py
~~~
