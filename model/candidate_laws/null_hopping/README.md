# Bare hopping

This control runs the contrast-memory seam at zero agreement energy on the single cell and the 20-, 80- and 320-cell federations. The exchange scale, seam register, carrier and port action are those of the reference law. Its full local Hamiltonian is the reference echo conjugate of the private exchange term at κ = 0. The degree-reference energy also vanishes. The comparison specifies free evolution and a first seam reading at a supplied dwell. [Reference specification](../contrast_memory_spherical_v1/SPEC.md) · [Registered operator](candidate.py)

## Sharp and coherent one-load preparations

With blank seam memories the local operator is half the identity minus half the port swap. Summing over seams gives the graph Laplacian divided by two in the one-load sector. The agreement and degree terms of the reference add exactly 6κ times the identity on that sector. The reference and the control therefore give identical one-load probabilities for every initial density matrix, dwell and common measurement. The same equality holds with an idle ancilla. Coherently comparing different load sectors or conditional evolution on two interferometer arms can expose their relative phase.

On one cell the levels of L/2 are 0, (5−√5)/2, 3 and (5+√5)/2, with multiplicities 1, 3, 5 and 3. The identity map intertwines the port basis with L/2. Its first nonconstant eigenspace is the image of G/4, where G is the Gram matrix of the twelve unit port directions. A triplet in this measurement therefore cannot distinguish κ = 8 from κ = 0. The executable tests check the matrix identity on all four carriers and the slow Gram subspace on the cell. [Operator and record controls](test_control.py)

## Two-load discriminator

Put one unit of load at each of ports 0 and 1 on a single cell, with every other port empty and all seams blank. After one unit of exchange time, each seam reads whether both ends hold a unit. The complete first reading has 30 possible seam addresses and one outcome for a separated pair. Each preparation has unit weight. The table gives calculated Born probabilities.

| Agreement energy κ | Pair on any seam | Separated pair |
| ---: | ---: | ---: |
| 0 | 0.540513658207 | 0.459486341793 |
| 8 | 0.992611885462 | 0.007388114538 |

On this sector the agreement energy, up to a constant, is −2κ times the projector onto adjacent pairs. A hop can change adjacency, so this term affects the free record distribution. The comparison distinguishes the laws on a prepared pair and a finite dwell. A stable particle identification requires additional preparations, lifetimes and transport records. [All outcome probabilities](results.json) · [Free release and seam reader](replay.py)

Run `python -m model run --law null_hopping --check` from the repository root. The independent regression constructs all pair swaps on the 66 unordered pairs directly and checks every first-read probability using dense matrix exponentiation.
