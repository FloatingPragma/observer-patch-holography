# A two-load defect on the 320-cell screen

The finite patch in cell 13 reported an adjacent pair on one of its 30 internal seams with probability **0.99098** after one native exposure on the 320-cell screen. Turning off the contrast interaction reduced that probability to **0.45000** at the same exposure. The native result exceeds the fixed 0.80 gate and the exchange-only control stays below its 0.65 ceiling even with the full numerical allowance. Both launched from the same first gateway-touching seam, (156,157). [Frozen protocol](PROTOCOL.md) · [All 31 probabilities and sampled counts per arm](results.json)

| Registered screen | Native patch probability | Exchange only | Native hits / 16,384 | Control hits / 16,384 |
|---|---:|---:|---:|---:|
| 20 cells | 0.990982505 | 0.450001676 | 16,227 | 7,375 |
| 320 cells | 0.990982502 | 0.450001733 | 16,241 | 7,449 |

The observer interrogates the load pair through each physical seam’s native echo, copies the agreed-pair bit and undoes the echo. With exactly two binary loads, a positive reading selects a unique address. Quiet includes every configuration outside the 30 measured local seam states. The 320-cell screen has 3,840 physical ports, 10,080 physical seams and 7,370,880 binary pair configurations. All original internal, external and degree terms act during flight. No criterion concerning long-time binding, charge, spin or electron identity was tested by this single-exposure comparison. [Local reader](../row_04/__main__.py) · [Other exposure and perturbation measurements](../row_04/README.md)

## Bound on the full-screen answer

The radius-six calculation retains 62 cells with 744 ports and 276,396 pair configurations. Its diagonal includes every contribution from the original full-screen degrees. To verify the boundary implementation, a strict radius-two compression on the smaller 20-cell screen was compared entry by entry to the principal submatrix of each complete Hamiltonian. Both discrepancies were exactly zero. On the full 20-cell screen a sparse exponential separately reproduced all 30 seam probabilities within 10⁻¹² of the Chebyshev recurrence for both interactions. [Hamiltonian and verification](experiment.py)

Split the generator into terms preserving each occupied cell and off-diagonal external hops V. Each occupied port meets at most one external seam. Every binary-pair row has at most two external off-diagonal entries of magnitude one half. Symmetry and the row-sum bound give ‖V‖ ≤ 1 for either interaction. In the interaction picture each external hop can increase one occupied cell's distance from cell 13 by at most one. Starting with both loads in cell 13, the first six orders of the full and projected Dyson series therefore agree. At exposure one, each tail has norm at most e/7!, giving state-norm error at most 2e/7!. For any projection, the probability error is at most twice the state error. The archived conservative floating-point allowance adds 10⁻⁹, including the Chebyshev truncation (below 10⁻¹²). The resulting 320-cell full-screen intervals are **[0.98882513, 0.99313987]** for native repair and **[0.44784437, 0.45215911]** for exchange only. The sampled frequencies and their binomial counting intervals are recorded separately. [Results](results.json)

Run the following from the reverse-engineering-reality checkout to reproduce and compare the entire receipt.

```sh
OPENBLAS_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_large.experiment
```
