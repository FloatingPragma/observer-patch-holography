# How many directions can an Echo observer use for the Einstein test?

A single echosahedral observer recovered the same rank-three spatial frame in **32 of 32 simulated count tables** at 65,536 preparations per source port. The source and detector labels were permuted separately for each table. The inverse received the count table and the uniform, reciprocal, weak-exchange model class. At 128 preparations per port it abstained on all 32 tables. The high-budget measurements determine a usable angular design before a stress or curvature sensor is chosen. [Measurement protocol](PROTOCOL.md) · [Archived design](results.json) · [Native spatial instrument](../../evaluation/patch_spatial_readout.py)

| Available directional readings | Subsets spanning nine tensor combinations | Subsets spanning eight |
| --- | ---: | ---: |
| 9 | 160 | 60 |
| 10 | 60 | 6 |
| 11 | 12 | 0 |
| 12 | 1 | 0 |

In the native port chart, **0, 1, 2, 3, 4, 5, 6, 8, 9** are one of the best conditioned nine-port choices. We used the native chart only to name one example among tied solutions after reconstructing the frame. The smallest singular value is 0.563 with unit spatial directions and the Euclidean Frobenius norm on tensor coefficients. The least sensitive full-rank nine-port choice gives 0.267. Averaging each direction with its opposite leaves only six independent readings. Losing the orientation suppresses all three time–space components. Keeping eleven directed readings guarantees full rank if any single receiver is unavailable. [All subset counts and sensitivity](results.json)

The calculation uses the normalized inner products of the measured slow modes. If those directions are realized as k = (1,n) in one orthonormal spacetime frame, the linear functional on a symmetric tensor is the contraction with k kᵀ. The Frobenius inner product of two such functionals is (1+n·m)². The rank of their Gram matrix gives their independent count without a choice of spatial axes. A second calculation gives the same ranks using explicit symmetric tensor coordinates. All null contractions vanish on the metric itself. Independent matter conservation and a geometric identity enter the Einstein paper's separate step for recovering the trace term. [Tensor closure and small-ball premises](../../../paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.tex)

G5 still needs measured null events in a regional frame and remainder bounds on successively smaller balls. This experiment supplies finite sampling and receiver-selection evidence for its directional clause. Gate exposure is supplied, the observation is on one cell, and both stress and curvature readings remain to be implemented in this frame. The earlier source test found a 61.64% mismatch when the gauge-current readings of a phase-conjugate pair were identified directly with null stress. [Scope of the current control](../../candidate_laws/interval_echo/prior_calibration.json)

From the RER root run `python -m model.experiments.interval_echo_angular` to reproduce the tables and all subset comparisons.
