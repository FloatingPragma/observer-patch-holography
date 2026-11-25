# The area response across two source preparations

A balanced material return has a nonlinear entropy slope as the strength of pairing varies. A single unflagged input has an affine density matrix, so the expected value of every fixed area apparatus is affine in that strength. The two slopes cannot cancel on an interval. Repeating the source supplies another resource. A bounded readout requires at least order \(1/\ell^2\) independent inputs for a fourth-order cancellation, and the product of two majority bits attains that scaling for the measured scalar with a chosen gain. [Proof](PROOF.md) · [Native return and collective readout](response.py) · [Kernel-checked endpoint bound](../../../../Lean/ObserverPatchHolography/EinsteinBranch/AreaReadoutSpread.lean)

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -m matched_einstein_limit.source_implication.area_resources.response
OPENBLAS_NUM_THREADS=1 python -m pytest -q code/matched_einstein_limit/source_implication/area_resources
(cd Lean && lake build ObserverPatchHolography.EinsteinBranch.AreaReadoutSpread)
```
