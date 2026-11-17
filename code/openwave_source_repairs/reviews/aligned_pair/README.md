# Full-pair spatial check

Run `python3 -B check.py /path/to/2026-11-16_codex_aligned_pair_spatial` with
Python 3 and NumPy. The default path selects the neighboring producer packet.
The script verifies the complete sealed producer before checking independent
geometry, source matrices, pump powers, all joint positions and error bounds.
It writes a fresh local `RESULT.json`; execute a copy if preserving this seal.
`REVIEW.md` records the scientific scope and same-group verdict.
