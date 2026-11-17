# Portable complete source flow

The original source packet is preserved at `workspace/openwave/shared/compact_logarithmic_flow/`, with its exact producer manifest and complete source dependencies at their original relative paths. The outer manifest authenticates this portable layout; it is not a replacement for the original source seal.

Run `python3 workspace/openwave/shared/compact_logarithmic_flow/headline.py --attempt stage2_002` for the bounded independent source/record check. It does not replay millions of processor pulses. The full isolated producer replay, with all scientific fields and ten artifact bytes matching, is recorded in the package-level `reviews/flow/`. No download or sibling checkout is needed. Python3 and NumPy are required.

For a new full execution, run the inner `run_flow.py --record --attempt NEW_UNIQUE_NAME`. The complete retained outcome generation takes about17–24minutes on the recorded machines.
