# What can the free end of a receiving wing measure?

Comparing two local current readings on a three-dimer receiving wing resolves all three previously measured one-quantum excitation energies. An observer gets these frequencies even on the collision sham. When the observer first identifies a single arriving quantum and waits in silence, the signed terminal current has zero expectation. The same material supports spectroscopy and supplies a direct test of whether the detector starts the oscillation. **OW-03 remains partial.** [Frozen protocol](PROTOCOL.md) · [Integer records](counts.npz) · [Decoded results](results.json) · [Official criterion](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L437)

The source uses the published OW-08 collision on a guarded four-dimer ring feeding two three-dimer wings. At exposure 32 the observer records the phase of an ordered ring, closes its exchange shutters and leaves its measured load in contact with the inner wing. An unordered ring is recorded as a miss. The left wing's supplied energy-sector analyzer reports “one” if the result lies in the span of its three previously calibrated single-quantum levels. The read preserves superpositions within that span. The far end has a physical two-load dimer, with two possible signed load orderings (0,2) and (2,0). Isolate it, rotate by a quarter exchange period, then measure the signed ordering to read its imaginary exchange coherence, Y. After a starting read, undo that basis rotation before reopening the wing. The second wing is unread. These operations use timers, shutters, read pointers and a supplied spectrometer. [Original local spectrum](../../row_08_annihilation/daughter_spectrum/README.md) · [Geometry and mode identification](../../row_08_annihilation/daughter_modes/README.md)

| Observer record | Native contrast 8 | Contrast 4 during waiting only |
| --- | ---: | ---: |
| Clock lines after a starting read on one-quantum arrivals | 31.80078, 32.51538, 33.21468 | 15.80876, 16.53027, 17.22164 |
| Same read on all ordered sham records | 31.80076, 32.51538, 33.21463 | 15.80871, 16.53027, 17.22163 |
| Held-out RMS error with three lines, then with the strongest individual candidate | 0.0100, 0.2508 | 0.0099, 0.2511 |
| Largest absolute terminal mean current after a silent one-quantum wait | 0.0381 | 0.0511 |
| Joint 99% sampling radius on these terminal means | 0.0593 | 0.0592 |

Frequencies are angular frequencies in inverse model exchange-time units with ℏ = 1. Timing is supplied. Each setting and appointment uses 262,144 fresh launches. Every ring miss, both ordered phases, both energy flags and every requested sign are in the counts. The native silent panel admitted 3.0745% of launches. On separate copies, after a starting current read only **4,037 of 8,115** one-quantum events retained that flag on an immediate reread. The exact survival probability is one half. With silent native evolution the same sector is preserved. Halving contrast only after the native herald preserves the preparation; the changed generator need not preserve the native energy sectors. Discovery uses even appointments, with odd appointments reserved for prediction. The frequencies agree with the previous energy read to 0.00005 on the native wall. There is no zero-wave-number mode in this finite wing. Obtaining independent inertia and an autonomous rest oscillation remains an OW-03 obligation. [Proof and reference conditions](PROOF.md)

From the RER root, decode the archive, or add `--reacquire` for a fresh same-environment simulation. Replay checks the original local energy calibration, all count totals, held appointments, the changed law, and seeded positive and zero signals.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -O -W error -m model.candidate_laws.interval_echo.openwave.row_03_clock.free_end.measure
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_03_clock/free_end
cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_03_clock/free_end/FreeEnd.lean
```
