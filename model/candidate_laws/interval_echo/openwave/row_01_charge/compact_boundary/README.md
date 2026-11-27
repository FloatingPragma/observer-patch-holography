# Can a remote material probe identify a compact charge?

An exterior receiver distinguished an enclosed vacancy, an excess and an undisturbed host when the gateway transmitted a mobile probe. Parking the gateway exchange gave exactly the same exterior outcome probabilities for all three sources. The native diagonal contact with the source cell was retained in both experiments. [Frozen protocol](PROTOCOL.md) · [Finite generator and exterior counters](measure.py) · [Complete simulated trial banks](counts.json)

An emitter in owner 0 prepared the existing compact spatial mode, with one load added or removed relative to the filled host. The full-board Hamiltonian preserves each preparation; its amplitude vanishes on every gateway. The probe was independently prepared in owner 6 on port 82. The aperture retains the two owners and their full-board diagonal terms. At each fixed appointment two sharp trit readers, calibrated against zero stock, reported the loads at ports 82 and 81. A simultaneous terminal census reported how many nonunit loads lay inside owner 0. There is no assignment of an original particle identity when the loads are identical. The total sample was 180 banks of 4,194,304 virtual starts. Every bank includes both readouts and every source outcome.

At native strength eight, with a vacancy probe and wait eight, both exterior ports read one on 2,738,967 starts with an enclosed vacancy, 2,146,182 with the uniform host and 3,509,004 with an enclosed excess. The two signed sources differ by 18.36 percentage points. The simultaneous 99% half-width for a two-bank difference is 0.274 percentage points. Removing bulk contrast preserves discrimination in the open geometry. These counts establish a material scattering response that also depends on exchange. The exterior detector reports no source-dependent flux in the parked geometry.

For the parked opening, the integer generator obeys an exact intertwining identity on the tensor product of the compact orbital with every exterior probe position. Its restriction is the exterior probe Hamiltonian plus a scalar depending on the source. The scalar contributes a common phase. The identity holds for either sign of source and probe and at both strengths. Exterior operations that preserve this encoding, including repeated number-preserving exterior reads, give identical record distributions for every finite history. An instrument able to populate the interior through exchange tests a different physical process. [Exact integer calculation](measure.py) · [Finite proof](QuietBoundary.lean) · [Integer support check at all times](quiet_support.py)

## What happens as the opening becomes weak?

A calibrated adjustable aperture multiplies just the bridge exchange by g, holding native contrast and exchange elsewhere fixed. At g=0 the first change of the amplitude transfers a defect between the owners. It has disjoint interior-number support from the undisturbed amplitude. Consequently the linear term of every exterior trit probability vanishes. This conclusion applies at each fixed finite wait.

Separate cohorts at g=1/16 and 1/8 tested a second-order expansion calculated about the parked apparatus. For each exterior event the absolute Taylor remainder is at most (g t)³/6. Indeed the bridge term has operator norm 1/2. Three derivatives of the evolved density in trace norm cost at most (2t‖V‖)³ by time-ordered unitary conjugation and the trace-norm commutator bound. Over 216 held probabilities, 48 moved beyond the sampling interval around the parked prediction. The quadratic prediction alone covered 204; including the analytic remainder covered 216. The largest central residual was 0.000173 and the largest predicted quadratic shift was 0.01285. This is a response to a deliberately opened material path. [Held protocol](WEAK_OPENING.md) · [Expansion and decoder](weak.py) · [Counts](weak_counts.json)

Run from the repository root with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_01_charge.compact_boundary
python -Bm "$p.measure"
python -Bm "$p.weak"
python -Bm pytest -q --pyargs "$p.test_boundary"
~~~

Each acquisition accepts `--acquire PATH` for a new simulated record. The primary preparation, clock, number aperture, two sharp readers and intermediate exchange settings are specified apparatus.
