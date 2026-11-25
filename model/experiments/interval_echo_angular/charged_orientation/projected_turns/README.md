# Projection can imitate a spinor order test

OpenWave records [OW-07 as untested](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L441); the full local criterion is open. A two-state spatial encoding of a native signed defect gives a negative perpendicular-half-turn order sign after projection back into its code. Each turn discards half the inputs. Its literal2π repeat has positive phase. This separates a projection-induced order sign from an intrinsic spinor rotation. [Protocol](PROTOCOL.md) · [All histories and controls](results.json) · [All-code bound](PROOF.md)

## A tight limit for coherent dressing

For any coherent encoding V, let P select its code and let R,S be commuting Hermitian full physical half-turns. The accepted maps A=V*RV and B=V*SV obey

||[A,B]|| ≤ 2 ||(I−P)RV|| ||(I−P)SV||.

If they are anticommuting unitaries scaled by input-independent successes pR,pS, this requires **pR+pS≤1**. Exact accepted spinor-like turns therefore cannot both succeed above one half in this class. The proof permits arbitrary finite material/memory dressing and untouched quantum references, while retaining the stated full physical actions. A genuinely noncommuting microscopic lift lies outside its hypotheses.

The constructed code saturates the limit: its accepted native turns are X/√2 and Z/√2. The physical operations are proper spatial π rotations of the actual12-port template with perpendicular axes, each implemented by18 guarded native π-seam contacts. Both full operators commute. Four physical ports carry the logical wavefunction. The coherent source isometry, ideal code/complement instrument and phase-reference interferometer are supplied; their native local compilation is unestablished.

## What the complete records say

Both material signs and all six Pauli-axis inputs are retained. After each filtered turn, the source may remain inside or leave the code. Rejected states continue into the second measurement and every history is recorded. A complete coherent-arm read compares opposite orders within each history.

| Protocol | Selected-branch phase | Probability per original start |
| --- | --- | ---: |
| Native opposite orders without filters | Positive | 1 |
| Native opposite orders with a filter after each turn | Negative | 1/4 |
| Same native half-turn repeated twice, filtered | Positive | 1/4 |
| Same native half-turn repeated four times, filtered | Positive | 1/16 |

The last two comparisons use an idle reference attenuated by an independent fair filter at each stage to match the selected amplitude. All its refusal histories are included. In the first vacancy/X-minus bank, the selected order gives **65,781 negative outputs** among262,144 starts;196,363 starts occupy the other histories. The literal full turn gives **65,444 positive outputs** with196,700 other-history starts. Thus the order sign alone would misidentify this apparatus.

The separately supplied spinor reference gives negative order and2π phases and positive4π phase through the same interference calculation. It installs the representation under test and earns no native credit. No two-core exchange is measured in this lane.

## Changed-law controls and calibration

Changing native contrast from8 to4 retains the same carrier, routing, ternary alphabet and symmetry. At these integer settings each π contact is exactly the same positive swap. The physical turn calibration errors are below3.6×10⁻¹⁵, and the projection-induced sign persists. The sign cannot identify the native contrast mechanism.

Halving exchange at the original contrast and pulse durations retains those same declared structures. The contact words then have distances1.99359 and1.99997 from the intended spatial half-turns; their full-turn calibration errors are1.98387 and1.88973. All order, full-turn and rejection records remain available. A phase from these failed spatial calibrations is not a tested rotation sign.

The result excludes obtaining a high-success intrinsic spinor solely by dressing or selecting these commuting physical rotations. Other memory-dependent dynamics, adiabatic paths, physical angular momentum, reconstructed spacetime rotations and exchange statistics require their own common-source tests.

## Replay

From the RER root:

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python -O -B -W error -m model.experiments.interval_echo_angular.charged_orientation.projected_turns.measure
python -B -m pytest -q model/experiments/interval_echo_angular/charged_orientation/projected_turns
~~~

The replay compares the complete archive without rewriting it. An explicit --output PATH creates a separate first record exclusively. Seven tests independently verify the physical axes, symbolic isometry, leakage identity, untouched-reference case, both rejection stages, literal full-turn control, optimized mismatch rejection and the actual optimized replay without altering its archive. Every bank contains262,144 synthetic starts. Raw probabilities and explicit complementary-bin rounding bounds are retained.
