# Interval Echo

The OpenWave evidence archive includes free particle persistence, a test of the same pair's colour response, and matter-wave experiments with blank and written seam memories. The research note gives the simulated patch records, controls and headless replay. [Rows 4 and 30](openwave/README.md) · [Evidence bundle](openwave/interval_echo_evidence.zip)

Interval Echo is the main candidate under test. Its elementary interaction is the κ = 8 contrast-memory seam. Its distinguishing rule is carried by an observer. The patch remembers how many of its own repair opportunities elapsed since the last verified handoff, then uses fresh local repairs to decide whether the closing mark advances its hand. Both patch owners read the two elapsed words on one signed seam at reunion. [Executable specification](SPEC.md) · [Reference interaction](../contrast_memory_spherical_v1/SPEC.md)

## What changes when the clock carries more memory

The three-bit clock accepts a closing primary unconditionally at interval length 8. Giving three physical fields one extra bit apiece extends the working range from 7 to 15. This costs one additional echosahedron at each travelling checkpoint. At length 16 the wider clock also overflows and releases the closing primary. The three-bit layout is the default for comparison with the registered radar runs.

| Completed handoffs | Three-bit traveller/home | Four-bit traveller/home | Predicted four-bit ratio |
| --- | ---: | ---: | ---: |
| Alternating 2 and 6 | 0.95535 ± 0.00190 | 0.95535 ± 0.00190 | 0.95602 |
| Alternating 6 and 2 | 0.95650 ± 0.00189 | 0.95650 ± 0.00189 | 0.95602 |
| Eight services per handoff | 0.99804 ± 0.00194 | 0.99005 ± 0.00193 | 0.99216 |

Each of six comparisons uses 8,192 independent paired preparations. Both widths see the same primary outcomes and the same random bits when they examine the same closing event. The directly paired change at length 8 was **−0.007994 ± 0.000173**, with 2,096 different closing decisions. The uncertainties are one standard error over independent reunion words. The stochastic cohort conditions on the completed transfer schedule. The distinct register run executes every physical transfer on the common 20-, 80- and 320-cell federations. A selected 20-cell witness with a positive eighth primary lets both owners read the consequences of the first return. [Material histories, finite predictions and cohort estimates](RESULTS.json) · [Individual reunion words](reunions.npz)

## What the other patches read

The research workspace has a searchable catalogue of nine defects, excitations and comparison preparations on the shared κ = 8 seam. The registered million-patch run contributes counts for signed repairs and mobile doubles. Separate pair experiments measure what a patch can detect after two loads approach or a scheduled repair writes a seam. Compact shapes with three to six loads have fixed-load calculations. All 17 physical Standard Model identifications remain open. [Local pair measurements](../../experiments/interval_echo_encounter/README.md) · [Particle identification criteria](../../evaluation/SCORECARD.md#standard-model-particle-identifications)

On the 320-cell sphere, 320 singleton observers and 960 overlapping three-cell observers read 23 signed repair records and tracked 34 conserved mobile doubles through 2,399 transfers. The same finite hands posted 41,108 positive marks and retained 40,812 quiet opportunities. All shared checkpoint readings agreed. The preparation probabilities, matching schedule, tag protection and gate windows belong to the experiment. An overlap comparison at one checkpoint and a unique outcome after differently ordered repairs are distinct measurements. The two-seam order test produced distinguishable terminal records. A prepared local priority resolved that fork in all six sharp preparation and priority combinations. [Bulk measurements](observations.json) · [Raw patch records](patch_records.npz) · [Consensus experiment](CONSENSUS.md)

The port-current measurement gives a direct test of the gauge premise. The protected menu computed from the rotations of the native cell has twelve independent responses and closes under every one of the 66 pairwise commutators. Allowing an additional bare port-load response gives seventeen visible directions. Identifying the complete menu available to a physical observer matters for the inference to a gauge algebra. [Port instrument and complete local outcome tables](CURRENTS.md)

The one-load sector passes the common rank-three cell test. The joint rest-energy and inertia requirement fails for the inherited free Hamiltonian. The rank-three transport of an observer’s full history, a universal relativistic clock comparison, a common stress and current measurement, and physical particle charges are open. [Scorecard](../../evaluation/SCORECARD.md)

## Measurements that distinguish the design choices

| Question | Experiment and result |
| --- | --- |
| Which repairs give a clock its rate? | A coherent seam oscillator and the registered hand agree in all three guided comparisons when they use accepted local ticks, including a fresh equal-leg replication. In the specified guided windows, counting all positive primaries raises the traveller/home ratio by 3.2–4.4 percentage points. A single trip followed by rest dilutes the fractional gap. At slower sustained motion the terminal oscillator can miss a gap that repair records resolve. [Duration and sensitivity](../../experiments/interval_echo_clock/DURATION.md) [Clock comparison](../../experiments/interval_echo_clock/README.md) |
| Does widening the hand fix the optical discrepancy? | On the same held 8/8 trip, the three-bit hand missed the optical prediction by +0.007606 ± 0.000248. The five-bit hand differed by −0.000035 ± 0.000261. Both uncertainties are one standard error. [Held links and speeds](../../experiments/interval_echo_optics/README.md) |
| Which local spatial responses propagate? | Three neighbouring patches jointly decoded 144 of 144 held load patterns. A local survey predicted 72 of 72 further patterns at held stations and exposures. The clock is idle and the dwell is supplied. [Intercell propagation](../../experiments/interval_echo_space/README.md) |
| Can observers recover three spatial directions? | In the programmed scalar instrument, eight anonymous stations yield three independent binary directions. The reconstructed chart predicts 24 held distributions within 0.001022 absolute error. Its prescribed refinement gives a flat 3+1D free scalar wave. [Clocked scalar measurements and theorem](../../experiments/seam_program/CLOCKED_FIELD.md) |
| Can observers range free propagation on the federation? | The flat fit misses the 5% tolerance on every held path. On 20, 80 and 320 cells the measured angular screen and calibration ranges also exclude every radial assignment in the round-screen, constant-lapse product. A shared autonomous geometry remains open. [Station measurements and controls](../../experiments/interval_echo_space/README.md#can-radar-supply-depth-on-a-round-angular-screen) |
| Does the first read prepare the pair? | A free ensemble of separated loads has stable-band weight 0.0029121. A sharp positive reading selects a branch with weight 0.9955693. Averaging over every result of the first reading also changes the band weight. A scheduled coherent write retained in the seam memory increased the chance of a later pair on disjoint seams from 0.3485% to 12.6725%. [Encounter and memory controls](../../experiments/interval_echo_encounter/README.md) |
| Can an accepted repair preserve phase? | A two-cell window requires absolute written contrast two. Giving a ready preferred seam the first commit yields the same complete endpoint for either request order and preserves a locally readable relative phase within that address. Interrogating both seams erases it. [Two-cell consensus](../../experiments/interval_echo_consensus/README.md) |

## Reproduce or inspect

Run from the RER repository root.

```sh
python -m model
python -m model --full
python -m model --full --seed 42 --output /tmp/echo-42
python -m model --full --seed 42 --verify /tmp/echo-42
python -m pytest -q model/candidate_laws/interval_echo
```

The following reads one individual reunion and one three-cell observer from the deposited data. Columns of the reunion are home, three-bit traveller and four-bit traveller. The four fields for each cell are signed repair in units of two, presence of a mobile double, last primary and hand position modulo five. The middle cell owns the triple’s reference hand. [Data producer](observations.py)

```python
from pathlib import Path
import numpy as np

folder = Path("model/candidate_laws/interval_echo")
with np.load(folder / "reunions.npz") as runs:
    print(runs["uniform_8"][0])
with np.load(folder / "patch_records.npz") as runs:
    print(runs["L2_triple_support"][0], runs["L2_triple_checkpoint"][0])
```
