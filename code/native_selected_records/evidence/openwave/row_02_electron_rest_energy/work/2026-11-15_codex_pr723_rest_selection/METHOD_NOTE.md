# Method and results

The literal endpoint event is

    j=(q_a-q_b-(2c-1)(abs(q_a-q_b) mod2))/2;
    q'=q-j e_a+j e_b.

The source-selected bond reference has

    P=[[3/4,1/4],[1/2,1/2]], pi=(2/3,1/3), Delta=log2,
    Pr(W/Delta=-1,0,1)=(1/6,1/2,1/3).

| Term/result | Executable definition | Recorded result |
| --- | --- | --- |
| Endpoint current and source | inputs/native.py::native; literal sha in CANDIDATE/PREREG | 14,520 census calls; 28,685 persistent native events |
| Compact ray under actual unread source | experiment.py current-group outer products, starting line61 | fidelity 2503/2520; outside mass1/315 |
| Source-directed bond selection | trajectories.py::Observer.consume and Trial.discover | all 24 signed starting preparations discover their own first-movement bond |
| Record-fed stationary reference | experiment.py complete HOLD columns; trajectories.py::Trial.epoch | all 120 bond/orientation/sign kernels P; feedback-disabled reference uniform |
| Perturbation recovery | experiment.py return census; trajectories.py::Trial.displacement_and_return | every start covered exactly; max32-attempt timeout33/2^32;24 trajectories return |
| Work distribution | experiment.py source work pulse after counted pi; trajectories.py::Trial.assay | exact population1/6,1/2,1/3;768 finite draws130/395/243 |
| Feedback reset | experiment.py reset columns and exact capped powers | exponential mean171/256, complete cap failure retained |
| Copy separation | experiment.py product/exclusive kernels; trajectories.py::pair_run | exact kernels;2048 persistent epochs each; independent4 versus exclusive3 native events per epoch |
| Independent verification | independent_check.py and trajectory_check.py | endpoint floor/ceiling, channel, reference, policy, tape and intervention checks pass |

Scientific preregistration: PREREG.txt and row-local LEDGER.md precede the first
census run. TRAJECTORY_PREREG.json separately precedes the persistent run after
root review requested it. Original reader and wrapper versions are retained in
amendments; checker-coverage and combined-replay updates have dated ledger hashes.
Local ordering is not upstream preregistration.

All finite kernels use exact fractions and integer branch counts. A graph-only
NumPy eigensolver provides the full 240 spectrum census. The independent reader
uses a different floor/ceiling endpoint formula and recomputes its own graph
Laplacian, importing no experiment or native module. The trajectory reader
reconstructs every observer decision from its public preparation and record tape.

Long kernel histories and joint Cartesian products are disclosed exact
compressions, not extra native calls. Persistent tapes store all actual events,
all gate/idle/selector coins, the displacement intervention, returns and every
work sample. Changed held-law control events are counted separately (384).
The equilibrium population work law uses pi from independently counted source
columns. Eight-epoch finite preparation has TV error at most 1/98304; no sample
is claimed exactly stationary. The finite sample mean 1553/1536 differs from1
within the disclosed diagnostic; the exact identity is only the population law.

Python3 and installed NumPy are required. `python3 headline.py` replays census and
trajectories, comparing all four output files byte for byte in the recorded
environment; the combined run took2.364 seconds. Portable independent readers
use declared float tolerances for spectrum/log diagnostics. Replay imports no
remote source and starts no subprocess. The optional mutation harness starts
only the already-read local checker files; thirteen mutations reject and two clean
controls pass.

The prospectively frozen quarter-gate alternative executes 6,240 additional
native calls and independently checks pi=(4/5,1/5), gap log4 and its work law.
The general relation gap=log(1/r) makes scheduler tunability explicit; the
localization support is likewise supplied by the address policy. Quarter epochs
use three native events, half-gate epochs two.

Not computed: physical electron identity, a common physical energy and vacuum
reference, controller/erasure-inclusive physical work, proper time, a sharp
nonzero rest-energy selected state or an unanchored measured mass/frequency.
The dynamically selected stationary state has transfer energy zero. Its modular
reference has two levels and nonzero variance. The original 511 keV anchor is used
only to reproduce its disclosed2.45084870234 keV separation failure.

**Scope:** positive finite record-fed selection/return and same-source modular
work, with retained failures of proposed mass dictionaries. No physical row
award, external-group confirmation, upstream pin change or publication.
