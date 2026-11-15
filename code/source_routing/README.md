# Readback along a seam path

This package executes a declared calibration schedule using only scalar pair
means on simple paths in the captured twelve-port W12 support at level three.
The observer at the destination retains its own initial value and its value
after each calibration sweep. From these local records alone it reconstructs
the initial values of every port on the path, including the remote source.
The observer-like system includes mutable port state, local readback,
protected records, canonical repair feedback, and a public event bundle.

`specification.json` fixes the paths' lengths, preparations, interventions and
stopping rule. `support_w12_l3.json` is a finite snapshot of the production
simulator's declared support. `capture_support.py` can refresh that snapshot
from the sibling `oph-physics-sim` checkout without executing a simulation
campaign; its local source file hashes identify the capture. The recorded base checkout
commit supplies context and does not assert that those local source files
are retrievable at that commit. Port
assignment is supplied. The independent verifier checks the complete mesh
incidence, slot uniqueness, connected support and every routed seam.

From the research repository root:

```sh
python3 code/source_routing/build_routing.py
python3 code/source_routing/verify_routing.py
python3 -m pytest -q code/source_routing/test_routing.py
cd Lean
lake env lean Geometry/SourceSeamPathTomography.lean
```

The producer writes `runtime/path_tomography_receipt.json`. The verifier
imports no producer and independently constructs and inverts the exact
observation matrix, replays every versioned operation, checks source and
nonsource interventions, and recomputes all event, record and rational-payload
costs. Hash-chain order is custody; read-from writers generate causal edges.
The receiver's decoder sees only its protected local samples. No remote
initial value or remote read log enters that decoder.

The theorem and its scope are stated in
`paper/tex_fragments/SOURCE_SEAM_PATH_TOMOGRAPHY.tex`, with the recursive
inverse, injectivity, attenuation and seam-count algebra in
`Lean/Geometry/SourceSeamPathTomography.lean`. Tests reject changed laws,
missing events, stale versions, remote side channels, optimistic noise
bounds and capacity undercounts. A forward sweep without calibration is
also tested as an actual information-loss control.

This is an exact finite, destructive readout construction. Initial values,
path, protected history and schedule are declared; no other seam touches
the path during the episode. All mean operations in the retained runs are active;
the general theorem also permits equal-endpoint mean evaluations that are
identity updates, without claiming strict repair progress or fairness.
The calibrated reconstruction can be badly
conditioned. Payload-bit counts describe these finite rational runs and
count written scalar register payloads, excluding event metadata (including
quadratic decrements), temporary decoding arithmetic, identifiers, hashes
and machine representation overhead; no bounded
physical precision, bandwidth or clock is inferred. Repeated old-version
reads, a full metric-neighbor routing/refinement, population selection and
count-volume attachment are separate requirements.

The analytic observation-map condition number for every path length is
`g_d = ((1+sqrt(2))^(d+1) + (1-sqrt(2))^(d+1))/2`, with
`g_0=1`, `g_1=3`, `g_d=2*g_(d-1)+g_(d-2)`. The inverse rows have alternating
delay coefficients, so this is the exact sharp maximum-norm amplification of
independently bounded sample errors under exact dynamics and decoding.
At fixed range and reconstruction tolerance, the extra sample precision
therefore grows linearly, approximately 1.271553 bits per seam. The large
gain alone establishes neither impracticality nor a physical precision bound.
Exact matrix tests cover every depth from zero through 24; the four embedded
receipt depths keep their complete operational custody.
