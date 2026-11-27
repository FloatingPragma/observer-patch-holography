# Local preparation and observer records

An observer can prepare a constrained quantum state by carrying a scratch
bit between neighboring centers, reading it locally, and retaining the
outcome. This experimental scalar construction executes that protocol and
connects its prepared state to finite source-on pulses and complete pointer
readouts. It is an alternative to the primary
[contrast-memory model](../contrast_memory_spherical_v1/README.md), with different dynamics.

The [proof](PROOF.md) gives the source map and the complete instrument,
including interrupted routes and quantum references. The
[construction register](OBLIGATIONS.json) distinguishes this finite
implementation from autonomous control and interacting matter.

To decode a one-body operator on a connected test graph, call
`model.candidate_laws.scalar_channel.one_body.read_one_body(graph, terms)`.
It returns the one-body matrix on three or more modes. On a two-mode edge
it returns the identifiable trace and sets `matrix` to `None`. Both results
include the maximum residual in the Pauli coefficients acting on the cycle-$+1$
code. Multiplying a term by a cycle stabilizer gives the same action there.
Terms that send the input outside the code enter the residual. The archived
`instruments.decode_paulis` assumes three or more modes and runs on 174
in the registered instrument. [Even-sector identification proof](PROOF.md)

## Reproduce the construction

The registered cell operator couples seam and wedge variables. Squaring its two-block form gives the identity ω²=m²+λ for its incidence Gram eigenvalues. The observed finite cage line uses its prepared clock and reference. The cell evaluation does not test a larger screened federation or an emergent light cone. The existing cage experiment already constructs signed seam–wedge incidence on the 20-, 80- and 320-cell boards with the same block rule. Its federation operator uses the [cage constructor](../../experiments/cage/construction.py); the registered cell-only compiler and instrument still need a federation implementation for a full end-to-end comparison. The registered current vanishes on its projected accessible code; enlarging the code would change that test.

Install `model/requirements.txt`. From the repository root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.candidate_laws.scalar_channel
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m pytest -q model/candidate_laws/scalar_channel model/candidate_laws/contrast_memory_spherical_v1/test_replay.py
```

On PowerShell, set `$env:OPENBLAS_NUM_THREADS='1'` and
`$env:OMP_NUM_THREADS='1'` before the Python commands. The default replay
checks the source freeze and regenerates every result without rewriting
evidence. The extended replay uses 256 ramp substeps per star:

```sh
python -B -O -m model.candidate_laws.scalar_channel --extended
```

`protocol.json` specifies the construction, controls and gates.
`FREEZE.json` binds its production code and inputs. `results.json` and
`extended.json` retain the complete finite outputs, including complementary
sectors and failure controls. Recording requires an explicit unused path
through `--record`; existing evidence cannot be overwritten by that command.

The construction is described in
[the construction manual](https://github.com/FloatingPragma/observer-patch-holography/blob/512b47075b0dd00a570b08cc9b576bec7559acb5/model/construction_manual/CONSTRUCTION_MANUAL.md), whose
section 11 says which of its parts this package executes. No sibling checkout
or network service is needed for replay.

## The cage

The subpackage [cage](../../experiments/cage/README.md) builds an object that this law holds on
a ring of four seams, and executes an observer's experiments on it: the swing
of its seam record at one line, the level difference of its two modes from a
fuel dial, and a detector that reads its excitation after five waits. It has
its own receipts and its own replay:

```sh
python -B -O -m model.experiments.cage
```

The experimental member [four-child register](../../experiments/refinement_spectrum/README.md) couples the hop of the
scalar incidence, without the level term, to an internal register with uniform hopping. Its clock and fuel
records measure rest frequencies and transfer between cells. The complete spectrum bounds the ratio of
successive rest frequencies at every tree depth.

## Local conversion of opposite charges

The experimental [pair-conversion law](../../experiments/radiative_repair/README.md) is a separate generator on the port
carrier with one common energy scale. Opposite charges prepared in its lowest band release energy into a
neutral field. In a second experiment, which starts from a pair localized on one seam, a meter in a
neighbouring cell registers two neutral arrivals and a converter there turns them back into a charged pair.

The [mobile-memory candidate](../../experiments/mobile_memory/README.md) has a mobile lowest charged band and a cage
link weighted by the native contrast square or a retained signed-current record. Its joint records
compare bound-pair presence, cage arrival, and the effect of reading the stored sign at a later event.
