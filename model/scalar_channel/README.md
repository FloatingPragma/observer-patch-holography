# Local preparation and observer records

An observer can prepare a constrained quantum state by carrying a scratch
bit between neighboring centers, reading it locally, and retaining the
outcome. This experimental scalar construction executes that protocol and
connects its prepared state to finite source-on pulses and complete pointer
readouts. It is an alternative to the primary
[contrast-memory model](../polynomial_sv/README.md), with different dynamics.

The [proof](PROOF.md) gives the source map and the complete instrument,
including interrupted routes and quantum references. The
[construction register](OBLIGATIONS.json) distinguishes this finite
implementation from autonomous control and interacting matter.

## Reproduce the construction

Install `model/requirements.txt`. From the repository root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -m model.scalar_channel
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m pytest -q model/scalar_channel model/test_contrast_replay.py
```

On PowerShell, set `$env:OPENBLAS_NUM_THREADS='1'` and
`$env:OMP_NUM_THREADS='1'` before the Python commands. The default replay
checks the source freeze and regenerates every result without rewriting
evidence. The extended replay uses 256 ramp substeps per star:

```sh
python -B -O -m model.scalar_channel --extended
```

`protocol.json` specifies the construction, controls and gates.
`FREEZE.json` binds its production code and inputs. `results.json` and
`extended.json` retain the complete finite outputs, including complementary
sectors and failure controls. Recording requires an explicit unused path
through `--record`; existing evidence cannot be overwritten by that command.

The construction is described in
[the construction manual](../construction_manual/CONSTRUCTION_MANUAL.md), whose
section 11 says which of its parts this package executes. No sibling checkout
or network service is needed for replay.

## The cage

The subpackage [cage](cage/README.md) builds an object that this law holds on
a ring of four seams, and executes an observer's experiments on it: the swing
of its seam record at one line, the level difference of its two modes from a
fuel dial, and a detector that reads its excitation after five waits. It has
its own receipts and its own replay:

```sh
python -B -O -m model.scalar_channel.cage
```
