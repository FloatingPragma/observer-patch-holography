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

## Reproduce the construction

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
[the construction manual](../../construction_manual/CONSTRUCTION_MANUAL.md), whose
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
