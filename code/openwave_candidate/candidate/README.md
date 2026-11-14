# Portable OPH construction

The active candidate uses one canonical port/transit geometry, modular half-seam
records and a local coarse C/E instrument. Its explicit uniform-record
preparation admits an exact reduced simulation of every native outcome; the
same dynamics also supports definite-record preparations with a certified
sparse observer. Physical charge/species and many-load completion remain open.
From this directory, install with `python -m pip install .`, then run:

```sh
oph-openwave uniform_native --output native-observer.json
oph-openwave compact_native --output compact-state.json
oph-openwave sparse_native --output definite-record-control.json
```

The [uniform-fiber theorem](methods/uniform_native.md) applies to the entire
local instrument, including per-seam continuous scheduling and external quantum
references. The [compact preparation](methods/compact_native.md) persists under
that same instrument. It is a flat transit mode, not an identified electron.
Uniform and definite records are distinct inputs. Native holonomy-preserving
operations cannot raise a definite record's initial zero-holonomy weight;
the uniform resource requires a supplied preparation service or an explicitly
accounted conditional protocol. The bounded integer `JointRecords` instrument
retains its fully monitored comparison domain.

The full 31-row dispatch and OpenWave task handoff live in `../integration/`
and `../platform_packet/`. Only executed supporting observables receive local
results; missing physical observables remain explicit prerequisites. M10 is
admitted to OpenWave; the local computations award no official cells.

## Earlier independent finite computations

The package recomputes five finite constructions from explicit equations:
charge weights, a gauge-singlet spin kernel, the lepton read register,
class-conjugate repair branches and a twelve-direction scalar symbol.
Each command retains its mathematical domain and tests altered inputs.

**Scope**: these are calculations on the stated finite constructions. No
maintainer-awarded cell, complete physical particle criterion or empirical
prediction is claimed. Installing both packages
tests dependency compatibility and headless execution; it does not add an
OpenWave model column or a Taichi field kernel.

## Clean reproduction

Use Python 3.12 or later. The public source is
[RER](https://github.com/FloatingPragma/observer-patch-holography).
From a clean clone of that repository:

```sh
python3 -m venv .venv-openwave
. .venv-openwave/bin/activate
python -m pip install ./code/openwave_candidate/candidate
cd /tmp
oph-openwave all --output oph-portable-results.json
```

The command works outside the checkout and imports no RER module.
For the exact OpenWave dependency-compatibility check, install the pinned
upstream source in the same environment before installing the candidate:

```sh
git clone https://github.com/openwave-labs/openwave.git
git -C openwave checkout ef55175d04737e73c92f33f9f9cf3d17c593da54
python -m pip install ./openwave
```

The tested Python and library versions are recorded in the
[portable receipt](../evidence/portable_run.json). Exact formulas and
declarations are in the method notes. No large dataset or private checkout
is required by these five commands:

| Command | Method |
| --- | --- |
| `oph-openwave charge` | [Charge lattice and its limits](methods/charge.md) |
| `oph-openwave spin` | [Finite spin and exterior exchange](methods/spin.md) |
| `oph-openwave leptons` | [Register and conditional mass comparison](methods/leptons.md) |
| `oph-openwave conjugation` | [Complete local conjugation](methods/conjugation.md) |
| `oph-openwave dispersion` | [Definition-only scalar expansion](methods/dispersion.md) |

Install the test extra with `python -m pip install './code/openwave_candidate/candidate[test]'`
and run `python -m pytest code/openwave_candidate/candidate/tests` from the
RER root. The receipt records each scientific payload, checks, mutation
controls, versions, elapsed time and candidate source hashes.

## Independent audits

The [charge/spin audit](audits/charge_spin_review.py) uses polynomial
generating functions, vector-times-spinor generators and a permutation
oracle. The [dispersion oracle](audits/dispersion_oracle.py) uses sign-paired
generating polynomials and independent high-precision evaluations. Both
dispersion implementations were completed from the definition packet
before their results were compared. This is a local mathematical
recomputation, not an external blind certification.

Source-parity audits and the exact quark replay additionally require a
separate public RER checkout at
`3bcacf1276ddb3ead35be2b9adcd63b6643e2d35`:

```sh
python candidate/audits/source_parity.py --rer /path/to/pinned-rer --output evidence/source_parity_audit.json
python candidate/audits/check_spin_source_chart.py --rer /path/to/pinned-rer
python candidate/tools/replay_native_quark.py --rer /path/to/pinned-rer --output evidence/native_quark.json
```

Run these from `code/openwave_candidate`. The exact quark replay rebuilds
the actual source matrices and selected spectral measures; it can take
several minutes. It retains success weight 1/2, failure weight 1/2 and the
equal up/down response of that source. This is not a physical quark-mass
prediction. The separate numerical quark audit uses Gram-orthonormal
matrices, direct diagonalization and resolvent solves.

Software inherits the Apache-2.0 terms of RER software. Authored method
notes retain the repository's documentation terms; no third-party paper
or private planning document is included.

The `construction` runner module integrates the audited finite incidence law:
`python -m oph_openwave.runner construction`. `Construction` uses the same
incidence matrix for cost and the massive port/transit Hamiltonian, evolves
independent loads in separate transit spaces, and applies literal native events
with exact signed integer current counters. The tests check the massive square,
additive eigenphases, independence, the mean native law and the Gauss identity.
These are local algebra checks. The coherent transit amplitudes are distinct
from integer current counters. The joined finite instrument below supplies an
explicit bridge; interacting-load statistics and the physical domain remain
open. Earlier row-specific controls remain available for comparison.

`python -m oph_openwave.runner quantum_records` runs the joined finite quantum
register protocol. `JointRecords` lifts the port/transit generator to integer
half-seam records, in the sector of a specified initial port. Its Gauss charge
includes transit occupancy and is exactly conserved by every half-hop. Integer
records are bounded by a declared cutoff; rejected boundary events are retained.
The monitored native instrument has rank-one Kraus operators labelled by source,
seam and coin, and includes noops. This makes the channel trace preserving even
when distinct sources have the same destination. `step` composes unitary evolution
with a supplied commit probability. It is a discrete protocol, with no asserted
fast-commit limit. The classical one-load C-current law is recovered at zero
coherent time away from the boundary. Transit events are explicitly retained as
noops. Full position/record monitoring, the cutoff and transit policy are supplied
choices. The instrument chooses C (coin,current) and E (coin,pre-repair contrast)
records with equal probability and retains both noops and each record value.
The separate offset-law dynamics, many-load joint register, and removal of the
integer cutoff remain open. This does not establish all 31 physical criteria.

## Local records and charge apparatus

`oph-openwave modular_native` executes the separately named
`LocalCoarseModular` candidate. It groups source amplitudes by the actual local
C/E outcome, retains coherence in degenerate outcomes, and uses Z6 half-seam
records. Its continuous generator assigns a supplied rate to each seam.
The complete small-sector comparison reproduces the change caused by global
monitoring. This is a distinct instrument and regulator, with Gauss residues
rather than an integer Gauss-vector theorem. See
[the local instrument](methods/modular_native.md).

`oph-openwave regional_charge` measures complete regional charge outcomes on
the bounded integer construction, with exact boundary and continuity accounting.
`oph-openwave general_cut` proves the diagonal cut identity on the canonical
ring without enumerating its quantum sector. The [apparatus](methods/regional_charge.md)
and [general-cut contract](methods/general_cut.md) distinguish fractional
expectations from sharp charge, integer circulation from complete-boundary
current, and operational load units from physical electric charge.

These commands share the same port/transit geometry and explicit record
conventions. Uniform-fiber reduction computes canonical-screen records on its
supplied code; direct full-sector enumeration is limited to small controls.
A terminal effect can be compressed even when measuring it leaves that code.
Continuing the reduced instrument requires its actual outcome maps. Sparse
definite-record certificates retain their separate preparation/error domain.

`oph-openwave sparse_native` executes the canonical-ring local-coarse channel
on a sparse modular path window. Its exact rational Taylor certificate retains
numerical norm defect and bounds complete outcome probabilities against the
full modular sector for one interval from the stated basis preparation.
See [the error certificate](methods/sparse_native.md). To retain every outcome:

```sh
python -m oph_openwave.sparse_native --full-ledger --output native-window.json
```
