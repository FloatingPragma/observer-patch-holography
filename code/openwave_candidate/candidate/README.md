# Portable finite OPH computations

The package recomputes five finite constructions from explicit equations:
charge weights, a gauge-singlet spin kernel, the lepton read register,
class-conjugate repair branches and a twelve-direction scalar symbol.
Each command retains its mathematical domain and tests altered inputs.

**Scope**: these are calculations on the stated finite constructions. No
OpenWave admission, maintainer-awarded cell, complete physical particle
criterion or empirical prediction is claimed. Installing both packages
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
