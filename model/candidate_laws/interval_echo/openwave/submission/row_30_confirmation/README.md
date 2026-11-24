# OW-30 candidate confirmation

The nominated object is a moving reversal in the orientations of loaded dimers on a commissioned strand. The proposal binds the state, memory convention, source schedule, local clock, records and numerical decoder to one Git revision. The canonical pin, source pointer and explicit acceptance of the declaration must agree on official OpenWave main before a public selector is derived. The local registration patch proposes those changes together. OpenWave's public decision remains separate. [Candidate](../../row_30_collective/confirmation/candidate.json) · [Terms](../../row_30_collective/confirmation/TASK.md) · [Canonical admission rule](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#L22-L27)

**Author-side result.** The frozen v3 source `92370140edfbdb0ac2be647b2b8166d821e5f1a2` passed on all three reserved local validation boards. Two local coupling measurements forecast 22 positive frequencies within 0.000575 on the lattice ruler. Fourteen lines with ordinary wave number at most 0.50 agreed within 0.003165, against a tolerance of 0.0048. All 170 count arrays were present. OpenWave has no landing for these terms or decision on this candidate. [Validation measurements and counts manifest](validation.json) · [Package revisions and checksums](PACKAGE.json) · [Version-bound review](REVIEW.md) · [Public row](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L467)

| v3 validation board | Strand dimers | Enrolled frequencies | Odd-hold residual with measured wavelengths | Pass ceiling | Same records with launch wavelengths |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 cells | 23 | 6 | 422.56 | 554.07 | 2,786,515 |
| 80 cells | 36 | 8 | 860.42 | 1137.74 | 792,450 |
| 320 cells | 30 | 8 | 786.58 | 963.37 | 1,794,431 |

The odd-hold spatial residual is a weighted sum of squares over all held receiver addresses and scans. The independently trained wave numbers are used directly for the dispersion forecasts. Blank and altered-coupling controls, held pair sums, ordinary and lattice wave numbers, and the matched lattice tangent all have separate gates. All three boards in each of the two exposed cohorts also passed. Changing flight coupling after calibration gives a complete scientific rejection on exposed boards. Removing an expected acquisition file gives an incomplete verdict. Corrupted counts give an invalid verdict. [Exposed pilot](pilot.json) · [Reanalysis under v3 of the v2 validation materials](repeat.json) · [Executable control tests](../../row_30_collective/confirmation/test_acquire.py)

## Why the wave-number reader changed

The v2 experiment used the wave numbers from the full-record, three-address recurrence in its dispersion gates. It also fitted the whole spatial profile on even holds and tested the result on odd holds. The three exposed geometries passed. On the three independently reserved geometries every frequency, control and transport gate passed, and substituting the launch labels failed on every board. Two of the recurrence estimates missed the odd-hold spatial cutoff. The all-board v2 result is **rejected**, with all 130 required count arrays present. [Frozen version and hashes](v2/PACKAGE.json) · [Complete validation manifest](v2/validation.json)

| v2 validation board | Cutoff | Interior-recurrence wave numbers | Spatial fit from even holds | Source launch labels |
| ---: | ---: | ---: | ---: | ---: |
| 20 cells | 190.89 | 120.58 | 92.44 | 20,289,334 |
| 80 cells | 416.58 | 551.03 | 303.44 | 2,540,128 |
| 320 cells | 1166.64 | 1252.79 | 954.34 | 692,904 |

Each column scores identical receiver counts. A lower residual fits better. On the 80-cell board the largest improvements come from the fourth and fifth measured wave numbers. On the 320-cell board, changing the third by about 0.00002255 radian per dimer is sufficient to bring the reported recurrence prediction below cutoff. Reapplying all wavelength-dependent gates with the even-hold fit to the six exposed boards passes on each. That retrospective diagnosis informed the v3 candidate. The v2 verdict is preserved. [Reanalysis code](diagnose.py) · [Record-level comparison](diagnosis.json)

## What the successor tests

V3 fits wave numbers from even-hold receiver profiles and uses those same measurements to forecast the odd-hold shape and to compare dispersion. Temporal frequencies come from separate counters. The foil substitutes the source wavelengths on those same odd holds. The cutoff, coupling calibration, shot budgets, geometries, enrollment rule, control amplitudes and tolerances use the published definitions in the terms. The old recurrence and its residuals remain as diagnostics. The changed decoder requires its own candidate identity and fresh confirmation.

Three fixed validation seeds, 785029, 514823 and 690271, occupy a different domain from the prospective public draw. The v2 boards remain enrolled in a separate repeat. Separate measurement output, a read-only inspection and an independent replay use the same bound source and draw. Acquisitions record physical rejections separately from incomplete or invalid evidence. [Reader](../../row_30_collective/confirmation/spatial.py) · [Acquisition and exits](../../row_30_collective/confirmation/acquire.py)

## Reproduce from the frozen checkout

The source bundle is a local, unsubmitted Git object with a public RER base. The OpenWave patch stages the exact source, port, canonical selection and proposed acceptance record for review. That local artifact cannot satisfy the official acceptance check. The working research checkout contains additional edits to the common candidate wrapper. The following commands keep the research checkout intact and use the version in the packet. [Package manifest](PACKAGE.json) · [Source Git bundle](SOURCE.bundle) · [Proposed OpenWave patch](OPENWAVE.patch)

```sh
project=$(git rev-parse --show-toplevel)
packet="$project/model/candidate_laws/interval_echo/openwave/submission/row_30_confirmation"
trial=$(mktemp -d)
git clone --shared --no-checkout "$project" "$trial/rer"
git -C "$trial/rer" fetch "$packet/SOURCE.bundle" refs/heads/ow30-confirmation-v3
git -C "$trial/rer" sparse-checkout set model/candidate_laws/interval_echo/openwave model/candidate_laws/contrast_memory_spherical_v1 model/executable_model
git -C "$trial/rer" checkout --detach FETCH_HEAD
runner="$trial/rer/model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation/run.py"
python="$project/.venv/bin/python"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
"$python" -I -B "$runner" check-sources
"$python" -I -B "$runner" acquire --mode validation --release "$trial/validation"
"$python" -I -B "$runner" replay --mode validation --release "$trial/validation"
```

The one-thread settings reproduce the archived validation counts and diagnostics. In a default-thread inspection on the review host, the 80-cell recurrence diagnostic differed from the archive by about 7.14 × 10⁻⁸, above the frozen absolute replay tolerance of 10⁻⁸. The verifier returned `invalid` (4). A fresh one-thread replay reproduced all 170 count arrays and passed every scientific gate.

Use the pinned versions in `model/requirements.txt`. The source and OpenWave entry points compare every input to the frozen bytes. The public entry point additionally verifies its own launcher and vendor files against official OpenWave Git blobs. With the word `public` as the mode, an unaccepted candidate exits with code 5 before writing a draw. An open registration requires publishing the frozen source for independent inspection and a maintainer decision on the proposed pin. [Complete status contract](../../row_30_collective/confirmation/README.md)
