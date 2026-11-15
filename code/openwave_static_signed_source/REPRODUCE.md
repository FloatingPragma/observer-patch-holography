# Reproduce the finite static experiments

Python 3.12 or 3.13, NumPy, SciPy and SymPy are required. From the repository root:

```sh
python3 -B code/openwave_static_signed_source/replay.py --verify-only
python3 -B code/openwave_static_signed_source/replay.py
```

The first command checks every published evidence file, closure of the file
list, and the original scientific input pins. The second runs the complete
compact-coordinate audit, localized branch reader, and common signed-source
experiment in a disposable copy. The original evidence is checked before and
after execution. The compact audit compares every retained result field except
elapsed time and artifact hashes, whose bytes are checked separately. The
common-source script compares its full regenerated receipt and branch records.

## Transport boundary

The historical export and preregistration include one CPython 3.13 bytecode
cache. Its source file is published and verified. The interpreter-dependent
cache is omitted from the public package. `PORTABILITY.json` identifies its
exact path and historical hash. `replay.py` removes only that cache entry from
the preregistration in its disposable execution copy. All other fields, source
files, observations, thresholds and expected results are unchanged. The
historical manifests are archived byte for byte; their raw checks require the
omitted cache and therefore are not the portable replay command.

The evidence folder contains 315 original scientific and protocol files plus
the original export manifest. `MANIFEST.json` binds exactly those 316 files.
The cache omission is a transport amendment, not a replacement scientific run.

The localized provider also pins 150 original paths whose identical bytes are
included under the export's `standalone/` prefixes. `PROVIDER_PATHS.json` maps
these aliases. The adapter checks each destination against that provider's
original pins and creates the aliases in scratch before running its reader.
There is no substituted theorem or changed source file.

That reader writes an elapsed-time receipt in scratch. Its regenerated fields
are compared with the frozen reader receipt, excluding elapsed time; the
original receipt bytes are restored before the common-source consumer checks
its input pins. The published receipt is never rewritten.

## Interpretation

The energy pointer certifies a selected eigenray and retains all complementary
outcomes as unresolved energy. The 511 keV reporting scale is a calibration.
The tensor/frame experiment uses a finite addressed-cell rotation. The
positive coherent memory law is a candidate primitive; passive repair alone
is not claimed to compile every body branch. Separate preparations assay one
source species. No physical clock, free lifetime, spin magnitude, exchange
statistics or global continuous rotation symmetry is inferred.

Archived documents use the provider's internal term `M6` for a native service;
it does not refer to OpenWave's M6 Ouroboros model. The source-pin and review
fields inside frozen provider records describe their original execution and
are not declarations of platform acceptance.
