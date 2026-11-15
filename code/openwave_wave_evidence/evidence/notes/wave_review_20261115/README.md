# Wave review package

Start with **[REVIEW.md](REVIEW.md)**. It records four concrete replay findings,
the independent finite checks, and the scope of the written-proof assessment.

- [Proof review](PROOF_REVIEW.md)
- [Draft platform task](SUBMISSION_TASK_DRAFT.md)
- [Review plan](PLAN.md) and [review-tool attempt history](REVIEWER_ATTEMPTS.md)
- `results/`: machine receipts, full captured subprocess output, profile
  diagnosis and the exact reviewer launcher used in each environment
- `MANIFEST.sha256`: review-package seal; original source seals are separate

## Reproduce

Requires Git, Python 3 and NumPy. Run from this repository's root:

```sh
python -B -X utf8 notes/wave_review_20261115/reproduce.py --label local --output /path/to/new-results
python -B -X utf8 notes/wave_review_20261115/independent.py --output /path/to/new-results/INDEPENDENT.json
python -B -X utf8 notes/wave_review_20261115/diagnose_profiles.py /path/to/new-results/local/REPRODUCTION.json --output /path/to/new-results/PROFILE_DIAGNOSTICS.json
python -B -X utf8 notes/wave_review_20261115/verify_review.py --check-seal
```

Choose a writable output path outside this sealed review directory. The
reproducer creates a unique scratch copy, extracts the exact `60519f4` Git
objects and checks their blob IDs. It executes the original bundle, then the
scientific leaves to diagnose wrapper failures, and finally optimized/mutated
verification probes. Scratch copies are retained for diagnosis.

On Windows, an optional `--wsl-python /absolute/linux/python` argument runs
the same subjects through Ubuntu WSL. This review installed NumPy 2.3.5 into
a separate temporary WSL virtual environment; its environment receipt records
the actual Python, NumPy, platform and zlib versions. No project dependency
or system Python installation was changed.

`reproduce.py` exits successfully when it has produced a complete **review
report**, even if that report records a failed original replay. Consult
`original_replay.exit_code`, each leaf's exit code and the byte comparisons.
The optimized runs are deliberate negative verification probes and must not
be used to claim a scientific pass.

`independent.py` fails explicitly on an arithmetic/record mismatch.
`verify_review.py --check-seal` verifies the saved review claims, log/source
bindings, comparator negative cases and review seal. It does not rerun the
scientific experiments or award a gate.

The post hoc diagnostic `1e-10` comparison classifies discrepancies only.
It is not a new acceptance policy. Original failed byte gates remain failed.
