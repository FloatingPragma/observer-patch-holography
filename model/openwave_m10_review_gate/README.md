# OPH's OpenWave column

OPH is admitted as M10. Its31 physics rows are untested. OpenWave's
independent arithmetic review reproduces seven dispersion coefficients;
the review gives no physics row for that calculation. The executable here
authenticates that status and reruns all three upstream implementations.

Run `python -m model.openwave_m10_review_gate --write --openwave-root <checkout>`
and the same command with `--check`. The checkout must be the exact upstream
revision in [contract.json](contract.json). [PROOF.md](PROOF.md) states the
native-process bridge required for a criterion claim.
