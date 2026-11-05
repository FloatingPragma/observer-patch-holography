# Reproducing the IPI dispersion note

The reader-facing source is `../OPH_FROZEN_DISPERSION_AUDIT_2026-11-02.md`.
The four-page PDF is beside the note at
`docs/comms/OPH_FROZEN_DISPERSION_AUDIT_2026-11-02.pdf`.

`audit_results.json` records the independent geometry, coefficient and negative
controls, original frozen-artifact hashes, and relevant working-tree input hashes.
`verification.txt` records 304 passing focused tests and four independent
verifier commands in the initial audit, plus 53 further tests for direct
carrier signatures and the existing field/clock link (357 total).
The intuitive spacetime revision adds 12 passing local-clock/motion checks
(369 targeted tests overall) and preserves the finite theorem premises.
`direct_tests_results.json` records 3,960 exact native branch checks and a
separate response-rank proposal, including the missing physical identification.
The audit does not claim an all-repository test run, a new
Lean kernel build, verified external timestamp anchoring, or an empirical test.
The PDF's hypothetical measurement and Planck-length example are illustrations.

From the repository root:

```sh
python3 docs/comms/dispersion-audit-2026-11-02/audit_dispersion.py
python3 docs/comms/dispersion-audit-2026-11-02/audit_direct_tests.py
python3 docs/comms/dispersion-audit-2026-11-02/build_pdf.py
pdftoppm -scale-to 1600 -png \
  docs/comms/OPH_FROZEN_DISPERSION_AUDIT_2026-11-02.pdf \
  tmp/pdfs/oph-dispersion-page
```

The audit needs NumPy, SymPy and mpmath. The PDF builder needs ReportLab,
Matplotlib, Pillow and pypdf, and uses macOS's Times New Roman fonts. It checks
for exactly four pages; inspect all four rendered pages after changing the note.
The original build used an isolated environment with those dependencies.
The source predictions and executable model were read without altering their
frozen artifacts or source files.

The retained JSON and verification log preserve the original audit on working-tree
inputs at baseline commit `fc2aab8c0e87279c26f45828ece2d6be79913f51`. Their
hashes are historical provenance, not assertions that current main has identical
source bytes. The scripts were relocated from the meta workspace by changing
repository/path resolution only. Historical command paths in `verification.txt`
refer to that original layout; use the commands above in this repository.
Rerunning an audit writes a new receipt for the current checkout.
