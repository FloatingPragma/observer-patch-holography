# A finite charged-energy detector

A twenty-qubit phase pointer measures the actual finite charged source
using its original Hamiltonian and an explicit neutral-energy reference
phase. A complete compensated CAR rotation prepares the input with
probability 1/4. Every phase outcome and preparation failure is retained.
The accepted source has ground fidelity above 0.99999999439512, including
a systematic rational-gate approximation bound. This is a finite boundary
source measurement, with no relativistic particle identification.

The [proof](PROOF.md) gives the source, complete apparatus, clock and
logical capacity. The source's complete 70928-state residual verifier is
replayed. The much larger source-pointer matrix is controlled analytically.
Exact finite coefficient and linear-solve services are declared; extra
gate-synthesis workspace and efficient local hardware are outside the
implementation claim.

```sh
python -m model.affine_gauss_charged_clock --write
python -m model.affine_gauss_charged_clock --check
python -m pytest -q -W error model/affine_gauss_charged_clock
```

The [source ancestry](PROVENANCE.md) identifies the frozen independent
calculation and the portable public dependency closure.
