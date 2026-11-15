# Actual charged ground and its field dressing

The complete70928-state one-charge sector of the fixed mobile Gauss source
has a unique ground energy in
`[1.8738097549651658,1.8738097549671659]`, in source clock units.
The retained charge compensator contributes exactly1. Subtracting that
known reference energy gives the source transition contribution near
`0.873809754966`. These are finite excitation energies, not GeV masses.

The verifier applies the entire Hamiltonian to an untrusted stored trial
vector, with outward source-coefficient enclosures and explicit arithmetic
error bounds. It admits only a residual at most`10^-12`. The independently
proved full-complement second-level bound then isolates the actual ground.
No photon states are removed beyond the source's declared total-occupation
cutoff three. The actual ground has a strictly positive transverse-field
occupation probability near`0.00007224654`.

Complete position records and a compensated preparation with retained
failure and energy-complement outcomes are included. [The proof](PROOF.md)
states every source, arithmetic and apparatus premise. The fixed boundary
compensator still prevents interpreting this result as a translation-invariant
relativistic rest mass or an identified Standard Model particle.

```sh
python -m model.affine_gauss_full_spectrum --check
python -m pytest -q -W error model/affine_gauss_full_spectrum
```

The verifier uses no numerical eigensolver. Reconstructing the trial is
optional: `python -m model.affine_gauss_full_spectrum.trial --output /tmp/charged-trial.npy`
generates an untrusted candidate without changing the committed receipt.
Only the complete residual certificate admits a trial.
