# Native nucleon current calculation

The calibrated common hadron source has exact vector charge normalization
`g_V=1`. Its native axial ratio lies in `[1.665734201267,1.667599132066]`,
which excludes the measured magnitude `1.2753 +/- 0.0013`. Calibration of
the masses does not calibrate the weak current. The complete colour-sector
current has nonzero excited-baryon leakage.

Run `python3 -m model.openwave_baryon_weak_formfactors` to regenerate the
exact integer local-current witness, rational finite-hopping bound and
floating current matrices. `PROOF.md` states the sector and comparison
assumptions. The empirical value in `inputs.json` is held out.
