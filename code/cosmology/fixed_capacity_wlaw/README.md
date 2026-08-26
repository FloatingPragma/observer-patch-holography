# Fixed-capacity equation-of-state comparisons

This package separates the exact conditional capacity law from empirical
comparisons. `FixedCapacityWLaw.lean` proves the algebraic implication

\[
w(a)=-1+\frac{1}{3}\frac{d\ln N}{d\ln a}
\]

after the density and continuity definitions are supplied. It does not derive
those definitions or a capacity history.

`official_desi_dr2_chain_audit.py` postprocesses the four official DESI DR2
Cobaya chains for each default CMB combination and checks every chain against
the collaboration's SHA-256 manifest. For CPL on `0 <= z <= 2`, the conditional
monotone-capacity subset is exactly

```text
w0 >= -1  and  w0 + (2/3) wa >= -1.
```

The same run also hash-checks the four official flat base-`Lambda`CDM
BAO+CMB chains and computes `Lambda*l_P^2` for every weighted sample directly
from the chain's paired `H0` and `omegal` columns. This sample-level posterior
replaces any assumed or hand-entered `H0`--density correlation. The displayed
conversion uses stated central SI constants and does not propagate their much
smaller uncertainties.

Download the public inputs and reproduce the receipt:

```bash
python3 code/cosmology/fixed_capacity_wlaw/download_official_desi_dr2_chains.py \
  --data-dir /path/to/desi-dr2-chains
python3 code/cosmology/fixed_capacity_wlaw/official_desi_dr2_chain_audit.py \
  --data-dir /path/to/desi-dr2-chains \
  --output code/cosmology/fixed_capacity_wlaw/runtime/official_desi_dr2_fz13_retrospective.json
python3 -m pytest -q \
  code/cosmology/fixed_capacity_wlaw/test_official_desi_dr2_chain_audit.py
```

The output is a retrospective, model- and prior-dependent posterior diagnostic.
The rare monotone-subset fractions are reported with per-chain tail counts and
resolution warnings; they are not a branch evidence ratio, a frequentist
exclusion, a direct capacity measurement, an OPH confirmation, or a frozen
FZ-13 score. Agreement with the fixed point is shared with LambdaCDM and earns
no OPH-specific credit.
