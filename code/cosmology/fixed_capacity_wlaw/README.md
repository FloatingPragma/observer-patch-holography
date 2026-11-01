# Fixed-capacity equation-of-state comparisons

This package separates the exact conditional capacity law from empirical
comparisons. `FixedCapacityWLaw.lean` proves the algebraic implication

\[
w(a)=-1+\frac{1}{3}\frac{d\ln N}{d\ln a}
\]

after the density and continuity definitions are supplied. It does not derive
those definitions or a capacity history.

`capacity_alpha_interval_certificate.py` upgrades the local declared
capacity-to-alpha calculation from a centered finite-difference diagnostic to
an outward-rounded interval certificate. On the committed comparison branch
and the certified domain
`|Delta log(alpha)| <= 1e-5`, it verifies the sign-definite implicit
denominators for the nested `m_Z` and pixel-closure roots, proves that the
selected mathematical branch is continuously differentiable, and encloses

```text
d log(N) / d log(alpha) in [-0.2141738647, -0.2061760309]
d log(alpha) / d log(N) in [-4.850224325, -4.669103775]
```

The reciprocal is formed only after the first interval excludes zero. An
independent replay evaluates the same tangent through the factorized chain
rule and verifies source hashes, denominator signs, and mean-value bounds:

```bash
python3 code/cosmology/fixed_capacity_wlaw/capacity_alpha_interval_certificate.py
python3 code/cosmology/fixed_capacity_wlaw/verify_capacity_alpha_interval_certificate.py
python3 -m pytest -q \
  code/cosmology/fixed_capacity_wlaw/test_capacity_alpha_interval_certificate.py \
  code/cosmology/fixed_capacity_wlaw/test_capacity_alpha_tangent.py
```

This closes only branch differentiability for the declared finite formulas.
B1, B2, and the physical branch-selection and readout part of B3 remain
undischarged. The certificate supplies no physical epoch evolution. The
retrospective wrapper now uses the interval upper bound for measurements inside
the certified domain and emits no mapped bound for wider measurements.

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
python3 code/cosmology/fixed_capacity_wlaw/verify_official_desi_dr2_independent.py \
  --data-dir /path/to/desi-dr2-chains
python3 -m pytest -q \
  code/cosmology/fixed_capacity_wlaw code/cosmology/postdiction_ledger -W error
```

The v3 receipt reads the published decimal columns exactly. Both posterior
families use one rational weighted-moment accumulator and one chain parser.
Means, population covariances, weight-concentration ESS, subset decisions and
quantile thresholds are independent of row order and chain merge grouping.
The base-LambdaCDM conversion uses the printed central SI decimal constants
before taking moments; it does not square an already rounded tiny display.
The CPL endpoint test is `w0 >= -1 and 3*(w0+1)+2*wa >= 0`, evaluated exactly.
The excluded subset is counted before conversion to a displayed probability.

Final roots use a private 90-digit numerical context. A nonzero binary64
summary field must retain relative rounding accuracy of `1e-12`; otherwise
the producer refuses that summary with a field-specific error. This is a
numerical reporting criterion, not an interval certificate or a bound on
posterior sampling error. Singular covariance leaves the valid chain summary
intact and marks the optional two-dimensional Gaussian diagnostic unavailable.
A positive Gaussian survival probability outside resolved binary64 range is
`null`, with an explicit status and the log tail retained when representable.
The normal quantile uses a small-displacement `expm1` calculation and, for
larger displacements, SciPy's inverse log-CDF. The diagnostic uses the exact
covariance, not the rounded covariance reconstructed from displayed fields.

The independent verifier imports no producer. It pins all twenty input files
separately, recomputes all 25 chain/combined summaries with centered Decimal
moments, and checks the four Gaussian diagnostics through elimination and an
independent `erfc` evaluation. It requires the full default receipt. CI runs
small complete fixtures and corrupt-evidence controls on Linux;
the public-chain replay uses the explicit download command above.

The review of this implementation found that numerical replay alone accepted
altered redshift ranges, formulas, provenance and physical interpretations.
The verifier now also binds the complete scientific context to the reviewed
v3 contract at `a14a3dac`, with a separate canonical digest. Computed statistics
and per-file records are excluded from that digest and must pass independent
replay. Unknown fields, duplicate JSON keys, non-finite JSON constants and
Boolean chain identities are refused. The verifier still executes without
importing the producer. Changes to the declared context require reviewing and
updating this explicit binding; passing arithmetic alone cannot authorize them.

The manifest URL is mutable. On 2026-10-07 it returned 226,892 bytes with SHA-256
`ec6614202355892f0d8b4c70934e73e215527d100faebc0ca0e5995d0b444910`
(HTTP Last-Modified: 2026-10-02); all twenty retained chain paths and hashes
agree with that fetched manifest. The receipt preserves the historical
manifest pin `df78872aa8b2d3473a9e8de78f498180efd7cbcbeb18211ce4787fac52067ee5`
recorded by the August producer. It is not the hash of today's URL contents.
Offline replay checks the retained chain identities directly; it neither
substitutes newer data nor claims to recover the earlier full manifest bytes.

The accounting audit reproduced erased narrow spreads, accepted negative
CPL endpoint offsets, lost excluded tails and quantile atoms, weight-scaling
overflow, and crashes on valid Gaussian tails. The original 30 controls had
22 failures on `de60560b`. Replaying the 416,469 published rows changes 188
numeric fields by at most `4.10e-10` relative; all subset masses, tail counts
and source hashes are unchanged. The live postdiction ledger is canonically
regenerated with the new parent binding and three refreshed citation anchors
whose formulas/statements remain in current main. Its verdicts and displayed
comparisons are unchanged. Papers, claim payloads, frozen registrations and
the pinned mandatory runner are unchanged. These repairs fall under
[#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).

The output is a retrospective, model- and prior-dependent posterior diagnostic.
The rare monotone-subset fractions are reported with per-chain tail counts and
resolution warnings; they are not a branch evidence ratio, a frequentist
exclusion, a direct capacity measurement, an OPH confirmation, or a frozen
FZ-13 score. Agreement with the fixed point is shared with LambdaCDM and earns
no OPH-specific credit.
