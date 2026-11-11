# Fixed lepton ratios and their comparison budget

At the declared phase `theta=2/9` and balance `c=sqrt(2)`, the positive
register eigenvalues have shape `1+c cos(theta+2 pi k/3)`. Their squares,
assigned to electron, muon and tau in modes 1, 2 and 0, give
`mu/e=206.7703159727293886` and `tau/e=3477.4728371045985323`.
The [complete-cycle proof](PROOF.md) derives the limiting effect rates on
the supplied apparatus. It supplies no interacting physical mass-pole map.

The comparison freezes the direct muon/electron ratio
`206.7682827 +/- 0.0000046` from the
[2022 CODATA table](https://physics.nist.gov/cuu/Constants/Table/allascii.txt),
the electron energy `0.51099895069 +/- 0.00000000016 MeV` from the
[same adjustment](https://physics.nist.gov/cuu/pdf/all.pdf), and the tau mass
`1776.93 +/- 0.09 MeV` from the
[2024 Particle Data Group average, page 1](https://pdg.lbl.gov/2024/listings/rpp2024-list-tau.pdf).
The direct muon ratio consumes its published adjustment uncertainty. The
tau reference is the ratio of the two stated energies. These editions are
fixed comparison inputs; no claim of a prospective blind test is made.

## Experimental uncertainty and empirical exposure

For `r=t/e`, write `a=u_t/e` and `b=t u_e/e^2`. Linear covariance
propagation gives `u_r^2=a^2+b^2-2 rho a b`. Positive semidefiniteness
requires `-1<=rho<=1`, hence the sharp range
`|a-b|<=u_r<=a+b`. The receipt retains both endpoints because the selected
tables give no cross-edition covariance. The algebraic box
`[(t-u_t)/(e+u_e),(t+u_t)/(e-u_e)]` is conditional on its input box;
standard uncertainties alone do not imply that box has a stated coverage.
No joint covariance of the two comparison ratios is supplied, so the
comparison reports no joint chi-square or probability of agreement.

Dividing residuals by experimental standard uncertainties gives
reference-only standardized discrepancies. These cannot be interpreted as
the significance of a complete physical theory: the register-to-pole map
and its errors are unspecified. The small fractional muon discrepancy is
large relative to the experimental uncertainty. Its cause is not derived
by the model. No finite-round correction is selected to cancel it.

The architecture and round-selection functions accept no comparison target.
The exposure record names the source independence, analyzer weights, phase,
balance, squared-rate interpretation and mode assignment. Particle masses
and Koide relations were known during construction of this architecture.
The optional MeV display uses the measured electron energy as an explicit
calibration; it does not derive an absolute energy scale.

## Finite execution and preparation sensitivity

The complete-cycle theorem bounds half diamond distance, with arbitrary
untouched references, by `14/(81N)` for the weak channel and `349/(1296N)`
for the mass channel at declared time one. An accuracy `10^-6` is fixed
independently of experimental residuals. Half is allocated to truncation,
half to the combined channel and independent fresh-preparation error.
The receipt computes the least integer `N` for that allocation and reports
the conditional per-cycle service allowance. This is an implementation
design, not a measurement of hardware accuracy. A bound on channel distance
does not give an additive error on physical masses or their ratios.

A separate preparation theorem applies to a fixed analyzer projector `A`.
For source states `rho,sigma` independent of the register/reference with
trace distance at most `eta`, the compressed effects satisfy
`||E_rho-E_sigma||<=eta`: for every unit register vector `v`, its source
effect `<v|A|v>` lies between zero and identity, and the variational
characterization of trace distance bounds its expectation difference.
Weyl's inequality then bounds each ordered effect eigenvalue by `eta`.
If `eta<e_e` and twice `eta` is below every spectral separation, the labels
are preserved and

```
(e_j-eta)/(e_e+eta) <= e'_j/e'_e <= (e_j+eta)/(e_e-eta).
```

These enclosures concern limiting effect-rate ratios with the analyzer
held fixed. They are not an independent physical uncertainty, nor a bound
for a changed analyzer. The receipt also retains the exact shape Jacobian
for separately declared phase/balance perturbations. Independence from the
continuing register/reference is essential; source marginals alone do not
satisfy the complete-cycle preparation theorem.

```sh
python3 -m model.lepton_cycle.compare --check
python3 -m pytest -q model/lepton_cycle/test_compare.py
```
