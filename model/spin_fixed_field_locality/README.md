# Preparing local marginals at positive selection field

For the actual native spin source and every fixed B>0, a uniform finite
spectral gap and source locality give an explicit volume modulus for
local ground-state marginals. The [proof](PROOF.md) compares a larger
source to the ground of an interior box by turning on only its crossing
bonds. Its error depends on the observation region, boundary area and
distance, with no exterior-volume or Hilbert-dimension factor.

The uniform gap B/8 is an explicit application of Rayudu--Takahashi's
[Lee-Yang theorem](https://arxiv.org/html/2607.10765v1). Gaussian approximate
spectral flow then gives the marginal bound directly. The producer checks
the actual source basis transformation and finite ground derivatives,
and chooses finite box sizes using rational upper bounds. It does not
allocate those large boxes or claim practical simulation costs.

The same bound controls a complete instrument confined to the observation
region, retaining every outcome. Independently attached references are
allowed. General correlated-source inputs and the zero-field ordered-state
modulus require separate estimates. Positive-field preparation does not
supply a particle mass or a relativistic dispersion.

```sh
python -W error -m pytest -q model/spin_fixed_field_locality
python -W error -m model.spin_fixed_field_locality --check model/spin_fixed_field_locality/results.json
```
