# Measuring geometry fluctuations

A complete coordinate measurement supplies a covariance estimate and its
energy cost. An exact finite-sample bound needs a number of independent
samples that does not grow with the geometry register count; each sample's
work grows with that count. Converting the same source's rest and photon
responses exposes an additional normalization factor. Physical metric
noise and the native-process realization are separate requirements.

The [proof](PROOF.md) states the full instrument and common-unit map.
Run `python3 -m model.openwave_geometry_covariance_protocol --write` and
`python3 -m model.openwave_geometry_covariance_protocol --check`.
