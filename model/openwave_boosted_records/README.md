# Moving particle records on the shared source

A positive spatial kinetic term gives the charged branches of the
longitudinal matter model moving, complete observer records. Their exact
rest values are preserved. A portable Taichi kernel executes the full
one-lepton reducing sector on 16, 32 and 64-site probes, including every
native particle/hole detector outcome and preparation failure. Its group
speed is bounded by one, and its continuum dispersion and complete record
errors have explicit quadratic refinement bounds.

The square-root kinetic is nonlocal. Full interacting Lorentz covariance,
physical spatial contraction and transverse photons are separate claims.
See [the proof](PROOF.md) and [the results](results.json).

    python -m model.openwave_boosted_records --write --openwave-root /path/to/openwave --python /path/to/taichi/python
    python -m model.openwave_boosted_records --check --openwave-root /path/to/openwave --python /path/to/taichi/python

The environment variables `OPH_OPENWAVE_ROOT` and `OPH_TAICHI_PYTHON`
provide the same paths for the ordinary `--check` command. The runner
requires the pinned OpenWave checkout and a Taichi interpreter; it has no
substitute NumPy-only replay. The recorded execution uses Taichi 1.7.4 on
Python 3.12.
