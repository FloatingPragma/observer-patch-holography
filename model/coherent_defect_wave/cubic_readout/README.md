# Cubic translation and frequency records

The shared-link occupation source runs on a declared three-dimensional
cubic circuit. Its exact one-wave encoding retains every changed local
defect value and all positive-Gauss degrees of freedom. No dense product
of gauge-pair qubits is required.

The [proof](PROOF.md) and [receipt](results.json) give an actual full tick
on side-four and side-eight cubes, with 9344 amplitudes and defect movement
probability `1/450`. The complete joint instruments have 32 and 256
translation/quadrature outcomes, respectively. Their first moments obey

`c(q)=(-29/30) product_i[(1-i)/2+(1+i)cos(q_i)/2]`.

This function is an expectation of the tick, not a particle dispersion.
A comparison with fixed links has the same first-moment frequency records
and zero defect movement. The separate face readout detects the difference.
The result identifies no physical particle or rest mass.

```sh
python -m model.coherent_defect_wave.cubic_readout --check
python -m pytest -q model/coherent_defect_wave/cubic_readout
```

`lattice.py` supplies the exact positive-Gauss one-occupation source API.
`fourier.py` compresses translation orbits by a proved isometry.
`readout.py` computes complete records without fitting a spectral curve.
The owning [particle path](../PARTICLE_PATH.md) states the physical attachment
requirements and distinguishes other programs and clocks.
