# Complete planar wave deflection

All961 planar wave modes plus the vacuum propagate on the central slice
of the native31³ lapse source. The complete position and velocity
instruments retain every outcome and admit arbitrary quantum references.
At source time16 the γ=1 packet bends toward the source by
1.9417333×10⁻⁵ sites relative to flat propagation, with transverse
velocity change −3.6701444×10⁻⁶. A full2209-mode side47 control keeps
the sign. Reflection, γ=0, γ=−1 and zero-profile controls are included.

The local second-order wave equation independently reproduces the
positive-frequency propagation to numerical roundoff. The analytic
Chebyshev tail is bounded separately; it does not include arithmetic
roundoff. Complete source data and instrument probabilities are in
`results.json`; `PROOF.md` gives the finite theorem and scope.

The principal-ray and finite-wave responses are explicitly compared.
γ=−1 has zero ray response but a small finite-wave residual. The declared
static scalar plane is not a three-dimensional photon proof or the
two-cell autonomous quantum geometry. Physical G and lensing angles are
not inferred from these dimensionless finite records.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_wave_deflection --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.openwave_wave_deflection --check
python -m unittest model.openwave_wave_deflection.test_waves
```
