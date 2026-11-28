# The recoil of a moving polarization defect

The independent fractional-polarization experiment fixes a ten-port belt and two mobile walls on it. One defect has three ordinary loads and polarization −⅓, the other four loads and polarization +⅔. On both preparations the dominant temporal lines resolve at every measured wavelength. **All twelve held comparisons at wavelengths shorter than the calibration exclude the massive continuum dispersion**, for every positive mass and wave speed, on the belt's measured perimeter coordinate. The nearest-step lattice version also misses all twelve. The two reversals at the calibration wavelength agree. [Polarization experiment](../../quark_winding/belt_walls/README.md) · [Design, enrollment and prior exposure](PROTOCOL.md)

Each source superposes a uniform wall amplitude with one rotating character. At every appointment a sharp load read gives the position of the anomalous gap or an outcome outside all ten wall words. The initial position counts identify the wavelength. Across appointments their signed spatial Fourier component measures the phase advance relative to the uniform character. The complete generator evolves all 120 or 210 load configurations. A separate spectral calculation checks whether a single energy carries at least 95% of each character's weight. It supplies a conservative tail bound for the frequency reader. On the native contrast law all audited weights exceed 98.69%. With contrast switched off they range from 29.30% to 58.88%, leaving all fourteen held comparisons unresolved. [Executable acquisition and count decoder](measure.py) · [All 448 count panels and source audits](counts.json)

| Harmonic on ten ports | Three loads, forward recoil | Four loads, forward recoil |
| --- | ---: | ---: |
| +1, calibration | 1.2058 to 1.2068 | 1.2175 to 1.2186 |
| +2 | 1.7010 to 1.7020 | 1.7353 to 1.7364 |
| +3 | 0.1926 to 0.1937 | 0.1678 to 0.1688 |
| +4 | 0.6175 to 0.6186 | 0.6271 to 0.6281 |

Intervals are rounded outward, in radians per supplied exchange-time unit. They include the fixed spectral-tail allowance and a simultaneous counting allowance of at least 99% conditional on the simulated laws. For measured first recoil $d>0$, every massive continuum curve predicts a recoil between $3d$ and $9d$ at three times the wave number. The lower ends of the two predicted ranges are 3.6174 and 3.6527. A shared calibration of the clock or a change of the energy zero preserves the discrepancy. All signs, individual frequency intervals, predicted envelopes, and unresolved results are accessible with the count decoder below.

To obtain the range, put $x=c|q_1|$, write $m$ for the rest frequency, and rationalize the difference of square roots. The ratio at $q=rq_1$ is $r^2(m+\sqrt{m^2+x^2})/(m+\sqrt{m^2+r^2x^2})$. For $r\ge 1$ the denominator lies between $m+\sqrt{m^2+x^2}$ and $r(m+\sqrt{m^2+x^2})$. Thus the ratio lies between $r$ and $r^2$, including for a band with the overall energy sign reversed.

The anomalous gap moves three perimeter ports in one elementary wall hop, as the independent transport task records. Within the subspace of the ten sharp walls the compressed generator has a constant diagonal and hopping amplitude −½ to each of the two addresses displaced by three. Its Fourier symbol therefore includes $-\cos(3q)$. Excursions outside that subspace change the measured curve. The ten-port belt's longest nonconstant wavelength has $q=2\pi/10$, or $3q\simeq1.885$ radians per direct wall hop. This apparatus cannot resolve the long-wavelength limit relative to that hop. The result excludes the tested dispersions on the enrolled harmonics, with the uniform character as the rest preparation. A test near that rest point requires a larger independently specified material and a clock that can resolve its smaller recoil. [Independent native-incidence propagation and projected hopping check](test_measure.py)

To regenerate all counts, decode the archive, or run the regression checks from the repository root, use

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_30_collective.belt_recoil
python -B -m "$pkg.measure"
python -B -m "$pkg.measure" --counts-only
python -B -m pytest -q --pyargs "$pkg.test_measure"
~~~
