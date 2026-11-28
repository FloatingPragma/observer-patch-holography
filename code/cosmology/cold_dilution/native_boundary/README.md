# One sharp record layer and the microwave background

If one layer of overlap records supplies all the early cold gravitating stress, its most favorable finite energy ratio is 17/260 of the ordinary input, before transport losses. Map that entire increment to an adiabatic pressureless component with the same gravitational coupling as baryons. CAMB gives a third-to-second temperature peak height ratio of **0.535**. The measured Planck PR3 peak-bin ratio is **0.980**, and the standard cold-density control gives **0.980**. This tests a specific cosmological identification of the native energy. The finite source itself specifies no recombination history, primordial spectrum or metric perturbation equations. [Native layer bound](../../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/README.md) · [Planck PR3 measured temperature bins](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_PowerSpect_CMB-TT-binned_R3.01.txt) · [Measured and predicted spectra](comparison.pdf)

| Spectrum | Matter–radiation equality redshift | First three peak multipoles | Third / first peak power | Third / second |
| --- | ---: | --- | ---: | ---: |
| Planck measured bins | | 225, 525, 825 | 0.435 | 0.980 |
| Standard cold density | 3,403 | 220, 536, 813 | 0.443 | 0.980 |
| Maximal one-layer density | 569 | 291, 759, 1,194 | 0.210 | 0.535 |
| Zero cold density | 534 | 294, 771, 1,213 | 0.205 | 0.527 |

Use the existing Planck control parameters, including baryonic density 0.02237, physical cold density 0.12 in the control, tilt 0.9649, flatness and the ordinary adiabatic initial law. The upper-bound model has physical cold density 0.001463. At this density radiation still dominates at standard recombination, changing both the driving of the oscillations and their angular separations. All three calculations include lensing from linear perturbations. CAMB's empirical nonlinear HMCode correction failed to integrate at the small cold density and was disabled uniformly. The baseline is calibrated with Planck, so its agreement is a control. [Transfer settings and bandpower calculation](measure.py) · [Planck cosmological parameters](https://doi.org/10.1051/0004-6361/201833910)

A shape control grants each theoretical spectrum a free uniform horizontal dilation and amplitude change to reproduce the first observed peak. At the third measured crest it then predicts **947 µK²** for one layer, **915 µK²** for zero cold density and **2,563 µK²** for standard cold density. Planck measured **2,522 ± 21 µK²** in that bin. The first-peak calibration uncertainty and inter-bin covariance are excluded from the quoted observational bin error. The fractional RMS across the 350–1500 effective multipoles is 39.0% for one layer after alignment and 7.5% for the control. With the unaltered control the RMS is 1.1%. The aligned control moves a smooth theoretical maximum to the centre of a finite measured bin, which can worsen its fit at other multipoles. Heights and residuals are descriptive comparisons at bin centres. A likelihood marginalized over cosmological parameters or perturbation dynamics would be a different test. [Selection and calibration rules](PROTOCOL.md) · [All bin-centre predictions and residuals](results.json)

Run from the repository root using CAMB 1.6.6. The default checks the deposit against a fresh Boltzmann run and verifies the SHA-256 of the archived public table.

~~~sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH=code:. python -Bm cosmology.cold_dilution.native_boundary.measure
python code/cosmology/cold_dilution/native_boundary/figure.py
~~~
