# Does a single record service supply enough cold matter for the acoustic polarization?

A single service on the tested blank federation can supply at most **3/7** of the native incident energy under external matching. Even loosening the bound to **1** leaves large discrepancies in Planck's temperature–polarization cross-spectrum (TE) and E-mode polarization power (EE) when the incident stock is identified with cosmological baryons and all the supplied energy behaves as cold matter. On 310 polarization bandpowers held out of the temperature calibration, the covariance-weighted residual is **17,641** at 3/7 and **11,410** at 1. A standard cold-to-baryon control at 5.36 gives **293**. [Native supply calculation](../../../model/experiments/interval_echo_shell_flux/product_address/green_response/material_recoil/overlap_source/native_boundary/population/results.json) · [Direct comparison figure](polarization.pdf) · [Numerical results](results.json)

A subsequent comparison uses repeated, paid services on three specified board sizes. Their exact cold/baryon inputs are 2.104, 5.134 and 10.520. Only the middle level reproduces the overall measured polarization error scale, with 307.30 on 310 bands. It requires 259.25 units of external work on an initial one-load energy of 50.5. [Three commissioned trails and their sky comparison](TRAILS.md)

## What goes into the prediction

Write \(q=\omega_c/\omega_b\), where \(\omega_i=\Omega_i h^2\). To convert finite-board energy into acoustic observables we grant the gained budget the entire cosmological response of standard adiabatic, pressureless, universally gravitating cold matter by recombination. We identify the initial native energy density with the entire baryon density. The stock, its formation before recombination, its perturbations and this energy conversion require separate OPH derivations. A primordial reservoir, repeated supply, a different gravitational response or a different perturbation law would change what can be inferred from the board bound. [Native temperature comparison](../cold_dilution/native_boundary/README.md)

CAMB 1.6.6 evolves flat models with \(H_0=70\) km/s/Mpc, \(\omega_b=0.02237\), optical depth 0.055, a 0.06 eV total neutrino mass, standard radiation and BBN helium. The scalar spectrum starts at amplitude \(2.1\times10^{-9}\) and index 1. We use lensed angular spectra computed from the linear matter transfer. The nominal baryon density and standard control use the published Planck TT,TE,EE+lowE+lensing benchmark, whose ancestry includes polarization. The source caps arise independently of Planck. The controls test the transfer and the scoring method. [CAMB source](https://github.com/cmbant/CAMB)

We used the whole-bin interval fixed before accessing the polarization bandpowers. It retains 155 TT, 155 TE and 155 EE bands over multipoles 100–1494. Each model uses **TT only** to calibrate a common amplitude \(A\), angular dilation \(d\) and tilt \(t\), mapping \(D_\ell\) to \(A(\ell/500)^tD_{\ell/d}\) in all three channels before binning. The nominal ranges are \(0.55\le d\le1.60\) and \(-0.20\le t\le0.20\). These nuisance parameters approximate smooth changes in angles and power. The table also reports a wider tilt range tried after the nominal fits touched its upper endpoint. [Protocol and amendments](PROTOCOL.md)

TE and EE observe the same sky as TT. For polarization \(P=(TE,EE)\), residuals \(r=y-\mu\), and the Planck bandpower covariance \(C\), we use
\[
r_{P|T}=r_P-C_{PT}C_{TT}^{-1}r_T,\quad
S=C_{PP}-C_{PT}C_{TT}^{-1}C_{TP},\quad
Q_{P|T}=r_{P|T}^{T}S^{-1}r_{P|T}.
\]
Correlated TE and EE contributions need their full joint covariance. Plik-lite marginalizes foregrounds and incorporates Planck's instrument and polarization calibration. Its covariance uses a fiducial sky. Scores far from the fiducial model are discrepancies on that error scale. A sampling significance for each extreme alternative would require its own covariance and nuisance treatment. [Planck high-multipole likelihood](https://arxiv.org/abs/1907.12875) · [ESA likelihood data](https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz)

## Temperature-calibrated prediction

The first score uses 155 training bins. The next two use 155 held bins each. The last two columns use the same 310 held bins. The broad-tilt sensitivity allows \(-1\le t\le1\). The first column lists literal model inputs, and \(z_{\rm eq}\) is CAMB's derived redshift of matter–radiation equality.

| Cold/baryon input | \(z_{\rm eq}\) | TT | TE given TT | EE given TT | Joint given TT | Joint, broad tilt |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 534 | 14,065 | 12,088 | 13,763 | 24,066 | 18,383 |
| 17/260, previous restricted stock | 569 | 13,232 | 11,531 | 13,140 | 22,976 | 17,668 |
| 17/82, pooled classical federation | 645 | 11,578 | 10,408 | 11,859 | 20,749 | 16,203 |
| 17/80, pooled coherent federation | 647 | 11,522 | 10,369 | 11,815 | 20,672 | 16,152 |
| 1/4, ordered sharp stock | 668 | 11,121 | 10,092 | 11,493 | 20,115 | 15,785 |
| 3/7, external matching on the federation | 763 | 9,383 | 8,867 | 10,059 | 17,641 | 14,138 |
| 1/2, sharp overwrite of classical inputs | 801 | 8,762 | 8,418 | 9,528 | 16,727 | 13,523 |
| 17/32, sharp overwrite of coherent inputs | 818 | 8,502 | 8,229 | 9,302 | 16,341 | 13,261 |
| 1, unrestricted coherent single-service bound | 1,069 | 5,393 | 5,839 | 6,419 | 11,410 | 9,849 |
| 2, control above both main caps | 1,604 | 2,166 | 2,745 | 2,677 | 5,008 | 4,977 |
| 0.12/0.02237, standard control | 3,403 | 154 | 154 | 138 | 293 | 293 |

These are fixed points on a grid. The prespecified windows show both the large errors and one accurate TE prediction at a cold/baryon ratio of 1. The signs matter for TE. Entries are means of \(D_\ell\) in µK². Window weights use inverse marginal variances, and errors use the full within-window covariance.

| Channel and multipole | Planck | 3/7 | 1 | Standard control |
|---|---:|---:|---:|---:|
| TE 120–170 | −44.73 ± 1.31 | −28.94 | −29.80 | −43.94 |
| TE 270–330 | 112.25 ± 1.79 | 152.05 | 141.41 | 111.71 |
| TE 400–480 | −57.58 ± 1.20 | −105.22 | −96.54 | −55.29 |
| TE 730–810 | −111.85 ± 1.44 | −107.98 | −111.86 | −115.65 |
| EE 270–330 | 7.70 ± 0.13 | 14.07 | 12.68 | 7.74 |
| EE 400–480 | 16.82 ± 0.21 | 25.25 | 23.35 | 16.72 |

At 3/7 and 1, respectively, the wide TT-only fits select \(t=0.424\) and 0.330. Propagating the three local TT calibration errors reduces the nominal joint scores to 13,742 and 9,225. This linearized treatment is approximate at the nominal tilt boundary. With polarization efficiency fixed in turn to 0.95, 1 or 1.05, the respective score ranges are 12,708–24,892 and 7,885–16,980. Recomputing the transfer at baryon densities 0.0215 and 0.0235 gives 16,535–18,603 and 10,352–12,340. These are separate sensitivities. [Full scores and calibration parameters](results.json)

## Published physical parameter fits

Planck VI Table 2 varies all six base ΛCDM parameters under the collaboration's priors. The TE+lowE analysis gives \(\omega_b=0.02249\pm0.00025\), \(\omega_c=0.1177\pm0.0020\). EE+lowE alone gives \(0.0240\pm0.0012\) and \(0.1158\pm0.0046\). Both include low-multipole E-mode polarization. Neither includes the high-multipole TT likelihood used for our calibration. Their data products and the comparison above share the same underlying observation. [Planck VI, sections 2.2.4–2.2.5 and Table 2](https://arxiv.org/abs/1807.06209)

We also computed the ratio at each joint posterior draw in all four released chains for each polarization combination. The release marks burn-in as removed. We retained the original multiplicities and calculated equal-tailed quantiles after marginalizing over the sampled nuisance and cosmological parameters.

| Deposited fit | Median cold/baryon | Central 68% | Central 95% | Smallest sampled | Chain rows |
|---|---:|---:|---:|---:|---:|
| TE+lowE | 5.235 | 5.107–5.359 | 4.986–5.484 | 4.732 | 29,884 |
| EE+lowE | 4.835 | 4.426–5.265 | 4.071–5.737 | 3.520 | 28,587 |

The finite chains contain zero weight below either cap. The quoted intervals apply to base ΛCDM under Planck's likelihood and priors. The declared density prior ranges are 0.005–0.1 for baryons and 0.001–0.99 for cold matter, which admit both 3/7 and 1 at the nominal baryon density. The authors' autocorrelation-corrected effective counts for the two density parameters range from 2,531 to 3,016. The four per-chain ratio medians span 5.226–5.239 for TE and 4.818–4.853 for EE. These separate fits share lowE and sky coverage. [Planck mission chain archive](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/) · [Exact joint-sample results](chain_results.json) · [Chain members and hashes](chain_provenance.json)

The supply shortfall has a quantitative consequence. If each service converts at most 3/7 of its incident energy, the whole retained output is additive and gravitates as this cold component, then the cumulative incident energy per baryon energy must be at least the inferred cold/baryon ratio divided by 3/7. At the lower endpoint of the EE central 95% interval this requires **9.50 baryon-energy units**. Regardless of whether the supply comes from repeated services or extra material, a cosmological model must specify its origin, timing and surviving perturbations to predict these spectra.

## Reproduction

Run from the repository root in a Python environment with CAMB, NumPy, SciPy, Matplotlib and pytest installed.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -m cosmology.cold_acoustic_polarization.measure
python -m cosmology.cold_acoustic_polarization.figure
python -m cosmology.cold_acoustic_polarization.chains
python -m pytest -q code/cosmology/cold_acoustic_polarization
```

The packaged arrays retain all 613 original bandpowers, band windows, covariance and the ESA test spectrum. To rebuild them from the primary release use `python -m cosmology.cold_acoustic_polarization.fetch`. Its seven input member hashes are recorded. Evaluating ESA's own complete test vector yields −292.286311296 against the official log likelihood −292.2863113. The tests also verify the band windows, the exact joint/conditional covariance decomposition, the separation of training data and a fresh 3/7 CAMB replay. [Data provenance](provenance.json)

The packaged chain columns contain the original multiplicity and the two physical densities, separated by chain. To rebuild them, download the 471 MiB base Plik archive and run `python -m cosmology.cold_acoustic_polarization.chains --archive /path/to/COM_CosmoParams_base-plikHM_R3.01.zip`. The replay compares the two weighted means and standard deviations with the archived marginals. [Primary ZIP](https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/cosmoparams/COM_CosmoParams_base-plikHM_R3.01.zip)
