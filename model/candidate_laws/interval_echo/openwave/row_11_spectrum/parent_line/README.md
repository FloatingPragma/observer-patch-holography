# Does the many-defect source prepare a heavy lepton line?

On the 320-cell federation, the 1,979-defect source has mean formation energy 120,325 and standard deviation 29.4364 in native exchange units. Its energy distribution gives at most **0.04512%** probability to a window of one unit on either side of the proposed heavy energy if the 35-unit compact vacancy is the light particle. Taking the lowest negative-unit state as the light reference instead gives at most **0.001783%** in the corresponding window. Both bounds concern this source. Other preparations in this large material sector can have different spectra. [Parent recipes](../excited_bounds.json) · [Moments and bounds](results.json)

## The two charge preparations

The source puts a balanced six-vacancy, six-double cut word in each selected owner and restores one double to an ordinary unit in the first owner. Owners are ordered by distance from owner zero, breaking ties by owner index. The two archived recipes select the first 10 and 165 owners. All federation seams open for free evolution, with blank memory and a filled host outside the preparation. The middle and heavy ratios 206.770316 and 3477.472837 come from the separate first-order register. Their use here sets the energies to test. [Preparation recipe](../excited_bounds.py) · [Ratio inputs](../../electron_field/lepton_family/arrival_spin/anchor_spectrum.json)

| Preparation | Unequal seams | Mean formation energy | Energy variance | Standard deviation / mean | Curvature of return probability |
| --- | ---: | ---: | ---: | ---: | ---: |
| 119 defects | 218 | 7,293 | 54.5 | 0.101226% | −109 |
| 1,979 defects | 3,466 | 120,325 | 866.5 | 0.0244641% | −1,733 |

A zero-load degree reference and the energy of the same filled host fix the energy convention. The one-vacancy compact vector on owner zero has integer amplitudes (0,0,0,0,2,−1,−1,2,−1,0,0,−1). The verifier checks every row of the 3,840-port equation $2Hv=70v$, including all three external gateways. The independently certified negative-unit ground energy lies in [32.59603971, 32.59603973]. The sharp degree-five vacancy used for the variational upper bound has mean 34.5 and variance 1.25. Each choice of a physical light branch is a separate hypothesis. [Exact lightest-state bounds](../lightest/README.md) · [Compact reference measurement](../../../../../experiments/interval_echo_lepton_register/pole_meter/README.md)

| Assumed light state | Middle target energy | Upper spectral weight within ±1 | Heavy target energy | Upper spectral weight within ±1 |
| --- | ---: | ---: | ---: | ---: |
| Compact vacancy at 35 | 7,236.961060 | 1.7674% | 121,711.549295 | 0.04512% |
| Negative-unit ground | 6,739.89341–6,739.89346 | 0.017877% | 113,351.84266–113,351.84278 | 0.001783% |

For the compact reference the window is centered on 35 times the six-decimal comparison ratio, and includes the unrounded comparison. For the ground reference the calculation multiplies the ground interval by the ratio interval rounded outward by half a unit in its sixth decimal place, then expands by one native unit on either side. The receipt also gives bounds for the point at the quoted compact-reference ratio and for the complete unexpanded ground-target interval. An energy shift proportional to conserved charge can shift all host-relative formation energies together. A physical assignment needs a reference fixed independently by an observer. Spectral weight in a prepared state alone determines neither particle rest mass nor lifetime.

## What observers can read

Prepare the sharp word, release it, then ask every owner to read its terminal loads. All owners matching the launch word is one outcome. Any departure is the other. The time-zero curvature of the first probability is $-2\sigma^2$. This prediction concerns the port records in exchange time calibrated by a native seam. Simultaneous terminal comparisons can be made on fresh preparations for each delay. The curvature is an exact local-in-time prediction. The numerical short-time check evolves all 81 four-rail words from the 99-state seam law. No large-board time series or heavy-particle counts enter the table. [Measurement code](measure.py) · [Independent small-system evolution](test_parent_line.py)

For a sharp input $x$, subtracting the filled-host energy shifts only the diagonal. Each seam whose endpoints have unequal loads takes $x$ to a distinct swapped word with amplitude $-1/2$ in $Hx$. The degree references contribute only to the diagonal. Thus $\sigma^2=\|(H-\langle H\rangle)x\|^2=N_{\rm unequal}/4$, and the raw second moment is $\langle H\rangle^2+\sigma^2$. Doubling only the native exchange strength multiplies the curvature by four on the same words, ports and clock. The filled sham returns with probability one at every delay. A measured zero curvature on either active source would falsify the predicted local exchange response.

For an energy window a distance $d>0$ to one side of the mean, its spectral weight is at most $\sigma^2/(\sigma^2+d^2)$. To see this, shift the centered energy toward the tail by $\sigma^2/d$ and apply the second-moment Markov bound to its square, reversing the sign for a lower tail. This bounds the whole tail, including any degenerate line in the window. The bound gives no positive lower weight there.

Run from the repository root with single-threaded BLAS.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.parent_line.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/parent_line
```
