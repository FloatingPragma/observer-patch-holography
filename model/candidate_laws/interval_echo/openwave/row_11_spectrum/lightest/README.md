# What is the lightest negative charge on a blank, filled federation?

On each of the three native boards the lowest-energy state with one net missing load is a **single vacancy**. Its ground energy is simple and the vacancy has a positive amplitude at every port. The independently prepared compact vacancy at 35 units sits about 2.4 units above it. These are formation energies relative to the blank filled host, with all registered contacts open and all native degree references present. [Integer certificates](results.json) · [Independent native-seam test](test_lightest.py) · [Compact reference and storage experiment](../../row_04_compact/README.md)

| Cells | Global ground energy in native units | Compact energy minus ground energy |
| ---: | ---: | ---: |
| 20 | 32.53708587–32.53708589 | 2.46291411–2.46291413 |
| 80 | 32.59650254–32.59650256 | 2.40349744–2.40349746 |
| 320 | 32.59603971–32.59603973 | 2.40396027–2.40396029 |

## Why one vacancy gives the global ground

Let $D$ count vacancies plus double occupancies. In the blank sector every free exchange preserves their counts. At net load $-1$, $D$ is odd and $D=1$ is exactly the one-vacancy subspace. The local-word certificate bounds the formation energies at $D=3,5,7$ below by 88, 112 and 144 respectively. For $D\ge9$ the same owner inequality and positive exchange give $H\ge6D-12\ge42$. A vacancy at any degree-five port has mean energy $32+5/2=34.5$. By the variational principle its sector minimum lies below every other composition minimum. The sharp-port mean serves only as a trial upper bound. [Local-word certificate](../bounds.py) · [Composition bounds](../bounds.json)

For one vacancy at port $i$ of degree $d_i$, each incident blank seam supplies the diagonal $8+1/2$ and a hopping amplitude $-1/2$. Each of the $6-d_i$ private reference contacts lowers the formation energy by eight relative to a filled port. Thus the full one-vacancy operator is
$$
H=8\,\mathrm{diag}(2d-6)+(\mathrm{diag}(d)-A)/2
$$
for the ordinary seam adjacency $A$. The measured boards are connected. The matrix $64I-H$ is nonnegative, irreducible and has positive diagonal entries. Perron–Frobenius gives a simple lowest eigenvalue of $H$ with a strictly positive normalized eigenvector $u$. Any symmetry that permutes physical ports and preserves the blank Hamiltonian must fix $u$, since it preserves positivity and normalization. More generally every unitary that preserves this line acts there as a scalar. Two such operations have order commutator $+1$. This one-dimensional ground space cannot carry the pair of anticommuting spatial half-turns of a spin-½ doublet. Degeneracy and exchange in other charged sectors require their own measurements. [Charged doublet in an applied field](../../../../../experiments/interval_echo_angular/charged_orientation/flux_doublet/README.md)

For a positive integer vector $v$, let $L=\min_i (Hv)_i/v_i$. Multiplication of $Hv\ge Lv$ by $u^T$ proves $E_0\ge L$. The Rayleigh quotient $v^THv/(v^Tv)$ gives the upper bound. Both numbers are evaluated exactly as rational numbers on the deposited integer witnesses. A floating-point shifted solve supplies candidate witnesses and plays no part in verifying an interval. [Positive integer vectors](positive_witnesses.npz) · [Exact verifier](measure.py)

## What can a charge detector infer?

A vacancy prepared at any sharp port has strictly positive spectral weight $|u_i|^2$ on this lowest line. A coherent return against the independently characterized, same-charge compact mode at 35 would contain the tabulated difference frequency. Its weight can be small. The separate storage experiment measures the compact mode. The ground energies in the table are operator certificates. No ground-line fit or finite-time resolution is attributed to the storage data. [Compact reference read](../../row_04_compact/README.md) · [Charged reference interferometer](../../../../../experiments/interval_echo_lepton_register/pole_meter/README.md)

Changing the origin of an entire charged sector moves its two formation energies by the same amount and preserves the measured frequency difference. Comparing the energy-35 line, the ground line, and a proposed heavy line therefore requires a common calibration. For three measured charged levels the ratio of energy differences $(E_h-E_l)/(E_m-E_l)$ also cancels a common scale. On the 320-cell witness the heavy sharp state's mean of 120325 gives a bound on its spectral reach. Its association with a rest mass requires spectral weight at that energy and a recoil measurement on the same state. [Excited charged source](../README.md#can-an-excited-charged-defect-evade-the-ground-energy-bound) · [Energy-reference ambiguity](../../row_03_clock/neutral_reference/README.md)

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.lightest.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/lightest
```

Passing `--write` to the verifier reruns the sparse discovery and replaces the integer witnesses and their bounds.
