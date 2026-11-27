# Six loads bind below the two-trimer threshold

A time-port tuned to the calibrated low-energy line admitted **66,953,343 of 67,108,864** oriented six-load caps on one blank cell. It admitted **10** preparations with the same six loads on triangular faces in two independent cells. Only 611 admitted caps occupied two disjoint faces. On the three registered federations the fully released sharp cap retains a connected occupied support with probability at least **7/8 at every time**. The measured axial quadrupole of the oriented material is negative about its equal-mass centroid. Independent proton and neutron identities, a polarized nuclear spin state and its electric moment are unestablished. The signed-source follow-up finds a positive material quadrupole on the rare launches where distinct regional probes record +1 and zero below their separated threshold. A cheaper competing fragmentation keeps the full deuteron criterion open. [Signed-source experiment](signed_pair/README.md) [Official criterion](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#simplest-test-per-criterion) · [Prospective protocol](PROTOCOL.md)

## Comparing energies and reading both shapes

Each launch starts with ordinary binary loads and blank seams, as in the [three-load face experiment](../row_13_baryon/README.md). The states use the same absolute native exchange and contrast energies, including the private reference cost at every port. The cell radius, polar axis, fresh zero references, coherent energy comparator, dwell and preparation are supplied. A signed echo from each occupied port reads +1 against its reference. Distinct readout cohorts contain sharp loads 0, 1 and 2. Their signed difference records give the enrolled answers on every calibration start. No fitted energy, length or electric standard from the deuteron entered the run. The energy reader used the six-load ground line rounded to 132.924 in the particle paper's binding table. The line location was known before acquisition. [Native readout and timing](apparatus.py)

| Source in the native absolute-energy filter | All starts | Bright | Bright and two disjoint faces | Bright and a cap |
| --- | ---: | ---: | ---: | ---: |
| Polar cap, then immediate terminal surveys | 67,108,864 | 66,953,343 | 611 | 66,798,388 |
| Polar cap, eight later exchange units | 67,108,864 | 66,953,190 | 618 | 66,798,461 |
| Two disjoint sharp faces in one cell, immediate survey | 67,108,864 | 210 | 82 | 118 |
| The same two-face source, eight later exchange units | 67,108,864 | 201 | 86 | 105 |
| Two isolated cell registers with one sharp face apiece | 67,108,864 | 10 | Separate registers | Separate registers |

The complete banks contain the port word on every bright launch and one dark refusal count. The separated run retains both labelled three-load words on a bright launch. The arms compare 8,192 coherent delays. Randomizing the delay coherence gives 8,346 bright in the native control. With exchange alone the common zero-energy comparison admits both the cap and separated stocks, 74,286 and 1,918 bright respectively. Their differing overlaps with the ground state survive the loss of the separation gap. The independently calculated sector ground energies are zero for all six occupancies in that control. Reversing contrast raises the six-load ground to −201.885 against −276.921 for two trimers and −288 for six independent loads. The reversed cap can persist geometrically at high energy. A terminal shape by itself would miss the loss of binding. [Raw counts, controls and spectral comparison](records.json.gz)

All counts are seeded draws from the calculated finite Born probabilities on the declared probability grid. At time zero, with simultaneous 99% exact binomial coverage over the five reported comparisons, the cap's bright probability is at least 0.997664, and the separated-trimer probability is at most 0.000000360. The fraction of bright caps with two faces lies between 0.00000802 and 0.00001033. The two-face input gives bright probability between 0.00000250 and 0.00000386, with two faces on 28.89–49.92% of bright outputs. These intervals describe the sampler. Rounding to the 40-bit probability grid changes any event probability by at most the number of categories divided by the grid size for these banks (below 4.5 × 10⁻⁸ for the largest bank). The six-load cap passes the prespecified separation and leakage gates. The rare events from the two-face source include leakage through the timing aperture. Their conditional shape distribution conflicts with an assignment of all bright records to the cap band. [Complete distributions](measure.py) · [Subspace and aperture bounds](PROOF.md)

## How much of the low-energy cluster can contain two triangular cores?

Each separated three-load trimer has energy at least 96. The calculated ground is 100.27234 on the single cell. The six-load ground is 132.92384. Integer edge counts and the positivity of exchange give a six-load variational upper bound of 134 on each open federation, against at least 176 for *every* proper separation of six loads and at least 192 for the two-trimer channel. These are native exchange units. The 20-, 80- and 320-cell bound concerns full free evolution. The sampled terminal banks concern one isolated six-load cell or two such independent trimer cells.

In the single cell a twelve-state band lies between 132.92384 and 132.92419. On every coherent or mixed state in that band an ideal position survey has probability at most **135/29,282**, or 0.4611%, of finding two disjoint native triangular faces. Direct diagonalization gives the sharper numerical maximum 0.001840%. The exact bound uses only the spectral gap of at least eleven above the cap compression and two integer contact-degree counts. The nearest two-face word needs at least two exchanges from any cap. A bound load aggregate need not retain two constituent three-load faces. [Integer certificate and proof](PROOF.md)

## Which moment does the polarized preparation carry?

Observers reconstruct the traceless quadrupole from all twelve terminal port readings, in material-charge times squared cell-radius units. The polar axis and cap orientation are fixed before the energy measurement. Sharp alignment standards give +2.4 and −2.4 about the cell origin. For the bright cap, the axial mean about that origin is within 0.000001 of zero. Complementing all loads and inverting the sphere fixes the cap, commutes with its evolution and energy filter, and reverses the cell-centred quadrupole on every axis. The exact zero applies to this oriented preparation. The two-face preparation gives +2.4 on the second surveyed axis at launch.

At the centroid of the six equal masses, the same bright cap words give **−3.48669** on the chosen axis at immediate arrest and **−3.48667** eight exchange units later. The independently aligned positive calibration gives +2.88240 at its mass centre. The measured sign on the cap disagrees with the positive deuteron value, and the two inputs are identical load triplets. These experiments have no independent proton or neutron rest-mass, spin-polarization or electric-charge calibration. [Sharp standards and full banks](records.json.gz) · [Geometry of the two centres](PROOF.md)

From the repository root, replay with

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_16_binding
.venv/bin/python -B -m "$pkg.measure" --output /tmp/row16-binding.json.gz
.venv/bin/python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_16_binding
```

The tests rebuild the six-load matrix from the 99-state seam, compare an independently factored binary delay tree to the terminal counts, and check the block-transfer inequalities by direct word enumeration.
