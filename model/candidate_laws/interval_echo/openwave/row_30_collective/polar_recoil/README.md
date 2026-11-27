# Does the neutral island propagate as a Klein–Gordon wave?

The neutral-particle experiment supplies a mobile gapped excitation on a passive material ring. Its resting cost and first measured recoil fix the free massive dispersion with no remaining speed fit. On two physical rings, the ordinary Klein–Gordon forecast passes one of fifteen further frequency comparisons. The lattice square-root relation passes three, as does a nearest-neighbour hopping band. The distinction between the two latter predictions is small enough to require long material phase records. [Prospective comparison and prior exposure](PROTOCOL.md) · [Complete counts](counts.json) · [Same-stock mass and charge experiment](../../../../../experiments/material_inertia/pair_spectrum/neutral_guide/polar_ring/README.md)

On the 80-cell board, observers prepare the uniform single upper dimer and one rotating character together. They read each possible single-upper word and one complementary outcome. The largest nonconstant component of the initial spatial histogram supplies the measured wavelength. Fresh copies at 16 nonzero waits determine the direction and phase advance. Both directions of all four measured harmonics pass the spatial identification check on the ten- and twelve-dimer rings. The preparation and calibration use specified local gates and a material clock. The rotating source is freely released.

Let $A_j(t)$ be the return amplitude of the single-upper Fourier word with character $j$. Translation symmetry gives the one-upper probability at site $x$ as $|A_0(t)+e^{2\pi ijx/n}A_j(t)|^2/(2n)$. The spatial Fourier component at the resolved nonzero harmonic therefore reads the relative phase of these two returns. The retained harmonics avoid the Nyquist bin. Complementary terminal words contribute zero to this moment and remain in the original denominator. We checked the translation reduction against direct evolution of every orientation word at two wavelengths. A normalized return with dominant line weight $w>1/2$ differs in phase from that line by at most $\arcsin[(1-w)/w]$. We propagate both line bounds and simultaneous 99 percent sampling bounds through every time step, requiring a unique allowed phase winding. A coarse existing band estimate supplies the first winding interval of one exchange unit in either direction. The final wait is 32,768 exchange times. [Acquisition and decoder](measure.py) · [Independent word-space check](test_measure.py)

The sharp-host clock gives a rest cost in [39.9131416, 39.9131972]. The first forward spatial recoil on ten dimers lies in [0.00945271, 0.00945392]. Each of the three forecasts uses these two readings. For wave number $q$ and calibration wave number $q_0$, set $r=(q/q_0)^2$ for the ordinary relation or $r=\sin^2(q/2)/\sin^2(q_0/2)$ for the lattice. If $m$ is the rest cost and $d$ the first recoil, the predicted recoil is $\sqrt{m^2+d(2m+d)r}-m$. The hopping forecast is $dr$ with the lattice value of $r$. The decoder propagates the full calibration intervals. At fixed positive $d,r$ the square-root forecast is monotone in $m$ because its derivative has the sign of $r-1$. At fixed $m,r$ it increases with $d$.

The greatest separation between the predicted and measured intervals was 0.061169 exchange units for the continuum Klein–Gordon forecast, 0.000347 for its lattice version and 0.000257 for nearest-neighbour hopping. All three use the same two measured calibration values.

The fifteen comparisons include the reverse of the calibration wavelength. They have shared calibration uncertainties. The result excludes these exact dispersions on the tested material and momentum domain. It leaves the near-rest effective mass and the dispersion of a different domain-wall preparation as separate measurements.

Run from the repository root with numerical thread counts set to one.

```sh
pkg=model.candidate_laws.interval_echo.openwave.row_30_collective.polar_recoil
.venv/bin/python -m "$pkg.measure"
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_collective/polar_recoil
```
