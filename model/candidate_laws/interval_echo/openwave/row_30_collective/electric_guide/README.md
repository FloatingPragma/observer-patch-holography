# Bulk matter waves on the electric-response guide

A single endpoint pulse on the 44-dimer polarization guide gives 43 distinct bulk frequencies. The 42-dimer guide gives 41. On both, the frequencies and spatial wavelengths recovered from separate counters satisfy the massive lattice dispersion across the whole band. The largest absolute discrepancy is below **0.000008 radians per exchange-clock unit**. The guide, its blank seam memories, charge inventory, endpoint settings and active guards are those of the neutral electric-response task. Its static solution already indicated the spectral form and motivated this transfer test. The temporal and spatial counts in this test are fresh. [Source material](../../row_20_electric/README.md) · [Source choices, comparisons and fixed tolerances](PROTOCOL.md) · [Complete count arrays](counts.npz)

| Polarizable dimers | Bulk waves with accepted held spatial profile | Wave-number interval from the material read | Lattice square root, all waves | Ordinary continuum, k ≤ 0.5 | Matched nearest-neighbour hopping, k ≤ 0.5 | Freely fitted hopping, held interior |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 44 | 43 | 0.07155 to 3.07015 | 43/43 | 5/6 | 6/6 | 4/41 |
| 42 | 41 | 0.07486 to 3.06682 | 41/41 | 5/6 | 6/6 | 4/39 |

Each ratio counts absolute errors at or below the fixed 0.001 tolerance. All 84 frequencies use the same two local calibrations for their own guide. The ordinary continuum and a parabola tangent at rest each pass the first five waves on each guide. Nearest-neighbour hopping with the same rest value and curvature passes all six long waves on each, giving an alternative explanation in that window. Across the full band its largest error is 0.06334. Even when the hopping offset and coefficient are fitted to the two measured band endpoints, the largest held error exceeds 0.01529 on both guides. Short waves distinguish the lattice predictions. The continuum approximation's maximum errors there exceed 2.52. Each guide is a fixed finite instance of the guarded material. [Frequency and wavelength reader, all four comparisons](measure.py)

### Source, clock and spatial read

The independent electric experiment selects the guide's actual native ports on the 320-owner board. Each dimer contains one double and one vacancy with their total held fixed. Both neutral endpoint loads equal one. The service shutters the outside exchanges and the guards prevent transfer between dimers. At the terminal read, X, Y and Z are the two transverse orientation quadratures and the sharp orientation of this pair of rails. Basis rotations and collecting individual dimer signs are part of the instrument. On a mixed material preparation, the first and final readout signs are multiplied. The source for the spatial scan measures just Y on the first dimer, with both first outcomes retained. At receiver i the observable is a signed product of orientation quadratures along the guide. Its support extends from the launch end to i. The spectrum uses a separate ensemble with matched first and last probes alternated uniformly among sites. No wavelength or normal-mode launch is supplied. The observer's time unit is the phase of an isolated matched dimer, and the appointments use a shared service schedule. Position uses distances along the prepared guide in dimer spacings. [Original stock, routings and shutters](../../row_20_electric/PROTOCOL.md) · [Neutral boundary settings](../../row_20_electric/neutral.py)

The return spectrum has a resolved non-oscillating component of 0.02273 on the longer guide and 0.02382 on the shorter. The frequency reader measures every resolved positive peak. The spatial reader infers each profile from a separate array of endpoint-release counts, estimates k on two of every three sites, and predicts the signed wave at the held sites. The worst held RMS discrepancy is 0.371% of peak amplitude, against a 5% limit. In the 44-dimer exchange-only control the reader reports 43 positive peaks. Forty-two pass the spatial check. The final peak at 1.49867 has a held RMS error of 19.7%. A separate spectrum check places two close lines there. All 42 resolved profiles match the control's own lattice forecast. Suppressing transverse exchange inside every dimer yields just one nonzero frequency, 32.5, with no accepted spatial profile. The changed law, positive tests and control share the same acquisition and decoder. The numerical tolerances specify a finite comparison. A joint confidence interval for all fitted spatial and temporal parameters is outside this analysis.

### Microscopic relation

The two guarded orientations at equal native degree have a common degree energy. In a staggered gauge write $z_i=(-1)^iZ_i$ and $y_i=(-1)^iY_i$. Neutral endpoints contribute constants. Up to a scalar the native guard-compressed generator is $H=-h\sum_i X_i-J\sum_i z_i z_{i+1}$. A separate dimer clock and a load-compensated guarded connector give $2h$ and $2J$ respectively. On the 44-dimer guide they give $m=2|J-h|=31.4999998$ and $c^2=4Jh=32.4999957$ in the common units.

Let $A_i=(\prod_{j<i}X_j)z_i$ and $B_i=(\prod_{j<i}X_j)y_i$. Their linear Heisenberg equations give $\ddot B=-C^TC B$, where $C=2(hI-JS)$ and $S$ shifts one site forward. The two-time mixed-state sign correlation is the corresponding element of $\cos(t\sqrt{C^TC})$ since the strings are orthonormal under the normalized trace. In the material interior this yields

$$\ddot B_i=-m^2 B_i+c^2(B_{i+1}-2B_i+B_{i-1}),\qquad
\omega^2=m^2+4c^2\sin^2(k/2).$$

At long wavelength the correction to $\omega^2=m^2+c^2k^2$ starts at $-c^2k^4/12$. This is a gapped collective orientation response on the guarded guide. The mode frequencies come from the time counts, the wave numbers from the spatial counts and the two coefficients from the separate bench. An independent direct four-dimer calculation forms the 70-word microscopic space, compresses it to the 16 permitted rail words, constructs every string as a dimer observable, and checks its cross-site evolution at four exposures. [Direct native-word regression](test_measure.py)

### Releasing the dimer guards

On four dimers, removing the dimer-number guards while keeping the two endpoint loads and exterior fixed permits all 70 placements of four doubles on eight free ports. At exposure 64, a separate stock check found one double per dimer on **99.25%** of launches. The signed endpoint return was **−0.18459** versus **+0.00592** in the guarded experiment. Invalid endpoint orientation reads give zero in the signed average. A high probability of finding the original stock at readout thus allows a measurable change of phase. The four-dimer release gives no long-guide wave-number spectrum. Extending the freely evolving material test across lengths, or implementing the guards on the candidate with an independent in-flight check, would address this dependence. [Complete physical release and count decoder](release.py)

**One passive connector cannot enforce this guard and retain its ordering energy.** On each of the four endpoint words 00, 02, 20, 22, the ten-dimensional written-memory block of the native seam acts on memory alone. For 02 and 20 the words lie outside the write image. On 00 and 22 the image comes only from equal input loads, with zero contrast and exchange energy. On the blank block the unlike word has diagonal cost $4\kappa+1/2$ and swap amplitude $-1/2$, while the equal word has zero cost. The swapped unlike word changes the total load in both adjacent dimers. The blank and written blocks are orthogonal reducing sectors of the free generator. A common connector state of blank weight $p$ therefore gives ordering coefficient $J_{\rm eff}=p(2\kappa+1/4)$ from its initial work contrast, and probability $p\sin^2(t/2)$ of a stock violation for an isolated unlike endpoint word at exposure $t$. Blank/written coherence cannot cancel the swap. At κ=8 retaining $J=16.25$ requires $p=1$, and exact free guard protection requires $p=0$. This result concerns one passive connector on these endpoint loads. More contacts, encoded dimers and time-dependent controls define other experiments. The exact symbolic check covers all 8,712 output coefficients for both generator parts on the 44 enrolled input columns. [Exact native-column check](connector_guard.py)

### Can a distant dimer receive a local signal?

A separate intervention chooses the sign of X on a single central dimer and releases it into the same static 44-dimer guide. Each receiver reads X on its own pair of rails. At exposure 12 the signed responses six dimers away are +0.03648 and +0.03663; ten dimers away they are +0.02002 and +0.01997. All four readings at exposure ¼ agree with zero, as do all four remote receivers at every appointment with the launch connectors cut or the transverse flips suppressed. The simultaneous 99% counting radius across all 1,584 signed averages is 0.000615 under the simulated law. The strictly local detector establishes transmission of the intervention in guarded material. Its bilinear response combines waves. The single-mode frequencies and wavelengths in the table come from the separate string instrument. [Enrolled local test](LOCAL.md) · [All receiver counts](local_counts.json) · [Local sign decoder](local_signal.py)

Run the measurements from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_30_collective.electric_guide
python -B -m "$pkg.measure" --counts-only
python -B -m "$pkg.measure"
python -B -m "$pkg.local_signal"
python -B -m "$pkg.connector_guard"
python -B -m pytest -q --pyargs "$pkg.test_measure"
~~~
