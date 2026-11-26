# A mobile material string with a measured breaking length

Two opposite orientation walls on a registered Interval Echo strand cost **7.996 exchange-energy units per additional dimer** over separations two through five. While only their two boundary bonds are held, the intervening dimers retain their native swaps. At separation six the ground changes to a screened configuration with two additional walls. The relaxed cost stays within 0.000002 of 79.950127 through separation twelve. On a separate release, with all wall pins removed, the first fourteen dimer readings resolve additional wall pairs in the same biased medium. [Frozen experiment](PROTOCOL.md) · [Readout bank](receivers.npz) · [Counts and calculated spectra](results.json)

| Wall separation in dimers | Relaxed energy | Level-0 intact line | Level-0 screened line | Weight of enrolled ground-to-ground transition |
|---:|---:|---:|---:|---:|
| 2 | 48.000000 | 48.000098 | 48.000134 | 1 |
| 3 | 55.995837 | 55.996042 | 80.016889 | 0.99941 |
| 4 | 63.991657 | 63.991828 | 79.950995 | 0.99882 |
| 5 | 71.987477 | 71.987670 | 79.950278 | 0.99823 |
| 6 | 79.950126 | 79.983500 | 79.950253 | 0.99129 (screened) |
| 7 | 79.950127 | 87.979394 | 79.950254 | 0.99102 |
| 8 | 79.950127 | 95.975254 | 79.950249 | 0.99071 |
| 10 | 79.950127 | 111.966938 | 79.950244 | 0.99009 |
| 12 | 79.950127 | 127.958613 | 79.950269 | 0.98947 |

The intact line at six has only 0.000042 weight on the ground transition. The enrolled screened line tracks the true minimum there. The 0.033 separation from the intact line is below one Fourier bin over the longest hold. Separate prepared clocks distinguish their dominant phases, while the calculated spectrum and weights identify the ground. By separation twelve the bright intact return still follows its original slope, even though its energy exceeds the minimum by 48.01. The below-crossover slope fitted on lengths two and three predicts the two reserved lines at four and five within 0.00027 on all three board deposits. Across the nine separations, the enrolled ground lines differ from direct diagonalization by at most 0.00021. Their Ramsey visibility exceeds 0.97. In the degree-off control only the competing screened preparation at separations four and higher falls below 0.75 visibility, reaching 0.00054 on the level-zero stock at twelve. Fits to that unenrolled fringe have no ground-energy interpretation. [Complete visibility, both candidate weights and both clocks](results.json)

The source has one conserved load on each blank dimer. Its orientation points toward either a degree-six gateway or its adjacent degree-five port. The latter pays for a native private reference. Two one-load dimers per owner continue this orientation order through eight owners. Before testing the string, local phase receivers on the three deposits give the connector cost 31.9824–32.0082, the swap gap 0.99957–1.00088 and the orientation cost 7.9964–8.0069. These independent readings forecast the crossover and all nine relaxed costs within 0.066 units. Parking the degree service on the same formed stock removes the bulk bias. Work in that comparison is exactly 32 at every separation. The two signed walls still exist. [Microscopic derivation and operating protocol](PROOF.md) · [Registered material preparation](commission.json)

The free flight starts with an eight-site upper interval in a lower background or the opposite preparation, then turns on all fourteen swaps. At the registered exposure π/20, the calculated probability that both original signed walls survive with an additional interior pair is **0.00916174** in the upper launch and **0.0000114445** in the reversed launch. All word probabilities and all terminal word counts are retained. The seeded receiver banks give the following joint-event counts out of 16,777,216 per preparation.

| Board | Upper | Reversed | Degree off, upper | Degree off, reversed |
|---:|---:|---:|---:|---:|
| 0 | 153617 | 209 | 33160 | 33358 |
| 1 | 153625 | 218 | 33261 | 33061 |
| 2 | 153379 | 190 | 32999 | 33272 |

The positive native excess is about 0.00914 on each board, above the predicted 0.001 threshold. A simultaneous 99% Hoeffding radius of 0.000669 for the nine primary differences quantifies the seeded sampler's counting resolution. Suspending the transverse swaps prevents every new pair at every enrolled time. At time 16 the upper preparation gives 461, 477 and 513 counted four-site opposite-phase islands with both original walls present. The reversed preparation gives zero at this sampler depth. The exposure dependence is nonmonotone, including a small reversal of the extra-pair contrast at time 1/2. Measurements on fresh copies establish nucleation during reversible flight. They do not track an individual pair through its later motion. The deposited Born laws come from native-sector evolution. The quadrature totals are seeded binomial draws, and the complete word totals are seeded multinomial draws. All boards realize the same effective Hamiltonian on independently compiled physical addresses. [Source and terminal-word scorer](material.py) · [Counts-only phase decoder and word tables](measure.py)

Along the active material one dimer step is two physical seam edges. The inferred short-distance slope is about 3.998 exchange-energy units per active seam edge. Across the full board, the two wall contacts can have the same graph separation at different material lengths. On level zero, for example, lengths four and five each admit contacts six seam edges apart, with relaxed energies 63.991657 and 71.987477. All 234 embeddings across nine lengths and three boards are included in the independent ruler audit. [Ruler enrollment](RULER_PROTOCOL.md) · [All midpoint distances](ruler.json)

**OW-22 remains partial.** This is native finite-distance tension, screening in the static ground, and new signed walls in the fully released material. When all interior states are allowed, the exact all-length ceiling for this source family is 80. The observed charge is a change of material orientation. Its identification with colour, the Coulomb contribution of a Cornell potential, and independent standards relating native energy and length to about 1 GeV/fm remain open. [Variational proof](PROOF.md) · [Machine-checked bound](MobileString.lean)

Run the following from the repository root to replay the material checks and receiver inference. Use the measured native stock and independent clock in the deposited bank. The optional acquisition command reproduces the archived samples from seed 2026112622.

```sh
pkg=model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string
folder=model/candidate_laws/interval_echo/openwave/row_22_confinement/mobile_string
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 15 .venv/bin/python -B -m pytest -q "$folder/test_mobile_string.py"
(cd Lean && nice -n 15 lake env lean "../$folder/MobileString.lean")
nice -n 15 .venv/bin/python -B -m "$pkg.measure" --acquire
```
