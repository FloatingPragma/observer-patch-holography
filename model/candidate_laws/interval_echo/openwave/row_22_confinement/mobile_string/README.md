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

The shared decoder accepts a line only when the observed quadratures support its winding. For each of the 13 holds, two independent binomial counts with N = 2²² draws give a complex quadrature vector. With component error 5/√N the enclosing disk has radius ε = 5√(2/N) = 0.003453. Hoeffding's bound and a union over 26 components give a probability of less than 0.00020 per preparation that a true quadrature lies outside its interval. For observed magnitude r > ε the phase uncertainty is at most arcsin(ε/r). The decoder checks the principal interval of the first hold and each nearest-turn boundary using both adjacent phase uncertainties and the ratio of their hold times. Its rate is conditional on the chosen first-hold interval and the sharp preparation. The separate eigenweights determine whether that rate can identify a ground transition.

Of 108 stored preparations, 99 have resolved windings. All nine unresolved readings are the degree-off screened preparation at separations 8, 10 and 12 on each board. The analysis gives them a null line and retains their attempted winding and affected holds. All native preparations and the enrolled degree-off intact readings pass. Thus every enrolled ground comparison and both held slope residuals on each board retain a reading. The decoder emits null slopes or residuals when a required line is unresolved. The original receiver bank and original numerical output are preserved. Run `python -m model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string.measure --check` for the current analysis. [Receiver bank](receivers.npz) · [Original numerical output](results_original.json) · [Visibility-aware results](results.json)

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
## Equal-degree route control

On each deposited board, drop the first port of the sixteen-dimer path and append the external mate of its last gateway. Pair consecutive ports on this shifted path. In order, the blank dimers alternate between two degree-five ports and two degree-six ports. Each connector joins adjacent dimers on a native seam. The full-board stock compiler leaves one load per dimer and prepares the same coherent connectors and idle spectator collars. Define the upper port on the right of even-indexed dimers and the left of odd-indexed dimers. This retains the alternating connector orientation of the original strand.

The hypothesis under test is that the bulk slope measures the difference in degree-reference energy between the two orientations. The previously calculated zero-bias prediction is 32 units of static work. This route was selected to test that prediction. Both degree-five ends retain their native private energy on the shifted route. With one load per dimer the private energy is orientation-independent. Inverting every free interior orientation takes the lower-pinned Hamiltonian into the upper-pinned one up to the two boundary walls, whose combined cost is the calibrated contact energy 4J. This predicts constant ground work, including on the held lengths.

Use the original two-copy Ramsey protocol, its 13 holds, both sharp preparations and both quadratures with 2²² simulated draws each. Sample a separate exchange clock and calibrate the contact, swap and private-port energy on each board. Use lengths 2 and 3 as short controls and 4, 5, 6, 7, 8, 10 and 12 as held checks. The intact line must be within 0.002 exchange-clock units of 32 on each held length. A changed-law control removes agreement energy during the witnessed evolution after the same stock preparation. Its prediction is zero work with a unit exchange gap. The routed medium, pins and clock remain supplied. These readings test the material energy of a prepared orientation wall. [Physical ports and acquisition](equal_degree.py) · [Receiver bank](equal_degree.npz)

All 21 held intact lines on the native law lie between 31.9997847 and 31.9999853. The smallest visibility over the holds is 0.9904784. On the original unequal-degree ports the short-work slope is about 8. Re-pairing ports therefore removes the measured bulk tension with the degree service active.

The zero-agreement control exposes a limitation of the phase reader. The interior consists of independent swaps with gap one, so the two-copy signal is real and equals cos(t/2) raised to twice the number of free dimers. Every measured quadrature across both preparations agrees with this zero-work law within 3.68 binomial standard deviations. All 21 held intact windings at κ = 0 are unresolved. The attempted rates differ from zero by more than 0.002 on 12 settings and happen to fall within 0.002 on nine. The same shared decoder accepts all 21 held intact lines at κ = 8. The deposits retain every hold and both ports.

Run `python -m model.candidate_laws.interval_echo.openwave.row_22_confinement.mobile_string.equal_degree` from the repository root to decode the deposit. Pass `--acquire /tmp/equal_degree.npz` to generate a separate bank with seed 220032. The regression checks the physical-port inversion identity and the positive unequal-degree control. [Tests](test_equal_degree.py)
