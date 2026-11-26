# A neutral material profile and its charged comparison in Interval Echo

**OW-14 remains open.** A prepared neutral four-defect state on a unit-filled cell has a reproducible positive core and negative distant shell around each of three preassigned gateways. Centering the detector on the defects themselves reverses both signs. Only 0.233% of its terminal records put all four defects in one connected cluster. The fixed-gateway pattern identifies a polarization caused by the cell's native site references. The measured quantity is signed material load relative to the filled host. Assigning a neutron or a proton still requires an independently calibrated electric detector, free particle transport and a rest-energy measurement on the same preparations. [Protocol](PROTOCOL.md) · [Complete simulated counts](records.json)

## Starting from the OW-13 trimer

The earlier freely released source puts three unit loads on an empty host. Swaps preserve that multiset. Each port's excess is nonnegative, so any nonnegative spatial average of its excess remains nonnegative. There is no negative shell or neutral total for this observable. On a uniform unit-filled host, exactly three nonhost ports carry three values in {−1,+1}, whose sum is odd. The earlier three-constituent colour witness does not supply fractional electric charges. An earlier signed written-junction comparison used a separate extension. [Released trimer](../row_13_baryon/open_release/README.md) · [Written-junction comparison](../quark_winding/written_junction/README.md)

This test keeps the registered exchange and contrast (κ = 8) and records the blank-memory sector. Cell 0 supplies the source and cell 3 supplies storage. Both have the same three degree-six gateways at labelled sites 1, 3 and 10. Their exterior contacts are shuttered. Every one of their 30 internal seams and nine private reference terms remains active. In each occupied cell the apparatus starts from an incoherent uniform distribution over the supplied multiset. Its coherent timing filter is calibrated to the sector floor on a 1/128 grid. Two successive central outputs select a narrow energy band. The timing apparatus, prepared inventory, shutter and spectral calibration are supplied. The program simulates their projective outcome probabilities, including rejections, through exact-sector diagonalization and seeded aggregate sampling. [Preparation and measurement](measure.py)

## What the observer reads

On the neutral preparation with two excesses and two vacancies, a terminal snapshot gives the signed port loads relative to the host. Its intrinsic center minimizes the sum of squared graph distances to the four occupied defect positions, using no signs. Ties are drawn uniformly. The independent landmark read uses each predeclared gateway as a center on a separate cohort.

| Center | At center | One link | Two links | Three links |
| --- | ---: | ---: | ---: | ---: |
| Intrinsic center | −0.001623 | −0.002420 | +0.004043 | 0 |
| Gateway 1 | +0.332564 | −0.234778 | +0.119671 | −0.217457 |
| Gateway 3 | +0.332734 | −0.234373 | +0.119131 | −0.217492 |
| Gateway 10 | +0.332415 | −0.235184 | +0.120048 | −0.217279 |

Each row reports mean **total** shell load. Around the intrinsic center, the simultaneous 99% count intervals are [−0.001759, −0.001487] at the center, [+0.001789, +0.006299] at distances at least two, and [0.2263%, 0.2404%] for a connected four-site record. Around every gateway the distant total is negative. The degree references pick the landmark pattern even though the preparation supplied no signed site labels. The filtered mixed state is a function of the guarded Hamiltonian, so the same probabilities apply after every guarded free delay. A translating excitation on the open board would require a new preparation and readout. [Record-only inference](readout.py)

Under sign conjugation, each deviation changes as $q_i=n_i-1\mapsto -q_i$. The exchange and $(q_i-q_j)^2$ terms are invariant. The private reference contributes $8\sum_i b_i(q_i^2+2q_i)$, where $b_i=6-d_i$. Only $16\sum_i b_i q_i$ changes sign. Shutting the private references restores sign conjugation on the entire neutral sector. The uniform seed, the timing filter and center selection then respect this symmetry, so the exact signed mean on every shell vanishes. The sampled shut-reference panel, exchange-only panel and unfiltered source all agree with zero at both kinds of center.

The balanced scattering probe on this same source resolves the underlying populations. On connected four-defect records, its intrinsic-center elastic probability is 0.51571, and calibrated fringe and intensity banks give roughly 26% excesses and 26% vacancies there. The whole neutral cell gives a dark forward port until the receiver records which site it probes. [Contact scattering and calibration](scattering/README.md)

Opening every bridge and internal contact sends positive signed material to the three neighboring gateway ports in a short free flight. Across the three federation sizes, the launch-centered signed second moment is positive at each of four clock appointments. [Unrestricted release and local terminal reads](free_release/README.md)

## An exact bound on connected preparations

The neutral sector contains $\binom{12}{4}\binom{4}{2}=2970$ words. Exactly 1260 have connected defect support. Compress the native Hamiltonian to those words. Its off-diagonal entries are nonpositive. Thus each Gershgorin lower edge is the corresponding row sum. Enumeration in doubled integer units gives a minimum of 310. Any state confined entirely to connected words has energy at least 155 in the original units. The legal basis word with excesses at 0, 1 and vacancies at 2, 4 has diagonal energy 144. The variational principle therefore separates every entirely connected state from the unrestricted ground by at least **11 exchange units**. This statement concerns support at one appointment and allows arbitrary coherent superpositions over all connected words. An independent sparse diagonalization gives floors 142.242600 on all words and 158.264213 on connected words. [Enumeration and variational certificate](compact.py)

## Comparing energies on one reference

The neutral data cell holds two excesses and two vacancies. The positive one holds two excesses and one vacancy. Independent cohorts sample the joint swap-parity and squared-contrast outcome on a random seam, and the relative squared load at a random private-reference port. The native formation energies relative to the *same* filled host are 142.243 ± 1.030 for neutral and 119.184 ± 0.601 for positive. The positive run also places one vacancy in the matching storage cell, where it costs 32.600 ± 0.099. Storage is filled in the neutral run. All energies use the exchange-clock unit.

The neutral data cell costs **23.059 ± 1.630** more locally. The *two-cell* endpoint difference, neutral plus filled storage minus positive plus stored vacancy, is **−9.541 ± 1.729**. Both endpoints contain the same total inventory. A uniform energy per excess $\mu\sum_i q_i$ shifts the local difference by $-\mu$ and the vacancy cost by $-\mu$, and cancels from the two-cell result. It serves only to diagnose the reference dependence. These measurements compare separately prepared endpoints. A neutron–proton rest-mass comparison would require free particle preparations on a calibrated common energy scale.

The Lean kernel checks the 30 labeled seams and degree five at every port, the 2,970-word census, connectivity on all 1,260 compact words, the doubled row bound of 310 and the doubled trial energy of 288. The connected matrix row sum bounds its minimum eigenvalue by Gershgorin’s theorem. The Rayleigh quotient of the unrestricted trial bounds the full minimum from above. The numerical eigensolver supplies neither bound. From the repository’s Lean directory, run `lake build NeutralClusterEnergy`. [Kernel-checked finite certificate](../../../../../Lean/Screen/NeutralClusterEnergy.lean)

## Reproduction

Run from the RER root with its NumPy/SciPy environment. The default full replay compares every archived count bin, preparation rejection and numerical panel. Each cohort starts with $2^{34}$ simulated preparations. The intrinsic neutral native cohort accepts 10,101,424 of them. Energy intervals use bounded-variable Hoeffding radii. Signed shell intervals sum forward-checked Clopper–Pearson bounds on the integer tail counts. With fewer than 2,500 interval components across the panels the union bound gives at least 99% simultaneous coverage conditional on the accepted sample sizes. The sparse enumeration and spectrum comparison are independent of these count intervals. [Verification tests](test_profile.py)

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -m model.candidate_laws.interval_echo.openwave.row_14_profile.measure
python -B -m model.candidate_laws.interval_echo.openwave.row_14_profile.readout
python -B -m model.candidate_laws.interval_echo.openwave.row_14_profile.compact
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_14_profile
```
