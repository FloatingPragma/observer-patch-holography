# Which signed inventories form bound aggregates on the same cell?

All five vacancy stocks and all five excess stocks of sizes two through six have a ground energy below **every** proper breakup on separate copies of the same native twelve-port cell. Each copy retains the thirty internal bonds, its nine private references and the same unit-filled host. Binding here is an energy difference in native exchange units. The experiment supplies a coherent energy clock and sharp stock preparations. [Frozen inventory survey](PROTOCOL.md) · [All 42 energy cohorts and sharp census distributions](counts.json) · [Decoder](measure.py)

| Defects on the cell | Vacancy binding interval | Excess binding interval |
| ---: | ---: | ---: |
| 2 | 11.65–21.69 | 0.0764860–0.0764862 |
| 3 | 27.16–39.20 | 2.10–16.53 |
| 4 | 26.71–39.68 | 18.46–32.89 |
| 5 | 26.12–40.14 | 2.00–16.53 |
| 6 | 34.67–48.73 | 25.58–40.99 |

The intervals give the cheapest allowed breakup energy minus the joined ground energy. Except for the two excesses, each interval uses an exact lower bound and a measured upper bound on each fragment or joined ground. Endpoints in the table are rounded outwards. They bound the finite energies; their width includes the choice of trial word. The two-excess interval uses an exact rational certificate checked with a separate clock measurement. No empirical nuclear mass sets any scale.

## How the inventory bounds are measured

On each sharp support an echo against a unit stock determines the signed load. The two write outcomes for the observer's hand agree on the load difference. On a one-sign word the diagonal contrast equals κ times the number of occupied-to-host bonds, plus the private-reference terms relative to the filled host. The census covers every support of each enrolled size. The remaining exchange is a sum of positive odd-swap projectors, so the minimum contrast bounds every energy from below.

To obtain an upper bound, the observer prepares a word minimizing its predicted sharp mean (contrast plus half the number of unequal bonds). Fresh copies receive a uniformly chosen local energy read from the 39 available terms. Multiplying its average by 39 estimates the chosen word's energy, which bounds the ground from above. All 42 cohorts agree with the signed-port and parity prediction within simultaneous 99% counting bounds. For each proper integer partition of the inventory the separate cells have additive ground energies. Taking the minimum of those fragment sums at both endpoints gives the displayed interval. This method certifies nine of ten native binding signs in its first pass. At half contrast it certifies seven. At zero contrast the uniform superposition of each stock inventory has exactly zero energy on every term, including every independent fragment.

## Resolving the two-excess threshold

The original two-excess interval crossed zero. The registered second test bounds the bottom energies using a strictly positive integral vector x on each sector. The off-diagonal exchange entries are nonpositive. Taking absolute amplitudes can only lower energy, so a nonnegative ground vector g exists and has positive overlap with x. The ratio ⟨g,Hx⟩/⟨g,x⟩ equals the ground energy and is a weighted mean of (Hx)ᵢ/xᵢ. The stored integer weights and twice-generator certify every endpoint with rational arithmetic. At native contrast the joined ground costs 84.72475194–84.72475198 units up to rounding, versus 84.80123804 for two independently stocked cells. Halving the contrast gives a binding energy of about 0.240734. The finite spectra do not require a universal increase of binding with contrast. [Refinement protocol](TWO_EXCESSES.md) · [Positive certificate vectors and clock records](two_counts.json) · [Calculation](two.py)

A uniformly sampled sharp stock gives every eigenmode of its inventory equal weight. A controlled Fourier clock with 131,072 equally spaced delays compared the predicted joined and separated lines. Its first zero lies 0.006136 exchange units from a line. At the forecast joined energy 84.72475 it returned 254,096 bright outcomes on the joined pair and 196 on separated single excesses, each in 16,777,216 simulated original starts. At zero contrast the ground line is zero for both stocks. Their bright counts there differ because their zero-energy spectral weights differ.

A blind Fourier scan of all 220 three-vacancy supports resolves at least eight distinct native energies below the cheapest breakup. Its clock bins were fixed before inspecting the excited energies. At half contrast six disjoint windows qualify, and at zero contrast none lie below breakup. [Subthreshold spectrum and leakage bound](LEVELS.md) · [All 65,536 clock bins](levels_counts.json.gz)

The sign of the finite energetic binding changes the interpretation of a material picture. Opposite inventories on a connected blank graph obey the separate sign-fragment bound, whereas identical defects bind on this cell. Identifying protons, neutrons or particular nuclear isotopes and measuring halo extent and lifetime would require other records. [Opposite-sign comparison](../sign_fragments/README.md) · [Trimer clock protocol](LEVELS.md)

At the weak two-excess line the relative position survey found 97.43% of pairs nonadjacent. Independently cooled excesses on matched cells, conditioned on distinct addresses, gave 99.21%. A repeated separation check gave different first-detection rates at two clock cadences, and far pairs could return. This experiment measures a spatial profile and a monitored crossing; a halo lifetime remains open. [Body-relative shell and return records](outer_pair/README.md)

From the repository root use its scientific Python environment with one numerical thread.

~~~sh
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_17_nuclei.like_stock
python -Bm "$p.measure" --acquire /tmp/oph-one-sign.json
python -Bm "$p.two" --acquire /tmp/oph-two-excess.json
python -Bm "$p.measure"
python -Bm "$p.two"
python -Bm pytest -q --pyargs "$p.test_stocks"
~~~

The leading contact-count estimate for vacancy breakup is twice the contrast per shared seam. It omits hopping and the diagonal energy of private references. With nine referenced ports on the tested cell, the sizes two through six have breakup costs about 14.86, 30.59, 32.07, 32.14 and 41.00 in native units. The six-vacancy estimate from contact counting alone is 48. These are finite breakup energies on the reference material.
