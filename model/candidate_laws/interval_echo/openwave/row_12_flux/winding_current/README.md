# What passes the observer when a wall crosses a cut?

On the ten-word cycle of ideal period-three walls, a gap travels around the belt three times while an ordinary material-current meter counts exactly −1 or +2 units. In free Interval Echo flight, additional routes through other words carry current. For the four-load wall these routes can return the original word after an odd number of crossings. The same departure from a rigid two-unit pump is visible in the energy spectrum of the isolated wall band. **OW-12 remains Partial.** The local instrument measures conserved port load. Identifying the carrier as a quark of fractional electric charge also requires a separately established electromagnetic coupling and electric standard. [Official OW-12 criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L114) · [Global observable ow-12-quarks](../../../../../../docs/registers/OBSERVABLES.json) · [Earlier flux-dispersion experiment](../README.md)

## A calibrated current and a freely returning wall

The preparation uses the twenty contacts on the ten-port belt of one icosahedral cell, blank seam memory, contrast eight and the degree compensation of the original belt experiment. A local two-position memory flips when an ordinary load crosses the chosen cut through any of the three traversing contacts. Its coupling is the conjugation of each such exchange by a controlled flip on one endpoint. The source, pointer coupling, synchronized choice of contacts, final occupation read and time unit are apparatus inputs. There is no intermediate wall measurement in a free-flight launch. The memory records crossing parity and the final port word on the same launch. Calibration with one ordinary load in both orientations, an empty pair, parked/half/full transfers and all twenty contacts at all ten cut positions gives 1,800 joint calibration histograms. A parked ordinary load gives zero flips. Every measured crossing gives one flip. [Protocol](PROTOCOL.md) · [Calibration counts](calibration.json)

The sharp walls form a connected cycle. Each directed exchange moves the gap centre by +3 positions and one microscopic load by −1 (three-load gap) or +2 (four-load gap). On a complete ten-exchange tour, each fixed cut registers three wall passages and signed material current −1 or +2. Removing reversed edge pairs from a closed walk leaves an integer multiple of a full tour. Thus the even parity prediction for any four-load return holds exactly on the sharp graph. At cut zero, only two of the three positive wall passages in a forward four-load tour coincide with actual load crossings. A single ideal gap passage already distinguishes the fractional *cycle average* from the response of an individual transferred ordinary unit. [Complete cut-zero tours and word decoding](results.json)

| Original word | Shortest same-word odd path in full law | Count of shortest paths | Odd returns at time 16, cut zero, per 1,073,741,824 launches |
| --- | ---: | ---: | ---: |
| Three ordinary loads at ports 0, 3, 6 | 5 hops | 18 | 40,096,218 |
| Four ordinary loads at ports 0, 1, 4, 7 | 5 hops | 32 | 630,665 |

For four loads at time 16 the ten detector positions gave 630,188–633,061 odd returns each. The simultaneous 99% exact-binomial intervals over the primary menu give a lower bound of 0.00058395 at this time on every cut. The completely disconnected meter, parked exchanges and restricted sharp cycle give zero. The effect resolves on every cut at all seven nonzero registered waits. Each shortest odd-return path has the same leading transition amplitude, since every exchange is −1/2. The coefficient of the time-to-the-fifth term is the number of paths times (i/2)⁵/5!. Additional diagonal energies only enter at higher order. This proves that the signal exists in the complete finite Hamiltonian independently of synthetic shot fluctuations. The three-load ideal cycle first allows odd return at ten hops. All clock times, changed-law outcomes, terminal words and detector parities are retained. Each bank in the archive has axes [law/measurement, cut, wait, word, pointer]. [All 800 bank histograms](counts.npz) · [Graph witnesses, bank indexing and bounds](results.json)

## Testing the isolated band

Could the odd event arise wholly from the small overlap of a sharp initial word with other bands? Holding the same counter in either of its two X eigenstates gives H(0) and H(π). The latter multiplies the exchange across the physical cut by −1. On the even-pump graph every closed path encounters an even number of such signs. A diagonal change of signs at word vertices then relates the two Hamiltonians, so their energy multisets must agree even when the states within a degenerate energy level change. In the four-load native top band the two spectra differ. The line at −9.15990442 lies 0.00889684 exchange units from the closest opposite line. The conservative certified lower bound is **0.00856873**, while both gaps below the top ten states exceed 13.96. Exact even-pump spectral equivalence therefore fails within the isolated band itself. This still permits an approximate fractional effective charge or a fractional part with additional integer transport. [Spectral test protocol](BAND_PROTOCOL.md) · [All spectra, enclosures and clock windows](band.json)

| Loads | Law or domain | Certified minimum separation for nominated line | Original clock launches in window, X+ / X− |
| ---: | --- | ---: | ---: |
| 3 | Registered contrast | 0.26996 | 208,625 / 7 |
| 3 | Sharp cycle | 0.27877 | 209,279 / 6 |
| 3 | Exchange only | 0.41858 | 68,866 / 348 |
| 4 | Registered contrast | 0.00857 | 206,994 / 6 |
| 4 | Sharp cycle | 0 | 101,800 / 102,198 |
| 4 | Exchange only | 0.19518 | 148,935 / 4 |

Each clock column uses 1,048,576 original inputs per preparation, with every off-window outcome retained. The zero-contrast controls have no positive gap underneath the top ten in at least one counter branch. They demonstrate sensitivity to the contrast-dependent isolation criterion. The choice of spectral window uses the calculated top ten energies and a source-overlap threshold of 0.005. Its counts are a separate numerical phase-estimation experiment with that nominated window, and do not constitute blind discovery. The clock uses 2²³ coherent times spaced 1/128 exchange units apart. Such a clock is a supplied readout service; its implementation inside the belt has not been demonstrated. The full histograms use sparse [signed bin, count] arrays with zeros implicit. [Clock outcomes](clock.npz)

The exact-integer enclosure uses the dyadic approximants V and D to the numerical eigenvectors and eigenvalues. If the Gram error has norm at most g < 1, the polar factor of V is unitary and differs from V in operator norm by at most g. The ordered energy error is consequently bounded by the residual norm plus (‖H‖+‖D‖)g. Integer row and column sums bound the required norms. The recorded bound additionally allows rounding of the nominated numerical eigenvalues to the certified dyadic grid. The independently simulated clock uses rounded numerical lines and source weights. The five-significant-figure spectral separation comes from the integer enclosure. [Certificate and clock sampler](band.py)

## Reproduce

From the repository root, with its Python dependencies installed, use a single numerical-library thread for the seeded receipts.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.candidate_laws.interval_echo.openwave.row_12_flux.winding_current
python -m "$pkg.calibrate"
python -m "$pkg.measure"
python -m "$pkg.band"
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_12_flux/winding_current
```

The tests independently assemble the doubled detector Hamiltonian contact by contact, propagate it in the original pointer basis, check two physical cut placements against the Fourier calculation, evaluate the first nonzero odd-return coefficient by integer matrix powers, and verify normalization of every archived histogram. The signs of signed current and wall advance use the displayed cyclic port order. All selections of the antipodal pair give an isomorphic deformed ring in this one-cell geometry. The calculation addresses this belt and preparation family.
