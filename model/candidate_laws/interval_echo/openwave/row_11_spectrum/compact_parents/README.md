# A measured charged formation line near the muon ratio

A preparation of 72 vacancies and 71 double occupancies on 143 distinct owners gives **7236.87525 native energy units** against a six-seam material store. The simultaneous 99% interval is **[7236.87446, 7236.87605]**. Its exact formation energy is $7257-9\sqrt5=7236.87538820$. With the 35-unit compact vacancy provisionally identified as the light state, the ratio is 206.76786823. The comparison selects the closest energy after applying an electron anchor to the compact vacancy. The applied phase impulse between neighboring cells leaves the surplus charge in a flat, coherently confined band. Local contributions account for the entire energy, with zero binding between owners. The source is a spectroscopically resolved collection of charged and neutral constituents. [Full counts and parent recipe](results.json) · [Physical circuit](measure.py) · [Recorded target inputs](../../electron_field/lepton_family/arrival_spin/anchor_spectrum.json)

The closest line was selected from 698,858 energies after the target ratio was known. Within ±20 native units of that target the family supplies 1,955 distinct lines, with mean spacing about 0.0205. The reported agreement therefore does not distinguish a muon from other additive compact preparations.

## How many exact lines does the native model produce?

Prepare at most one compact defect in each occupied cell of the blank, unit-filled federation. The complete energy support in this family can be enumerated before supplying any particle ratios. All normal and external contacts are open for free evolution. The spectrum includes adjacent occupied owners.

| Cells | Owners with additional compact frequencies | Distinct net-negative energies | Largest native formation energy |
| ---: | ---: | ---: | ---: |
| 20 | 0 | 10 | 953.0000 |
| 80 | 24 | 10,478 | 4027.8328 |
| 320 | 104 | 698,858 | 16317.2755 |

Write $p$ for the number of double occupancies. There are $p+1$ vacancies. A compact state vanishes at its owner's gateways, so the only occupied positions have degree five. Relative to the filled host their native contrast potentials are 32 for a vacancy and 64 for a double. With $L$ the local Laplacian, the kinetic generator is $L/2$. On each actual gateway chart the verifier computes the largest invariant kernel of the escape rows by exact rational elimination. Its spectrum for $L$ is either $6$ with multiplicity two, or $6$ with multiplicity three and $5\pm\sqrt5$ once each. The code constructs ordinary integer eigenvectors and algebraic vectors $u\pm\sqrt5 v$ for both additional lines. It checks each row and the zeros at all three gateways. [Exact local and federation enumeration](spectrum.py) · [Underlying compact-space computation](../../../../../experiments/interval_echo_encounter/confinement.py)

Every inter-owner seam sees the ordinary word $(1,1)$ throughout the compact product, so each such seam annihilates the state relative to the host. The local generator cannot leave a vacancy or a double at a gateway, including after the kinetic terms arriving there are summed. Thus each owner's local eigenvalue adds exactly for arbitrarily placed occupied owners. This statement also holds for coherent superpositions inside their retained compact spaces. Moving an entire occupied orbital to another compatible empty owner costs zero energy. Any coherent superposition of such placements has the same exact eigenvalue. A phase impulse depending only on the signed loads of whole owners changes their relative phases within that degenerate space, so the spatially resolved signed density stays fixed during subsequent free motion. The independent test applies the complete native diagonal and every off-diagonal swap to three nearby defects with exact coefficients in $\mathbb Z[\sqrt5]$. Changing one amplitude breaks the eigenvalue equation. [Three-defect word test](test_compact_parents.py) · [Lean product eigenvalue proof](CompactEnergySum.lean)

If $h$ owners use the two additional branches and $d$ is the number on the upper branch minus the number on the lower, then
$$
  2E=70+204p-h+d\sqrt5.
$$
Here $0\le p\le\lfloor(N-1)/2\rfloor$, $0\le h\le\min(2p+1,N_{\mathrm{special}})$, and $d=-h,-h+2,\ldots,h$. Every permitted pair of coefficients occurs by placing the factors on the enumerated owners. Distinct pairs give distinct energies by the irrationality of $\sqrt5$. Because $h\le104<204$, the rational coefficient also uniquely determines $p,h$. The table counts distinct frequencies, independent of degeneracy due to placements or local modes.

## Comparison with the empirical and register ratios

The targets enter after construction of the full spectrum. When 35 is the light reference, the closest of the 698,858 lines to the measured muon/electron ratio uses $p=71$, $h=40$, $d=-18$. The nearest to the distinct register ratio 206.77031597 uses $h=51$, $d=-13$ and gives 206.77044452. The first table entry below is the state selected for the physical experiment. Changing the light reference changes which many-body product wins. The nearest-line certificate uses interval arithmetic for $\sqrt5$, the decimal targets and the independently enclosed ground state. The reported standard uncertainties are those of the measured ratio of 206.7682827(46). [CODATA, Eq. 187](https://physics.nist.gov/cuu/pdf/RevModPhys.97.025002.pdf#page=34) · [Lowest charged line](../lightest/README.md)

| Assumed light branch, 320 cells | Closest heavy/light ratio to the muon central value | Difference in reported standard uncertainties |
| --- | ---: | ---: |
| Ordinary compact vacancy, 35 | 206.76786823 | −90.10 |
| Lowest negative state, approximately 32.59604 | 206.76805275 | −49.99 |
| Lower compact vacancy, approximately 33.38197 | 206.76787296 | −89.07 |
| Upper compact vacancy, approximately 35.61803 | 206.76834367 | +13.25 |

Even dividing the largest product energy by the lowest net-negative energy gives a ratio of only 500.591. Every line in this family is below the heavy register target of 3477.473 under all four light calibrations. Other compositions with multiple interacting defects in an owner or states with amplitude at the gateways have different energy supports. The large number of independently placeable constituents, the exposure of the separate signed loads, and the flat spatial response distinguish these formation lines from evidence for a single freely propagating muon.

## What the observers measured

For the chosen 143-defect parent, prepare positive and negative loads by native partial exchanges and signed private-contact phases; the inventory is supplied at preparation. There are 196 external seams between occupied owners. At each of the 143 participating owners in turn, use fresh parent preparations and a common distant reference port of degree six. Split that one defect between a local orbital and the store. Full native swaps and their inverse route this operation through arbitrary intervening words; the longest route is 64 seams. Close the six contacts incident to the storage port for the exposure, leaving its energy at the declared zero relative to the filled host. Each occupied neighbor is a gateway node of its orbital. The other 142 factors acquire the same phase in both arms. Restore the local preparation, combine the paths and measure the two quadratures. The third bin includes every other final port. This measures each occupied factor against the same store, which suffices to reconstruct the formation energy because additivity has been established independently. [Storage protocol and phase unwrapping](metrology.py)

| Branch | Number of vacancies | Vacancy frequency | Number of doubles | Double frequency |
| --- | ---: | ---: | ---: | ---: |
| Lower | 15 | 33.38196611 | 14 | 65.38196588 |
| Ordinary | 52 | 34.99999994 | 51 | 67.00000004 |
| Upper | 5 | 35.61803409 | 6 | 67.61803381 |

Each occupied owner contributes 262,144 simulated starts at each quadrature and wait. Pool counts only for distinct owners in the same sign and orbital branch. An independent native two-port dimer supplies the common time unit, with 268,435,456 simulated starts per setting. Nineteen doubling waits run from 1/64 to 4096. The shortest wait fixes the unambiguous phase; each doubling unwraps the next phase from the preceding observation. Simultaneous Hoeffding bounds on the two quadratures of all seven groups give a phase error below $\pi/3$ at every wait, ensuring unique continuation. The conservative interval includes all six line errors and the shared clock error. The older 20-cell baselines with just a single-time 320-cell differential read give the much coarser reconstruction 7233.51. The direct 320-cell scan above uses fresh reads of every participating owner.

For the independent free-flight test, split the surplus vacancy between adjacent owners at ports 96 and 1162 using their actual external seam. With ordinary compact orbitals on both arms, leave all seams open and apply a relative phase of zero, $\pi/2$ or $\pi$ to the second owner's charge. After inverse preparation and recombination the zero setting returns all 262,144 starts to the bright port at every tested wait. The $\pi$ setting returns all to the other port. At one time unit the measured signed gap is −0.000229, with exact value zero. Changing just one arm to the lower or upper internal branch yields measured frequencies −1.62057 and +0.61584 against predictions −1.61803 and +0.61803. Both branches stay in their own homes. A phase applied site by site within a single orbital can break its confinement and requires a separate propagation experiment. [Free-flight panels](results.json) · [Protocol](PROTOCOL.md)

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.compact_parents.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/compact_parents
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_11_spectrum/compact_parents/CompactEnergySum.lean)
```
