# Can the low-energy three-cell state contain three resolved walls?

A ten-port belt supports a phase slip of the period-three density pattern. Depending on whether three or four ports are occupied, this wall has material polarization −1/3 or +2/3. Connecting three belts allows one to compare three resolved walls with three spatially separated clusters of the same conserved loads. In any state spectrally supported at or below the separated preparation's mean energy, the probability of resolving all three walls is at most **0.181% for net polarization zero** and **0.149% for net polarization one**. The result includes coherent superpositions, virtual hops and all free exposures on this guide. [Single-wall and three-wall experiments](../quark_winding/belt_walls/README.md) · [Protocol](PROTOCOL.md) · [Counts and certificates](results.json)

Use three connected cells of the 20-cell federation, ten belt ports per cell, all 60 belt contacts and both intercell bridges. The preparation supplies the guides and their degree references. Each wall is read from the spacings around a belt. The comparison puts three consecutive clusters on separate belts away from the bridges. Both states use ten or eleven ordinary loads. The comparison state supplies a variational upper bound on the minimum energy and can be prepared with the same conserved load counts. Its energy is a finite comparison, without an identification as an asymptotic breakup threshold.

| Total loads | Basis configurations | Three-wall words | Separated mean energy | Three-wall source mean | Three-wall probability ceiling from the mean energy alone | Ceiling with spectral support below the comparison |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 30,045,015 | 3,000 | −167 | 3 | 56.25% | 0.1801% |
| 11 | 54,627,300 | 3,000 | −199 | −12 | 50.37% | 0.1482% |

All energies subtract the same constant within each sector. On each seam an occupied pair and the exchange-odd outcome are orthogonal alternatives. A sharp unequal pair gives either exchange parity with probability one half. With 16,384 independent simulated preparations per seam, the separated energy estimates were −166.995 and −199.003; the three-wall estimates were 2.997 and −11.997. Reversing the contrast while preserving all routing reverses the ordering. [Source words and seam counts](results.json)

## Why a bound on the whole spectrum is possible

Write the shifted Hamiltonian as $H=K-16A$, where the exchange operator $K$ is a sum of positive seam projectors and $A$ counts occupied pairs on seams. Let $P$ project onto all sharp three-wall words and $Q=I-P$. Within $P$ at most three or four seams respectively have two occupied endpoints. Grouping all configurations by each belt's load count and its bridge bits gives full-sector maxima of 20 and 21. Positivity of $K$ yields the mean-energy ceilings in the table. They alone allow a majority of the probability on the wall words.

Write $F$ for the spectral projector of $H$ up to the separated comparison energy $E$. On the wall space $PHP$ is bounded below by $L=-48$ or $-64$. Each matrix entry connecting a wall word to its complement is $-1/2$. Exhausting the 3,000 words gives at most 34 or 36 exits per word, and at most three wall parents for any outside word. The row and column sums give $\|PHQ\|\le\sqrt{34\cdot3}/2$ or $\sqrt{36\cdot3}/2$. On the respective subspaces the operator $X=PF$ satisfies
\[
  (PHP)X-X(H|_{F})=-PHQ\,QF .
\]
The integral solution with factors $e^{-tPHP}$ and $e^{t(H|_F)}$ bounds $\|X\|$ by $\|PHQ\|/(L-E)$. Squaring gives the last column for every state in $F$. Shifts by a common sector energy cancel. Upper-band composites, dressed walls giving other records, and alternative guides require separate tests. The three-wall read measures material polarization; this comparison performs no colour or electric-charge measurement.

From the repository root run the measurement and the independent checks.

```sh
python -m model.candidate_laws.interval_echo.openwave.row_13_walls.energy
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_13_walls
```
