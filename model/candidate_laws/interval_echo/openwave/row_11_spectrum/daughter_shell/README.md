# Can the charged triplet decay into a lighter charge and a neutral pair?

Observers on the 20-cell federation can see a vacancy cross into another cell while a neutral vacancy–excess pair stays behind. At one eighth of an exchange time the two owners jointly record the chosen three addresses in **11,359 of 4,194,304** simulated releases. The charged parent’s upper band carries at least **11.67857 more energy units** than *any* independent vacancy and neutral pair together can carry away. The energy gap rules out that two-body decay of this band. Other parent bands, additional products, and preloaded material have different energy budgets. [Preparation and protocol](PROTOCOL.md) · [All counts and bounds](results.json) · [Executable](measure.py)

## A charge coincidence and its energy budget

One unit fills each port and all collars start blank. The sharp parent word has vacancies at 1 and 3 and an excess at 8. All 630 native seams evolve with their degree references active. A successful terminal read finds the vacancy at 82 in owner 6 and the neutral pair at 3 and 8 in owner 0. A signed echo against a zero-load reference distinguishes all three port occupancies. Conservation of the individual defect counts guarantees that this read accounts for the entire word. Closing the three bridges of the charged receiver keeps its charge at zero and gives zero coincidences at every exposure. Starting with the output word gives unit acceptance at zero exposure. [Charge reader](../../row_01_charge/measure.py) · [Charged-triplet preparation and band proof](../../row_02_dressed/README.md)

| Autonomous exposure | Native coincidence probability, whole-board enclosure | Simulated coincidences / 4,194,304 | Exchange-only control |
| --- | ---: | ---: | ---: |
| 0 | 0 | 0 | 0 |
| 1/16 | 0.0892924–0.0892930% | 3,588 | 4,037 |
| 1/8 | 0.2743162–0.2748075% | 11,359 | 15,441 |

The two preparations differ only by moving one vacancy across the 1–82 bridge. A native balanced exchange and a phase accumulated on the degree-five port 74 prepare their sum and difference. After free evolution a second balanced bridge exchange compares the two output ports. Terminal words outside the two chosen branches count as a third outcome. From the two signed count differences, observers infer a parent-minus-separated energy of **15.322** at time 1/64 and **15.524** at 1/128. The exact finite-time slopes, before sampling, are 15.3364 and 15.4590. The initial slope is 15.5. Removing the contrast during free evolution gives −0.498 and −0.475 against an initial slope of −0.5 with the same preparation and decoder. These are readings of sharp-state *mean* energy differences. [Native bridge calibration and individual outcomes](results.json)

For the two arms $p,d$ write $X=|p\rangle\langle d|+|d\rangle\langle p|$ and $Y=-i|p\rangle\langle d|+i|d\rangle\langle p|$. Subtracting the two input densities gives $X$. The final balanced contact reads $Y$ from the daughter-minus-parent port record. Thus the half-difference of the readings is $m(t)=\mathrm{Tr}(Y e^{-itH}Xe^{itH})/2$ and $m'(0)=H_{pp}-H_{dd}$. The sharp energies are 199.5 and 184. The latter splits into 51 on the charged port and 133 on the neutral pair. For every charged location in owner 6 and every pair in owner 0 away from its common gateway, the two supports have no connecting edge. In this 1,320-dimensional chamber the compression of the full three-defect Hamiltonian equals the compression of the independent daughter sum entry by entry. The equality describes the energy form on these supports. Free motion can leave the chamber.

## A bound on every independent outgoing energy

For a vacancy at a port of degree $d$, the host-relative diagonal is $8(2d-6)+d/2$ and each neighbor receives amplitude $-1/2$ in the Hamiltonian column. Give degree-six ports weight eight and degree-five ports weight one. The largest absolute weighted off-diagonal row sum plus the diagonal is **209/4** on each federation. Weighted Gershgorin bounds every single-vacancy energy above by this number.

For a neutral pair let $d_h,d_x$ be the degrees of the vacancy and excess and let $a$ be one when their ports touch. Their host-relative potential is $V=8[12+2(d_h-d_x+a)]$. There are $n=d_h+d_x-a$ allowed different-load swaps, each with diagonal contribution $1/2$ and off-diagonal amplitude $-1/2$. The largest potential, 128, requires a touching degree-six vacancy and degree-five excess. Give these words weight seven and every other word weight one. A hop creating a high-potential word must pass through a shared neighbor or interchange the two defects. Each shared neighbor qualifies for at most one of the two possible hops because the final vacancy and excess require different degrees. Enumerating physical two-hop paths gives at most two shared neighbors for any two distinct ports. In a high word the interchange reverses the degrees and leaves the high class. Its weighted row bound is at most $128+5+(2+8/7)/2=942/7$. Every low word has $V\le112$, $n\le12$ and at most three high neighbors. Its bound is at most $112+12+3(7-1)/2=133$. The independent test also checks the weighted row inequality on all 57,360 ordered neutral-pair words of the 20-cell board and the native 99-state seam on each realized local degree/contact class. [Microscopic and full-pair checks](test_shell.py)

| Cells | Upper charged-band rank | Charged parent energy at least | All independent charged + neutral energies at most | Minimum deficit |
| ---: | ---: | ---: | ---: | ---: |
| 20 | 104 | 198.5 | 5231/28 | 327/28 |
| 80 | 384 | 199 | 5231/28 | 341/28 |
| 320 | 1,504 | 199 | 5231/28 | 341/28 |

Each output energy of $H_C\otimes I+I\otimes H_N$ is a sum of daughter energies, including any degeneracies or relative motion on these finite boards. Its spectrum therefore ends at or below $209/4+942/7=5231/28$. On parent and daughter energy eigenvectors the time-translation intertwining equation $H_{\rm out}T=TH_{\rm parent}$ gives $(E_{\rm out}-E_{\rm parent})T_{fi}=0$. All these coefficients vanish for an upper-band parent. The charge convention drops out of the deficit because both sides carry total charge −1 and use the same filled host. A tighter experimental resolution than the listed deficit would separate the energy ranges. Interaction energy, residual material elsewhere, or work at a detector would enter a different channel. The formal proof takes the spectral bounds and the stationarity equation as explicit hypotheses. [Finite energy and transfer lemmas](EnergyShell.lean)

The sharp source has 15 different-load swaps, each leading to a distinct word. Its energy variance is 15/4. The complementary spectrum ends at 192, so the one-sided variance bound puts at most $(15/4)/(15/4+(199.5-192)^2)=1/16$ of the sharp preparation outside the upper band. The flight counts enroll the unfiltered source. They establish a separation of material loads. Assigning energy to an individual coincidence requires a daughter energy read on that trial. A charged-lepton interpretation also requires rest energies, neutral flavors and a decay width. [Prepared upper band](../../row_02_dressed/README.md) · [Two-neutral-product inventory](../neutral_ejecta/README.md)

## Error bounds and replay

The primary finite-time coincidence calculation restricts the 6,825,840-dimensional sector to all 1,721 words within four swaps of the source, retaining their full diagonals. An integer recurrence evaluates the degree-64 Taylor polynomial of the centered Hamiltonian between the two words. Taylor's integral remainder on the imaginary axis is at most $(t\|H-cI\|)^{65}/65!$. In the interaction picture the hopping norm is at most nine, since at most eighteen swaps leave any word. A path omitted by the restriction must exit at distance five and return to distance one, requiring at least nine swaps. The full and restricted amplitudes therefore differ by at most $2(9t)^9/[9!(1-9t/10)]$. Probabilities in the table include both errors and outward rounding. The bridge energy read uses sparse numerical evolution on the radius-three union around both words. Omitted paths there require at least eight hops. A radius-five sparse calculation independently checks the primary intervals. The counts are seeded samples from the computed terminal distributions.

From the repository root run

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.daughter_shell.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/daughter_shell
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_11_spectrum/daughter_shell/EnergyShell.lean)
```
