# Which charged formation energies can passive seam marks produce?

A single vacancy on a filled Interval Echo board can have a much smaller positive formation energy when some seams hold native zero-contrast marks. With any choice of these marks, the ratio between **any two positive one-vacancy formation energies** is at most **122.289** above the corresponding marked, filled host. This excludes both the electron–muon and electron–tau mass ratios on this entire family of preparations, including excited levels. The corresponding single-excess ceiling is $69/8=8.625$. The energy convention fixes the clean private contacts and assigns the filled host zero. The result applies across the deposited federations and refinements with the same five/six-degree cell architecture. [Protocol](PROTOCOL.md) · [Operator and observer records](results.json)

## Which contact does the mark turn off?

Prepare one material load on each port. On any selected seam apply the native echo to the filled blank pair. Its output has unchanged loads and lies entirely in the two zero-contrast coin labels. Both parts of the free generator annihilate the eighteen local columns with those labels, for all nine choices of endpoint loads; the remaining blank columns carry the registered exchange and contrast. The experiment checks the native 99-state matrices using integer products of the write pattern. Thus every selected seam stays silent during free release, including when the charge reaches it. Every unselected seam stays live. The mask can include external seams. In a coherent superposition or mixture of these mask sectors the Hamiltonian is block diagonal, and the bounds apply to their combined formation spectrum. A source with nonzero written contrast, extra mobile pairs or altered private references has a different generator. [Native computation](measure.py)

Let $d$ count live incident seams and let $p=6-d_{\mathrm{physical}}$ be one on a five-seam port and zero on a gateway. With the same filled and marked host subtracted, the vacancy and excess generators on the graph of live edges are
$$
 H_- = L/2+8\,\mathrm{diag}(d-p),\qquad H_+ = L/2+8\,\mathrm{diag}(d+3p).
$$
The exchange form is one half the sum of squared edge differences. Its upper bound $L/2\leq\mathrm{diag}(d)$ gives $H_-\leq54I$ and $H_+\leq69I$. Nonisolated excess states cost at least eight. An isolated excess has energy zero at a gateway and 24 at a five-seam port. An isolated vacancy has energy zero or $-8$ at those same port types.

## Sharp positive vacancy floor

Put $a=(19-\sqrt{297})/4=0.4415780151$ and $b=a/(1-2a)$. On each edge from a leaf with amplitude $u$ to a non-leaf with amplitude $v$, completing the square gives
$$
\tfrac12|u-v|^2\;\geq\;a|u|^2-b|v|^2.
$$
The leaf potential is nonnegative. A core vertex with $d\geq2$ has potential at least $8(d-1)\geq a+db$, since $8-2b=a$ and $b<8$. It pays for every neighboring leaf and its own $a$ contribution. All other edge forms are positive. This proves $H_-\geq aI$ on every connected component of size at least three. A two-port component with two five-seam ends has energies 0 and 1. If a gateway is present, apply the leaf inequality to the other end and use the gateway potential of eight to cover $a+b$. Together with singleton components, the entire vacancy spectrum is contained in $\{-8,0\}\cup[a,54]$. Three live five-seam ports in a path attain $a$. [Lean quadratic-form bounds](../../../../../../Lean/SourceLaw/SilentSeamVacancy.lean) · [Independent small-graph and conjugation tests](test_silent.py)

On the same component, $H_+ - H_- = 32\,\mathrm{diag}(p)$. A connected finite component has a strictly positive ground vector for the exchange generator with any diagonal potential, by the Perron–Frobenius theorem applied after a scalar shift and sign reversal. If the component has a five-seam port, the excess ground state therefore has strictly larger energy than the vacancy ground state. If all ports are gateways, both generators coincide. The all-gateway graphs on the 20-, 80- and 320-cell federations are forests of maximum degree two, with largest components of 6, 6 and 12 ports. Removing additional edges gives paths and isolated ports. On each nontrivial all-gateway path both signs have energies in $[8,18]$. This compares conjugates on identical supports. Coincident energies on different components require their own source assignment.

## What did the released defects actually do?

In the first cell leave just the two seams of the path $5$–$0$–$7$ active on the whole federation. These three ports each have one clean private contact. The vacancy energies above the filled reference are **0.441578, 0.5 and 9.058422**. The excess energies are higher by exactly 32 at all three levels. A sharp vacancy or excess at port 5 propagates freely. The terminal echo read against a zero-load reference resolves which port carries the signed defect. Each of the four banks includes 32,768 starts at each of 1,024 exposure settings and retains a missing-charge outcome. The three-port generator is identical on the three deposited boards, so a single count bank represents each law/sign setting.

| Contact contrast | Defect | Fast frequency, hub occupation | Slow frequency, leaf imbalance |
| --- | --- | ---: | ---: |
| 8 | Vacancy | 8.616865 | 0.058425 |
| 8 | Excess | 8.616855 | 0.058424 |
| 4 | Vacancy | 4.717016 | 0.108491 |
| 4 | Excess | 4.716993 | 0.108498 |

The equal native frequencies for opposite signs follow from the common energy shift on these three ports. Their nominal spectral gaps are 8.616844 and 0.058422. The signed position records resolve those gaps. Absolute formation energy additionally requires a material reference. Varying the contrast at fixed incidence and reader changes both frequencies. [Decoded count tests](test_silent.py) · [Complete compressed counts](results.json)

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.silent_seams.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/silent_seams
```
