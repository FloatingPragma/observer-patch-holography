# Does a locally written medium carry a Klein–Gordon mode?

A uniform local commissioning service opens energy gaps on all six tested boards. The selected band can carry a sharp load several cells under free evolution. The proposed spatial waves each have multiple frequencies in that band, even when the source is confined to ports with exactly one write. These preparations leave OW‑30 open. [Protocol](PROTOCOL.md) · [Native experiment](experiment.py) · [Complete construction, clock and receiver records](results.json)

## What forms and what propagates

Start with empty ports and blank collars. Each seam receives two independent fair bits from a native half-contact on a fresh dimer followed by a signed port read. Two positive reads trigger the native echo on that empty seam. The echo writes a zero-contrast collar. The experiment specifies this commissioning rule uniformly, without a path, repeated geometry, or target gap. After commissioning, the unmodified registered generator acts on every seam and degree reference. The formation rule, source and measuring apparatus are supplied. This tests material sustained in free evolution after commissioning. The free law on a completely blank initial state preserves that sector. [Registered law](../../SPEC.md) · [All-blank control](../row_30_route/BLANK.md) · [Earlier route-dependent gap](../row_30_route/README.md)

In the one-load sector the empty-written contact annihilates the endpoint vacuum and both endpoint load states under each local generator term. Blank contacts exchange the load. Write \(r_v\) for the number of empty writes incident on port \(v\), and \(L\) for the Laplacian of the contacts still conducting. The blank private degree references restore missing degree up to six. In native exchange units the exact one-load Hamiltonian is
\[
 H=8\,\mathrm{diag}(6-r_v)+L/2.
\]
Every board has maximum physical degree six. Thus \(0\leq L/2\leq6I\). By the min–max principle, adding it to the diagonal moves each ordered eigenvalue upward by at most six. The closed intervals \([48-8r,54-8r]\) are disjoint. Each contains exactly as many eigenvalues as there are ports with \(r\) writes. This counts all components. The experiment checks the complete spectrum of its largest conducting component on every board, and an independent test diagonalizes the full 20-cell boards. The separate exchange-only control retains precisely the same deposited contacts and gives \(H=L/2\). All five enrolled spatial waves then have a single frequency, with free-generator residual below \(6\times10^{-15}\).

One write is the most frequent count on each largest component. A sharp launch from its first one-write port has over 99% weight in the predicted interval [40,46] on all six boards. A 256-port timing interferometer tags frequencies [39,47] against a zero-energy stored load on an isolated six-seam port. The timing arms, Fourier recombination and native spatial source pulses are controlled measurement services. The archive includes all tag failures. After the tag, the load propagates freely for the indicated additional exposure.

| Commissioning seed | Cells | Cells on largest component | At least two cell hops, exposure 0 | Exposure 8 | Exposure 64 | Exposure 512 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 178121 | 20 | 19 | 155 | 170 | 249 | 333 |
| 178121 | 80 | 78 | 157,217 | 223,905 | 204,042 | 199,212 |
| 178121 | 320 | 297 | 1,508 | 4,038 | 2,262 | 4,406 |
| 178122 | 20 | 18 | 0 | 0 | 0 | 0 |
| 178122 | 80 | 65 | 1 | 2 | 0 | 0 |
| 178122 | 320 | 307 | 0 | 0 | 16 | 291 |

Each entry counts marked arrivals among 1,048,576 launches. Distance is the shortest distance in the original cell graph. The finite timing instrument allows travel before exposure zero. On the first 80-cell board the remote fraction among tagged launches rises from 15.16% to 21.58% after eight additional units. The independent 80-cell board gives two such detections at the same exposure. The two layouts have very different transport from the prescribed launch port.

## Can one wavelength identify one frequency?

Use the commissioning records to construct the spatial Laplacian of transmitting seams. Prepare its constant wave and four lowest nonconstant waves on the largest component, using local splitters and degree-reference phases. The independent spatial read samples transmitting bonds with a phase and half-contact, giving the antisymmetric probability \(\langle L\rangle/(2|E|)\). Frequencies come from fresh copies and an independent clock. The delay step \(2\pi/64\) covers the entire exact spectrum [0,54] without aliasing. Each single frequency in [40,46] gives over 80% probability in its two nearest timing bins. The exact Fourier sidelobes allow at most 4.76% leakage into [39,47] for an eigenline in any other coordination band.

The strongest pair of adjacent timing bins receives 34.23–50.07% of tagged returns for the six constant waves. Within the exact band the two 320-cell constant waves have standard deviations 0.51222 and 0.52667 in native frequency units. All thirty component-wide spatial preparations exclude a single frequency in that band at simultaneous 99% confidence under the registered sampling model. The smallest excess over the rejection threshold is 0.0404 of all launches.

The follow-up prepares the uniform wave on the largest connected set of one-write sites on each board. Its internal bonds read zero spatial variation. The timing measurement retains over 97% of this source, so mixing widely separated coordination bands makes only a small contribution. On fresh copies we undo the native source after the time read and count recovery of the exact incoming spatial wave by the launch-port detector.

| Seed | Cells | Ports in one-write domain | Tagged launches | Incoming wave recovered given tag |
| --- | ---: | ---: | ---: | ---: |
| 178121 | 20 | 14 | 98.01% | 44.86% |
| 178121 | 80 | 16 | 97.67% | 58.32% |
| 178121 | 320 | 18 | 97.37% | 43.34% |
| 178122 | 20 | 15 | 97.76% | 53.93% |
| 178122 | 80 | 23 | 97.99% | 25.16% |
| 178122 | 320 | 25 | 98.15% | 38.90% |

For the joint test, let \(p\) be the tag probability, \(q\) the joint probability of a tag and recovery, and \(\epsilon\) the maximum foreign-band leakage. Suppose all in-band spectral weight \(w\) of the source lies at one frequency. Write \(A\) for its vector of accepted complex timing amplitudes. Its contribution to \(p\) is \(w\|A\|^2\ge p-\epsilon\), so \(w\ge p-\epsilon\). In the spatial-return amplitude the same contribution has norm \(w\|A\|\ge p-\epsilon\), since \(\|A\|\le1\). The total foreign amplitude in those bins has norm at most \((1-w)\sqrt\epsilon\). The reverse triangle inequality then gives

\[
 q \;\ge\; \big[\max(0,\,p-\epsilon-(1-p+\epsilon)\sqrt\epsilon)\big]^2.
\]

The lower bound uses the measured lower confidence limit on \(p\). The upper confidence limit on the observed \(q\) falls below it on all six one-write domains, by 0.257–0.587 of all launches. The bounds use Hoeffding's inequality simultaneously across both independent banks, all three board sizes, six spatial preparations per board, and the clock and return tests. They quantify simulated counting fluctuations. A narrower source or a spatial mode that includes amplitudes on neighboring coordination classes could still select a single branch. A candidate for OW‑30 would need to select that mode independently of its frequency and predict the held energy–momentum readings under free propagation.

## Replay

From the RER root, use the model's NumPy and SciPy environment.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_formed.experiment
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_formed/test_formation.py
```

Replay regenerates the native repair words, spatial detector reads, clock bins and radial counts. Use `--output NEW.json` to save a separate receipt. [Independent checks](test_formation.py)
