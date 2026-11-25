# Counting the boundary of a material product-address ball

Three doubles moving on independently surveyed 12-port guides give a triple of position records. A parallel arrival clock measures the distance from a central source chamber; an independent single-coordinate survey counts how many ways leave each ball. The five boundaries have **24, 96, 216, 384 and 600** outward edges. The ordinary unbiased walk on these same guides conserves total flux outside the source. It carries nearly three times as much current per edge through the middle of a face of the second ball as it carries through a corner of that face. [Protocol](PROTOCOL.md) · [Native computation](measure.py) · [Complete acquisition](records.json)

## Source, boards and rulers

On the 20-owner federation the three paths occupy all twelve ports of owners 0, 1 and 2 respectively. The addressed seams are the first lexicographic Hamilton path of a native cell; every consecutive port pair shares a native seam. On owner zero its port sequence is (0,1,5,4,2,3,9,8,6,7,10,11). The other two paths use the identical local ports. The unaddressed contacts are shuttered. The graph of address changes has three independent coordinates because each register is separately prepared, actuated and read. These twelve-site guides are linear in their ordered address record. They bend in the icosahedral cell.

A fourth double on an independent seam of owner three distinguishes a ready chamber from the released states. At ready the guide doubles are parked at indices (5,5,5). One source appointment moves the flag to its fired position and independently chooses index 5 or 6 in each of the three guides. The eight tags are publicly recorded. Each is a different reachable active vertex; the tag choosing (5,5,5) also fires the flag. Reset after escape belongs to external preparation. The admitted graph has exactly eight ready-to-active edges, all internal to every positive-index ball. This construction keeps the nine source and nearest addresses distinct at fixed total material load.

First survey every offered local seam in both directions with a sharp double and a vacancy. A full native π pulse exchanges their addresses, including their complex phases, by the blank 99-state generator. With exchange suppressed at the same exposure the double stays put. The sharp position reader supplies the public outcome. All 33 rail links gave arrival on all 32 enrolled launches per direction; parked exposures gave zero. The source-flag seam has a separate transfer check. The receipt groups deterministic sharp histograms after checking the operator columns. A new clock survey starts at ready, tries all eight initial corners, advances each guide one surveyed step per subsequent appointment, and parks each marker once its detector fires. For all 1,728 target triples it retains all eight arrival times. The measured radius is the minimum of the eight complete times. At most one link per coordinate can be crossed in one range appointment; independent simultaneous steps attain that bound. The resulting distance to an active address $j$ is
\[
n(j)=1+\max_{i=1,2,3}\min(|j_i-5|,|j_i-6|).
\]
A separate breadth-first test enumerates all permitted parallel moves from each initial corner and reproduces the entire timing panel. Parallel ranging and the single-coordinate transport menu use different appointments. Their operations and the central dispatch are specified apparatus. [Range and aperture checks](test_product.py)

## The exact finite boundary count

With coordinates translated by $-5.5$, the active ball of radius $n$ is the cube $|x_i|\le n-1/2$ on the half-integer grid. Each coordinate therefore has $2n$ possible inside positions. An outward *single-coordinate* read changes one axis at just one face of the cube. Choose the axis, the outward sign, and the two tangential positions on that face. These choices determine both ends of exactly one native single-exchange edge; every outward edge determines these choices uniquely. At an edge or corner of the cube, exiting along different axes gives different measured destinations. The count is $3\cdot2\cdot(2n)^2=24n^2$ for every $1\le n\le5$ on the board. Layer six contains the outer detector addresses. The same argument applies to longer paths through any depth leaving room for an outer step. [Lean cardinality of the face records](ProductShell.lean) · [Generic finite-shell theorem](../../../../Lean/Screen/LayeredDiscreteGauss.lean)

Enumerating the 4,752 undirected active product edges and the eight distinct chamber releases by native transfers gives the five counts, independently of the formula. Keeping the original timing calibration but parking the 6-to-7 link on the first axis eliminates 144 successful address-edge transfers and reduces the first shell count from 24 to 20. Parking the 5-to-6 link instead removes another set of 144 address edges without altering any of the five cut counts. The cut census can therefore detect a damaged boundary and distinguish it from a loss of internal links.

The same source permits independent changes of energy and three-guide entropy, and a pointwise test of the far current from measured step rates. [Corner preparations](corner_source/README.md) · [Current and its angular dependence](green_response/README.md)

## Does flux spread evenly?

After each release the ordinary walker selects one of the six signed single-coordinate moves with probability 1/8 each and otherwise idles. A collector absorbs the triple on first reaching layer six. Every interior coordinate move is unbiased. Faces, edges and corners receive no compensating interventions. On a finite box let $V$ vanish at the collector and solve
\[
6V(x)-\sum_{y\sim x}V(y)=\tfrac18\mathbf1_{x\in\{5,6\}^3}.
\]
For one launched walker the expected number of visits to $x$ is $8V(x)$, by the visit-balance equation for the lazy walk. Each individual oriented edge therefore carries expected net crossings $V(x)-V(y)$. Give the chamber-to-corner edges flux $1/8$. The source at ready has total outflow one; every other vertex through layer five has zero divergence, including the eight injection corners after counting their incoming release. All five shell totals equal one. Repeating the operation with reinjection converts per-launch fluxes into steady rates with one common renewal-time factor. The per-launch normalization uses directly counted releases.

The Green function is obtained exactly over the rationals after identifying the 35 interior classes under coordinate permutations and sign flips. Substitution into the full unsymmetrized graph verifies all 1,000 interior vertex balances. For a fixed face of the second ball the four central crossings (both tangential offsets zero from the central rail pair) carry **0.0172204164** each. Its four corner crossings (both offsets one) carry **0.0059800403** each. The complete exact fractions are stored in the receipt. Cubic symmetry relates like edges on different faces but permits these two different values.

For each independent launch we counted the net signed crossings on the 24 middle edges minus those on the 24 corner edges of shell two. All 65,536 original launches reached their collector, the longest taking 391 walk appointments. The measured contrast was **0.271103** per release, against **0.269769** predicted. To bound the uncertainty despite correlations within an episode, let $g(x,y)$ be that signed single-edge reward. With zero collector values, the first moment $u(x)$ solves $(6I-A)u=\sum_{y\sim x}g(x,y)$ and the second moment $s(x)$ solves $(6I-A)s=\sum_{y\sim x}[g(x,y)^2+2g(x,y)u(y)]$. The source's corner-average variance is 0.988725. Chebyshev's inequality across the independent complete launches gives a 99% half-width below 0.03885. The flagged source, occupation records and exact shell count thus supply the finite cardinality and source-balance premises on this chosen graph. Equal sharing fails for the ordinary walk on the very same graph.

The measured radius counts parallel arrival appointments on three guided markers. Interpreting it as a freely measured spatial radius for one body, calibrating its physical units, or identifying these material currents with electric or gravitational flux requires additional evidence. The icosahedral screen's subdivision frequencies do not enter either measurement.

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -O -W error -m model.experiments.interval_echo_shell_flux.product_address.measure
.venv/bin/python -B -W error -m pytest -q model/experiments/interval_echo_shell_flux/product_address
```
