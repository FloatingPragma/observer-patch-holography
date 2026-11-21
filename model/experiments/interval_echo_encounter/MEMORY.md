# A scheduled native repair with its working memory retained

At duration one, cell 0 applies Interval Echo’s native 99-state echo U to seam (0, 1). The seam begins blank. Its action on that input is the contrast write

$$
W|a,b\rangle = \frac{1}{\sqrt 2} \sum_{c=0,1} |c,a-b\rangle\,|y_c(a+b)\rangle .
$$

The first arm keeps that coherent working memory and lets the federation evolve for four more units. The matched erase arm applies U† immediately after U and then evolves for the same duration. The bypass arm uses the same separated preparation and free evolution. All 36 sharp, initially nonadjacent preparations receive the service regardless of the loads present on (0, 1). The only observer reading in the retained arm is the final seam/quiet measurement on cell 0. Its 21 internal seams disjoint from ports 0 and 1 provide the primary outcome. The free flight, the local service and the reading use the fixed sphere and the registered operators. The service appointment, waiting durations and reader are supplied. [Protocol](MEMORY_PROTOCOL.md) · [Native operators](../../candidate_laws/contrast_memory_spherical_v1/law.py)

## Observer results

| Measurement, all 36 separated preparations | Retained write | Write followed by inverse |
| --- | ---: | ---: |
| Final pair on an internal seam disjoint from the repair endpoints | 0.12672518 | 0.00348463 |
| Weight in the reference bound band immediately after the service, global diagnostic | 0.05241507 | 0.00291210 |

Retaining the repair increases the disjoint-seam yield by 12.324 percentage points and the immediate bound-band weight by 4.950 percentage points. Erasing the write restores the bypass probabilities on every local seam and on quiet, within 10⁻¹⁰. The disjoint terminal detector and unconditional enrollment exclude a coincidence with the written seam and selection on a first detection. The comparison changes both the coherent load-memory state prepared by the service and the subsequent response of the tagged seam. It measures their combined effect. The schedule gives the load-memory system an external intervention. A first record of autonomous encounter and a later repeated reading of that object are separate requirements. [All terminal probabilities](memory.json)

## Evolution with a written seam

The other seams have blank working memory and preserve their blank subspaces. On the addressed seam the matrix is UGU† throughout the flight. The code removes the blank term for that seam from the two-load Hamiltonian, acts with the remaining terms on every memory component, and inserts the native nonblank term on its two actual endpoints. This accounts for the physical effect of leaving that seam written during the flight. [Implementation](memory.py)

Binary input loads have contrast 0 or ±1. The write occupies the symmetric coin combination for zero contrast and two coin states for each parity of the ±1 contrast labels. The calculation checks directly in all 99 local basis states that the write maps the incoming states into this five-dimensional memory space times the binary endpoint occupations, and that the full tagged term maps each of the three component spaces into itself. Blank terms on the other seams preserve binary loads by exchange. A separate one-cell calculation evolves the full eleven-state memory and all 78 possible two-load configurations, including all twelve doubles, and compares the final probabilities and the load-band diagnostic. [Independent full-basis check](test_retained_memory.py)

The load-band diagnostic immediately after the service uses the blank-sector bound projector on the loads and identity on the working memory. This determines which load components the write creates or removes in the reference band at the instant of the service. A subsequent joint Hamiltonian acts on the retained memory. The terminal seam read is available to the patch. The projected spectral quantity uses the state of the full federation.

## Reproduce

From the RER repository root with the model dependencies installed, run

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -m model.experiments.interval_echo_encounter.memory
python -m pytest -q -W error model/experiments/interval_echo_encounter
```
