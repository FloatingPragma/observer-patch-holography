# Can a negative material defect emit two independently detectable neutral products?

**A single vacancy in a filled Interval Echo host cannot produce a charged daughter plus a neutral click in a separate, initially pristine receiver.** For two separated neutral receivers to click together, the parent must initially contain at least two doubles and three vacancies. The threshold is sharp for material records. A source at exactly this inventory has a nonzero joint arrival amplitude in two distinct neighboring owners on all three open federations. The charged parent in that construction is a prepared load arrangement. An independently identified mass, outgoing flavor and lifetime are still needed for the physical decay assignment. [Observer records and exact coefficients](results.json) · [Executable](measure.py)

OpenWave lists muon decay as untested for M10 and asks whether a muon relaxes to an electron plus neutral ejecta. Its wording specifies neutral ejecta without requiring two independent detectors. The two-detector version here tests the more specific leptonic decay interpretation. For a single neutral receiver the corresponding lower bound is one double and two vacancies. [Official criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#simplest-test-per-criterion)

## What the receiving observers can tell

Prepare one load per receiver port and blank receiving collars. Each receiver records its terminal port loads against a supplied zero-load reference. A terminal contrast echo gives the exact load 0, 1 or 2 on a sharp input. The read can be done on disjoint fresh reference rails after free flight. A neutral click requires **both** a signed sum of zero relative to the filled host and a word different from the original all-ones word. A separate group of ports reads the remaining negative charge. Zero charge on its own always accepts the unperturbed receiver. This is a material charge and an ideal calibrated terminal read. [Signed-load reader](../../row_01_charge/measure.py)

More generally, take any local positive detector effect supported in the receiver's zero-charge sector and having zero probability on its filled, blank state. Such an effect also tests for a departure that could be detected coherently. With all collars blank there is only one all-ones basis vector. Positivity implies that the effect annihilates this vector and has no matrix elements to it. Initial correlations or a superposition in the source therefore cannot hide a joint click in the forbidden sector.

## Why there is a minimum inventory

1. Each freely evolving seam preserves the unordered counts of 0, 1 and 2 at its ends. It also preserves the seam's own collar class, blank or written with absolute contrast 0, 1 or 2. In the native write matrix put \(A=\sqrt{2}W\), \(B=2I-AA^T\) and write \(G\) for the generator on written memory and its endpoints. The integer certificate verifies \(A^T(2G)B=0\) separately for the exchange and contrast parts. The two diagonal blocks of eight times the physical term are \(2A^T(2G)A\) and \(B(2G)B\). Across the two blocks, all 180 nonzero exchange entries and all 90 nonzero contrast entries preserve both quantities. All 99 input columns participate, including initially written source collars. Other seams leave a particular collar untouched, and the registered degree terms are diagonal on port loads. The same conservation laws therefore hold for the full generator and its unitary at every time. [Native seam](../../../../contrast_memory_spherical_v1/law.py) · [Earlier histogram result](../../row_01_charge/README.md#can-the-law-select-a-unique-elementary-amount)
2. In any accepted neutral receiver word the numbers of doubles and vacancies agree. A nontrivial word must therefore contain at least one of each. Disjoint clicks require distinct doubles. If the charged daughter retains \(p_d\) doubles and there are \(r\) disjoint neutral clicks, the parent has \(p\ge p_d+r\) doubles. Its net charge −1 requires \(p+1\) vacancies. Applying the bound to every word proves the same zero probability for coherent sources and neutral effects as above. The Lean lemmas prove the per-receiver implication and finite disjoint counting bound, with arbitrary unobserved ports. [Formal proof](../../../../../../Lean/SourceLaw/NeutralEjecta.lean)
3. A two-stage leptonic cascade would use at least four doubles beyond the final charged daughter's inventory **if** all four neutral products are resolved in separate pristine patches and retained without recycling. A detector in an initially written medium can carry an internal excitation with all loads equal to one. Such a detector requires its own source and sham comparison. A coincidence measurement between load and energy on a preloaded detector has a different support. [Neutral detector with an internal frequency](../../row_03_clock/written_coin/README.md) · [Neutral products measured on prepared wings](../../row_09_neutral/README.md)

## Does the registered law allow the threshold to be reached?

In owner zero, place two doubles at its first two outgoing gateways, a vacancy just inside each, and a surplus vacancy at a third site. All other ports and collars start in their ordinary filled, blank state. In the terminal word the surplus vacancy is still at its start address. Each receiving owner holds one vacancy at its gateway and one double on the next internal port. The rest of the federation reads its host word. This is a terminal word with two separately measurable neutral records. Whole-board reconstruction for this certificate uses the record of every owner. The inclusive coincidence sums probabilities of this and other accepted words.

The verifier enumerates all different-load swaps within four exchanges of the initial and final words, counting all shortest paths where the two enumerations meet. The shortest length is eight on each of the 20-, 80- and 320-cell federations. There are 280 shortest paths in each case. An individual arm can transport its adjacent pair in four exchanges in two orders. The two arms contribute \(2^2 {8\choose4}=280\) combined orders. The full-board enumeration checks for other shorter routes. Each off-diagonal swap has amplitude −1/2 in the native generator. Its diagonal terms first contribute at a higher order. Consequently, for this **specific sharp word** in native time units,

\[
P(t)=\left(\frac{280}{2^8 8!}\right)^2t^{16}+O(t^{18})
     =\frac{t^{16}}{1{,}358{,}954{,}496}+O(t^{18}).
\]

Reality of the finite Hamiltonian and of both words makes the probability even in time. The nonzero leading coefficient proves support at sufficiently short nonzero times. The Taylor coefficient controls the onset. A probability at a finite appointment also depends on higher orders. It is independent of the diagonal contrast strength on the blank sector. The source has neither an observed heavy mass nor a chosen internal frequency. Redistribution of already prepared neutral pairs suffices to pass this two-click test. [All addresses and counts](results.json)

## Finite receiver counts and a false-positive control

A smaller bench keeps the two four-port paths and the actual degree references. All other contacts are closed. The surplus vacancy stays on a separate isolated port. Each path evolves with its own three native seam terms on all 81 load words. The bench records neither, left only, right only, or both at times 1/4, 1, 4, 16 and 64. Each setting has 1,048,576 synthetic launches. The table gives the exposure-1 panel. [Full distributions and all exposures](results.json)

| Law and preparation | Predicted double-click probability | Simulated double clicks |
| --- | ---: | ---: |
| Native, zero or one seeded pair | 0 | 0 |
| Native, two seeded pairs | 4.55 × 10⁻²¹ | 0 |
| Native, two seeded pairs, both bridges cut | 0 | 0 |
| Pair-converting control, no seeded pair | 0.00173448 | 1,832 |
| Pair-converting control, one seeded pair | 0.02058272 | 21,511 |
| Pair-converting control, two seeded pairs | 0.10913788 | 114,874 |

The control adds the previously tested symmetric conversion between \(|11\rangle\) and \((|02\rangle+|20\rangle)/\sqrt{2}\), with strength eight on every enabled seam. The alphabet, routes, load conservation and terminal detectors are fixed. Its unseeded receiving medium can produce a coincident signal of its own. The count banks distinguish the predicted zeros from a missing read. Their sensitivity is too low to observe the small positive two-pair probability on the native bench. [Changed generator](../../row_01_charge/measure.py)

## What the rest-energy ratios would predict in the receivers

Suppose the light particle is a single vacancy and a parent in one fixed negative-unit composition is measured above the same blank, unit-filled host. At a parent/light ratio of at least 206.7682827 or 3477, the native formation-energy bounds require at least 83 or 1,387 defect ports in the parent. Its signed load −1 then fixes the inventory to at least **41 or 693 doubles**, with one more vacancy than double. This arithmetic uses the prior spectrum certificate. It applies to formation energies with that reference and to these input ratios. [Spectral proof](../README.md#can-an-excited-charged-defect-evade-the-ground-energy-bound) · [Spectral receipt](../excited_bounds.json)

If the final charged daughter is that single vacancy, it carries zero doubles. If the rest of the board returns to the filled host apart from two separately observed neutral patches, every double must appear in those patches, accompanied there by a vacancy. The bound holds in every accepted terminal word of a fixed-composition source.

| Parent/light ratio threshold | Minimum pairs across both neutral patches | Minimum pairs in the larger patch | Minimum ports in that patch |
| --- | ---: | ---: | ---: |
| 206.7682827 | 41 | 21 | 42 |
| 3477 | 693 | 347 | 694 |

A single 12-port owner reading zero net load can hold at most six pairs. If each receiver is just one owner as in the short-time test above, **at least 29 or 681 doubles must remain elsewhere** after the single-vacancy daughter is detected. Observers can check that residual inventory directly. Larger detector patches can instead measure how many pairs each neutral patch contains without assuming a compact neutral product. A dressed charged daughter can retain some of the doubles. A prior preparation of written material calls for its own formation-energy bound. Testing the predicted counts requires a separate rest-energy measurement and a decay preparation. [Executable cross-check](measure.py) · [Two-receiver counting lemmas](../../../../../../Lean/SourceLaw/NeutralEjecta.lean)

## What this says about the charged-lepton triple

On the seven-cell comparison with the lepton mass ratios, the installed mass operator and the free charge detector conserve each of the three family characters. On the matched two-cell written-field comparison, supplying contact energy can yield neutral excitation and a lower charged mode, while all six measured separated channels miss the isolated parent energy. This experiment tests a third, native possibility and fixes the minimum source content for passive neutral detection in a pristine material host. Identifying a muon would require selecting the **same** charged parent by its energy and independently measuring a lighter charged daughter and neutral ejecta in free evolution. Identifying its weak flavor also requires an interaction tag. The short-time coefficient predicts a specific material word without making either particle assignment. [Guided character measurements](../../electron_field/lepton_family/identity_audit/README.md) · [Written-field conversion](../../electron_field/lepton_family/decay/README.md) · [Charged-current selection rule](../../../../../../Lean/Screen/ChargedCurrentDictionary.lean)

From the repository root, replay the certificate and synthetic experiments with

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.neutral_ejecta.measure
(cd Lean && lake env lean -j1 SourceLaw/NeutralEjecta.lean)
```
