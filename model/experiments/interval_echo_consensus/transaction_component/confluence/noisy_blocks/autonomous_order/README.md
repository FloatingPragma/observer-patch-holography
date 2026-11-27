# Repairs choose their next offer from measured material

Two observing stations find **20 different orders of eight repairs from one initial state**. All 6,561 possible sequences of native fault outcomes retire every label once. The material that remains after each repair helps select the next label. From one initial mismatch, the receiver's expected terminal count is **1905/4096**, including all outcomes. [Prospective first experiment](PROTOCOL.md) · [Physical selector follow-up](SCANNER_PROTOCOL.md) · [Sharp branch and mask records](scan_records.json.gz)

## What selects the offer

Eight balanced dimers on owners 18 and 19 record which of the eight offers remain. Each offer covers a pair or triple of positions on one of two repair triangles. A ninth dimer marks the block boundary. Label 7 covers the second triangle's triple. The three remaining dimers on these two owners operate the selector. Each has native orientations 20 and 02. A full seam exchange interchanges the orientations. Two independently routed pointer stations resolve each rail and restore their paths. [Stock and routing](dispatcher.py) · [Independent pulse and pointer checks](test_dispatcher.py)

The first control dimer starts at zero. On each triangle the observers read which count-indexed material word is active. They compare each active position directly with its protected boundary at both stations. Every measured mismatch triggers a swap of the parity dimer. Its two exits arm the ascending or descending scan's seek dimer. The address circuit for both scans is fixed. An armed seeker examines the next material label. A live label retires itself and disarms that seeker with the same detector-controlled pulse signal. This selects the least label on even parity and the greatest on odd parity. Both heads verify retirement. The three control dimers return to zero. Only then is the fresh noise clock interrogated. [Native selector and controls](scanner.py)

The registered launch has flag word 8, with one mismatch on the second triangle. Its first offer is 7. A correcting first outcome clears the mismatch, so the next offer is 0. A preserving outcome leaves one mismatch, so the next offer is 6. Across all 512 combinations of initial material mask and parity, an available label is retired exactly once per successful activation. Empty initial stock produces no offer. Parking retirement keeps the mask, blocks repair, leaves the triangle unchanged and consumes no fault draw. Parking the six parity swaps leaves the same entry mask and six measured flags but changes the first physical offer from 7 to 0. [Complete two-head tapes](scan_records.json.gz) · [Independent decoder](test_dispatcher.py)

The final offer leaves the remaining-label shelf empty. After its repair, the observed empty shelf triggers eight refill pulses and toggles the wrap marker. An independently routed scalar run with sixteen native appointments read the successive orders 7,0,1,2,3,4,5,6 and 7,6,5,4,3,2,1,0. Its 3,532 single-pointer outcomes verify the two wraps and the fault settings selected for that trajectory. The three spare control dimers were separately tested on every branch of the single-block experiment. [Scalar quantum readout](scalar.py) · [Baseline and scalar tapes](records.json.gz)

## Contraction for adaptive order

The conditional law for a newly prepared noise clock is correction with probability 12/16, preservation with 3/16 and corruption with 1/16. All sites in an offer see that same outcome. Write $x$ for one site's zero-or-one mismatch and $r$ for the number of still-unused labels covering it. Before the next independently drawn fault, consider the local potential

$$
M(r,x)=(3/16)^r(x-1/13).
$$

If this site is untouched, both arguments stay fixed. If touched, its conditional next mismatch mean is $(3x+1)/16$ and its remaining exposure count drops by one. In either case the conditional expected potential equals its current value. Each site starts with three remaining exposures and finishes with none. Summing the six terminal means yields

$$
\mathbb{E}[D_{\rm end}\mid\mathcal F_{\rm start}]
  =(27D_{\rm start}+1878)/4096.
$$

This remains valid when the dispatcher uses every prior public record to choose the next unused offer. The finite tree proof permits completely different subsequent scheduling decisions on the three current noise outcomes. Randomized nonanticipating choices also satisfy the conditional identity by averaging over the current selection. With repeated complete blocks the limiting mean is $6/13$. The all-64-input integer recursion checks the exact mean separately for this material policy. Holding a single fault setting across a whole block from six initial mismatches instead gives $3/2$. That control has the same isolated command probabilities and violates the fresh conditional law. [Lean decision-tree and reserve-budget proofs](Adaptive.lean) · [Independent integer recursion](test_dispatcher.py)

At an offer of size $k$ containing $b$ current mismatches, the conditional expected change in total mismatch is $(k-13b)/16$. The eight labels have a total of 18 site exposures. In a strict prefix at least one size-two-or-larger label remains, so the total exposed thus far is at most 16 on each trajectory. Summing the conditional drift bounds gives $\mathbb E[D_t\mid\mathcal F_{\rm start}]\le D_{\rm start}+1$. This gives $L=8$, $A=1$, $\beta=1$, $\lambda=27/4096$ and $\epsilon=939/2048$. The prefix bound is conservative. The previous fixed-order sharp-prefix calculation uses different scheduling hypotheses. [Prior registered process and clock](../README.md)

The four normal forms are the two protected boundary choices on each triangle with all six repair flags clear. The distance counts repair mismatches. A listener can retain the public remaining shelf and wrap marker alongside the repair quotient. The future of the flag word can depend on that shelf before the next block boundary. The bounds condition on the full past. A wrap may move the controller marker even on a normal repair word.

The apparatus produces the order and its block boundaries. Supplied resources include preparation of the initial balanced stock, ongoing isolated activation appointments, within-appointment routing and pulse sequencing, and the projected and restored fresh noise clock. Boundaries are measured in completed offers. There is no demonstrated bound on the elapsed time between successive external activations. This establishes PR-121 for the declared finite material process. The global derivable premise remains Partial. The archived outcomes are exact branches of a simulated native instrument with Born weights. The 6,561 paths and 20 orders are counts of alternatives, with highly unequal weights. [All initial states and complete joint terminal law](measure.py)

From the repository root, replay both acquired banks, the independent tests and the proof with one numerical thread.

~~~sh
pkg=model.experiments.interval_echo_consensus.transaction_component.confluence.noisy_blocks.autonomous_order
folder=model/experiments/interval_echo_consensus/transaction_component/confluence/noisy_blocks/autonomous_order
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
nice -n 15 .venv/bin/python -B -m "$pkg.measure"
nice -n 15 .venv/bin/python -B -m "$pkg.scanner"
nice -n 15 .venv/bin/python -B -m pytest -q "$folder"
(cd Lean && nice -n 15 lake env lean -j1 -DwarningAsError=true "../$folder/Adaptive.lean")
~~~
