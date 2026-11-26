# Measuring the quadratic load of a twelve-port word

A receiving observer can infer the sum of the squared loads on someone else's twelve ports without reading the source list. Each trial routes the load of one source port, and all twelve addresses receive equal numbers of trials. A blank contact splits its load between two arms, a private reference acts on the first, and a second contact recombines them. The receiver counts nonempty readings on the second arm. Repeating this for each address costs the same number of contact activations for a six-double source and an all-single source. The inferred square sums are 24 and 12. [Frozen design](PROTOCOL.md) · [Every detector count and signed transaction](records.json)

## Which physical reading gives the square

Two independent runs exchange the roles of cells 0 and 13 on the 20-cell board. A full addressed exchange on a blank tree edge swaps the endpoint loads and the return path restores the intermediate ports. Fresh copies put the selected source load against an empty receiver port. All other couplings are parked. On a load of size q the private zero reference contributes the phase κq²t. A half exchange and its inverse therefore put the load on the other arm with probability (1 − cos(κq²t))/2. The final signed echo against a third, empty port resolves zero, one and two. The decoder pools the two nonzero outcomes. Both readout edges retain blank heads throughout preparation. [Native acquisition and readout](measure.py)

At the first appointment, π/κ, only a single can click. At the second, π/(4κ), a double always clicks and a single clicks with independently calibrated probability a. For pooled rates A and B the single fraction is A and the double fraction is B − aA. The quadratic reading is therefore 12[A + 4(B − aA)]. This decoder takes two receiver rates and a unit calibration. A third appointment at 3π/(8κ) tests a frequency it did not use. A linear private phase and a reference with contrast removed use the same gates and appointment times. [Finite decoding and gap proof](Readback.lean)

## What the two observers compared

Start with six doubles and six vacancies. One equal-load offer gives a refusal. Six successful writes then transfer one unit at a time from a double to a vacancy. Each change has a record giving the event, two addresses and port loads before and after. A directly prepared all-single word provides an independent endpoint. The complete word and repair counts are withheld from the receiver decoder. The unit and empty calibration sources have totals twelve and zero, respectively.

Both receivers measured 24.000 on the initial word. After the sixth success their readings were 12.043 and 11.943, with simultaneous radius 0.423. Each resolved all six drops of two. The refused offer left the initial reading at 24.000. The largest discrepancy on the held appointment was 0.00110 against an allowance of 0.00948. Replacing the private square by a linear phase on the initial word gave 12.004 and 12.026. Removing contrast gave zero. Every reported address and appointment used 65,536 simulated starts. The conservative bounds apply Hoeffding's inequality to independent, equally weighted address samples and include the uncertainty of the unit standard. [Regression and decoder checks](test_readback.py)

For twelve ternary loads with total twelve, each double requires one vacancy. Consequently the exact square sum is 12 plus twice the number of doubles. The all-single word alone reaches 12, and the next allowed value is 14. Actual gate use has another meaning. The selected address needs 25 to 31 routing contacts in these runs, plus two splitter contacts, the same private dwell at a given appointment, and a terminal echo and undo. This resource bill depends on the route. It takes exactly the same values for every source word at that address.

PR-117 remains Partial. The experiment realizes a quadratic phase reading on the ternary, total-twelve sector. The certified lineage also declares a full signed integer charge fiber and a normalized central readback cost. The measured phase scale does not establish a universal energetic or time cost of all possible readouts of that fiber. [Lineage and its distinct accounting conventions](../../../../../paper/screen_microphysics_and_observer_synchronization.pdf)

Run these commands from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_consensus.routed_counts.phase_cost.measure --write
python -m pytest -q model/experiments/interval_echo_consensus/routed_counts/phase_cost
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_consensus/routed_counts/phase_cost/Readback.lean)
~~~
