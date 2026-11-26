# Which earlier records did a delivery use?

Two observers on a twenty-owner Interval Echo board independently read the decision, receiving bank, delivered bit, selected continuation and final response. Across three preparation histories and every six-bit input, they agree on all 96 accepted deliveries and all 96 refusals. They also read which earlier writer supplied each necessary input. Rewriting the selected payload with the same value changes the recorded parent in 32 comparisons while preserving every public output. These counts enumerate deterministic simulated measurements on the declared apparatus.

## Apparatus and measurements

The [protocol](PROTOCOL.md) was frozen in `00359861c`. Each input occupies four balanced rail pairs, one for its value and three for its last-writer stamp. Six unconditional appointments set permit, bank, left payload, right payload, onward and spare. In two histories a seventh appointment writes either payload again with the same initial command operand. The native comparison echo reads both rails relative to a reference. A missing or doubled response refuses decoding. The dispatcher consults all six registers. On acceptance, the addressed receiver records presence and reads its selected payload. A continuation contact records its address. The dispatcher stamps this composite record eight. The selected continuation station reads it, writes payload XOR onward and stamps the response nine. Observers at owners 16 and 19 compare the result. The deposited archive retains each sealed mutation, source and continuation readings, the two receiving tapes and the producer’s register/value/writer declarations. All source appointments have empty semantic support because their commanded operands are fixed independently of the board. Local actuation still reads the rails.

At a launch, a subset of registers **forces** the public outcome when every completion agreeing on those registers gives the same outcome. For each of 64 inputs we test all 64 subsets against all 64 physically prepared completions within its history. At every refusal, permit alone forces the outcome. At every acceptance, permit, bank, the chosen payload and onward form its unique least forcing subset. Flipping any one of those bits separately changes a measured outcome. Every larger set containing them is sufficient. Adding spare fails minimality in all 96 acceptances. Dropping onward loses sufficiency in all 96. If the entire diagnostic tape is included as a public result, all six inputs become necessary.

The generic proof in [Support.lean](Support.lean) derives the least-certificate property from sufficiency and individual pivotal witnesses. It proves both hypotheses for the dispatch function at every input, and proves that the continuation depends on its one composite register. The independently decoded [measurements](records.json.gz) check that the operated device realizes the corresponding truth table. The support definition is specific to this observation menu. It allows every other input to vary jointly. For example, if enrollment instead required the left and right payloads to agree, either payload register could determine the delivered bit. The two minimal certificates in that restricted table are incomparable.

## Writer authentication

Each sealed event supplies its written register and the versions it used. A second decoder checks every affected writer stamp and the unchanged values and stamps outside the write. For each certified register, the parent must have written it, the child must name that writer, and its received value must equal the parent's deposited value. Exactly 480 direct parent incidences pass in the 192 nominal launches. On a held equal-payload launch, rejecting an uncertified register, a false write declaration, an obsolete equal-valued writer and a discontinuous value each exercises a different condition of the [formal parent relation](../../../../Lean/ObserverPatchHolography/Provenance/SemanticEventProvenance.lean).

Changing the continuation circuit to onward XOR spare while retaining the original declaration exposes the missing spare in all 32 accepted altered launches. Parking the dispatcher’s comparison echo leaves 24 invalid pair readings and produces a quiet receiver and continuation. These controlled registrations establish PR-81 for this finite producer. Identifying the resulting order with physical signal propagation or spacetime needs additional premises.

Run from the repository root:

```sh
python -O -W error -m model.experiments.interval_echo_consensus.versioned_dispatch.measure
python -m pytest -q model/experiments/interval_echo_consensus/versioned_dispatch
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_consensus/versioned_dispatch/Support.lean)
```
