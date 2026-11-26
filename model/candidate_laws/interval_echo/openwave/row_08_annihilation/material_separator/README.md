# Does the collision signal cross an immobile material dimer?

At exchange time 64, two opposite walls prepared on the four-dimer source ring produced **26,854 double clicks in 67,108,864 launches** on the two uninterrupted calorimeters. The ordered sham produced zero. Freezing just the middle material dimer of either guide removes all source dependence from the receiver beyond it. With both frozen, the complete two-receiver distribution is independent of the source. The calorimeters still responded through both material pins to an explicitly added neutral field. [Frozen design](PROTOCOL.md) · [Complete word probabilities and integer counts](counts.npz) · [Decoded results](results.json)

OpenWave lists pair annihilation to radiation as its criterion and positronium branching as a follow-on. M10 has the upstream development mark 🚧. The local collision exhibits opposite phase defects, recovery of source order and transport of energy through polarizable matter. The complete criterion is Partial. Its proposed conjugate inputs give opposite dipole readings and zero common-mode material charge. The transported daughter modes have an independently derived 31.5-unit gap in the continued native guide. The material-isolation test addresses whether an independent outgoing channel reaches these same receivers. [Official criterion and decisions](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#simplest-test-per-criterion) · [Incoming and outgoing charges](../wall_charge/README.md) · [Daughter spectrum](../daughter_modes/README.md)

## A continuous collision with both receivers attached

Each wing has three dimers. Shuttering only the middle dimer's own exchange fixes its orientation. Both adjacent links, the source attachment, the far dimer and the calorimeter coupling remain active. The pin starts without an orientation measurement. The source, guides and calorimeters evolve without intermediate observation. At the final read each launch yields four ring bits, the two middle orientations and both meter bits. The rows below are wall-minus-sham differences from seeded draws of the calculated finite probabilities. Each sham and wall cohort contains $2^{26}$ copies. A separate dimer hand read time 64.00056. [Native apparatus](apparatus.py)

| Frozen middle dimer | Excess left clicks | Excess right clicks | Excess double clicks | Ordered ring on wall launches, calculated |
| --- | ---: | ---: | ---: | ---: |
| Neither | 0.026423 | 0.026474 | 0.0004002 | 6.0358% |
| Left | −0.0000010 | 0.0001488 | 0 | 0.2531% |
| Right | 0.0001470 | 0.0000006 | 0 | 0.2531% |
| Both | −0.0000005 | −0.0000020 | 0 | 0.0348% |

The unblocked source loaded the two meters with 1.720 excess energy units at this appointment. Its initial wall surplus was 65 on all four fixtures. The pinned-side contrasts and every bilateral joint-bin contrast passed the 0.001 gate. The **remaining open guide missed its 0.001 gate on both unilateral fixtures**. Its calculated marginal contrast was 0.00014612. Changing one guide also changes the collision at the ring, as its order probabilities show. Thus subtracting the open and pinned signals cannot assign a numerical fraction of the original energy to individual trajectories.

For each fixed middle orientation the far dimer and meter form a separate Hamiltonian factor. Summing the retained orientations gives the protected marginal for any source-local intervention without outcome selection. This also explains why a selected final ring word is a different comparison. At any vertex separator on the full native contact graph, continuously pinning the individual ternary loads factors the compressed generator, including on nonuniform material with coherent written memories. The source mapping inspects all 33 full local columns for either endpoint at each fixed load. A quiet port at the terminal appointment alone cannot establish that it stayed pinned. [Proof, native columns and finite-sector theorem](PROOF.md)

## A neutral field with the same distal material and meter

The sensitivity control adds one neutral two-state field register at each of the six seam addresses along each physical route. Its reciprocal XX chain has a gapless particle–hole branch in its half-filled uniform continuation. Each source flips the field register nearest the ring. The far field register couples to the original far-dimer swap. The original far dimer, calorimeter gap 32.51538, material–meter coupling and detector bit remain active. The field hopping is tuned to put one finite mode in the detector band. The two sources and their energies are supplied independently. This arm starts in the conditional far-dimer ground with the middle orientation fixed. [Finite field and detector operator](neutral_field.py)

At time 64, across the four initial pairs of middle orientations, the source-minus-sham probability for both meters to click was **0.005834–0.005846**. At time 128 it was **0.070507–0.070564**. For the two zero orientations, interrupting only the left field at 64 left marginal excesses of 0.0000001 on the left and 0.076471 on the right. Interrupting only the right gave 0.076477 and 0.0000005. Cutting both field routes gives exact source independence. The field pulses do 22.7147 units of work per arm. This demonstrates sensitivity to energy carried independently of the middle material exchange. It supplies neither an emission mechanism on the native collision nor an electromagnetic identification for the comparison field.

The integer banks are seeded draws from the calculated native and comparison probabilities. The simultaneous 99% bounded-difference allowance for the stated source-minus-sham proportions is 0.000427 under that sampler. The exact zeros follow from the source algebra. An independent sparse propagation of both unilateral native fixtures reproduced all primary word probabilities within $6.0\times10^{-12}$. [Calculation and decoder](measure.py) · [Independent checks](test_separator.py)

From the repository root, replay the archived counts and tests with

~~~sh
python -m model.candidate_laws.interval_echo.openwave.row_08_annihilation.material_separator.measure
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_08_annihilation/material_separator
(cd Lean && lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_08_annihilation/material_separator/PinnedRead.lean)
~~~

Use the acquisition option to regenerate the probabilities and sample banks with the installed NumPy sampler.
