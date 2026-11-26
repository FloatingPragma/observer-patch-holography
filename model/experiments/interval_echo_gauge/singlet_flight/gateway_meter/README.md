# What do the three gateway clocks and the rooted loop measure?

All **1,260 gateway pairs** on the 20-, 80- and 320-owner Interval Echo boards return the predicted signed triplet overlaps. The largest error is 0.00113. The same three clocks on the half-exchange material miss a pair overlap by 0.44250. A separate measurement on all **25,200 rooted internal faces** distinguishes the odd and even exchange preparations perfectly at native contrast. These two tests and the exact finite identities below fulfill PR-150 on its specified boards and ordered material sector. [Frozen protocol](PROTOCOL.md) · [Seeded detector counts](counts.json.gz) · [Simulator and decoder](measure.py)

## Three time comparisons select the triplet

Each owner opens its thirty internal blank contacts and parks bridges and private reference holds. On one load their sum is five times the contrast coefficient plus half of the degree-five cell Laplacian. The actual isolated-cell generator includes one additional, spatially constant private-reference energy. Thus both preparations have the same energy projectors. Shutting only the bridges on a live federation can leave a nonuniform gateway potential. The address preparation, phase reference, gated time superpositions and terminal clocks here are commissioned apparatus. [Isolated-cell and mobile measurements](../README.md)

At each appointment a balanced two-path clock compares a free exposure of duration τ with a closed hold, using the selected frequency as reference. Its bright output applies the average of the identity and the phase-referenced free evolution. The four prior adjacency eigenvalues are 5, √5, −1 and −√5. Select the √5 level. For each of the other three energies use τ equal to π divided by its distance from the selected energy. The three bright outcomes pass the selected level with certainty, and each unwanted level has one dark appointment with certainty. The tape keeps the first dark outcome and the three-bright outcome. A sharp port gives its diagonal matrix entry. The difference in three-bright frequencies for the normalized sum and difference of two ports, divided by two, gives their signed real entry.

For the 12-port adjacency matrix A and its antipodal permutation J, write B = I − J and C = BA. Integer arithmetic on the registered cell verifies B² = 2B, BC = CB = 2C, C² = 10B and AC = 5B. The symmetric matrix P = (B + C/√5)/4 therefore obeys P² = P and AP = √5 P. Its trace is three. On the odd space A² = 5I. On the even space (I+J)(A+I) is the all-ones matrix, which gives eigenvalues 5 and −1. The positive and negative odd eigenspaces each have dimension three by the two traces. P is consequently the entire first nonconstant spectral projector. Its port-column Gram matrix equals P. Every stored internal seam has been checked, including owners with antipodal gateways. [Exact product identity](GatewayRead.lean)

Every sharp-port prediction was within 0.00131. The 1,260 distinct pairs divide into 128 antipodal entries of −1/4, 136 entries of +1/(4√5), and 996 entries of −1/(4√5). A simultaneous 99% Hoeffding allowance for 7,560 measured bright fractions is 0.00260 rounded up, and also bounds each signed half-difference. Halving every internal exchange while retaining its contrast energy, endpoints and clock settings gives a maximum diagonal error of 0.36210. These outcomes distinguish frequency selection from a detector that labels geometry independently of the flight law.

## A loop outcome bounded by its vacancy source

On a blank addressed seam a native exposure π exchanges the material exactly. The contrast cost contributes an integer multiple of 2π and the antisymmetric exchange amplitude changes sign. Three such exposures around (r,a,b,r) permute a and b and restore r. The source reader first tests whether the vacancy is at r. On success the dark port comparing the closed triangle with an idle path selects the antisymmetric excess-address combination. Every other excess address is fixed by the triangle. For amplitudes u and v at a and b the outcome weight is |u−v|²/2. This is at most |u|²+|v|²; summing over any spectator reference gives the quadratic-form bound 0 ≤ Aᵣ,ₐᵦ ≤ Pᵣ on the complete ordered sector. The root gate is a recorded part of the instrument. [Reference-inclusive inequality](GatewayRead.lean) · [Protected exchanges](../../loop_transport/HOLONOMY.md)

The even source gave zero dark outcomes, and the odd source gave 26,424,115,200 dark outcomes on that many trials across all rooted faces. A sharp excess gave dark frequency 0.500002. A vacancy away from the root was always refused. Under the half-exchange control the even, odd and sharp sources all give dark probability 5/8. Each count is simulated and carries its original launch denominator. This result identifies the local loop measurement used in the finite vacancy-variance theorem. A propagating independent colour field remains a separate question.

Run from the repository root with numerical threads set to one.

~~~sh
python -m model.experiments.interval_echo_gauge.singlet_flight.gateway_meter.measure
python -m pytest -q model/experiments/interval_echo_gauge/singlet_flight/gateway_meter
~~~
