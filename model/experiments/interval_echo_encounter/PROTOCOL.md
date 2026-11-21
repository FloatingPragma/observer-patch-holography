# Protocol for the separated-pair readout comparison

Use the level-zero spherical federation, one observer on cell 0, and the same free two-load sector used in the [two-read encounter](../../evaluation/observer_encounter.py). The preparation is the uniform mixture over all 36 sharp nonadjacent pairs inside that cell. Evolve under the registered contrast interaction at κ = 8 and the exchange-only control at κ = 0. The clock program is idle. The local read windows use supplied free-evolution durations.

At duration one, enroll the lexicographically first internal seam. Compare three services on independently prepared copies.

- Bypass the seam reader.
- Read the two-outcome seam projector and retain both outcomes, including every quiet trial.
- Select the sharp service on half the trials with an independent fair bit. Bypass it on the rest. Keep the service choice and all outcomes.

After four more free units, the cell resolves every internal adjacent pair and the complementary quiet outcome. Retain the unconditional late detection probability on all 36 starts, the sharp yes/quiet contributions, the conditional probabilities, and the difference from bypass. As an audit of the coherences removed, apply a local sign flip to the chosen seam on an independently random fraction of preparations. Fractions one half and one quarter must reproduce the unread sharp and half-serviced terminal distributions, respectively, to absolute error below 10⁻¹⁰ per recorded outcome. The sign comparison is a channel identity on the local binary predicate. The randomized service uses ordinary sharp reads and a supplied fair bit.

For κ = 8, also compute the total weight in the low-energy bound band on the full two-load sector, before the read, conditioned on each sharp outcome and averaged across both. The gap theorem for the shifted Hamiltonian separates energies ≤ −10 from energies ≥ 0. Use the spectral projection below −5. As this projector commutes with free evolution, any change in its ensemble mean across the service measures redistribution by that service. It is a simulator diagnostic. The finite observer can read local seam occupations.

The registered comparisons are as follows.

1. Test whether sharp versus bypass changes the all-start local late detection probability by more than 10⁻⁶, retaining the sign.
2. Test whether the separated ensemble carries less than 5% dressed bound-band weight and the selected click branch exceeds 95%. Test whether the absolute change in the all-outcome ensemble exceeds 10⁻⁶, retaining the sign. Bound the projector approximation error by 10⁻¹⁰ in operator norm and include that bound in the decision.
3. Report whether the unconditional bypass yield for κ = 8 exceeds κ = 0. Local adjacency in the bypass arm is an outcome at one time. A pre-existing durable capture record would require a separate history marker.

Retain failed comparisons. Each distribution contains all 30 internal seams plus quiet. Cross-check the sharp arm against the deposited two-read receipt to 10⁻¹⁰ and independently propagate at least one start through its projected quiet branch. The input preparations, visits, waits, interaction control and significance thresholds are fixed for this experiment.
