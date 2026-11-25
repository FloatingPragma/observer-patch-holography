# Finite theorem checks for the dark-sector assessment

Run the following two commands from the RER Lean directory with its registered toolchain.

```sh
lake env lean ObserverPatchHolography/EinsteinBranch/DarkSector.lean
lake env lean ../openwave/interval_echo/row_19/check_theorems.lean
```

The first elaborates the source of the dark-sector results. The second inspects the axioms of the finite split, collar localization, explicit finite expectation, recovery switch-off, total-source response, comoving dilution and deep-profile mass law. It also checks the exact fourth-power normalization equivalence, the ten-observer decoder, the finite-radius tensor bound and convergence of the tensor residual. The mathematical dependencies and the physical assumptions are different questions. The complete definitions and theorem statements retain the declared boost, bulk Gibbs condition where used, recovery bound, universal geometric coupling where used, conserved charge for cold scaling, and the profile hypotheses. The RER register lists the physical realization of anomalous collar stress at the ordinary coupling as open. [Dark-sector paper](../../../cosmology/oph_dark_matter_paper.tex) · [Finite split](../../../Lean/ObserverPatchHolography/EinsteinBranch/FiniteCapGeneratorSplit.lean) · [Decision and observations](DECISION.md)

Both commands completed successfully. Lean reported only propext, Classical.choice and Quot.sound for all eleven inspected declarations. The signed-area and disk readbacks were also replayed. Premise PR-100, physical gravitating collar stress at the ordinary coupling, is classified Open in the status register. [Lean output](theorem_receipt.txt) · [Status register](../../../docs/registers/OBSERVABLES.json)

The quantitative Einstein theorem measures matter, entropy and fixed-volume area on common proper balls. Its calibrated limiting coefficient gives the ordinary Newton coupling. The tensor bound controls the bulk, area and stationarity errors divided by the fourth power of the same proper radius. The signed collar transfer must supply those measurements or a physically identified directional consequence. [Matched small-ball theorem](../../../Lean/ObserverPatchHolography/EinsteinBranch/QuantitativeSmallBall.lean) · [Einstein tower dependency index](../../../Lean/docs/EINSTEIN_BRANCH_INDEX.md)
