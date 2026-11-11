# Three finite observation laws

The executable gives a complete joint screen-reference interference readout,
an exact four-read counting law on all 28,680 two-occupation configurations,
and a certified measurement/feedback response from a freshly executed
interacting source. Each result retains its failures and original records.
The untouched reference marginal is unchanged; only joint correlations
distinguish the two repair laws.

These are theorems and executions of the declared finite constructions.
Their source services and clocks are mathematical inputs, without a physical
universe identification. The [proof](PROOF.md) specifies the domains and
the [receipt](results.json) contains the complete records and proposed
dashboard wording.

```sh
python3 -m unittest model.finite_observable_closures.test_closures
python3 -m model.finite_observable_closures --write
python3 -m model.finite_observable_closures --check
```
