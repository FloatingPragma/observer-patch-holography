# Native selected excitation, clock and empty-cell transfer

Official row2 is🚧 under the2026-11-15 PR723 decision. This new private N+
candidate demonstrates source-dependent local mode selection, an actual native
phase clock, and native transport of that clock-bearing source into an empty
cell. A separately frozen blank-ticket refinement gives a bounded frequency
limit and addresses the fixed-word phase alias. No pin or award is changed.

Run the small, entirely offline checks with Python3 and NumPy:

```
python3 -B headline.py
python3 -B check.py
```

`headline.py` executes the fully inspectable261-line fixed-clock source,
167-line refinement and45-line changed-law source and checks their saved
outputs. It fetches nothing. `check.py` imports no producer module and uses
floor/ceiling endpoint arithmetic, independent Fourier formulas, full raw
padding checks and literal work-trace replay. The original fixed run failed
only on JSON conversion of a NumPy boolean; its source/freeze/error remain.
`experiment_serialized.py` changes only that conversion.

The localized basis input is filtered by the actual source transfer, reaching
compact-mode fidelity0.99885947 with complete success probability5.64866e-6.
The same geometry and observer under always-swap repair retain only4/7.
Every failure and all positive/negative phase outcomes remain. The positive
clock state has at least0.99611 fidelity and retains that overlap under holding.
A21-native-event word carries its complete state into an empty identical region,
with the work coin returned zero. The refined record-only reader recovers both
frequency bands; one511keV anchor leaves a limiting pair of1022keV for independent separated copies
and625.844629keV for the second internal source frequency.

Source enumeration, complete finite population composition and sampled histories
are distinguished. There is no Monte Carlo claim of witnessing the rare selected
branch. Triangle reservation and coherent memory services are supplied architecture.
The native endpoint rule determines the transfer, selected mode, phase and transport.
The derived energy generator is not an installed update. The empty vacuum has
zero phase; the charged constant mode also does, so no general charged mass floor
is inferred. Laboratory time, mechanical inertia and ambient confinement remain
outside this finite claim. Scientific/adversarial review and full other-group G6
are still required for an award request.

Start with DESIGN.md, refinement/DESIGN.md, both PROOF.md files and
RESPONSE_TO_REVIEW.md. INTERFACE.md gives the exact source/arc map for a common
candidate; it does not assume that another row already runs on the same state.
ALTERNATIVES.md and pilots/ preserve the failed and impractical routes. The old
2.45085keV excess and80 lower modes retain their original disposition.
