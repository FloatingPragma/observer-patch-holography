# Finite metric response and capacity histories

A matter disturbance changes the lapse and causal-cone readouts of the
specified reciprocal scalar source. Its exact lapse response is
`-99/200` times its area response, with the complete thermal update and all
rejected proposals retained. A conditional regional Einstein theorem also
turns entropy/area error budgets directly into exterior Newton-force error
bounds. The capacity density law supplies normalized expansion histories
without choosing an absolute cosmic capacity.

The [proof](PROOF.md) separates these three source domains. The finite
metric response has no spatial gravitational field; the force bound requires
the Einstein theorem's matched regional premises, and the cosmological
history requires its density and Friedmann identifications. No measured
gravitational constant or numerical cosmological constant is derived.

```sh
python3 -m model.openwave_gravity_capacity --write
python3 -m model.openwave_gravity_capacity --check
python3 -m pytest -q -W error model/openwave_gravity_capacity
cd Lean
lake env lean ../model/openwave_gravity_capacity/Transfer.lean
```

The receipt regenerates exact source spectra and response intervals. The
tests independently evolve the complete source transition matrix, compare
proper-time propagators, and integrate the radial error source. The Lean
file checks the scalar transfer identities using the existing regional and
capacity theorem modules; its kernel command is separate from the receipt.
