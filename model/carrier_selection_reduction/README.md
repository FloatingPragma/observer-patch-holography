# Selecting twelve directions

Six directional channels are the minimum needed to read every quadratic
response in a supplied three-dimensional inner-product space. Among those
readouts, maximizing the weakest measured response forces equal channel
weights and the six icosahedral axes. Exposing both signs as inverse ports
gives twelve directions with degree five, thirty edges and twenty triangular
faces. A finite tolerance in the information criterion suffices for exactly
the same incidence. Native selection of the response space and the
measurement objective requires a separate attachment.

The [proof](PROOF.md) states the weighted theorem, its stability bound and
the operational premises. The [exact receipt](results.json) records the
attaining configuration, complete sign classification and premise controls.
The certificate uses SymPy exact algebra and rational inequalities; it does
not use measured particle data or a numerical search.

Run from the repository root:

```sh
.venv/bin/python -m model.carrier_selection_reduction.certificate --check
.venv/bin/python -m unittest model.carrier_selection_reduction.test_certificate
```

Running the certificate without `--check` writes its receipt. The producer
has no dependency on a frozen campaign or private continuation package.

The comparison with the determinant selector in
[the screen paper](../../paper/screen_microphysics_and_observer_synchronization.tex)
is in the proof. Its shared six-axis uniqueness step is also checked by
[the carrier certificate](../../code/a5_closure/a5_selection_certificate.py).
This package supplies the weighted minimum-channel reduction and the
information-tolerance implication, with explicit controls for dimension,
weights, completeness, optimality, minimality and inverse ports.
