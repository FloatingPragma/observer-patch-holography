# Measuring electric and magnetic loop strengths

Complete color-loop records determine the ratio of electric and magnetic
strengths in the specified finite Higgs system. The calculation evolves all
438 and 1,596 Gauss states and retains every Higgs outcome. The Higgs
contribution cancels from a useful short-time coefficient, while its
interaction changes the sixth-order response. Six finite clock settings
give rigorous bounds around the source's bare ratio one. Connecting this
controlled loop measurement to a growing physical color field requires a
native detector construction and a continuum comparison.

[Proof and premises](PROOF.md) · [exact source and readout](source.py) ·
[reproducible certificate](results.json)

From the repository root:

```sh
python -m model.color_loop_readout --write model/color_loop_readout/results.json
python -m model.color_loop_readout --check model/color_loop_readout/results.json
python -m pytest -q -W error model/color_loop_readout
```

The measured finite coefficient is e/|m|. Its value one and the source's
color-Casimir coefficient ratio two are not estimates of the physical
strong coupling. The receipt gives every declared interface, exact
probability enclosure and a finite sampling bound.
