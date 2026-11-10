# Mobile paired flux

An executable boundary-changing extension of the coherent source. It adds a
flat link sector and moves complete protected inside/exterior flux pairs
between carriers while preserving Gauss constraints and total pair count.

Run from the repository root:

```sh
python -m pytest model/coherent_defect_wave/mobile_flux
python -m model.coherent_defect_wave.mobile_flux --write
python -m model.coherent_defect_wave.mobile_flux --check
```

See [PROOF.md](PROOF.md) for the literal link-domain change, preparation,
complete records, clock and scope. The separate
[protected-object routing construction](ROUTING_COMPLETION.md) gives a local
two-lane program whose one-object code has a controlled free Dirac window.
Its value `m=2` follows from the chosen unit-cost face-count coin; it is not a
GeV prediction or the original source clock's mass. `results.json` records
exact probabilities and source hashes. The source identifier is
`coherent-mobile-paired-flux-v1`; previous fixed-boundary receipts do not
automatically describe it.
