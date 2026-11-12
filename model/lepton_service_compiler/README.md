# Preparing and reading a complete lepton cycle

The sixty-address lepton source has an explicit implementation in pulses of
the six-level native processor. Its preparation, analyzer, inverse tick,
four return outcomes and source-preserving reset are compiled together.
One, two and three cycles retain the complete quantum archives and every
return history. An exact conjugation identity makes the large analyzer
cheap to evaluate without replacing its compiled source preparation.

[Proof](PROOF.md), [executable](service.py) and [receipt](results.json)
specify the supplied coherent loading, blank reset, routing and readout
services. Their autonomous implementation and physical precision require
separate source constructions.

Run `python3 -m model.lepton_service_compiler --check` from the repository
root. The full source replay takes several minutes.
