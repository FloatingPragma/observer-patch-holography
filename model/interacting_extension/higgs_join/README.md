# Interacting leptons and a Higgs field

The literal charged-lepton Higgs vertex acts with mobile fermions, two native
60-state gauge links, electric energy and a colored loop on complete finite
Gauss spaces of dimensions 438 and 1,596. The interaction changes the propagator
and the response of six reference-compensated addition probes. The graph,
scalar cap, family couplings and prepared state are declared controls; no
physical mass or continuum limit is identified.

[Proof and scope](PROOF.md) · [constructor](interacting_carrier.py) ·
[complete response](run_joined_response.py)

From the repository root:

```sh
python -m model.interacting_extension.higgs_join.interacting_carrier --out /tmp/oph-joined-operators
python -m model.interacting_extension.higgs_join.run_joined_response --output /tmp/oph-joined-response
python -m model.interacting_extension.higgs_join.precision_audit --response-root /tmp/oph-joined-response
python -m model.interacting_extension.higgs_join.audit_exact_join --output /tmp/oph-joined-audit
python -m model.interacting_extension.higgs_join.certificate --check
python -m pytest -q model/interacting_extension/higgs_join
```

The constructor checks immutable source bytes from r2062, so a descendant
checkout with the same inputs is admissible. Large eigensystem files are
produced on demand. The retained compact packet contains every exact basis
vector and operator numerator, together with sparse normalized operators.
The certificate replays the exact source construction and the complete response.
Its numerical acceptance windows are toleranced checks, not interval proofs.
