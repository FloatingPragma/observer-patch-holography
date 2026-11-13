# Complete Fock and field transport

An explicit product of 24 even two-mode implementers transports the complete
192-mode native source, its actual currents, local charges and longitudinal
Gauss dressing. All products, pair terms, quantum references and work records
remain. The inverse factor checks avoid enumerating 2^192 occupations.

A separate full-Fock rotor control verifies the field-dependent connection,
its pair terms, derivative and square, all odd sectors and the entire integer
flux complement. It quantifies why this independent rotor is not obtained by
relabeling the finite longitudinal source. No GeV calibration or pole is added.

Run `python -W error -m model.openwave_fock_field_transport --write`, then
`python -W error -m model.openwave_fock_field_transport --check` separately.
Focused tests are in `test_transport.py`; `PROOF.md` gives the exact domains
and source-equivalence obstruction.
