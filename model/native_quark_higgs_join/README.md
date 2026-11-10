# Quark and Higgs response

The native quantum Higgs operator and an elementary quark removal have a
strictly positive shared contribution in an energy bin of the same
interacting eight-cell source. Complete interference records retain
preparation failure, pump failure and the outside-bin outcome. The
[proof](PROOF.md) derives the exact operators, every local line and the
full-source error bounds. Energies use the source clock; the result does
not identify physical quark masses or asymptotic colored particles.

Replay with `python -B -W error -m model.native_quark_higgs_join --check`.
Run the eight tests with
`python -B -W error -m unittest model.native_quark_higgs_join.test_join`.

After energy-bin acceptance, the actual prepared elementary and composite
states have common-history distance below `5.193e-6`, including their
retained preparation environment. Their original success rates differ;
the comparison excludes reading the pump-setting label.
