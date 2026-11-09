# Populated spatial and energy records

A coherent phase scan of the actual eight-cell pair detector has positive
original energy-band weight in every scalar Fourier record. The actual
Wilson matrices control how its band bound varies with the phase.
The complete order-32 phase configuration has a reported joint
spatial/energy union probability greater than one sixteenth, with all
failures and the idle quantum reference retained.
Parseval bounds the eight-row scalar effect by one quarter of the identity;
its original final-excitation union consequently has weight above 15/256.

The [proof](PROOF.md) derives these bounds on the unchanged source state,
fixes complete clock and phase errors, and distinguishes the finite
analyzer from conserved momentum or a rest-mass observable. Its primary
GeV mass table is empty. The independent six-state physical CAR/link
control checks the complete instrument and reference logic; it is not
the primary spatial source.

Run `python -m model.particle_primary_populated_momentum --write`,
`python -m model.particle_primary_populated_momentum --check`, and
`python -m pytest -W error model/particle_primary_populated_momentum`.
