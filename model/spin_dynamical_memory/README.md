# Dynamical response from the spin source

The selected spin source determines an exact response equation with a memory
term that retains every complementary excitation. A complete 256-state
certificate measures that term and a noncentral time history using the source's
own ground and clock. Both conserved spin-charge channels have nonzero memory;
the full response retains memory after finite field removal. Bounded resolvent
filters transfer the response equation
to the same canonical continuum histories without an unbounded-moment shortcut.
The [proof](PROOF.md) distinguishes these source results from Lorentz covariance
and surviving continuum scattering, whose dynamical limits are unproved.

```sh
python -W error -m pytest -q model/spin_dynamical_memory
python -W error -m model.spin_dynamical_memory --write
python -W error -m model.spin_dynamical_memory --check
```
