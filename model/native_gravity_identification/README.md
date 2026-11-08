# Entropy and gravity on the fair source

The growing fair source has an exact regional entropy split and a controlled
modular first law, using every allowed local occupation sector. Its canonical
reference also has a strict finite entropy maximum. The [proof](PROOF.md)
shows that this reference has trivial modular flow on number-preserving
local observables, while native repair changes those observables. A fixed
native cut also fails the required edge-entropy variation. A different
state and responsive geometry require their own source attachment before
these finite results can imply Einstein dynamics.

```sh
python -m model.native_gravity_identification.certificate --check
python -m pytest -q -W error model/native_gravity_identification
```

The [certificate](results.json) evaluates local modular data at source levels
0, 6 and 200 without constructing the global Hilbert space. Its direct source
hashes and finite controls do not claim a replay of every imported producer.
The local preparation error enters a modular-flow estimate only with a
positive support-eigenvalue bound; no archive factorization or physical
clock conversion is inferred.
