# A charged object with its gauge field

The existing native gauge Hamiltonian moves a fermion together with its
link field on a tree. Exact Gauss constraints leave one dressed state at
every vertex, including the reference position. The same Hamiltonian and
complete position detector determine its propagation. The
[proof](PROOF.md) exhausts the invariant sector and gives a volume bound.
The tree, static reference and evolution parameters are specified inputs;
these reference-bound states have no identified physical particle mass.

Reproduce the exact calculations:

```sh
python -m model.native_gauge_strings.certificate --check
python -m pytest -q -W error model/native_gauge_strings
```
