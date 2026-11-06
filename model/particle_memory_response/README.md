# Memory corrections to family dynamics

The native record memory has an exact sum over every temporal correlation.
It produces a coherent family potential together with positive diffusion.
Keeping both terms gives a corrected reduced semigroup with a rational
finite-time error bound, including an untouched quantum reference. Exact
three-ray reduction replaces long lag scans and large saved transfer arrays.
The result uses stationary product preparation and an unread output; physical
masses and complete observation tapes require separate attachments.

[The proof](PROOF.md) states the result and its domain.
[The receipt](results.json) binds exact source identities and numerical
convergence checks.

```sh
python3 -m model.particle_memory_response.certificate --check
python3 -m pytest -q model/particle_memory_response
```
