This package compiles a conditional local classical action from the registered
fundamental, dual and singlet representations. It uses a supplied flat Spin
chart, connection, Higgs doublet and minimal kinetic/potential/Yukawa laws.
Fermion values, conjugates and first derivatives are exterior-algebra
variables. Gauge parameters retain their first and symmetric second jets.

Run from the repository root:

```sh
python3 code/sm_local_action/build_local_action.py --write
python3 code/sm_local_action/verify_local_action.py
python3 -m pytest -q code/sm_local_action/test_local_action.py
```

The producer needs Python's standard library; the independent verifier also
needs SymPy. Tests need pytest. The verifier uses a separate bit-mask exterior
algebra over SymPy Gaussian rationals and never imports the producer. It
reconstructs the action, all local Ward directions and variational currents.
`validate_custody(packet)` checks strict schema, declared scope and input hashes;
its result explicitly does **not** certify mathematical coefficient correctness.
`verify(packet)` performs the full independent mathematical replay.

The receipt stores exact Gaussian coefficients as `[real, imaginary]` canonical
rational strings and exterior words as strictly increasing generator indices.
All 65 coefficient sectors are independent: three gauge normalizations, six
minimal kinetic sectors, two scalar potential coefficients and 54 real/imaginary
Yukawa entries. The 900 generators cover three declared families, two spin
components, field/conjugate pairs and five jet slots. Bosonic test jets are
rational algebraic fixtures. No fitted or measured parameters are consumed.

The local current is the connection derivative of this same matter action. Its
covariant conservation follows on the matter equations, as stated in the owning
paper fragment. Neither a source-emitted action/current nor a quantum functional
measure, physical Spin identification or calibrated Standard Model is asserted.
The kernel-checked Lean module proves generic connection/curvature identities
and a Grassmann sign control; the full action theorem is analytic, with separate
exact finite executable controls. Gauge fixing and ghosts are outside this
classical assembly.

The intended observer attachment is to bounded patches with local state, ports,
readback, records and repair, certified by public execution bundles. This
package supplies the conditional local action for that attachment; it does
not supply those source histories by identifying arbitrary jets with records.
