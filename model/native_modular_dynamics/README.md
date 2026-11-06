# Modular histories from native counts

The source-counted complete-attempt transfer `P` determines the faithful
reference `rho_L=P^L/Tr(P^L)` and the generator `H=-log P`. This package
extends its static count reconstruction to ordered modular correlation
functions. Exact rational polynomials approximate `P^(it)` with a certified
binomial tail; each resulting correlation is a finite weighted sum of native
closed-word probabilities. The common clock normalization is `t=L s`.
Here `U(t)=P^(it)`; conventional Hamiltonian Heisenberg time is `-L s`.
Every rejection, zero event and normalizing return probability is retained.

The [proof](PROOF.md) gives both approximation and data-error bounds,
including coefficient amplification and conditioning. The
[certificate](certificate.py) attaches those bounds to the existing
current-only observer at a finite attempt cost. These are reconstructed
statistics, not a claim that the native process executes modular unitaries.
Finite-sample estimation error and physical time calibration are separate.

Coordinate sectors inherit faithful conditional references with exact nested
conditioning. A proper nonempty coordinate corner cannot intertwine the
original modular flow. Completing the reference with a specified record
partition gives compatible flows within its sectors, while changing the
global reference. This obstruction does not cover arbitrary tensor regional
algebras.

The [interaction test](interaction.py) strengthens the finite field
diagnosis. A four-configuration contrast vanishes for every quadratic
operator in the native occupation dictionary, even with arbitrary hopping
and on-site coefficients. It is nonzero for both `I-P` and the full
`-log P`, with an exact operator-norm separation. All 66 states are retained.
The construction does not establish a chiral or relativistic continuum.

```sh
python3 -m model.native_modular_dynamics.certificate --check
python3 -m pytest -q -W error model/native_modular_dynamics
```
