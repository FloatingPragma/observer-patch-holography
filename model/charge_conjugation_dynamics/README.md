# Symmetry of particle and antiparticle class counts

The finite native repair law has an exact conjugation symmetry: complement
each load about five, reverse the class and integer-lift signs, and complement
every repair coin. An unbiased, conjugation-symmetric ensemble therefore has
zero mean signed class imbalance at every finite time. Individual histories
can have an imbalance. The complete paired-history calculation cancels that
sampling fluctuation exactly.

The certificate checks all 2,592 local branches and executes paired histories
on the actual closed federation, including no-operation and rejected records.
A generic [Lean theorem](../../Lean/Dynamics/ChargeConjugation.lean) proves
that an equivariant finite kernel preserves even distributions and gives zero
mean for every odd observable. The executable supplies its native premises.
The [proof](PROOF.md) also bounds finite bias and symmetric postselection.

Class conjugation is a property of the declared classical construction. These
counts do not identify baryon number or determine a baryon-to-photon ratio.
Physical baryogenesis requires a source, dictionary and asymmetry mechanism
with their own complete records.

```sh
python3 -m model.charge_conjugation_dynamics --check
python3 -m pytest -q model/charge_conjugation_dynamics
cd Lean && lake env lean Dynamics/ChargeConjugation.lean
```
