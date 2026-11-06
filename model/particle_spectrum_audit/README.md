# Particle masses and comparison conventions

The declared circulant and six-power branch yields twelve mass parameters.
Reading the charged coefficients under an explicit pole-mass premise uses
the native values directly; the running-to-pole conversion is a separate
declared comparison attachment. All twelve values
lie outside the printed PDG marginal intervals. The charged ratios nearly
agree, while the light-quark ratios exclude a common normalization repair.
The branch supplies no particle abundance distribution or proof of physical
spectrum identity. Its architecture is target-informed.

[The proof](PROOF.md) gives the inversion and comparison rules.
[The receipt](results.json) contains all mass values, conventions and residuals.
[The reference data](pdg2026_reference.json) identify the PDG primary tables,
uncertainties and distinct mass observables.

From the repository root:

```sh
python3 -m model.particle_spectrum_audit.derive --check
python3 -m pytest -q model/particle_spectrum_audit
```

The default check independently evaluates the three-register charged law and
comparison arithmetic. The input packet and multiloop outputs are explicit
frozen inputs. A separate executable replays all five light-quark outputs,
the top, Higgs, W and Z, the two criticality conditions, and the reference
substitution controls against a built SMDR v1.3 source tree:

```sh
python3 -m model.particle_spectrum_audit.replay_smdr --smdr-dir /path/to/SMDR
```

The source pin checks every C/header file before compilation. The supplied
tree must contain `libsmdr.a`, `lib3vil.a` and `libtsil.a`. The retained
perturbative approximations have no certified total theory-error bound.
The electroweak decoupling adapter includes both corrections to the pinned
implementation: the threshold-function operator and the placement of the
strong-coupling assignment after the electroweak terms. The adapter retains
the upstream copyright notice and [GPL license](LICENSE.SMDR.txt).

The neutral gap-ratio comparison consumes the hypotheses of the
[finite response bound](../particle_finite_response/PROOF.md). The absolute
neutrino eigenmasses are unmeasured. Zero masses in the minimal declared
action give zero oscillation gaps, contrary to the nonzero measured gaps.
