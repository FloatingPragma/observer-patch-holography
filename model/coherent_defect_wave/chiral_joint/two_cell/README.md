# Complete quantum two-cell matter and scalar response

Both cells retain their full six-mode lepton algebra and five-state quantum
Higgs memory. The connecting quantum link keeps all six active states.
The source has 614400 states and one specified interacting Hamiltonian;
no scalar expectation or classical link replaces either quantum register.

The [exact Gauss projection](BASIS_PROOF.md) gives 1053 neutral physical
states and the complete charged boundary sectors. The
[full-spectrum certificate](SPECTRUM_PROOF.md) proves a unique neutral
ground near -4.554826635 source units and a gap exceeding 0.394804412.
The ground belongs to L=2, so a trajectory beginning in the empty L=0
state cannot substitute for its calculation or preparation.

The [symmetry construction](SYMMETRIES.md) proves exact odd-sector doublets
and maps right addition to reflected left removal with the compensating
boundary. A charged sector's absolute ground can be dark to the selected
probe. The complete compensated right probe has visible addition and
removal lines near 1.05204385 and 0.90979863 source energy units, with
certified positive record weights. Its absolute charged minimum and first
populated detector line are distinct.
The same proof constructs a separate two-link periodic Hamiltonian with
an exact antiperiodic translation and complete momentum-phase records.
Open-edge energies are not imported into that periodic source.

The [complete readout](READOUT_PROOF.md) retains 460 energy-window outcomes
conditioned on ground preparation, or 461 outcomes including its failure.
The [propagation control](PROPAGATION_PROOF.md) separately prepares a
physical pair/scalar state and certifies nonzero motion through the quantum
link, with both arrival and no-arrival records. Removing hopping makes
that arrival exactly zero.

The [refinement construction](REFINEMENT_PROOF.md) translates the periodic
charged source and its retained compensator together, giving physical joint
phase records. A blank second link fails the attempted open-to-periodic
transfer: a compressed zero residual hides an actual omitted-state residual
of squared norm 3/16. Each source therefore needs its own energy calculation
or a controlled full-operator comparison. The six finite dependencies and
physical acceptance tests are recorded in the [particle path](../../PARTICLE_PATH.md#six-finite-dependencies-and-the-physical-refinement).

These are finite quantum energies and observer records. Three-dimensional
relativistic dispersion, physical particle identity and a GeV reference
require their own controlled attachment.

```sh
python -m model.coherent_defect_wave.chiral_joint.two_cell --check
python -m pytest -q model/coherent_defect_wave/chiral_joint/two_cell
```

Replay checks the stored rational eigenbasis proposals with exact arithmetic.
Use `--generate-witnesses` explicitly to replace those proposals, followed
by `--write` to regenerate the receipt.
