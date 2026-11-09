# Source tests for particle calibrations

Complete source calculations test five inputs to the particle atlas. A
common unit cancels from source energy ratios, and an event-cost identity
selects a finite rank-three band without a dimension window. The phase and
attenuation tests expose the additional dynamics needed for a physical
mass prediction. [Proofs and reproducible calculations](results.json) retain
every complement and distinguish selected coefficients from particles.

A fixed construction is sufficient: the carrier, repair code, memory,
scheduler, preparation and observer controls are design choices. No theorem
forcing those choices or proving the uniqueness of OPH is a prerequisite.
The calculation must derive the complete observer output of that declared
program. A setting chosen from a measured mass is a calibration for that
mass; an independent output can be tested against a held-out measurement.
Here “source-derived” describes that computation, not uniqueness of its
hardware or software.

| Atlas input | Exact result | Observer output or empirical connection to establish |
| --- | --- | --- |
| Coherent phase | The free zero-momentum sea is phase-flat. A supplied nonzero momentum can favor an interior phase. The complete conserved sector of the finite Higgs vertex favors zero phase under energy minimization. | Specify the phase register and its coupling/preparation, then compute the common-source output. Removing a mass calibration requires an independently fixed setting or a derived response. |
| Attenuation | Free-sea minimization at every finite positive-weighted momentum menu selects r=1. The same Higgs-sector test selects r=1. | The retained loss/controller dynamics must produce a nontrivial amplitude without inserting three measured masses. |
| Capacity response | The literal reset C-duty controller has an exact finite-horizon response. At every n>=2 its ordered coefficient ratios are below 1171058/6633 and 8616/11. | A physical capacity-to-program law and physical observable attachment. Changing this controller's archive budget alone cannot give charged-lepton ratios under the linear mass dictionary. |
| Clock and units | The three exact native occupation rates have a common H=-log P. Common time conversion and block length cancel from their ratios. | Attach the charged action and its identified particles to a common operational clock; determine the energy conversion for GeV values. |
| Families and chirality | Per-event cost selects the smallest positive rank-three band. Full-Fock gauge-invariant dressing gives exact canonical fermions on a 240-state quantum-link carrier. | Fix the cooling and matter program, then establish its physical family output, spatial chirality, mirror decoupling and controlled interacting limit. |

The phase and attenuation are alternative inputs to overlapping charged
programs; the ten atlas cards are not ten independent scalar fit parameters.
None of these tests promotes the calibrated Standard Model comparison to
a source-derived physical spectrum. No physical mass prediction is exported.

The proofs are [phase and attenuation](PROOF_PHASE_ATTENUATION.md),
[capacity response](PROOF_CAPACITY_RESPONSE.md),
[clock and units](PROOF_CLOCK_UNITS.md), and
[family cost and canonical dressing](PROOF_FAMILIES_CHIRALITY.md).
Each states the existing finite results it consumes and the actual source
domain. The finite constructions use different declared sources; their
conjunction is not an intertwiner between them.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_calibration_selection.run --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_calibration_selection.run --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/particle_calibration_selection
```

The result file binds imported scientific inputs and package code/proofs.
Its decimal gap values are displays of rational enclosures. Measured
particle values enter no source selector; the capacity comparison accepts
them only after the independent source inequality has been computed.
