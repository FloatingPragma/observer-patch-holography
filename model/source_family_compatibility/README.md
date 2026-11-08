# A common clock for spin and fermionic source sectors

The full-ground strong-coupling theorem imposes an exact compatibility
constraint on particle reconstruction. Throughout its stated domain,
keeping the ordered spin field's ballistic clock nonzero sends every
occupation-charged source excitation to infinite energy. A chemical
reference shift cannot retain both conjugate charges on that clock.
This applies to the specified projection family, not to every OPH source.

The [proof](PROOF.md) gives the complete energy-weight bound, source grading,
an exact finite countercontrol with both spin and odd CAR excitations,
and the precise finite gauge, chirality and statistics attachments.
It also replays an existing gauge-invariant experiment whose two admitted
operation parameters give different complete records. No new Hamiltonian
or fitted mass is introduced.

The [repair-clock construction](REPAIR_TRAJECTORY.md) supplies an exact
hole-position detector on the full source while preserving its internal
spin information. Its complete repair records give calculable dephasing
and diffusion; a growing global-attempt schedule bounds their disturbance
by `4/n` while retaining the entire archive. Two executed native programs also have equal record capacity
and resource counts but different coherent responses. A specified program
has a definite response; these capacity data alone do not select it.

The [full interacting defect budget](DEFECT_BUDGET.md) removes the need for
strict one-hole invariance in the repair comparison. It bounds the complete
source/reference/archive error by occupation-defect expectations along the
same ideal Hubbard evolution. The positive interaction energy controls
those expectations, including actual doublon--hole creation, without a
global strong-coupling projection premise. Energy-changing controls require
their own work budget. The action, initial state and physical clock remain
separate selection premises.

At fixed couplings, the [local full-source construction](LOCAL_MATTER_READOUT.md)
retains odd fermion probes and interacting spin observables on one centered
ground and clock. Complete local records have an explicit repair and boundary
error for the original Weyl action with its bare repair; distant raw repair
values are inaccessible to that protocol. No alignment-sector invariance
or frame covariance of that repair is assumed.
[Original charged energy bands](CHARGED_SPECTRAL_BANDS.md) have positive
volume-independent weight for both addition and removal. The exact dimer
also generates its own composite response clock. These are source energies
and interacting local responses, without a relativistic particle assignment.

The separate [aligned-hole sector](ALIGNED_HOLE_DISPERSION.md) supplies an
exact vacancy-position detector, magnetic translation momentum and a
controlled free nonrelativistic continuum. Its reference gap and kinetic
curvature follow from the signed source action. This aligned reference is
not the interacting full ground; its gap is not a physical rest mass.

```sh
python -W error -m pytest -q model/source_family_compatibility
python -W error -m model.source_family_compatibility --check model/source_family_compatibility/results.json
```

An interacting relativistic limit requires a source-selected scaling and
controlled spatial spectral branches on this common field. The fixed-spacing
construction does not identify its centered state with the ordered spin
limit or derive the Standard Model's chiral gauge representations.
