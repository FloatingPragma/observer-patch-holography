# From a twist to a particle

The [explainer atlas](../../docs/research/OPH_PARTICLE_ATLAS_V3.svg) follows
an illustrative topological defect in nine numbered steps: the
compare, repair and record loop; a twist owned by a whole loop; why no local
operation removes it; its core, charge, spin and antiparticle; motion as
re-establishment one place over; position as a written record; velocity as
the slope of successive records; the proposed connection between standing
repair cost and rest energy; and the field description observers build from
agreed records. The chain combines distinct finite examples and physical
interpretation targets. It is not an established derivation of the Standard
Model spectrum. A candidate dictionary, measured roster and the ten-node
[dependency map](../particle_joint_reconstruction/DEPENDENCIES.md) close the page.

Every step carries a two-part line. The teal half states what the finite
model computes; the amber half states the link the program has to derive.
The [horizontal atlas](../../docs/research/OPH_PARTICLE_ATLAS_V2.html) keeps
every family chain, calibration path and proof link.

The carrier drawings use the certified icosahedral coordinates and the
recorded creation and transport run of
[holonomy transport](../particle_observer_outputs/HOLONOMY_TRANSPORT.md).
The detector inset reads the two-sector click law of the
[defect attachment](../particle_defect_attachment/PROOF.md). Node statuses
come from the dependency register. Reference masses are CODATA 2022 and
PDG 2026 comparison data; the builder rejects any input that declares them
predictions. The evidence panel reads the
[candidate correspondence](../standard_model_source_map/candidates/CANDIDATES.md)
and isolated gauge-transition receipt, checks all their producer hashes, and
distinguishes a rare primary-source transition and static fusion from an
executed contact-interaction control in a separate supplied action. The GeV
mass table is empty; the reference masses do not supply source parameters.

```sh
python -m model.particle_atlas_explainer.build --write
python -m model.particle_atlas_explainer.build --check
python -m pytest -q model/particle_atlas_explainer
```

`--write` generates both SVG and PNG using `rsvg-convert`; `--check` compares
both outputs. The PNG embeds its source SVG digest. `--slices DIR` writes
1000-unit viewBox slices for visual review. The SVG is standalone and dark;
its animations stop under reduced-motion settings. The renderer consumes
`build.py`, `canvas.py`, `parts.py`, `scenes_a.py` and `scenes_b.py`. The
separate `card.py` is not part of this atlas or its source closure.
