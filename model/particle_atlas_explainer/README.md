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

The opening repair panel binds the
[coherent seam-repair construction](../coherent_seam_repair/PROOF.md): a reversible
seam correction stores exchange disagreement in coherent memory; its declared
address reflection produces a finite wave. Preparation and optional detector
operations use the same walk. A five-level encoded transport changes a complete
finite return response in a static link background. The panel distinguishes
this source from the dissipative committed-record pump and from the illustrative
protected textures. It assigns no physical species or GeV mass to either.

The accepted observer-chart comparison transports the detector with the state
and links. A separate 30-bit link-memory control transports regional flux while
holding its boundary fixed. The panel shows why that boundary condition protects
the regional charge and why it does not supply a joint dynamical gauge field for
the coherent walk. Its [boundary proof](../coherent_seam_repair/PROTECTION.md)
also bounds the finite signed-graph gap without identifying it with rest mass.

The [joint defect-wave construction](../coherent_defect_wave/PROOF.md) couples
those registers in one fixed unitary tick. Its panel shows both directions of
the interaction, the detector change caused by recording a coherent branch,
and the complete same-clock preparation on the full register and the admitted
interior gauge sector. The return apparatus uses the retained link frame and
preserves that constraint. Its finite 72000-dimensional source has no assigned
physical particle species or rest masses.

Every step carries a two-part line. The teal half states what the finite
model computes; the amber half states the link the program has to derive.
The [source catalogue](../particle_candidate_confirmation/CATALOGUE.md) gives
each protected object its label, family, source domain and proof. The six
[confirmation tests](../particle_candidate_confirmation/CONFIRMATION_STEPS.md)
specify the physical evidence required for each named candidate.

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
The [candidate identity tests](../particle_candidate_confirmation/RESULTS.md)
bind the six finite gauge multiplets, the seven-record composite counter and
the dressed charged pair's positive energy window. These detector results
constrain the particle dictionary without identifying a named species.

```sh
python -m model.particle_atlas_explainer.build --write
python -m model.particle_atlas_explainer.build --check
python -m pytest -q model/particle_atlas_explainer
```

`--write` generates both SVG and PNG using `rsvg-convert`; `--check` compares
both outputs. The PNG embeds its source SVG digest. `--slices DIR` writes
1000-unit viewBox slices for visual review. The SVG is standalone and dark;
its animations stop under reduced-motion settings. The renderer consumes
`build.py`, `canvas.py`, `parts.py`, `scenes_a.py`, `scenes_b.py`, `families.py`
and `repair_panels.py`, `joint_panels.py`. The
separate `card.py` is not part of this atlas or its source closure.

The public [candidate atlas](../../docs/research/OPH_PARTICLE_CANDIDATE_ATLAS.svg)
is owned by the leaf renderer below. It retains the numerical atlas produced by
`model.particle_candidate_atlas` with unchanged rendered contents inside its enlarged canvas and
adds the authenticated coherent repair panel, source-family catalogue and confirmation steps. The
immutable numerical producer and its audit retain their scientific domain;
its standalone SVG comparison command does not own this augmented public view.
The leaf audit binds both the numerical rendering and the candidate receipt.

```sh
python -m model.particle_atlas_explainer.candidate_view --write
python -m model.particle_atlas_explainer.candidate_view --check
```
