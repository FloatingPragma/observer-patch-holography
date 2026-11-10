# Coherent defect-wave observer

Executable source programs connect protected link flux, wave propagation
and complete observer records. Their Hamiltonians, preparations and
occupation rules are specified separately. The [particle path](PARTICLE_PATH.md)
identifies which results can be joined and the physical conditions needed
for a measured-particle comparison.

## Massive free particles

The [graded flux field](many_body_route/PROOF.md) encodes every fermionic
occupation in literal Gauss-positive link slots. Its local Wilson clock
has a positive filled-sea reference, particle and antiparticle partners,
and a controlled massive free Dirac limit. The
[complete mass instrument](mass_instrument/PROOF.md) reads rest energy two
in source units, retaining preparation failures and apparatus work.
Grading, regulator and face-count coefficient are architecture choices;
interacting Standard Model identification and an independent GeV scale
require additional source attachments.

## Chiral fields and quantum Higgs response

The [sitewise chiral source](chiral_source/PROOF.md) carries every occupation
of the actual fifteen native charge states with two left-Weyl components.
Its spectral regulator controls a free chiral limit with explicit finite
nonlocality. The [quantum Higgs calculation](chiral_higgs/PROOF.md) preserves
all twelve currents and executes complete scalar/matter energy records.
The [conditional field dictionary](chiral_particle_map/PROOF.md) gives seven
color-resolved Dirac fields and one neutral Weyl field with masses `|y v|`.
The quantum scalar reference does not select that background, the Yukawa
coefficients or measured species; an interacting chiral continuum requires
its own proof.

## An interacting physical reference

The [periodic quantum source](chiral_joint/periodic_spectrum/README.md)
keeps every active lepton occupation and both quantum scalar/link memories.
Its unique physical ground supports equal charged addition/removal lines
at 0.97009613 source units, with all 2627 preparation/readout alternatives
retained. An [actual Higgs composite](chiral_joint/higgs_composites/PROOF.md)
excites the same line without a fitted scalar background. The
[particle path](PARTICLE_PATH.md) separates this finite result from physical
mass: source-relative dispersion, the three-dimensional interacting
reference and an independent energy scale require their own control.

## Native interactions and charges

The [full quantum-link source](chiral_joint/peter_weyl_links/README.md)
transports neutral packets and all five native matter representations while
retaining their references and complete detector records. Four matter hops
also produce a contractible gauge-loop record, without an inserted magnetic
term. The [complete local ground](chiral_joint/local_reference/README.md)
is a spin-one-half doublet with gap about 0.08955533 source units. The
[particle path](PARTICLE_PATH.md) specifies the many-cell reference and
dispersion needed to interpret such finite responses as physical masses.

The [native-routing join](native_routing_join/PROOF.md) couples flux and
native update on one clock. The
[propagating completion](propagating_native_join/PROOF.md) moves native waves
with controlled free kinetic windows. Its contact bound concerns one flux
and one native occupation with a preserved finite native momentum window;
its density vertex fails the scalar-boost test. An
[actual native occupation cloud](field_identity_audit/PROOF.md) carries the
fifteen fermion-family weights and admits current-preserving bound transport.
The idle repair fixes this cloud, giving no native mass correction.
Propagation and Higgs interactions in the chiral source are declared
additional laws.

## Moving the protected object

The [mobile-flux source](mobile_flux/PROOF.md) adds a flat orbit and transports
complete inside/exterior flux pairs between cells. Five face locations do
not identify five observed species. A
[two-lane route](mobile_flux/ROUTING_COMPLETION.md) gives the same object a
controlled one-object Dirac window while retaining every collision state.
The full graded field above specifies a different occupation law and
positive reference; their interaction and state data cannot be transferred
without an explicit intertwiner.

## Finite energy and collision controls

A [positive Hamiltonian](energy/PROOF.md) built from native repair and hopping
has constrained grounds and an exact six-state two-wave collision with
complete records. Its [populated triplet bands](dispersion/PROOF.md) fail the
common relativistic dispersion test. The
[cubic instrument](cubic_readout/PROOF.md) executes complete records at 8 and
64 momenta on a different circuit clock. These controls identify what the
specified laws produce, without treating an expectation value as a mass.

The original joined observer also has a
[prepared pure background](background/PROOF.md),
[complete same-U readout](spectroscopy/PROOF.md),
[connected spatial circuit](spatial/PROOF.md) and
[shared-link occupation law](species/PROOF.md). The source contract below
refers to that fixed-boundary program.

```sh
python -m model.coherent_defect_wave --check
python -m model.coherent_defect_wave.background --check
python -m model.coherent_defect_wave.spectroscopy --check
python -m model.coherent_defect_wave.spatial --check
python -m model.coherent_defect_wave.species --check
python -m model.coherent_defect_wave.energy --check
python -m model.coherent_defect_wave.dispersion --check
python -m model.coherent_defect_wave.mobile_flux --check
python -m model.coherent_defect_wave.many_body_route --check
python -m model.coherent_defect_wave.mass_instrument --check
python -m model.coherent_defect_wave.propagating_native_join --check
python -m model.coherent_defect_wave.field_identity_audit --check
python -m model.coherent_defect_wave.chiral_source --check
python -m model.coherent_defect_wave.chiral_regulator --check
python -m model.coherent_defect_wave.chiral_higgs --check
python -m model.coherent_defect_wave.chiral_particle_map --check
python -m pytest -q model/coherent_defect_wave
```

## Source, state, reference and clock

The source is the pinned twelve-port carrier with thirty link bits, a five-level
internal code, a two-state coherent repair record and sixty oriented seam
addresses. The fixed-boundary, one-regional-defect sector has ten link words.
The joint Hilbert dimension is 72000. A sparse basis key is
`(link_word, port, internal, parity, address)`; coefficients are exact rational or
algebraic numbers divided by a shared square-root normalization.

One tick is `U=C4 C3 C2 C1 C0 B`: a link-dependent wave advance B, followed by
five occupation-controlled local defect transports in sorted carrier-edge order.
The wave step reads every link sign. A transport moves the defect only if its
incident flux permits movement and the coherent wave occupies an endpoint.
The link memory carries the alternatives coherently. Idle execution does not
measure or archive a moved/held flag; doing so would replace U by a channel.

The [proof](PROOF.md) establishes the unitary inverse, conserved boundary and
internal charges, and reciprocal interaction. The [same-U preparation](PREPARATION.md)
uses a complete heralded filter to prepare a faithful stationary reduced
reference with success probability `281/1080`. This reference is distinct from
the pure operational controls used to exhibit interaction. The retained
preparation archive determines a larger, conditionable state.

The [interior constraint](GAUSS.md) has a 36000-dimensional positive sector.
The clock and the complete regional-flux and wave-basis detectors preserve it.
A raw link-word detector is also available as an ambient-register apparatus;
it dephases gauge partners and is outside that constrained readout grammar.
No claim transfers this one central constraint to a complete continuous gauge
theory or to all boundary charges.

## Executing the controls

```python
from model.coherent_defect_wave import JointObserver
observer = JointObserver()
state = observer.inverse_wave(observer.embedded_ports(1, (5, 6)))
output = observer.step(state)
for record, branch in observer.region_face_measurement(output):
    print(record, branch.norm_squared())
```

The example yields two flux positions with probabilities one half. An inverse
controlled-transport apparatus can recombine the alternatives. Its complete
wave-return detector succeeds with probability one; inserting the complete
flux-face measurement first lowers that success probability to one half.
This return apparatus belongs to the ambient register experiment: its
port-superposition projector need not preserve the interior Gauss constraint.
A separate conditional-frame return apparatus in `gauss.py` preserves that
constraint and reproduces the same one-versus-one-half contrast on a normalized
constrained preparation. It reads no raw link label and retains both outcomes.
A wave-basis detector also distinguishes two different link inputs with the
same prepared wave, with probabilities `1/120` and zero. These are finite
operational interaction and coherence controls, not particle-scattering data.

[The exact receipt](results.json) binds the implementation and its source
imports, records every detector complement, and contains the separate declared
architecture and physical-scope boundaries. Regenerate it with `--write` after
a source change. No empirical particle inputs enter the construction.
