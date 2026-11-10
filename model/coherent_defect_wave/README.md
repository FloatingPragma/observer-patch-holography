# Coherent defect-wave observer

A finite unitary clock couples a movable regional flux to a coherent carrier
wave in both directions. The program uses all thirty central-link signs,
preserves a fixed boundary charge and an interior Gauss constraint, and prepares
its own stationary reference. Complete detectors distinguish the reciprocal
responses and preserve interference until a recorded measurement intervenes.
The construction does not identify a physical particle or supply a GeV mass.

The [particle path](PARTICLE_PATH.md) connects these finite results to the
common dependency map. Its executable extensions include a
[prepared pure stationary background](background/PROOF.md),
[complete same-U frequency readout](spectroscopy/PROOF.md),
[connected spatial circuit](spatial/PROOF.md) and
[shared-link occupation law](species/PROOF.md). Their distinct observation
domains are explicit. The tested spatial spectator band is quadratic at
small momentum; the primitive triplet copies have identical spectra.
Neither result establishes a relativistic species or three lepton masses.

```sh
python -m model.coherent_defect_wave --check
python -m model.coherent_defect_wave.background --check
python -m model.coherent_defect_wave.spectroscopy --check
python -m model.coherent_defect_wave.spatial --check
python -m model.coherent_defect_wave.species --check
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
