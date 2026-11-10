# Coherent seam repair

A seam can correct disagreement while keeping the removed information in a
coherent memory. This finite observer combines that reversible repair with an
address reflection, giving a wave recurrence on the actual twelve-port carrier.
The same walk prepares a faithful stationary reference. Optional detectors read
complex amplitudes, and a five-level fibre carries a static gauge-loop signal.
The construction supplies no physical particle identification or mass prediction.

The named program is `coherent-seam-repair-v1`. Its declared complex-state,
coherent-memory and apparatus operations differ from the committed coin/current
observer in `model/observer/native.py`. The latter keeps its original semantics.
The [proof](PROOF.md) fixes each source, preparation, reference and clock; the
[receipt](results.json) regenerates its exact finite checks from public sources.

```sh
python -m model.coherent_seam_repair --check
python -m model.observer.coherent_repair --check
python -m pytest -q model/coherent_seam_repair/test_repair.py
```

Use `--write` to regenerate `results.json`. SymPy and NumPy are required. The
implementation keeps the 1440- and 7200-component amplitudes as finite vectors;
its proofs use local matrices and connection blocks rather than dense walk
matrices.

`Session(observer, state, max_events=...)` retains every operation and detector
outcome up to its declared record capacity. It returns both measurement branches
and reports capacity exhaustion. Its wave-step counter counts idle W operations;
apparatus and preparation requests have their own records.

```python
from fractions import Fraction
from model.observer.coherent_repair import Observer

observer = Observer()
state = observer.embedded_port(0)
state = observer.step(observer.step(state))
assert observer.source_populations(state)[0] == Fraction(643, 900)
outcomes = observer.embedded_return_measurement(state, 0)
# Both labelled, unnormalized detector branches are retained.
```

`Observer(fibre=5, central_edge=(0,1))` executes the declared central-link
background. The code does not evolve that background or derive a many-particle
interaction. Antisymmetric seam amplitudes do not by themselves establish
fermionic statistics. The [physical boundary](PROOF.md#protection-and-physical-boundaries)
distinguishes conserved internal charge, removable flux and finite wave gaps.

The [finite chart-agreement construction](PREMISES.md) transports the apparatus
with the source and proves equality of complete adaptive history probabilities.
A separate [thirty-bit link program](PROTECTION.md) realizes a movable regional
flux protected by a fixed boundary. Its evolving links have no composed
stationary-reference construction with the static-link walk.
