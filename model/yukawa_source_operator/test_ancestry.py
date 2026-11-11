"""Empirical-input isolation controls for the complete native calculation."""
import json
from pathlib import Path

from . import certificate, operator, three_point


def test_target_mutation_does_not_change_source_output(monkeypatch, tmp_path):
    target_names = {'reference_values.json', 'experimental_inputs.json',
                    'calibrated_matter_parameters.json', 'calibrated_neutrino_parameters.json'}
    original = Path.read_text
    reads = []
    payload = {'electron_mass': 1e99, 'top_Yukawa': -7, 'muon_mass': 0}
    def guarded_read(path, *args, **kwargs):
        reads.append(path.name)
        if path.name in target_names:
            return json.dumps(payload)
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'read_text', guarded_read)
    # Recompute this package's exact contraction and trace objects each time.
    operator.tensor_lines.cache_clear()
    three_point.source_data.cache_clear()
    three_point.source_moments.cache_clear()
    first = certificate.run()
    payload.update(electron_mass=-1e99, top_Yukawa=1e30, muon_mass=1e-99)
    operator.tensor_lines.cache_clear()
    three_point.source_data.cache_clear()
    three_point.source_moments.cache_clear()
    second = certificate.run()
    assert first == second
    assert not set(reads) & target_names
    assert first['empirical_mass_or_Yukawa_inputs'] == []
