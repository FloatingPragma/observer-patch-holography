from fractions import Fraction as F
import subprocess
import sys

from . import packet


def test_complete_detector_receipt_and_transitive_source_custody():
    value = packet._produce()
    packet.validate(value)
    assert value['formal_native_pulses_per_preparation_and_detection'] == 65822721
    assert F(value['correct_native_label_probability_lower']) > F(99999999, 10**8)
    assert value['scope']['operational_model_Weyl_flavors_identified'] == 8
    assert 'model/particle_species_detector/physical_targets.json' in value['source_bindings']
    assert 'model/certified_reference/certificate.py' in value['source_bindings']
    assert 'model/observer/source_scalar/compiler.py' in value['source_bindings']


def test_replay_rejects_supplied_bad_receipt(tmp_path):
    path = tmp_path/'bad.json'
    path.write_text('{"schema":"wrong"}')
    result = subprocess.run([sys.executable, '-W', 'error', '-m', 'model.particle_species_detector', '--check', str(path)],
                            cwd=packet.ROOT, capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert 'Receipt keys mismatch' in result.stderr
