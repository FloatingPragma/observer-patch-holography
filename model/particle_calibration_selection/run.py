"""Replay the first five source-selection tests and bind their evidence."""
import argparse
import hashlib
import json
from pathlib import Path

from . import phase_attenuation, capacity_response, clock_units, families_chirality

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RESULT = HERE/'results.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inputs():
    retained = json.loads((ROOT/'model/particle_candidate_atlas/audit.json').read_text())
    pins = dict(retained['input_sha256'])
    for path, digest in pins.items():
        if sha(ROOT/path) != digest:
            raise ValueError('changed retained scientific source: '+path)
    extras = set(phase_attenuation.SOURCE_INPUTS)|set(families_chirality.SOURCE_INPUTS)
    extras.update(str(p.relative_to(ROOT)) for p in HERE.iterdir()
                  if p.suffix in ('.py', '.md'))
    extras.update(('model/native_modular_dynamics/PROOF.md',
                   'model/native_quantum_reconstruction/certificate.py',
                   'model/native_clock_comparison/PROOF.md',
                   'Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean'))
    for path in extras:
        pins[path] = sha(ROOT/path)
    return dict(sorted(pins.items()))


def build():
    pins = inputs()
    value = {
        'schema': 'oph.particle_calibration_selection.v1',
        'phase_and_attenuation': phase_attenuation.certificate(),
        'capacity_response': capacity_response.certificate(),
        'clock_units': clock_units.certificate(),
        'families_chirality': families_chirality.certificate(),
        'source_selected_physical_masses': [],
        'physical_calibrations_eliminated': [],
        'restricted_inputs_eliminated': [
            'common clock/energy unit and integer block length from same-source gap ratios',
            'dimension window from the smallest-positive-band calculation for the declared near-rotation cost'],
        'input_sha256': pins,
    }
    if any(sha(ROOT/path) != digest for path, digest in pins.items()):
        raise ValueError('scientific inputs changed during computation')
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = build()
    if args.write:
        RESULT.write_text(result)
    elif RESULT.read_text() != result:
        raise ValueError('source-selection results differ from their complete replay')
    print('PASS: five source-input tests; physical mass predictions remain empty')


if __name__ == '__main__':
    main()
