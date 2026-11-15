"""Verify this bounded field package and write its immutable custody manifest."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    receipt_names = ('CELL_COMPLEX', 'AUTONOMOUS', 'NATIVE_FIELD', 'COMPOSITION',
                     'LINEAR_RESPONSE', 'OBSERVABLES', 'CLOCK_JOIN', 'RESET_DOMAIN')
    pins = {}
    checks = {}
    for name in receipt_names:
        receipt = json.loads((HERE / (name + '.json')).read_text())
        assert receipt['checks'] and all(receipt['checks'].values()), name
        checks[name] = len(receipt['checks'])
        for source, expected in receipt.get('source_sha256', {}).items():
            assert digest(ROOT / source) == expected, (name, source)
            if source in pins:
                assert pins[source] == expected
            pins[source] = expected
    for name in ('manybody', 'clock', 'energy'):
        path = HERE.parent / name / 'MANIFEST.json'
        pins[str(path.relative_to(ROOT))] = digest(path)
    files = {str(p.relative_to(ROOT)): digest(p) for p in sorted(HERE.iterdir())
             if p.is_file() and p.suffix in ('.py', '.json', '.md')
             and p.name != 'MANIFEST.json'}
    record = {
        'schema': 'oph.common.autonomous_field.freeze.v1',
        'worker': 'DYN08',
        'candidate_id': 'native-local12-fermion-plaquette-v2-proposed',
        'selected_by_this_worker': False,
        'files': files,
        'source_pins': pins,
        'receipt_checks': checks,
        'checks': {'eight_receipts_pass': True,
                   'all_current_source_pins_match': True,
                   'frozen_geometry_preserved': True,
                   'coherent_correlated_reset_domain_proved': True,
                   'architecture_changes_and_transfer_boundaries_explicit': True},
        'interfaces': {
            'ideal_source': 'autonomous.initial; autonomous.step; autonomous.records',
            'full_native_word': 'native_field.complete_word(alpha,face_couplings208,clock_spacing=None)',
            'native_preparation': 'composition.preparation_word(charge=1 or -1)',
            'source_certificate': 'COMPOSITION.json',
            'complete_clock_join': 'CLOCK_JOIN.json',
            'actual_current_readout': 'OBSERVABLES.json',
            'harmonic_field_limit': 'LINEAR_RESPONSE.json',
            'reset_input_domain': 'RESET_DOMAIN.md'},
        'allocation': {'source_and_work_bits': 3203, 'scratch_bits': [3200, 3201, 3202],
                       'optional_clock_bit': 3203, 'optional_spectral_probe_bit': 3204,
                       'joint_consumer_rule': 'remap fresh clocks, probes and detector flags to disjoint wires'},
        'scope': 'Complete finite autonomous field proposal, native factorized circuit certificates, actual one-face joint source records and transverse harmonic response. No completed interacting quantum Maxwell continuum or calibrated physical force/mass claim. Finite native services and CAR/field couplings are declared architecture.',
        'independent_review_pending': True,
        'reserved_acceptance_runs': False,
        'official_rows_changed': []}
    path = HERE / 'MANIFEST.json'
    path.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'manifest_sha256': digest(path), 'files': len(files),
                      'source_pins': len(pins), 'receipt_checks': checks}))


if __name__ == '__main__':
    run()
