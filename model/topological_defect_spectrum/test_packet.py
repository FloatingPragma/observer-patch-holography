from copy import deepcopy
import hashlib

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_cold_receipt_retains_all_spectra_and_complete_source_custody(replay):
    stored = strict_json.load(packet.OUTPUT)
    packet.validate(stored, replay)
    expected = {str(path.relative_to(packet.ROOT)) for path in packet.OUTPUT.parent.glob('*.py')}
    expected.update(packet.PACKAGE+'/'+name for name in
                    ('README.md', 'PROOF.md', 'CORE_PROOF.md', 'spectral_witness.json'))
    expected.add(packet.PRIOR)
    inherited = strict_json.load(packet.ROOT/packet.PRIOR)['source_bindings']
    assert expected | set(inherited) <= set(stored['source_bindings'])
    for path, digest in stored['source_bindings'].items():
        assert hashlib.sha256((packet.ROOT/path).read_bytes()).hexdigest() == digest
    assert stored['core']['controlled_27_site_cube']['full_Fock_dimension'] == 2**54
    assert stored['scope']['two_to_54_Fock_matrix_executed'] is False


def test_conditional_scope_and_spectral_mutations_are_rejected(replay):
    for key, value in replay['scope'].items():
        wrong = deepcopy(replay)
        wrong['scope'][key] = not value
        with pytest.raises(ValueError):
            packet.validate(wrong, replay)
    for section, key, changed in (
            ('algebra', 'linear_texture_zero_multiplicity', 1),
            ('bulk', 'Coulomb_general_signed_multiplicity', 'n^2'),
            ('lattice', 'full_exact_lattice_zero_mode_dimension', 4),
            ('sampling', 'additional_lattice_gap_levels_excluded', True),
            ('core', 'virtual_nonzero_modes_retained_in_error_bound', False)):
        wrong = deepcopy(replay)
        wrong[section][key] = changed
        with pytest.raises(ValueError):
            packet.validate(wrong, replay)
    wrong = deepcopy(replay)
    wrong['lattice_bounds']['positive_levels'].pop()
    with pytest.raises(ValueError):
        packet.validate(wrong, replay)
    for key, changed in (('all_pass', 1), ('physical_comparison_count_added', False)):
        wrong = deepcopy(replay)
        wrong[key] = changed
        with pytest.raises(ValueError):
            packet.validate(wrong, replay)


def test_receipt_custody_requires_proofs_and_exact_rational_witness(replay):
    for name in ('PROOF.md', 'CORE_PROOF.md', 'spectral_witness.json'):
        wrong = deepcopy(replay)
        del wrong['source_bindings'][packet.PACKAGE+'/'+name]
        with pytest.raises(ValueError, match='manifest'):
            packet.validate(wrong)
    wrong = deepcopy(replay)
    wrong['source_bindings'][packet.PRIOR] = '0'*64
    with pytest.raises(ValueError, match='manifest'):
        packet.validate(wrong)


def test_inherited_manifest_verified_without_reexecuting_parent(monkeypatch):
    def forbidden():
        raise AssertionError('an inherited proof was reexecuted')
    monkeypatch.setattr(packet.parent_packet, 'report', forbidden)
    monkeypatch.setattr(packet.parent_packet, '_produce', forbidden)
    assert packet.inherited_custody()['every_inherited_source_digest_verified'] is True
    prior = strict_json.load(packet.ROOT/packet.PRIOR)
    prior['source_bindings'][next(iter(prior['source_bindings']))] = '0'*64
    monkeypatch.setattr(packet.strict_json, 'load', lambda path: prior)
    with pytest.raises(ValueError, match='manifest'):
        packet.inherited_custody()
