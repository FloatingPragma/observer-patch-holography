from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import linalg

from ..centered_particle_spectrum import packet as ground_packet
from . import certificate, packet, source


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_complete_search_keeps_weak_lines_and_source_energy_sum_rules(result):
    spectrum = result['spectrum']
    assert spectrum['complete_block_lines'] == spectrum['strictly_positive_lines'] == 138
    assert spectrum['weight_positivity_unresolved_lines'] == 0
    assert spectrum['weight_outside_reducing_block'] == '0'
    assert min(F(row['weight_interval'][0]) for row in spectrum['rows']) > 0
    assert min(F(row['weight_interval'][1]) for row in spectrum['rows']) < F(1, 10**12)
    check = spectrum['independent_sum_rules']
    for left, right in (('spectral_weight_interval', 'direct_AstarA_interval'),
                        ('spectral_first_energy_moment_interval', 'direct_Astar_commutator_interval')):
        a, b = map(F, check[left])
        c, d = map(F, check[right])
        assert max(a, c) <= min(b, d)
    assert check['source_rows'] == [4900, 1568]
    assert spectrum['source_selection']['all_linear_CAR_responses_orthogonal'] is True


def test_replay_accepts_only_integer_proposals_without_numerical_solver(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('cold verifier used an eigensolver or QR proposer')
    monkeypatch.setattr(linalg, 'eigh', forbidden)
    monkeypatch.setattr(linalg, 'qr', forbidden)
    source._block.cache_clear()
    source.parent.reduced.cache_clear()
    value = certificate.verify(packet.load_witness())
    assert value['strictly_positive_lines'] == 138


@pytest.mark.parametrize('mutation', ['energy', 'vector', 'omitted', 'float'])
def test_wrong_or_incomplete_spectral_proposals_fail(mutation):
    witness = packet.load_witness()
    if mutation == 'energy':
        witness['energies'][0] -= 1 << certificate.ENERGY_BITS
    elif mutation == 'vector':
        witness['vectors'][:, 0] = 0
    elif mutation == 'omitted':
        witness['vectors'] = witness['vectors'][:, :-1]
    else:
        witness['energies'] = witness['energies'].astype(float)
    with pytest.raises((ArithmeticError, ValueError)):
        certificate.eigensystem(witness)


@pytest.mark.parametrize('columns', [[], (True,), (0,)*138])
def test_bad_projector_columns_are_rejected_after_cache_warmup(columns):
    witness = packet.load_witness()
    source.block(tuple(map(int, witness['columns'])))
    with pytest.raises(ValueError):
        source.block(columns)


def test_zero_probe_is_rejected_instead_of_claiming_particles(monkeypatch):
    zero = 0*source.probe()
    monkeypatch.setattr(source, 'probe', lambda: zero)
    witness = packet.load_witness()
    old = ground_packet.load_witness()
    g = source.parent.reduced('ground', tuple(map(int, old['ground_columns'])))
    b = source.block(tuple(map(int, witness['columns'])))
    with pytest.raises(ArithmeticError, match='identically zero'):
        source.selection_certificate(g, b)


def test_direct_configuration_space_moments_match_numerical_control(result):
    old = ground_packet.load_witness()
    b = source.parent.reduced('ground', tuple(map(int, old['ground_columns'])))
    g = b['full_basis'] @ (old['ground_vectors'][:, 0].astype(float)/(2.0**60))
    g /= np.linalg.norm(g)
    a = source.probe()
    q = a @ g
    first = q @ ((source.parent.centered(5, 2) @ q -
                   a @ (source.parent.centered(4, 4) @ g))/16)
    for interval, value in ((result['spectrum']['total_weight_interval'], q @ q),
                             (result['spectrum']['first_energy_moment_interval'], first)):
        lo, hi = map(float, interval)
        assert lo <= value <= hi


@pytest.mark.parametrize('flag', [
    'all_dark_charged_states_exhausted', 'physical_particle_species_identified',
    'finite_gaps_identified_as_relativistic_rest_masses', 'native_detector_executed',
    'autonomous_localized_defect_proved', 'continuum_or_gravity_conclusion'])
def test_unproved_promotions_fail(result, flag):
    altered = deepcopy(result)
    altered['spectrum'][flag] = True
    with pytest.raises(ValueError, match='beyond proved scope'):
        packet.validate(altered)


def test_missing_weak_line_and_deleted_energy_control_fail(result):
    altered = deepcopy(result)
    altered['spectrum']['rows'].pop()
    with pytest.raises(ValueError, match='coverage'):
        packet.validate(altered)
    altered = deepcopy(result)
    altered['spectrum']['independent_sum_rules']['independent_of_composite_eigenvectors'] = False
    with pytest.raises(ValueError, match='energy moment'):
        packet.validate(altered)


def test_retained_receipt_and_table_match_exact_replay(result):
    retained = packet.strict_json.loads(packet.OUTPUT.read_text())
    packet.validate(retained)
    assert packet.canonical(result) == packet.canonical(retained)
    assert packet.table(result) == packet.TABLE.read_text()
