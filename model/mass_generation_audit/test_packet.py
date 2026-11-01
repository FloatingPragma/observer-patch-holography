from copy import deepcopy

import pytest

from . import packet


def test_packet_exact_scopes_and_transitive_sources():
    value = packet._produce()
    assert value['scope']['positive_canonical_charge_gap_certified']
    pins = value['source_bindings']
    assert 'model/periodic_momentum/momentum.py' in pins
    assert 'model/common_source/mass_response/gap_certificate.py' in pins
    assert 'paper/tex_fragments/MASS_GENERATION_AUDIT.tex' in pins


@pytest.mark.parametrize('key', ['massive_three_dimensional_species_identified',
                                 'native_ground_preparation_executed',
                                 'complete_repair_ground_invariance_proved',
                                 'empirical_mass_fit_used'])
def test_no_receipt_promotion(key):
    value = deepcopy(packet._produce())
    value['scope'][key] = True
    with pytest.raises(ValueError): packet.validate(value)


def test_import_time_custody_change_is_rejected(monkeypatch):
    pins = packet.source_bindings()
    pins['model/mass_generation_audit/gap.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: pins)
    with pytest.raises(ValueError, match='source changed since import'):
        packet._produce()
