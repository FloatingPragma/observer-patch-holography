"""The sealed supplement keeps its parent class, adoption and claim scope."""
from copy import deepcopy
from unittest import mock

import pytest

from tools import check_field_record_addenda_claims as claims


@pytest.fixture
def data():
    return claims.strict_json.load(claims.ROOT/'claims/addenda/field_record_addenda.json')


def test_current_class_is_equivalent_without_rewriting_the_seal(data):
    before = deepcopy(data)
    claims.check_parent_class(data)
    assert data == before
    assert data['sealed_program_class_sha256'] == claims.SEALED_CLASS_SHA256


@pytest.mark.parametrize('field', ('sealed_parent_commit', 'sealed_program_class_sha256'))
def test_replacing_the_sealed_authority_is_rejected(data, field):
    data[field] = '0'*len(data[field])
    with pytest.raises(ValueError, match='Sealed parent authority'):
        claims.check_parent_class(data)


def test_changed_adoption_semantics_are_rejected(data):
    klass = claims.strict_json.load(claims.ROOT/claims.program_class_v2.PATH)
    klass['adoption']['status'] = 'PROPOSED_FOR_INDEPENDENT_REVIEW'
    with mock.patch.object(claims.strict_json, 'load', side_effect=lambda path: (
            deepcopy(klass) if str(path).endswith('program_class_v2.json') else claims.strict_json.loads(path.read_text()))):
        with pytest.raises(ValueError, match='adoption semantics'):
            claims.check_parent_class(data)


@pytest.mark.parametrize('field', ('changes_existing_physical_statuses', 'physical_particle_identifications',
                                   'source_derived_masses_GeV'))
def test_physical_promotion_is_rejected(data, field):
    data[field] = True
    with mock.patch.object(claims.strict_json, 'load', return_value=data):
        with pytest.raises(ValueError, match='physical promotion'):
            claims.check()


def test_current_supplement_passes(capsys):
    claims.check()
    assert 'PASS' in capsys.readouterr().out
