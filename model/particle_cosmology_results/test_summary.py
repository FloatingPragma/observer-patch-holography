import json
from pathlib import Path

import pytest

from .summary import authenticate, certificate, digest, require_false


def test_actual_source_join_and_mass_ancestry():
    result = certificate()
    assert len(result['tasks']) == 10
    assert len(result['conditional_GeV_targets']) == 17
    assert result['physical_premise_promotions'] == []


def test_modified_dependency_and_forged_success_rejected(tmp_path: Path):
    directory = tmp_path / 'model' / 'control'
    directory.mkdir(parents=True)
    source = directory / 'source.py'
    source.write_text('fixed_source = 1\n')
    record = {'schema': 'oph.control.v1', 'all_checks_pass': True,
              'source_bindings': {'model/control/source.py': digest(source)}}
    path = directory / 'results.json'
    path.write_text(json.dumps(record))
    authenticate('control', tmp_path)
    source.write_text('fixed_source = 2\n')
    with pytest.raises(ValueError, match='changed input'):
        authenticate('control', tmp_path)
    record['all_checks_pass'] = False
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='failed scientific checks'):
        authenticate('control', tmp_path)
    del record['all_checks_pass']
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='failed scientific checks'):
        authenticate('control', tmp_path)


def test_physical_promotion_cannot_follow_from_receipt_success():
    require_false(False, 'control')
    for forged in (True, None, 0, 'false'):
        with pytest.raises(ValueError, match='Unsupported physical promotion'):
            require_false(forged, 'control')
