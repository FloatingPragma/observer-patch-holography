import hashlib
import json

import pytest

from .verify import verify_manifest


def fixture(tmp_path):
    source = tmp_path / 'input.json'
    source.write_text('{"sample": [1, 2, 3]}\n')
    raw = source.read_bytes()
    row = {'path': 'input.json', 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
    manifest = tmp_path / 'input_manifest.json'
    manifest.write_text(json.dumps({'schema': 'oph.scientific_comparison_inputs.v1', 'files': [row]}))
    return source, manifest, row


def test_content_change_rejected(tmp_path):
    source, manifest, row = fixture(tmp_path)
    assert verify_manifest(manifest)['files'] == 1
    source.write_text('{"sample": [1, 2, 4]}\n')
    with pytest.raises(ValueError, match='bytes differ'):
        verify_manifest(manifest)


@pytest.mark.parametrize('name', ['../input.json', '/input.json'])
def test_escape_rejected(tmp_path, name):
    source, manifest, row = fixture(tmp_path)
    row['path'] = name
    manifest.write_text(json.dumps({'schema': 'oph.scientific_comparison_inputs.v1', 'files': [row]}))
    with pytest.raises(ValueError, match='escaping'):
        verify_manifest(manifest)


def test_duplicate_rejected(tmp_path):
    source, manifest, row = fixture(tmp_path)
    manifest.write_text(json.dumps({'schema': 'oph.scientific_comparison_inputs.v1', 'files': [row, row]}))
    with pytest.raises(ValueError, match='Duplicate'):
        verify_manifest(manifest)


def test_nonfinite_input_requires_explicit_ineligibility(tmp_path):
    source, manifest, row = fixture(tmp_path)
    source.write_text('{"unresolved_value": NaN}\n')
    row.update(bytes=source.stat().st_size, sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    value = {'schema': 'oph.scientific_comparison_inputs.v1', 'files': [row]}
    manifest.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='ineligible numerical classification'):
        verify_manifest(manifest)
    row.update(nonfinite_constant_count=1, numerical_status='raw_nonfinite_values',
               eligible_for_automatic_prediction_or_passing_numeric_receipt=False)
    manifest.write_text(json.dumps(value))
    report = verify_manifest(manifest)
    assert report['nonstandard_json_constants_retained'] == [
        {'path': 'input.json', 'nonstandard_constants': 1,
         'eligible_for_prediction_or_passing_numeric_receipt': False}]
    assert b'NaN' in source.read_bytes()
