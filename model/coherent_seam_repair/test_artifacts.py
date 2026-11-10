"""Custody for the complete finite construction and its proof artifacts."""
import hashlib
import json
from pathlib import Path


def test_package_manifest_binds_every_artifact():
    root = Path(__file__).resolve().parent
    expected = json.loads((root / 'manifest.json').read_text())['files']
    actual = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob('*')
        if path.is_file() and path.name != 'manifest.json'
        and not {'__pycache__', '.pytest_cache'}.intersection(path.relative_to(root).parts)
    }
    assert expected == actual
