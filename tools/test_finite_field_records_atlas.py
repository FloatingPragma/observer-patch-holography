"""A PNG check must detect stale pixels without silently rewriting them."""
import base64
import hashlib
from types import SimpleNamespace

import pytest

from tools import build_finite_field_records_atlas as atlas


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="1" height="1"/>\n'
    monkeypatch.setattr(atlas, 'OUT', tmp_path / 'atlas.svg')
    monkeypatch.setattr(atlas, 'render', lambda: svg)
    raw = base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=')

    def renderer(command, **kwargs):
        assert command[:2] == ['rsvg-convert', '--format=png']
        assert kwargs['input'] == svg.encode()
        return SimpleNamespace(stdout=raw)

    monkeypatch.setattr(atlas.subprocess, 'run', renderer)
    return svg


def test_written_png_binds_svg_and_check_is_read_only(prepared):
    atlas.main(['--png'])
    png = atlas.OUT.with_suffix('.png')
    expected = png.read_bytes()
    assert hashlib.sha256(prepared.encode()).hexdigest().encode() in expected
    before = (atlas.OUT.stat().st_mtime_ns, png.stat().st_mtime_ns)
    atlas.main(['--check', '--png'])
    assert png.read_bytes() == expected
    assert before == (atlas.OUT.stat().st_mtime_ns, png.stat().st_mtime_ns)


@pytest.mark.parametrize('mutation', ['missing', 'stale'])
def test_png_check_rejects_stale_artifact_without_repair(prepared, mutation):
    atlas.main(['--png'])
    png = atlas.OUT.with_suffix('.png')
    if mutation == 'missing':
        png.unlink()
    else:
        png.write_bytes(b'stale raster')
    with pytest.raises(SystemExit, match='PNG is stale'):
        atlas.main(['--check', '--png'])
    assert (not png.exists()) if mutation == 'missing' else png.read_bytes() == b'stale raster'


def test_svg_check_rejects_stale_source_without_repair(prepared):
    atlas.OUT.write_text('stale vector')
    with pytest.raises(SystemExit, match='atlas is stale'):
        atlas.main(['--check', '--png'])
    assert atlas.OUT.read_text() == 'stale vector'
