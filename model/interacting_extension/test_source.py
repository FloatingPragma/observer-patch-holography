"""A commit label cannot authenticate modified source bytes or a stale cache."""
from pathlib import Path
import shutil

import pytest

from ._source import ROOT, verify_sources
from .dependencies import source_reference

NAME='model/observer_particles/results/unified_yukawa_operators.json'


def copy_input(tmp_path):
    p=tmp_path/NAME;p.parent.mkdir(parents=True)
    shutil.copyfile(ROOT/NAME,p)
    return p


def test_source_contract_works_without_git_and_rejects_dirty_bytes(tmp_path):
    p=copy_input(tmp_path)
    assert verify_sources(tmp_path,[NAME])==verify_sources(ROOT,[NAME])
    p.write_bytes(p.read_bytes()+b' ')
    with pytest.raises(ValueError,match='source bytes'):verify_sources(tmp_path,[NAME])


def test_cached_exact_data_cannot_bypass_revalidation(tmp_path,monkeypatch):
    p=copy_input(tmp_path);sentinel=object()
    monkeypatch.setattr(source_reference,'SOURCE_FILES',(NAME,))
    monkeypatch.setattr(source_reference,'_exact_data',lambda path:sentinel)
    assert source_reference.exact_source_packet(tmp_path) is sentinel
    p.write_bytes(p.read_bytes()+b' ')
    with pytest.raises(ValueError,match='source bytes'):source_reference.exact_source_packet(tmp_path)


def test_missing_unknown_empty_and_escaping_sources_rejected(tmp_path):
    for names in ([],[NAME],['../outside']):
        with pytest.raises(ValueError):verify_sources(tmp_path,names)
    p=tmp_path/NAME;p.parent.mkdir(parents=True);p.symlink_to(ROOT/NAME)
    with pytest.raises(ValueError,match='escaping'):verify_sources(tmp_path,[NAME])
