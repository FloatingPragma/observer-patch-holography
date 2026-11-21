"""Bind the implemented instrument and its independent forecast before opening records."""
from hashlib import sha256

from tools.strict_json import load
from model.candidate_laws.interval_echo import candidate, clock
from model.candidate_laws.contrast_memory_spherical_v1 import law
from .gate import ROOT, require
from .preparation import HERE


def sources():
    files = [HERE/name for name in ("clocked_protocol.json","clocked_field.py","clocked_apparatus.py",
             "clocked_run.py","receiver_sources.json","preparation_sources.json")]
    files.append(ROOT/"model/candidate_laws/interval_echo/clock.py")
    return {p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n",b"\n")).hexdigest()
            for p in files}


def verify_sources():
    require(candidate.KAPPA == clock.KAPPA == 8 and clock.law is candidate.native_law is law,
            "clock and data need the same native generator")
    require(sources()==load(HERE/"clocked_sources.json"), "changed clocked field sources")
